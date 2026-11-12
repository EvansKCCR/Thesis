"""Retrieve strictly reciprocal one-to-one validation negative homologues."""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys

from pipeline_paths import resource_path as _layout_resource
from .datasets import sha, table
from .reports import read

ROOT = Path(__file__).resolve().parents[2]
HOMOLOGY = ROOT.parent / 'homology_retrieval' / 'reciprocal_blast.py'


def fasta_records(path):
    header = None
    chunks = []
    path=_layout_resource(path)
    opener=gzip.open if path.suffix.casefold()=='.gz' else open
    with opener(path,'rt',encoding='utf-8-sig') as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(chunks)
                header, chunks = line[1:], []
            else:
                if header is None:
                    raise ValueError('Sequence before FASTA header: '+str(path))
                chunks.append(re.sub(r'\s+', '', line).upper())
    if header is not None:
        yield header, ''.join(chunks)


def discover_benchmark_input(reference_root, explicit=None):
    if explicit is not None:
        path = _layout_resource(Path(explicit)).resolve()
        if not path.is_dir():
            raise ValueError('Benchmark protein directory does not exist: '+str(path))
        return path
    root = _layout_resource(Path(reference_root)).resolve()
    candidates = [p for p in root.rglob('*') if p.is_dir() and p.name.casefold() == 'input'
                  and p.parent.name.casefold() in ('benchmark', 'benchmarks')]
    if len(candidates) > 1:
        raise ValueError('Multiple Benchmark/Input directories found; set the validation benchmark input path explicitly: '
                         +', '.join(str(p) for p in sorted(candidates)))
    return candidates[0] if candidates else None


def load_rows(path):
    with _layout_resource(path).open(encoding='utf-8-sig', newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def proteome_label(path):
    name=path.name
    if name.casefold().endswith('.gz'):
        name=name[:-3]
    for ext in ('.fasta','.fsa','.faa','.fas','.pep','.aa','.fa'):
        if name.casefold().endswith(ext):
            name=name[:-len(ext)]
            break
    return re.sub(r'[^A-Za-z0-9_.-]','_',name)


def run(negative_controls, reference_root, out, *, threads=1, benchmark_input=None, qc=None):
    """Run the existing RBH tool and append only fully checked 1:1 RBH controls.

    ``negative_controls`` are already labelled human controls from the frozen
    validation input. The resulting panel is additive and carries its exact
    source IDs and RBH audit fields; ambiguous/one-to-many matches are omitted.
    """
    out = _layout_resource(Path(out)).resolve()
    panel = discover_benchmark_input(reference_root, benchmark_input)
    if panel is None:
        return [], dict(status='unavailable', reason='no Benchmark/Input directory under the configured raw reference root',
                        confirmed_one_to_one=0)
    queries = [r for r in negative_controls if r.get('control_class') != 'positive']
    if not queries:
        return [], dict(status='not_applicable', reason='validation input has no negative controls',
                        confirmed_one_to_one=0)
    if not HOMOLOGY.is_file():
        raise ValueError('Missing homology_retrieval pipeline: '+str(HOMOLOGY))
    query_dir = _layout_resource(out/'query')
    query_dir.mkdir(parents=True, exist_ok=False)
    query_fasta = _layout_resource(query_dir/'homo_sapiens_negative_controls.faa')
    grouped_queries={}
    for row in queries:
        grouped_queries.setdefault(row['sequence_sha256'],[]).append(row)
    unique_queries=[rows[0] for rows in grouped_queries.values()]
    query_rows=[]
    with query_fasta.open('w',encoding='ascii',newline='\n') as handle:
        for index,row in enumerate(unique_queries,1):
            query_id=f'candidate|NEG{index:07d}'
            handle.write(f'>{query_id}\n{row["sequence"]}\n')
            source_ids=[item['original_id'] for item in grouped_queries[row['sequence_sha256']]]
            query_rows.append(dict(query_id=query_id,source_control_id=row['original_id'],
                source_control_ids=';'.join(source_ids),
                source_accession=row.get('accession',''),sequence_sha256=row['sequence_sha256'],
                control_class=row.get('control_class','negative')))
    source_map=[]
    for row in query_rows:
        source_map.append(dict(reference_file=query_fasta.name,reference_id=row['query_id'],
            reference_proteome=query_fasta.name,source_control_id=row['source_control_id'],
            source_control_ids=row['source_control_ids'],source_accession=row['source_accession'],
            sequence_sha256=row['sequence_sha256']))
    table(_layout_resource(query_dir/'query_source_map.tsv'),source_map,
          ['reference_file','reference_id','reference_proteome','source_control_id',
           'source_control_ids','source_accession','sequence_sha256'])
    output = _layout_resource(out/'homology_run')
    command=[sys.executable,'-B',str(HOMOLOGY),'--reference-dir',str(query_dir),
             '--target-dir',str(panel),'--background-dir',str(query_dir),
             '--source-map',str(query_dir/'query_source_map.tsv'),
             '--outdir',str(output),'--threads',str(max(1,int(threads)))]
    completed=subprocess.run(command,cwd=ROOT,check=False)
    if completed.returncode:
        raise ValueError('Negative-control homology retrieval failed; inspect '+str(output/'logs'))
    query_by_id={r['query_id']:r for r in query_rows}
    qc=qc or {'min_length':30,'max_ambiguous_fraction':0.05,'max_ambiguous_run':20}
    min_length=int(qc.get('min_length',30));max_ambiguous=float(qc.get('max_ambiguous_fraction',0.05))
    max_run=int(qc.get('max_ambiguous_run',20));ambiguous=set('BXZJUO')

    def quality_failure(sequence):
        invalid=set(sequence)-set('ACDEFGHIKLMNPQRSTVWYBXZJUO')
        if invalid:return 'invalid_amino_acids:'+''.join(sorted(invalid))
        if len(sequence)<min_length:return 'below_min_length'
        fraction=sum(residue in ambiguous for residue in sequence)/len(sequence)
        if fraction>max_ambiguous:return 'ambiguous_fraction_exceeded'
        run=longest=0
        for residue in sequence:
            run=run+1 if residue in ambiguous else 0
            longest=max(longest,run)
        if longest>max_run:return 'ambiguous_stretch_exceeded'
        return ''
    summaries=load_rows(_layout_resource(output/'targets_summary.tsv'))
    target_files={proteome_label(path):path for path in panel.rglob('*')
                  if path.is_file() and path.name.casefold().endswith((
                      '.fa','.faa','.fasta','.fas','.fsa','.pep','.aa','.fa.gz','.faa.gz','.fasta.gz','.fas.gz','.fsa.gz','.pep.gz','.aa.gz'))}
    if len(target_files)!=len([p for p in panel.rglob('*') if p.is_file() and p.name.casefold().endswith((
            '.fa','.faa','.fasta','.fas','.fsa','.pep','.aa','.fa.gz','.faa.gz','.fasta.gz','.fas.gz','.fsa.gz','.pep.gz','.aa.gz'))]):
        raise ValueError('Benchmark proteome filenames do not map uniquely to homology species labels')
    if {row['species'] for row in summaries}!=set(target_files):
        raise ValueError('Homology tool target files differ from the selected Benchmark/Input panel')
    target_sequences={};target_sequence_counts={}
    audit=load_rows(_layout_resource(output/'homology_audit.tsv'))
    eligible=[]
    for row in audit:
        query=query_by_id.get(row.get('reference_id',''))
        row['source_control_id']=query['source_control_id'] if query else ''
        row['source_control_ids']=query['source_control_ids'] if query else ''
        row['source_control_accession']=query['source_accession'] if query else ''
        row['source_sequence_sha256']=query['sequence_sha256'] if query else ''
        same_source=(row.get('mapping_method')=='identical_sequence'
                     and row.get('source_assignment')=='explicit_source_map'
                     and 'homo_sapiens' in row.get('reference_proteome','').casefold())
        nonhuman='homo_sapiens' not in row.get('species','').casefold()
        eligible.append(row.get('status')=='RBH_candidate' and same_source and nonhuman)
        row['one_to_one_candidate']=False
        row['one_to_one_decision']='not_verified_reciprocal_or_exact_source_mapping'
    candidate_indexes=[i for i,ok in enumerate(eligible) if ok]
    target_counts={}
    query_counts={}
    for i in candidate_indexes:
        row=audit[i]
        target_counts[(row['species'],row['best_target_id'])]=target_counts.get((row['species'],row['best_target_id']),0)+1
        query_counts[(row['species'],row['reference_id'])]=query_counts.get((row['species'],row['reference_id']),0)+1
    added=[]
    for i in candidate_indexes:
        row=audit[i]
        if target_counts[(row['species'],row['best_target_id'])]!=1 or query_counts[(row['species'],row['reference_id'])]!=1:
            row['one_to_one_decision']='reciprocal_hit_not_one_to_one_in_panel'
            continue
        target_path=target_files.get(row['species'])
        if target_path is None:
            raise ValueError('Homology audit target is absent from targets_summary.tsv: '+row['species'])
        if target_path not in target_sequences:
            target_sequences[target_path]={header.split()[0]:sequence for header,sequence in fasta_records(target_path)}
            counts={}
            for sequence in target_sequences[target_path].values():
                digest=hashlib.sha256(sequence.encode('ascii')).hexdigest()
                counts[digest]=counts.get(digest,0)+1
            target_sequence_counts[target_path]=counts
        target_id=row['best_target_id']
        sequence=target_sequences[target_path].get(target_id)
        if not sequence:
            raise ValueError('Confirmed RBH target sequence not found: '+row['species']+'/'+target_id)
        sequence_digest=hashlib.sha256(sequence.encode('ascii')).hexdigest()
        failure=quality_failure(sequence)
        if failure:
            row['one_to_one_decision']='target_sequence_failed_quality_control:'+failure
            continue
        if target_sequence_counts[target_path].get(sequence_digest)!=1:
            row['one_to_one_decision']='target_sequence_not_unique_in_reference_proteome'
            continue
        existing=[control for control in negative_controls if control.get('sequence_sha256')==sequence_digest]
        if existing:
            row['one_to_one_decision']='identical_sequence_already_in_validation_control_panel'
            continue
        token=hashlib.sha256((row['species']+'\0'+target_id).encode()).hexdigest()[:16]
        original=f'BenchmarkNegative|{token}'
        row['one_to_one_candidate']=True
        row['one_to_one_decision']='confirmed_unique_reciprocal_best_hit_to_exact_human_source'
        added.append(dict(original_id=original,species=row['species'].replace('_',' '),control_class='negative',
            accession='',genes='',protein_name='',trusted_family='',sequence=sequence,
            sequence_sha256=hashlib.sha256(sequence.encode('ascii')).hexdigest(),length=len(sequence),
            negative_homology_source='human_negative_control_panel',
            negative_homology_source_id=row['source_control_id'],
            negative_homology_source_ids=row['source_control_ids'],
            negative_homology_source_accession=row['source_control_accession'],
            negative_homology_reference_query_id=row['reference_id'],
            negative_homology_target_id=target_id,negative_homology_species=row['species'],
            negative_homology_status=row['status'],negative_homology_one_to_one=True,
            negative_homology_evidence=';'.join(f'{k}={row.get(k,"")}' for k in (
                'status','reference_proteome','mapped_background_id','best_target_id','reverse_best_id',
                'forward_pident','forward_qcov','forward_scov','reverse_pident','reverse_qcov','reverse_scov'))))
    table(_layout_resource(out/'homology_matches.tsv'),audit)
    table(_layout_resource(out/'confirmed_one_to_one_negative_homologues.tsv'),
          [row for row in audit if row['one_to_one_candidate']],
          list(audit[0]) if audit else ['reference_id','species','status','one_to_one_candidate','one_to_one_decision'])
    decisions={}
    for row in audit:
        key=row.get('one_to_one_decision','')
        decisions[key]=decisions.get(key,0)+1
    summary=dict(status='complete',benchmark_input=str(panel),query_controls=len(queries),
        distinct_query_sequences=len(query_rows),duplicate_source_controls_collapsed=len(queries)-len(query_rows),
        all_reciprocal_and_nonambiguous=sum(eligible),confirmed_one_to_one=len(added),
        rejected_or_unverified=len(audit)-sum(bool(r['one_to_one_candidate']) for r in audit),
        decisions=decisions,quality_control=qc,
        source_mapping='exact identity to each labelled Homo sapiens negative-control sequence in the generated source background',
        confirmation='RBH_candidate + exact source-sequence mapping + unique query/target pair + unique target sequence + configured sequence QC',
        unverified_matches_used=False)
    _layout_resource(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
    return added,summary
