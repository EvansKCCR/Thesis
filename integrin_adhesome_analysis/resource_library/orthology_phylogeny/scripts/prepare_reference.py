"""Prepare one complete reference proteome using explicit protein-to-gene mapping."""

# Resolve source resources independently of the caller working directory.
import sys as _layout_sys
from pathlib import Path as _LayoutPath
_LAYOUT_BUNDLE_ROOT = next(p for p in _LayoutPath(__file__).resolve().parents
                          if (p/'config/layout_paths.json').is_file() and (p/'pipeline_paths.py').is_file())
if str(_LAYOUT_BUNDLE_ROOT) not in _layout_sys.path:
    _layout_sys.path.insert(0, str(_LAYOUT_BUNDLE_ROOT))
from pipeline_paths import activate as _layout_activate, resource_path as _layout_resource
_layout_activate()

import argparse
import hashlib
import json
import re
import shutil
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path
from pipeline import records, digest, read_tsv, write_tsv, AA


BENCHMARK_SEQUENCE_POLICY = 'remove terminal stop; mask internal stop markers as X; apply configured sequence QC; retain ID provenance'
DEFAULT_QC = dict(min_length=30, max_ambiguous_fraction=0.05, max_ambiguous_run=20)
AMBIGUOUS = set('XBZJ')


def normalize_sequence(seq, protein_id, *, internal_stop_policy='mask'):
    """Normalize stops, leaving invalid-residue decisions to the QC report."""
    if internal_stop_policy not in ('reject', 'mask'):
        raise ValueError('Unknown internal stop policy: '+internal_stop_policy)
    seq = seq.upper().removesuffix('*')
    if internal_stop_policy == 'mask':
        seq = seq.replace('*', 'X')
    return seq


def infer_gene(header, protein_id, gene_regex=None):
    """Use only explicit header tags or conservative common isoform suffixes."""
    if gene_regex:
        match = re.search(gene_regex, header)
        if match and match.groups() and match.group(1):
            return match.group(1), 'configured_gene_regex', None
    match = re.search(r'(?:^|[\s|;])gene(?:_id)?[=:]([^|\s;]+)', header, re.I)
    if match:
        return match.group(1), 'header_gene_tag', None
    patterns = (r'^(.*?)[._-](?:isoform|transcript|trans|variant|t|i)(\d+)$',
                r'^(.*?)[._-](\d+)$', r'^(.*?)-R([A-Z]+)$')
    for pattern in patterns:
        match = re.match(pattern, protein_id, re.I)
        if match and match.group(1):
            return match.group(1), 'suffix_inference', protein_id
    return None, 'sequence_cluster_fallback', None


def _mapping_table(mapping):
    if mapping is None:
        return {}
    fields, rows = read_tsv(mapping)
    if not {'protein_id', 'gene_id'} <= set(fields or []):
        raise ValueError('protein_gene.tsv requires protein_id and gene_id columns')
    result = {}
    for row in rows:
        pid, gene = (row.get('protein_id') or '').strip(), (row.get('gene_id') or '').strip()
        if not pid or not gene or pid in result or not re.fullmatch(r'[A-Za-z0-9_.-]+', gene):
            raise ValueError(f'Duplicate or invalid protein_gene.tsv row: {pid}')
        result[pid] = gene
    return result


def _cdhit_groups(species, fallback, out, executable):
    if not fallback:
        return {}
    command_path = shutil.which(executable)
    if not command_path:
        raise ValueError(f'CD-HIT is required for unmapped proteins without usable header/suffix gene IDs; install `{executable}` or supply protein_gene.tsv')
    with tempfile.TemporaryDirectory(prefix=f'.{species}_cdhit_', dir=out) as temp_name:
        temp = Path(temp_name)
        source, target = temp/'fallback.faa', temp/'clusters.faa'
        with source.open('w', encoding='ascii') as stream:
            for pid, seq in sorted(fallback.items()):
                stream.write(f'>{pid}\n{seq}\n')
        command = [command_path, '-i', str(source), '-o', str(target), '-c', '0.98', '-aS', '0.90', '-d', '0', '-T', '1', '-M', '0']
        result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        (out/'cdhit_fallback.log').write_text('COMMAND: '+' '.join(command)+'\n'+result.stdout, encoding='utf-8')
        if result.returncode:
            raise ValueError(f'CD-HIT fallback failed with exit status {result.returncode}; inspect {out / "cdhit_fallback.log"}')
        clstr = Path(str(target)+'.clstr')
        groups, cluster_index = {}, -1
        for line in clstr.read_text(encoding='utf-8').splitlines():
            if line.startswith('>Cluster '):
                cluster_index += 1
            else:
                member = re.search(r'>(.*?)\.\.\.', line)
                if member:
                    groups[member.group(1)] = f'cdhit_{cluster_index:06d}'
        if set(groups) != set(fallback):
            raise ValueError('CD-HIT cluster report does not account for every fallback sequence')
        shutil.copyfile(clstr, out/'cdhit_fallback.clstr')
        return groups


def prepare(species, source, mapping=None, out=None, *, internal_stop_policy='mask', gene_regex=None,
            qc=None, cdhit='cd-hit'):
    if not re.fullmatch(r'[A-Za-z][A-Za-z0-9_]*', species) or '__' in species:
        raise ValueError('Use a species code starting with a letter and no double underscores')
    if out is None:
        raise ValueError('Output directory is required')
    qc = dict(DEFAULT_QC, **(qc or {}))
    if qc['min_length'] < 1 or not 0 <= qc['max_ambiguous_fraction'] <= 1 or qc['max_ambiguous_run'] < 0:
        raise ValueError('Invalid sequence QC thresholds')
    out.mkdir(parents=True, exist_ok=False)
    mapping_rows = _mapping_table(mapping)
    raw_by_id = defaultdict(list)
    for index, (header, raw_seq) in enumerate(records(source), 1):
        pid = header.split()[0]
        seq = normalize_sequence(raw_seq, pid, internal_stop_policy=internal_stop_policy)
        raw_by_id[pid].append(dict(protein_id=pid,header=header,raw_sequence=raw_seq,sequence=seq,index=index))
    if not raw_by_id:
        raise ValueError('Empty reference proteome')
    unique, qc_rows = {}, []
    for pid, copies in raw_by_id.items():
        ordered = sorted(copies, key=lambda row:(-len(row['sequence']),row['sequence'],row['header'],row['index']))
        chosen = ordered[0]
        invalid = sorted(set(chosen['sequence'])-AA)
        ambiguous = sum(char in AMBIGUOUS for char in chosen['sequence'])
        fraction = ambiguous/len(chosen['sequence']) if chosen['sequence'] else 1.0
        longest_run = max((len(match.group()) for match in re.finditer(r'[XBZJ]+',chosen['sequence'])),default=0)
        reasons=[]
        if invalid: reasons.append('invalid_residue:'+''.join(invalid))
        if len(chosen['sequence']) < qc['min_length']: reasons.append('below_minimum_length')
        if fraction > qc['max_ambiguous_fraction']: reasons.append('ambiguous_fraction')
        if longest_run > qc['max_ambiguous_run']: reasons.append('ambiguous_run')
        identical_duplicates=sum(copy['sequence']==chosen['sequence'] for copy in copies)-1
        conflict_duplicates=len(copies)-1-identical_duplicates
        if conflict_duplicates: reasons.append('duplicate_id_conflict_resolved_deterministically')
        decision=';'.join(reasons) if reasons else 'qc_pass'
        row=dict(**chosen,ambiguous_fraction=fraction,longest_ambiguous_run=longest_run,
                 invalid_residues=''.join(invalid),duplicate_record_count=len(copies),
                 identical_duplicate_count=identical_duplicates,
                 terminal_stop_removed=int(chosen['raw_sequence'].endswith('*')),
                 internal_stops_masked=chosen['raw_sequence'].removesuffix('*').count('*'),decision=decision)
        unique[pid]=row
        qc_rows.append(row)
    eligible={pid:row for pid,row in unique.items() if not any(reason in row['decision'] for reason in
        ('invalid_residue','below_minimum_length','ambiguous_fraction','ambiguous_run'))}
    explicit = {}
    fallback = {}
    suffix_candidates = defaultdict(list)
    sources = {}
    for pid,row in eligible.items():
        if mapping is not None:
            gene=mapping_rows.get(pid)
            if gene is None:
                raise ValueError(f'{pid} is absent from the explicitly specified protein_gene.tsv')
            explicit[pid]=gene;sources[pid]='protein_gene.tsv'
            continue
        gene,method,candidate=infer_gene(row['header'],pid,gene_regex)
        if method=='suffix_inference' and candidate:
            suffix_candidates[gene].append(pid)
        elif gene:
            explicit[pid]=gene;sources[pid]=method
        else:
            fallback[pid]=row['sequence']
    genes={}
    for candidate,pids in suffix_candidates.items():
        if len(pids)>1:
            for pid in pids: genes[pid]=candidate;sources[pid]='shared_isoform_suffix'
        else:
            pid=pids[0];fallback[pid]=eligible[pid]['sequence']
    cluster_map=_cdhit_groups(species,fallback,out,cdhit)
    for pid,gene in explicit.items(): genes[pid]=gene
    for pid,cluster in cluster_map.items(): genes[pid]=cluster;sources[pid]='cd-hit_98pct_short_coverage_90pct'
    groups=defaultdict(list)
    for pid,row in eligible.items(): groups[genes[pid]].append((pid,row['sequence']))
    if not groups:
        raise ValueError('No proteins passed sequence-level QC')
    target = _layout_resource(out/f'{species}.faa')
    audit = []
    selected={}
    with _layout_resource(target).open('w', encoding='ascii') as f:
        for gene, variants in sorted(groups.items()):
            chosen = sorted(variants, key=lambda x:(-len(x[1]),x[0]))[0]
            prepared = species+'__'+gene
            f.write(f'>{prepared}\n{chosen[1]}\n')
            selected[gene]=(chosen[0],prepared)
    for pid,row in unique.items():
        gene=genes.get(pid,'')
        representative,prepared=selected[gene] if gene in selected else ('','')
        decision='retained' if representative==pid else row['decision'] if row['decision']!='qc_pass' else 'not_selected_longest_isoform_or_cluster_representative'
        audit.append(dict(protein_id=pid,original_id=pid,original_header=row['header'],gene_id=gene,
            prepared_id=prepared,representative_id=representative,selected=int(representative==pid),
            length=len(row['sequence']),ambiguous_fraction=row['ambiguous_fraction'],
            longest_ambiguous_run=row['longest_ambiguous_run'],invalid_residues=row['invalid_residues'],
            terminal_stop_removed=row['terminal_stop_removed'],internal_stops_masked=row['internal_stops_masked'],
            duplicate_record_count=row['duplicate_record_count'],identical_duplicate_count=row['identical_duplicate_count'],
            gene_assignment_method=sources.get(pid,''),decision=decision))
    fields=list(audit[0])
    write_tsv(_layout_resource(out/'identifier_map.tsv'), fields, audit)
    write_tsv(_layout_resource(out/'sequence_qc.tsv'), fields, audit)
    (_layout_resource(out/'preparation.json')).write_text(json.dumps(dict(pipeline_version='1.0',species=species,
        input_protein_ids=len(raw_by_id),input_fasta_records=sum(len(x) for x in raw_by_id.values()),
        genes=len(groups),sequence_qc_thresholds=qc,
        policy='QC before deterministic representative selection; mapped gene IDs when specified; otherwise header/suffix inference and CD-HIT fallback',
        internal_stop_policy=internal_stop_policy,
        sequence_policy=BENCHMARK_SEQUENCE_POLICY,
        sequence_qc=dict(retained=sum(r['decision']=='retained' for r in audit),
            excluded=sum(r['decision']!='retained' for r in audit),
            invalid_protein_ids=sum(bool(r['invalid_residues']) for r in audit),
            internal_stop_proteins=sum(r['internal_stops_masked']>0 for r in audit),
            internal_stop_sites=sum(r['internal_stops_masked'] for r in audit),
            duplicate_ids=sum(r['duplicate_record_count']>1 for r in audit),
            cdhit_fallback_proteins=sum(r['gene_assignment_method']=='cd-hit_98pct_short_coverage_90pct' for r in audit),
            report_sha256=digest(_layout_resource(out/'sequence_qc.tsv'))),
        input_fasta=str(_layout_resource(source).resolve()),input_sha256=digest(source),
        mapping=str(_layout_resource(mapping).resolve()) if mapping is not None else None,
        mapping_sha256=digest(mapping) if mapping is not None else None,
        gene_regex=gene_regex, cdhit_thresholds=dict(identity=0.98,shorter_sequence_coverage=0.90) if fallback else None,
        output_sha256=digest(target),unused_mapping_rows=len(set(mapping_rows)-set(raw_by_id))),indent=2)+'\n')
    return target


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--version',action='version',version='adhesome-gene-preparation 1.0')
    p.add_argument('--species',required=True)
    p.add_argument('--fasta',type=Path,required=True)
    p.add_argument('--map',dest='mapping',type=Path,help='Optional protein_gene.tsv; omit to infer from headers/suffixes then use CD-HIT fallback')
    p.add_argument('--gene-regex',help='Optional regex with first capture group as gene ID')
    p.add_argument('--min-length',type=int,default=DEFAULT_QC['min_length'])
    p.add_argument('--max-ambiguous-fraction',type=float,default=DEFAULT_QC['max_ambiguous_fraction'])
    p.add_argument('--max-ambiguous-run',type=int,default=DEFAULT_QC['max_ambiguous_run'])
    p.add_argument('--cdhit',default='cd-hit')
    p.add_argument('--out',type=Path,required=True,help='New output directory')
    a=p.parse_args()
    try:
        print(prepare(a.species,a.fasta,a.mapping,a.out,gene_regex=a.gene_regex,
            qc=dict(min_length=a.min_length,max_ambiguous_fraction=a.max_ambiguous_fraction,max_ambiguous_run=a.max_ambiguous_run),cdhit=a.cdhit))
    except (ValueError,OSError,KeyError) as exc:
        p.exit(1,f'ERROR: {exc}\n')


if __name__=='__main__':
    main()
