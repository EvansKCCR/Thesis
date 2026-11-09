#!/usr/bin/env python3
"""Prepare/evaluate a reference-only orthogroup reconstruction test.

The inference directory contains reference proteomes only. Positive controls are
exported to a separate query FASTA and are read only after OrthoFinder finishes.
Evaluation uses exact full-length amino-acid sequence identity, never IDs.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
from pathlib import Path
import re
import shutil
import sys


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda:handle.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def fasta(path):
    records={};identifier=None;parts=[]
    def save():
        seq=''.join(parts).upper().removesuffix('*')
        if not identifier or not seq or identifier in records:
            raise ValueError(f'Empty or duplicate FASTA record in {path}: {identifier}')
        try:seq.encode('ascii')
        except UnicodeEncodeError as exc:raise ValueError(f'Non-ASCII sequence in {path}: {identifier}') from exc
        records[identifier]=seq
    with Path(path).open(encoding='utf-8-sig') as handle:
        for raw in handle:
            line=raw.strip()
            if line.startswith('>'):
                if identifier is not None:save()
                words=line[1:].split();identifier=words[0] if words else '';parts=[]
            elif line:
                if identifier is None:raise ValueError(f'Sequence before first FASTA header: {path}')
                parts.append(line)
        if identifier is not None:save()
    if not records:raise ValueError(f'No FASTA records: {path}')
    return records


def digest(sequence):return hashlib.sha256(sequence.encode('ascii')).hexdigest()


def table(path, delimiter=None):
    with Path(path).open(encoding='utf-8-sig',newline='') as handle:
        reader=csv.DictReader(handle,delimiter=delimiter or ('\t' if str(path).lower().endswith('.tsv') else ','))
        if not reader.fieldnames or any(not f for f in reader.fieldnames) or len(reader.fieldnames)!=len(set(reader.fieldnames)):
            raise ValueError(f'Missing or duplicate table headers: {path}')
        rows=list(reader)
    if any(None in row for row in rows):raise ValueError(f'Malformed row in {path}')
    return rows


def write_tsv(path, rows, fields):
    with Path(path).open('w',encoding='utf-8',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=fields,delimiter='\t',extrasaction='ignore',lineterminator='\n')
        writer.writeheader();writer.writerows(rows)


def prepare(args):
    out=args.out.resolve()
    if out.exists():raise ValueError(f'Output already exists; choose a new blind-run directory: {out}')
    panel=args.panel_proteomes.resolve()
    targets=json.loads(args.target_species.read_text(encoding='utf-8-sig'))
    if isinstance(targets,dict):targets=targets.get('targets',targets.get('target_species'))
    if not isinstance(targets,list) or not targets:raise ValueError('--target-species must contain a nonempty JSON list')
    target_names={str(x).strip() for x in targets}
    if any(not re.fullmatch(r'[A-Za-z0-9_.-]+',name) or name in {'.','..'} for name in target_names):
        raise ValueError('Unsafe target species name in target list')
    controls=table(args.controls)
    required={'Specie','protein_id','original_id','control_class','sequence_sha256'}
    if not controls or not required<=set(controls[0]):raise ValueError('Controls manifest lacks required labelled sequence fields')
    positives=[r for r in controls if r['control_class'].strip().lower().replace(' ','_') in {'positive','positive_control'}]
    if not positives:raise ValueError('No positive controls found in the controls manifest')
    if any(r['Specie'] not in target_names for r in positives):raise ValueError('A positive-control species is absent from target-species JSON')
    if any(not re.fullmatch(r'[0-9a-f]{64}',r['sequence_sha256'] or '') for r in positives):raise ValueError('Invalid positive-control sequence checksum')
    species_files={p.stem:p for p in panel.glob('*.faa') if p.is_file()}
    missing=target_names-species_files.keys()
    if missing:raise ValueError('Target FASTA missing from panel: '+', '.join(sorted(missing)))
    reference_files={name:path for name,path in species_files.items() if name not in target_names}
    if len(reference_files)<2:raise ValueError('Reference-only OrthoFinder inference requires at least two reference proteomes')
    target_sequences={name:fasta(species_files[name]) for name in target_names}
    seen=set();query_rows=[];query_records=[]
    for index,row in enumerate(positives,1):
        key=(row['Specie'],row['protein_id'])
        if not all(key) or key in seen:raise ValueError('Missing or duplicate positive-control target key: '+str(key))
        seen.add(key)
        sequence=target_sequences[key[0]].get(key[1])
        if sequence is None:raise ValueError('Positive control absent from saved target proteome: '+str(key))
        actual=digest(sequence)
        if actual!=row['sequence_sha256']:raise ValueError('Positive control sequence differs from its frozen manifest: '+str(key))
        query_id=f'BlindControl{index:06d}'
        query_records.append((query_id,sequence))
        query_rows.append(dict(query_id=query_id,species=row['species'],target_column=row['Specie'],
            protein_id=row['protein_id'],original_id=row['original_id'],sequence_sha256=actual))
    out.mkdir(parents=True)
    refs=out/'reference_proteomes';refs.mkdir()
    for name,path in sorted(reference_files.items()):shutil.copyfile(path,refs/(name+'.faa'))
    with (out/'withheld_positive_controls.faa').open('w',encoding='utf-8',newline='\n') as handle:
        for query_id,sequence in query_records:
            handle.write('>'+query_id+'\n')
            for i in range(0,len(sequence),60):handle.write(sequence[i:i+60]+'\n')
    write_tsv(out/'withheld_positive_controls.tsv',query_rows,
        ['query_id','species','target_column','protein_id','original_id','sequence_sha256'])
    manifest=dict(version='1.0',status='prepared_reference_only_blind_test',inference_inputs='reference proteomes only',
        withheld_positive_controls_excluded_from_inference=True,matching='exact_full_length_amino_acid_sequence_identity',
        target_species=sorted(target_names),reference_proteomes={name:{'file':name+'.faa','sha256':sha(path)} for name,path in sorted(reference_files.items())},
        controls_manifest_sha256=sha(args.controls),target_species_sha256=sha(args.target_species),
        query_count=len(query_rows),query_fasta_sha256=sha(out/'withheld_positive_controls.faa'))
    (out/'blind_test_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
    print('Prepared reference-only OrthoFinder input: '+str(refs))
    print('Withheld positive controls (do not add to inference input): '+str(out/'withheld_positive_controls.faa'))
    return 0


def evaluate(args):
    root=args.prepared.resolve();groups=args.orthogroups.resolve();out=args.out.resolve()
    if out.exists():raise ValueError(f'Output already exists; choose a new evaluation directory: {out}')
    manifest=json.loads((root/'blind_test_manifest.json').read_text(encoding='utf-8'))
    if manifest.get('status')!='prepared_reference_only_blind_test' or not manifest.get('withheld_positive_controls_excluded_from_inference'):
        raise ValueError('Not a valid reference-only blind-test preparation')
    refs=root/'reference_proteomes';query_fasta=root/'withheld_positive_controls.faa'
    if sha(query_fasta)!=manifest['query_fasta_sha256']:raise ValueError('Withheld-control FASTA changed after preparation')
    ref_sequences={name:fasta(refs/(name+'.faa')) for name in manifest['reference_proteomes']}
    ref_by_hash=defaultdict(list)
    for species,records in ref_sequences.items():
        expected=manifest['reference_proteomes'][species]['sha256']
        if sha(refs/(species+'.faa'))!=expected:raise ValueError(f'Reference FASTA changed: {species}')
        for protein,sequence in records.items():ref_by_hash[digest(sequence)].append((species,protein))
    with groups.open(encoding='utf-8-sig',newline='') as handle:
        reader=csv.DictReader(handle,delimiter='\t');headers=reader.fieldnames or [];rows=list(reader)
    if not headers or headers[0]!='Orthogroup' or len(headers)!=len(set(headers)):raise ValueError('Expected OrthoFinder Orthogroups.tsv')
    expected_columns=set(ref_sequences)
    if set(headers[1:])!=expected_columns:
        raise ValueError('Orthogroup species columns do not exactly match prepared reference FASTAs; expected '+', '.join(sorted(expected_columns)))
    protein_group={};group_members={}
    for row in rows:
        group=row.get('Orthogroup','')
        if not group:raise ValueError('Orthogroups.tsv contains an empty group ID')
        for species in expected_columns:
            for protein in (p.strip() for p in row.get(species,'').split(',') if p.strip()):
                key=(species,protein)
                if key in protein_group:raise ValueError('Reference protein occurs in multiple groups: '+str(key))
                protein_group[key]=group
    queries=fasta(query_fasta);controls=table(root/'withheld_positive_controls.tsv')
    if set(queries)!={r['query_id'] for r in controls}:raise ValueError('Query FASTA and withheld-control table do not match')
    details=[]
    for row in controls:
        sequence=queries[row['query_id']]
        if digest(sequence)!=row['sequence_sha256']:raise ValueError('Withheld control sequence hash mismatch: '+row['query_id'])
        hits=ref_by_hash.get(row['sequence_sha256'],[])
        assigned=sorted({protein_group[(species,protein)] for species,protein in hits if (species,protein) in protein_group})
        if not hits:kind='no_identical_sequence_in_reference_panel'
        elif not assigned:kind='identical_reference_sequence_not_assigned_to_group'
        elif len(assigned)==1:kind='identical_sequence_maps_to_one_orthogroup'
        else:kind='identical_sequence_maps_to_multiple_orthogroups'
        details.append(dict(species=row['species'],target_column=row['target_column'],protein_id=row['protein_id'],
            original_id=row['original_id'],query_id=row['query_id'],match_type=kind,
            identical_reference_sequences=len(hits),matched_reference_proteomes=';'.join(sorted({s for s,_ in hits})),
            matched_reference_ids=';'.join(f'{s}:{p}' for s,p in hits),orthogroups=';'.join(assigned),
            reconstructed=kind=='identical_sequence_maps_to_one_orthogroup'))
    by_species=defaultdict(list)
    for row in details:by_species[row['species']].append(row)
    summary=[]
    for species,items in sorted(by_species.items()):
        counts=Counter(r['match_type'] for r in items)
        assessable=sum(r['identical_reference_sequences']>0 for r in items)
        recovered=counts['identical_sequence_maps_to_one_orthogroup']
        summary.append(dict(species=species,positive_controls=len(items),identity_assessable_controls=assessable,
            identical_sequence_maps_to_one_orthogroup=counts['identical_sequence_maps_to_one_orthogroup'],
            identical_sequence_maps_to_multiple_orthogroups=counts['identical_sequence_maps_to_multiple_orthogroups'],
            identical_reference_sequence_not_assigned_to_group=counts['identical_reference_sequence_not_assigned_to_group'],
            no_identical_sequence_in_reference_panel=counts['no_identical_sequence_in_reference_panel'],
            recovered_fraction_all_controls=recovered/len(items),
            recovered_fraction_identity_assessable=recovered/assessable if assessable else 'not_estimable'))
    if details:
        assessable=sum(r['identical_reference_sequences']>0 for r in details)
        recovered=sum(r['reconstructed'] for r in details)
        summary.append(dict(species='All positive controls',positive_controls=len(details),identity_assessable_controls=assessable,
            identical_sequence_maps_to_one_orthogroup=sum(r['match_type']=='identical_sequence_maps_to_one_orthogroup' for r in details),
            identical_sequence_maps_to_multiple_orthogroups=sum(r['match_type']=='identical_sequence_maps_to_multiple_orthogroups' for r in details),
            identical_reference_sequence_not_assigned_to_group=sum(r['match_type']=='identical_reference_sequence_not_assigned_to_group' for r in details),
            no_identical_sequence_in_reference_panel=sum(r['match_type']=='no_identical_sequence_in_reference_panel' for r in details),
            recovered_fraction_all_controls=recovered/len(details),
            recovered_fraction_identity_assessable=recovered/assessable if assessable else 'not_estimable'))
    out.mkdir(parents=True)
    write_tsv(out/'blind_reconstruction_by_species.tsv',summary,list(summary[0]) if summary else ['species'])
    write_tsv(out/'blind_reconstruction_controls.tsv',details,
        ['species','target_column','protein_id','original_id','query_id','match_type','identical_reference_sequences',
         'matched_reference_proteomes','matched_reference_ids','orthogroups','reconstructed'])
    receipt=dict(version='1.0',status='complete',inference_scope='reference proteomes only',
        control_labels_used_in_inference=False,control_sequences_used_in_inference=False,
        identity_criterion='exact_full_length_amino_acid_sequence_identity',orthogroups_sha256=sha(groups),
        reference_only_manifest_sha256=sha(root/'blind_test_manifest.json'),positive_controls=len(details),
        recovered=sum(r['reconstructed'] for r in details),summary=summary)
    (out/'blind_test_results.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
    print('Blind-test results: '+str(out))
    return 0


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest='command',required=True)
    p=sub.add_parser('prepare',help='Separate reference-only OrthoFinder inputs from withheld positives')
    p.add_argument('--panel-proteomes',type=Path,required=True)
    p.add_argument('--target-species',type=Path,required=True,help='JSON list from the saved panel/targets.json')
    p.add_argument('--controls',type=Path,required=True,help='Branch inputs/sequence_manifest.tsv')
    p.add_argument('--out',type=Path,required=True,help='New, nonexistent blind-test directory')
    p.set_defaults(func=prepare)
    p=sub.add_parser('evaluate',help='Map withheld controls to reference-only OrthoFinder groups by exact sequence identity')
    p.add_argument('--prepared',type=Path,required=True)
    p.add_argument('--orthogroups',type=Path,required=True,help='Reference-only Orthogroups.tsv')
    p.add_argument('--out',type=Path,required=True,help='New, nonexistent evaluation directory')
    p.set_defaults(func=evaluate)
    args=parser.parse_args(argv)
    try:return args.func(args)
    except (OSError,ValueError,KeyError,json.JSONDecodeError) as exc:
        print('ERROR: '+str(exc),file=sys.stderr);return 1


if __name__=='__main__':sys.exit(main())
