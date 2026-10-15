"""Validate the current canonical rule, HMM and motif resources without inference."""
import argparse,hashlib,json,re,sys
from pathlib import Path

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'01_detection'))
sys.path.insert(0,str(ROOT/'03_functional_motif_annotation'))
from adhesome_components import load_crosswalk
from adhesome_core import digest, table, sequence_dict
import master_library


def validate():
    detection=ROOT/'01_detection'
    rules=json.loads((detection/'config/family_rules.json').read_text())
    models=json.loads((detection/'databases/pfam/adhesome_v1.0.hmm.metadata.json').read_text())
    library=detection/'databases/pfam/adhesome_v1.0.hmm'
    if digest(library)!=models['library_sha256']:raise ValueError('HMM library checksum mismatch')
    for kind in ('ga','fallback'):
        if digest(str(library)+'.'+kind+'.hmm')!=models['split_sha256'][kind]:raise ValueError('HMM split checksum mismatch')
    required=set(re.findall(r'PF\d{5}',json.dumps(rules['families'])))
    if required-set(models['models']):raise ValueError('Family rules require unbundled Pfam models: '+str(sorted(required-set(models['models']))))
    for name,rule in rules['families'].items():
        if not rule.get('domain_alternatives'):raise ValueError('No architecture for '+name)
        if rule.get('parent_family') and rule['parent_family'] not in rules['families']:raise ValueError('Missing family parent for '+name)
        seen={name};parent=rule.get('parent_family')
        while parent:
            if parent in seen:raise ValueError('Family hierarchy cycle: '+name)
            seen.add(parent);parent=rules['families'][parent].get('parent_family')
    manifest=detection/'config/adhesome_pfam_manifest.csv'
    if digest(manifest)!=models['manifest_sha256']:raise ValueError('Rebuild HMM library after manifest changes')
    crosswalk_path=detection/'config/adhesome_component_crosswalk.csv'
    crosswalk=load_crosswalk(crosswalk_path,manifest)
    if digest(crosswalk_path)!=models['crosswalk_sha256']:raise ValueError('Crosswalk/library metadata mismatch')
    panel=detection/'references/adhesome_v1.0'
    lock=json.loads((panel/'reference_lock.json').read_text())
    for name,expected in lock['files'].items():
        if digest(panel/name)!=expected:raise ValueError('Reference panel checksum mismatch: '+name)
    component_refs=table(panel/'component_catalogue.csv')
    if {row['component_symbol'] for row in component_refs}!={row['Official Symbol'] for row in crosswalk} or len(component_refs)!=len(crosswalk):
        raise ValueError('Reference panel does not cover every crosswalk component exactly once')
    if lock.get('crosswalk_sha256')!=digest(crosswalk_path):raise ValueError('Rebuild component references after crosswalk changes')
    expansion=json.loads((panel/'component_expansion.json').read_text())
    if expansion.get('rules_sha256')!=digest(detection/'config/family_rules.json'):raise ValueError('Rebuild component family labels after rule changes')
    forward=table(panel/'search_catalogue.csv');reverse=table(panel/'metazoan_catalogue.csv')
    for catalogue_rows,fasta_name in ((forward,'search_references.faa'),(reverse,'metazoan_references.faa')):
        sequences=sequence_dict(panel/fasta_name)
        if len(catalogue_rows)!=len(sequences) or {r['reference_id'] for r in catalogue_rows}!=set(sequences):
            raise ValueError('Reference FASTA/catalogue mismatch: '+fasta_name)
        for row in catalogue_rows:
            if hashlib.sha256(sequences[row['reference_id']].encode()).hexdigest()!=row['sequence_sha256']:
                raise ValueError('Reference sequence checksum mismatch: '+row['reference_id'])
    reverse_by_id={r['reference_id']:r for r in reverse}
    for row in forward:
        if row['reference_id'] not in reverse_by_id or any(row[k]!=reverse_by_id[row['reference_id']][k] for k in ('accession','family','sequence_sha256')):
            raise ValueError('Forward reference missing from reverse competition: '+row['reference_id'])
    motif=ROOT/'03_functional_motif_annotation/Motif_library'
    core=master_library.load(motif/'adhesion_core.json')
    full=master_library.load(motif/'motif_master.json','master')
    core_rules={r['rule_id']:r['regex_python'] for r in core['rules']}
    full_rules={r['rule_id']:r['regex_python'] for r in full['rules']}
    if any(full_rules.get(k)!=v for k,v in core_rules.items()):raise ValueError('Adhesion/master rule mismatch')
    return dict(pipeline_version='1.0',family_rules=len(rules['families']),pfam_models=len(models['models']),
                components=len(crosswalk),component_references=len(component_refs),forward_references=len(forward),
                reciprocal_references=len(reverse),adhesion_motifs=len(core_rules),master_motifs=len(full_rules))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--version',action='version',version='adhesome-resource-validator 1.0')
    parser.parse_args();print(json.dumps(validate(),indent=2))


if __name__=='__main__':main()
