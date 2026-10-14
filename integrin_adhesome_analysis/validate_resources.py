"""Validate the current canonical rule, HMM and motif resources without inference."""
import argparse,hashlib,json,re,sys
from pathlib import Path

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'01_detection'))
sys.path.insert(0,str(ROOT/'03_functional_motif_annotation'))
from adhesome_components import load_crosswalk
from adhesome_core import digest
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
    motif=ROOT/'03_functional_motif_annotation/Motif_library'
    core=master_library.load(motif/'adhesion_core.json')
    full=master_library.load(motif/'motif_master.json','master')
    core_rules={r['rule_id']:r['regex_python'] for r in core['rules']}
    full_rules={r['rule_id']:r['regex_python'] for r in full['rules']}
    if any(full_rules.get(k)!=v for k,v in core_rules.items()):raise ValueError('Adhesion/master rule mismatch')
    return dict(pipeline_version='1.0',family_rules=len(rules['families']),pfam_models=len(models['models']),
                components=len(crosswalk),adhesion_motifs=len(core_rules),master_motifs=len(full_rules))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--version',action='version',version='adhesome-resource-validator 1.0')
    parser.parse_args();print(json.dumps(validate(),indent=2))


if __name__=='__main__':main()
