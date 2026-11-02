"""Validate the current canonical rule, HMM and motif resources without inference."""

# Resolve source resources independently of the caller working directory.
import sys as _layout_sys
from pathlib import Path as _LayoutPath
_LAYOUT_BUNDLE_ROOT = next(p for p in _LayoutPath(__file__).resolve().parents
                          if (p/'config/layout_paths.json').is_file() and (p/'pipeline_paths.py').is_file())
if str(_LAYOUT_BUNDLE_ROOT) not in _layout_sys.path:
    _layout_sys.path.insert(0, str(_LAYOUT_BUNDLE_ROOT))
from pipeline_paths import activate as _layout_activate, resource_path as _layout_resource
_layout_activate()



# Resolve shared pipeline imports from the installed script, not the caller cwd.
import sys
from pathlib import Path
sys.path.insert(0,str(_layout_resource(Path(__file__)).resolve().parents[1]))
import argparse,hashlib,json,re,sys
from pathlib import Path

ROOT=_LAYOUT_BUNDLE_ROOT
sys.path.insert(0,str(_layout_resource(ROOT/'detection')))
sys.path.insert(0,str(_layout_resource(ROOT/'motif_annotation')))
from adhesome_components import load_crosswalk
from adhesome_core import digest, table, sequence_dict
import master_library
from pipeline_paths import recorded_resource


def validate(allow_stale_calibration=False,audit_resources_only=False):
    detection=_layout_resource(ROOT/'detection')
    rules=json.loads((_layout_resource(detection/'config/family_rules.json')).read_text())
    models=json.loads((_layout_resource(detection/'databases/pfam/adhesome_v1.0.hmm.metadata.json')).read_text())
    library=_layout_resource(detection/'databases/pfam/adhesome_v1.0.hmm')
    if digest(library)!=models['library_sha256']:raise ValueError('HMM library checksum mismatch')
    for kind in ('ga','fallback'):
        if digest(str(library)+'.'+kind+'.hmm')!=models['split_sha256'][kind]:raise ValueError('HMM split checksum mismatch')
    required=set(re.findall(r'PF\d{5}',json.dumps(rules['families'])))
    missing_models=sorted(required-set(models['models']))
    if missing_models and not audit_resources_only:raise ValueError('Family rules require unbundled Pfam models: '+str(missing_models))
    for name,rule in rules['families'].items():
        if not rule.get('domain_alternatives'):raise ValueError('No architecture for '+name)
        if rule.get('parent_family') and rule['parent_family'] not in rules['families']:raise ValueError('Missing family parent for '+name)
        seen={name};parent=rule.get('parent_family')
        while parent:
            if parent in seen:raise ValueError('Family hierarchy cycle: '+name)
            seen.add(parent);parent=rules['families'][parent].get('parent_family')
    manifest=_layout_resource(detection/'config/adhesome_pfam_manifest.csv')
    if digest(manifest)!=models['manifest_sha256']:raise ValueError('Rebuild HMM library after manifest changes')
    crosswalk_path=_layout_resource(detection/'config/adhesome_component_crosswalk.csv')
    crosswalk=load_crosswalk(crosswalk_path,manifest)
    if digest(crosswalk_path)!=models['crosswalk_sha256']:raise ValueError('Crosswalk/library metadata mismatch')
    panel=_layout_resource(detection/'references/adhesome_v1.0')
    lock=json.loads((_layout_resource(panel/'reference_lock.json')).read_text())
    for name,expected in lock['files'].items():
        if digest(panel/name)!=expected:raise ValueError('Reference panel checksum mismatch: '+name)
    component_refs=table(_layout_resource(panel/'component_catalogue.csv'))
    if {row['component_symbol'] for row in component_refs}!={row['Official Symbol'] for row in crosswalk} or len(component_refs)!=len(crosswalk):
        raise ValueError('Reference panel does not cover every crosswalk component exactly once')
    if lock.get('crosswalk_sha256')!=digest(crosswalk_path):raise ValueError('Rebuild component references after crosswalk changes')
    expansion=json.loads((_layout_resource(panel/'component_expansion.json')).read_text())
    reference_labels_current=expansion.get('rules_sha256')==digest(_layout_resource(detection/'config/family_rules.json'))
    if not reference_labels_current and not audit_resources_only:raise ValueError('Rebuild component family labels after rule changes')
    forward=table(_layout_resource(panel/'search_catalogue.csv'));reverse=table(_layout_resource(panel/'metazoan_catalogue.csv'))
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
    motif=_layout_resource(ROOT/'motif_annotation/Motif_library')
    core=master_library.load(_layout_resource(motif/'adhesion_core.json'))
    full=master_library.load(_layout_resource(motif/'motif_master.json'),'master')
    core_rules={r['rule_id']:r['regex_python'] for r in core['rules']}
    full_rules={r['rule_id']:r['regex_python'] for r in full['rules']}
    if any(full_rules.get(k)!=v for k,v in core_rules.items()):raise ValueError('Adhesion/master rule mismatch')
    identity_dir=_layout_resource(ROOT/'annotation')
    identities=json.loads((_layout_resource(identity_dir/'family_identity_rules.json')).read_text())
    if identities['pipeline_version']!='1.0':raise ValueError('Identity library version differs')
    if digest(identity_dir/identities['source_document'])!=identities['source_example_sha256']:
        raise ValueError('Family architecture example checksum mismatch')
    identifiers=[rule['id'] for rule in identities['rules']]
    if len(set(identifiers))!=len(identifiers):raise ValueError('Duplicate family identity rule IDs')
    # A returned, observed CDD model may complement a domain absent from the
    # bundled detection HMM subset; this does not add an HMM or infer a region.
    allowed=set(models['models'])|set(identities['additional_cdd_domains'])|set(identities['pfam_cdd_aliases'])
    allowed|={m for r in identities['rules'] for g in r['required_groups'] for m in g['models'] if re.fullmatch(r'PF\d{5}',m)}
    scopes=set(rules['families'])|set(json.loads((_layout_resource(motif/'context_rules.json')).read_text())['families'])
    from stage2_resolution import BROAD_LABELS
    scopes|=set(BROAD_LABELS)|{rule['family'] for rule in identities['rules']}|{'component_reference_candidate','unassigned_domain_candidate'}
    for rule in identities['rules']:
        if set(rule['parents'])-scopes:raise ValueError('Unknown parent scope: '+rule['id'])
        if not re.fullmatch(r'[A-Za-z0-9_.-]+',rule['family']):raise ValueError('Unsafe identity family name')
        re.compile(rule['specific_cdd_pattern'],re.I)
        for group in rule['required_groups']:
            if group['count']<1 or set(group['models'])-allowed:raise ValueError('Invalid observed identity core: '+rule['id'])
    for pattern in list(identities['pfam_cdd_aliases'].values())+list(identities['additional_cdd_domains'].values()):re.compile(pattern,re.I)
    from control_exclusion import controls
    from identity_probability import DEFAULT_MODEL,fingerprint
    bank=controls();calibration_status='not_fitted';training_controls=0
    if _layout_resource(DEFAULT_MODEL).exists():
        calibration=json.loads(_layout_resource(DEFAULT_MODEL).read_text(encoding='utf-8'))
        if not calibration['context'].get('positive_control_only'):raise ValueError('Refit identity probabilities using positive controls only')
        stale=calibration['code_fingerprint']!=fingerprint()
        if stale and not allow_stale_calibration:raise ValueError('Refit identity probabilities after scoring code/policy changes; --allow-stale-calibration validates resources with numerical probabilities disabled')
        for name,expected in calibration['input_hashes'].items():
            path=recorded_resource(name)
            if not _layout_resource(path).is_file() or digest(path)!=expected:
                if not allow_stale_calibration:raise ValueError('Identity calibration input changed: '+name)
                stale=True
        calibration_status='stale_disabled_pending_local_refit' if stale else 'positive_control_only';training_controls=len(calibration['samples'])
    return dict(pipeline_version='1.0',validation_scope='audit_resources_only' if audit_resources_only else 'full_resources',
                stage1_ready=not missing_models and reference_labels_current,unbundled_stage1_pfam_models=missing_models,
                reference_family_labels_current=reference_labels_current,
                family_rules=len(rules['families']),pfam_models=len(models['models']),
                components=len(crosswalk),component_references=len(component_refs),forward_references=len(forward),
                reciprocal_references=len(reverse),adhesion_motifs=len(core_rules),master_motifs=len(full_rules),
                scoped_identity_rules=len(identifiers),negative_exclusion_sequences=len(bank['negative_by_hash']),
                identity_calibration_status=calibration_status,positive_identity_training_controls=training_controls)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--version',action='version',version='adhesome-resource-validator 1.0')
    parser.add_argument('--allow-stale-calibration',action='store_true',help='Check resources without fitting; report stale numerical probability calibration as disabled')
    parser.add_argument('--audit-resources-only',action='store_true',help='Verify existing audit resources while explicitly reporting Stage 1 model/reference-label gaps')
    args=parser.parse_args();print(json.dumps(validate(args.allow_stale_calibration,args.audit_resources_only),indent=2))


if __name__=='__main__':main()
