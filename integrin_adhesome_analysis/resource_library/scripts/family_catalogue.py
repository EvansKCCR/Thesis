"""Reconcile biological family sets and programmatic dispatch without truth leakage."""

# Resolve source resources independently of the caller working directory.
import sys as _layout_sys
from pathlib import Path as _LayoutPath
_LAYOUT_BUNDLE_ROOT = next(p for p in _LayoutPath(__file__).resolve().parents
                          if (p/'config/layout_paths.json').is_file() and (p/'pipeline_paths.py').is_file())
if str(_LAYOUT_BUNDLE_ROOT) not in _layout_sys.path:
    _layout_sys.path.insert(0, str(_LAYOUT_BUNDLE_ROOT))
from pipeline_paths import activate as _layout_activate, resource_path as _layout_resource
_layout_activate()

import hashlib
import json
from pathlib import Path
import re

ROOT=_LAYOUT_BUNDLE_ROOT
ANNOTATION=_layout_resource(ROOT/'annotation')
FAMILY_EQUIVALENTS={'src':'Src','Src':'Src','SRC':'Src'}


def canonical_family(name):
    """Explicit user-confirmed aliases only; never generic case folding of biology."""
    return FAMILY_EQUIVALENTS.get(name,name)
METADATA_FIELDS={'family_set','calibration_positive_accessions','calibration_positive_count',
    'calibration_user_inspection_required','calibration_review_flags','calibration_reference_accessions',
    'rule_evidence_origin','rule_status','benchmark_status','review_reason','calibration_test_status'}
GATE_FIELDS={'no_auto_promotion','user_inspection_required','require_ncbi_curated_specific_hit',
    'candidate_rescue_only','required_cdd_repeat_indices'}


def registry_path():
    reconciled=_layout_resource(ANNOTATION/'family_set_registry_reconciled.json')
    return reconciled if _layout_resource(reconciled).is_file() else _layout_resource(ANNOTATION/'family_set_registry.json')


def pair():
    """The supplied canonical pair supersedes the removed *_expanded files."""
    return _layout_resource(ANNOTATION/'family_rules.json'),_layout_resource(ANNOTATION/'family_identity_rules.json')


def resources():
    return [*pair(),registry_path(),_layout_resource(ANNOTATION/'family_set_review_gate.py'),
            _layout_resource(ANNOTATION/'family_set_rule_coverage_summary.tsv')]


def report():
    rules,identities=(json.loads(_layout_resource(p).read_text()) for p in pair())
    reg=json.loads(_layout_resource(registry_path()).read_text());sets=reg['family_sets'];entries=identities['rules']
    stage1=json.loads((_layout_resource(ROOT/'detection/config/family_rules.json')).read_text())['families']
    models=json.loads((_layout_resource(ROOT/'detection/databases/pfam/adhesome_v1.0.hmm.metadata.json')).read_text())['models']
    policies=[]
    for family,rule in sorted(stage1.items()):
        routes=[r for r in entries if family in r['parents'] or r['family']==family or r.get('family_set')==family]
        policies.append(dict(family=family,family_set=rule.get('family_set',family),
            stage2_configured=family in rules['families'],identity_rule_ids=[r['id'] for r in routes],
            handling='screening_hypothesis_only' if rule.get('screening_only') else
            'review_proposal_no_automatic_promotion' if rule.get('no_auto_promotion') or rule.get('user_inspection_required') else
            'diagnostic_identity_rules_and_architecture' if routes else 'architecture_family_only_no_paralog_transfer'))
    return dict(pipeline_version='1.0',authority=str(registry_path().relative_to(ROOT)),
        reconciled_registry_supplied=(_layout_resource(ANNOTATION/'family_set_registry_reconciled.json')).is_file(),
        authority_policy='Prefer explicit reconciled registry; otherwise use supplied registry and report reconciliation gaps',
        family_sets=len(sets),stage1_families=len(stage1),stage2_families=len(rules['families']),
        identity_rules=len(entries),stage1_policies=policies,explicit_family_aliases=FAMILY_EQUIVALENTS,
        stage1_stage2_conflicts=[f for f in stage1.keys()&rules['families'].keys() if stage1[f]!=rules['families'][f]],
        missing_registry_routes=[f for f in sets if not any(r.get('family_set',r['family'])==f for r in entries)],
        duplicate_rule_ids=sorted({r['id'] for r in entries if sum(e['id']==r['id'] for e in entries)>1}),
        proposed_rules=[r['id'] for r in entries if r.get('no_auto_promotion') or r.get('user_inspection_required')],
        coverage_is_not_biological_validation=True,dispatch_is_not_ancestry=True,
        unbundled_stage1_pfam_models=sorted(set(re.findall(r'PF\d{5}',json.dumps(stage1)))-set(models)),
        architectural_alternatives={name:rule['domain_alternatives'] for name,rule in rules['families'].items()
            if len(rule.get('domain_alternatives',[]))>1},
        resource_sha256={str(p.relative_to(ROOT)):hashlib.sha256(_layout_resource(p).read_bytes()).hexdigest() for p in resources()})
