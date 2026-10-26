"""Version 1.0 evidence recruitment and family confirmation, without gold labels."""

# Resolve the organized bundle independently of the caller working directory.
import sys as _bundle_sys
from pathlib import Path as _BundlePath
_BUNDLE_ROOT=next(p for p in _BundlePath(__file__).resolve().parents if (p/'pipeline_paths.py').is_file())
if str(_BUNDLE_ROOT) not in _bundle_sys.path:_bundle_sys.path.insert(0,str(_BUNDLE_ROOT))
from pipeline_paths import activate as _activate_bundle_paths
_activate_bundle_paths()
from contextlib import contextmanager
from functools import lru_cache
import copy
import json
from pathlib import Path
import re

import family_identity as identity

CONFIG=Path(__file__).resolve().parents[1]/'annotation/stage2_evidence_rules.json'


@lru_cache(maxsize=1)
def configuration():return json.loads(CONFIG.read_text())


def default_profile():return configuration()['default_profile']


def profile(rules):
    name=rules.get('_stage2_profile') or configuration()['default_profile']
    if name not in configuration()['profiles']:raise ValueError('Unknown Stage 2 profile: '+name)
    return name,configuration()['profiles'][name]


@contextmanager
def annotation_context(rules):
    name,options=profile(rules)
    if not options['annotation']:
        yield;return
    library=copy.deepcopy(identity.library())
    library.update(import_observed_pfam_accessions=True,accept_cdd_pfam_anchor=True)
    existing={rule['id'] for rule in library['rules']}
    additions=[r for r in configuration()['family_rules'] if r['family'] not in options.get('excluded_new_rule_families',[])]
    if existing & {r['id'] for r in additions}:raise ValueError('Duplicate enhanced identity rule')
    library['rules'].extend(copy.deepcopy(additions))
    with identity.using_configuration(library):yield


def candidate_hypotheses(row,existing,rules,domain_parser):
    """Recruit new hypotheses inside the sealed Stage 1 pool, never new IDs.

    Coordinate-bearing observed domain evidence plus a family-specific CDD
    diagnosis is required. Orthogroups alone, target names and labels cannot route.
    """
    name,options=profile(rules)
    if not options['recruitment']:return []
    hits=identity.combined_domains(row,domain_parser(row));proposals=[]
    for rule in configuration()['family_rules']:
        family=rule['family']
        if family in options.get('excluded_new_rule_families',[]):continue
        if family in existing or not identity.diagnostic_cdd(row,rule['specific_cdd_pattern']):continue
        observed=[len(identity.distinct_footprints(hits,set(group['models']))) for group in rule['required_groups']]
        required=[group['count'] for group in rule['required_groups']]
        recovered=sum(min(a,b) for a,b in zip(observed,required))/sum(required)
        if recovered<rule.get('provisional_min_domain_fraction',1):continue
        core={model for group in rule['required_groups'] for model in group['models']}
        anchors=[h for h in hits if h['model'] in core and h['start'] is not None]
        if not anchors:continue
        proposals.append(dict(row,family=family,assigned_family='',row_origin='stage2_recovery',
            discovery_rule='diagnostic_domain_recruitment:'+rule['id'],recovered_from_families='; '.join(sorted(existing)),
            recruitment_evidence=json.dumps(dict(rule=rule['id'],observed_groups=observed,required_groups=required,
                recovered_fraction=recovered,specific_cdd=True,orthogroup_used_for_recruitment=False),sort_keys=True)))
    return proposals


def partial_confirmation(row,family,result,rules):
    name,options=profile(rules)
    if not options['annotation'] or result['score']==4:return result
    # Partial support names only the tested family; arbitrary shared PH/LIM/FERM
    # fragments must not acquire the name of a different adhesome component.
    partial=[]
    for trace in result.get('rule_trace',[]):
        if trace['family']!=family or not trace['diagnostic_cdd_matched'] or not trace['pfam_anchor_passed']:continue
        if not all(trace[key] for key in ('coordinates_passed','length_passed','maximum_counts_passed')):continue
        groups=trace['domain_groups'];fraction=sum(min(g['observed'],g['required']) for g in groups)/sum(g['required'] for g in groups)
        rule=next((r for r in identity.library()['rules'] if r['id']==trace['rule_id']),None)
        minimum=(rule or {}).get('provisional_min_domain_fraction',1)
        if minimum<=fraction<1:partial.append((fraction,trace['rule_id']))
    if partial:
        fraction,rule_id=max(partial)
        return dict(result,score=2,fraction=fraction,family=family,complete=False,cdd=True,identity_rule=None,
            note='Partial observed '+family+' core and diagnostic CDD; incomplete '+rule_id+'; provisional only')
    if 0<result['score']<4 and not result.get('cdd'):
        return dict(result,score=0,complete=False,identity_rule=None,
            note=result['note']+'; shared partial domains lack family-specific confirmation',
            baseline_gate=dict(configured=True,complete=False,reason_codes=['PROVISIONAL_SPECIFICITY_UNCONFIRMED']))
    return result
