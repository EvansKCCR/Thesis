"""Empirical correctness calibration; never transfers a control's protein name."""

# Resolve source resources independently of the caller working directory.
import sys as _layout_sys
from pathlib import Path as _LayoutPath
_LAYOUT_BUNDLE_ROOT = next(p for p in _LayoutPath(__file__).resolve().parents
                          if (p/'config/layout_paths.json').is_file() and (p/'pipeline_paths.py').is_file())
if str(_LAYOUT_BUNDLE_ROOT) not in _layout_sys.path:
    _layout_sys.path.insert(0, str(_LAYOUT_BUNDLE_ROOT))
from pipeline_paths import activate as _layout_activate, resource_path as _layout_resource
_layout_activate()


from collections import defaultdict
import hashlib
import json
import math
from pathlib import Path
from pipeline_paths import recorded_resource

ROOT=_LAYOUT_BUNDLE_ROOT
DEFAULT_MODEL=_layout_resource(ROOT/'validation/identity_probability_model.json')


def sha(path):
    with _layout_resource(Path(path)).open('rb') as handle:return hashlib.file_digest(handle,'sha256').hexdigest()


def fingerprint():
    from family_identity import LIBRARY
    files=['scripts/cdd_evidence.py','scripts/family_catalogue.py','scripts/stage2_evidence.py','scripts/identity_probability.py','scripts/family_identity.py','scripts/stage2_resolution.py','scripts/control_exclusion.py','scripts/calibrate_identity.py','scripts/label_ontology.py','detection/config/family_rules.json',
        'annotation/integrated_adhesome_pipeline.py',
        'annotation/audit_policy.json',
        'annotation/family_identity_rules.json']
    return {name:sha(LIBRARY if name.endswith('/family_identity_rules.json') else ROOT/name) for name in files}


def resource_hashes():
    # Negative exclusion gates must invalidate recovery even without a fitted model.
    result={name:sha(ROOT/name) for name in ('scripts/identity_probability.py','scripts/stage2_resolution.py','scripts/control_exclusion.py','scripts/cdd_evidence.py','scripts/family_catalogue.py','scripts/phylogenetic_identity.py')}
    for name in ('calibration_set.xlsx','calibration.faa'):
        path=_layout_resource(ROOT/'validation')/name
        result[path.relative_to(ROOT).as_posix()]=sha(path)
    if _layout_resource(DEFAULT_MODEL).exists():
        result[str(DEFAULT_MODEL.relative_to(ROOT)).replace('\\','/')]=sha(DEFAULT_MODEL)
        model=json.loads(_layout_resource(DEFAULT_MODEL).read_text(encoding='utf-8'))
        for name in model['input_hashes']:
            path=recorded_resource(name);result[name]=sha(path) if _layout_resource(path).exists() else None
    return result


def features(row):
    return dict(score=float(row['total_score_0_10']),
        stratum=str(row.get('family_confidence') or row.get('assignment_confidence') or '')+'|'+str(row.get('protein_name_scope') or ''),
        sequence_sha256=str(row.get('sequence_sha256') or ''),orthology_panel_sha256=str(row.get('orthology_panel_sha256') or ''))


def blocks(samples,min_bin=8):
    """Pool-adjacent-violators with contiguous score bins; no fixed probabilities."""
    grouped=defaultdict(list)
    for sample in samples:grouped[sample['score']].append(sample['correct'])
    bins=[];pending=None
    for score,values in sorted(grouped.items()):
        if pending is None:pending=dict(lo=score,hi=score,n=0,correct=0)
        pending['hi']=score;pending['n']+=len(values);pending['correct']+=sum(values)
        if pending['n']>=min_bin:bins.append(pending);pending=None
    if pending:
        if bins:
            bins[-1]['hi']=pending['hi'];bins[-1]['n']+=pending['n'];bins[-1]['correct']+=pending['correct']
        else:bins.append(pending)
    result=[]
    for item in bins:
        result.append(item)
        while len(result)>1 and result[-2]['correct']/result[-2]['n']>result[-1]['correct']/result[-1]['n']:
            right=result.pop();left=result.pop()
            result.append(dict(lo=left['lo'],hi=right['hi'],n=left['n']+right['n'],correct=left['correct']+right['correct']))
    return result


def estimate(samples,query,minimum=20):
    selected=[s for s in samples if s['stratum']==query['stratum']]
    scope='evidence_stratum'
    if len(selected)<minimum:selected=samples;scope='pooled_evidence'
    if len(selected)<minimum:return dict(status='insufficient_calibration_controls')
    if not min(s['score'] for s in selected)<=query['score']<=max(s['score'] for s in selected):
        return dict(status='outside_observed_calibration_score_range')
    fitted=blocks(selected);chosen=fitted[-1]
    for left,right in zip(fitted,fitted[1:]):
        if query['score']<=(left['hi']+right['lo'])/2:chosen=left;break
    n=chosen['n'];prob=chosen['correct']/n;z=1.959963984540054
    centre=(prob+z*z/(2*n))/(1+z*z/n)
    delta=z*math.sqrt(prob*(1-prob)/n+z*z/(4*n*n))/(1+z*z/n)
    return dict(status='empirical_control_estimate',probability=round(prob,6),lower_95=round(max(0,centre-delta),6),
        upper_95=round(min(1,centre+delta),6),n=n,scope=scope)


def annotate(rows,model_path=DEFAULT_MODEL,family_rules_sha256=None,orthogroup_context_weight=None):
    """A matching training sequence is excluded in full before its probability is fit."""
    if orthogroup_context_weight is not None:
        from orthogroup_context import validate_weight
        validate_weight(orthogroup_context_weight)
        for row in rows:
            for key in ('assignment_probability','assignment_probability_lower_95','assignment_probability_upper_95',
                        'probability_calibration_n','probability_calibration_scope'):row.pop(key,None)
            row['probability_method']='unavailable_for_heuristic_context_configuration'
        return dict(status='probability_model_not_calibrated_for_heuristic_context_configuration',probability_interpretation=False)
    if not _layout_resource(model_path).exists():return dict(status='no_labelled_probability_model')
    model=json.loads(_layout_resource(model_path).read_text(encoding='utf-8'))
    if family_rules_sha256 and model['code_fingerprint'].get('detection/config/family_rules.json')!=family_rules_sha256:
        return dict(status='probability_model_incompatible_with_selected_family_rules')
    if model['code_fingerprint']!=fingerprint():return dict(status='probability_model_stale_after_scoring_change')
    if 'input_hashes' not in model:return dict(status='probability_model_incomplete')
    for name,expected in model['input_hashes'].items():
        path=recorded_resource(name)
        if not _layout_resource(path).exists() or sha(path)!=expected:return dict(status='probability_model_inputs_changed')
    if not {'samples','context','minimum_controls','cross_validation'}<=model.keys():return dict(status='probability_model_incomplete')
    samples=model['samples'];applied=0
    for row in rows:
        if row['assigned_family']=='unassigned':continue
        positive_scope=model['context'].get('positive_family_scope',[])
        if positive_scope and row['assigned_family'] not in positive_scope:
            row['probability_method']='unavailable_outside_positive_control_family_scope';continue
        query=features(row)
        panels=model['context'].get('reference_panels',[])
        if panels and query['orthology_panel_sha256'] not in panels:
            row['probability_method']='unavailable_different_or_missing_reference_panel';continue
        if not query['sequence_sha256']:row['probability_method']='unavailable_missing_sequence_hash';continue
        overlap=any(s['sequence_sha256']==query['sequence_sha256'] for s in samples)
        train=[s for s in samples if s['sequence_sha256']!=query['sequence_sha256']]
        predicted=estimate(train,query,model['minimum_controls'])
        row['probability_method']=('leave_one_sequence_out|' if overlap else 'pooled_control_transfer|')+predicted['status']
        if 'probability' not in predicted:continue
        row.update(assignment_probability=predicted['probability'],assignment_probability_lower_95=predicted['lower_95'],
            assignment_probability_upper_95=predicted['upper_95'],probability_calibration_n=predicted['n'],
            probability_calibration_scope=predicted['scope'])
        applied+=1
    return dict(status='empirical_labelled_control_calibration',applied=applied,model_sha256=sha(model_path),
        context=model['context'],training_controls=len(samples),validation=model['cross_validation'])
