"""CDD hit-category interpretation; a model name is never a specific-hit test."""

# Resolve source resources independently of the caller working directory.
import sys as _layout_sys
from pathlib import Path as _LayoutPath
_LAYOUT_BUNDLE_ROOT = next(p for p in _LayoutPath(__file__).resolve().parents
                          if (p/'config/layout_paths.json').is_file() and (p/'pipeline_paths.py').is_file())
if str(_LAYOUT_BUNDLE_ROOT) not in _layout_sys.path:
    _layout_sys.path.insert(0, str(_LAYOUT_BUNDLE_ROOT))
from pipeline_paths import activate as _layout_activate, resource_path as _layout_resource
_layout_activate()

import math
import re

SOURCE='https://www.ncbi.nlm.nih.gov/Structure/cdd/cdd_help.shtml'


def number(value):
    try:
        value=float(value)
        return value if math.isfinite(value) else None
    except (TypeError,ValueError):return None


def interpretation(hit):
    category=str(hit.get('hit_type') or hit.get('hit_category') or '').strip().lower().replace('_','-')
    source='reported_category'
    # The old parser wrote only specific rows to this explicitly labelled export.
    # This exact marker supports migration; arbitrary coordinate/name strings do not.
    if not category and hit.get('interpretation')=='CDD_specific_domain_prediction_not_family_or_motif_validation':
        category='specific';source='legacy_specific_only_export'
    accession=str(hit.get('accession') or '')
    curated=bool(re.fullmatch(r'cd\d+',accession,re.I))
    evalue=number(hit.get('evalue'));score=number(hit.get('bitscore'))
    threshold=number(hit.get('specific_bitscore_cutoff'))
    cutoff=number(hit.get('specific_evalue_cutoff'))
    reasons=[]
    if category!='specific':reasons.append('CDD_SPECIFIC_CATEGORY_NOT_REPORTED')
    if not accession:reasons.append('CDD_ACCESSION_MISSING')
    if re.fullmatch(r'cl\d+',accession,re.I):reasons.append('CDD_SUPERFAMILY_NOT_FAMILY_DIAGNOSIS')
    if hit.get('evalue') not in (None,'') and evalue is None:reasons.append('CDD_INVALID_EVALUE')
    if hit.get('bitscore') not in (None,'') and score is None:reasons.append('CDD_INVALID_BITSCORE')
    if evalue is not None and evalue<0:reasons.append('CDD_INVALID_EVALUE')
    if threshold is not None and (score is None or score<threshold):reasons.append('CDD_MODEL_SCORE_CUTOFF_FAILED')
    if cutoff is not None and (evalue is None or evalue>cutoff):reasons.append('CDD_MODEL_EVALUE_CUTOFF_FAILED')
    if hit.get('specific_threshold_met') in (False,'false','False','0'):reasons.append('CDD_REPORTED_THRESHOLD_FAILED')
    coverage=number(hit.get('model_coverage'))
    if coverage is not None and not 0<=coverage<=1:reasons.append('CDD_INVALID_MODEL_COVERAGE')
    return dict(category=category or 'unreported',category_source=source,curated_model=curated,
        specific=not reasons,reason_codes=reasons,evalue=evalue,bitscore=score,
        model_coverage=coverage,threshold_source='model_cutoff' if threshold is not None or cutoff is not None else
        'reported_specific_category' if category=='specific' else 'not_available')


def diagnostic(hit,pattern,curated=False):
    state=interpretation(hit)
    return bool(state['specific'] and (not curated or state['curated_model']) and
                re.search(pattern,str(hit.get('short_name') or ''),re.I))
