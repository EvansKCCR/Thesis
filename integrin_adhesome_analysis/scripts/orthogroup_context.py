"""Version 1.0: optional heuristic scaling of orthogroup context, never a probability."""

# Resolve the organized bundle independently of the caller working directory.
import sys as _bundle_sys
from pathlib import Path as _BundlePath
_BUNDLE_ROOT=next(p for p in _BundlePath(__file__).resolve().parents if (p/'pipeline_paths.py').is_file())
if str(_BUNDLE_ROOT) not in _bundle_sys.path:_bundle_sys.path.insert(0,str(_BUNDLE_ROOT))
from pipeline_paths import activate as _activate_bundle_paths
_activate_bundle_paths()
import math

HEURISTIC_WEIGHT=0.541


def validate_weight(weight):
    if weight is None:return None
    if isinstance(weight,bool) or not isinstance(weight,(float,int)) or not math.isfinite(weight) or not 0<=weight<=1:
        raise ValueError('orthogroup_context_weight must be a finite number between zero and one')
    return float(weight)


def configuration(weight=None):
    weight=validate_weight(weight)
    return dict(mode='baseline' if weight is None else 'disabled' if weight==0 else 'heuristic',
        orthogroup_context_weight=weight,probability_interpretation=False,
        application='One coefficient applied to the raw orthogroup-context score; replaces, never multiplies, any benchmark coefficient',
        independent_evidence='Sequence, Pfam/CDD, motif, structure, localization and family/exclusion gates are unchanged')
