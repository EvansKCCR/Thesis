"""Version 1.0: optional heuristic scaling of orthogroup context, never a probability."""

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
