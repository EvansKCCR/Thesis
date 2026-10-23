"""Version 1.0: optional heuristic scaling of orthogroup context, never a probability."""
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
