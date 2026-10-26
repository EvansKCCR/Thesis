"""Use the current audit's independent identity gate during detection."""
from functools import lru_cache
import importlib.util
from pathlib import Path

IDENTITY_FIELDS=('reference_subfamily','subfamily_support_method','reciprocal_reference_subfamily',
                 'phylogeny_concordant','subfamily_profile_calibrated','top_profile_subfamily',
                 'profile_margin_pass','reference_accession','reference_dataset_version','profile_domain_gate_verified')


@lru_cache(maxsize=1)
def gate():
    path=Path(__file__).resolve().parents[1]/'annotation/integrated_adhesome_pipeline.py'
    spec=importlib.util.spec_from_file_location('current_adhesome_identity',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module.identity


def identity_evidence(support,rule):
    verified=[r for r in support if gate()(r,rule)[0]==3]
    details={}
    # A conflict remains visible and cannot supply a resolving identity.
    for field in IDENTITY_FIELDS:
        values={str(r[field]) for r in verified if r.get(field) not in (None,'')}
        if len(values)>1:return False,{}
        if values:details[field]=values.pop()
    return bool(verified),details
