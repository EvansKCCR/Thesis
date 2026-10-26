#!/usr/bin/env python3
"""Explicit calibration-only exploration; never tune using independent validation."""

# Resolve the organized bundle independently of the caller working directory.
import sys as _bundle_sys
from pathlib import Path as _BundlePath
_BUNDLE_ROOT=next(p for p in _BundlePath(__file__).resolve().parents if (p/'pipeline_paths.py').is_file())
if str(_BUNDLE_ROOT) not in _bundle_sys.path:_bundle_sys.path.insert(0,str(_BUNDLE_ROOT))
from pipeline_paths import activate as _activate_bundle_paths
_activate_bundle_paths()

# Resolve shared pipeline imports from the installed script, not the caller cwd.
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import argparse
import copy
import importlib.util
import json
from pathlib import Path
import sys
from unittest.mock import patch

from schisto_adhesome.datasets import ROOT, dump, sha, table
from schisto_adhesome.workflows import verify_frozen
from schisto_adhesome.reports import read, evaluate_branch, summarize_predictions
from orthogroup_context import validate_weight
import family_identity
import adhesome_workflow as core


def run(args):
    run=args.run.resolve(); contract=json.loads((run/'run_contract.json').read_text())
    verify_frozen(contract)
    branch=run/'calibration'; audit=branch/'downstream/final'
    state=json.loads((branch/'branch_status.json').read_text())
    if state['status']!='complete':raise ValueError('Complete the isolated calibration self-test first')
    core.verify(audit,'audit')
    evidence=audit/'joined_evidence.tsv';rows=core.read(evidence)
    allowed={r['protein_id'] for r in read(branch/'inputs/sequence_manifest.tsv')}
    if any(r['protein_id'] not in allowed for r in rows):
        raise ValueError('Ablation input includes non-calibration query evidence')
    source=core.AUDIT
    spec=importlib.util.spec_from_file_location('isolated_calibration_classifier',source)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    family,identity=core.audit_configuration(contract['frozen_settings']['family_coverage'])
    arms=[('baseline',contract['frozen_settings']['orthogroup_context_weight']),
          ('original_benchmark_policy',None),('heuristic_0541',.541),('ablation_0',0.)]
    for weight in args.additional_weight:
        weight=validate_weight(weight)
        if weight not in [v for k,v in arms]:arms.append(('weight_'+str(weight).replace('.','_'),weight))
    out=args.out.resolve();out.mkdir(parents=True,exist_ok=False)
    manifest=read(branch/'inputs/sequence_manifest.tsv')
    cal_hashes={r['sequence_sha256'] for r in manifest}
    references={r['sequence_sha256'] for r in core.read(ROOT/'detection/references/adhesome_v1.0/search_catalogue.csv',',')}
    model=json.loads((ROOT/'validation/identity_probability_model.json').read_text())
    training={r['sequence_sha256'] for r in model.get('samples',[])}
    outputs=[]; evidence_invariant=None
    from tools.evaluate_orthogroup_ablation import invariants
    for name,weight in arms:
        print('Calibration-only classification: '+name,flush=True)
        rules=json.loads(family.read_text());rules['_sha256']=sha(family)
        rules['_audit_policy']=json.loads(core.AUDIT_POLICY.read_text());rules['_orthogroup_context_weight']=weight
        with family_identity.using_library(identity):
            audited,retained=module.classify_rows(copy.deepcopy(rows),rules)
        fingerprint=invariants(audited)
        if evidence_invariant is None:evidence_invariant=fingerprint
        elif fingerprint!=evidence_invariant:raise ValueError('Independent evidence changed between ablation arms')
        before={r['protein_id']:r for r in module.collapse_rows(audited)}
        final={r['protein_id']:r for r in retained}
        excluded={k:r for k,r in before.items() if r.get('stage2_excluded')}
        if set(final)&set(excluded) or set(final)|set(excluded)!=set(before):
            raise ValueError('Incomplete ablation candidate accounting')
        def read_arm(unused,code):
            return tuple({k:r for k,r in group.items() if r.get('Specie')==code} for group in (before,final,excluded))
        with patch('evaluate_validation.read_audit',side_effect=read_arm):
            predictions=evaluate_branch(branch,manifest,audit,cal_hashes,references,training)
        folder=out/name
        summary=summarize_predictions(predictions,folder,'calibration')
        table(folder/'candidate_assignments.tsv',list(final.values()))
        outputs.append(dict(arm=name,weight=weight,**summary['supported']))
    table(out/'ablation_comparison.tsv',outputs)
    dump(out/'study_manifest.json',dict(version='1.0',dataset='calibration_only',
        parent_contract_sha256=contract['contract_sha256'],evidence_sha256=sha(evidence),
        independent_evidence_sha256=evidence_invariant,arms=outputs,
        validation_labels_used=False,configuration_modified=False,automatic_selection=False,
        recommendation='Keep approved baseline until a justified configuration is explicitly frozen; validate only that frozen configuration.',
        heuristic_interpretation='0.541 is a context coefficient, never individual orthology probability'))
    print('Calibration exploration saved: '+str(out))
    return 0


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run',type=Path,required=True,help='Isolated run containing a completed calibration branch')
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--additional-weight',type=float,action='append',default=[])
    try:return run(p.parse_args(argv))
    except (OSError,ValueError,KeyError) as exc:print('ERROR: '+str(exc),file=sys.stderr);return 1


if __name__=='__main__':sys.exit(main())
