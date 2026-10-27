#!/usr/bin/env python3
"""Calibration-first cached Stage 2 ablation; validation cannot select a profile."""

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
from collections import Counter, defaultdict
import copy
from datetime import datetime, timezone
import gc
import importlib.util
import json
import math
from pathlib import Path
import sys
from unittest.mock import patch

import family_identity
from schisto_adhesome.datasets import ROOT, dump, sha, table
from schisto_adhesome.reports import read, evaluate_branch, summarize_predictions, metrics
from schisto_adhesome.validation_review import family_metrics
import stage2_evidence

ARMS=[('original','original',1.),('recruitment_only','recruitment_only',1.),
      ('annotation_only','annotation_only',1.),('combined','evidence_v1',1.),
      ('combined_without_orthogroup','evidence_v1',0.),('combined_context_0541','evidence_v1',.541),
      ('combined_calibrated','calibrated_v1',1.),('combined_calibrated_without_orthogroup','calibrated_v1',0.),
      ('combined_calibrated_context_0541','calibrated_v1',.541)]


def local(value):
    from evaluate_family_expansion import local_path
    return local_path(str(value))


def load_module():
    path=ROOT/'annotation/integrated_adhesome_pipeline.py'
    spec=importlib.util.spec_from_file_location('stage2_experiment_classifier',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def cohort(run,kind):
    branch=run/kind;audit=branch/'downstream/final'
    path=audit/'joined_evidence.tsv';seal=json.loads((audit/'workflow_manifest.json').read_text())
    if sha(path)!=seal['files']['joined_evidence.tsv']:raise ValueError('Cached evidence changed')
    rows=read(path);manifest=read(branch/'inputs/sequence_manifest.tsv')
    allowed={r['protein_id'] for r in manifest}
    if any(row['protein_id'] not in allowed for row in rows):raise ValueError('Query cohorts are mixed')
    return branch,audit,rows,manifest


def diagnostics(run,kind,out):
    branch,audit,rows,manifest=cohort(run,kind)
    predictions=read(branch/'reports/predictions.tsv')
    grouped=defaultdict(list)
    for row in rows:grouped[row['protein_id']].append(row)
    traces=defaultdict(list)
    missing={r['protein_id']:r for r in predictions if r['control_class']=='positive' and r['supported']=='False'}
    import csv
    with (audit/'stage2_rule_gate_trace.tsv').open(encoding='utf-8-sig') as handle:
        for row in csv.DictReader(handle,delimiter='\t'):
            p=missing.get(row['protein_id'])
            if p and row['family'] in p['expected_families'].split(';'):traces[row['protein_id']].append(row)
    result=[]
    for pid,p in missing.items():
        evidence=grouped.get(pid,[]);expected=traces[pid]
        best=min(expected,key=lambda r:(r['diagnostic_cdd_matched']!='True',len(r['missing_domains']),r['reason_codes'].count(';'))) if expected else {}
        codes=json.loads(best.get('reason_codes') or '[]')
        if not evidence:codes=['NO_STAGE1_RECRUITMENT']
        elif not p['expected_families']:codes=['EXPECTED_FAMILY_UNMAPPED']
        elif not expected:codes=['NO_CONFIGURED_EXPECTED_IDENTITY_RULE']
        elif best.get('diagnostic_cdd_matched')=='True' and 'MISSING_REQUIRED_DOMAIN' in codes:
            codes.append('CHECK_OBSERVED_CDD_PFAM_ACCESSION_IMPORT')
        ev=evidence[0] if evidence else {}
        hypotheses=';'.join(sorted({r['family'] for r in evidence}))
        action=('Review Stage 1 separately; no external-evidence or Stage 2 fabrication' if not evidence else
                'Curate trusted identity mapping separately; do not copy predictions into truth' if not p['expected_families'] else
                'Review diagnostic family model and compatible measured domain core' if not expected else
                'Reconcile exact returned CDD Pfam accessions, routing and core completeness; retain partial cores as provisional')
        result.append(dict(protein_id=pid,original_id=p['original_id'],species=p['species'],expected_family=p['expected_families'],
            stage1_recruited=bool(evidence),stage2_hypotheses=hypotheses,predicted_family=p['predicted_family'],grade=p['grade'],
            Pfam_architecture=ev.get('architecture',''),CDD_architecture=ev.get('NCBI_CDD_architecture',''),
            diagnostic_CDD_models=best.get('matched_cdd_models',''),required_domain_groups=best.get('domain_groups',''),
            orthology_reference_count=ev.get('Orthology_direct_reference_count',''),
            similarity=ev.get('similarity',''),similarity_reason=ev.get('similarity_reason',''),
            motif_context_supported=ev.get('MotifPrediction_context_supported_count',''),
            rejection_reason=';'.join(codes),recommended_action=action,
            unverified_causes='Sequence fragmentation and alignment-based coverage loss cannot be inferred from a failed domain gate alone'))
    out.mkdir(parents=True,exist_ok=True);table(out/'false_negative_audit.tsv',result)
    dump(out/'false_negative_summary.json',dict(kind=kind,missed_positives=len(result),reasons=dict(Counter(code for r in result for code in r['rejection_reason'].split(';'))),
        input_evidence_sha256=sha(audit/'joined_evidence.tsv'),labels_used_for_diagnostics_only=True))
    return result


def run_arm(module,run,kind,out,name,profile,weight,family_pair=None,ablation=None):
    branch,audit,rows,manifest=cohort(run,kind)
    rules_path=audit/'configuration/family_rules.json';identity_path=audit/'configuration/family_identity_rules.json'
    if family_pair:rules_path,identity_path=family_pair
    rules=json.loads(rules_path.read_text());rules['_sha256']=sha(rules_path)
    rules['_audit_policy']=json.loads((ROOT/'annotation/audit_policy.json').read_text())
    rules['_orthogroup_context_weight']=weight;rules['_stage2_profile']=profile
    rules['_identity_ablations']=ablation or {}
    print(f'{kind}: {name} (cached evidence only)',flush=True)
    from contextlib import ExitStack
    architecture=module.architecture
    def diagnostic_only(row,family,rules):
        result=architecture(row,family,rules)
        structural={'MISSING_REQUIRED_DOMAIN','DOMAIN_ORDER_CONFLICT','DOMAIN_COORDINATES_UNCONFIRMED',
            'PFAM_ANCHOR_MISSING','LENGTH_OUTSIDE_RULE_RANGE','EXCESS_DOMAIN_COPIES'}
        eligible={trace['family'] for trace in result.get('rule_trace',[]) if trace.get('reached') and
            trace.get('diagnostic_cdd_matched') and not set(trace['reason_codes'])-structural}
        from family_catalogue import canonical_family
        eligible={canonical_family(f) for f in eligible}
        if len(eligible)==1:
            result=dict(result,score=4,complete=True,family=next(iter(eligible)),cdd=True,
                note='EXPERIMENTAL architecture ablation: reported diagnostic CDD and dispatch retained; full core/order/length/copy gates removed; not a production assignment')
        return result
    with ExitStack() as stack:
        stack.enter_context(family_identity.using_library(identity_path))
        if (ablation or {}).get('architecture') is False:stack.enter_context(patch.object(module,'architecture',side_effect=diagnostic_only))
        audited,resolved=module.classify_rows(rows,rules)
    before={r['protein_id']:r for r in module.collapse_rows(audited)};final={r['protein_id']:r for r in resolved}
    excluded={pid:r for pid,r in before.items() if r.get('stage2_excluded')}
    if set(final)&set(excluded) or set(final)|set(excluded)!=set(before):raise ValueError('Candidate accounting incomplete')
    cal_manifest=read(run/'calibration/inputs/sequence_manifest.tsv')
    references={r['sequence_sha256'] for r in read(branch/'reference_overlap.tsv')}
    for filename in ('search_catalogue.csv','metazoan_catalogue.csv'):
        import csv
        with (ROOT/'detection/references/adhesome_v1.0'/filename).open(encoding='utf-8-sig') as handle:
            references.update(r['sequence_sha256'] for r in csv.DictReader(handle))
    model=json.loads((ROOT/'validation/identity_probability_model.json').read_text())
    def read_arm(unused,code):return tuple({pid:r for pid,r in group.items() if r.get('Specie')==code} for group in (before,final,excluded))
    with patch('evaluate_validation.read_audit',side_effect=read_arm):
        pred=evaluate_branch(branch,manifest,audit,{r['sequence_sha256'] for r in cal_manifest},references,
                             {r['sequence_sha256'] for r in model.get('samples',[])})
    folder=out/kind/name;folder.mkdir(parents=True,exist_ok=False)
    summary=summarize_predictions(pred,folder,kind)
    confusion=Counter((r['expected_families'] or r['control_class'],r['predicted_family'] if r['supported'] else 'abstained') for r in pred)
    table(folder/'family_confusion_matrix.tsv',[dict(expected=a,predicted=b,count=n) for (a,b),n in sorted(confusion.items())])
    table(folder/'identity_outcomes.tsv',[dict(identity_outcome=k,count=v) for k,v in Counter(r.get('identity_outcome','legacy') for r in resolved).items()])
    from audit_report import prepare_sheet
    exported,fields,removed=prepare_sheet('Audited_family_assignment',resolved)
    table(folder/'candidate_assignments.tsv',[{key:row.get(key) for key in fields} for row in exported])
    recruitment=[dict(protein_id=r['protein_id'],family=r['family'],discovery_rule=r.get('discovery_rule'),
        recruitment_evidence=r.get('recruitment_evidence',''),domain_score=r['domain_score_0_4'],grade=r['grade'],
        assigned_family=r['assigned_family']) for r in audited if r.get('row_origin')=='stage2_recovery']
    table(folder/'recruited_hypotheses.tsv',recruitment)
    result=dict(arm=name,profile=profile,context_weight=weight,kind=kind,**summary['supported'],
        family=family_metrics(pred,'supported'),hard_negative_FP=sum(r['control_class']=='hard_negative' and r['supported'] for r in pred),
        background_negative_FP=sum(r['control_class']=='background_negative' and r['supported'] for r in pred),
        recovery_hypotheses=len(recruitment),source_evidence_sha256=sha(audit/'joined_evidence.tsv'))
    dump(folder/'arm_summary.json',result)
    del audited,resolved,before,final,rows;gc.collect()
    return result,pred


def freeze_selection(out,results,predictions):
    baseline=predictions['original'];base={r['protein_id']:r for r in baseline}
    eligible=[]
    for result in results:
        if result['arm']=='original' or result['context_weight']!=1.:continue
        pred=predictions[result['arm']]
        baseline_stats=results[0]
        regressions=[r for r in pred if base[r['protein_id']]['supported'] and base[r['protein_id']]['control_class']=='positive'
                    and base[r['protein_id']]['predicted_family'] in base[r['protein_id']]['expected_families'].split(';')
                    and not (r['supported'] and r['predicted_family'] in r['expected_families'].split(';'))]
        okay=not regressions and result['FP']<=baseline_stats['FP'] and result['hard_negative_FP']<=baseline_stats['hard_negative_FP'] and result['background_negative_FP']<=baseline_stats['background_negative_FP'] and result['family']['FP']<=baseline_stats['family']['FP']
        result['eligible_by_calibration']=okay;result['correct_family_regressions']=len(regressions)
        if okay:eligible.append(result)
    selected=max(eligible,key=lambda r:(r['family']['F1'] or 0,r['F1'] or 0,r['TP'])) if eligible else results[0]
    files=[ROOT/'scripts/stage2_evidence.py',stage2_evidence.CONFIG,ROOT/'scripts/family_identity.py',ROOT/'annotation/integrated_adhesome_pipeline.py']
    frozen=dict(version='1.0',selected_arm=selected['arm'],selected_profile=selected['profile'],context_weight=selected['context_weight'],
        selected_at_utc=datetime.now(timezone.utc).isoformat(),selection_dataset='calibration_only',validation_labels_used=False,
        code_files={str(p):sha(p) for p in files},calibration_results=results,
        interpretation='Context weights are coefficients, not individual orthology correctness probabilities; 0.541 is a historical sensitivity arm.')
    dump(out/'frozen_selection.json',frozen);return frozen


def validation_started(out):
    return (out/'validation_started.json').exists() or any((out/'validation'/name).exists() for name,_,_ in ARMS[1:])


def finalize_freeze(module,run,out):
    """Verify production integration on calibration before sealing validation."""
    if validation_started(out):raise ValueError('Cannot change a freeze after validation starts')
    path=out/'frozen_selection.json';frozen=json.loads(path.read_text());prior=sha(path)
    result,pred=run_arm(module,run,'calibration',out,'integration_verification',frozen['selected_profile'],frozen['context_weight'])
    selected=read(out/'calibration'/frozen['selected_arm']/'predictions.tsv')
    old={r['protein_id']:r for r in selected}
    for row in pred:
        if any(str(row[key])!=old[row['protein_id']][key] for key in ('supported','predicted_family','grade','provisional')):
            raise ValueError('Production integration changed selected calibration calls')
    files=[ROOT/'scripts/stage2_evidence.py',stage2_evidence.CONFIG,ROOT/'scripts/family_identity.py',
        ROOT/'annotation/integrated_adhesome_pipeline.py',
        ROOT/'adhesome_workflow.py',ROOT/'schisto_adhesome/branch.py',ROOT/'schisto_adhesome/workflows.py',
        ROOT/'configs/isolated_workflow.json']
    frozen.update(code_files={str(p):sha(p) for p in files},validation_frozen=True,
        production_integration=dict(finalized_at_utc=datetime.now(timezone.utc).isoformat(),previous_selection_sha256=prior,
            selection_unchanged=True,calibration_predictions_identical=True,verification=result))
    dump(path,frozen)
    print('Frozen for validation: '+frozen['selected_profile'],flush=True)


def phase(args):
    run=args.run.resolve();out=args.out.resolve();out.mkdir(parents=True,exist_ok=True)
    module=load_module()
    if args.phase=='calibration':
        prior=out/'frozen_selection.json'
        if validation_started(out) or prior.exists() and (not args.extend_calibration or json.loads(prior.read_text()).get('validation_frozen')):
            raise ValueError('Calibration selection is already frozen; no extension after validation')
        prior_hash=sha(prior) if prior.exists() else None
        diagnostics(run,'calibration',out/'calibration')
        results=[];predictions={}
        for name,profile,weight in ARMS:
            folder=out/'calibration'/name
            if args.extend_calibration and folder.exists():
                result=json.loads((folder/'arm_summary.json').read_text());pred=read(folder/'predictions.tsv')
                for row in pred:
                    for key in ('supported','stage1_recruited','provisional','including_provisional','before_exclusion'):
                        row[key]=row[key]=='True'
            else:result,pred=run_arm(module,run,'calibration',out,name,profile,weight)
            results.append(result);predictions[name]=pred
            if name=='original':
                saved=read(run/'calibration/reports/predictions.tsv');saved={r['protein_id']:r for r in saved}
                for row in pred:
                    old=saved[row['protein_id']]
                    if any(str(row[key])!=old[key] for key in ('supported','predicted_family','grade')):raise ValueError('Original arm does not reproduce saved baseline')
        frozen=freeze_selection(out,results,predictions)
        if prior_hash:
            frozen['calibration_extension']=dict(previous_selection_sha256=prior_hash,
                reason='Precision-preserving calibration-only variant excludes new LIMK-family rescue after same-family positive/negative conflict',
                validation_evaluated=False)
            dump(out/'frozen_selection.json',frozen)
        print(json.dumps({k:frozen[k] for k in ('selected_arm','selected_profile','context_weight')},indent=2))
    elif args.phase=='freeze':
        finalize_freeze(module,run,out)
    else:
        frozen_path=out/'frozen_selection.json';frozen=json.loads(frozen_path.read_text());before=sha(frozen_path)
        if not frozen.get('validation_frozen'):raise ValueError('Finalize the production freeze before validation')
        if validation_started(out):raise ValueError('Validation already started; do not repeat selection or overwrite results')
        for path,digest in frozen['code_files'].items():
            if sha(Path(path))!=digest:raise ValueError('Frozen experimental configuration changed')
        dump(out/'validation_started.json',dict(started_at_utc=datetime.now(timezone.utc).isoformat(),frozen_selection_sha256=before))
        diagnostics(run,'validation',out/'validation')
        results=[]
        # The measured baseline is reused; independent validation never selects an arm.
        baseline=json.loads((run/'validation/reports/summary.json').read_text())
        dump(out/'validation/original_saved_summary.json',baseline)
        for name,profile,weight in ARMS[1:]:
            result,pred=run_arm(module,run,'validation',out,name,profile,weight);results.append(result)
        dump(out/'validation/ablation_results.json',dict(frozen_selection_sha256=before,selected_arm=frozen['selected_arm'],arms=results,
            validation_used_for_selection=False,baseline=baseline))
        if sha(frozen_path)!=before:raise ValueError('Selection changed during validation')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run',type=Path,required=True);parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--phase',choices=('calibration','freeze','validation'),required=True)
    parser.add_argument('--extend-calibration',action='store_true',help='Add unrun calibration variants only, before any validation evaluation')
    args=parser.parse_args()
    try:phase(args);return 0
    except (OSError,ValueError,KeyError) as exc:print('ERROR: '+str(exc),file=sys.stderr);return 1


if __name__=='__main__':raise SystemExit(main())
