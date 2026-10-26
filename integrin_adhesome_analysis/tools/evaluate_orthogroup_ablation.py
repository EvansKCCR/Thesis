#!/usr/bin/env python3
"""Version 1.0: paired classification replay with fixed orthogroup-context weights.

No recruitment, motif search, orthology inference, phylogeny, fitting or tuning.
Gold labels are used only after each arm's classifier has finished.
"""

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
from collections import Counter
import csv
import hashlib
import html
import importlib.util
import json
from pathlib import Path
import platform
import re
import sys
import time

import evaluate_validation as evaluation
from evaluate_family_expansion import local_path,source_context,predictions,validate_pair
import family_identity
import identity_probability
from orthogroup_context import HEURISTIC_WEIGHT,configuration,validate_weight
from validation_reporting import reports,macro_family_metrics,populations,policy_expected_families,policy_prediction,yes

ROOT=Path(__file__).resolve().parents[1]
AUDIT=ROOT/'annotation'
ARMS=[('baseline',None),('benchmark_weighted',HEURISTIC_WEIGHT),('orthogroup_ablation',0.0),
      ('sensitivity_025',.25),('sensitivity_040',.40),('sensitivity_070',.70),('sensitivity_100',1.0)]
POLICIES=('dual_gate_supported','dual_gate_strong','audit_including_provisional')
VIEW_FIELDS=['Specie','protein_id','assigned_family','pre_exclusion_family','retained','excluded','grade','family_confidence',
    'protein_name','protein_name_scope','family_support_level','annotation_resolution','domain_gate_passed',
    'domain_score_0_4','orthology_raw_score_0_3','orthogroup_context_weight','orthogroup_context_policy','orthology_score_0_3',
    'motif_score_0_2','localization_score_0_1','total_score_0_10','equally_ranked_hypotheses',
    'alternative_family_hypotheses','architecture','architecture_evidence','identity_rule','identity_gate_trace']


def digest(value):return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':')).encode()).hexdigest()


def dump(path,value):path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n',encoding='utf-8')


def table(path,rows,fields=None):
    rows=list(rows);fields=fields or list(dict.fromkeys(key for row in rows for key in row)) or ['status']
    cleaned=[{key:json.dumps(value,sort_keys=True) if isinstance(value,(dict,list)) else value for key,value in row.items()} for row in rows]
    evaluation.write_table(path,cleaned,fields)


def sealed_file(path,hashes):
    manifest=path.parent/'workflow_manifest.json'
    data=json.loads(manifest.read_text())
    expected=data['files'].get(path.name)
    if not expected or evaluation.sha(path)!=expected:raise ValueError('Unsealed/changed source: '+str(path))
    for item in (path,manifest):hashes[str(item.resolve())]=evaluation.sha(item)
    return data


def unique(rows):
    result={}
    for row in rows:
        key=(row['Specie'],row['protein_id'])
        if key in result:raise ValueError('Duplicate resolved candidate: '+str(key))
        result[key]=row
    return result


def invariants(audited):
    fields=('audit_family','domain_score_0_4','motif_score_0_2','localization_score_0_1','architecture',
            'NCBI_CDD_architecture','architecture_evidence','reference_subfamily',
            'subfamily_support_method','identity_rule')
    rows=[]
    for row in audited:
        rows.append(dict(Specie=row['Specie'],protein_id=row['protein_id'],origin=row['row_origin'],
            original_row=row.get('original_hypothesis_row'),family=row['preaudit_family'],
            independent_evidence={field:row.get(field) for field in fields}))
    return digest(rows)


def candidate_view(before,final):
    result=[]
    for key,row in sorted(before.items()):
        item={field:row.get(field) for field in VIEW_FIELDS}
        item.update(pre_exclusion_family=row['assigned_family'],retained=key in final,
                    excluded=bool(row.get('stage2_excluded')),assigned_family=final[key]['assigned_family'] if key in final else 'unassigned')
        result.append(item)
    return result


def label_summary(pred):
    result=[];families=[]
    for population,cohort in populations([row for row in pred if row['evaluation_status']=='scored']):
        for policy in POLICIES:
            subset=[row for row in cohort if row.get(policy) is not None]
            tp=fp=fn=unmapped=0
            for row in subset:
                expected=policy_expected_families(row,policy);called=yes(row.get(policy));predicted=policy_prediction(row,policy)
                if row['control_class']=='positive':
                    if not expected:unmapped+=1;continue
                    correct=called and predicted in expected
                    tp+=correct;fn+=not correct;fp+=called and not correct
                else:fp+=called
            macro,per_group=macro_family_metrics(subset,policy)
            result.append(dict(population=population,policy=policy,TP=tp,FP=fp,FN=fn,
                micro_precision=evaluation.rate(tp,tp+fp),micro_recall=evaluation.rate(tp,tp+fn),
                micro_F1=evaluation.rate(2*tp,2*tp+fp+fn),unmapped_positive_controls=unmapped,
                independent_holdout=False,**macro))
            families.extend(dict(population=population,**row) for row in per_group)
    return result,families


def label_errors(pred):
    rows=[]
    for row in pred:
        if row['evaluation_status']!='scored':continue
        for policy in POLICIES:
            if row.get(policy) is None:continue
            expected=policy_expected_families(row,policy);called=yes(row[policy]);predicted=policy_prediction(row,policy)
            if row['control_class']=='positive' and not expected:continue
            wrong=called and (row['control_class']!='positive' or predicted not in expected)
            missed=row['control_class']=='positive' and not (called and predicted in expected)
            if wrong or missed:
                rows.append(dict(policy=policy,accession=row['accession'],protein_id=row['protein_id'],
                    control_class=row['control_class'],expected_families=';'.join(sorted(expected)),predicted_family=predicted,
                    error='FP;FN' if wrong and missed else 'FP' if wrong else 'FN',grade=row['grade'],
                    stage1_recruited=row['stage1_recruitment'],negative_exclusion=row['negative_exclusion']))
    return rows


def run_arm(name,weight,out,module,rows,family_path,identity_path,context,gold_context):
    folder=out/name;folder.mkdir(exist_ok=True)
    rules=json.loads(family_path.read_text());rules['_sha256']=evaluation.sha(family_path)
    rules['_audit_policy']=json.loads((AUDIT/'audit_policy.json').read_text())
    rules['_orthogroup_context_weight']=validate_weight(weight)
    original=digest(rows);started=time.monotonic()
    print('CLASSIFY '+name+' ('+('unchanged baseline' if weight is None else str(weight))+')',flush=True)
    with family_identity.using_library(identity_path):audited,retained=module.classify_rows(rows,rules)
    if digest(rows)!=original:raise ValueError('Classifier changed the immutable joined evidence')
    before=unique(module.collapse_rows(audited));final=unique(retained)
    expected={(row['Specie'],row['protein_id']) for row in rows}
    excluded={key:row for key,row in before.items() if row.get('stage2_excluded')}
    if set(before)!=expected or set(final)&set(excluded) or set(final)|set(excluded)!=expected:
        raise ValueError('Candidate accounting failed in '+name)
    view=candidate_view(before,final)
    if weight is not None and {row['orthogroup_context_weight'] for row in view}!={weight}:
        raise ValueError('Context coefficient is inconsistent across candidates/species')
    table(folder/'candidate_assignments.tsv',view,VIEW_FIELDS)
    contributions=['Specie','protein_id','original_hypothesis_row','row_origin','preaudit_family','audit_family','assigned_family',
        'domain_score_0_4','orthology_raw_score_0_3','orthogroup_context_weight','orthogroup_context_policy','orthology_score_0_3',
        'motif_score_0_2','localization_score_0_1','total_score_0_10','grade','domain_gate_passed','stage2_excluded']
    table(folder/'hypothesis_contributions.tsv',audited,contributions)
    target=next(iter({row['Specie'] for row in context[4]}))
    prior={pid:row for (sp,pid),row in before.items() if sp==target}
    accepted={pid:row for (sp,pid),row in final.items() if sp==target}
    negatives={pid:row for (sp,pid),row in excluded.items() if sp==target}
    # The frozen gold labels are first consulted AFTER classification is complete.
    pred=predictions(*context[:4],*gold_context,prior,accepted,negatives,{})
    metrics,confusions,family_metrics,cells=evaluation.aggregate(pred)
    diagnostics=reports(pred,metrics)
    labels,per_family=label_summary(pred)
    for filename,data in [('control_predictions.tsv',pred),('metrics.tsv',metrics),('confusion_matrices.tsv',confusions),
        ('family_metrics.tsv',family_metrics),('family_confusion_matrix.tsv',cells),('family_performance_by_population.tsv',labels),
        ('per_family_by_population.tsv',per_family),('false_assignments.tsv',label_errors(pred)),*diagnostics.items()]:
        table(folder/filename,data)
    species=[]
    for sp in sorted({key[0] for key in expected}):
        selected=[row for row in view if row['Specie']==sp]
        species.append(dict(Specie=sp,input_candidates=len(selected),retained=sum(row['retained'] for row in selected),
            excluded=sum(row['excluded'] for row in selected),**{'grade_'+g:sum(row['grade']==g for row in selected) for g in 'ABCD'}))
    table(folder/'species_summary.tsv',species)
    summary=dict(arm=name,configuration=configuration(weight),input_hypotheses=len(rows),input_candidates=len(expected),
        retained_candidates=len(final),excluded_candidates=len(excluded),independent_evidence_sha256=invariants(audited),
        input_rows_sha256=original,elapsed_seconds=round(time.monotonic()-started,3),species=species,
        family_performance=[r for r in labels if r['population']=='all_controls'],
        candidate_accounting_complete=True,independent_holdout_available=False,model_refit=False)
    dump(folder/'summary.json',summary)
    dump(folder/'candidate_records.json',view)
    output_hashes={path.name:evaluation.sha(path) for path in folder.iterdir() if path.is_file()}
    dump(folder/'arm_receipt.json',dict(files=output_hashes))
    print('DONE '+name+': '+str(len(final))+' retained; '+str(round(time.monotonic()-started,1))+' s',flush=True)
    return view,summary


def compare(baseline,arm,name,target):
    base=unique(baseline);other=unique(arm)
    if set(base)!=set(other):raise ValueError('Arm changed the candidate universe')
    result=[]
    for key,old in sorted(base.items()):
        new=other[key]
        result.append(dict(arm=name,Specie=key[0],protein_id=key[1],population='validation_control' if key[0]==target else 'target_species',
            baseline_family=old['assigned_family'],arm_family=new['assigned_family'],family_changed=old['assigned_family']!=new['assigned_family'],
            baseline_grade=old['grade'],arm_grade=new['grade'],grade_changed=old['grade']!=new['grade'],
            baseline_confidence=old['family_confidence'],arm_confidence=new['family_confidence'],
            confidence_changed=old['family_confidence']!=new['family_confidence'],
            baseline_retained=old['retained'],arm_retained=new['retained'],
            baseline_weight=old['orthogroup_context_weight'],arm_weight=new['orthogroup_context_weight'],
            baseline_context=old['orthology_score_0_3'],arm_context=new['orthology_score_0_3'],
            baseline_total=old['total_score_0_10'],arm_total=new['total_score_0_10'],
            baseline_provisional=old['grade']=='C',baseline_ambiguous=old['family_confidence']=='ambiguous',
            multi_domain=len(re.findall(r'PF\d{5}',str(old.get('architecture') or '')))>1,
            diagnostic_complete=old['domain_gate_passed'],baseline_name=old['protein_name'],arm_name=new['protein_name'],
            source_reference=old['identity_rule']))
    return result


def changes_summary(changes):
    rows=[]
    for name in sorted({row['arm'] for row in changes}):
        selected=[row for row in changes if row['arm']==name]
        groups=[('all_candidates',selected),('target_species',[r for r in selected if r['population']=='target_species'])]
        groups += [(key,[r for r in selected if r[key]]) for key in ('baseline_provisional','baseline_ambiguous','multi_domain','diagnostic_complete')]
        groups += [(sp,[r for r in selected if r['Specie']==sp]) for sp in sorted({r['Specie'] for r in selected})]
        for population,cohort in groups:
            n=len(cohort);families=sum(r['family_changed'] for r in cohort)
            rows.append(dict(arm=name,population=population,candidates=n,family_changes=families,
                family_change_fraction=evaluation.rate(families,n),grade_changes=sum(r['grade_changed'] for r in cohort),
                confidence_changes=sum(r['confidence_changed'] for r in cohort),
                lost_retention=sum(r['baseline_retained'] and not r['arm_retained'] for r in cohort)))
    return rows


def historical_check(path,view):
    old=unique(evaluation.read_table(path));base=unique(view);changes=[]
    retained={key for key,row in base.items() if row['retained']}
    for key in sorted(set(old)|retained):
        previous=old.get(key,{});current=base.get(key,{})
        if key not in retained or key not in old or previous.get('putative_audited_family')!=current.get('assigned_family') or previous.get('grade')!=current.get('grade'):
            changes.append(dict(Specie=key[0],protein_id=key[1],saved_family=previous.get('putative_audited_family'),
                replay_family=current.get('assigned_family'),saved_grade=previous.get('grade'),replay_grade=current.get('grade')))
    return changes


def family_deltas(rows):
    keyed={(row['arm'],row['population'],row['policy'],row['expected_admissible_families']):row for row in rows}
    result=[]
    for arm in sorted({row['arm'] for row in rows}-{'baseline'}):
        keys={(p,policy,fam) for a,p,policy,fam in keyed if a in ('baseline',arm)}
        for pop,policy,family in sorted(keys):
            base=keyed.get(('baseline',pop,policy,family),{});other=keyed.get((arm,pop,policy,family),{})
            entry=dict(arm=arm,population=pop,policy=policy,expected_admissible_families=family,
                       baseline_group_kind=base.get('group_kind','not_present'),arm_group_kind=other.get('group_kind','not_present'))
            for field in ('TP','FP','FN','precision','recall','F1'):
                old=float(base[field]) if base.get(field) not in ('',None) else None
                new=float(other[field]) if other.get(field) not in ('',None) else None
                entry.update({f'baseline_{field}':old,f'arm_{field}':new,f'delta_{field}':new-old if old is not None and new is not None else None})
            result.append(entry)
    return result


def report(out,summary,performance,changes,group_changes,family_changes):
    def fmt(value):return 'n.a.' if value is None else f'{value:.3f}' if isinstance(value,float) else str(value)
    rows=[r for r in performance if r['population']=='all_controls' and r['policy']=='dual_gate_supported']
    lines=['# Orthogroup-context weighting ablation', '',
        '**Recommendation: retain the expanded baseline. No independent held-out family-label set is available.**', '',
        'The F1-derived coefficient 0.541 is a fixed heuristic evidence-strength adjustment, not an individual orthology correctness probability. The orthology benchmark is not used as validation of downstream family accuracy.', '',
        '## Formulation and controlled comparison', '',
        'Baseline retains the existing benchmark-policy coefficient. Explicit arms use `context = raw_context × w`, replacing any existing benchmark coefficient so it is applied once. The total remains domain + context + motif + localization. All seven target species and the validation target use the same coefficient.', '',
        'The grade rule is nonlinear: complete architecture and positive context gives A; complete architecture without positive context gives B; partial architecture and positive context gives C; other cases give D. Positive fractional coefficients therefore need not change A/B/C grades or family rankings. Zero is a genuine ablation. Architecture/identity, direct sequence, motif, topology and exclusion gates are unchanged. Numeric identity calibration is not refitted or interpreted from the heuristic.', '',
        '## Labelled control self-test', '',
        'The existing labels and reference/calibration-overlap strata are reused. These controls informed family rules and overlap reference/calibration resources. Even a no-exact-overlap stratum is not an independently separated rule-development holdout. Positive controls evaluate family calls; negatives evaluate exclusion gates. Unrecruited positives remain false negatives. Wrong accepted positive-family calls count as both a false positive and a false negative.', '',
        '| Arm | Micro precision | Micro recall | Micro F1 | Macro F1 | FP | FN |',
        '|---|---:|---:|---:|---:|---:|---:|']
    for row in rows:lines.append('| '+' | '.join(fmt(row[key]) for key in ('arm','micro_precision','micro_recall','micro_F1','macro_F1','FP','FN'))+' |')
    lines += ['', '## Candidate effects', '', '| Arm | Target candidates | Family changes | Grade changes |', '|---|---:|---:|---:|']
    for row in group_changes:
        if row['population']=='target_species':lines.append(f"| {row['arm']} | {row['candidates']} | {row['family_changes']} | {row['grade_changes']} |")
    lines += ['', '## Interpretation and limitations', '',
        'Counts and confidence shifts are behavior, not evidence of improved accuracy. The report separates supported A/B calls, strong A calls, and exploratory provisional calls. Sensitivity coefficients were fixed before replay; no best weight was fitted or selected on these controls. Seven-species assignments have no independent ground truth here. Potential benefits or deterioration must be tested on separately curated proteins before changing the default.', '',
        'Family-level changes, false assignments, confidence transitions, ambiguous/multi-domain/provisional subsets, and individual score contributions are exported in the accompanying TSV files and per-arm directories. Missing denominators remain unavailable, not zero. Integrin, talin, kindlin, FAK, ILK and vinculin behavior is included in `biological_family_review.tsv`; unchanged labels alone do not establish their biological correctness.', '',
        '## Reproducibility', '',
        'Frozen joined evidence and the same expanded family/identity pair, policies, labels and exclusions were used in all arms. Only the context coefficient changed. Source hashes and verification are in `study_manifest.json`; every arm records its configuration and per-candidate raw/weighted context. Classification is deterministic, so a random seed is not applicable. No recruitment, search, inference, phylogeny, numeric model training, production-output replacement or backup was performed.', '',
        f"Input candidates: {summary['input_candidates']}. Input hypotheses: {summary['input_hypotheses']}. Independent-evidence invariance: {summary['independent_evidence_unchanged']}. Historical baseline replay: {summary['historical_baseline_reproduced']}."]
    text='\n'.join(lines)+'\n';(out/'report.md').write_text(text,encoding='utf-8')
    # A self-contained, printable scientific summary without an additional workbook.
    def html_table(data,fields):
        return '<table><thead><tr>'+''.join('<th>'+html.escape(k.replace('_',' '))+'</th>' for k in fields)+'</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+html.escape(fmt(row.get(k)))+'</td>' for k in fields)+'</tr>' for row in data)+'</tbody></table>'
    body='<h1>Orthogroup-context weighting ablation</h1><p><strong>Retain the expanded baseline.</strong> No independent held-out family labels are available; results below are exploratory control self-tests.</p>'
    body+='<h2>Family-classification performance</h2>'+html_table(rows,['arm','micro_precision','micro_recall','micro_F1','macro_F1','FP','FN'])
    body+='<h2>Changes across the seven target species</h2>'+html_table([r for r in group_changes if r['population']=='target_species'],['arm','candidates','family_changes','grade_changes','family_change_fraction'])
    body+='<h2>Method and limitations</h2><p>'+html.escape(lines[8])+'</p><p>Positive fractional weights need not change the nonlinear A/B/C gates. Zero removes orthogroup context. Domain, direct sequence, motif, topology and negative-exclusion evidence remain fixed. The coefficient is not a probability. Controls overlap rule-development and reference resources. No independent accuracy improvement is claimed, and reduced counts do not imply higher accuracy.</p><p>Complete methods, denominators, source hashes and candidate-level results accompany this summary.</p>'
    style='body{font:11pt Arial,sans-serif;color:#18232b;max-width:1100px;margin:36px auto;padding:0 20px}h1{font-size:22pt}h2{font-size:14pt;margin-top:28px}p{line-height:1.5}table{border-collapse:collapse;width:100%;font-size:10pt}th{background:#253949;color:white;text-align:left}th,td{padding:8px;border-bottom:1px solid #ccd3d8}td:not(:first-child){text-align:right}@media print{body{margin:12mm;max-width:none}table{page-break-inside:avoid}}'
    (out/'publication_summary.html').write_text('<!doctype html><html lang="en"><meta charset="utf-8"><title>Orthogroup-context ablation</title><style>'+style+'</style><body>'+body+'</body></html>',encoding='utf-8')


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--handoff',type=Path,required=True);parser.add_argument('--evidence',type=Path,required=True)
    parser.add_argument('--benchmark-report',type=Path,required=True);parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--validation-target',default='validation');parser.add_argument('--resume',action='store_true')
    args=parser.parse_args(argv)
    handoff=args.handoff.resolve();evidence=args.evidence.resolve();out=args.out.resolve();benchmark=args.benchmark_report.resolve()
    hashes={};source=sealed_file(evidence,hashes);sealed_file(benchmark,hashes)
    benchmark_data=json.loads(benchmark.read_text());measured=float(benchmark_data['metrics']['f1'])
    if round(measured,3)!=HEURISTIC_WEIGHT:raise ValueError('The fixed 0.541 heuristic does not match this source benchmark; do not silently optimize or substitute a weight')
    context=source_context(handoff,evidence,args.validation_target);hashes.update(context[-1])
    csv.field_size_limit(100_000_000)
    rows=evaluation.read_table(evidence)
    family_path=AUDIT/'family_rules_expanded.json';identity_path=AUDIT/'family_identity_rules_expanded.json'
    validate_pair(json.loads(family_path.read_text()),json.loads(identity_path.read_text()))
    if source.get('rules_sha256')!=evaluation.sha(family_path) or source.get('identity_library_sha256')!=evaluation.sha(identity_path):
        raise ValueError('Frozen baseline uses a different family/identity pair; select the matching expanded baseline')
    for name in ['tools/evaluate_orthogroup_ablation.py','scripts/evaluate_family_expansion.py','scripts/orthogroup_context.py','scripts/family_identity.py',
        'scripts/stage2_resolution.py','scripts/identity_probability.py','scripts/control_exclusion.py','scripts/classification_contract.py','scripts/audit_report.py',
        'scripts/evaluate_validation.py','scripts/validation_reporting.py','scripts/calibrate_identity.py','scripts/label_ontology.py',
        'annotation/integrated_adhesome_pipeline.py',
        'annotation/audit_policy.json']:
        path=ROOT/name;hashes[str(path.resolve())]=evaluation.sha(path)
    for path in (family_path,identity_path,evaluation.LABELS,evaluation.CONTROL_FASTA,identity_probability.DEFAULT_MODEL):
        hashes[str(path.resolve())]=evaluation.sha(path)
    for name in ('search_catalogue.csv','metazoan_catalogue.csv'):
        path=ROOT/'detection/references/adhesome_v1.0'/name;hashes[str(path.resolve())]=evaluation.sha(path)
    historical=evidence.parent/'audited_family_assignments.tsv';sealed_file(historical,hashes)
    gold_context=(*evaluation.load_gold(evaluation.LABELS,evaluation.CONTROL_FASTA,family_path,identity_path),
        {row['sequence_sha256'] for row in json.loads(identity_probability.DEFAULT_MODEL.read_text()).get('samples',[])})
    contract=dict(pipeline_version='1.0',inputs=hashes,evidence=str(evidence),handoff=str(handoff),arms=ARMS,
        validation_target=args.validation_target,random_seed=None,independent_holdout_available=False)
    if out.exists():
        if not args.resume or json.loads((out/'study_contract.json').read_text())!=json.loads(json.dumps(contract)):
            raise ValueError('Choose a new study directory, or resume with identical inputs/configuration')
    else:out.mkdir(parents=True);dump(out/'study_contract.json',contract)
    spec=importlib.util.spec_from_file_location('context_ablation_classifier',AUDIT/'integrated_adhesome_pipeline.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    views={};summaries={};performance=[];family_performance=[]
    for name,weight in ARMS:
        folder=out/name
        if (folder/'arm_receipt.json').is_file():
            receipt=json.loads((folder/'arm_receipt.json').read_text())
            for filename,sha in receipt['files'].items():
                if evaluation.sha(folder/filename)!=sha:raise ValueError('Completed study arm changed: '+str(folder/filename))
            # Candidate dictionaries must retain typed scores/booleans on resume.
            views[name]=json.loads((folder/'candidate_records.json').read_text())
            summaries[name]=json.loads((folder/'summary.json').read_text());print('REUSE '+name,flush=True)
        else:
            if folder.exists() and not args.resume:raise ValueError('Incomplete arm directory; use --resume with identical study inputs: '+str(folder))
            views[name],summaries[name]=run_arm(name,weight,out,module,rows,family_path,identity_path,context,gold_context)
        if name=='baseline':
            differences=historical_check(historical,views[name]);table(out/'historical_baseline_differences.tsv',differences)
            if differences:raise ValueError('Baseline replay differs from saved calls; inspect historical_baseline_differences.tsv before interpreting ablation')
        performance.extend(dict(arm=name,**row) for row in evaluation.read_table(folder/'family_performance_by_population.tsv'))
        family_performance.extend(dict(arm=name,**row) for row in evaluation.read_table(folder/'per_family_by_population.tsv'))
    # Restore numeric metrics after TSV serialization; missing denominators stay None.
    for row in performance:
        for key in ('micro_precision','micro_recall','micro_F1','macro_precision','macro_recall','macro_F1'):
            row[key]=float(row[key]) if row.get(key) not in ('',None) else None
    if len({summary['independent_evidence_sha256'] for summary in summaries.values()})!=1:
        raise ValueError('Non-context evidence changed across study arms')
    base=unique(views['baseline'])
    changes=[row for name in views if name!='baseline' for row in compare(views['baseline'],views[name],name,args.validation_target)]
    group_changes=changes_summary(changes)
    table(out/'candidate_comparisons.tsv',changes);table(out/'classification_change_summary.tsv',group_changes)
    table(out/'family_performance.tsv',performance);table(out/'per_family_performance.tsv',family_performance)
    table(out/'family_performance_changes.tsv',family_deltas(family_performance))
    transitions=Counter((row['arm'],row['Specie'],row['baseline_grade'],row['arm_grade']) for row in changes)
    table(out/'confidence_transitions.tsv',[dict(arm=a,Specie=s,baseline_grade=b,arm_grade=c,count=n) for (a,s,b,c),n in sorted(transitions.items())])
    watched=('integrin','talin','kindlin','FAK','ILK','vinculin','fibronectin')
    table(out/'biological_family_review.tsv',[dict(row,review_reason='named adhesome family; observational consistency only') for row in changes
        if any(name.lower() in (str(row['baseline_family'])+' '+str(row['arm_family'])).lower() for name in watched)])
    summary=dict(input_hypotheses=len(rows),input_candidates=len(base),arms=summaries,independent_evidence_unchanged=True,
        historical_baseline_reproduced=True,independent_validation_available=False,default='baseline',
        recommendation='Retain baseline; no independent labels are available and count/confidence changes alone do not establish accuracy gains',
        benchmark_source_f1=measured,heuristic_weight=HEURISTIC_WEIGHT,probability_interpretation=False,
        model_refit=False,random_seed=None,python=platform.python_version())
    for filename,sha in hashes.items():
        if evaluation.sha(Path(filename))!=sha:raise ValueError('Immutable source changed during study: '+filename)
    dump(out/'study_summary.json',summary)
    report(out,summary,performance,changes,group_changes,family_performance)
    dump(out/'study_manifest.json',dict(**contract,status='complete',source_files_unchanged=True,
        files={str(p.relative_to(out)):evaluation.sha(p) for p in out.rglob('*') if p.is_file()}))
    print('Study complete: '+str(out/'report.md'),flush=True)
    return 0


if __name__=='__main__':sys.exit(main())
