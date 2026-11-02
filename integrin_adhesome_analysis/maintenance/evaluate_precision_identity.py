#!/usr/bin/env python3
"""Calibration-first cached identity comparisons; never launch searches or inference."""

# Resolve source resources independently of the caller working directory.
import sys as _layout_sys
from pathlib import Path as _LayoutPath
_LAYOUT_BUNDLE_ROOT = next(p for p in _LayoutPath(__file__).resolve().parents
                          if (p/'config/layout_paths.json').is_file() and (p/'pipeline_paths.py').is_file())
if str(_LAYOUT_BUNDLE_ROOT) not in _layout_sys.path:
    _layout_sys.path.insert(0, str(_LAYOUT_BUNDLE_ROOT))
from pipeline_paths import activate as _layout_activate, resource_path as _layout_resource
_layout_activate()

import sys
from pathlib import Path
ROOT=_LAYOUT_BUNDLE_ROOT;sys.path.insert(0,str(ROOT))
from pipeline_paths import activate
activate()
import argparse
from collections import Counter
import json
import csv
import optimize_stage2 as experiment
from family_catalogue import pair,report
from schisto_adhesome.datasets import dump,sha,table
from schisto_adhesome.workflows import algorithm_inventory


def reported_scope_metrics(folder):
    """Correct broad identities are abstentions at finer scope, not wrong paralogs."""
    from calibrate_identity import expected_broad_families
    csv.field_size_limit(100_000_000)
    with (_layout_resource(folder/'candidate_assignments.tsv')).open(encoding='utf-8-sig') as h:
        calls={r['protein_id']:r for r in csv.DictReader(h,delimiter='\t')}
    with (_layout_resource(folder/'predictions.tsv')).open(encoding='utf-8-sig') as h:pred=list(csv.DictReader(h,delimiter='\t'))
    called=correct=unresolved_broad=wrong_specific=0
    for row in pred:
        if row['supported']!='True':continue
        expected=set(row['expected_families'].split(';'))-{''}
        if row['control_class']=='positive' and not expected:continue
        called+=1;call=calls.get(row['protein_id'],{})
        broad=call.get('family_support_level')=='broad_supported'
        accepted=expected|expected_broad_families(row,expected) if broad else expected
        okay=row['control_class']=='positive' and row['predicted_family'] in accepted
        correct+=okay;unresolved_broad+=bool(okay and broad)
        wrong_specific+=bool(not okay and not broad)
    result=dict(calls_with_evaluable_truth=called,correct_at_reported_scope=correct,
        precision_at_reported_scope=correct/called if called else None,
        supported_broad_with_unresolved_specific_identity=unresolved_broad,incorrect_specific_family_calls=wrong_specific,
        interpretation='Specific exact-label metrics remain in family_summary.tsv; broad ancestors are accepted only for explicitly broad-scope predictions. This does not score one-to-one orthology or exact paralog names.')
    dump(_layout_resource(folder/'reported_scope_family_metrics.json'),result);return result


def checked_labels(run,kind):
    from report_audits import checked_cohort,local,read
    from schisto_adhesome.datasets import load_dataset
    contract=json.loads((_layout_resource(run/'run_contract.json')).read_text());reference=contract['datasets'][kind]
    try:
        checked_cohort(run,kind,contract)
        return dict(status='byte_identical_to_frozen_sources')
    except ValueError as original:
        labels=local(reference['labels_file']);proteins=local(reference['protein_file'])
        if sha(proteins)!=reference['fasta_sha256']:raise original
        trusted,summary=load_dataset(kind,proteins,labels,select_labelled=True)
        manifest=read(_layout_resource(run/kind/'inputs/sequence_manifest.tsv'));by_id={r['original_id']:r for r in trusted}
        if len(manifest)!=len(by_id) or {r['original_id'] for r in manifest}!=set(by_id):raise original
        fields=('original_id','species','control_class','accession','genes','protein_name','trusted_family','sequence_sha256','length')
        if any(str(row[key])!=str(by_id[row['original_id']][key]) for row in manifest for key in fields):raise original
        return dict(status='byte_changed_but_all_frozen_evaluation_fields_identical',fields_verified=list(fields),
            records=len(manifest),frozen_labels_sha256=reference['labels_sha256'],current_labels_sha256=sha(labels),
            fasta_sha256=sha(proteins),original_run_modified=False)


def evaluate(run,out,phase,comparison_set='library'):
    out=_layout_resource(Path(out));kind=phase;label_check=checked_labels(_layout_resource(Path(run)),kind)
    if phase=='calibration':
        if _layout_resource(out).exists() and any(p.name!='frozen_identity_configuration.json' for p in out.iterdir()):
            raise ValueError('Comparison outputs already exist; choose a new output directory')
        out.mkdir(parents=True,exist_ok=True)
        dump(_layout_resource(out/'frozen_identity_configuration.json'),dict(code=algorithm_inventory(),family_resources=report(),
            comparison_set=comparison_set,profiles=['calibrated_v1'] if comparison_set=='library' else ['calibrated_v1','precision_v1'],
            selection='Previous calibrated_v1 algorithm maintained; current canonical library compared with the saved library' if comparison_set=='library' else 'Optional strict development comparison; not the main algorithm',
            validation_not_used_for_rule_selection=True))
    else:
        frozen=json.loads((_layout_resource(out/'frozen_identity_configuration.json')).read_text())
        if frozen['code']!=algorithm_inventory():raise ValueError('Code/libraries changed after calibration freeze')
        if frozen.get('comparison_set','development')!=comparison_set:raise ValueError('Comparison set differs from frozen calibration')
    dump(out/(kind+'_label_verification.json'),label_check)
    module=experiment.load_module();results=[];predictions={}
    arms=[('saved_baseline','calibrated_v1',None,None),('expanded_library','calibrated_v1',pair(),None)] if comparison_set=='library' else [
          ('saved_baseline','calibrated_v1',None,None),('canonical_calibrated','calibrated_v1',pair(),None),
          ('precision','precision_v1',pair(),None),('without_redispatch','precision_v1',pair(),{'redispatch':False}),
          ('without_competition_review','precision_v1',pair(),{'competition':False}),
          ('without_architecture_EXPERIMENT_ONLY','precision_v1',pair(),{'architecture':False})]
    for name,profile,resources,ablation in arms:
        result,pred=experiment.run_arm(module,_layout_resource(Path(run)),kind,out,name,profile,1.,resources,ablation)
        results.append(result);predictions[name]={r['protein_id']:r for r in pred}
        reported_scope_metrics(out/kind/name)
    table(_layout_resource(out/kind/'comparison.tsv'),[{k:v for k,v in result.items() if not isinstance(v,(dict,list))} for result in results])
    baseline=predictions['saved_baseline'];changed=[]
    selected_name='expanded_library' if comparison_set=='library' else 'precision'
    for pid,row in predictions[selected_name].items():
        previous=baseline[pid]
        changed.append(dict(protein_id=pid,control_class=row['control_class'],stage1_recruited=row['stage1_recruited'],
            baseline_family=previous['predicted_family'],updated_family=row['predicted_family'],
            baseline_supported=previous['supported'],updated_supported=row['supported'],
            expected_families=row['expected_families'],change='unchanged' if previous['predicted_family']==row['predicted_family'] and previous['supported']==row['supported'] else
            'upgraded' if row['supported'] and not previous['supported'] else 'downgraded_or_abstained' if previous['supported'] and not row['supported'] else 'reassigned'))
    table(_layout_resource(out/kind/'candidate_annotation_changes.tsv'),changed)
    precision=next(r for r in results if r['arm']==selected_name);baseline_result=results[0]
    dump(_layout_resource(out/kind/'precision_review.json'),dict(results=results,
        candidate_change_counts=dict(Counter(r['change'] for r in changed)),
        precision_not_lower=(precision['PPV'] or 0)>=(baseline_result['PPV'] or 0),
        broad_and_specific_annotation_performance='See family_metrics.tsv and family_summary.tsv; these are separate from candidate retention',
        phylogenetic_comparison='Not independent accuracy: no independent adjudicated gene-tree reference supplied',
        clade_holdout='Not estimated: no independently held-out clade labels supplied',
        pr_curve='Not estimated: categorical gate decisions are not a calibrated continuous family/protein score'))
    return results


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--run',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--phase',choices=['calibration','validation'],required=True)
    p.add_argument('--comparison-set',choices=['library','development'],default='library',help='Default: expand the library while preserving calibrated_v1. Strict ablations require development explicitly.')
    a=p.parse_args();print(json.dumps(evaluate(a.run,a.out,a.phase,a.comparison_set),indent=2))


if __name__=='__main__':main()
