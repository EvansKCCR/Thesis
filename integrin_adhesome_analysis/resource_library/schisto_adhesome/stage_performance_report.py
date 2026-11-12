"""Compare labeled performance at workflow stages and post hoc proposals."""
from __future__ import annotations

from collections import Counter
import csv
import html
import json
from pathlib import Path
import hashlib


def read_tsv(path):
    with Path(path).open(encoding='utf-8-sig',newline='') as handle:
        return list(csv.DictReader(handle,delimiter='\t'))


def truthy(value):
    return str(value).strip().casefold() in ('true','1','yes')


def families(value):
    if isinstance(value,(list,tuple,set)):
        return {str(item).strip() for item in value if str(item).strip()}
    return {item.strip() for item in str(value or '').split(';') if item.strip()}


def rate(n,d):
    return n/d if d else None


def score(rows, called, predicted_family):
    counts=Counter();family=Counter();unmapped=0
    for row in rows:
        positive=row.get('control_class')=='positive'
        call=bool(called(row));expected=families(row.get('expected_families'))
        counts['TP' if positive and call else 'FN' if positive else 'FP' if call else 'TN']+=1
        if positive and not expected:
            unmapped+=1
            continue
        label=predicted_family(row) if call else ''
        correct=positive and call and bool(expected & families(label))
        family['TP']+=bool(correct)
        family['FP']+=call and not correct
        family['FN']+=positive and not correct
    tp,fp,tn,fn=(counts[k] for k in ('TP','FP','TN','FN'))
    return dict(n=len(rows),TP=tp,FP=fp,TN=tn,FN=fn,
        precision=rate(tp,tp+fp),recall=rate(tp,tp+fn),F1=rate(2*tp,2*tp+fp+fn),
        family_TP=family['TP'],family_FP=family['FP'],family_FN=family['FN'],
        family_precision=rate(family['TP'],family['TP']+family['FP']),
        family_recall=rate(family['TP'],family['TP']+family['FN']),
        family_F1=rate(2*family['TP'],2*family['TP']+family['FP']+family['FN']),
        positive_family_unmapped=unmapped)


def render(path, kind, metric_rows, change_rows, summary):
    esc=lambda x:html.escape(str(x if x is not None else 'n.a.'))
    def table(rows):
        if not rows:return '<p>No rows.</p>'
        fields=list(rows[0])
        return '<table><thead><tr>'+''.join('<th>'+esc(k)+'</th>' for k in fields)+'</tr></thead><tbody>'+''.join(
            '<tr>'+''.join('<td>'+esc(row.get(k,''))+'</td>' for k in fields)+'</tr>' for row in rows)+'</tbody></table>'
    content=('<!doctype html><meta charset="utf-8"><title>'+esc(kind)+' stage performance</title>'
        '<style>body{font:15px system-ui;margin:30px}table{border-collapse:collapse;margin:18px 0;width:100%}th,td{border:1px solid #ccd;padding:6px;text-align:left}th{background:#eef}</style>'
        '<h1>'+esc(kind.capitalize())+' stage performance</h1><p>'+esc(summary['interpretation'])+'</p>'
        '<p>Every stage uses the same labelled cohort (n='+str(summary['denominator'])+'). Family metrics exclude positive records without mapped expected families.</p>'
        '<h2>Stage metrics</h2>'+table(metric_rows)+'<h2>Change from the prior stage</h2>'+table(change_rows))
    Path(path).write_text(content,encoding='utf-8')


def generate(predictions_path, posthoc_report, refinement, out, kind):
    predictions=read_tsv(predictions_path);posthoc=read_tsv(Path(posthoc_report)/'predictions.tsv')
    key=lambda r:(r.get('Specie',''),r.get('protein_id',''))
    base={key(row):row for row in predictions};proposed={key(row):row for row in posthoc}
    if len(base)!=len(predictions) or len(proposed)!=len(posthoc) or set(base)!=set(proposed):
        raise ValueError('Stage and post hoc reports do not cover the same unique labelled cohort')
    refinement=Path(refinement)
    raw_refinement=read_tsv(refinement/'Stage_2_posthoc_refinement.tsv')
    refinement_manifest=json.loads((refinement/'run_manifest.json').read_text(encoding='utf-8'))
    columns=refinement_manifest.get('column_mapping',{})
    raw_key=lambda row:(str(row.get(columns.get('Specie','Specie'),'')).strip(),
                        str(row.get(columns.get('protein_id','protein_id'),'')).strip())
    raw_refinement_by_key={raw_key(row):row for row in raw_refinement}
    # The all-labelled prediction cohort also contains proteins that Stage 1
    # did not recruit. They correctly remain no-calls and have no row in the
    # audit-derived refinement table. Require unique evidence rows that are a
    # subset of the labelled cohort; the post hoc performance report itself
    # still has to cover the full labelled cohort above.
    if len(raw_refinement_by_key)!=len(raw_refinement) or not set(raw_refinement_by_key).issubset(base):
        raise ValueError('Post hoc evidence table contains duplicate or out-of-cohort records')
    rows=[]
    for protein_key,source in base.items():
        post=proposed[protein_key]
        for field in ('control_class','expected_families'):
            if str(source.get(field,''))!=str(post.get(field,'')):
                raise ValueError('Stage/post hoc label mismatch for '+str(protein_key))
        evidence=raw_refinement_by_key.get(protein_key,{})
        rows.append(dict(source,posthoc_AB=truthy(post.get('posthoc_AB')),
            posthoc_family=post.get('posthoc_evaluation_family',''),
            posthoc_proposed_family=evidence.get('posthoc_proposed_family_or_class',''),
            posthoc_family_correct=truthy(post.get('posthoc_AB_family_correct')),
            original_supported=truthy(post.get('original_supported')),
            posthoc_evidence_rationale=evidence.get('posthoc_evidence_rationale',''),
            posthoc_status=post.get('posthoc_status',''),
            posthoc_tier=post.get('posthoc_tier') or evidence.get('posthoc_evidence_tier',''),
        ))
    stages=[
        ('Stage 1 detection',lambda r:truthy(r.get('stage1_recruited')),
         lambda r:r.get('stage1_family_hypotheses','').split(';')),
        ('Stage 2 annotation',lambda r:truthy(r.get('stage2_annotated')),
         lambda r:r.get('stage2_family','')),
        ('Stage 3 family assignment',lambda r:truthy(r.get('supported')),
         lambda r:r.get('predicted_family','')),
        ('Post hoc A/B proposals only',lambda r:truthy(r.get('posthoc_AB')),
         lambda r:r.get('posthoc_family','')),
        ('Stage 3 plus post hoc A/B proposals',lambda r:truthy(r.get('supported')) or truthy(r.get('posthoc_AB')),
         lambda r:r.get('predicted_family','') if truthy(r.get('supported')) else r.get('posthoc_family',''))]
    metrics=[]
    for name,call,family in stages:
        metrics.append(dict(stage=name,**score(rows,call,family)))
    changes=[]
    for prior,current in zip(metrics,metrics[1:]):
        changes.append(dict(transition=prior['stage']+' → '+current['stage'],
            TP_change=current['TP']-prior['TP'],FP_change=current['FP']-prior['FP'],
            FN_change=current['FN']-prior['FN'],precision_change=(current['precision']-prior['precision']) if current['precision'] is not None and prior['precision'] is not None else None,
            recall_change=(current['recall']-prior['recall']) if current['recall'] is not None and prior['recall'] is not None else None,
            F1_change=(current['F1']-prior['F1']) if current['F1'] is not None and prior['F1'] is not None else None,
            family_FP_change=current['family_FP']-prior['family_FP'],family_FN_change=current['family_FN']-prior['family_FN']))
    out=Path(out);out.mkdir(parents=True,exist_ok=False)
    # Candidate-level changes keep original audit calls intact and identify
    # post hoc assignments as unpromoted proposals.
    change_rows=[]
    for row in rows:
        proposal=bool(str(row.get('posthoc_proposed_family','')).strip())
        if proposal:
            ab_proposal=truthy(row.get('posthoc_AB'))
            rescued=ab_proposal and not truthy(row.get('supported'))
            change_rows.append(dict(Specie=row.get('Specie'),protein_id=row.get('protein_id'),
                control_class=row.get('control_class'),change=('posthoc_AB_rescue_proposal' if rescued else
                    'posthoc_AB_assignment_review' if ab_proposal else 'posthoc_descriptive_C_D_proposal'),
                original_family=row.get('predicted_family'),proposed_family=row.get('posthoc_proposed_family'),
                evaluation_family=row.get('posthoc_family'),
                proposal_tier=row.get('posthoc_tier'),proposal_family_correct=row.get('posthoc_family_correct'),
                evidence=row.get('posthoc_evidence_rationale','')))
    false_positive_path=Path(refinement)/'potential_false_positive_candidates.tsv'
    false_positives=read_tsv(false_positive_path) if false_positive_path.is_file() else []
    for row in false_positives:
        change_rows.append(dict(Specie=row.get('Specie'),protein_id=row.get('protein_id'),
            control_class='audit_assignment_review',change='potential_false_positive_flag',
            original_family=row.get('audited_family'),proposed_family=row.get('domain_proposed_family'),
            proposal_tier=row.get('evidence_tier'),proposal_family_correct='',evidence=row.get('domain_evidence','')))
    summary=dict(version='1.0',kind=kind,status='complete',denominator=len(rows),
        posthoc_candidates_rescued=len([r for r in rows if truthy(r.get('posthoc_AB')) and not truthy(r.get('supported'))]),
        posthoc_AB_proposals=len([r for r in rows if truthy(r.get('posthoc_AB'))]),
        posthoc_all_tier_proposals=sum(bool(str(r.get('posthoc_proposed_family','')).strip()) for r in rows),
        posthoc_rescued_positive_family_correct=sum(r.get('control_class')=='positive' and truthy(r.get('posthoc_AB_family_correct')) and not truthy(r.get('supported')) for r in rows),
        supported_assignments_flagged_as_potential_false_positive=len(false_positives),
        posthoc_assignment_changes_proposed=sum(bool(str(r.get('posthoc_proposed_family','')).strip()) and r.get('posthoc_family')!=r.get('predicted_family') for r in rows),
        posthoc_proposals_are_not_automatically_promoted=True,
        original_audit_preserved=True,
        source_predictions_sha256=hashlib.sha256(Path(predictions_path).read_bytes()).hexdigest(),
        source_posthoc_predictions_sha256=hashlib.sha256((Path(posthoc_report)/'predictions.tsv').read_bytes()).hexdigest(),
        interpretation='Stage 2 is the saved pre-exclusion family annotation; Stage 3 is the final supported A/B audit. Post hoc A/B proposals are shown separately and as a counterfactual cumulative proposal policy; they do not overwrite or promote original calls.')
    with (out/'stage_performance.tsv').open('w',encoding='utf-8',newline='') as h:
        writer=csv.DictWriter(h,fieldnames=list(metrics[0]),delimiter='\t');writer.writeheader();writer.writerows(metrics)
    with (out/'stage_changes.tsv').open('w',encoding='utf-8',newline='') as h:
        writer=csv.DictWriter(h,fieldnames=list(changes[0]) if changes else ['transition','TP_change','FP_change','FN_change','precision_change','recall_change','F1_change','family_FP_change','family_FN_change'],delimiter='\t');writer.writeheader();writer.writerows(changes)
    with (out/'candidate_changes.tsv').open('w',encoding='utf-8',newline='') as h:
        fields=['Specie','protein_id','control_class','change','original_family','proposed_family','evaluation_family','proposal_tier','proposal_family_correct','evidence']
        writer=csv.DictWriter(h,fieldnames=fields,delimiter='\t',extrasaction='ignore');writer.writeheader();writer.writerows(change_rows)
    render(out/'stage_performance.html',kind,metrics,changes,summary)
    summary['output_sha256']={name:hashlib.sha256((out/name).read_bytes()).hexdigest()
        for name in ('stage_performance.tsv','stage_changes.tsv','candidate_changes.tsv','stage_performance.html')}
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
    return summary
