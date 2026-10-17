"""Version 1.0: read-only, post-audit evaluation of the labelled validation target."""
import argparse
from collections import Counter,defaultdict
import csv
import hashlib
import html
import importlib.util
import json
from pathlib import Path
import re

from calibrate_identity import expected_families
from classification_contract import general_family

ROOT=Path(__file__).resolve().parent
LABELS=ROOT/'02_calibration_set/calibration_set.xlsx'
CONTROL_FASTA=ROOT/'02_calibration_set/calibration.faa'
POLICIES={
    'stage1_recruitment':'Exported stage 1 hypothesis, including provisional follow-up',
    'audit_before_negative_exclusion':'Stage 1 plus unambiguous A/B call passing the complete family gate, before negative exclusions',
    'dual_gate_supported':'Stage 1 plus unambiguous final A/B call passing the complete family gate, after negative exclusions',
    'dual_gate_strong':'Stage 1 plus a final strong grade A call passing the complete family gate, after negative exclusions',
    'audit_including_provisional':'Stage 1 plus any final A/B/C named family; exploratory coverage, not a passed full family gate'}


def sha(path):
    with Path(path).open('rb') as handle:return hashlib.file_digest(handle,'sha256').hexdigest()


def read_table(path,delimiter='\t'):
    with Path(path).open(encoding='utf-8-sig',newline='') as handle:return list(csv.DictReader(handle,delimiter=delimiter))


def write_table(path,rows,fields):
    with Path(path).open('w',encoding='utf-8',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=fields,delimiter='\t',extrasaction='ignore')
        writer.writeheader();writer.writerows(rows)


def records(path):
    header=None;parts=[]
    with Path(path).open(encoding='utf-8-sig') as handle:
        for line in handle:
            if line.startswith('>'):
                if header is not None:yield header,''.join(parts).upper().removesuffix('*')
                header=line[1:].strip();parts=[]
            elif line.strip():
                if header is None:raise ValueError('Sequence before FASTA header: '+str(path))
                parts.append(''.join(line.split()))
    if header is not None:yield header,''.join(parts).upper().removesuffix('*')


def accession(identifier):
    token=identifier.split()[0];pieces=token.split('|')
    return pieces[1] if len(pieces)>1 else token


def sequence_hash(sequence):return hashlib.sha256(sequence.encode()).hexdigest()
def truthy(value):return str(value).strip().lower() in ('true','1','yes')
def rate(numerator,denominator):return numerator/denominator if denominator else None


def binary_metrics(tp,fp,tn,fn):
    recall=rate(tp,tp+fn);specificity=rate(tn,tn+fp)
    return dict(TP=tp,FP=fp,TN=tn,FN=fn,n=tp+fp+tn+fn,PPV=rate(tp,tp+fp),recall=recall,
        F1=rate(2*tp,2*tp+fp+fn),specificity=specificity,NPV=rate(tn,tn+fn),
        accuracy=rate(tp+tn,tp+fp+tn+fn),
        balanced_accuracy=(recall+specificity)/2 if recall is not None and specificity is not None else None)


def control_class(value):
    key=re.sub('[^a-z]','',str(value).lower())
    return {'positivecontrol':'positive','hardnegativecontrol':'hard_negative',
        'backgroundnegativecontrol':'background_negative'}.get(key)


def resource_hashes():
    paths=[Path(__file__),LABELS,CONTROL_FASTA,ROOT/'calibrate_identity.py',ROOT/'classification_contract.py',
        ROOT/'01_detection/config/family_rules.json',ROOT/'01_detection/references/adhesome_v1.0/search_catalogue.csv',
        ROOT/'01_detection/references/adhesome_v1.0/metazoan_catalogue.csv']
    return {p.relative_to(ROOT).as_posix():sha(p) if p.exists() else None for p in paths}


def load_gold(labels,control_fasta):
    from openpyxl import load_workbook
    book=load_workbook(labels,read_only=True,data_only=True)
    try:
        values=iter(book['control_identity'].values);fields=next(values)
        rows=[dict(zip(fields,row)) for row in values if any(v is not None for v in row)]
    finally:book.close()
    gold={}
    for row in rows:
        key=str(row.get('accession') or '').strip()
        if not key or key in gold:raise ValueError('Missing/duplicate control accession: '+key)
        category=control_class(row.get('class'))
        if not category:raise ValueError('Unknown control class for '+key)
        gold[key]=dict(row,control_class=category)
    trusted={}
    for header,sequence in records(control_fasta):
        key=accession(header)
        if key in trusted:raise ValueError('Duplicate trusted control sequence: '+key)
        trusted[key]=sequence_hash(sequence)
    categories=defaultdict(set)
    for key,row in gold.items():
        if key in trusted:categories[trusted[key]].add(row['control_class']=='positive')
    if any(len(values)>1 for values in categories.values()):raise ValueError('Conflicting positive/negative labels on the same control sequence')
    config=json.loads((ROOT/'01_detection/config/family_rules.json').read_text());catalogue=defaultdict(set);references=set();reverse=set()
    for name,digests in [('search_catalogue.csv',references),('metazoan_catalogue.csv',reverse)]:
        for row in read_table(ROOT/'01_detection/references/adhesome_v1.0'/name,','):
            digests.add(row['sequence_sha256'])
            if row['family'] in config['families']:
                for gene in re.split('[;,\\s]+',row.get('gene','').upper()):catalogue[gene].add(general_family(row['family'],config['families']))
    for row in gold.values():
        explicit=str(row.get('trusted_family') or '').strip()
        row['expected_families']=set(explicit.split(';')) if explicit else expected_families(row,catalogue)
        row['expected_families'].discard('other_named_identity_outside_current_family_catalogue')
    return gold,trusted,references,reverse


def read_audit(audit_dir,target):
    from openpyxl import load_workbook
    path=audit_dir/'Stage_2_output.xlsx'
    if not path.is_file():raise ValueError('A completed version 1.0 Stage_2_output.xlsx is required')
    book=load_workbook(path,read_only=True,data_only=True)
    try:
        sheet=book['Screening_hypothesis'];values=iter(sheet.values);fields=next(values);screen=[]
        for values_row in values:
            row=dict(zip(fields,values_row));pid=str(row.get('protein_id') or '')
            species=row.get('Specie') or pid.split('__',1)[0]
            if species!=target or not pid:continue
            row['assigned_family']=row.get('audit_family') if row.get('grade') in ('A','B','C') else 'unassigned'
            screen.append(row)
    finally:book.close()
    module_path=ROOT/'05_conservative_adhesome_family_assignment_audit/integrated_adhesome_pipeline.py'
    spec=importlib.util.spec_from_file_location('validation_audit_contract',module_path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    before={r['protein_id']:r for r in module.collapse_rows(screen)}
    final_path=audit_dir/'audited_family_assignments.tsv'
    if not final_path.exists():raise ValueError('The compact audited_family_assignments.tsv is required')
    final={}
    for row in read_table(final_path):
        pid=row['protein_id']
        if (row.get('Specie') or pid.split('__',1)[0])!=target:continue
        if pid in final:raise ValueError('Duplicate final audited candidate: '+pid)
        final[pid]=row
    exclusions={}
    for row in read_table(audit_dir/'negative_control_exclusions.tsv'):
        if row.get('Specie')==target:
            if row['protein_id'] in exclusions:raise ValueError('Duplicate excluded validation candidate')
            exclusions[row['protein_id']]=row
    if set(final)&set(exclusions):raise ValueError('An excluded negative control remains in final assignments')
    if set(final)|set(exclusions)!=set(before):raise ValueError('Audited validation protein accounting is incomplete')
    return before,final,exclusions


def supported(row,strong=False):
    if not row:return False
    family=row.get('putative_audited_family') or row.get('assigned_family') or 'unassigned'
    confidence=row.get('assignment_confidence') or row.get('family_confidence') or ''
    complete=float(row.get('domain_score_0_4') or 0)==4
    return family not in ('unassigned','unresolved','') and complete and confidence!='ambiguous' and (
        row.get('grade')=='A' and confidence=='strong' if strong else row.get('grade') in ('A','B'))


def aggregate(predictions):
    scored=[row for row in predictions if row['evaluation_status']=='scored']
    metrics=[];confusions=[];family_metrics=[];family_cells=Counter()
    for policy in POLICIES:
        for population,allowed in [('all_controls',{'positive','hard_negative','background_negative'}),
                ('positive_vs_hard_negative',{'positive','hard_negative'}),('positive_vs_background_negative',{'positive','background_negative'}),
                ('sequence_verified_controls',{'positive','hard_negative','background_negative'})]:
            cohort=[r for r in scored if r['control_class'] in allowed and
                (population!='sequence_verified_controls' or r['label_sequence_verification']=='trusted_fasta_match')];counts=Counter()
            for row in cohort:
                actual=row['control_class']=='positive';called=bool(row[policy])
                counts['TP' if actual and called else 'FN' if actual else 'FP' if called else 'TN']+=1
            result=dict(policy=policy,population=population,**binary_metrics(*(counts[k] for k in ('TP','FP','TN','FN'))))
            metrics.append(result)
            for key,actual,predicted in [('TP','positive','accepted'),('FN','positive','not_accepted'),('FP','negative','accepted'),('TN','negative','not_accepted')]:
                confusions.append(dict(policy=policy,population=population,actual=actual,predicted=predicted,count=counts[key]))
        if policy=='stage1_recruitment':continue
        tp=fp=fn=unmapped=0
        for row in scored:
            actual_positive=row['control_class']=='positive';expected=set(row['expected_families'].split(';'))-{''}
            called=bool(row[policy]);predicted=row['pre_exclusion_family'] if policy=='audit_before_negative_exclusion' else row['audited_family']
            if actual_positive and not expected:unmapped+=1;continue
            if actual_positive:
                correct=called and predicted in expected
                tp+=correct;fn+=not correct;fp+=called and not correct
                actual=';'.join(sorted(expected))
            else:fp+=called;actual='__'+row['control_class']+'__'
            family_cells[policy,actual,predicted if called else '__not_accepted__']+=1
        family_metrics.append(dict(policy=policy,TP=tp,FP=fp,FN=fn,PPV=rate(tp,tp+fp),recall=rate(tp,tp+fn),
            F1=rate(2*tp,2*tp+fp+fn),positive_controls_without_mapped_family=unmapped))
    families=[dict(policy=p,actual_admissible_families=a,predicted_family=b,count=n) for (p,a,b),n in sorted(family_cells.items())]
    return metrics,confusions,family_metrics,families


def write_html(path,summary,metrics,confusions,family_metrics):
    esc=lambda value:html.escape(str(value))
    def table(rows,fields):
        return '<table><thead><tr>'+''.join('<th>'+esc(f)+'</th>' for f in fields)+'</tr></thead><tbody>'+''.join(
            '<tr>'+''.join('<td>'+('n.a.' if row.get(f) is None else f'{row[f]:.4f}' if isinstance(row.get(f),float) else esc(row.get(f,'')))+'</td>' for f in fields)+'</tr>' for row in rows)+'</tbody></table>'
    matrices=[]
    for policy in POLICIES:
        cells={(r['actual'],r['predicted']):r['count'] for r in confusions if r['policy']==policy and r['population']=='all_controls'}
        rows=[dict(actual=label,accepted=cells[label,'accepted'],not_accepted=cells[label,'not_accepted']) for label in ('positive','negative')]
        matrices.append('<h3>'+esc(policy)+'</h3><p>'+esc(POLICIES[policy])+'</p>'+table(rows,['actual','accepted','not_accepted']))
    body='<h1>Validation target performance</h1><p>'+esc(summary['interpretation'])+'</p>'
    body+='<p>Scored controls: '+str(summary['scored_controls'])+' / input sequences: '+str(summary['input_sequences'])+'. Unscored reasons: '+esc(summary['unscored_reasons'])+'.</p>'
    body+='<p>Sequence-verified labels: '+str(summary['sequence_verified_controls'])+'. Accession labels without a duplicate trusted sequence: '+str(summary['accession_label_only_controls'])+'. The sequence-verified comparison is included in metrics.tsv.</p>'
    body+='<h2>Membership performance</h2><p>PPV is precision. F1 = 2TP / (2TP + FP + FN). Rates are fractions from 0 to 1.</p>'+table([r for r in metrics if r['population']=='all_controls'],['policy','TP','FP','TN','FN','PPV','recall','F1','specificity','accuracy'])
    body+='<h2>Confusion matrices</h2><p>Rows are actual classes; columns are predicted decisions. A correct membership decision does not prove a correct family name.</p>'+''.join(matrices)
    body+='<h2>Family correctness at the reported scope</h2>'+table(family_metrics,['policy','TP','FP','FN','PPV','recall','F1','positive_controls_without_mapped_family'])
    body+='<p>A wrong positive-family call contributes one FP and one FN here. Admissible family labels come only from the trusted labels and ontology. This does not evaluate exact gene, isoform or paralogue names.</p>'
    body+='<h2>Interpretation and files</h2><p>Reference overlap: '+esc(summary['exact_forward_reference_overlap'])+' scored controls. Exact control exclusions and calibration overlap are reported per protein. Negative controls remain exclusion gates; post-exclusion specificity includes these configured rules.</p>'
    body+='<p>Undefined denominators are n.a. Reports never change assignments, confidence or evidence weights. See metrics.tsv, confusion_matrices.tsv, family_metrics.tsv, family_confusion_matrix.tsv and control_predictions.tsv for full data.</p>'
    style='body{font:15px system-ui;color:#15324a;margin:36px;max-width:1400px}table{border-collapse:collapse;margin:18px 0 32px;width:100%}th,td{border:1px solid #cedae4;padding:9px;text-align:right;overflow-wrap:anywhere}th:first-child,td:first-child{text-align:left}th{background:#e8f0f6}h2{margin-top:38px}p{line-height:1.5}'
    path.write_text('<!doctype html><html><head><meta charset="utf-8"><title>Validation performance</title><style>'+style+'</style></head><body>'+body+'</body></html>',encoding='utf-8')


def evaluate(audit_dir,handoff,out,target=None,labels=LABELS,control_fasta=CONTROL_FASTA):
    audit_dir=Path(audit_dir).resolve();handoff=Path(handoff).resolve();out=Path(out).resolve()
    if out.exists():raise FileExistsError('Choose a new validation report directory: '+str(out))
    hand=json.loads((handoff/'workflow_manifest.json').read_text(encoding='utf-8'))
    targets=set(hand.get('target_species',[]))
    if not hand.get('detection'):
        if target is None:return dict(status='validation_target_not_present')
        raise ValueError('The handoff has no linked detection batch')
    detection=Path(hand['detection'])
    if target is None:
        matches=[]
        for code in sorted(targets):
            config_path=detection/code/'proteome.json'
            if not config_path.exists():continue
            settings=json.loads(config_path.read_text(encoding='utf-8'));entry=settings['proteomes'][0]
            if code.lower()=='validation' or Path(entry['path']).stem.lower()=='validation' or entry.get('species','').lower()=='pipeline validation':matches.append(code)
        if not matches:return dict(status='validation_target_not_present')
        if len(matches)>1:raise ValueError('Multiple validation targets; specify --target')
        target=matches[0]
    if target not in targets:raise ValueError('Validation target is not in the audited handoff: '+target)
    config_path=detection/target/'proteome.json';settings=json.loads(config_path.read_text(encoding='utf-8'))
    raw=Path(settings['proteomes'][0]['path'])
    if not raw.is_absolute():raw=(config_path.parent/raw).resolve()
    prepared_dir=detection/target/'prepared';preparation=prepared_dir/'preparation.json'
    receipt=json.loads(preparation.read_text(encoding='utf-8'))['species'][target]
    prepared=prepared_dir/(target+'.faa');mapping_path=prepared_dir/(target+'.id_map.csv')
    if sha(raw)!=receipt['input_sha256'] or sha(prepared)!=receipt['output_sha256']:raise ValueError('Validation input/prepared sequences changed after detection')
    for path in (prepared,detection/target/'candidates/candidates.csv'):
        expected=hand.get('sources',{}).get(str(path.resolve()))
        if expected is None or sha(path)!=expected:raise ValueError('Unverified validation detection source: '+str(path))
    gold,trusted,references,reverse=load_gold(Path(labels),Path(control_fasta))
    mapping=read_table(mapping_path,',');by_original={r['original_id']:r for r in mapping}
    if len(by_original)!=len(mapping):raise ValueError('Duplicate validation original-ID mappings')
    sequences={h.split()[0]:s for h,s in records(prepared)}
    selected={r['protein_id'] for r in read_table(handoff/'hypotheses.tsv') if (r.get('Specie') or r['protein_id'].split('__',1)[0])==target}
    if not selected<=set(sequences):raise ValueError('A recruited validation candidate is absent from its prepared proteome')
    before,final,excluded=read_audit(audit_dir,target)
    if set(before)!=selected:raise ValueError('Validation audit does not account for every exported stage 1 candidate')
    model=ROOT/'02_calibration_set/identity_probability_model.json'
    training={r['sequence_sha256'] for r in json.loads(model.read_text()).get('samples',[])} if model.exists() else set()
    predictions=[];seen=set();input_ids=set()
    for header,sequence in records(raw):
        original=header.split()[0];key=accession(original);digest=sequence_hash(sequence)
        if key in seen:raise ValueError('Duplicate validation control accession: '+key)
        seen.add(key);input_ids.add(original)
        if original not in by_original:raise ValueError('Validation ID missing from preparation map: '+original)
        unit=by_original[original];pid=unit['protein_id'];retained=unit['decision']=='retained'
        if retained and sequences.get(pid)!=sequence:raise ValueError('Prepared validation sequence differs from original control: '+pid)
        actual=gold.get(key,{})
        status='unlabelled_control' if not actual else 'labelled_sequence_mismatch' if key in trusted and trusted[key]!=digest else 'scored'
        verification='not_labelled' if not actual else 'accession_label_only' if key not in trusted else 'trusted_fasta_match' if trusted[key]==digest else 'trusted_sequence_mismatch'
        current=final.get(pid,{}) if retained else {};prior=before.get(pid,{}) if retained else {}
        gate1=retained and pid in selected
        final_family=current.get('putative_audited_family') or current.get('assigned_family') or 'unassigned'
        expected=actual.get('expected_families',set())
        row=dict(accession=key,protein_id=pid,original_id=original,sequence_sha256=digest,evaluation_status=status,
            label_sequence_verification=verification,
            control_class=actual.get('control_class',''),trusted_protein_name=actual.get('protein_name',''),
            expected_families=';'.join(sorted(expected)),preparation_decision=unit['decision'],
            stage1_recruitment=gate1,audit_before_negative_exclusion=gate1 and supported(prior),
            dual_gate_supported=gate1 and supported(current),dual_gate_strong=gate1 and supported(current,strong=True),
            audit_including_provisional=gate1 and current.get('grade') in ('A','B','C') and final_family not in ('','unassigned','unresolved'),
            audited_family=final_family,pre_exclusion_family=prior.get('assigned_family','unassigned'),grade=current.get('grade',''),
            assignment_confidence=current.get('assignment_confidence',''),protein_name=current.get('protein_name',''),
            protein_name_scope=current.get('protein_name_scope',''),negative_exclusion=retained and pid in excluded,
            exclusion_gate=excluded.get(pid,{}).get('exclusion_gate','') if retained else '',
            exact_forward_reference_overlap=digest in references,exact_reverse_reference_overlap=digest in reverse,
            identity_calibration_training_overlap=digest in training)
        predictions.append(row)
    if input_ids!=set(by_original):raise ValueError('Validation input/preparation-map membership differs')
    scored=[r for r in predictions if r['evaluation_status']=='scored']
    if not scored:raise ValueError('No sequence-verified labelled controls in the validation target')
    metrics,confusions,family_metrics,families=aggregate(predictions)
    main=next(r for r in metrics if r['policy']=='dual_gate_supported' and r['population']=='all_controls')
    summary=dict(pipeline_version='1.0',status='complete',validation_target=target,input_sequences=len(predictions),
        labelled_input_controls=sum(r['accession'] in gold for r in predictions),scored_controls=len(scored),
        unique_scored_sequences=len({r['sequence_sha256'] for r in scored}),
        sequence_verified_controls=sum(r['label_sequence_verification']=='trusted_fasta_match' for r in scored),
        accession_label_only_controls=sum(r['label_sequence_verification']=='accession_label_only' for r in scored),
        unscored_reasons=dict(Counter(r['evaluation_status'] for r in predictions if r['evaluation_status']!='scored')),
        scored_control_classes=dict(Counter(r['control_class'] for r in scored)),
        labels_absent_from_validation_input=sorted(set(gold)-seen),exact_forward_reference_overlap=sum(r['exact_forward_reference_overlap'] for r in scored),
        identity_calibration_training_overlap=sum(r['identity_calibration_training_overlap'] for r in scored),
        negative_control_exclusions=sum(r['negative_exclusion'] for r in scored),primary_policy='dual_gate_supported',primary_metrics=main,
        policies=POLICIES,undefined_metric_policy='null in JSON, blank in TSV, n.a. in HTML; never fabricated zero',
        interpretation='Reference-panel control self-test with exact reference overlap, calibrated controls and configured negative exclusions. This is not independent held-out or novel-lineage performance. Membership and family correctness are reported separately; evaluation never changes calls or weights.')
    paths=[raw,prepared,mapping_path,preparation,config_path,Path(labels),Path(control_fasta),handoff/'workflow_manifest.json',handoff/'hypotheses.tsv',
        detection/target/'candidates/candidates.csv',audit_dir/'Stage_2_output.xlsx',audit_dir/'audited_family_assignments.tsv',audit_dir/'negative_control_exclusions.tsv']
    summary['input_sha256']={str(p):sha(p) for p in paths};summary['evaluation_resources']=resource_hashes()
    out.mkdir(parents=True,exist_ok=False)
    for name,rows in [('metrics.tsv',metrics),('confusion_matrices.tsv',confusions),('family_metrics.tsv',family_metrics),('family_confusion_matrix.tsv',families),('control_predictions.tsv',predictions)]:
        write_table(out/name,rows,list(rows[0]) if rows else ['policy','actual_admissible_families','predicted_family','count'])
    write_html(out/'performance_report.html',summary,metrics,confusions,family_metrics)
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
    return dict(status='complete',directory=str(out),summary_sha256=sha(out/'summary.json'),scored_controls=len(scored),
        primary_policy='dual_gate_supported',F1=main['F1'],PPV=main['PPV'],recall=main['recall'])


def after_audit(audit_dir,handoff):
    """Evaluation errors remain visible without invalidating completed audit calls."""
    out=Path(audit_dir)/'validation_evaluation'
    if out.exists():return dict(status='existing_validation_report_retained',directory=str(out))
    try:return evaluate(audit_dir,handoff,out)
    except Exception as error:
        import traceback
        out.mkdir(parents=True,exist_ok=True)
        result=dict(status='failed',error=str(error),directory=str(out))
        (out/'summary.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
        (out/'error.log').write_text(traceback.format_exc(),encoding='utf-8')
        return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--audit',type=Path,required=True);parser.add_argument('--handoff',type=Path)
    parser.add_argument('--out',type=Path,required=True);parser.add_argument('--target')
    parser.add_argument('--labels',type=Path,default=LABELS);parser.add_argument('--control-fasta',type=Path,default=CONTROL_FASTA)
    parser.add_argument('--version',action='version',version='adhesome-validation 1.0');args=parser.parse_args()
    from adhesome_workflow import verify
    audit=verify(args.audit,'audit')
    handoff=args.handoff or Path(json.loads((Path(audit['annotation'])/'workflow_manifest.json').read_text())['handoff'])
    verify(handoff,'handoff')
    if args.out.resolve()==args.audit.resolve() or args.audit.resolve() in args.out.resolve().parents:
        parser.error('Choose a separate report directory; standalone evaluation must not change a sealed audit')
    print(json.dumps(evaluate(args.audit,handoff,args.out,args.target,args.labels,args.control_fasta),indent=2))


if __name__=='__main__':main()
