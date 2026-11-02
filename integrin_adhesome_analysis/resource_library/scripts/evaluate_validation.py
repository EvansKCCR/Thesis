"""Version 1.0: read-only, post-audit evaluation of the labelled validation target."""

# Resolve source resources independently of the caller working directory.
import sys as _layout_sys
from pathlib import Path as _LayoutPath
_LAYOUT_BUNDLE_ROOT = next(p for p in _LayoutPath(__file__).resolve().parents
                          if (p/'config/layout_paths.json').is_file() and (p/'pipeline_paths.py').is_file())
if str(_LAYOUT_BUNDLE_ROOT) not in _layout_sys.path:
    _layout_sys.path.insert(0, str(_LAYOUT_BUNDLE_ROOT))
from pipeline_paths import activate as _layout_activate, resource_path as _layout_resource
_layout_activate()


import argparse
from collections import Counter,defaultdict
import csv
import hashlib
import html
import importlib.util
import json
from pathlib import Path
import re

from calibrate_identity import expected_family_resolution,expected_broad_families
from classification_contract import general_family

ROOT=_LAYOUT_BUNDLE_ROOT
LABELS=_layout_resource(ROOT/'validation/calibration_set.xlsx')
CONTROL_FASTA=_layout_resource(ROOT/'validation/calibration.faa')
POLICIES={
    'stage1_recruitment':'Exported stage 1 hypothesis, including provisional follow-up',
    'stage2_candidate_retention':'Stage 1 candidate present in the final audit after exclusions, including unresolved candidates; membership retention only, not family correctness',
    'audit_before_negative_exclusion':'Stage 1 plus unambiguous A/B call passing the complete family gate, before negative exclusions',
    'audit_without_exact_negative_exclusion':'Complete unambiguous A/B family gate with diagnostic negative exclusions retained and exact-sequence exclusions removed; counterfactual audit only',
    'dual_gate_supported':'Stage 1 plus unambiguous final A/B call passing the complete family gate, after negative exclusions',
    'dual_gate_strong':'Stage 1 plus a final strong grade A call passing the complete family gate, after negative exclusions',
    'broad_family_supported':'Stage 1 plus a supported broad or distinct named family at the declared resolution; no provisional or ambiguous promotion',
    'specific_family_supported':'Stage 1 plus a supported distinct named family at the declared resolution; does not establish an exact gene or paralogue; legacy unresolved resolution is unevaluable',
    'provisional_annotation':'Stage 1 plus a named provisional annotation, reported separately from supported calls',
    'ambiguous_annotation':'Stage 1 with competing family annotations, reported separately from supported and provisional calls',
    'audit_including_provisional':'Stage 1 plus any final A/B/C named family; exploratory coverage, not a passed full family gate'}


def sha(path):
    with _layout_resource(Path(path)).open('rb') as handle:return hashlib.file_digest(handle,'sha256').hexdigest()


def read_table(path,delimiter='\t'):
    with _layout_resource(Path(path)).open(encoding='utf-8-sig',newline='') as handle:return list(csv.DictReader(handle,delimiter=delimiter))


def write_table(path,rows,fields):
    with _layout_resource(Path(path)).open('w',encoding='utf-8',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=fields,delimiter='\t',extrasaction='ignore')
        writer.writeheader();writer.writerows(rows)


def records(path):
    header=None;parts=[]
    with _layout_resource(Path(path)).open(encoding='utf-8-sig') as handle:
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


def resource_hashes(family_rules=None,identity_rules=None):
    paths=[Path(__file__),LABELS,CONTROL_FASTA,_layout_resource(ROOT/'scripts/calibrate_identity.py'),_layout_resource(ROOT/'scripts/label_ontology.py'),_layout_resource(ROOT/'scripts/validation_reporting.py'),_layout_resource(ROOT/'scripts/classification_contract.py'),
        _layout_resource(ROOT/'annotation/family_identity_rules.json'),
        _layout_resource(ROOT/'validation/identity_probability_model.json'),
        _layout_resource(ROOT/'detection/config/family_rules.json'),_layout_resource(ROOT/'detection/references/adhesome_v1.0/search_catalogue.csv'),
        _layout_resource(ROOT/'detection/references/adhesome_v1.0/metazoan_catalogue.csv')]
    result={p.relative_to(ROOT).as_posix():sha(p) if _layout_resource(p).exists() else None for p in paths}
    if family_rules:result['detection/config/family_rules.json']=sha(family_rules)
    if identity_rules:result['annotation/family_identity_rules.json']=sha(identity_rules)
    return result


def load_gold(labels,control_fasta,family_rules=None,identity_rules=None):
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
    config=json.loads(_layout_resource(Path(family_rules or ROOT/'detection/config/family_rules.json')).read_text());catalogue=defaultdict(set);references=set();reverse=set()
    configured=set(config['families'])|{r['family'] for r in json.loads(_layout_resource(Path(identity_rules or ROOT/'annotation/family_identity_rules.json')).read_text())['rules']}
    for name,digests in [('search_catalogue.csv',references),('metazoan_catalogue.csv',reverse)]:
        for row in read_table(_layout_resource(ROOT/'detection/references/adhesome_v1.0')/name,','):
            digests.add(row['sequence_sha256'])
            if row['family'] in config['families']:
                for gene in re.split(r'[;,\s]+',row.get('gene','').upper()):catalogue[gene].add(general_family(row['family'],config['families']))
    for row in gold.values():
        resolution=expected_family_resolution(row,catalogue)
        row.update({key:value for key,value in resolution.items() if key!='families'})
        row['expected_families']=resolution['families']-{'other_named_identity_outside_current_family_catalogue'}
        row['expected_broad_families']=expected_broad_families(row,row['expected_families'])
        # Alternatives from an explicit label remain accepted at their declared
        # scope. A parent cannot substitute for a supplied finer label when the
        # distinct-family policy is being evaluated.
        from stage2_resolution import BROAD_LABELS
        specific=set(row['expected_families'])-set(BROAD_LABELS)
        row['expected_specific_families']=specific
        row['family_rule_coverage']='configured' if row['expected_families']&configured else 'outside_configured_family_gates' if row['expected_families'] else 'needs_mapping_review'
    return gold,trusted,references,reverse


def audit_rule_paths(audit_dir):
    """Evaluation uses the pair recorded by the completed Stage 2 invocation."""
    receipt=_layout_resource(Path(audit_dir)/'Stage_2_output.manifest.json')
    if not _layout_resource(receipt).exists():return None,None
    meta=json.loads(_layout_resource(receipt).read_text())
    if not meta.get('family_rules_path'):return None,None  # legacy/direct API output
    if not meta.get('identity_library_path'):
        raise ValueError('Audit has a family-rule path without its paired identity library')
    pair=_layout_resource(Path(meta['family_rules_path'])),_layout_resource(Path(meta['identity_library_path']))
    for path,expected in zip(pair,(meta.get('rule_sha256'),meta.get('identity_library_sha256'))):
        if not expected or sha(path)!=expected:
            raise ValueError('Audit configuration differs from its Stage 2 receipt: '+str(path))
    return pair


def read_audit(audit_dir,target):
    from openpyxl import load_workbook
    path=_layout_resource(audit_dir/'Stage_2_output.xlsx')
    if not _layout_resource(path).is_file():raise ValueError('A completed version 1.0 Stage_2_output.xlsx is required')
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
    module_path=_layout_resource(ROOT/'annotation/integrated_adhesome_pipeline.py')
    spec=importlib.util.spec_from_file_location('validation_audit_contract',module_path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    before={r['protein_id']:r for r in module.collapse_rows(screen)}
    final_path=_layout_resource(audit_dir/'audited_family_assignments.tsv')
    if not _layout_resource(final_path).exists():raise ValueError('The compact audited_family_assignments.tsv is required')
    final={}
    for row in read_table(final_path):
        pid=row['protein_id']
        if (row.get('Specie') or pid.split('__',1)[0])!=target:continue
        if pid in final:raise ValueError('Duplicate final audited candidate: '+pid)
        final[pid]=row
    exclusions={}
    for row in read_table(_layout_resource(audit_dir/'negative_control_exclusions.tsv')):
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


def resolution_decisions(row):
    """Reporting tiers never relax the full family gate or promote grade C."""
    if not row:return dict(broad_family_supported=False,specific_family_supported=False,
        provisional_annotation=False,ambiguous_annotation=False)
    level=str(row.get('family_support_level') or '').strip()
    resolution=str(row.get('annotation_resolution') or '').strip()
    family=row.get('putative_audited_family') or row.get('assigned_family') or 'unassigned'
    confidence=row.get('assignment_confidence') or row.get('family_confidence') or ''
    ambiguous=level=='ambiguous' or confidence=='ambiguous'
    complete=supported(row) and not ambiguous
    broad=complete and (not level or level in ('broad_supported','specific_supported'))
    if level:specific=complete and level=='specific_supported'
    elif not complete:specific=False
    elif resolution:specific=resolution in ('family','family_lineage','specific_family')
    else:specific=None
    provisional=(row.get('grade')=='C' or level=='provisional') and not ambiguous and family not in ('','unassigned','unresolved')
    return dict(broad_family_supported=broad,specific_family_supported=specific,
        provisional_annotation=provisional,ambiguous_annotation=ambiguous)


def aggregate(predictions):
    from validation_reporting import populations,policy_expected_families,policy_prediction
    scored=[row for row in predictions if row['evaluation_status']=='scored']
    metrics=[];confusions=[];family_metrics=[];family_cells=Counter()
    for policy in POLICIES:
        for population,cohort in populations(scored):
            unavailable=sum(row.get(policy) is None for row in cohort)
            cohort=[row for row in cohort if row.get(policy) is not None]
            counts=Counter()
            for row in cohort:
                actual=row['control_class']=='positive';called=truthy(row[policy])
                counts['TP' if actual and called else 'FN' if actual else 'FP' if called else 'TN']+=1
            result=dict(policy=policy,population=population,unevaluable_policy_controls=unavailable,**binary_metrics(*(counts[k] for k in ('TP','FP','TN','FN'))))
            metrics.append(result)
            for key,actual,predicted in [('TP','positive','accepted'),('FN','positive','not_accepted'),('FP','negative','accepted'),('TN','negative','not_accepted')]:
                confusions.append(dict(policy=policy,population=population,actual=actual,predicted=predicted,count=counts[key]))
        if policy in ('stage1_recruitment','stage2_candidate_retention','ambiguous_annotation'):continue
        tp=fp=fn=unmapped=0
        for row in scored:
            if row.get(policy) is None:continue
            actual_positive=row['control_class']=='positive';expected=policy_expected_families(row,policy)
            called=truthy(row[policy]);predicted=policy_prediction(row,policy)
            if actual_positive and not expected:unmapped+=1;continue
            if actual_positive:
                correct=called and predicted in expected
                tp+=correct;fn+=not correct;fp+=called and not correct
                actual=';'.join(sorted(expected))
            else:fp+=called;actual='__'+row['control_class']+'__'
            family_cells[policy,actual,predicted if called else '__not_accepted__']+=1
        from validation_reporting import macro_family_metrics
        macro,_=macro_family_metrics(scored,policy)
        family_metrics.append(dict(policy=policy,TP=tp,FP=fp,FN=fn,PPV=rate(tp,tp+fp),recall=rate(tp,tp+fn),
            F1=rate(2*tp,2*tp+fp+fn),positive_controls_without_mapped_family=unmapped,**macro))
    families=[dict(policy=p,actual_admissible_families=a,predicted_family=b,count=n) for (p,a,b),n in sorted(family_cells.items())]
    return metrics,confusions,family_metrics,families


def write_html(path,summary,metrics,confusions,family_metrics,reports=None):
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
    body+='<h2>Membership and disposition coverage</h2><p>PPV is precision. F1 = 2TP / (2TP + FP + FN). Rates are fractions from 0 to 1. Candidate retention, provisional annotation and ambiguity are coverage dispositions, not supported family acceptance. Specific supported means a distinct family at its declared scope, never an exact gene or paralogue.</p>'+table([r for r in metrics if r['population']=='all_controls'],['policy','TP','FP','TN','FN','PPV','recall','F1','specificity','accuracy','unevaluable_policy_controls'])
    body+='<h2>Confusion matrices</h2><p>Rows are actual classes; columns indicate the policy decision. For retention, provisional and ambiguous policies, accepted means that disposition was reported. A correct membership decision does not prove a correct family name.</p>'+''.join(matrices)
    body+='<h2>Family correctness at the reported scope</h2>'+table(family_metrics,['policy','TP','FP','FN','PPV','recall','F1','positive_controls_without_mapped_family','macro_precision','macro_recall','macro_F1'])
    body+='<p>A wrong positive-family call contributes one FP and one FN here. Admissible family labels come only from the trusted labels and ontology. This does not evaluate exact gene, isoform or paralogue names.</p>'
    body+='<p>Broad-family policy alone accepts the explicitly declared actin, Src/Csk and PET–LIM ontology ancestors. These are ontology roll-ups rather than independently curated new gold labels. Specific-family policy retains distinct gold labels. Macro scores give equal weight to each admissible gold-family set, keeping alternative labels together, plus predicted-only wrong-family groups. Undefined group denominators are omitted from that mean and their counts are exported in family_macro_metrics.tsv; per-group decisions are in family_per_group_metrics.tsv.</p>'
    if reports:
        body+='<h2>Negative-exclusion ablation</h2><p>Counterfactuals use the same saved family decisions. Removing exact matches retains diagnostic negative gates. Older single-gate records cannot establish this counterfactual and are counted as unevaluable. Removing all exclusions measures the saved classifier before those gates. Neither establishes independent specificity.</p>'+table(reports['exclusion_ablation.tsv'],['population','policy','n','unevaluable_policy_controls','TP','FP','TN','FN','PPV','recall','F1','specificity'])
        body+='<h2>Calibration and reference overlap</h2><p>Strata describe overlap with the saved training and reference resources. No-overlap strata are not independently retrained holdouts. Empty denominators remain n.a.</p>'+table([r for r in metrics if r['policy']=='dual_gate_supported' and r['population'] in ('identity_calibration_overlap','no_identity_calibration_overlap','no_forward_reference_overlap','no_reference_or_calibration_overlap')],['population','TP','FP','TN','FN','PPV','recall','F1'])
        body+='<h2>Distinct reporting tiers</h2>'+table(reports['confidence_tiers.tsv'],['control_class','tier','count'])
        body+='<h2>Family-label coverage</h2>'+table(reports['family_label_coverage.tsv'],['family_rule_coverage','positive_controls'])
        body+=table(reports['family_label_coverage_summary.tsv'],['positive_controls','mapped_expected_family_controls','unmapped_expected_family_controls','mapped_fraction','mapped_specific_family_controls','broad_ancestor_rollup_controls'])
        body+='<p>Trusted-name mapping can extend beyond the configured family gates. It does not add classifier rules or turn crosswalk components into supported calls. See family_mapping_review.tsv, missed_positive_review.tsv, family_misclassifications.tsv and per_expected_family.tsv.</p>'
        body+='<h2>Stage 2 gate diagnostics</h2><p>stage2_gate_diagnostics.tsv preserves the saved per-protein trace and reports only established reason codes. Missing or invalid traces are explicit; missing CDD support does not establish an unrecognized alias. Multiple limitation codes can apply to one candidate.</p>'+table(reports['stage2_reason_code_counts.tsv'],['control_class','reason_code','count'])
    body+='<h2>Interpretation and files</h2><p>Reference overlap: '+esc(summary['exact_forward_reference_overlap'])+' scored controls. Exact control exclusions and calibration overlap are reported per protein. Negative controls remain exclusion gates; post-exclusion specificity includes these configured rules.</p>'
    body+='<p>Undefined denominators are n.a. Reports never change assignments, confidence or evidence weights. See metrics.tsv, confusion_matrices.tsv, family_metrics.tsv, family_confusion_matrix.tsv and control_predictions.tsv for full data.</p>'
    style='body{font:15px system-ui;color:#15324a;margin:36px;max-width:1400px}table{border-collapse:collapse;margin:18px 0 32px;width:100%}th,td{border:1px solid #cedae4;padding:9px;text-align:right;overflow-wrap:anywhere}th:first-child,td:first-child{text-align:left}th{background:#e8f0f6}h2{margin-top:38px}p{line-height:1.5}'
    path.write_text('<!doctype html><html><head><meta charset="utf-8"><title>Validation performance</title><style>'+style+'</style></head><body>'+body+'</body></html>',encoding='utf-8')


def evaluate(audit_dir,handoff,out,target=None,labels=LABELS,control_fasta=CONTROL_FASTA):
    audit_dir=_layout_resource(Path(audit_dir)).resolve();handoff=_layout_resource(Path(handoff)).resolve();out=_layout_resource(Path(out)).resolve()
    if _layout_resource(out).exists():raise FileExistsError('Choose a new validation report directory: '+str(out))
    hand=json.loads((_layout_resource(handoff/'workflow_manifest.json')).read_text(encoding='utf-8'))
    targets=set(hand.get('target_species',[]))
    if not hand.get('detection'):
        if target is None:return dict(status='validation_target_not_present')
        raise ValueError('The handoff has no linked detection batch')
    detection=_layout_resource(Path(hand['detection']))
    if target is None:
        matches=[]
        for code in sorted(targets):
            config_path=_layout_resource(detection/code/'proteome.json')
            if not _layout_resource(config_path).exists():continue
            settings=json.loads(_layout_resource(config_path).read_text(encoding='utf-8'));entry=settings['proteomes'][0]
            if code.lower()=='validation' or _layout_resource(Path(entry['path'])).stem.lower()=='validation' or entry.get('species','').lower()=='pipeline validation':matches.append(code)
        if not matches:return dict(status='validation_target_not_present')
        if len(matches)>1:raise ValueError('Multiple validation targets; specify --target')
        target=matches[0]
    if target not in targets:raise ValueError('Validation target is not in the audited handoff: '+target)
    config_path=_layout_resource(detection/target/'proteome.json');settings=json.loads(_layout_resource(config_path).read_text(encoding='utf-8'))
    raw=_layout_resource(Path(settings['proteomes'][0]['path']))
    if not raw.is_absolute():raw=(_layout_resource(config_path.parent/raw)).resolve()
    prepared_dir=_layout_resource(detection/target/'prepared');preparation=_layout_resource(prepared_dir/'preparation.json')
    receipt=json.loads(_layout_resource(preparation).read_text(encoding='utf-8'))['species'][target]
    prepared=prepared_dir/(target+'.faa');mapping_path=prepared_dir/(target+'.id_map.csv')
    if sha(raw)!=receipt['input_sha256'] or sha(prepared)!=receipt['output_sha256']:raise ValueError('Validation input/prepared sequences changed after detection')
    for path in (prepared,_layout_resource(detection/target/'candidates/candidates.csv')):
        expected=hand.get('sources',{}).get(str(_layout_resource(path).resolve()))
        if expected is None or sha(path)!=expected:raise ValueError('Unverified validation detection source: '+str(path))
    selected_rules,selected_identity=audit_rule_paths(audit_dir)
    gold,trusted,references,reverse=load_gold(_layout_resource(Path(labels)),_layout_resource(Path(control_fasta)),selected_rules,selected_identity)
    mapping=read_table(mapping_path,',');by_original={r['original_id']:r for r in mapping}
    if len(by_original)!=len(mapping):raise ValueError('Duplicate validation original-ID mappings')
    sequences={h.split()[0]:s for h,s in records(prepared)}
    initial_hypotheses=defaultdict(set)
    for r in read_table(_layout_resource(handoff/'hypotheses.tsv')):
        if (r.get('Specie') or r['protein_id'].split('__',1)[0])==target:initial_hypotheses[r['protein_id']].add(r['family'])
    selected=set(initial_hypotheses)
    if not selected<=set(sequences):raise ValueError('A recruited validation candidate is absent from its prepared proteome')
    before,final,excluded=read_audit(audit_dir,target)
    if set(before)!=selected:raise ValueError('Validation audit does not account for every exported stage 1 candidate')
    model=_layout_resource(ROOT/'validation/identity_probability_model.json')
    training={r['sequence_sha256'] for r in json.loads(_layout_resource(model).read_text()).get('samples',[])} if _layout_resource(model).exists() else set()
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
        exact_excluded=excluded.get(pid,{}).get('exclusion_gate')=='exact_negative_control_sequence' if retained else False
        all_gates=excluded.get(pid,{}).get('exclusion_gates','') if retained else ''
        diagnostic_excluded='diagnostic_negative_subfamily_and_complete_architecture' in all_gates.split(';')
        exact_ablation_known=not exact_excluded or bool(all_gates)
        without_exact=prior if exact_excluded and not diagnostic_excluded else current
        row=dict(accession=key,protein_id=pid,original_id=original,sequence_sha256=digest,evaluation_status=status,
            label_sequence_verification=verification,
            control_class=actual.get('control_class',''),trusted_protein_name=actual.get('protein_name',''),
            expected_families=';'.join(sorted(expected)),preparation_decision=unit['decision'],
            expected_broad_families=';'.join(sorted(actual.get('expected_broad_families',set()))),
            expected_specific_families=';'.join(sorted(actual.get('expected_specific_families',set()))),
            trusted_genes=actual.get('genes',''),family_label_source=actual.get('family_label_source',''),family_rule_coverage=actual.get('family_rule_coverage',''),
            family_label_review_required=actual.get('family_label_review_required',''),
            family_label_review_reason=actual.get('family_label_review_reason',''),
            legacy_catalogue_families_ignored=actual.get('legacy_catalogue_families_ignored',''),
            stage1_recruitment=gate1,audit_before_negative_exclusion=gate1 and supported(prior),
            stage2_candidate_retention=gate1 and bool(current),
            audit_without_exact_negative_exclusion=(gate1 and supported(without_exact)) if exact_ablation_known else None,
            exact_exclusion_ablation_status='complete_gate_provenance' if exact_ablation_known else 'unavailable_legacy_single_gate_record',
            without_exact_exclusion_family=(without_exact.get('putative_audited_family') or without_exact.get('assigned_family') or 'unassigned'),
            dual_gate_supported=gate1 and supported(current),dual_gate_strong=gate1 and supported(current,strong=True),
            audit_including_provisional=gate1 and current.get('grade') in ('A','B','C') and final_family not in ('','unassigned','unresolved') and current.get('assignment_confidence')!='ambiguous',
            audited_family=final_family,pre_exclusion_family=prior.get('assigned_family','unassigned'),grade=current.get('grade',''),
            assignment_confidence=current.get('assignment_confidence',''),protein_name=current.get('protein_name',''),
            protein_name_scope=current.get('protein_name_scope',''),negative_exclusion=retained and pid in excluded,
            exclusion_gate=excluded.get(pid,{}).get('exclusion_gate','') if retained else '',
            exclusion_gates=all_gates,
            exact_forward_reference_overlap=digest in references,exact_reverse_reference_overlap=digest in reverse,
            identity_calibration_training_overlap=digest in training)
        row.update({name:(gate1 and decision) if decision is not None else None
                    for name,decision in resolution_decisions(current).items()})
        for field in ('identity_gate_trace','family_support_level','annotation_resolution'):
            row[field]=current.get(field) or prior.get(field) or ''
        for field in ('architecture_evidence','domain_score_0_4','orthology_score_0_3','identity_rule','identity_competitors','equally_ranked_hypotheses'):
            row['review_'+field]=prior.get(field,'')
        row['review_assignment_confidence']=prior.get('family_confidence','')
        row['stage1_family_hypotheses']=';'.join(sorted(initial_hypotheses.get(pid,set())))
        predictions.append(row)
    if input_ids!=set(by_original):raise ValueError('Validation input/preparation-map membership differs')
    scored=[r for r in predictions if r['evaluation_status']=='scored']
    if not scored:raise ValueError('No sequence-verified labelled controls in the validation target')
    metrics,confusions,family_metrics,families=aggregate(predictions)
    from validation_reporting import reports as diagnostic_reports
    reports=diagnostic_reports(predictions,metrics)
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
        family_label_coverage=reports['family_label_coverage.tsv'],
        family_label_coverage_summary=reports['family_label_coverage_summary.tsv'][0],
        gate_trace_coverage=dict(Counter(r['gate_trace_status'] for r in reports['stage2_gate_diagnostics.tsv'])),
        exact_exclusion_ablation_unavailable_controls=sum(r['audit_without_exact_negative_exclusion'] is None for r in scored),
        robustness_scope='Saved-output exclusion ablation, exact-sequence deduplication and reference/calibration overlap strata. No model retraining, reference removal or independent holdout inference performed.',
        calibration_training_overlap_is_acceptance_gate=False,
        identity_calibration_overlap_basis='Exact sequences in the saved training model, including historical overlap if its numerical estimates are stale/disabled; not a classification gate',
        negative_control_exclusions=sum(r['negative_exclusion'] for r in scored),primary_policy='dual_gate_supported',primary_metrics=main,
        policies=POLICIES,undefined_metric_policy='null in JSON, blank in TSV, n.a. in HTML; never fabricated zero',
        interpretation='Reference-panel control self-test with exact reference overlap, calibrated controls and configured negative exclusions. This is not independent held-out or novel-lineage performance. Membership and family correctness are reported separately; evaluation never changes calls or weights.')
    paths=[raw,prepared,mapping_path,preparation,config_path,_layout_resource(Path(labels)),_layout_resource(Path(control_fasta)),_layout_resource(handoff/'workflow_manifest.json'),_layout_resource(handoff/'hypotheses.tsv'),
        _layout_resource(detection/target/'candidates/candidates.csv'),_layout_resource(audit_dir/'Stage_2_output.xlsx'),_layout_resource(audit_dir/'audited_family_assignments.tsv'),_layout_resource(audit_dir/'negative_control_exclusions.tsv')]
    paths.extend(p for p in (selected_rules,selected_identity) if p is not None)
    summary['input_sha256']={str(p):sha(p) for p in paths};summary['evaluation_resources']=resource_hashes(selected_rules,selected_identity)
    out.mkdir(parents=True,exist_ok=False)
    for name,rows in [('metrics.tsv',metrics),('confusion_matrices.tsv',confusions),('family_metrics.tsv',family_metrics),('family_confusion_matrix.tsv',families),('control_predictions.tsv',predictions)]:
        write_table(out/name,rows,list(rows[0]) if rows else ['policy','actual_admissible_families','predicted_family','count'])
    empty_fields={'family_misclassifications.tsv':[*predictions[0],'review_reason'],
                  'missed_positive_review.tsv':[*predictions[0],'review_reason'],
                  'sequence_group_conflicts.tsv':['sequence_sha256','accessions','reason'],
                  'family_macro_metrics.tsv':['policy','macro_precision','macro_recall','macro_F1','macro_family_groups'],
                  'family_per_group_metrics.tsv':['policy','expected_admissible_families','TP','FP','FN','precision','recall','F1'],
                  'stage2_reason_code_counts.tsv':['control_class','reason_code','count'],
                  'exclusion_impact.tsv':['exclusion_gate','control_class','excluded','passed_family_gate_before_exclusion']}
    for name,rows in reports.items():write_table(out/name,rows,list(rows[0]) if rows else empty_fields.get(name,['status']))
    write_html(_layout_resource(out/'performance_report.html'),summary,metrics,confusions,family_metrics,reports)
    (_layout_resource(out/'summary.json')).write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
    return dict(status='complete',directory=str(out),summary_sha256=sha(_layout_resource(out/'summary.json')),scored_controls=len(scored),
        primary_policy='dual_gate_supported',F1=main['F1'],PPV=main['PPV'],recall=main['recall'])


def after_audit(audit_dir,handoff):
    """Evaluation errors remain visible without invalidating completed audit calls."""
    out=_layout_resource(Path(audit_dir)/'validation_evaluation')
    if _layout_resource(out).exists():return dict(status='existing_validation_report_retained',directory=str(out))
    try:return evaluate(audit_dir,handoff,out)
    except Exception as error:
        import traceback
        out.mkdir(parents=True,exist_ok=True)
        result=dict(status='failed',error=str(error),directory=str(out))
        (_layout_resource(out/'summary.json')).write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
        (_layout_resource(out/'error.log')).write_text(traceback.format_exc(),encoding='utf-8')
        return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--audit',type=Path,required=True);parser.add_argument('--handoff',type=Path)
    parser.add_argument('--out',type=Path,required=True);parser.add_argument('--target')
    parser.add_argument('--labels',type=Path,default=LABELS);parser.add_argument('--control-fasta',type=Path,default=CONTROL_FASTA)
    parser.add_argument('--version',action='version',version='adhesome-validation 1.0');args=parser.parse_args()
    from adhesome_workflow import verify
    audit=verify(args.audit,'audit')
    handoff=args.handoff or _layout_resource(Path(json.loads((Path(audit['annotation'])/'workflow_manifest.json').read_text())['handoff']))
    verify(handoff,'handoff')
    if _layout_resource(args.out).resolve()==_layout_resource(args.audit).resolve() or _layout_resource(args.audit).resolve() in _layout_resource(args.out).resolve().parents:
        parser.error('Choose a separate report directory; standalone evaluation must not change a sealed audit')
    print(json.dumps(evaluate(args.audit,handoff,args.out,args.target,args.labels,args.control_fasta),indent=2))


if __name__=='__main__':main()
