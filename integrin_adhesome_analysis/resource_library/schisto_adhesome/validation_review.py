"""Automatic terminal validation review. Reads saved calls; never runs inference."""

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
from collections import Counter, defaultdict
from datetime import datetime, timezone
import html
import json
import math
import os
from pathlib import Path
import re
import sys

from .datasets import category, dump, label_rows, sha
from .reports import read, metrics

POLICIES = {'stage1_recruited':'Stage 1 recruitment', 'supported':'Supported A/B',
            'provisional':'Provisional C only', 'including_provisional':'Supported plus provisional',
            'before_exclusion':'Supported before exclusions'}
FLAGS = tuple(POLICIES) + ('calibration_overlap','reference_overlap','training_overlap')


def local(value):
    value = str(value)
    if os.name == 'nt' and re.match(r'^/mnt/[a-z]/', value):
        return _layout_resource(Path(value[5].upper()+':/'+value[7:]))
    return _layout_resource(Path(value))


def boolean(value):
    if value in (True, 'True'):return True
    if value in (False, 'False'):return False
    raise ValueError('Invalid saved reporting boolean: '+str(value))


def predictions(path):
    rows = read(path)
    for row in rows:
        for key in FLAGS:
            if key in ('training_overlap','calibration_overlap') and row[key]=='unknown':continue
            row[key] = boolean(row[key])
    return rows


def family_metrics(rows, policy):
    tp=fp=fn=unmapped=0
    for row in rows:
        positive=row['control_class']=='positive'
        expected=set(row['expected_families'].split(';'))-{''}
        if positive and not expected:
            unmapped+=1;continue
        correct=positive and row[policy] and row['predicted_family'] in expected
        tp+=bool(correct);fn+=positive and not correct;fp+=row[policy] and not correct
    div=lambda a,b:a/b if b else None
    return dict(TP=tp,FP=fp,FN=fn,PPV=div(tp,tp+fp),recall=div(tp,tp+fn),
                F1=div(2*tp,2*tp+fp+fn),unmapped_positive_labels=unmapped)


def same(observed, expected):
    return observed is None and expected is None or (
        observed is not None and expected is not None and
        math.isclose(float(observed),float(expected),rel_tol=1e-10,abs_tol=1e-12))


def verify_reference(rows, labels):
    truth={str(row.get('fasta_id') or ''):row for row in labels}
    ids=[row['original_id'] for row in rows]
    if not rows or len(truth)!=len(labels) or len(ids)!=len(set(ids)) or set(ids)!=set(truth):
        raise ValueError('Validation review requires a unique, exact join of reference labels and predictions')
    for row in rows:
        reference=truth[row['original_id']]
        pairs=[(row['control_class'],category(reference.get('class'))),
               (int(row['length']),int(reference['sequence_length'])),
               (row.get('accession',''),str(reference.get('accession') or reference.get('mapped_accession') or '')),
               (row.get('genes',''),str(reference.get('genes') or '')),
               (row.get('protein_name',''),str(reference.get('protein_name') or reference.get('Protein_name') or ''))]
        if any(a!=b for a,b in pairs):
            raise ValueError('Validation reference differs from saved prediction: '+row['original_id'])


def verify_saved_metrics(rows, folder, summary):
    for saved in read(_layout_resource(folder/'metrics.tsv')):
        subset=rows if saved['population']=='all' else [r for r in rows if str(r[saved['population']])==saved['value']]
        calculated=metrics(subset,saved['policy'])
        if any(not same(value,None if saved[key]=='' else saved[key]) for key,value in calculated.items()):
            raise ValueError('Saved validation metric does not reconcile: '+saved['policy'])
    supported=metrics(rows,'supported');families=family_metrics(rows,'supported')
    if summary['records']!=len(rows) or any(not same(supported[k],v) for k,v in summary['supported'].items()):
        raise ValueError('Saved validation summary does not reconcile with predictions')
    for key in ('TP','FP','FN'):
        if summary['family'][key]!=families[key]:
            raise ValueError('Saved validation family summary does not reconcile')


class Document:
    def __init__(self):self.blocks=[]
    def heading(self,text):self.blocks.append(('heading',text))
    def paragraph(self,text):self.blocks.append(('paragraph',text))
    def table(self,headers,rows):
        rows=list(rows)
        if any(len(row)!=len(headers) for row in rows):raise ValueError('Invalid review table')
        self.blocks.append(('table',(headers,rows)))
    def links(self,items):self.blocks.append(('links',items))
    def save(self,folder,*,title='Validation report — workflow 1.0'):
        h=[];md=[];esc=lambda v:html.escape(str(v));cell=lambda v:str(v).replace('|','\\|').replace('\n',' ')
        for kind,value in self.blocks:
            if kind=='heading':
                level=1 if not h else 2
                h.append(f'<h{level}>{esc(value)}</h{level}>');md.append('#'*level+' '+value+'\n')
            elif kind=='paragraph':h.append('<p>'+esc(value)+'</p>');md.append(value+'\n')
            elif kind=='table':
                headers,rows=value
                h.append('<div class="table-wrap"><table><thead><tr>'+''.join('<th>'+esc(c)+'</th>' for c in headers)+'</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+esc(c)+'</td>' for c in row)+'</tr>' for row in rows)+'</tbody></table></div>')
                md.append('| '+' | '.join(map(cell,headers))+' |\n| '+' | '.join(['---']*len(headers))+' |\n'+'\n'.join('| '+' | '.join(map(cell,row))+' |' for row in rows)+'\n')
            else:
                h.append('<ul>'+''.join('<li><a href="'+html.escape(url,quote=True)+'">'+esc(label)+'</a></li>' for label,url in value)+'</ul>')
                md.append('\n'.join('- ['+label+']('+url+')' for label,url in value)+'\n')
        css='body{font:16px/1.55 system-ui,sans-serif;color:#183043;background:#f4f7fa;margin:0}main{max-width:1160px;margin:28px auto;background:white;padding:38px 44px;border-top:7px solid #176c78}h1{font-size:30px;color:#154a58}h2{font-size:22px;margin:34px 0 14px;border-top:1px solid #d9e4e9;padding-top:10px}p{max-width:100ch}a{color:#086a88}.table-wrap{overflow-x:auto;margin:18px 0}table{border-collapse:collapse;width:100%;font-size:14px}td,th{padding:10px 12px;border-bottom:1px solid #dce5ea;text-align:left;vertical-align:top}th{background:#e9f2f5}tr:nth-child(even){background:#f7f9fb}td:first-child{overflow-wrap:anywhere}@media(max-width:720px){main{margin:0;padding:22px 16px}}@media print{body{background:white}main{margin:0;padding:10px}h2{break-after:avoid}tr{break-inside:avoid}.table-wrap{overflow:visible}}'
        text='<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>'+esc(title)+'</title><style>'+css+'</style></head><body><main>'+''.join(h)+'</main></body></html>'
        (_layout_resource(folder/'validation_review.html')).write_text(text,encoding='utf-8')
        (_layout_resource(folder/'validation_review.md')).write_text('\n'.join(md),encoding='utf-8')


def generate(run, *, report_dir=None, audit_dir=None, calibration_reports=None, labels_file=None):
    """Generate from this run's reference, not a hard-coded validation cohort.

    The branch calls this after its audit. The scheduler calls it again once both
    parallel branches return, so a late calibration report is included safely.
    """
    run=_layout_resource(Path(run)).resolve();branch=_layout_resource(run/'validation')
    folder=_layout_resource(Path(report_dir)).resolve() if report_dir is not None else _layout_resource(branch/'reports')
    audit_dir=_layout_resource(Path(audit_dir)).resolve() if audit_dir is not None else _layout_resource(branch/'downstream/final')
    contract_path=_layout_resource(run/'run_contract.json');contract=json.loads(_layout_resource(contract_path).read_text())
    reference=contract['datasets']['validation'];label_file=local(labels_file or reference['labels_file'])
    if sha(label_file)!=reference['labels_sha256']:
        raise ValueError('Validation reference changed since the frozen run')
    sources=dict(reference=label_file,contract=contract_path,predictions=_layout_resource(folder/'predictions.tsv'),
                 metrics=_layout_resource(folder/'metrics.tsv'),summary=_layout_resource(folder/'summary.json'),
                 families=_layout_resource(folder/'family_metrics.tsv'),audit=_layout_resource(audit_dir/'audited_family_assignments.tsv'))
    for key,path in sources.items():
        if not _layout_resource(path).is_file():raise ValueError('Missing completed validation reporting input: '+str(path))
    before={key:sha(path) for key,path in sources.items()}
    rows=predictions(sources['predictions']);verify_reference(rows,label_rows(label_file))
    if reference.get('labelled_records',len(rows))!=len(rows) or (
            'counts' in reference and reference['counts']!=dict(Counter(r['control_class'] for r in rows))):
        raise ValueError('Validation cohort counts differ from the run contract')
    summary=json.loads(_layout_resource(sources['summary']).read_text());verify_saved_metrics(rows,folder,summary)
    audit_rows=read(sources['audit']);audit={row['protein_id']:row for row in audit_rows}
    exclusion_path=_layout_resource(audit_dir/'negative_control_exclusions.tsv')
    exclusion_rows=read(exclusion_path) if _layout_resource(exclusion_path).exists() else []
    excluded={row['protein_id']:row for row in exclusion_rows}
    if _layout_resource(exclusion_path).exists():sources['exclusions']=exclusion_path;before['exclusions']=sha(exclusion_path)
    recruited_ids={row['protein_id'] for row in rows if row['stage1_recruited']}
    if (len(audit)!=len(audit_rows) or len(excluded)!=len(exclusion_rows) or set(audit)&set(excluded)
            or set(audit)|set(excluded)!=recruited_ids):
        raise ValueError('Validation review audit candidate accounting differs from predictions')
    from evaluate_validation import supported
    for row in rows:
        actual=audit.get(row['protein_id'],{})
        name=actual.get('putative_audited_family') or actual.get('assigned_family') or 'unassigned'
        if supported(actual)!=row['supported'] or name!=row['predicted_family']:
            raise ValueError('Validation review decisions differ from saved audit: '+row['original_id'])
    positives=[r for r in rows if r['control_class']=='positive'];negatives=[r for r in rows if r['control_class']!='positive']
    binary={p:metrics(rows,p) for p in POLICIES};fm={p:family_metrics(rows,p) for p in ('supported','provisional','including_provisional')}
    mapped=[r for r in positives if r['expected_families']];unmapped=[r for r in positives if not r['expected_families']]
    review=[]
    for row in rows:
        if row['provisional'] and (row['control_class']!='positive' or row['expected_families'] and row['predicted_family'] not in row['expected_families'].split(';')):
            review.append(dict(original_id=row['original_id'],protein_id=row['protein_id'],species=row['species'],
                control_class=row['control_class'],genes=row['genes'],reference_protein_name=row['protein_name'],
                expected_families=row['expected_families'],predicted_family=row['predicted_family'],grade=row['grade'],
                reason='negative_control_provisionally_called' if row['control_class']!='positive' else 'reported_family_does_not_match_mapped_expected_family'))
    provisional_positive=[r for r in positives if r['provisional']]
    mismatches=sum(r['control_class']=='positive' for r in review)
    blind_root=_layout_resource(branch/'downstream/blind_reconstruction')
    blind_path=_layout_resource(blind_root/'results/blind_test_results.json')
    blind_rows_path=_layout_resource(blind_root/'results/blind_reconstruction_by_species.tsv')
    blind_rows=read(blind_rows_path) if _layout_resource(blind_rows_path).is_file() else []
    blind=json.loads(_layout_resource(blind_path).read_text()) if _layout_resource(blind_path).is_file() else dict(status='unavailable')
    if _layout_resource(blind_path).is_file():sources['blind_reference_reconstruction']=blind_path;before['blind_reference_reconstruction']=sha(blind_path)
    if _layout_resource(blind_rows_path).is_file():sources['blind_reference_reconstruction_by_species']=blind_rows_path;before['blind_reference_reconstruction_by_species']=sha(blind_rows_path)
    cal_folder=_layout_resource(Path(calibration_reports)).resolve() if calibration_reports is not None else _layout_resource(run/'calibration/reports')
    cal_state=_layout_resource(run/'calibration/branch_status.json');cal=None
    cal_ready=calibration_reports is not None or _layout_resource(cal_state).exists() and json.loads(_layout_resource(cal_state).read_text()).get('status')=='complete'
    if cal_ready and (_layout_resource(cal_folder/'predictions.tsv')).exists():
        sources['calibration_predictions']=_layout_resource(cal_folder/'predictions.tsv');before['calibration_predictions']=sha(sources['calibration_predictions'])
        cal=predictions(sources['calibration_predictions'])
    cal_status='available' if cal is not None else 'pending_or_unavailable'
    pct=lambda v:'n.a.' if v is None else f'{v:.1%}'
    dec=lambda v:'n.a.' if v is None else f'{v:.3f}'
    doc=Document();doc.heading('Validation report — workflow 1.0')
    doc.paragraph(f'Automatic evaluation against {label_file.name}: {len(rows):,} labelled proteins, '
                  f'{len(positives):,} positives and {len(negatives):,} negatives. This report reads saved calls only. '
                  'No inference, fitting or classification changes are performed by reporting.')
    doc.heading('Main findings')
    m=binary['supported'];recruited=binary['stage1_recruited']
    doc.paragraph(f'Stage 1 recruited {recruited["TP"]} of {len(positives)} positives ({pct(recruited["recall"])} recall) '
                  f'and {recruited["FP"]} negatives. Supported Stage 2 assignments retained {m["TP"]} positives and '
                  f'{m["FP"]} negatives: observed PPV {pct(m["PPV"])}, recall {pct(m["recall"])}, F1 {dec(m["F1"])}. '
                  'These are panel observations, not proof of precision on unseen proteins.')
    doc.paragraph(f'Of {len(provisional_positive)} positive provisional calls, {mismatches} have reported family labels '
                  f'that differ from their mapped expected family; {sum(not r["expected_families"] for r in provisional_positive)} '
                  'lack a mapped expected family. Membership recovery and family correctness are evaluated separately.')
    doc.heading('Reference population and accounting')
    grades=Counter(r['grade'] for r in rows)
    doc.table(['Item','Records'],[('Labelled evaluation',len(rows)),('Positive controls',len(positives)),
        ('Negative controls',len(negatives)),('Input FASTA',reference.get('fasta_records','n.a.')),
        ('Unlabelled input records excluded',len(reference.get('excluded_unlabelled_ids',[]))),
        ('Stage 1 candidates audited',len(recruited_ids)),('Candidates removed by exclusion gates',len(excluded))]+
        [(f'Grade {key}',value) for key,value in sorted(grades.items())])
    doc.paragraph('Every reference ID joins exactly once to saved predictions. Class, length, accession, genes and '
                  'reference names were checked; the label checksum matches the run contract. Unlabelled input '
                  'records are excluded from evaluation and are not counted as true negatives.')
    doc.heading('Tiered detection and classification performance')
    doc.table(['Policy','TP','FP','TN','FN','PPV','Recall','F1','Specificity','Accuracy'],
        [(POLICIES[p],m['TP'],m['FP'],m['TN'],m['FN'],pct(m['PPV']),pct(m['recall']),dec(m['F1']),
          pct(m['specificity']),pct(m['accuracy'])) for p,m in binary.items()])
    doc.paragraph('TP = called positive, FP = called negative, TN = negative not called, FN = positive not called. '
                  'PPV = TP/(TP+FP), recall = TP/(TP+FN), F1 = 2TP/(2TP+FP+FN). Undefined rates are n.a. '
                  'Supported means a named family, grade A/B, domain score 4/4 and a non-ambiguous assignment. '
                  'Provisional means grade C with a named family. Reference labels are used for terminal evaluation only.')
    m=binary['supported'];doc.table(['Actual class','Supported accepted','Not accepted'],[('Positive',m['TP'],m['FN']),('Negative',m['FP'],m['TN'])])
    doc.heading('Family assignment correctness')
    doc.table(['Family scoring policy','Correct (TP)','Incorrect (FP)','Unrecovered (FN)','Micro PPV','Micro recall','Micro F1'],
        [(POLICIES[p],m['TP'],m['FP'],m['FN'],pct(m['PPV']),pct(m['recall']),dec(m['F1'])) for p,m in fm.items()])
    doc.paragraph(f'{len(mapped)} positives have mapped expected families; {len(unmapped)} do not. '
                  f'{sum(r["supported"] for r in unmapped)} supported positives have no mapped family truth. '
                  'Those unmapped positives remain in binary metrics but not family scoring. A wrong called family '
                  'on a mapped positive contributes one FP and one FN; a called negative contributes one FP. '
                  'Alternative expected families are one admissible label set, not duplicated cases. Exact-family '
                  'mismatches can include differences in granularity and require review before a biological claim.')
    family=summary['family']
    doc.paragraph(f'Supported family macro precision is {pct(family.get("macro_precision"))} over '
                  f'{family.get("macro_precision_denominator",0)} defined groups; macro recall '
                  f'{pct(family.get("macro_recall"))} and macro F1 {dec(family.get("macro_F1"))} cover '
                  f'{family.get("macro_recall_denominator",0)} groups. Zero-call groups have undefined precision. '
                  'Expected families come from the run’s saved trusted-identity mapping, not predicted recruitment '
                  'families. These metrics do not independently verify specific protein names or probability calibration.')
    families=read(sources['families']);worst=sorted(families,key=lambda r:(-int(r['FN']),r['family']))[:15]
    doc.table(['Expected family/group','Correct supported','Missed/incorrect','Recall'],
              [(r['family'],r['TP'],r['FN'],pct(float(r['recall'])) if r['recall'] else 'n.a.') for r in worst])
    grouped=defaultdict(list)
    for row in unmapped:grouped[(row['accession'],row['protein_name'])].append(row)
    doc.table(['Unmapped reference identity','Accession','Positive records','Supported'],
              [(name,acc,len(rs),sum(r['supported'] for r in rs)) for (acc,name),rs in grouped.items()])
    doc.heading('Comparison with calibration')
    if cal is None:doc.paragraph('Calibration comparison is pending or unavailable. The scheduler refreshes this report '
        'after the parallel calibration branch finishes; validation reporting does not launch or refit calibration.')
    else:
        doc.table(['Policy','Calibration PPV','Validation PPV','Calibration recall','Validation recall','Calibration F1','Validation F1'],
            [(POLICIES[p],pct(metrics(cal,p)['PPV']),pct(binary[p]['PPV']),pct(metrics(cal,p)['recall']),pct(binary[p]['recall']),
              dec(metrics(cal,p)['F1']),dec(binary[p]['F1'])) for p in ('stage1_recruited','supported','including_provisional')])
        doc.paragraph(f'Calibration contains {len(cal)} controls, including {sum(r["control_class"]=="positive" for r in cal)} positives. '
                      'Different populations, negative classes and species make this a descriptive comparison, not a paired improvement estimate.')
    doc.heading('Performance by species')
    species_rows=[];mixed_species=False
    for species in sorted({r['species'] for r in rows}):
        subset=[r for r in rows if r['species']==species];m=metrics(subset,'supported')
        pos=sum(r['control_class']=='positive' for r in subset);neg=len(subset)-pos;mixed_species|=bool(pos and neg)
        species_rows.append((species,pos,neg,sum(r['stage1_recruited'] for r in subset),m['TP'],m['FP'],m['FN'],pct(m['recall']),pct(m['specificity'])))
    doc.table(['Species','Positives','Negatives','Recruited','TP','FP','FN','Recall','Specificity'],species_rows)
    if not mixed_species:doc.paragraph('No species has both positive and negative controls. Species and class are confounded; '
        'sensitivity and specificity cannot both be estimated within any one species from this panel.')
    doc.heading('Exact-sequence overlap strata')
    strata=[]
    for label,key in [('Calibration','calibration_overlap'),('Forward references','reference_overlap'),('Identity training','training_overlap')]:
        for value in (False,True,'unknown'):
            subset=[r for r in rows if r[key]==value];m=metrics(subset,'supported')
            if value=='unknown' and not subset:continue
            strata.append((label,'Unknown: comparison input unavailable' if value=='unknown' else 'Exact overlap' if value else 'No exact overlap',sum(r['control_class']=='positive' for r in subset),
                sum(r['control_class']!='positive' for r in subset),m['TP'],m['FP'],pct(m['recall']),dec(m['F1'])))
    clean=[r for r in rows if not any(r[key] for key in FLAGS[-3:])];clean_metrics=metrics(clean,'supported')
    strata.append(('All three','No exact overlap',sum(r['control_class']=='positive' for r in clean),sum(r['control_class']!='positive' for r in clean),
                   clean_metrics['TP'],clean_metrics['FP'],pct(clean_metrics['recall']),dec(clean_metrics['F1'])))
    doc.table(['Comparison','Stratum','Positives','Negatives','TP','FP','Recall','F1'],strata)
    doc.paragraph('These strata exclude exact sequence matches for the stated comparison. They do not exclude '
                  'reference homologues or establish independence from rule development. The combined row excludes all three exact-overlap types.')
    doc.heading('Missed positives, false positives and exclusion ablation')
    missed=[r for r in positives if not r['supported']];missed_grades=Counter(r['grade'] for r in missed)
    doc.table(['Unsupported positive disposition','Records'],sorted(missed_grades.items()))
    doc.table(['Positive missed by Stage 1','Reference name'],[(r['original_id'],r['protein_name']) for r in positives if not r['stage1_recruited']])
    doc.table(['Recruited negative','Predicted family','Grade','Supported','Provisional'],[(r['original_id'],r['predicted_family'],r['grade'],r['supported'],r['provisional']) for r in negatives if r['stage1_recruited']])
    before_m=binary['before_exclusion'];after_m=binary['supported']
    doc.paragraph(f'Supported calls before exclusions: TP {before_m["TP"]}, FP {before_m["FP"]}; after exclusions: '
                  f'TP {after_m["TP"]}, FP {after_m["FP"]}. Exclusion gates were recorded for '
                  f'{sum(bool(r.get("exclusion_gate")) for r in rows)} validation records. This ablation reports saved decisions, not a new inference run.')
    regression=[]
    for gene in ('PTPRF','NF2','TRIO'):
        subset=[r for r in positives if gene in re.split(r'[,;| /]+',r['genes'].upper())]
        if subset:regression.append((gene,len(subset),sum(r['supported'] for r in subset),sum(r['provisional'] for r in subset),
            '; '.join(f'{name} ({count})' for name,count in sorted(Counter(r['predicted_family'] for r in subset).items()))))
    if regression:doc.table(['Previously flagged identity','Positives','Supported','Provisional','Predicted family counts'],regression)
    doc.heading('Blind reference-only orthogroup reconstruction')
    doc.paragraph('Assessment status: '+str(blind.get('status','unavailable'))+'. OrthoFinder inference used only the prepared reference proteomes; '
                  'positive controls were withheld from inference and matched afterward by exact full-length sequence identity. '
                  'This measures reconstruction availability, not family-classification accuracy.')
    if blind_rows:
        doc.table(['Species','Positive controls','Identity assessable','One orthogroup','Multiple orthogroups',
                   'Identical reference, ungrouped','No identical reference'],
            [(r.get('species'),r.get('positive_controls'),r.get('identity_assessable_controls'),
              r.get('identical_sequence_maps_to_one_orthogroup'),r.get('identical_sequence_maps_to_multiple_orthogroups'),
              r.get('identical_reference_sequence_not_assigned_to_group'),r.get('no_identical_sequence_in_reference_panel'))
             for r in blind_rows])
        doc.paragraph(f"Reference-only orthogroup table hash: {blind.get('orthogroups_sha256','unavailable')}. "
                      f"Reconstructed {blind.get('recovered','n.a.')} of {blind.get('positive_controls','n.a.')} positive controls into exactly one orthogroup.")
    else:
        doc.paragraph('The blind assessment did not complete. Audit classification and its performance metrics remain available without it.')
    doc.heading('Interpretation and follow-up')
    doc.paragraph('Review unsupported positives and the provisional-family mismatch table separately. Recruitment '
                  'misses need a recruitment-stage investigation; recovered membership does not establish the reported '
                  'family identity. Curate unmapped expected families without copying predictions into truth. '
                  'Retain distinct confidence tiers and do not relax gates solely to improve these scores.')
    doc.paragraph('Exact reference overlap, species/class composition and incomplete family truth '
                  'limit transferability claims. No experimental-evidence requirement or suppression is introduced by this report.')
    doc.heading('Detailed outputs and provenance')
    items=[('All predictions',_layout_resource(folder/'predictions.tsv')),('Metric strata',_layout_resource(folder/'metrics.tsv')),('Per-family metrics',_layout_resource(folder/'family_metrics.tsv')),
        ('Unsupported positives and errors',_layout_resource(folder/'errors.tsv')),('Pipeline summary report',_layout_resource(folder/'report.html')),
        ('Validation reference',label_file),('Final audit workbook',_layout_resource(audit_dir/'Stage_2_output.xlsx'))]
    if cal is not None:items.append(('Calibration report',_layout_resource(cal_folder/'report.html')))
    if _layout_resource(blind_rows_path).is_file():items.append(('Blind reconstruction by species',blind_rows_path))
    if _layout_resource(blind_root/'results/blind_reconstruction_controls.tsv').is_file():items.append(('Blind reconstruction controls',_layout_resource(blind_root/'results/blind_reconstruction_controls.tsv')))
    links=[]
    for label,path in items:
        if _layout_resource(path).is_file():
            try:url=_layout_resource(Path(os.path.relpath(path,folder))).as_posix()
            except ValueError:url=path.as_uri()
            links.append((label,url))
    links.extend([('Provisional family review','validation_provisional_family_review.tsv'),('Reporting provenance','validation_review_provenance.json')])
    doc.links(links)
    doc.paragraph('Original analytical outputs remain unchanged. This terminal report uses the recorded label reference and '
                  'saved predictions. It creates no backups and launches no pipeline stages.')
    for key,path in sources.items():
        if sha(path)!=before[key]:raise ValueError('Reporting source changed while rendering: '+str(path))
    doc.save(folder)
    fields=['original_id','protein_id','species','control_class','genes','reference_protein_name','expected_families','predicted_family','grade','reason']
    # Keep a valid header even if no provisional errors are present.
    with (_layout_resource(folder/'validation_provisional_family_review.tsv')).open('w',encoding='utf-8',newline='') as handle:
        import csv
        writer=csv.DictWriter(handle,fields,delimiter='\t');writer.writeheader();writer.writerows(review)
    provenance=dict(pipeline_version='1.0',created_at_utc=datetime.now(timezone.utc).isoformat(),
        reference=str(label_file),labelled_records=len(rows),positive_records=len(positives),negative_records=len(negatives),
        source_files={key:dict(path=str(path),sha256=before[key]) for key,path in sources.items()},
        checks=dict(unique_exact_label_join=True,labels_match_frozen_contract=True,labels_lengths_accessions_names_match=True,
                    final_audit_decisions_match_predictions=True,all_metric_strata_recomputed=True,sources_unchanged=True),
        calibration_comparison=cal_status,supported=binary['supported'],family_by_policy=fm,
        combined_no_exact_overlap=clean_metrics,provisional_family_review_records=len(review),
        blind_reference_reconstruction_status=blind.get('status','unavailable'),analytical_pipeline_rerun=False,backups_created=False)
    dump(_layout_resource(folder/'validation_review_provenance.json'),provenance)
    return dict(report=str(_layout_resource(folder/'validation_review.html')),markdown=str(_layout_resource(folder/'validation_review.md')),
                provisional_review=str(_layout_resource(folder/'validation_provisional_family_review.tsv')),calibration_comparison=cal_status)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run',type=Path,required=True)
    args=parser.parse_args()
    try:print(json.dumps(generate(args.run),indent=2));return 0
    except (OSError,ValueError,KeyError) as exc:
        print('ERROR: '+str(exc),file=sys.stderr);return 1


if __name__=='__main__':raise SystemExit(main())
