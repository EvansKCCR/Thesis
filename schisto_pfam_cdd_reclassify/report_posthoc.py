#!/usr/bin/env python3
"""Report saved post hoc proposals against reference labels; never run classification."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
import csv
import html
import json
from pathlib import Path
import re
import sys
import zipfile
import xml.etree.ElementTree as ET

from schisto_pfam_cdd_reclassify import load_xlsx, sha256, write_tsv
from stage2_io import SUPPORTED, REVIEW, excluded, integer_length, token
import orthogroup_reconstruction as orthogroup

ROOT = Path(__file__).resolve().parent
POLICIES = ('original_supported', 'posthoc_A', 'posthoc_AB')


def read_table(path, sheet='control_identity'):
    path = Path(path)
    if path.suffix.lower()=='.xlsx':
        fields, values = load_xlsx(path,sheet)
        rows = [dict(zip(fields,row)) for row in values if any(row)]
    else:
        with path.open(encoding='utf-8-sig',newline='') as handle:
            reader=csv.DictReader(handle,delimiter='\t' if path.suffix.lower()=='.tsv' else ',')
            fields, rows = reader.fieldnames or [], list(reader)
    if not fields or any(not f for f in fields) or len(fields)!=len(set(fields)):
        raise ValueError('Missing/duplicate input headers: '+str(path))
    if any(None in row for row in rows):
        raise ValueError('Row extends beyond input header: '+str(path))
    return rows, fields


def category(value):
    key=token(value).replace('_','')
    result={'positive':'positive','positivecontrol':'positive',
        'hardnegative':'hard_negative','hardnegativecontrol':'hard_negative',
        'backgroundnegative':'background_negative','backgroundnegativecontrol':'background_negative',
        'negative':'negative','negativecontrol':'negative'}.get(key)
    if result is None:raise ValueError('Unknown reference control class: '+str(value))
    return result


def ontology():
    """Use the existing trusted-name ontology for evaluation only; no model fit."""
    main=ROOT.parent/'integrin_adhesome_analysis'
    sys.path.insert(0,str(main))
    from pipeline_paths import activate
    activate()
    from calibrate_identity import expected_family_resolution, expected_broad_families
    scripts=main/'resource_library/scripts'
    return expected_family_resolution,expected_broad_families,{
        p.name:sha256(p) for p in [scripts/'calibrate_identity.py',scripts/'label_ontology.py']}


def normalized(name, config):
    value=str(name or '').strip()
    # Explicit equivalences only. SCR is deliberately absent.
    return config['family_aliases'].get(value,value)


def labels_and_mapping(labels, mapping_path, kind, label_sheet, config, resolver, broad_resolver):
    labels, _ = read_table(labels,label_sheet)
    if not labels:raise ValueError('Reference labels are empty')
    mappings, _ = read_table(mapping_path)
    records={};by_pair={};used=set()
    for original in labels:
        identifier=str(original.get('fasta_id') or original.get('protein_id') or
                       (original.get('accession') if kind=='calibration' else '') or '').strip()
        if not identifier or identifier in records:
            raise ValueError('Missing/duplicate reference label ID: '+identifier)
        cat=category(original.get('class') or original.get('control_class'))
        record=dict(original,protein_name=original.get('protein_name') or original.get('Protein_name') or '')
        # Expected identity comes only from the reference label, never predictions/mapping names.
        resolution=resolver(record) if cat=='positive' else {'families':set(),'family_label_source':'negative_control'}
        families={normalized(f,config) for f in resolution['families']}
        broad={normalized(f,config) for f in broad_resolver(record,resolution['families'])} if cat=='positive' else set()
        records[identifier]=dict(label_id=identifier,control_class=cat,expected=families,
            expected_broad=broad,truth_source=resolution['family_label_source'],reference=original)
    for row in mappings:
        pair=(str(row.get('Specie') or '').strip(),str(row.get('protein_id') or '').strip())
        if not all(pair) or pair in by_pair:
            raise ValueError('Missing/duplicate species/protein key in sequence manifest: '+str(pair))
        identity=str(row.get('original_id') or '').strip()
        if identity not in records and kind=='calibration':
            # Calibration table is keyed by accession; the explicit manifest binds it to the FASTA ID.
            identity=str(row.get('accession') or '').strip()
        if identity not in records:
            raise ValueError('Sequence manifest record has no exact reference label: '+str(pair))
        if identity in used:raise ValueError('A reference label maps to multiple prediction keys: '+identity)
        used.add(identity)
        gold=records[identity];reference=gold['reference']
        if row.get('control_class') and category(row['control_class'])!=gold['control_class']:
            raise ValueError('Control class differs between labels and manifest: '+identity)
        length=reference.get('sequence_length') or reference.get('length')
        if length is not None and str(length).strip() and row.get('length'):
            if integer_length(length)!=integer_length(row['length']):
                raise ValueError('Label/manifest length mismatch: '+identity)
        gold['pair']=pair
        gold['sequence_manifest_row']=dict(row,Specie=pair[0],protein_id=pair[1])
        by_pair[pair]=identity
    if used!=set(records):
        raise ValueError('Sequence manifest must cover every labelled ID; missing: '+', '.join(sorted(set(records)-used)[:8]))
    return records,by_pair


def load_refinement(folder):
    folder=Path(folder).resolve()
    receipt=folder/'run_manifest.json'
    manifest=json.loads(receipt.read_text(encoding='utf-8'))
    if manifest.get('status')!='complete' or not manifest.get('source_record_accounting_passed'):
        raise ValueError('A completed, accounted post hoc run is required')
    path=folder/'Stage_2_posthoc_refinement.tsv'
    if manifest.get('output_sha256',{}).get(path.name)!=sha256(path):
        raise ValueError('Refinement output differs from its saved SHA-256')
    rows, fields=read_table(path)
    needed={'posthoc_status','posthoc_proposed_family_or_class','posthoc_evidence_tier'}
    if not needed<=set(fields):raise ValueError('Input is not the all-record post hoc refinement table')
    mapping=manifest['column_mapping']
    result={}
    for original in rows:
        row=dict(original)
        for name in ('Specie','protein_id','family_support_level','putative_audited_family','length'):
            row[name]=original.get(mapping.get(name),'') if mapping.get(name) else ''
        pair=(row['Specie'].strip(),row['protein_id'].strip())
        if not all(pair) or pair in result:raise ValueError('Missing/duplicate refinement key: '+str(pair))
        level=token(row['family_support_level']);status=row['posthoc_status']
        if status not in {'reviewed','supported_preserved','negative_excluded'}:
            raise ValueError('Unknown post hoc status: '+status)
        if status=='supported_preserved' and (level not in SUPPORTED or excluded(row)):
            raise ValueError('Inconsistent supported audit scope: '+str(pair))
        if status=='reviewed' and (level not in REVIEW or excluded(row)):
            raise ValueError('Inconsistent refinement review scope: '+str(pair))
        proposal=row['posthoc_proposed_family_or_class']
        tier=row['posthoc_evidence_tier'][:1]
        if status!='reviewed' and (proposal or tier):
            raise ValueError('A preserved/excluded record cannot have a proposal: '+str(pair))
        if status=='reviewed' and tier not in {'A','B','C','D'}:
            raise ValueError('Missing/invalid proposal tier: '+str(pair))
        row['_tier']=tier
        result[pair]=row
    expected=manifest.get('unique_species_protein_pairs')
    if len(result)!=expected:raise ValueError('Post hoc record count differs from its manifest')
    return result, manifest, {'refinement_table':path,'refinement_manifest':receipt}


def rate(a,b):return a/b if b else None


def metrics(rows, policy):
    tp=fp=tn=fn=0;family_tp=family_fp=family_fn=unmapped=0
    for r in rows:
        positive=r['control_class']=='positive';called=r[policy]
        tp+=positive and called;fp+=not positive and called
        tn+=not positive and not called;fn+=positive and not called
        if positive and not r['expected_families']:
            unmapped+=1
            continue
        correct=positive and called and r[policy+'_family_correct']
        family_tp+=bool(correct);family_fp+=called and not correct
        family_fn+=positive and not correct
    binary=dict(TP=tp,FP=fp,TN=tn,FN=fn,n=len(rows),PPV=rate(tp,tp+fp),recall=rate(tp,tp+fn),
        F1=rate(2*tp,2*tp+fp+fn),specificity=rate(tn,tn+fp),accuracy=rate(tp+tn,len(rows)))
    family=dict(TP=family_tp,FP=family_fp,FN=family_fn,
        PPV=rate(family_tp,family_tp+family_fp),recall=rate(family_tp,family_tp+family_fn),
        F1=rate(2*family_tp,2*family_tp+family_fp+family_fn),
        evaluable_positive_labels=sum(r['control_class']=='positive' and bool(r['expected_families']) for r in rows),
        unmapped_positive_labels=unmapped)
    return binary,family


def evaluate(gold, predictions, config):
    results=[]
    for identifier,truth in gold.items():
        pair=truth['pair'];row=predictions.get(pair,{})
        status=row.get('posthoc_status','absent_from_audit')
        proposal=row.get('posthoc_proposed_family_or_class','')
        original_family=normalized(row.get('putative_audited_family'),config)
        proposal_family=normalized(config['proposal_labels'].get(proposal,''),config)
        supported=(status=='supported_preserved' and original_family not in ('','unassigned','unresolved'))
        tier=row.get('_tier','')
        eligible=status=='reviewed'
        result=dict(label_id=identifier,Specie=pair[0],protein_id=pair[1],control_class=truth['control_class'],
            expected_families=';'.join(sorted(truth['expected'])),expected_broad_families=';'.join(sorted(truth['expected_broad'])),
            truth_source=truth['truth_source'],posthoc_status=status,posthoc_eligible=eligible,
            original_audited_family=row.get('putative_audited_family',''),
            original_protein_name=row.get('protein_name',''),original_support=row.get('family_support_level',''),
            original_grade=row.get('grade',''),original_supported=bool(supported),
            original_evaluation_family=original_family,
            posthoc_proposal=proposal,posthoc_evaluation_family=proposal_family,posthoc_tier=tier,
            posthoc_A=eligible and tier=='A',posthoc_AB=eligible and tier in ('A','B'),
            proposal_label_mapping_status='mapped' if proposal_family else 'unmapped' if proposal else 'no_proposal')
        # Broad audit support is evaluated at its explicitly declared scope.
        original_truth=truth['expected'] | (truth['expected_broad'] if token(row.get('family_support_level'))=='broad_supported' else set())
        result['original_supported_family_correct']=supported and original_family in original_truth
        for policy in ('posthoc_A','posthoc_AB'):
            result[policy+'_family_correct']=result[policy] and bool(proposal_family) and proposal_family in truth['expected']
        results.append(result)
    return results


def display(value):
    if value is None:return 'not estimable'
    return f'{value:.1%}' if isinstance(value,float) else str(value)


def render(folder, kind, summary, matrices, family_matrices, recovery=None):
    title=kind.capitalize()+' post hoc performance — workflow 1.0'
    paragraphs=[
        'Original supported audit calls and post hoc proposals are scored separately. No combined or promoted audit result is created.',
        'Original-supported results use the full labelled cohort. Proposal results are shown both for that cohort and for the eligible unresolved/provisional/ambiguous subset. Missing audit candidates remain no-calls in the full-cohort denominator.',
        'Binary PPV/recall/F1 describe issuing a call or proposal for positive controls; family PPV/recall/F1 additionally require a correct mapped expected family. A wrong accepted family contributes both a family FP and FN. Positive labels without an expected family remain in binary metrics and are explicitly excluded from family metrics.',
        'Tier A and A+B proposals are exploratory domain-rule hypotheses, not supported audit assignments, calibrated probabilities, exact protein identities or demonstrated independent improvement. C/D proposals remain descriptive and are not scored as family calls.',
        'Expected families come only from trusted reference labels and the existing evaluation ontology. Proposal label equivalences are explicit. An unmapped proposal counts as incorrect when an evaluable family label exists. Negative labels are evaluation/exclusion controls and do not train a model.',
        'No model fitting, rule optimization, classification, orthology or phylogeny is run by this report. Calibration and previously inspected validation labels do not establish unseen-sequence accuracy.'
    ]
    sections=[]
    sections.append(('Record accounting',['Item','Count'],[[k,str(v)] for k,v in summary.items()]))
    for name, rows, fields in [('Binary performance',matrices,['population','policy','n','TP','FP','TN','FN','PPV','recall','F1','specificity']),
                              ('Family assignment performance',family_matrices,['population','policy','TP','FP','FN','PPV','recall','F1','evaluable_positive_labels','unmapped_positive_labels'])]:
        sections.append((name,fields,[[display(r.get(f)) for f in fields] for r in rows]))
    if recovery is not None:
        paragraphs.append('Orthogroup recovery requires an identical full-length amino acid sequence from the same-species reference proteome in the saved group. Protein identifiers only locate records; their agreement does not affect recovery. Identical sequences available in the panel but missing from the control group are counted separately. Controls without a mapped reference are not assessable.')
        paragraphs.append('Positive-control labels are applied after inference for this report. Because control sequences were supplied as targets, this checks recovery of known sequences in the saved groups. It does not establish reconstruction of unseen sequences or one-to-one orthology.')
        fields = ['species','positive_controls','assigned_orthogroup','controls_with_reference',
                  'any_same_species_member','identical_sequence_available','identical_sequence_recovered',
                  'identical_sequence_available_not_recovered','other_same_species_member',
                  'no_same_species_member','reference_not_available','recovery_among_with_reference',
                  'recovery_among_available_matches']
        headers = ['Species','Positive controls','In an orthogroup','With reference','Any same-species member',
                   'Identical sequence in panel','Identical sequence recovered','Available but not recovered',
                   'Other same-species member only','No same-species member or identical sequence in panel',
                   'No mapped reference','Recovered / with reference','Recovered / identical sequence in panel']
        sections.append(('Positive-control orthogroup reconstruction',headers,
                         [[display(r[f]) for f in fields] for r in recovery]))
    esc=lambda x:html.escape(str(x),quote=True)
    body='<h1>'+esc(title)+'</h1>'+''.join('<p>'+esc(p)+'</p>' for p in paragraphs)
    md=['# '+title,'',*paragraphs]
    for name,fields,rows in sections:
        body+='<h2>'+esc(name)+'</h2><div class="table"><table><thead><tr>'+''.join('<th>'+esc(f)+'</th>' for f in fields)+'</tr></thead><tbody>'
        body+=''.join('<tr>'+''.join('<td>'+esc(v)+'</td>' for v in row)+'</tr>' for row in rows)+'</tbody></table></div>'
        md+=['','## '+name,'','| '+' | '.join(fields)+' |','| '+' | '.join(['---']*len(fields))+' |']
        md+=['| '+' | '.join(str(v).replace('|','\\|') for v in row)+' |' for row in rows]
    links=['predictions.tsv','binary_metrics.tsv','family_metrics.tsv','family_confusion.tsv',
           'unmapped_labels.tsv','unmapped_proposals.tsv','family_misclassifications.tsv','unmatched_refinement_records.tsv','report_manifest.json']
    if recovery is not None:links.extend(['orthogroup_reconstruction.tsv','orthogroup_reconstruction_controls.tsv'])
    body+='<h2>Detailed tables and provenance</h2><ul>'+''.join('<li><a href="'+esc(p)+'">'+esc(p)+'</a></li>' for p in links)+'</ul>'
    css='body{font:15px/1.6 system-ui;background:#f5f7f9;color:#183043;margin:0}main{max-width:1250px;margin:24px auto;padding:32px;background:white;border-top:6px solid #176c78}h1,h2{color:#154a58}.table{overflow-x:auto}table{border-collapse:collapse;width:100%;font-size:13px}th,td{padding:8px;text-align:left;border-bottom:1px solid #dce5ea}th{background:#e9f2f5}tr:nth-child(even){background:#f7f9fb}@media print{main{margin:0;padding:0}.table{overflow:visible}}'
    (folder/'performance_report.html').write_text('<!doctype html><html lang="en"><meta charset="utf-8"><title>'+esc(title)+'</title><style>'+css+'</style><main>'+body+'</main></html>',encoding='utf-8')
    (folder/'performance_report.md').write_text('\n\n'.join(md)+'\n',encoding='utf-8')


def generate(args):
    folder=args.refinement.resolve();destination=args.out.resolve()
    if destination.exists():raise ValueError('Choose a new report output directory: '+str(destination))
    protected=[folder,args.sequence_manifest.resolve().parent]
    if bool(args.orthogroups)!=bool(args.panel_proteomes) or (args.reference_column and not args.orthogroups):
        raise ValueError('Orthogroup reconstruction requires both --orthogroups and --panel-proteomes')
    if args.orthogroups:
        protected.extend([args.orthogroups.resolve().parent,args.panel_proteomes.resolve()])
    if any(destination==p or p in destination.parents for p in protected):
        raise ValueError('Keep reports outside refinement and frozen-input directories')
    config_path=args.aliases.resolve()
    config=json.loads(config_path.read_text(encoding='utf-8'))
    if config.get('version')!='1.0' or not all(isinstance(config.get(k),dict) for k in ('family_aliases','proposal_labels')):
        raise ValueError('Invalid evaluation label configuration')
    if any(not isinstance(k,str) or not k.strip() or not isinstance(v,str) or not v.strip()
           for name in ('family_aliases','proposal_labels') for k,v in config[name].items()):
        raise ValueError('Evaluation label equivalences must be nonempty strings')
    sources={'labels':args.labels.resolve(),'sequence_manifest':args.sequence_manifest.resolve(),'evaluation_aliases':config_path}
    before={k:sha256(p) for k,p in sources.items()}
    predictions, manifest, refinement_sources=load_refinement(folder)
    sources.update(refinement_sources)
    before.update({k:sha256(p) for k,p in refinement_sources.items()})
    resolver,broad,ontology_hashes=ontology()
    gold,by_pair=labels_and_mapping(args.labels,args.sequence_manifest,args.kind,args.label_sheet,config,resolver,broad)
    recovery=None
    if args.orthogroups:
        orthogroups=args.orthogroups.resolve()
        controls=[r['sequence_manifest_row'] for r in gold.values() if r['control_class']=='positive']
        for r in controls:
            if not r.get('species') or not r.get('original_id') or not re.fullmatch(r'[0-9a-f]{64}',r.get('sequence_sha256') or ''):
                raise ValueError('Reconstruction requires species, original_id and sequence_sha256 in the positive-control manifest')
        reference_mapping, panel_files=orthogroup.panel_inputs(orthogroups,args.panel_proteomes.resolve(),controls,
            {pair[0] for pair in by_pair},args.reference_column)
        recovery_sources={'orthogroups':orthogroups,**{'panel:'+c:p for c,p in panel_files.items()}}
        sources.update(recovery_sources)
        before.update({k:sha256(p) for k,p in recovery_sources.items()})
        recovery, recovery_details=orthogroup.reconstruction(orthogroups,controls,reference_mapping,panel_files)
    if predictions and not predictions.keys() & by_pair.keys():
        raise ValueError('No refinement keys match this reference panel/sequence manifest')
    for pair,row in predictions.items():
        if pair not in by_pair or not row['length']:continue
        label=gold[by_pair[pair]]['reference']
        length=label.get('sequence_length') or label.get('length')
        if length is not None and str(length).strip() and integer_length(length)!=integer_length(row['length']):
            raise ValueError('Audit/reference length mismatch: '+str(pair))
    rows=evaluate(gold,predictions,config)
    unmatched=[dict(Specie=p[0],protein_id=p[1],reason='No labelled reference key; excluded from performance denominators')
               for p in predictions if p not in by_pair]
    summary=dict(labelled_records=len(rows),positive_controls=sum(r['control_class']=='positive' for r in rows),
        negative_controls=sum(r['control_class']!='positive' for r in rows),
        present_in_refinement=sum(r['posthoc_status']!='absent_from_audit' for r in rows),
        absent_from_original_audit=sum(r['posthoc_status']=='absent_from_audit' for r in rows),
        supported_audit_calls=sum(r['original_supported'] for r in rows),
        eligible_for_posthoc_review=sum(r['posthoc_eligible'] for r in rows),
        explicit_exclusions=sum(r['posthoc_status']=='negative_excluded' for r in rows),
        unmapped_positive_family_labels=sum(r['control_class']=='positive' and not r['expected_families'] for r in rows),
        unlabelled_refinement_records=len(unmatched))
    matrices=[];families=[];confusion=[];misclassified=[]
    for population,subset in [('full_labelled_cohort',rows),('posthoc_eligible',[r for r in rows if r['posthoc_eligible']])]:
        for policy in POLICIES:
            if population=='posthoc_eligible' and policy=='original_supported':continue
            binary,family=metrics(subset,policy)
            matrices.append(dict(population=population,policy=policy,**binary))
            families.append(dict(population=population,policy=policy,**family))
            if population!='full_labelled_cohort':continue
            counts=Counter()
            for r in subset:
                expected=';'.join(r['expected_families'].split(';')) if r['control_class']=='positive' else '__negative__'
                if not expected:expected='__unmapped_truth__'
                family_key='original_evaluation_family' if policy=='original_supported' else 'posthoc_evaluation_family'
                predicted=r[family_key] or '__unmapped_proposal__' if r[policy] else '__no_call__'
                counts[(expected,predicted)]+=1
                if r[policy] and not r[policy+'_family_correct'] and (r['control_class']!='positive' or r['expected_families']):
                    misclassified.append(dict(policy=policy,**r))
            confusion += [dict(policy=policy,expected_family=e,predicted_family=p,records=n) for (e,p),n in sorted(counts.items())]
    if any(sha256(p)!=before[k] for k,p in sources.items()):
        raise ValueError('A reporting input changed; no outputs written')
    destination.mkdir(parents=True,exist_ok=False)
    write_tsv(destination/'predictions.tsv',rows,list(rows[0]))
    write_tsv(destination/'binary_metrics.tsv',matrices,list(matrices[0]))
    write_tsv(destination/'family_metrics.tsv',families,list(families[0]))
    write_tsv(destination/'family_confusion.tsv',confusion,['policy','expected_family','predicted_family','records'])
    write_tsv(destination/'unmapped_labels.tsv',[r for r in rows if r['control_class']=='positive' and not r['expected_families']],list(rows[0]))
    write_tsv(destination/'unmapped_proposals.tsv',[r for r in rows if r['posthoc_AB'] and not r['posthoc_evaluation_family']],list(rows[0]))
    write_tsv(destination/'family_misclassifications.tsv',misclassified,['policy',*list(rows[0])])
    write_tsv(destination/'unmatched_refinement_records.tsv',unmatched,['Specie','protein_id','reason'])
    if recovery is not None:
        write_tsv(destination/'orthogroup_reconstruction.tsv',recovery,list(recovery[0]))
        write_tsv(destination/'orthogroup_reconstruction_controls.tsv',recovery_details,orthogroup.DETAIL_FIELDS)
    render(destination,args.kind,summary,matrices,families,recovery)
    receipt=dict(version='1.0',kind=args.kind,status='complete',created_utc=datetime.now(timezone.utc).isoformat(),
        source_files={k:str(p) for k,p in sources.items()},source_sha256=before,
        source_audit_sha256=manifest['input_sha256'],source_audit_path=manifest['input'],
        report_code_sha256=sha256(Path(__file__)),ontology_sha256=ontology_hashes,
        summary=summary,original_audit_and_proposals_separate=True,combined_policy_generated=False,
        classification_rerun=False,model_fitted=False,exact_protein_identity_scored=False,
        negative_controls_used_for_training=False,
        output_sha256={p.name:sha256(p) for p in sorted(destination.iterdir()) if p.is_file()})
    if recovery is not None:
        receipt['orthogroup_reconstruction']=dict(match_criterion='identical_full_length_amino_acid_sequence',
            reference_columns=reference_mapping,totals=recovery[-1],
            code_sha256=sha256(Path(orthogroup.__file__)))
    (destination/'report_manifest.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
    print('Report: '+str(destination/'performance_report.html'))
    return 0


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--refinement',type=Path,required=True,help='Completed post hoc output directory')
    parser.add_argument('--labels',type=Path,required=True,help='Reference calibration XLSX or validation CSV/TSV')
    parser.add_argument('--sequence-manifest',type=Path,required=True,help='Matching frozen branch inputs/sequence_manifest.tsv')
    parser.add_argument('--kind',choices=['calibration','validation'],required=True)
    parser.add_argument('--label-sheet',default='control_identity')
    parser.add_argument('--aliases',type=Path,default=ROOT/'config/evaluation_labels.json')
    parser.add_argument('--orthogroups',type=Path,help='Saved OrthoFinder Orthogroups.tsv for optional reconstruction counts')
    parser.add_argument('--panel-proteomes',type=Path,help='Matching saved panel/proteomes directory')
    parser.add_argument('--reference-column',action='append',default=[],metavar='SPECIES=COLUMN',
        help='Override the matching reference column; empty value marks a species not assessable')
    parser.add_argument('--out',type=Path,required=True,help='New report directory outside saved input/output directories')
    args=parser.parse_args(argv)
    try:return generate(args)
    except (OSError,ValueError,KeyError,zipfile.BadZipFile,ET.ParseError) as exc:
        print('ERROR: '+str(exc),file=sys.stderr)
        return 1


if __name__=='__main__':
    sys.exit(main())
