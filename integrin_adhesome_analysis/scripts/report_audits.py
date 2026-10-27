#!/usr/bin/env python3
"""Report completed calibration/validation audits against their original labels.

Reads saved calls only: no inference, fitting, backups or changes to source runs.
"""

# Resolve the organized bundle independently of the caller working directory.
import sys as _bundle_sys
from pathlib import Path as _BundlePath
_BUNDLE_ROOT=next(p for p in _BundlePath(__file__).resolve().parents if (p/'pipeline_paths.py').is_file())
if str(_BUNDLE_ROOT) not in _bundle_sys.path:_bundle_sys.path.insert(0,str(_BUNDLE_ROOT))
from pipeline_paths import activate as _activate_bundle_paths
_activate_bundle_paths()
import argparse
from collections import Counter
from datetime import datetime,timezone
import json
from pathlib import Path
import sys

import adhesome_workflow as core
from schisto_adhesome.datasets import ROOT,dump,load_dataset,sha,table
from schisto_adhesome.reports import read,evaluate_branch,summarize_predictions
from schisto_adhesome.validation_review import Document,family_metrics,generate as validation_review,local


def checked_cohort(run,kind,contract,labels_override=None,allow_equivalent_calibration_label_bytes=False):
    """Recheck original truth and sequence hashes without preparing new inputs."""
    reference=contract['datasets'][kind]
    labels=local(labels_override or reference['labels_file']).resolve()
    protein_file=local(reference['protein_file']).resolve()
    if sha(labels)!=reference['labels_sha256'] and not (kind=='calibration' and allow_equivalent_calibration_label_bytes):
        raise ValueError(kind+': reference labels differ from the saved run')
    if sha(protein_file)!=reference['fasta_sha256']:raise ValueError(kind+': reference FASTA differs from the saved run')
    trusted,summary=load_dataset(kind,protein_file,labels,select_labelled=True)
    if (summary['counts']!=reference['counts'] or summary['labelled_records']!=reference['labelled_records']
            or summary['excluded_unlabelled_ids']!=reference.get('excluded_unlabelled_ids',[])):
        raise ValueError(kind+': reference population differs from the saved run')
    path=run/kind/'inputs/sequence_manifest.tsv';manifest=read(path)
    by_id={r['original_id']:r for r in trusted}
    if len(manifest)!=len(by_id) or len({r['protein_id'] for r in manifest})!=len(manifest) or set(by_id)!={r['original_id'] for r in manifest}:
        raise ValueError(kind+': sequence manifest must join every labelled ID exactly once')
    fields=('original_id','species','control_class','accession','genes','protein_name','trusted_family','sequence_sha256','length')
    for row in manifest:
        if any(str(row[key])!=str(by_id[row['original_id']][key]) for key in fields):
            raise ValueError(kind+': sequence/label metadata differ for '+row['original_id'])
    return manifest,labels,protein_file,path


def checked_audit(run,kind,audit):
    """Bind the standalone audit to the same annotation and handoff as its cohort."""
    meta=core.verify(audit,'audit')
    annotation=run/kind/'annotation';handoff=run/kind/'handoff'
    annotation_meta=core.verify(annotation,'annotation',required_files=[])
    if (meta.get('annotation_manifest')!=sha(annotation/'workflow_manifest.json')
            or local(meta['annotation']).resolve()!=annotation.resolve()
            or local(annotation_meta['handoff']).resolve()!=handoff.resolve()
            or annotation_meta.get('handoff_manifest')!=sha(handoff/'workflow_manifest.json')):
        raise ValueError(kind+': audit belongs to a different annotation or Stage 1 handoff')
    core.verify(handoff,'handoff',required_files=['hypotheses.tsv'])
    return meta


def check_output(out,run,audits):
    for source in (run,*audits):
        if out==source or source in out.parents or out in source.parents:
            raise ValueError('Report output must be separate from source run and audit directories')
    if out.exists():raise ValueError('Choose a new report directory: '+str(out))


def generate(args):
    run=args.run.resolve();out=args.out.resolve()
    audits={kind:getattr(args,kind+'_audit').resolve() for kind in ('calibration','validation')}
    check_output(out,run,list(audits.values()))
    contract_path=run/'run_contract.json';contract=json.loads(contract_path.read_text())
    sources={str(contract_path):sha(contract_path)};cohorts={};audit_meta={};label_checks={}
    def record(path):sources[str(path)]=sha(path)
    # All data/provenance checks precede creating reports. Do not validate the old
    # algorithm inventory: newer completed audits deliberately use newer rules.
    for kind in ('calibration','validation'):
        manifest,labels,fasta,manifest_path=checked_cohort(run,kind,contract,getattr(args,kind+'_labels',None),
            getattr(args,'allow_equivalent_calibration_label_bytes',False))
        label_checks[kind]=dict(frozen_labels_sha256=contract['datasets'][kind]['labels_sha256'],
            current_labels_sha256=sha(labels),byte_identical=sha(labels)==contract['datasets'][kind]['labels_sha256'],
            all_evaluation_fields_match_frozen_manifest=True,reference_fasta_byte_identical=True)
        audit_meta[kind]=checked_audit(run,kind,audits[kind]);cohorts[kind]=(manifest,labels)
        for path in (labels,fasta,manifest_path,audits[kind]/'workflow_manifest.json',
                     run/kind/'annotation/workflow_manifest.json',run/kind/'handoff/workflow_manifest.json',
                     run/kind/'handoff/hypotheses.tsv',run/kind/'reference_overlap.tsv'):
            record(path)
    references=set()
    for name in ('search_catalogue.csv','metazoan_catalogue.csv'):
        path=ROOT/'detection/references/adhesome_v1.0'/name;record(path)
        references.update(r['sequence_sha256'] for r in core.read(path,','))
    model_path=ROOT/'validation/identity_probability_model.json';record(model_path)
    training={r['sequence_sha256'] for r in json.loads(model_path.read_text()).get('samples',[])}
    for path in (ROOT/'scripts/calibrate_identity.py',ROOT/'scripts/label_ontology.py',ROOT/'schisto_adhesome/reports.py',Path(__file__).resolve()):record(path)
    overlap={r['sequence_sha256'] for r in cohorts['calibration'][0]}
    out.mkdir(parents=True,exist_ok=False);results={};predictions={};metrics_rows=[]
    for kind in ('calibration','validation'):
        branch=run/kind;manifest,labels=cohorts[kind]
        ref_hashes=references|{r['sequence_sha256'] for r in read(branch/'reference_overlap.tsv')}
        rows=evaluate_branch(branch,manifest,audits[kind],overlap,ref_hashes,training)
        folder=out/kind;summary=summarize_predictions(rows,folder,kind)
        # Reports describe fixed saved audits; reporting does not select or fit.
        summary.update(interpretation='Terminal evaluation of a completed audit using reference labels verified against the frozen cohort; no fitting or inference.',
            audit=str(audits[kind]),reference_labels=str(labels),stage2_profile=audit_meta[kind].get('stage2_profile','original'),
            labels_sha256=sha(labels),excluded_unlabelled_records=len(contract['datasets'][kind].get('excluded_unlabelled_ids',[])))
        dump(folder/'summary.json',summary)
        from schisto_adhesome.reports import render
        render(folder/'report.html',summary,read(folder/'metrics.tsv'))
        results[kind]=summary;predictions[kind]=rows
        metrics_rows.append(dict(dataset=kind,**summary['supported'],**{'family_'+k:v for k,v in family_metrics(rows,'supported').items()}))
    review=validation_review(run,report_dir=out/'validation',audit_dir=audits['validation'],
        calibration_reports=out/'calibration',labels_file=cohorts['validation'][1])
    table(out/'performance_comparison.tsv',metrics_rows)
    # Retain a human-readable calibration review, including negative subgroups,
    # separate family correctness and grade tiers, rather than only JSON metrics.
    doc=Document();doc.heading('Calibration report — completed optimized audit')
    doc.paragraph('Evaluation against '+cohorts['calibration'][1].name+'. Labels were joined through the saved sequence manifest. '
        'No biological analysis, fitting, threshold selection or change to the audit was performed. Calibration is an internal self-test, '
        'not evidence of independent accuracy. Negative controls are evaluated as exclusions, never protein families.')
    doc.table(['Control class','Records'],sorted(Counter(r['control_class'] for r in predictions['calibration']).items()))
    for title,path,fields in [('Supported / provisional / recruitment tiers',out/'calibration/metrics.tsv',
        ['policy','TP','FP','TN','FN','PPV','recall','F1','specificity']),
        ('Family correctness',out/'calibration/family_summary.tsv',
        ['TP','FP','FN','micro_precision','micro_recall','micro_F1','unmapped_positive_labels'])]:
        data=read(path)
        if title.startswith('Supported'):data=[r for r in data if r['population']=='all']
        doc.heading(title);doc.table(fields,[[r[k] for k in fields] for r in data])
    doc.heading('Negative-control performance')
    data=[r for r in read(out/'calibration/metrics.tsv') if r['population']=='control_class' and r['policy']=='supported']
    fields=['value','TP','FP','TN','FN','specificity'];doc.table(fields,[[r[k] for k in fields] for r in data])
    doc.links([('Full metric strata','metrics.tsv'),('Family metrics','family_metrics.tsv'),('Candidate predictions','predictions.tsv'),
        ('Errors and missed positives','errors.tsv'),('Report provenance','../reporting_provenance.json')])
    doc.save(out/'calibration',title='Calibration report — workflow 1.0')
    for ext in ('html','md'):(out/'calibration'/('validation_review.'+ext)).replace(out/'calibration'/('calibration_review.'+ext))
    index=Document();index.heading('Optimized calibration and validation performance')
    index.paragraph('These reports evaluate the specified completed audits against the original labelled populations. '
        'Stage 1 misses remain in the denominator; unlabelled input records remain excluded. '
        'Existing run outputs and dashboards are preserved. Reference/calibration overlaps and prior inspection limit independence claims.')
    fields=['dataset','TP','FP','TN','FN','PPV','recall','F1','specificity','family_TP','family_FP','family_FN','family_F1']
    index.table(fields,[[r[k] for k in fields] for r in metrics_rows])
    index.links([('Calibration report','calibration/calibration_review.html'),('Validation report','validation/validation_review.html'),
        ('Performance comparison','performance_comparison.tsv'),('Reporting provenance','reporting_provenance.json')])
    index.save(out,title='Optimized audit performance — workflow 1.0')
    for ext in ('html','md'):(out/('validation_review.'+ext)).replace(out/('report.'+ext))
    for path,digest in sources.items():
        if sha(Path(path))!=digest:raise ValueError('Reporting input changed: '+path)
    for kind in audits:
        if sha(audits[kind]/'workflow_manifest.json')!=sources[str(audits[kind]/'workflow_manifest.json')]:
            raise ValueError('Completed audit changed during reporting')
        core.verify(audits[kind],'audit',required_files=['Stage_2_output.xlsx','audited_family_assignments.tsv','negative_control_exclusions.tsv'])
    dump(out/'reporting_provenance.json',dict(version='1.0',created_at_utc=datetime.now(timezone.utc).isoformat(),
        original_run=str(run),audits={k:str(v) for k,v in audits.items()},
        source_hashes=sources,reference_labels={k:str(v[1]) for k,v in cohorts.items()},
        reference_label_checks=label_checks,
        checks=dict(exact_label_manifest_join=True,reference_sequences_verified=True,completed_audits_verified=True,
                    same_branch_annotation_and_handoff=True,sources_unchanged=True),
        analytical_pipeline_rerun=False,automatic_recalibration=False,backups_created=False))
    return dict(report=str(out/'report.html'),calibration_report=str(out/'calibration/calibration_review.html'),validation_report=review['report'])


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run',type=Path,required=True,help='Original isolated run containing both cohort manifests and handoffs')
    for kind in ('calibration','validation'):
        parser.add_argument('--'+kind+'-audit',type=Path,required=True,help='Completed standalone audit output')
        parser.add_argument('--'+kind+'-labels',type=Path,help='Optional matching label file; default is the original run contract')
    parser.add_argument('--out',type=Path,required=True,help='New report directory, separate from audits and original run')
    parser.add_argument('--allow-equivalent-calibration-label-bytes',action='store_true',
        help='Allow a repackaged calibration workbook only when all labels, cohort IDs, sequences and evaluation fields still match the frozen run exactly; validation label checks remain strict')
    try:print(json.dumps(generate(parser.parse_args(argv)),indent=2));return 0
    except (OSError,ValueError,KeyError) as exc:print('ERROR: '+str(exc),file=sys.stderr);return 1


if __name__=='__main__':raise SystemExit(main())
