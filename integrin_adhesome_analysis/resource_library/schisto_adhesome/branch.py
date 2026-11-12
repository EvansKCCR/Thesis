"""Process-isolated adapter to the existing detection/annotation/orthology core."""

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
import hashlib
import html
import json
import os
from pathlib import Path
import sys
import traceback
import time
from datetime import datetime, timezone

from .datasets import ROOT, dump, sha, table
from .workflows import verify_frozen
from .reports import read, evaluate_branch, summarize_predictions, target_report
from .validation_review import generate as generate_validation_review
import adhesome_workflow as core


def _validation_negative_panel(branch, contract, stages, threads):
    """Add only exact-source, unique RBH benchmark homologues to validation."""
    from .datasets import materialize, label_rows, reference_matches
    from .negative_homology import run as run_negative_homology, fasta_records
    branch=_layout_resource(branch);inputs=_layout_resource(branch/'inputs')
    output=_layout_resource(branch/'negative_homology');receipt=_layout_resource(output/'augmentation.json')
    manifest_path=_layout_resource(inputs/'sequence_manifest.tsv');config_path=_layout_resource(inputs/'targets.json')
    dataset=contract['datasets']['validation']
    source_hashes={str(_layout_resource(Path(dataset['protein_file']))):dataset['fasta_sha256'],
                   str(_layout_resource(Path(dataset['labels_file']))):dataset['labels_sha256']}
    if _layout_resource(receipt).is_file():
        saved=json.loads(_layout_resource(receipt).read_text(encoding='utf-8'))
        for source,expected in saved.get('source_hashes',{}).items():
            if not _layout_resource(Path(source)).is_file() or sha(_layout_resource(Path(source)))!=expected:
                raise ValueError('Original validation source changed after negative-homology augmentation: '+source)
        if saved.get('materialized_hashes',{}).get('sequence_manifest.tsv')!=sha(manifest_path) or saved.get('materialized_hashes',{}).get('targets.json')!=sha(config_path):
            raise ValueError('Saved validation negative-homology panel was prepared from different inputs')
        labels_checkpoint=_layout_resource(Path(saved.get('labels_path','')))
        if saved.get('labels_sha256') and (not _layout_resource(labels_checkpoint).is_file() or sha(labels_checkpoint)!=saved['labels_sha256']):
            raise ValueError('Saved validation label panel changed after negative-homology augmentation')
        for name,expected in saved.get('output_hashes',{}).items():
            path=_layout_resource(output/name)
            if not _layout_resource(path).is_file() or sha(path)!=expected:
                raise ValueError('Saved validation negative-homology evidence changed: '+str(path))
        stages['negative_homology']=dict(status=saved['summary']['status'],execution='reused',
            confirmed_one_to_one=saved['summary'].get('confirmed_one_to_one',0),manifest_sha256=sha(receipt))
        print(f'{datetime.now(timezone.utc).isoformat()} | STAGE NEGATIVE_HOMOLOGY REUSE; verified exact source and unique RBH checkpoint',flush=True)
        return saved
    output.mkdir(parents=True,exist_ok=False)
    print(f'{datetime.now(timezone.utc).isoformat()} | STAGE NEGATIVE_HOMOLOGY START',flush=True)
    entries=json.loads(_layout_resource(config_path).read_text(encoding='utf-8'))['proteomes']
    by_pid={}
    for entry in entries:
        path=_layout_resource(inputs/(entry['code']+'.faa'))
        for header,sequence in fasta_records(path):
            unit=header.split()[0];pid=entry['code']+'__'+unit
            by_pid[pid]=sequence
    original_rows=[];negative=[]
    for row in read(manifest_path):
        sequence=by_pid.get(row['protein_id'])
        if not sequence or hashlib.sha256(sequence.encode('ascii')).hexdigest()!=row.get('sequence_sha256'):
            raise ValueError('Validation protein manifest/FASTA sequence mismatch: '+row['protein_id'])
        record=dict(row,sequence=sequence)
        original_rows.append(record)
        if str(row.get('control_class','')).casefold() not in ('positive','positive_control','positivecontrol'):
            negative.append(record)
    reference_root=_layout_resource(Path(contract['references']['root']))
    qc=json.loads(_layout_resource(config_path).read_text(encoding='utf-8'))['qc']
    added,summary=run_negative_homology(original_rows,reference_root,output/'run',threads=threads,qc=qc)
    summary['base_validation_records']=len(original_rows)
    summary['base_negative_controls']=len(negative)
    summary['augmented_validation_records']=len(original_rows)+len(added)
    dump(_layout_resource(output/'homology_summary.json'),summary)
    original_label_path=_layout_resource(Path(contract['datasets']['validation']['labels_file']))
    if added:
        merged=original_rows+added
        materialize('validation',merged,inputs,qc)
        input_summary=json.loads(_layout_resource(branch/'input_manifest.json').read_text(encoding='utf-8'))
        input_summary['homology_added_negative_controls']=len(added)
        input_summary['labelled_records']=len(merged)
        counts={}
        for row in merged:counts[row['control_class']]=counts.get(row['control_class'],0)+1
        input_summary['counts']=counts
        input_summary['species']={}
        for row in merged:input_summary['species'][row['species']]=input_summary['species'].get(row['species'],0)+1
        dump(_layout_resource(branch/'input_manifest.json'),input_summary)
        selected=contract['references'].get('selection')
        selected_files=[reference_root/name for name in selected] if selected is not None else None
        overlap=reference_matches(merged,reference_root,selected_files)
        table(_layout_resource(branch/'reference_overlap.tsv'),overlap)
        # Keep original labels as supplied and add generated negative labels in a
        # run-local table; no source label workbook/CSV is edited.
        labels=label_rows(original_label_path)
        fields=list(dict.fromkeys(k for row in labels for k in row))
        additions=[dict(fasta_id=row['original_id'],**{'class':'negative_control'},sequence_length=row['length'],
            accession='',genes='',protein_name='') for row in added]
        for item in additions:
            for key in item:
                if key not in fields:fields.append(key)
        table(_layout_resource(inputs/'validation_labels_with_homologues.tsv'),labels+additions,fields)
    evidence=_layout_resource(output/'run/homology_matches.tsv')
    confirmed=_layout_resource(output/'run/confirmed_one_to_one_negative_homologues.tsv')
    evidence_path=evidence if _layout_resource(evidence).is_file() else confirmed if _layout_resource(confirmed).is_file() else None
    output_files=[]
    for path in _layout_resource(output).rglob('*'):
        if _layout_resource(path).is_file() and path.name!='augmentation.json':
            output_files.append((path.relative_to(output).as_posix(),sha(path)))
    # Include the final sequence inputs: their hashes prove the augmentation
    # was applied before Stage 1 detection and is stable during a resume.
    input_outputs={name:sha(_layout_resource(inputs/name)) for name in ('sequence_manifest.tsv','targets.json')}
    if evidence_path is not None:
        input_outputs['negative_homology_evidence']=sha(evidence_path)
    saved=dict(summary=summary,source_hashes=source_hashes,output_hashes=dict(output_files),
        materialized_hashes=input_outputs,evidence_path=str(evidence_path) if evidence_path else '',
        labels_path=str(_layout_resource(inputs/'validation_labels_with_homologues.tsv')) if added else
            str(original_label_path))
    saved['labels_sha256']=sha(_layout_resource(Path(saved['labels_path'])))
    dump(receipt,saved)
    stages['negative_homology']=dict(status=summary['status'],execution='recomputed',
        confirmed_one_to_one=len(added),manifest_sha256=sha(receipt))
    print(f'{datetime.now(timezone.utc).isoformat()} | STAGE NEGATIVE_HOMOLOGY COMPLETE; confirmed {len(added)} one-to-one homologues; unverified matches excluded',flush=True)
    return saved


def _posthoc_refinement(branch, audit, kind, inputs, stages):
    """Run proposal-only Pfam/CDD refinement and labeled reporting after Stage 3."""
    branch=_layout_resource(branch);audit=_layout_resource(audit)
    review=_layout_resource(branch/'posthoc_refinement')
    receipt=_layout_resource(review/'run_manifest.json')
    reused=False
    if _layout_resource(receipt).is_file():
        reused=True
        saved=json.loads(_layout_resource(receipt).read_text(encoding='utf-8'))
        if saved.get('status')!='complete' or saved.get('input_sha256')!=sha(_layout_resource(audit/'Stage_2_output.xlsx')):
            raise ValueError('Existing post hoc result is incomplete or belongs to a different Stage 3 audit')
        for name,expected in saved.get('output_sha256',{}).items():
            path=_layout_resource(review/name)
            if not _layout_resource(path).is_file() or sha(path)!=expected:
                raise ValueError('Saved post hoc output changed: '+str(path))
        print(f'{datetime.now(timezone.utc).isoformat()} | STAGE POSTHOC_REFINEMENT REUSE; source audit and outputs verified',flush=True)
    elif _layout_resource(review).exists():
        raise ValueError('Incomplete post hoc output exists; inspect and move it before resuming: '+str(review))
    else:
        print(f'{datetime.now(timezone.utc).isoformat()} | STAGE POSTHOC_REFINEMENT START',flush=True)
        core.refine_unresolved(argparse.Namespace(input=audit,sheet='Audited_family_assignment',kind=kind if kind in ('calibration','validation') else None,
            out=review,check_inputs=False,compare=None,strict=False))
        saved=json.loads(_layout_resource(receipt).read_text(encoding='utf-8'))
        if saved.get('status')!='complete' or not saved.get('source_record_accounting_passed'):
            raise ValueError('Post hoc refinement did not produce a complete accounting receipt')
        print(f'{datetime.now(timezone.utc).isoformat()} | STAGE POSTHOC_REFINEMENT COMPLETE',flush=True)
    stages['posthoc_refinement']=dict(status='complete',execution='reused' if reused else 'recomputed',
        rescued_or_reviewed=saved.get('unassigned_review_records',0),
        potential_false_positive_candidates=saved.get('potential_false_positive_candidates',0),
        source_sha256=saved.get('input_sha256'))
    if kind in ('calibration','validation'):
        contract=json.loads(_layout_resource(branch.parent/'run_contract.json').read_text(encoding='utf-8'))
        validation_checkpoint=_layout_resource(branch/'negative_homology/augmentation.json')
        label_path=Path(contract['datasets'][kind]['labels_file'])
        if kind=='validation' and _layout_resource(validation_checkpoint).is_file():
            homology=json.loads(_layout_resource(validation_checkpoint).read_text(encoding='utf-8'))
            label_path=Path(homology.get('labels_path') or label_path)
        report_dir=_layout_resource(branch/'reports/posthoc_performance')
        report_receipt=_layout_resource(report_dir/'report_manifest.json')
        if not _layout_resource(report_receipt).is_file():
            if _layout_resource(report_dir).exists():
                raise ValueError('Incomplete post hoc performance report exists: '+str(report_dir))
            core.report_refinement(argparse.Namespace(refinement=review,labels=label_path,
                sequence_manifest=_layout_resource(inputs/'sequence_manifest.tsv'),kind=kind,
                label_sheet='control_identity',aliases=None,orthogroups=None,panel_proteomes=None,
                reference_column=[],out=report_dir))
        else:
            report=json.loads(_layout_resource(report_receipt).read_text(encoding='utf-8'))
            if report.get('status')!='complete' or report.get('source_audit_sha256')!=saved.get('input_sha256'):
                raise ValueError('Saved post hoc performance report does not match the refinement source')
        stages['posthoc_performance_report']=dict(status='complete',directory=str(report_dir))
    return saved


def _stage_performance_report(branch, kind, stages):
    from .stage_performance_report import generate as generate_stage_performance
    branch=_layout_resource(branch);stage_report=_layout_resource(branch/'reports/stage_performance')
    predictions=_layout_resource(branch/'reports/predictions.tsv')
    posthoc_report=_layout_resource(branch/'reports/posthoc_performance')
    refinement=_layout_resource(branch/'posthoc_refinement')
    manifest=_layout_resource(stage_report/'summary.json')
    if _layout_resource(manifest).is_file():
        current=json.loads(_layout_resource(manifest).read_text(encoding='utf-8'))
        if (current.get('status')!='complete' or current.get('source_predictions_sha256')!=sha(predictions)
                or current.get('source_posthoc_predictions_sha256')!=sha(_layout_resource(posthoc_report/'predictions.tsv'))):
            raise ValueError('Saved stage performance report does not match current audit/post hoc predictions')
        for name,expected in current.get('output_sha256',{}).items():
            path=_layout_resource(stage_report/name)
            if not _layout_resource(path).is_file() or sha(path)!=expected:
                raise ValueError('Saved stage performance output changed: '+str(path))
    else:
        if _layout_resource(stage_report).exists():
            raise ValueError('Incomplete stage performance report exists: '+str(stage_report))
        summary=generate_stage_performance(predictions,posthoc_report,refinement,stage_report,kind)
    main_report=_layout_resource(branch/'reports/report.html')
    if _layout_resource(main_report).is_file():
        content=_layout_resource(main_report).read_text(encoding='utf-8')
        if 'stage_performance/stage_performance.html' not in content:
            summary=json.loads(_layout_resource(manifest).read_text(encoding='utf-8'))
            metric_rows=read(_layout_resource(stage_report/'stage_performance.tsv'))
            fields=('stage','n','TP','FP','TN','FN','precision','recall','F1')
            esc=lambda value:html.escape(str(value if value not in (None,'') else 'n.a.'))
            cell=lambda row,field:(f'{float(row[field]):.1%}' if field in ('precision','recall','F1') and row.get(field) not in ('',None) else row.get(field))
            table_html='<table><thead><tr>'+''.join('<th>'+esc(field)+'</th>' for field in fields)+'</tr></thead><tbody>'
            table_html+=''.join('<tr>'+''.join('<td>'+esc(cell(row,field))+'</td>' for field in fields)+'</tr>' for row in metric_rows)+'</tbody></table>'
            with _layout_resource(main_report).open('a',encoding='utf-8') as handle:
                handle.write('<h2>Stage-by-stage and post hoc performance</h2><p>All stages use the same labelled cohort (n='+str(summary.get('denominator','n.a.'))+'). Stage 2 is pre-exclusion family annotation; Stage 3 is the supported final audit. Post hoc proposals are shown for comparison and are not promoted. Rescued candidates: '+str(summary.get('posthoc_candidates_rescued',0))+'; potential false-positive flags: '+str(summary.get('supported_assignments_flagged_as_potential_false_positive',0))+'.</p>'+table_html+'<p>Incremental changes and candidate evidence: <a href="stage_performance/stage_performance.html">detailed stage report</a>.</p>\n')
    stages['stage_performance_report']=dict(status='complete',directory=str(stage_report))


def reuse_detection(folder, config):
    """Import completed selected targets only after checking recorded inputs/resources."""
    from evaluate_family_expansion import local_path
    for entry in config['proteomes']:
        code=entry['code']; base=folder/code
        prior=json.loads((_layout_resource(base/'prepared/preparation.json')).read_text())['species'][code]
        raw=_layout_resource(Path(entry['path']))
        if not raw.is_absolute():
            raise ValueError('Reuse requires absolute target FASTA paths in --targets-config')
        if sha(raw)!=prior['input_sha256'] or sha(_layout_resource(base/'prepared'/f'{code}.faa'))!=prior['output_sha256']:
            raise ValueError('Target sequence differs from detection receipt: '+code)
        classification=json.loads((_layout_resource(base/'candidates/classification.json')).read_text())
        for name,digest in classification['checksums'].items():
            path=local_path(name)
            if sha(path)!=digest:
                raise ValueError('Detection input changed: '+str(path))
        for path in (_layout_resource(ROOT/'detection/config/family_rules.json'),
                     _layout_resource(ROOT/'detection/config/adhesome_pfam_manifest.csv'),
                     _layout_resource(ROOT/'detection/config/adhesome_component_crosswalk.csv')):
            matches=[v for k,v in classification['checksums'].items() if _layout_resource(Path(k)).name==path.name]
            if matches!=[sha(path)]:
                raise ValueError('Detection configuration is not current: '+str(path))
        scan=json.loads((_layout_resource(base/'domains/scan.json')).read_text())
        if scan['metadata_sha256']!=sha(_layout_resource(ROOT/'detection/databases/pfam/adhesome_v1.0.hmm.metadata.json')):
            raise ValueError('Detection domain library changed')


def execute(args):
    flags={'SCHISTO_NO_ARCHIVES':'1','SCHISTO_ISOLATED_BRANCH':args.kind}
    previous={key:os.environ.get(key) for key in flags}
    try:
        os.environ.update(flags)
        return _execute(args)
    finally:
        for key,value in previous.items():
            if value is None:os.environ.pop(key,None)
            else:os.environ[key]=value


def _execute(args):
    run=_layout_resource(args.run).resolve(); branch=run/args.kind
    contract=json.loads((_layout_resource(run/'run_contract.json')).read_text()); settings=contract['frozen_settings']
    verify_frozen(contract)
    # These flags prohibit the legacy helper from archiving partial outputs and
    # prevent its old single-target evaluator from scoring a branch against itself.
    status_file=_layout_resource(branch/'branch_status.json')
    stages={}; state_path=_layout_resource(branch/'stages.json')
    if _layout_resource(state_path).exists():stages=json.loads(_layout_resource(state_path).read_text())
    def stage(name,folder,action,kind):
        receipt=_layout_resource(folder/'workflow_manifest.json')
        if _layout_resource(receipt).exists():
            core.verify(folder,kind)
            stages[name]=dict(status='complete',execution='reused',manifest_sha256=sha(receipt))
            print(f'{datetime.now(timezone.utc).isoformat()} | STAGE {name.upper()} REUSE; checkpoint verified',flush=True)
        else:
            started=time.perf_counter();started_utc=datetime.now(timezone.utc).isoformat()
            stages[name]=dict(status='running',execution='recomputed',started_utc=started_utc);dump(state_path,stages)
            print(f'{started_utc} | STAGE {name.upper()} START',flush=True)
            try:
                action()
                core.verify(folder,kind)
                elapsed=time.perf_counter()-started
                stages[name]=dict(status='complete',execution='recomputed',manifest_sha256=sha(receipt),
                    started_utc=started_utc,completed_utc=datetime.now(timezone.utc).isoformat(),elapsed_seconds=round(elapsed,3))
                print(f'{datetime.now(timezone.utc).isoformat()} | STAGE {name.upper()} COMPLETE; elapsed {elapsed:.1f}s',flush=True)
            except Exception as exc:
                elapsed=time.perf_counter()-started
                stages[name]=dict(status='failed',execution='recomputed',error=str(exc),started_utc=started_utc,
                    failed_utc=datetime.now(timezone.utc).isoformat(),elapsed_seconds=round(elapsed,3))
                dump(state_path,stages)
                print(f'{datetime.now(timezone.utc).isoformat()} | STAGE {name.upper()} FAILED; elapsed {elapsed:.1f}s; {exc}',flush=True)
                # Keep the originating file and line in branch.log; the short
                # status error alone is insufficient to diagnose NameError and
                # other exceptions from imported workflow modules.
                traceback.print_exc(file=sys.stdout)
                raise
        dump(state_path,stages)
    try:
        if args.kind=='target':
            inputs=_layout_resource(Path(contract['datasets']['target']['input_root'])); config=_layout_resource(branch/'targets.json')
        else:
            inputs=_layout_resource(branch/'inputs');config=_layout_resource(inputs/'targets.json')
        if args.kind=='validation':
            _validation_negative_panel(branch,contract,stages,args.threads)
            dump(state_path,stages)
        if args.reuse_detection:
            reuse_detection(args.reuse_detection,json.loads(_layout_resource(config).read_text()))
        cmd=['start','--input-dir',str(inputs),'--targets-config',str(config),
             '--out',str(_layout_resource(branch/'handoff')),'--threads',str(args.threads)]
        if args.reuse_detection:cmd+=['--detection',str(args.reuse_detection)]
        elif (_layout_resource(branch/'handoff/detection/batch.json')).exists():cmd+=['--resume']
        def invoke(command):
            if core.main(command):raise ValueError('Core stage failed: '+' '.join(command[:2]))
        stage('detection',_layout_resource(branch/'handoff'),lambda:invoke(cmd),'handoff')
        # A first run remains useful without externally returned topology/CDD.
        if not args.context_dir and not args.allow_missing_context and not (_layout_resource(branch/'annotation/workflow_manifest.json')).exists():
            result=dict(status='awaiting_external_evidence',
                handoff=str(_layout_resource(branch/'handoff')),candidate_export=str(_layout_resource(branch/'handoff/external_candidates.faa')),
                hypotheses=str(_layout_resource(branch/'handoff/hypotheses.tsv')),stages=stages)
            dump(status_file,result);return 2
        cmd=['resume','--handoff',str(_layout_resource(branch/'handoff')),'--out',str(_layout_resource(branch/'annotation')),
             '--annotation-only','--threads',str(args.threads)]
        if args.context_dir:cmd+=['--context-dir',str(args.context_dir)]
        if args.allow_missing_context:cmd+=['--allow-missing-context']
        if (_layout_resource(branch/'annotation')).exists() and not (_layout_resource(branch/'annotation/workflow_manifest.json')).exists():cmd+=['--retry-failed']
        stage('annotation',_layout_resource(branch/'annotation'),lambda:invoke(cmd),'annotation')
        downstream=_layout_resource(branch/'downstream')
        cmd=['continue','--annotation',str(_layout_resource(branch/'annotation')),'--out',str(downstream),
             '--reference-raw',contract['references']['root'],'--threads',str(args.threads),
             '--analysis-threads',str(args.threads),'--family-coverage',settings['family_coverage'],
             '--stage2-profile',settings.get('stage2_profile','original'),
             '--benchmark-scope',settings.get('benchmark_scope','all')]
        for filename in settings.get('reference_files',[]):cmd+=['--reference-file',filename]
        if settings.get('reference_panel'):cmd+=['--reference-panel',settings['reference_panel']]
        for key in ('orthogroup_context_weight','orthofinder_search','orthofinder_method',
                    'orthofinder_inflation','benchmark_selection'):
            if settings.get(key) is not None:cmd+=['--'+key.replace('_','-'),str(settings[key])]
        if args.kind=='target' and getattr(args,'phylogeny',False):
            cmd+=['--phylogeny','--phylogeny-grouping',getattr(args,'phylogeny_grouping','orthogroup')]
        if (_layout_resource(downstream/'recovery.json')).exists():
            # Recovery restores the original frozen settings. Never adjust weights.
            cmd=['recover','--run',str(downstream),'--threads',str(args.threads),'--analysis-threads',str(args.threads)]
            if args.kind=='target' and getattr(args,'phylogeny',False):
                cmd+=['--phylogeny','--phylogeny-grouping',getattr(args,'phylogeny_grouping','orthogroup')]
        if settings.get('phylogenetic_identity_resolution',False):
            cmd+=['--phylogenetic-identity-resolution']
            phylo=contract.get('evolutionary_evidence',{}).get(args.kind)
            if phylo:
                source=_layout_resource(Path(phylo['path']))
                if _layout_resource(source).is_file() and sha(source)==phylo['sha256']:
                    cmd+=['--phylogenetic-identity-evidence',phylo['path']]
                else:
                    dump(_layout_resource(branch/'optional_phylogenetic_evidence_warning.json'),dict(status='unavailable',reason='Frozen optional phylogenetic evidence missing or changed; family audit will complete without it'))
        stage('orthology_audit',downstream,lambda:invoke(cmd),'completed_workflow')
        audit=_layout_resource(downstream/'final');reports=_layout_resource(branch/'reports')
        _posthoc_refinement(branch,audit,args.kind,inputs,stages)
        dump(state_path,stages)
        if args.kind=='target':
            summary=target_report(audit,reports)
        else:
            manifest=read(_layout_resource(branch/'inputs/sequence_manifest.tsv'))
            cal_path=_layout_resource(run/'calibration/inputs/sequence_manifest.tsv')
            cal=read(cal_path) if cal_path.is_file() else None
            references=set()
            references.update(r['sequence_sha256'] for r in read(_layout_resource(branch/'reference_overlap.tsv')))
            for path in (_layout_resource(ROOT/'detection/references/adhesome_v1.0/search_catalogue.csv'),
                         _layout_resource(ROOT/'detection/references/adhesome_v1.0/metazoan_catalogue.csv')):
                references.update(r['sequence_sha256'] for r in core.read(path,','))
            model_path=_layout_resource(ROOT/'validation/identity_probability_model.json')
            model=json.loads(model_path.read_text()) if model_path.is_file() else None
            training={r['sequence_sha256'] for r in model.get('samples',[])} if model is not None else None
            predictions=evaluate_branch(branch,manifest,audit,{r['sequence_sha256'] for r in cal} if cal is not None else None,references,training)
            summary=summarize_predictions(predictions,reports,args.kind)
            _stage_performance_report(branch,args.kind,stages)
            dump(state_path,stages)
            if args.kind=='calibration':
                drift=[]
                if args.compare_report:
                    old={r['original_id']:r for r in read(_layout_resource(args.compare_report/'predictions.tsv'))}
                    for r in predictions:
                        prior=old.get(r['original_id'])
                        if prior and any(str(prior.get(k))!=str(r.get(k)) for k in ('predicted_family','grade','supported')):
                            drift.append(dict(original_id=r['original_id'],previous_family=prior.get('predicted_family'),
                                current_family=r['predicted_family'],previous_grade=prior.get('grade'),current_grade=r['grade']))
                    dump(_layout_resource(reports/'drift_summary.json'),dict(status='compared',changed_assignments=len(drift),
                        shared_controls=sum(r['original_id'] in old for r in predictions),
                        new_controls=sum(r['original_id'] not in old for r in predictions),
                        previous_report_sha256=sha(_layout_resource(args.compare_report/'predictions.tsv'))))
                else:dump(_layout_resource(reports/'drift_summary.json'),dict(status='no_previous_report_supplied'))
                table(_layout_resource(reports/'assignment_drift.tsv'),drift)
        review=generate_validation_review(run) if args.kind=='validation' else None
        verify_frozen(contract)
        result=dict(status='complete',report=review['report'] if review else str(_layout_resource(reports/'report.html')),summary=summary,stages=stages,
                    reference_version=contract['references']['version'],automatic_recalibration=False)
        if review:result['validation_review']=review
        dump(status_file,result);return 0
    except Exception as exc:
        dump(status_file,dict(status='failed',error=str(exc),stages=stages));raise


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--kind',choices=['calibration','validation','target'],required=True)
    parser.add_argument('--run',type=Path,required=True);parser.add_argument('--threads',type=int,required=True)
    parser.add_argument('--phylogeny',action='store_true',help='Target branch: request optional candidate phylogenies')
    parser.add_argument('--phylogeny-grouping',choices=['orthogroup','family'],default='orthogroup')
    parser.add_argument('--context-dir',type=Path);parser.add_argument('--allow-missing-context',action='store_true')
    parser.add_argument('--reuse-detection',type=Path)
    parser.add_argument('--compare-report',type=Path)
    try:return execute(parser.parse_args())
    except Exception as exc:
        print('ERROR: '+str(exc),file=sys.stderr);return 1


if __name__=='__main__':sys.exit(main())
