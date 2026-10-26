"""Process-isolated adapter to the existing detection/annotation/orthology core."""

# Resolve the organized bundle independently of the caller working directory.
import sys as _bundle_sys
from pathlib import Path as _BundlePath
_BUNDLE_ROOT=next(p for p in _BundlePath(__file__).resolve().parents if (p/'pipeline_paths.py').is_file())
if str(_BUNDLE_ROOT) not in _bundle_sys.path:_bundle_sys.path.insert(0,str(_BUNDLE_ROOT))
from pipeline_paths import activate as _activate_bundle_paths
_activate_bundle_paths()
import argparse
import json
import os
from pathlib import Path
import sys

from .datasets import ROOT, dump, sha, table
from .workflows import verify_frozen
from .reports import read, evaluate_branch, summarize_predictions, target_report
from .validation_review import generate as generate_validation_review
import adhesome_workflow as core


def reuse_detection(folder, config):
    """Import completed selected targets only after checking recorded inputs/resources."""
    from evaluate_family_expansion import local_path
    for entry in config['proteomes']:
        code=entry['code']; base=folder/code
        prior=json.loads((base/'prepared/preparation.json').read_text())['species'][code]
        raw=Path(entry['path'])
        if not raw.is_absolute():
            raise ValueError('Reuse requires absolute target FASTA paths in --targets-config')
        if sha(raw)!=prior['input_sha256'] or sha(base/'prepared'/f'{code}.faa')!=prior['output_sha256']:
            raise ValueError('Target sequence differs from detection receipt: '+code)
        classification=json.loads((base/'candidates/classification.json').read_text())
        for name,digest in classification['checksums'].items():
            path=local_path(name)
            if sha(path)!=digest:
                raise ValueError('Detection input changed: '+str(path))
        for path in (ROOT/'detection/config/family_rules.json',
                     ROOT/'detection/config/adhesome_pfam_manifest.csv',
                     ROOT/'detection/config/adhesome_component_crosswalk.csv'):
            matches=[v for k,v in classification['checksums'].items() if Path(k).name==path.name]
            if matches!=[sha(path)]:
                raise ValueError('Detection configuration is not current: '+str(path))
        scan=json.loads((base/'domains/scan.json').read_text())
        if scan['metadata_sha256']!=sha(ROOT/'detection/databases/pfam/adhesome_v1.0.hmm.metadata.json'):
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
    run=args.run.resolve(); branch=run/args.kind
    contract=json.loads((run/'run_contract.json').read_text()); settings=contract['frozen_settings']
    verify_frozen(contract)
    # These flags prohibit the legacy helper from archiving partial outputs and
    # prevent its old single-target evaluator from scoring a branch against itself.
    status_file=branch/'branch_status.json'
    stages={}; state_path=branch/'stages.json'
    if state_path.exists():stages=json.loads(state_path.read_text())
    def stage(name,folder,action,kind):
        receipt=folder/'workflow_manifest.json'
        if receipt.exists():
            core.verify(folder,kind)
            stages[name]=dict(status='complete',execution='reused',manifest_sha256=sha(receipt))
        else:
            stages[name]=dict(status='running',execution='recomputed');dump(state_path,stages)
            action()
            core.verify(folder,kind)
            stages[name]=dict(status='complete',execution='recomputed',manifest_sha256=sha(receipt))
        dump(state_path,stages)
    try:
        if args.kind=='target':
            inputs=Path(contract['datasets']['target']['input_root']); config=branch/'targets.json'
        else:
            inputs=branch/'inputs';config=inputs/'targets.json'
        if args.reuse_detection:
            reuse_detection(args.reuse_detection,json.loads(config.read_text()))
        cmd=['start','--input-dir',str(inputs),'--targets-config',str(config),
             '--out',str(branch/'handoff'),'--threads',str(args.threads)]
        if args.reuse_detection:cmd+=['--detection',str(args.reuse_detection)]
        elif (branch/'handoff/detection/batch.json').exists():cmd+=['--resume']
        def invoke(command):
            if core.main(command):raise ValueError('Core stage failed: '+' '.join(command[:2]))
        stage('detection',branch/'handoff',lambda:invoke(cmd),'handoff')
        # A first run remains useful without externally returned topology/CDD.
        if not args.context_dir and not args.allow_missing_context and not (branch/'annotation/workflow_manifest.json').exists():
            result=dict(status='awaiting_external_evidence',
                handoff=str(branch/'handoff'),candidate_export=str(branch/'handoff/external_candidates.faa'),
                hypotheses=str(branch/'handoff/hypotheses.tsv'),stages=stages)
            dump(status_file,result);return 2
        cmd=['resume','--handoff',str(branch/'handoff'),'--out',str(branch/'annotation'),
             '--annotation-only','--threads',str(args.threads)]
        if args.context_dir:cmd+=['--context-dir',str(args.context_dir)]
        if args.allow_missing_context:cmd+=['--allow-missing-context']
        if (branch/'annotation').exists() and not (branch/'annotation/workflow_manifest.json').exists():cmd+=['--retry-failed']
        stage('annotation',branch/'annotation',lambda:invoke(cmd),'annotation')
        downstream=branch/'downstream'
        cmd=['continue','--annotation',str(branch/'annotation'),'--out',str(downstream),
             '--reference-raw',contract['references']['root'],'--threads',str(args.threads),
             '--analysis-threads',str(args.threads),'--family-coverage',settings['family_coverage'],
             '--stage2-profile',settings.get('stage2_profile','original'),
             '--benchmark-scope',settings.get('benchmark_scope','all')]
        for key in ('orthogroup_context_weight','orthofinder_search','orthofinder_method',
                    'orthofinder_inflation','benchmark_selection'):
            if settings.get(key) is not None:cmd+=['--'+key.replace('_','-'),str(settings[key])]
        if (downstream/'recovery.json').exists():
            # Recovery restores the original frozen settings. Never adjust weights.
            cmd=['recover','--run',str(downstream),'--threads',str(args.threads),'--analysis-threads',str(args.threads)]
        stage('orthology_audit',downstream,lambda:invoke(cmd),'completed_workflow')
        audit=downstream/'final';reports=branch/'reports'
        if args.kind=='target':
            summary=target_report(audit,reports)
        else:
            manifest=read(branch/'inputs/sequence_manifest.tsv')
            cal=read(run/'calibration/inputs/sequence_manifest.tsv')
            references=set()
            references.update(r['sequence_sha256'] for r in read(branch/'reference_overlap.tsv'))
            for path in (ROOT/'detection/references/adhesome_v1.0/search_catalogue.csv',
                         ROOT/'detection/references/adhesome_v1.0/metazoan_catalogue.csv'):
                references.update(r['sequence_sha256'] for r in core.read(path,','))
            model=json.loads((ROOT/'validation/identity_probability_model.json').read_text())
            training={r['sequence_sha256'] for r in model.get('samples',[])}
            predictions=evaluate_branch(branch,manifest,audit,{r['sequence_sha256'] for r in cal},references,training)
            summary=summarize_predictions(predictions,reports,args.kind)
            if args.kind=='calibration':
                drift=[]
                if args.compare_report:
                    old={r['original_id']:r for r in read(args.compare_report/'predictions.tsv')}
                    for r in predictions:
                        prior=old.get(r['original_id'])
                        if prior and any(str(prior.get(k))!=str(r.get(k)) for k in ('predicted_family','grade','supported')):
                            drift.append(dict(original_id=r['original_id'],previous_family=prior.get('predicted_family'),
                                current_family=r['predicted_family'],previous_grade=prior.get('grade'),current_grade=r['grade']))
                    dump(reports/'drift_summary.json',dict(status='compared',changed_assignments=len(drift),
                        shared_controls=sum(r['original_id'] in old for r in predictions),
                        new_controls=sum(r['original_id'] not in old for r in predictions),
                        previous_report_sha256=sha(args.compare_report/'predictions.tsv')))
                else:dump(reports/'drift_summary.json',dict(status='no_previous_report_supplied'))
                table(reports/'assignment_drift.tsv',drift)
        review=generate_validation_review(run) if args.kind=='validation' else None
        verify_frozen(contract)
        result=dict(status='complete',report=review['report'] if review else str(reports/'report.html'),summary=summary,stages=stages,
                    reference_version=contract['references']['version'],automatic_recalibration=False)
        if review:result['validation_review']=review
        dump(status_file,result);return 0
    except Exception as exc:
        dump(status_file,dict(status='failed',error=str(exc),stages=stages));raise


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--kind',choices=['calibration','validation','target'],required=True)
    parser.add_argument('--run',type=Path,required=True);parser.add_argument('--threads',type=int,required=True)
    parser.add_argument('--context-dir',type=Path);parser.add_argument('--allow-missing-context',action='store_true')
    parser.add_argument('--reuse-detection',type=Path)
    parser.add_argument('--compare-report',type=Path)
    try:return execute(parser.parse_args())
    except Exception as exc:
        print('ERROR: '+str(exc),file=sys.stderr);return 1


if __name__=='__main__':sys.exit(main())
