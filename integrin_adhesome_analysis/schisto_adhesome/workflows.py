"""Independent branch scheduling and sealed, dependency-aware execution."""
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys

from .datasets import ROOT, sha, dump, table, load_dataset, overlaps, materialize, reference_inventory, reference_matches, fasta


def algorithm_inventory():
    """Only analytical code/libraries, excluding datasets and generated results."""
    paths = set(ROOT.glob('*.py')) | set((ROOT/'schisto_adhesome').glob('*.py'))
    for name in ('01_detection','03_functional_motif_annotation','04_orthology_phylogenetic_classification',
                 '05_conservative_adhesome_family_assignment_audit','templates'):
        for p in (ROOT/name).rglob('*'):
            if p.is_file() and not {'__pycache__','results','runs','tests'} & set(p.relative_to(ROOT).parts) and (
                    p.suffix in ('.py','.json','.csv','.hmm','.faa') or p.name.endswith('.hmm.metadata.json')):
                paths.add(p)
    paths.add(ROOT/'02_calibration_set/identity_probability_model.json')
    return {p.relative_to(ROOT).as_posix():sha(p) for p in sorted(paths) if p.exists()}


def verify_frozen(contract):
    for name, digest in contract['algorithm_files'].items():
        if sha(ROOT/name) != digest:
            raise ValueError('Frozen algorithm/library changed: '+name)
    reference = reference_inventory(Path(contract['references']['root']))
    if reference != contract['references']:
        raise ValueError('Frozen reference composition changed; use a new run')


def prepare(args):
    settings = json.loads(args.config.read_text())
    allowed = {'version','family_coverage','orthogroup_context_weight','qc','orthofinder_search',
               'orthofinder_method','orthofinder_inflation','benchmark_scope','benchmark_selection'}
    if set(settings)-allowed or settings.get('version') != '1.0':
        raise ValueError('Unknown configuration fields or configuration version != 1.0')
    from orthogroup_context import validate_weight
    if validate_weight(settings.get('orthogroup_context_weight')) is None:
        raise ValueError('Isolated routine runs require a fixed context coefficient; use 1.0 for neutral expanded-baseline context')
    if settings.get('family_coverage') not in ('expanded','original') or 'qc' not in settings:
        raise ValueError('Select paired family coverage and explicit sequence QC')
    calibration, cal_summary = load_dataset('calibration',args.calibration_fasta,args.calibration_labels)
    datasets = {'calibration':(calibration,cal_summary)}
    if args.mode == 'validation':
        validation,val_summary = load_dataset('validation',args.validation_fasta,args.validation_labels,
                                               select_labelled=args.select_labelled_validation)
        datasets['validation'] = (validation,val_summary)
    else:
        if not args.targets_config or not args.input_dir:
            raise ValueError('Target mode requires --input-dir and --targets-config')
        config = json.loads(args.targets_config.read_text())
        if config.get('qc') != settings['qc']:
            raise ValueError('Target sequence QC differs from frozen branch configuration')
        sys.path.insert(0,str(ROOT/'01_detection'))
        try:
            from adhesome_pipeline import check_inputs
            sources,ignored,summary = check_inputs(args.input_dir,config)
        finally:
            sys.path.pop(0)
        reserved = re.compile(r'^(?:Cal_|Val_|validation$|calibration$)',re.I)
        controls = {r['original_id'] for r in calibration}
        # Inspect the validation IDs solely for a production leakage gate, never run it.
        if args.validation_labels.exists():
            from .datasets import label_rows
            controls |= {str(r.get('fasta_id')) for r in label_rows(args.validation_labels)}
        for code, path in sources.items():
            if reserved.search(code) or controls & set(fasta(path)):
                raise ValueError('Evaluation control input cannot enter production target panel: '+code)
        target_manifest = dict(kind='target',targets_config=str(args.targets_config.resolve()),
            targets_config_sha256=sha(args.targets_config),input_root=str(args.input_dir.resolve()),
            files={str(p):sha(p) for p in sources.values()}, species=summary,ignored_fastas=ignored)
        datasets['target'] = (config,target_manifest)
    control_ids={r['original_id'] for r in calibration}
    if 'validation' in datasets:control_ids.update(r['original_id'] for r in datasets['validation'][0])
    for path in args.reference_raw.iterdir():
        if path.is_file() and path.suffix.lower() in ('.faa','.fa','.fasta') and control_ids & set(fasta(path)):
            raise ValueError('Evaluation query identifiers present in raw reference input: '+str(path))
    contract = dict(version='1.0', mode=args.mode, frozen_settings=settings,
        algorithm_files=algorithm_inventory(),references=reference_inventory(args.reference_raw),
        datasets={kind:summary for kind,(_,summary) in datasets.items()},
        automatic_recalibration=False, orthology_isolation='independent inference per branch',
        shared_reference_queries_modify_production=False)
    encoded = json.dumps(contract,sort_keys=True).encode()
    contract['contract_sha256'] = hashlib.sha256(encoded).hexdigest()
    out = args.out.resolve()
    if (out/'run_contract.json').exists():
        if not args.resume:
            raise ValueError('Run exists; use --resume with unchanged frozen inputs/configuration')
        old = json.loads((out/'run_contract.json').read_text())
        if contract != old:
            raise ValueError('Input/configuration/code/reference contract changed; this run cannot be relabelled or resumed')
    else:
        if out.exists() and any(out.iterdir()):
            raise ValueError('Output directory is not empty')
        out.mkdir(parents=True,exist_ok=True)
        dump(out/'run_contract.json',contract)
    preflight = dict(status='passed',datasets=contract['datasets'],overlap_warning=None,
        no_length_filter_added=True,reference_version=contract['references']['version'])
    for kind,(data,summary) in datasets.items():
        branch = out/kind; branch.mkdir(exist_ok=True)
        dump(branch/'input_manifest.json',summary)
        if kind != 'target':
            materialize(kind,data,branch/'inputs',settings['qc'])
            table(branch/'excluded_unlabelled.tsv',[dict(original_id=i,reason='not_in_label_table') for i in summary['excluded_unlabelled_ids']])
            matches=reference_matches(data,args.reference_raw)
            table(branch/'reference_overlap.tsv',matches)
            preflight.setdefault('reference_overlap',{})[kind]=len({r['original_id'] for r in matches})
        else:
            dump(branch/'targets.json',data)
    if args.mode == 'validation':
        overlap = overlaps(calibration,datasets['validation'][0])
        table(out/'calibration_validation_overlap.tsv',overlap)
        preflight['overlap_warning'] = dict(exact_sequence_overlaps=len(overlap),
            interpretation='Retain all labelled records; report overlap-free strata separately. Do not claim all records are independent.')
    dump(out/'preflight.json',preflight)
    return contract


def launch_branches(kinds, worker, workers):
    """Always schedule both gates; one failure must not cancel the other gate."""
    results = {}
    with ThreadPoolExecutor(max_workers=workers) as pool:
        pending = {pool.submit(worker,kind):kind for kind in kinds}
        for future in as_completed(pending):
            kind = pending[future]
            try:
                results[kind] = future.result()
            except Exception as exc:
                results[kind] = dict(status='failed',error=str(exc))
    return results


def run(args):
    contract = prepare(args); out=args.out.resolve()
    if args.preflight_only:
        print('Preflight passed: '+str(out/'preflight.json')); return 0
    workers = 2 if args.schedule=='parallel' and args.threads>=2 else 1
    per_branch = max(1,args.threads//workers)
    kinds = ['calibration',args.mode]
    started = datetime.now(timezone.utc).isoformat()
    def worker(kind):
        cmd = [sys.executable,'-B','-m','schisto_adhesome.branch','--kind',kind,'--run',str(out),
               '--threads',str(per_branch)]
        context = getattr(args,kind+'_context',None)
        if context:cmd += ['--context-dir',str(context.resolve())]
        if args.allow_missing_context:cmd += ['--allow-missing-context']
        if kind=='calibration' and args.compare_calibration:
            cmd += ['--compare-report',str(args.compare_calibration.resolve())]
        if kind=='target' and args.reuse_target_detection:
            cmd += ['--reuse-detection',str(args.reuse_target_detection.resolve())]
        log = out/kind/'branch.log'
        with log.open('a',encoding='utf-8') as handle:
            proc = subprocess.run(cmd,cwd=ROOT,stdout=handle,stderr=subprocess.STDOUT)
        status_file = out/kind/'branch_status.json'
        status = json.loads(status_file.read_text()) if status_file.exists() else dict(status='failed',exit_code=proc.returncode)
        if proc.returncode not in (0,2):status=dict(status='failed',exit_code=proc.returncode,error=status.get('error','Branch failed before completion'))
        status['log'] = str(log)
        return status
    gates = launch_branches(kinds,worker,workers)
    verify_frozen(contract)
    status = 'complete' if all(g.get('status')=='complete' for g in gates.values()) else (
        'failed' if any(g.get('status')=='failed' for g in gates.values()) else 'awaiting_external_evidence')
    summary = dict(version='1.0',selected_mode=args.mode,status=status,timestamp=started,
        completed_at=datetime.now(timezone.utc).isoformat(),scheduler=dict(workers=workers,threads_per_branch=per_branch),
        frozen_contract_sha256=contract['contract_sha256'],weights=contract['frozen_settings'].get('orthogroup_context_weight'),
        automatic_recalibration=False,gates=gates)
    dump(out/'run_summary.json',summary)
    links = '\n'.join(f'- {kind}: {gate.get("status")} — [{kind} report]({kind}/reports/report.html), [{kind} log]({kind}/branch.log)' for kind,gate in gates.items())
    (out/'README.md').write_text(f'# Workflow 1.0: {args.mode}\n\nStatus: {status}\n\n{links}\n\nCalibration is a frozen internal self-test. No fitting occurs.\n',encoding='utf-8')
    print(json.dumps(summary,indent=2))
    return 1 if status=='failed' else 2 if status!='complete' else 0
