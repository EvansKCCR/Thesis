"""Independent branch scheduling and sealed, dependency-aware execution."""

# Resolve source resources independently of the caller working directory.
import sys as _layout_sys
from pathlib import Path as _LayoutPath
_LAYOUT_BUNDLE_ROOT = next(p for p in _LayoutPath(__file__).resolve().parents
                          if (p/'config/layout_paths.json').is_file() and (p/'pipeline_paths.py').is_file())
if str(_LAYOUT_BUNDLE_ROOT) not in _layout_sys.path:
    _layout_sys.path.insert(0, str(_LAYOUT_BUNDLE_ROOT))
from pipeline_paths import activate as _layout_activate, resource_path as _layout_resource
_layout_activate()


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
from input_dependencies import InputRegistry


def algorithm_inventory():
    """Only analytical code/libraries, excluding datasets and generated results."""
    paths = set()
    for name in ('main', 'maintenance', 'resource_library', 'config'):
        for p in (ROOT/name).rglob('*'):
            if not p.is_file() or {'__pycache__','results','runs','tests'} & set(p.relative_to(ROOT).parts):
                continue
            if name == 'maintenance' and p.suffix != '.py':
                continue  # Maintenance receipts are history, not algorithm inputs.
            if p.suffix in ('.py','.json','.csv','.hmm','.faa','.yml','.txt','.tsv'):
                paths.add(p)
    paths.add(ROOT/'validation/identity_probability_model.json')
    return {p.relative_to(ROOT).as_posix():sha(p) for p in sorted(paths) if p.exists()}


def verify_frozen(contract):
    from pipeline_paths import recorded_digest_matches
    for name, digest in contract['algorithm_files'].items():
        path = _layout_resource(ROOT/name)
        exact = path.is_file() and sha(path) == digest
        relocated = ROOT == _LAYOUT_BUNDLE_ROOT and recorded_digest_matches(name, digest)
        if not exact and not relocated:
            raise ValueError('Frozen algorithm/library changed: '+name)
    raw=_layout_resource(Path(contract['references']['root']))
    selected=[raw/name for name in contract['references']['selection']] if 'selection' in contract['references'] else None
    reference = reference_inventory(raw,selected)
    if reference != contract['references']:
        raise ValueError('Frozen reference composition changed; use a new run')


def _verified_handoff(folder):
    """Require an intact detection checkpoint before adopting newer code."""
    manifest_path=_layout_resource(folder/'workflow_manifest.json')
    if not _layout_resource(manifest_path).is_file():return False
    try:manifest=json.loads(_layout_resource(manifest_path).read_text(encoding='utf-8'))
    except (OSError,ValueError):return False
    if manifest.get('kind')!='handoff':return False
    for name,expected in manifest.get('files',{}).items():
        path=_layout_resource(folder/name)
        if not _layout_resource(path).is_file() or sha(path)!=expected:return False
    return True


def _verified_folder(folder,kind):
    """Verify a sealed branch artifact before adopting an audit-only code fix."""
    manifest_path=_layout_resource(folder/'workflow_manifest.json')
    if not _layout_resource(manifest_path).is_file():return False
    try:manifest=json.loads(_layout_resource(manifest_path).read_text(encoding='utf-8'))
    except (OSError,ValueError):return False
    if manifest.get('kind')!=kind:return False
    for name,expected in manifest.get('files',{}).items():
        path=_layout_resource(folder/name)
        if not _layout_resource(path).is_file() or sha(path)!=expected:return False
    return True


def _awaiting_context_checkpoint(out,mode):
    """Only sealed detection-only gate runs can adopt changed implementation code."""
    summary_path=_layout_resource(out/'run_summary.json')
    if not _layout_resource(summary_path).is_file():return False
    try:summary=json.loads(_layout_resource(summary_path).read_text(encoding='utf-8'))
    except (OSError,ValueError):return False
    if summary.get('status')!='awaiting_external_evidence':return False
    kinds=['calibration',mode] if mode in ('validation','target') else []
    if not kinds or any(kind not in summary.get('gates',{}) for kind in kinds):return False
    for kind in kinds:
        gate=summary['gates'][kind]
        if gate.get('status')!='awaiting_external_evidence':return False
        if gate.get('stages',{}).get('detection',{}).get('status')!='complete':return False
        branch=_layout_resource(out/kind)
        status_path=_layout_resource(branch/'branch_status.json')
        if not _layout_resource(status_path).is_file():return False
        try:status=json.loads(_layout_resource(status_path).read_text(encoding='utf-8'))
        except (OSError,ValueError):return False
        if status.get('status')!='awaiting_external_evidence' or not _verified_handoff(_layout_resource(branch/'handoff')):
            return False
    return True


def _audit_join_resume_checkpoint(out,mode):
    """Permit an audit-join-only fix after inference, before any audit is sealed."""
    kinds=list(dict.fromkeys(['calibration',mode])) if mode in ('validation','target') else []
    if not kinds:return False
    for kind in kinds:
        branch=_layout_resource(out/kind)
        handoff=_layout_resource(branch/'handoff')
        annotation=_layout_resource(branch/'annotation')
        downstream=_layout_resource(branch/'downstream')
        if not _verified_handoff(handoff) or not _verified_folder(annotation,'annotation'):
            return False
        checkpoint=_layout_resource(downstream/'recovery.json')
        if not _layout_resource(checkpoint).is_file():return False
        try:state=json.loads(_layout_resource(checkpoint).read_text(encoding='utf-8'))
        except (OSError,ValueError):return False
        if state.get('annotation_sha256')!=sha(_layout_resource(annotation/'workflow_manifest.json')):
            return False
        stages=state.get('stages',{})
        if stages.get('audit',{}).get('status')=='complete':return False
        if _layout_resource(downstream/'workflow_manifest.json').exists() or _layout_resource(downstream/'final/workflow_manifest.json').exists():
            return False
        for stage_name in ('references','panel'):
            receipt=stages.get(stage_name,{})
            if receipt.get('status')!='complete':return False
            for name,expected in receipt.get('files',{}).items():
                path=_layout_resource(downstream/name)
                if not _layout_resource(path).is_file() or digest(path)!=expected:return False
        orthology=stages.get('orthology',{})
        if orthology.get('status')=='complete':
            for name,expected in orthology.get('files',{}).items():
                path=_layout_resource(downstream/name)
                if not _layout_resource(path).is_file() or digest(path)!=expected:return False
    return True


def reconcile_resume_contract(out,contract):
    """Resume unchanged inputs; adopt code-only updates at sealed handoff checkpoints."""
    out=_layout_resource(Path(out)).resolve();path=_layout_resource(out/'run_contract.json')
    if not _layout_resource(path).is_file():return contract
    old=json.loads(_layout_resource(path).read_text(encoding='utf-8'))
    if old==contract:return contract
    stable=lambda value:{k:v for k,v in value.items() if k not in ('algorithm_files','contract_sha256')}
    if stable(old)!=stable(contract):
        changed=sorted(k for k in set(old)|set(contract) if k not in ('algorithm_files','contract_sha256') and old.get(k)!=contract.get(k))
        raise ValueError('Input/configuration/code/reference contract changed outside the algorithm inventory: '+', '.join(changed)+'. This run cannot be relabelled or resumed.')
    old_files=old.get('algorithm_files',{});new_files=contract.get('algorithm_files',{})
    changed_files=sorted(name for name in set(old_files)|set(new_files) if old_files.get(name)!=new_files.get(name))
    if not changed_files:return contract
    unsafe=[name for name in changed_files if not (name.endswith(('.py','.yml','.yaml')) or name=='config/layout_paths.json')]
    if unsafe:
        raise ValueError('Frozen scientific rule/model resources changed: '+', '.join(unsafe)+'. Keep the original resources or start a new output directory.')
    audit_join_only=(set(changed_files)<= {'resource_library/scripts/audit_join.py',
                         'resource_library/schisto_adhesome/workflows.py'}
                     and 'resource_library/scripts/audit_join.py' in changed_files
                     and _audit_join_resume_checkpoint(out,contract.get('mode')))
    if not _awaiting_context_checkpoint(out,contract.get('mode')) and not audit_join_only:
        raise ValueError('Workflow code changed, and this run has no sealed detection-only external-evidence checkpoint. Keep the original code or start a new output directory.')
    log_path=_layout_resource(out/'run_contract_adoptions.json')
    history=json.loads(_layout_resource(log_path).read_text(encoding='utf-8')) if _layout_resource(log_path).is_file() else dict(adoptions=[])
    history.setdefault('adoptions',[]).append(dict(from_contract_sha256=old.get('contract_sha256'),
        to_contract_sha256=contract.get('contract_sha256'),changed_algorithm_files=changed_files,
        completed_gate_checkpoints_verified=True,
        adoption_scope='audit_join_only_before_audit' if audit_join_only else 'sealed_detection_handoff',
        adopted_utc=datetime.now(timezone.utc).isoformat()))
    dump(log_path,history);dump(path,contract)
    print('RESUME CONTRACT: adopted code-only update; sealed calibration/selected-gate detection handoffs verified; '
          f'{len(changed_files)} changed algorithm files recorded in {log_path}',flush=True)
    return contract


def prepare(args):
    settings = json.loads(_layout_resource(args.config).read_text())
    allowed = {'version','family_coverage','orthogroup_context_weight','qc','orthofinder_search',
               'orthofinder_method','orthofinder_inflation','benchmark_scope','benchmark_selection','stage2_profile','phylogenetic_identity_resolution','reference_files','reference_panel'}
    if set(settings)-allowed or settings.get('version') != '1.0':
        raise ValueError('Unknown configuration fields or configuration version != 1.0')
    settings['phylogenetic_identity_resolution']=bool(settings.get('phylogenetic_identity_resolution',False) or getattr(args,'phylogenetic_identity_resolution',False))
    for key in ('reference_files','reference_panel'):
        if getattr(args,key,None):settings[key]=getattr(args,key)
    from stage2_evidence import configuration as stage2_configuration
    settings.setdefault('stage2_profile',stage2_configuration()['default_profile'])
    if settings['stage2_profile'] not in stage2_configuration()['profiles']:
        raise ValueError('Unknown frozen Stage 2 profile')
    from orthogroup_context import validate_weight
    if validate_weight(settings.get('orthogroup_context_weight')) is None:
        raise ValueError('Isolated routine runs require a fixed context coefficient; use 1.0 for neutral expanded-baseline context')
    if settings.get('family_coverage') not in ('expanded','original') or 'qc' not in settings:
        raise ValueError('Select paired family coverage and explicit sequence QC')
    registry=InputRegistry(ROOT)
    decisions=[registry.assess('calibration_'+kind,'calibration.inputs',mandatory=False,
                              override=getattr(args,'calibration_'+kind),execution=True) for kind in ('fasta','labels')]
    calibration=[]
    datasets={}
    if all(event['available'] for event in decisions):
        calibration,cal_summary=load_dataset('calibration',args.calibration_fasta,args.calibration_labels)
        datasets['calibration']=(calibration,cal_summary)
    if args.mode == 'validation':
        for kind in ('fasta','labels'):
            registry.require('validation_'+kind,'validation.inputs',getattr(args,'validation_'+kind))
        validation,val_summary = load_dataset('validation',args.validation_fasta,args.validation_labels,
                                               select_labelled=args.select_labelled_validation)
        datasets['validation'] = (validation,val_summary)
    else:
        if not args.targets_config or not args.input_dir:
            raise ValueError('Target mode requires --input-dir and --targets-config')
        registry.require('target_configuration','target.configuration',args.targets_config)
        config = json.loads(_layout_resource(args.targets_config).read_text())
        if config.get('qc') != settings['qc']:
            raise ValueError('Target sequence QC differs from frozen branch configuration')
        registry.require('target_raw','target.inputs',args.input_dir)
        target_root=Path(args.input_dir).resolve()
        controls_paths=[getattr(args,name,None) or registry.path_for(name) for name in
                        ('calibration_fasta','calibration_labels','validation_fasta','validation_labels')]
        controls_paths += [registry.path_for(name) for name in ('negative_exclusion_fasta','negative_exclusion_labels')]
        for control_path in controls_paths:
            if target_root in Path(control_path).resolve().parents:
                raise ValueError('Evaluation controls must be outside the target-analysis input root')
        if Path(args.reference_raw).resolve()==target_root:
            raise ValueError('Reference and target input roots must be separate')
        sys.path.insert(0,str(_layout_resource(ROOT/'detection')))
        try:
            from adhesome_pipeline import check_inputs
            sources,ignored,summary = check_inputs(args.input_dir,config)
        finally:
            sys.path.pop(0)
        reserved = re.compile(r'^(?:Cal_|Val_|validation$|calibration$)',re.I)
        controls = {r['original_id'] for r in calibration}
        # Inspect the validation IDs solely for a production leakage gate, never run it.
        for label_path in (args.calibration_labels,args.validation_labels):
            if not _layout_resource(label_path).exists():continue
            from .datasets import label_rows
            controls |= {str(r.get('fasta_id') or r.get('accession')) for r in label_rows(label_path) if r.get('fasta_id') or r.get('accession')}
        for code, path in sources.items():
            if reserved.search(code) or controls & set(fasta(path)):
                raise ValueError('Evaluation control input cannot enter production target panel: '+code)
        target_manifest = dict(kind='target',targets_config=str(_layout_resource(args.targets_config).resolve()),
            targets_config_sha256=sha(args.targets_config),input_root=str(_layout_resource(args.input_dir).resolve()),
            files={str(p):sha(p) for p in sources.values()}, species=summary,ignored_fastas=ignored)
        datasets['target'] = (config,target_manifest)
    control_ids={r['original_id'] for r in calibration}
    if not calibration and Path(args.calibration_labels).is_file():
        from .datasets import label_rows
        control_ids.update(str(r.get('fasta_id') or r.get('accession')) for r in label_rows(args.calibration_labels) if r.get('fasta_id') or r.get('accession'))
    if 'validation' in datasets:control_ids.update(r['original_id'] for r in datasets['validation'][0])
    registry.require('reference_raw','orthology.references',args.reference_raw)
    selected=registry.select_references(panel=settings.get('reference_panel'),files=settings.get('reference_files'),raw=args.reference_raw)
    for path in selected if selected is not None else args.reference_raw.iterdir():
        if _layout_resource(path).is_file() and path.suffix.lower() in ('.faa','.fa','.fasta') and control_ids & set(fasta(path)):
            raise ValueError('Evaluation query identifiers present in raw reference input: '+str(path))
    evolutionary_evidence={}
    if settings['phylogenetic_identity_resolution']:
        for kind in datasets:
            path=getattr(args,kind+'_phylogenetic_evidence',None)
            if path:evolutionary_evidence[kind]=dict(path=str(_layout_resource(path).resolve()),sha256=sha(path) if _layout_resource(path).is_file() else None)
    contract = dict(evolutionary_evidence=evolutionary_evidence,version='1.0', mode=args.mode, frozen_settings=settings,
        algorithm_files=algorithm_inventory(),references=reference_inventory(args.reference_raw,selected),
        datasets={kind:summary for kind,(_,summary) in datasets.items()},
        automatic_recalibration=False, orthology_isolation='independent inference per branch',
        shared_reference_queries_modify_production=False,
        input_decisions=decisions, skipped_gates={} if 'calibration' in datasets else {
            'calibration':dict(status='skipped',reason='Optional calibration FASTA/labels unavailable',input_decisions=decisions)})
    encoded = json.dumps(contract,sort_keys=True).encode()
    contract['contract_sha256'] = hashlib.sha256(encoded).hexdigest()
    out = _layout_resource(args.out).resolve()
    if (_layout_resource(out/'run_contract.json')).exists():
        if not args.resume:
            raise ValueError('Run exists; use --resume with unchanged frozen inputs/configuration')
        contract=reconcile_resume_contract(out,contract)
    else:
        if _layout_resource(out).exists() and any(out.iterdir()):
            raise ValueError('Output directory is not empty')
        out.mkdir(parents=True,exist_ok=True)
        dump(_layout_resource(out/'run_contract.json'),contract)
    preflight = dict(status='passed',datasets=contract['datasets'],overlap_warning=None,
        no_length_filter_added=True,reference_version=contract['references']['version'])
    for kind,(data,summary) in datasets.items():
        branch = out/kind; branch.mkdir(exist_ok=True)
        dump(_layout_resource(branch/'input_manifest.json'),summary)
        if kind != 'target':
            materialize(kind,data,_layout_resource(branch/'inputs'),settings['qc'])
            table(_layout_resource(branch/'excluded_unlabelled.tsv'),[dict(original_id=i,reason='not_in_label_table') for i in summary['excluded_unlabelled_ids']])
            matches=reference_matches(data,args.reference_raw,selected)
            table(_layout_resource(branch/'reference_overlap.tsv'),matches)
            preflight.setdefault('reference_overlap',{})[kind]=len({r['original_id'] for r in matches})
        else:
            dump(_layout_resource(branch/'targets.json'),data)
    if args.mode == 'validation' and 'calibration' in datasets:
        overlap = overlaps(calibration,datasets['validation'][0])
        table(_layout_resource(out/'calibration_validation_overlap.tsv'),overlap)
        preflight['overlap_warning'] = dict(exact_sequence_overlaps=len(overlap),
            interpretation='Retain all labelled records; report overlap-free strata separately. Do not claim all records are independent.')
    dump(_layout_resource(out/'preflight.json'),preflight)
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
    contract = prepare(args); out=_layout_resource(args.out).resolve()
    if args.preflight_only:
        print('Preflight passed: '+str(_layout_resource(out/'preflight.json'))); return 0
    kinds = [kind for kind in ('calibration',args.mode) if kind in contract['datasets']]
    workers = 2 if len(kinds)>1 and args.schedule=='parallel' and args.threads>=2 else 1
    per_branch = max(1,args.threads//workers)
    started = datetime.now(timezone.utc).isoformat()
    def worker(kind):
        cmd = [sys.executable,'-B','-m','schisto_adhesome.branch','--kind',kind,'--run',str(out),
               '--threads',str(per_branch)]
        context = getattr(args,kind+'_context',None)
        if context:cmd += ['--context-dir',str(_layout_resource(context).resolve())]
        if args.allow_missing_context:cmd += ['--allow-missing-context']
        if kind=='calibration' and args.compare_calibration:
            cmd += ['--compare-report',str(_layout_resource(args.compare_calibration).resolve())]
        if kind=='target' and args.reuse_target_detection:
            cmd += ['--reuse-detection',str(_layout_resource(args.reuse_target_detection).resolve())]
        log = _layout_resource(out/kind/'branch.log')
        branch_started=datetime.now(timezone.utc)
        print(f'{branch_started.isoformat()} | BRANCH {kind.upper()} START; log: {log}',flush=True)
        with _layout_resource(log).open('a',encoding='utf-8') as handle:
            env=os.environ.copy()
            env['PYTHONPATH']=os.pathsep.join([str(ROOT/'resource_library'),str(ROOT),str(ROOT/'resource_library/scripts'),env.get('PYTHONPATH','')])
            proc=subprocess.Popen(cmd,cwd=ROOT,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,env=env,
                                  text=True,encoding='utf-8',errors='replace',bufsize=1)
            assert proc.stdout is not None
            for line in proc.stdout:
                handle.write(line);handle.flush()
                print(f'[{kind}] {line.rstrip()}',flush=True)
            proc.wait()
        elapsed=(datetime.now(timezone.utc)-branch_started).total_seconds()
        status_file = _layout_resource(out/kind/'branch_status.json')
        status = json.loads(_layout_resource(status_file).read_text()) if _layout_resource(status_file).exists() else dict(status='failed',exit_code=proc.returncode)
        if proc.returncode not in (0,2):status=dict(status='failed',exit_code=proc.returncode,error=status.get('error','Branch failed before completion'))
        status['log'] = str(log)
        status['elapsed_seconds']=round(elapsed,3)
        status['completed_utc']=datetime.now(timezone.utc).isoformat()
        print(f'{datetime.now(timezone.utc).isoformat()} | BRANCH {kind.upper()} {status.get("status", "finished").upper()}; elapsed {elapsed:.1f}s',flush=True)
        return status
    gates = launch_branches(kinds,worker,workers)
    if args.mode=='validation' and gates.get('validation',{}).get('status')=='complete':
        # Calibration can finish after validation in parallel mode. Refresh only
        # the terminal report after both workers return; never rerun either gate.
        from .validation_review import generate
        gate=gates['validation']
        try:
            review=generate(out)
            gate.update(report=review['report'],validation_review=review)
        except (OSError,ValueError,KeyError) as exc:
            gate.update(status='failed',error='Validation reporting failed: '+str(exc))
        dump(_layout_resource(out/'validation/branch_status.json'),gate)
    verify_frozen(contract)
    status = 'complete' if all(g.get('status')=='complete' for g in gates.values()) else (
        'failed' if any(g.get('status')=='failed' for g in gates.values()) else 'awaiting_external_evidence')
    summary = dict(version='1.0',selected_mode=args.mode,status=status,timestamp=started,
        completed_at=datetime.now(timezone.utc).isoformat(),scheduler=dict(workers=workers,threads_per_branch=per_branch),
        frozen_contract_sha256=contract['contract_sha256'],weights=contract['frozen_settings'].get('orthogroup_context_weight'),
        automatic_recalibration=False,gates=gates,skipped_gates=contract.get('skipped_gates',{}))
    dump(_layout_resource(out/'run_summary.json'),summary)
    links = '\n'.join(f'- {kind}: {gate.get("status")} — [{kind} report]({kind}/reports/'+
        ('validation_review.html' if kind=='validation' and gate.get('validation_review') else 'report.html')+
        f'), [{kind} log]({kind}/branch.log)' for kind,gate in gates.items())
    links += ''.join(f"\n- {kind}: skipped — {gate['reason']}" for kind,gate in contract.get('skipped_gates',{}).items())
    (_layout_resource(out/'README.md')).write_text(f'# Workflow 1.0: {args.mode}\n\nStatus: {status}\n\n{links}\n\nCalibration is a frozen internal self-test. No fitting occurs.\n',encoding='utf-8')
    print(json.dumps(summary,indent=2))
    return 1 if status=='failed' else 2 if status!='complete' else 0
