#!/usr/bin/env python3
"""Repair a returned-context checkpoint; run analysis ONLY with --resume.

This maintenance adapter lives outside the frozen analytical-code inventory.
It never launches the calibration branch or changes analytical rules/references.
--prepare queues a context refresh, preserving the old sealed annotation until
the user explicitly resumes. No backups or archived attempts are created.
"""

# Resolve the organized bundle independently of the caller working directory.
import sys as _bundle_sys
from pathlib import Path as _BundlePath
_BUNDLE_ROOT=next(p for p in _BundlePath(__file__).resolve().parents if (p/'pipeline_paths.py').is_file())
if str(_BUNDLE_ROOT) not in _bundle_sys.path:_bundle_sys.path.insert(0,str(_BUNDLE_ROOT))
from pipeline_paths import activate as _activate_bundle_paths
_activate_bundle_paths()
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

TOPOLOGY = 'deeptmhmm_predicted_topologies.3line'
PLAN = 'validation/context_refresh.json'


def local(value):
    value = str(value)
    if os.name == 'nt' and re.match(r'^/mnt/[a-z]/', value):
        return Path(value[5].upper() + ':/' + value[7:])
    return Path(value)


def sha(path):
    with Path(path).open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def load(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def dump(path, value):
    Path(path).write_text(json.dumps(value, indent=2) + '\n', encoding='utf-8')


def now():
    return datetime.now(timezone.utc).isoformat()


def check_files(base, files):
    for name, expected in files.items():
        path = base / name.replace('\\', '/')
        if not path.is_file() or sha(path) != expected:
            raise ValueError('Protected/checkpoint file changed: ' + str(path))


def verify_seal(folder, kind):
    data = load(folder / 'workflow_manifest.json')
    if data['kind'] != kind:
        raise ValueError('Unexpected sealed stage: ' + str(folder))
    check_files(folder, data['files'])
    return data


def fasta(path):
    records = {}; identifier = None; parts = []
    for line in Path(path).read_text(encoding='utf-8-sig').splitlines() + ['>']:
        if line.startswith('>'):
            if identifier is not None:
                if identifier in records:
                    raise ValueError('Duplicate FASTA identifier: ' + identifier)
                records[identifier] = ''.join(parts)
            identifier = line[1:].split()[0] if line[1:].strip() else None
            parts = []
        elif line.strip():
            parts.append(line.strip())
    return records


def topologies(path, sequences):
    lines = [line.strip() for line in Path(path).read_text(encoding='utf-8-sig').splitlines() if line.strip()]
    if len(lines) % 3:
        raise ValueError('Incomplete 3line topology record')
    records = {}
    for index in range(0, len(lines), 3):
        header, sequence, topology = lines[index:index + 3]
        identifier = header[1:].split()[0] if header.startswith('>') else ''
        if (identifier in records or sequences.get(identifier) != sequence or
                len(topology) != len(sequence) or set(topology) - set('IOMSBX')):
            raise ValueError('Duplicate, mismatched or invalid topology: ' + identifier)
        records[identifier] = topology
    return records


def intact_search(path):
    """Require a complete gzip stream and complete BLAST6 rows, including EOF."""
    if not path.is_file() or not path.stat().st_size:
        return False
    try:
        with gzip.open(path, 'rt', encoding='ascii') as handle:
            for line in handle:
                fields = line.rstrip('\n').split('\t')
                if not line.endswith('\n') or len(fields) != 12:
                    return False
                float(fields[11])
        return True
    except (OSError, EOFError, ValueError):
        return False


def bind_search_panel(wd, panel):
    species = {}
    for line in (wd / 'SpeciesIDs.txt').read_text().splitlines():
        key, name = line.split(': ', 1); species[int(key)] = name
    ids = dict(line.split(': ', 1) for line in (wd / 'SequenceIDs.txt').read_text().splitlines())
    for number, name in species.items():
        expected = fasta(panel / 'proteomes' / name)
        actual = {ids[k]: v for k, v in fasta(wd / f'Species{number}.fa').items()}
        if actual != expected:
            raise ValueError('Prepared search sequences differ from panel: ' + name)
    return species


def prepare(run, apply=False):
    branch = run / 'validation'; annotation = branch / 'annotation'
    if (run / PLAN).exists():
        raise ValueError('A refresh plan already exists; use --resume with that plan')
    if load(run / 'calibration/branch_status.json')['status'] != 'complete':
        raise ValueError('This adapter requires completed calibration')
    check_files(run.parents[1], load(run / 'run_contract.json')['algorithm_files'])
    verify_seal(branch / 'handoff', 'handoff')
    meta = verify_seal(annotation, 'annotation')
    recovery = load(branch / 'downstream/recovery.json')
    if (branch / 'downstream/workflow_manifest.json').exists() or (branch / 'downstream/final').exists():
        raise ValueError('Final audit already exists; review its dependent outputs separately')
    if recovery['annotation_sha256'] != sha(annotation / 'workflow_manifest.json'):
        raise ValueError('Downstream checkpoint uses another annotation')
    if meta['handoff_manifest'] != sha(branch / 'handoff/workflow_manifest.json'):
        raise ValueError('Annotation/handoff binding changed')
    sources = meta['sources']
    keys = [key for key in sources if key.replace('\\', '/').endswith('/' + TOPOLOGY)]
    if len(keys) != 1:
        raise ValueError('Expected one topology source')
    key = keys[0]; source = local(key)
    for name, expected in sources.items():
        if name != key and sha(local(name)) != expected:
            raise ValueError('Other external evidence changed: ' + name)
    sequences = fasta(branch / 'handoff/external_candidates.faa')
    old = topologies(annotation / 'context' / TOPOLOGY, sequences)
    updated = topologies(source, sequences)
    if set(updated) != set(sequences):
        raise ValueError('The topology return must cover every exported candidate')
    if any(updated.get(pid) != value for pid, value in old.items()):
        raise ValueError('Previously accepted topology changed; this adapter accepts additions only')
    added = sorted(set(updated) - set(old))
    if not added:
        raise ValueError('No newly returned topology records')
    protected = {}
    for folder in ('calibration',):
        # All calibration analytical seals plus its reports and branch receipts.
        for path in (run / folder).rglob('*'):
            if path.is_file() and (path.name == 'workflow_manifest.json' or 'reports' in path.parts
                                   or path.name in ('branch_status.json', 'stages.json')):
                protected[path.relative_to(run).as_posix()] = sha(path)
    for name in ('run_contract.json', 'validation/handoff/workflow_manifest.json',
                 'validation/handoff/external_candidates.faa'):
        protected[name] = sha(run / name)
    for name in ('references', 'panel'):
        stage = recovery['stages'][name]
        if stage['status'] != 'complete':
            raise ValueError('Reference/panel preparation must be complete')
        check_files(branch / 'downstream', stage['files'])
        protected.update({'validation/downstream/' + k.replace('\\', '/'): v for k, v in stage['files'].items()})
    orth = branch / 'downstream/orthology'
    searches = {}; interrupted = []
    if recovery['stages'].get('orthology', {}).get('status') == 'complete':
        check_files(branch / 'downstream', recovery['stages']['orthology']['files'])
        protected.update({'validation/downstream/' + k.replace('\\', '/'): v
                          for k, v in recovery['stages']['orthology']['files'].items()})
    elif not (orth / 'inference_complete.json').exists():
        results = list((orth / 'inference').glob('Results_*'))
        if len(results) != 1 or set(p.name for p in results[0].iterdir()) != {'Log.txt', 'WorkingDirectory'}:
            raise ValueError('Expected an interrupted search-only OrthoFinder attempt')
        wd = results[0] / 'WorkingDirectory'
        log = (results[0] / 'Log.txt').read_text()
        if 'Started OrthoFinder version 3.1.5' not in log or 'Search program: diamond' not in log:
            raise ValueError('Partial-search adapter is restricted to the inspected OrthoFinder 3.1.5 / DIAMOND run')
        if not all(value in log for value in ('Scoring matrix: BLOSUM62', 'Gap extend: 1', 'Gap open: 11',
                                              'MSA program: famsa', 'Tree program: fasttree')):
            raise ValueError('Search or tree configuration differs from the inspected checkpoint')
        if any(recovery['settings'].get('orthofinder_' + name) is not None
               for name in ('inflation', 'method', 'search')):
            raise ValueError('Custom inference settings require a separate restart review')
        failed_log = '\n'.join(p.read_text(errors='replace') for p in (orth / 'logs').glob('command_*.log'))
        interrupted = sorted(set(re.findall(r'Blast(\d+_\d+)\.txt', failed_log)))
        if 'KeyboardInterrupt' not in failed_log:
            raise ValueError('Partial search interruption was not verified')
        for path in wd.iterdir():
            if path.is_file() and not path.name.startswith('Blast'):
                protected[path.relative_to(run).as_posix()] = sha(path)
        for path in wd.glob('Blast*.txt.gz'):
            pair = path.name.removeprefix('Blast').removesuffix('.txt.gz')
            if pair not in interrupted and intact_search(path):
                searches[pair] = sha(path)
                protected[path.relative_to(run).as_posix()] = searches[pair]
        if (wd / 'blast_commands.txt').exists():
            raise ValueError('Unexpected queued search commands; review before preparing a restart')
        bind_search_panel(wd, branch / 'downstream/panel')
    plan = dict(version='1.0', status='pending_local_resume', prepared_at=now(),
                topology_source=key, topology_sha256=sha(source),
                previous_topology_sha256=sources[key], added_ids=added,
                candidate_count=len(sequences), previous_topology_count=len(old),
                annotation_manifest_before=sha(annotation / 'workflow_manifest.json'),
                protected_files=protected, completed_searches=searches,
                interrupted_searches=interrupted, pipeline_executed=False, backups_created=False,
                phase='queued', script_sha256=sha(Path(__file__)))
    if not apply:
        print(json.dumps({k: v for k, v in plan.items() if k != 'protected_files'} |
                         dict(status='inspection_only', protected_files_checked=len(protected)), indent=2))
        return
    dump(run / PLAN, plan)
    recovery['pending_annotation_refresh'] = dict(plan='../context_refresh.json',
        source_sha256=plan['topology_sha256'], added_candidates=len(added))
    dump(branch / 'downstream/recovery.json', recovery)
    stages = load(branch / 'stages.json')
    stages['annotation'].update(status='pending_context_refresh', execution='pending_local_resume')
    if stages.get('orthology_audit', {}).get('status') == 'running':
        stages['orthology_audit']['status'] = 'interrupted'
    dump(branch / 'stages.json', stages)
    dump(branch / 'branch_status.json', dict(status='ready_to_resume', stages=stages,
        reason='All topology records returned; selective context refresh queued. Calibration is untouched.'))
    summary = load(run / 'run_summary.json')
    summary['status'] = 'ready_to_resume'; summary['gates']['validation'] = load(branch / 'branch_status.json')
    summary['context_refresh_prepared_at'] = now()
    dump(run / 'run_summary.json', summary)
    helper = 'maintenance/resume_validation_context.py'
    (run / 'README.md').write_text('# Workflow 1.0: validation\n\nStatus: ready to resume locally.\n\n'
        f'All {len(sequences)} candidates now have exact-sequence topology returns; {len(added)} were added.\n'
        'The old annotation is retained as a sealed prior result until the local refresh.\n'
        'Calibration is complete. Validation references/panel are complete. Validation orthology was interrupted during searches.\n\n'
        'From the pipeline installation directory, run:\n\n```bash\n'
        f'python -B {helper} --run runs/{run.name} --resume --threads 4\n```\n\n'
        'This helper skips calibration entirely, refreshes annotation, reuses completed orthology, '
        'and completes only missing/interrupted DIAMOND searches when validation inference is unfinished. '
        'It then runs the remaining benchmark/audit/reports. No backups or archived attempts are made.\n'
        'Use this command for the pending context refresh rather than the ordinary resume command.\n', encoding='utf-8')
    check_files(run, protected)
    print(json.dumps({k: v for k, v in plan.items() if k != 'protected_files'}, indent=2))


def remove_stale(folder, parent):
    """Discard only obsolete annotation products, inside the explicitly checked scope."""
    if folder.is_symlink() or not folder.resolve().is_relative_to(parent.resolve()):
        raise ValueError('Refusing to remove an out-of-scope annotation product')
    if folder.is_dir():
        if any(p.is_symlink() for p in folder.rglob('*')):
            raise ValueError('Symlink in obsolete annotation output')
        shutil.rmtree(folder)
    elif folder.exists():
        folder.unlink()


def resume(run, threads, root):
    if os.name == 'nt':
        raise ValueError('Run --resume inside the existing WSL conda environment')
    plan = load(run / PLAN)
    if sha(Path(__file__)) != plan['script_sha256']:
        raise ValueError('Maintenance helper changed since checkpoint preparation')
    check_files(run, plan['protected_files'])
    if sha(local(plan['topology_source'])) != plan['topology_sha256']:
        raise ValueError('Topology return changed since preparation')
    sys.path.insert(0, str(root)); os.chdir(root)
    import adhesome_workflow as core
    from schisto_adhesome import branch as adapter
    from schisto_adhesome.workflows import verify_frozen
    verify_frozen(load(run / 'run_contract.json'))
    if plan['status'] == 'complete':
        print('Validation already completed; nothing to resume')
        return
    os.environ['SCHISTO_NO_ARCHIVES'] = '1'
    os.environ['SCHISTO_ISOLATED_BRANCH'] = 'validation'
    annotation = run / 'validation/annotation'; downstream = run / 'validation/downstream'
    if plan['phase'] in ('queued', 'refreshing_annotation'):
        if plan['phase'] == 'queued':
            if sha(annotation / 'workflow_manifest.json') != plan['annotation_manifest_before']:
                raise ValueError('Annotation changed since context-refresh preparation')
            verify_seal(annotation, 'annotation')
            plan['phase'] = 'refreshing_annotation'; dump(run / PLAN, plan)
        elif (annotation / 'workflow_manifest.json').exists():
            # Annotation finished before a process interruption; do not run it twice.
            meta = verify_seal(annotation, 'annotation')
            if meta['sources'].get(plan['topology_source']) == plan['topology_sha256']:
                plan['phase'] = 'annotation_refreshed'; dump(run / PLAN, plan)
        if plan['phase'] == 'refreshing_annotation':
            remove_stale(annotation / 'motifs', annotation)
            remove_stale(annotation / 'workflow_manifest.json', annotation)
            state = load(annotation / 'annotation_recovery.json')
            state['status'] = 'failed'; state['contract']['sources'][plan['topology_source']] = plan['topology_sha256']
            dump(annotation / 'annotation_recovery.json', state)
            code = core.main(['resume', '--handoff', str(run / 'validation/handoff'),
                '--out', str(annotation), '--context-dir', str(local(plan['topology_source']).parent),
                '--annotation-only', '--retry-failed', '--threads', str(threads)])
            if code:
                raise ValueError('Annotation refresh failed; inspect its latest command log, then retry this helper')
            plan['phase'] = 'annotation_refreshed'; dump(run / PLAN, plan)
    if plan['phase'] == 'annotation_refreshed':
        meta = verify_seal(annotation, 'annotation')
        if meta['sources'].get(plan['topology_source']) != plan['topology_sha256']:
            raise ValueError('Refreshed annotation did not adopt the returned topology')
        state = load(downstream / 'recovery.json')
        state['annotation_sha256'] = sha(annotation / 'workflow_manifest.json')
        state.pop('pending_annotation_refresh', None)
        state['context_refresh'] = dict(source_sha256=plan['topology_sha256'], added_ids=plan['added_ids'], completed_at=now())
        dump(downstream / 'recovery.json', state)
        plan['phase'] = 'downstream'; dump(run / PLAN, plan)
    original_orthology = core.orthology

    def preserved_orthology(args):
        out = args.out.resolve(); receipt = out / 'inference_complete.json'
        if receipt.exists():
            saved = load(receipt)
            if saved.get('maintenance_resume'):
                if saved['panel_manifest'] != sha(args.panel / 'workflow_manifest.json'):
                    raise ValueError('Inference panel changed')
                check_files(out, saved['files'])
                results = Path(saved['results'])
            else:
                return original_orthology(args)
        else:
            results = finish_searches(core, out, args.panel, plan, run, threads)
            saved = dict(maintenance_resume=True, panel_manifest=sha(args.panel / 'workflow_manifest.json'),
                         results=str(results), files={p.relative_to(out).as_posix(): sha(p)
                         for p in (out / 'inference').rglob('*') if p.is_file()}, completed_at=now())
            dump(receipt, saved)
        args.results = results; args.input_manifest = args.panel / 'workflow_manifest.json'
        # Candidate extraction has not completed if this action is reached.
        # Never archive or remove its partial output automatically.
        original_orthology(args)

    core.orthology = preserved_orthology
    args = argparse.Namespace(kind='validation', run=run, threads=threads,
        context_dir=local(plan['topology_source']).parent, allow_missing_context=False,
        reuse_detection=None, compare_report=None)
    print('REUSE calibration: branch will not be launched', flush=True)
    code = adapter.execute(args)
    check_files(run, plan['protected_files'])
    if code:
        raise ValueError('Validation branch did not finish')
    plan.update(status='complete', phase='complete', completed_at=now(), pipeline_executed=True)
    dump(run / PLAN, plan)
    summary = load(run / 'run_summary.json')
    summary['gates']['validation'] = load(run / 'validation/branch_status.json')
    summary['status'] = 'complete'; summary['completed_at'] = now()
    dump(run / 'run_summary.json', summary)
    (run / 'README.md').write_text('# Workflow 1.0: validation\n\nStatus: complete\n\n'
        '- Calibration: reused without launching its branch.\n'
        '- Validation: [report](validation/reports/report.html).\n'
        '- Topology: all returned candidates imported; [refresh receipt](validation/context_refresh.json).\n', encoding='utf-8')


def orthofinder_version(text):
    """Read both 'OrthoFinder:v3.1.5' and the ordinary version banner."""
    plain = re.sub(r'\x1b\[[0-?]*[ -/]*[@-~]', '', text)
    match = re.search(r'\bOrthoFinder(?:\s+version)?\s*:?\s*v?\s*'
                      r'(\d+\.\d+\.\d+)(?![\d.])', plain, flags=re.IGNORECASE)
    return match.group(1) if match else None


def finish_searches(core, out, panel, plan, run, threads):
    """Only the inspected 3.1.5 search checkpoint; no inference parameters are tuned."""
    roots = sorted((out / 'inference').glob('Results_*'))
    finished = [p for p in roots if (p / 'Orthogroups/Orthogroups.tsv').exists()
                and any((p / 'Orthologues').glob('**/*.tsv'))
                and 'OrthoFinder run completed' in (p / 'Log.txt').read_text(errors='replace')]
    if len(finished) == 1:
        return finished[0]
    # Additional partial post-search inference needs explicit review; do not rerun it.
    if len(roots) != 1:
        raise ValueError('A subsequent OrthoFinder attempt is partial; inspect logs before selecting its restart stage')
    version = subprocess.run(['orthofinder', '-v'], capture_output=True, text=True)
    detected = orthofinder_version(version.stdout + '\n' + version.stderr)
    if version.returncode:
        detail = re.sub(r'\x1b\[[0-?]*[ -/]*[@-~]', '', version.stderr or version.stdout).strip()
        raise ValueError(f'OrthoFinder version query failed (exit {version.returncode}): {detail[-1500:]}')
    if detected != '3.1.5':
        raise ValueError('Keep OrthoFinder 3.1.5 for this saved inference checkpoint; '
                         f'detected: {detected or "unrecognized version output"}')
    old = roots[0]; wd = old / 'WorkingDirectory'
    if (wd / 'blast_commands.txt').exists():
        raise ValueError('Unexpected queued OrthoFinder search commands; refusing possible duplicate searches')
    # Bind working sequences to the sealed panel, not only species labels.
    species = bind_search_panel(wd, panel)
    progress_file = out / 'search_resume.json'
    progress = load(progress_file) if progress_file.exists() else dict(completed=plan['completed_searches'].copy())
    missing = []
    for query in species:
        for target in species:
            pair = f'{query}_{target}'; path = wd / f'Blast{pair}.txt.gz'
            if pair in progress['completed']:
                if sha(path) != progress['completed'][pair]:
                    raise ValueError('Completed search changed: ' + pair)
            else:
                missing.append((query, target))
    print(f'REUSE {len(progress["completed"])} search pairs; finish {len(missing)} missing/interrupted pairs', flush=True)
    logs = out / 'search_resume_logs'; logs.mkdir(exist_ok=True)

    def search(pair):
        query, target = pair; key = f'{query}_{target}'; path = wd / f'Blast{key}.txt.gz'
        command = ['diamond', 'blastp', '--ignore-warnings', '-d', str(wd / f'diamondDBSpecies{target}'),
            '-q', str(wd / f'Species{query}.fa'), '-o', str(wd / f'Blast{key}.txt'),
            '--matrix', 'BLOSUM62', '--gapopen', '11', '--gapextend', '1', '--more-sensitive',
            '-p', '1', '--quiet', '-e', '0.001', '--compress', '1']
        with (logs / f'Blast{key}.log').open('a') as handle:
            proc = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT)
        if proc.returncode or not intact_search(path):
            raise ValueError('Search failed: ' + key + '; inspect ' + str(logs))
        return key, sha(path)

    with ThreadPoolExecutor(max_workers=threads) as pool:
        jobs = [pool.submit(search, pair) for pair in missing]
        for job in as_completed(jobs):
            key, value = job.result(); progress['completed'][key] = value
            dump(progress_file, progress)
    if len(progress['completed']) != len(species) ** 2:
        raise ValueError('Missing search results; OrthoFinder continuation refused')
    command = ['orthofinder', '-b', str(old), '-t', str(threads), '-a', str(threads)]
    core.execute(command, out / 'logs')
    candidates = [p for p in (out / 'inference').glob('Results_*') if p != old]
    if (len(candidates) != 1 or not (candidates[0] / 'Orthogroups/Orthogroups.tsv').exists()
            or not any((candidates[0] / 'Orthologues').glob('**/*.tsv'))
            or 'OrthoFinder run completed' not in (candidates[0] / 'Log.txt').read_text(errors='replace')):
        raise ValueError('Cannot identify complete OrthoFinder continuation results')
    return candidates[0]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--pipeline-root', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--threads', type=int, default=4)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--inspect', action='store_true', help='Read-only input/checkpoint integrity inspection; no files changed')
    mode.add_argument('--prepare', action='store_true', help='Validate returned inputs and queue the checkpoint only; no analysis')
    mode.add_argument('--resume', action='store_true', help='Explicitly run the pending local validation work; skip calibration')
    args = parser.parse_args()
    if args.threads < 1:
        parser.error('--threads must be positive')
    try:
        if args.prepare or args.inspect:
            prepare(args.run.resolve(), apply=args.prepare)
        else:
            resume(args.run.resolve(), args.threads, args.pipeline_root.resolve())
    except (OSError, ValueError, KeyError) as exc:
        parser.exit(1, 'ERROR: ' + str(exc) + '\n')


if __name__ == '__main__':
    main()
