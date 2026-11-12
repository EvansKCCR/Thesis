"""Run OrthoFinder on the benchmark core, then assign the remaining panel taxa.

The final OrthoFinder result remains the full predicted orthogroup set for
benchmark scoring. A separate core-consensus table records the groups present
in at least the configured number of benchmark species.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def read_table(path: Path):
    with path.open(encoding='utf-8-sig', newline='') as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        if not reader.fieldnames:
            raise ValueError(f'Empty OrthoFinder group table: {path}')
        return reader.fieldnames, list(reader)


def species_columns(fields):
    metadata = {'OG', 'HOG', 'Orthogroup', 'Gene Tree Parent Clade', 'Gene_Tree_Parent_Clade'}
    return [field for field in fields if field not in metadata]


def group_table(result: Path) -> Path:
    candidates = [result / 'Orthogroups' / 'Orthogroups.tsv',
                  result / 'Phylogenetic_Hierarchical_Orthogroups' / 'N0.tsv',
                  result / 'Orthogroups' / 'Orthogroups.txt']
    for path in candidates:
        if path.is_file():
            return path
    raise ValueError(f'No supported OrthoFinder orthogroup table in {result}')


def read_group_table(path: Path, species: set[str]):
    """Read OrthoFinder's tabular output or its native OG: member-list file."""
    if path.name != 'Orthogroups.txt':
        return read_table(path)
    observed = set()
    fields = ['OG', *sorted(species)]
    rows = []
    with path.open(encoding='utf-8-sig') as stream:
        for line_number, line in enumerate(stream, 1):
            line = line.strip()
            if not line:
                continue
            group, separator, members_text = line.partition(':')
            if not separator or not group.strip().startswith('OG'):
                raise ValueError(f'Malformed OrthoFinder group line {line_number} in {path}')
            row = {field: '' for field in fields}
            row['OG'] = group.strip()
            members = members_text.split()
            if not members:
                raise ValueError(f'OrthoFinder group has no members at line {line_number} in {path}')
            for member in members:
                matches = [name for name in species if member.startswith(name + '__')]
                if len(matches) != 1:
                    raise ValueError(f'Cannot map OrthoFinder member {member!r} to one panel species in {path}')
                name = matches[0]
                observed.add(name)
                row[name] = (row[name] + ' ' + member).strip()
            rows.append(row)
    if not rows:
        raise ValueError(f'No orthogroups found in {path}')
    fields = ['OG', *sorted(observed)]
    return fields, rows


def find_results(root: Path, species: set[str]) -> Path:
    matches = []
    seen = set()
    for table in [*root.rglob('Orthogroups.tsv'), *root.rglob('Orthogroups.txt')]:
        result = table.parent.parent
        if result in seen:
            continue
        seen.add(result)
        try:
            fields, _ = read_group_table(table, species)
        except (OSError, ValueError):
            continue
        if set(species_columns(fields)) == species:
            matches.append(result)
    if len(matches) != 1:
        raise ValueError(f'Expected one OrthoFinder result with {len(species)} species; found {len(matches)}')
    return matches[0]


def find_results_optional(root: Path, species: set[str]):
    matches = []
    seen = set()
    for table in [*root.rglob('Orthogroups.tsv'), *root.rglob('Orthogroups.txt')]:
        result = table.parent.parent
        if result in seen:
            continue
        seen.add(result)
        try:
            fields, _ = read_group_table(table, species)
        except (OSError, ValueError):
            continue
        if set(species_columns(fields)) == species:
            matches.append(result)
    return matches


def run_command(command, logfile, label):
    logfile.parent.mkdir(parents=True, exist_ok=True)
    print(f'STAGE {label} START', flush=True)
    with logfile.open('w', encoding='utf-8') as stream:
        stream.write('COMMAND: ' + subprocess.list2cmdline([str(x) for x in command]) + '\n\n')
        stream.flush()
        result = subprocess.run([str(x) for x in command], stdout=stream, stderr=subprocess.STDOUT)
    if result.returncode:
        raise RuntimeError(f'{label} failed with exit status {result.returncode}; inspect {logfile}')
    log_text = logfile.read_text(encoding='utf-8', errors='replace')
    reported_errors = [line.strip() for line in log_text.splitlines() if line.lstrip().startswith('ERROR:')]
    if reported_errors:
        raise RuntimeError(f'{label} reported an error despite exit status 0: {reported_errors[0]}; inspect {logfile}')
    traceback_logged = 'Traceback (most recent call last):' in log_text
    unequal_alignment_logged = 'Sequences must all be the same length' in log_text
    if traceback_logged or unequal_alignment_logged:
        detail = next(
            (line.strip() for line in reversed(log_text.splitlines())
             if line.lstrip().startswith(('AttributeError:', 'ImportError:', 'ModuleNotFoundError:', 'TypeError:', 'ValueError:'))),
            'Sequences must all be the same length' if unequal_alignment_logged
            else 'Python traceback found in command log',
        )
        raise RuntimeError(f'{label} logged an internal error despite exit status 0: {detail}; inspect {logfile}')
    print(f'STAGE {label} COMPLETE', flush=True)


def consensus_table(core_result: Path, core_species: list[str], minimum: int, output: Path):
    path = group_table(core_result)
    fields, rows = read_group_table(path, set(core_species))
    absent = set(core_species) - set(species_columns(fields))
    if absent:
        raise ValueError('Core OrthoFinder table lacks benchmark species: ' + ', '.join(sorted(absent)))
    group_field = next((name for name in ('OG', 'HOG', 'Orthogroup') if name in fields), None)
    if group_field is None:
        raise ValueError(f'Cannot identify orthogroup ID column in {path}')
    selected = []
    for row in rows:
        present = [sp for sp in core_species if (row.get(sp) or '').strip()]
        if len(present) >= minimum:
            selected.append(dict(orthogroup=row[group_field], benchmark_species_present=len(present),
                                 benchmark_species=';'.join(present),
                                 core_members=';'.join(f'{sp}:{(row.get(sp) or "").strip()}' for sp in present)))
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=['orthogroup','benchmark_species_present','benchmark_species','core_members'], delimiter='\t')
        writer.writeheader()
        writer.writerows(selected)
    return path, len(rows), len(selected)


def panel_quality_report(entries, inference: Path):
    bm = [entry for entry in entries if entry.get('orthology_role') == 'benchmark_core']
    rows = []
    for entry in bm:
        try:
            provenance = json.loads(entry.get('source') or '{}')
        except json.JSONDecodeError:
            provenance = {}
        scope = provenance.get('reference_scope') or {}
        level = scope.get('assembly_level') or scope.get('assemblyLevel') or ''
        completeness = scope.get('completeness_percent') or scope.get('busco_completeness') or ''
        rows.append(dict(species=entry['species'],source_fasta=Path(provenance.get('fasta','')).name,
                         assembly_level=level,completeness=completeness,
                         assembly_metadata_status='provided' if level and completeness != '' else 'unverified',
                         assembly_category_status='eligible' if str(level).casefold() in ('chromosome','scaffold','chromosome-level','scaffold-level') else 'unverified_or_ineligible'))
    report = inference / 'benchmark_panel_quality_review.tsv'
    with report.open('w', encoding='utf-8', newline='') as stream:
        fields = ['species','source_fasta','assembly_level','completeness','assembly_metadata_status','assembly_category_status']
        writer=csv.DictWriter(stream,fieldnames=fields,delimiter='\t');writer.writeheader();writer.writerows(rows)
    targets=[entry for entry in entries if entry.get('orthology_role')=='target']
    outgroups=[entry for entry in entries if entry.get('orthology_role')=='reference_outgroup']
    summary = dict(status='metadata_unverified' if not all(row['assembly_metadata_status']=='provided' for row in rows) else 'metadata_present',
                   benchmark_species=len(bm),target_species=len(targets),reference_species=len(outgroups),
                   assembly_metadata_complete=all(row['assembly_metadata_status']=='provided' for row in rows),
                   eligible_assembly_count=sum(row['assembly_category_status']=='eligible' for row in rows),
                   report=str(report))
    (inference / 'benchmark_panel_quality_review.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
    if not summary['assembly_metadata_complete']:
        print('PANEL QUALITY REVIEW: assembly metadata incomplete; see benchmark_panel_quality_review.tsv', flush=True)
    return summary


def fasta_layout(path: Path):
    """Read sequence identifiers and lengths without loading large FASTAs in memory."""
    identifiers=[];lengths=[];identifier=None;length=0
    with path.open(encoding='utf-8',errors='replace') as stream:
        for raw in stream:
            line=raw.strip()
            if line.startswith('>'):
                if identifier is not None:lengths.append(length)
                words=line[1:].split();identifier=words[0] if words else ''
                if not identifier or identifier in identifiers:
                    raise ValueError(f'Empty or duplicate sequence ID in {path}: {identifier!r}')
                identifiers.append(identifier);length=0
            elif line:
                if identifier is None:raise ValueError(f'Sequence before first FASTA header in {path}')
                length+=len(line)
    if identifier is not None:lengths.append(length)
    if not identifiers or any(length<1 for length in lengths):
        raise ValueError(f'No sequences or an empty sequence in {path}')
    return identifiers,lengths


def repair_core_alignments(core_result: Path, threads: int, report: Path):
    """Preflight cached core MSAs; deterministically realign malformed files from source sequences."""
    working=core_result/'WorkingDirectory'
    alignments=working/'Alignments_ids';sequences=working/'Sequences_ids'
    if not alignments.is_dir() or not sequences.is_dir():
        raise ValueError(f'Core OrthoFinder result is missing reusable MSA/source directories: {working}')
    malformed=[]
    for aligned in sorted(alignments.glob('*.fa')):
        ids,lengths=fasta_layout(aligned)
        if len(ids)>1 and len(set(lengths))!=1:
            malformed.append((aligned,ids,lengths))
    if not malformed:
        report.parent.mkdir(parents=True,exist_ok=True)
        with report.open('w',encoding='utf-8',newline='') as stream:
            writer=csv.DictWriter(stream,fieldnames=['orthogroup','status','sequence_count','alignment_length','source_sha256','previous_alignment_sha256','aligned_sha256'],delimiter='\t')
            writer.writeheader()
        return []
    executable=shutil.which('famsa')
    if not executable:
        names=', '.join(path.name for path,_,_ in malformed[:10])
        raise ValueError(f'{len(malformed)} cached core alignments have unequal sequence lengths ({names}); install FAMSA to repair them before assignment')
    repairs=[]
    for aligned,expected_ids,old_lengths in malformed:
        source=sequences/aligned.name
        if not source.is_file():
            raise ValueError(f'Cannot repair malformed core alignment without its original sequences: {source}')
        source_ids,_=fasta_layout(source)
        if set(source_ids)!=set(expected_ids):
            raise ValueError(f'Cannot safely repair {aligned.name}: alignment/source sequence membership differs')
        temporary=aligned.with_name(aligned.stem+'.repaired.tmp.fa')
        result=subprocess.run([executable,'-t',str(max(1,threads)),str(source),str(temporary)],
            stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,encoding='utf-8',errors='replace')
        if result.returncode:
            temporary.unlink(missing_ok=True)
            raise RuntimeError(f'FAMSA could not repair {aligned.name} (exit {result.returncode}): {result.stdout[-2000:]}')
        repaired_ids,repaired_lengths=fasta_layout(temporary)
        if set(repaired_ids)!=set(expected_ids) or (len(repaired_ids)>1 and len(set(repaired_lengths))!=1):
            temporary.unlink(missing_ok=True)
            raise ValueError(f'FAMSA output failed ID/equal-length validation for {aligned.name}')
        source_hash=sha256(source);previous_hash=sha256(aligned);aligned_hash=sha256(temporary)
        temporary.replace(aligned)
        repairs.append(dict(orthogroup=aligned.stem,status='realigned_with_famsa',sequence_count=len(repaired_ids),
            alignment_length=repaired_lengths[0],source_sha256=source_hash,previous_alignment_sha256=previous_hash,
            aligned_sha256=aligned_hash))
    report.parent.mkdir(parents=True,exist_ok=True)
    with report.open('w',encoding='utf-8',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=['orthogroup','status','sequence_count','alignment_length','source_sha256','previous_alignment_sha256','aligned_sha256'],delimiter='\t')
        writer.writeheader();writer.writerows(repairs)
    print(f'CORE MSA PREFLIGHT: repaired {len(repairs)} malformed alignments from saved source sequences; see {report}',flush=True)
    return repairs


def execute(args):
    panel = args.panel.resolve()
    inference = args.out.resolve()
    snapshot_path = panel / 'panel.snapshot.json'
    snapshot = json.loads(snapshot_path.read_text(encoding='utf-8'))
    entries = snapshot.get('proteomes', [])
    core = [entry for entry in entries if entry.get('orthology_role') == 'benchmark_core']
    added = [entry for entry in entries if entry.get('orthology_role') != 'benchmark_core']
    if len(core) != 12:
        raise ValueError(f'Core assignment requires exactly 12 benchmark proteomes; found {len(core)}')
    targets = [entry for entry in entries if entry.get('orthology_role') == 'target']
    outgroups = [entry for entry in entries if entry.get('orthology_role') == 'reference_outgroup']
    if not targets or len(added) != len(targets) + len(outgroups) or not added:
        raise ValueError(f'Panel requires at least one user-selected target and valid target/reference roles; found {len(targets)} targets and {len(outgroups)} references')
    if args.minimum_core_presence < 1 or args.minimum_core_presence > len(core):
        raise ValueError('Core-presence threshold must be between 1 and 12')
    core_codes = [entry['species'] for entry in core]
    all_codes = [entry['species'] for entry in entries]
    if len(set(all_codes)) != len(all_codes):
        raise ValueError('Species codes in panel are not unique')
    input_hashes = {}
    for entry in entries:
        fasta = (panel / entry['fasta']).resolve() if not Path(entry['fasta']).is_absolute() else Path(entry['fasta'])
        if not fasta.is_file():
            raise ValueError(f"Prepared proteome missing for {entry['species']}: {fasta}")
        input_hashes[entry['species']] = sha256(fasta)
    settings = dict(threads=args.threads, analysis_threads=args.analysis_threads,
                    orthofinder_options=args.option, minimum_core_presence=args.minimum_core_presence,
                    core_species=core_codes, all_species=all_codes, input_hashes=input_hashes,
                    panel_snapshot_sha256=sha256(snapshot_path))
    state_path = inference / 'core_assignment_state.json'
    if inference.exists() and not args.resume:
        raise ValueError(f'Inference output already exists; use --resume only for a matching interrupted core-assignment run: {inference}')
    inference.mkdir(parents=True, exist_ok=True)
    quality=panel_quality_report(entries,inference)
    if state_path.exists():
        state = json.loads(state_path.read_text(encoding='utf-8'))
        if state.get('settings') != settings:
            raise ValueError('Core-assignment checkpoint settings or panel inputs changed; use a new output directory')
    else:
        state = dict(settings=settings, status='preparing')
    core_input = inference / 'core_input'
    added_input = inference / 'added_input'
    if not state.get('inputs_prepared'):
        for folder, group in ((core_input, core), (added_input, added)):
            folder.mkdir(parents=True, exist_ok=True)
            for entry in group:
                source = (panel / entry['fasta']).resolve() if not Path(entry['fasta']).is_absolute() else Path(entry['fasta'])
                target = folder / f"{entry['species']}.faa"
                if target.exists():
                    if sha256(target) != sha256(source):
                        raise ValueError(f'Interrupted input staging contains a changed file: {target}')
                else:
                    shutil.copyfile(source, target)
        state['inputs_prepared'] = True
        state['status'] = 'inputs_prepared'
        state_path.write_text(json.dumps(state, indent=2) + '\n', encoding='utf-8')
    help_file = inference / 'orthofinder_help.txt'
    help_result = subprocess.run(['orthofinder', '-h'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    if help_result.returncode or '--assign' not in help_result.stdout or '--core' not in help_result.stdout:
        raise ValueError('Installed OrthoFinder must support --assign and --core for benchmark-core inference')
    help_file.write_text(help_result.stdout, encoding='utf-8')
    common = []
    if args.option:
        common.extend(args.option)
    core_result = Path(state['core_result']) if state.get('core_result') else None
    if core_result is None or not core_result.is_dir() or not group_table(core_result).is_file():
        core_out = inference / 'core_run'
        recovered_core = find_results_optional(core_out, set(core_codes)) if args.resume and core_out.exists() else []
        if len(recovered_core) > 1:
            raise ValueError('Multiple core OrthoFinder results found; inspect the checkpoint and use a new output if needed')
        if recovered_core:
            core_result = recovered_core[0]
            print(f'STAGE ORTHOLOGY CORE REUSE: {core_result}', flush=True)
        else:
            command = ['orthofinder', '-f', core_input, '-t', str(args.threads), '-a', str(args.analysis_threads), '-o', core_out, *common]
            run_command(command, inference / 'core.log', 'ORTHOLOGY CORE (12 benchmark species)')
            core_result = find_results(core_out, set(core_codes))
        state.update(core_result=str(core_result), core_group_table_sha256=sha256(group_table(core_result)), status='core_complete')
        state_path.write_text(json.dumps(state, indent=2) + '\n', encoding='utf-8')
    elif state.get('core_group_table_sha256') != sha256(group_table(core_result)):
        raise ValueError('Core orthogroup table changed since checkpoint')
    _, core_group_count, consensus_count = consensus_table(core_result, core_codes, args.minimum_core_presence,
                                                            inference / 'core_consensus.tsv')
    state.update(core_group_count=core_group_count, consensus_group_count=consensus_count,
                 minimum_core_presence=args.minimum_core_presence, panel_quality=quality)
    state_path.write_text(json.dumps(state, indent=2) + '\n', encoding='utf-8')
    final_result = Path(state['final_result']) if state.get('final_result') else None
    if final_result is not None and final_result.is_dir():
        try:
            final_fields, _ = read_group_table(group_table(final_result), set(all_codes))
            final_species = set(species_columns(final_fields))
            complete_panel = final_species == set(all_codes)
        except (OSError, ValueError):
            complete_panel = False
        if not complete_panel:
            # Earlier checkpoints could mistake an Orthogroups.txt core file for
            # a full-panel result because their parser filled absent species columns.
            state.pop('final_result', None)
            state.pop('final_group_table_sha256', None)
            state['status'] = 'core_complete'
            state_path.write_text(json.dumps(state, indent=2) + '\n', encoding='utf-8')
            final_result = None
    if final_result is None or not final_result.is_dir():
        recovered = find_results_optional(inference, set(all_codes)) if args.resume else []
        if len(recovered) > 1:
            raise ValueError('Multiple completed full-panel OrthoFinder results found; select a clean new output directory')
        if recovered:
            final_result = recovered[0]
            state.update(final_result=str(final_result),
                         final_group_table_sha256=sha256(group_table(final_result)), status='complete')
            state_path.write_text(json.dumps(state, indent=2) + '\n', encoding='utf-8')
    if final_result is None or not final_result.is_dir():
        repair_core_alignments(core_result,args.threads,inference/'alignment_repairs.tsv')
        command = ['orthofinder', '--assign', added_input, '--core', core_result,
                   '-t', str(args.threads), '-a', str(args.analysis_threads), *common]
        run_command(command, inference / 'assign.log', f'ORTHOLOGY ASSIGN (remaining {len(added)} species)')
        final_result = find_results(inference, set(all_codes))
        state.update(final_result=str(final_result),
                     final_group_table_sha256=sha256(group_table(final_result)), status='complete')
        state_path.write_text(json.dumps(state, indent=2) + '\n', encoding='utf-8')
    (inference / 'final_results.txt').write_text(str(final_result.resolve()) + '\n', encoding='utf-8')
    print(f'Core consensus: {consensus_count}/{core_group_count} orthogroups present in >= {args.minimum_core_presence}/12 benchmark taxa', flush=True)
    print(f'Full {len(all_codes)}-species OrthoFinder result: {final_result}', flush=True)
    return final_result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--panel', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--threads', type=int, default=4)
    parser.add_argument('--analysis-threads', type=int, default=4)
    parser.add_argument('--minimum-core-presence', type=int, default=10)
    parser.add_argument('--option', action='append', default=[], help='One OrthoFinder option token; repeat for key and value')
    parser.add_argument('--resume', action='store_true')
    args = parser.parse_args()
    try:
        if args.threads < 1 or args.analysis_threads < 1:
            raise ValueError('Thread counts must be positive')
        execute(args)
    except (OSError, ValueError, RuntimeError, KeyError) as exc:
        parser.exit(1, f'ERROR: {exc}\n')


if __name__ == '__main__':
    main()
