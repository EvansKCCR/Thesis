"""Build consensus single-copy gene trees and species trees from saved inference.

This is separate from the optional adhesome-family tree workflow. It uses the
benchmark-core occupancy table, MAFFT L-INS-i, trimAl gappyout, IQ-TREE ML gene
trees, a concatenated IQ-TREE analysis, and either ASTRAL-IV or ASTRAL-III.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def read_tsv(path):
    with Path(path).open(encoding='utf-8-sig', newline='') as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        if not reader.fieldnames:
            raise ValueError(f'Empty table: {path}')
        return reader.fieldnames, list(reader)


def read_fasta(path):
    records = {}
    name, chunks = None, []
    with Path(path).open(encoding='ascii') as stream:
        for line in stream:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if name is not None:
                    if name in records:
                        raise ValueError(f'Duplicate FASTA ID {name} in {path}')
                    records[name] = ''.join(chunks)
                name, chunks = line[1:].split()[0], []
            elif name is None:
                raise ValueError(f'Sequence before FASTA header in {path}')
            else:
                chunks.append(line.upper())
    if name is not None:
        records[name] = ''.join(chunks)
    if not records:
        raise ValueError(f'No sequences in {path}')
    return records


def write_fasta(path, records):
    with Path(path).open('w', encoding='ascii', newline='\n') as stream:
        for name, seq in records.items():
            stream.write(f'>{name}\n')
            for start in range(0, len(seq), 80):
                stream.write(seq[start:start + 80] + '\n')


def parse_alignment(path):
    return read_fasta(path)


def run_logged(command, logfile):
    with Path(logfile).open('w', encoding='utf-8') as stream:
        stream.write('COMMAND: ' + subprocess.list2cmdline([str(x) for x in command]) + '\n\n')
        stream.flush()
        result = subprocess.run([str(x) for x in command], stdout=stream, stderr=subprocess.STDOUT)
    if result.returncode:
        raise RuntimeError(f'Command exited {result.returncode}; inspect {logfile}')


def executable(name):
    path = shutil.which(name)
    if not path:
        raise ValueError(f'Required executable not found on PATH: {name}')
    return path


def prepare(args):
    orthology = args.orthology.resolve()
    panel = args.panel.resolve()
    inference = orthology / 'inference'
    if not (inference / 'core_consensus.tsv').is_file():
        raise ValueError('No core consensus table; orthology must use benchmark-core assignment')
    if not (inference / 'final_results.txt').is_file():
        raise ValueError('No final OrthoFinder result pointer in saved orthology run')
    result_path = Path((inference / 'final_results.txt').read_text().strip())
    group_table_path = result_path / 'Orthogroups' / 'Orthogroups.tsv'
    if not group_table_path.is_file():
        raise ValueError(f'Full predicted orthogroup table is missing: {group_table_path}')
    _, core_rows = read_tsv(inference / 'core_consensus.tsv')
    consensus = {row['orthogroup'] for row in core_rows}
    fields, rows = read_tsv(group_table_path)
    group_field = next((key for key in ('OG', 'HOG', 'Orthogroup') if key in fields), None)
    if not group_field:
        raise ValueError('Cannot identify orthogroup ID field in final OrthoFinder table')
    snapshot = json.loads((panel / 'panel.snapshot.json').read_text(encoding='utf-8'))
    entries = snapshot['proteomes']
    core_species = [entry['species'] for entry in entries if entry.get('orthology_role') == 'benchmark_core']
    target_species = [entry['species'] for entry in entries if entry.get('orthology_role') == 'target']
    outgroups = [entry['species'] for entry in entries if entry.get('orthology_role') == 'reference_outgroup']
    if len(core_species) != 12 or not target_species:
        raise ValueError(f'Expected 12 benchmark core proteomes and at least one configured target; found {len(core_species)} core and {len(target_species)} targets')
    all_species = core_species + target_species + outgroups
    if set(species_columns(fields)) != set(all_species):
        raise ValueError('Final orthogroup table species do not match the sealed three-tier panel')
    minimum_bench = args.minimum_benchmark_presence
    minimum_target = args.minimum_target_presence if args.minimum_target_presence is not None else min(5,len(target_species))
    minimum_total = args.minimum_total_taxa if args.minimum_total_taxa is not None else min(15,minimum_bench+minimum_target)
    selected = []
    for row in rows:
        og = row[group_field]
        if og not in consensus:
            continue
        membership = {sp: [member.strip() for member in (row.get(sp) or '').split(',') if member.strip()]
                      for sp in all_species}
        if any(len(membership[sp]) > 1 for sp in core_species + target_species):
            continue
        if sum(bool(membership[sp]) for sp in core_species) < minimum_bench:
            continue
        if sum(bool(membership[sp]) for sp in target_species) < minimum_target:
            continue
        if sum(bool(membership[sp]) for sp in core_species + target_species) < minimum_total:
            continue
        # Outgroups are included as tips when their inferred group has at most
        # one member; their occupancy does not affect the single-copy gate.
        if any(len(membership[sp]) > 1 for sp in outgroups):
            continue
        selected.append((og, membership))
    if not selected:
        raise ValueError('No orthogroups pass consensus occupancy and strict single-copy criteria')
    out = args.out.resolve()
    if out.exists() and not args.resume:
        raise ValueError(f'Output exists; use --resume to continue this species-tree project: {out}')
    (out / 'groups').mkdir(parents=True, exist_ok=True)
    sequence_maps = {entry['species']: read_fasta(panel / 'proteomes' / f"{entry['species']}.faa") for entry in entries}
    report = []
    for og, membership in selected:
        records = {}
        for species, members in membership.items():
            if members:
                pid = members[0]
                if pid not in sequence_maps[species]:
                    raise ValueError(f'{og}: OrthoFinder member {pid} missing from {species} prepared proteome')
                records[pid] = sequence_maps[species][pid]
        folder = out / 'groups' / og
        folder.mkdir(exist_ok=True)
        raw = folder / 'input.faa'
        write_fasta(raw, records)
        report.append(dict(orthogroup=og, benchmark_occupancy=sum(bool(membership[sp]) for sp in core_species),
                           target_occupancy=sum(bool(membership[sp]) for sp in target_species),
                           outgroup_occupancy=sum(bool(membership[sp]) for sp in outgroups),
                           taxa=len(records), input=str(raw)))
    with (out / 'groups.tsv').open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(report[0]), delimiter='\t')
        writer.writeheader(); writer.writerows(report)
    provenance = dict(pipeline_version='1.0',status='prepared', orthology_result=str(result_path), group_table_sha256=digest(group_table_path),
                      core_consensus_sha256=digest(inference / 'core_consensus.tsv'),
                      panel_snapshot_sha256=digest(panel / 'panel.snapshot.json'),
                      group_inputs_sha256={row['orthogroup']:digest(row['input']) for row in report},
                      benchmark_species=core_species, target_species=target_species, outgroups=outgroups,
                      gates=dict(minimum_benchmark_presence=minimum_bench, minimum_target_presence=minimum_target,
                                 minimum_total_taxa=minimum_total, single_copy='at most one protein per included species'),
                      selected_groups=len(report), groups_sha256=digest(out / 'groups.tsv'))
    (out / 'manifest.json').write_text(json.dumps(provenance, indent=2) + '\n', encoding='utf-8')
    print(f'Selected {len(report)} strict single-copy consensus orthogroups -> {out}', flush=True)


def tree_group(folder, args):
    raw = folder / 'input.faa'
    aligned = folder / 'aligned.faa'
    trimmed = folder / 'trimmed.faa'
    tree = folder / 'gene_tree.treefile'
    state_path = folder / 'tree_state.json'
    input_sha = digest(raw)
    if args.resume and state_path.is_file():
        state = json.loads(state_path.read_text(encoding='utf-8'))
        if state.get('input_sha256') != input_sha:
            raise ValueError(f'Input changed since the saved tree run for {folder.name}')
        outputs = {'aligned.faa':aligned,'trimmed.faa':trimmed,'gene_tree.treefile':tree}
        if all(path.is_file() and digest(path)==state.get('outputs_sha256',{}).get(name) for name,path in outputs.items()):
            return folder.name, 'reused'
        raise ValueError(f'Saved tree outputs failed checksum validation for {folder.name}; inspect before retrying')
    result = subprocess.run([args.mafft, '--localpair', '--maxiterate', '1000', str(raw)],
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    (folder / 'mafft.log').write_text('COMMAND: ' + subprocess.list2cmdline([args.mafft,'--localpair','--maxiterate','1000',str(raw)]) + '\n\n' + result.stderr, encoding='utf-8')
    if result.returncode:
        (folder / 'mafft.log').write_text(result.stderr, encoding='utf-8')
        raise RuntimeError(f'MAFFT failed for {folder.name}; inspect {folder / "mafft.log"}')
    aligned.write_text(result.stdout, encoding='ascii')
    run_logged([args.trimal, '-in', aligned, '-out', trimmed, '-gappyout'], folder / 'trimal.log')
    if not trimmed.is_file() or not trimmed.stat().st_size:
        raise RuntimeError(f'trimAl produced no alignment for {folder.name}')
    run_logged([args.iqtree, '-s', trimmed, '-m', 'MFP', '-B', '1000', '-T', '1', '--prefix', folder / 'gene_tree'],
               folder / 'iqtree.log')
    if not tree.is_file():
        raise RuntimeError(f'IQ-TREE did not produce {tree}')
    state=dict(input_sha256=input_sha,outputs_sha256={'aligned.faa':digest(aligned),'trimmed.faa':digest(trimmed),
                'gene_tree.treefile':digest(tree)})
    state_path.write_text(json.dumps(state,indent=2)+'\n',encoding='utf-8')
    return folder.name, 'complete'


def concatenate(project):
    rows = []
    _, groups = read_tsv(project / 'groups.tsv')
    species_order = []
    for row in groups:
        alignment = row['orthogroup']
        seqs = parse_alignment(project / 'groups' / alignment / 'trimmed.faa')
        if not seqs:
            continue
        lengths = {len(value) for value in seqs.values()}
        if len(lengths) != 1:
            raise ValueError(f'Unequal sequence lengths in trimmed alignment {alignment}')
        species_order.extend(name.split('__', 1)[0] for name in seqs if name.split('__',1)[0] not in species_order)
        rows.append((alignment, seqs, lengths.pop()))
    concatenated = {sp: [] for sp in species_order}
    for _, seqs, width in rows:
        by_species = {name.split('__',1)[0]: seq for name,seq in seqs.items()}
        for sp in species_order:
            concatenated[sp].append(by_species.get(sp, '-' * width))
    output = {sp: ''.join(parts) for sp,parts in concatenated.items()}
    write_fasta(project / 'concatenated_trimmed.faa', output)
    with (project / 'concatenation.tsv').open('w', encoding='utf-8', newline='') as stream:
        writer=csv.writer(stream,delimiter='\t');writer.writerow(['orthogroup','alignment_length','offset_start_1based','offset_end_1based'])
        offset=1
        for og,_,width in rows:
            writer.writerow([og,width,offset,offset+width-1]);offset+=width


def run(args):
    project = args.project.resolve()
    manifest_path = project / 'manifest.json'
    if not manifest_path.is_file():
        raise ValueError('Run prepare before run')
    manifest = json.loads(manifest_path.read_text(encoding='utf-8'))
    if manifest.get('groups_sha256') != digest(project / 'groups.tsv'):
        raise ValueError('Prepared group membership changed; prepare a new project')
    for og, expected in manifest.get('group_inputs_sha256',{}).items():
        if digest(project/'groups'/og/'input.faa') != expected:
            raise ValueError(f'Prepared sequence input changed for {og}; prepare a new project')
    if not args.resume and manifest.get('status') == 'complete':
        raise ValueError('Species-tree run is already complete; use a new output directory')
    mafft, trimal, iqtree = executable(args.mafft), executable(args.trimal), executable(args.iqtree)
    astral4 = shutil.which(args.astral4) if args.astral4 else None
    astral3 = Path(args.astral3_jar).resolve() if args.astral3_jar else None
    if not astral4 and not (astral3 and astral3.is_file() and shutil.which(args.java)):
        raise ValueError('Install ASTRAL-IV (`astral4`) or supply --astral3-jar before starting this computationally intensive run')
    args.mafft, args.trimal, args.iqtree = mafft, trimal, iqtree
    folders = sorted((project / 'groups').iterdir())
    manifest['status']='running'
    manifest_path.write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
    statuses=[]
    with ThreadPoolExecutor(max_workers=args.threads) as pool:
        futures={pool.submit(tree_group,folder,args):folder.name for folder in folders}
        for index,future in enumerate(as_completed(futures),1):
            og=futures[future]
            try:
                _,status=future.result()
                statuses.append(dict(orthogroup=og,status=status,error=''))
            except Exception as exc:
                statuses.append(dict(orthogroup=og,status='failed',error=str(exc)))
            print(f'STAGE GENE_TREES {index}/{len(folders)} {og}: {statuses[-1]["status"]}',flush=True)
    with (project/'gene_tree_status.tsv').open('w',encoding='utf-8',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=['orthogroup','status','error'],delimiter='\t');writer.writeheader();writer.writerows(statuses)
    failed=[row for row in statuses if row['status']=='failed']
    if failed:
        raise ValueError(f'{len(failed)} gene-tree jobs failed; inspect gene_tree_status.tsv and resume after review')
    concatenate(project)
    outgroups=manifest['outgroups']
    concat_command=[iqtree,'-s',project/'concatenated_trimmed.faa','-m','MFP','-B','1000','-T',str(args.iqtree_threads),'--prefix',project/'concatenated_tree']
    if outgroups:
        concat_command.extend(['-o',','.join(outgroups)])
    run_logged(concat_command,project/'concatenated_iqtree.log')
    gene_trees=project/'gene_trees.tre'
    gene_trees.write_text(''.join((project/'groups'/og/'gene_tree.treefile').read_text(encoding='utf-8').strip()+'\n' for og in sorted(folder.name for folder in folders)),encoding='utf-8')
    if astral4:
        command=[astral4,'-i',gene_trees,'-o',project/'species_tree.astral4.tre']
        run_logged(command,project/'astral4.log')
        species_tree='species_tree.astral4.tre'
        method='ASTRAL-IV'
    else:
        command=[args.java,'-jar',astral3,'-i',gene_trees,'-o',project/'species_tree.astral3.tre']
        run_logged(command,project/'astral3.log')
        species_tree='species_tree.astral3.tre'
        method='ASTRAL-III'
    if not (project/species_tree).is_file():
        raise ValueError(f'{method} did not produce the expected species tree: {project/species_tree}')
    manifest.update(status='complete',gene_tree_count=len(folders),failed_gene_trees=0,
                    gene_trees_sha256=digest(gene_trees),concatenated_alignment_sha256=digest(project/'concatenated_trimmed.faa'),
                    coalescent_method=method,coalescent_tree=species_tree,
                    coalescent_tree_rooting='unrooted; outgroup tips are retained for downstream rooting',
                    concatenated_tree='concatenated_tree.treefile',group_status_sha256=digest(project/'gene_tree_status.tsv'))
    manifest_path.write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
    print(f'STAGE SPECIES_TREES COMPLETE: {project}',flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest='command',required=True)
    p=sub.add_parser('prepare')
    for name in ('orthology','panel','out'):
        p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--minimum-benchmark-presence',type=int,default=10)
    p.add_argument('--minimum-target-presence',type=int,default=None,help='Default is min(5, number of configured target species)')
    p.add_argument('--minimum-total-taxa',type=int,default=None,help='Default is min(15, benchmark plus minimum target presence)')
    p.set_defaults(func=prepare)
    p=sub.add_parser('run')
    p.add_argument('--project',type=Path,required=True)
    p.add_argument('--threads',type=int,default=4)
    p.add_argument('--iqtree-threads',type=int,default=2)
    p.add_argument('--mafft',default='mafft')
    p.add_argument('--trimal',default='trimal')
    p.add_argument('--iqtree',default='iqtree2')
    p.add_argument('--astral4',default='astral4')
    p.add_argument('--astral3-jar',type=Path)
    p.add_argument('--java',default='java')
    p.add_argument('--resume',action='store_true')
    p.set_defaults(func=run)
    args=parser.parse_args()
    try:
        if args.command=='prepare':
            prepare(args)
        else:
            if args.threads<1 or args.iqtree_threads<1:
                raise ValueError('Thread counts must be positive')
            run(args)
    except (OSError,ValueError,RuntimeError,KeyError) as exc:
        parser.exit(1,f'ERROR: {exc}\n')


def species_columns(fields):
    metadata={'OG','HOG','Orthogroup','Gene Tree Parent Clade','Gene_Tree_Parent_Clade'}
    return [field for field in fields if field not in metadata]


if __name__=='__main__':
    main()
