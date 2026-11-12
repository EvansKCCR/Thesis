#!/usr/bin/env python3
"""Eight-species preparation, guarded OrthoFinder execution, and candidate extraction."""

# Resolve source resources independently of the caller working directory.
import sys as _layout_sys
from pathlib import Path as _LayoutPath
_LAYOUT_BUNDLE_ROOT = next(p for p in _LayoutPath(__file__).resolve().parents
                          if (p/'config/layout_paths.json').is_file() and (p/'pipeline_paths.py').is_file())
if str(_LAYOUT_BUNDLE_ROOT) not in _layout_sys.path:
    _layout_sys.path.insert(0, str(_LAYOUT_BUNDLE_ROOT))
from pipeline_paths import (activate as _layout_activate, resource_path as _layout_resource,
                            default_thread_count as _default_thread_count,
                            warn_slow_wsl_mounts as _warn_slow_wsl_mounts)
_layout_activate()

import argparse
import csv
import gzip
import hashlib
import json
import re
import subprocess
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

ROOT = _layout_resource(Path(__file__)).resolve().parents[1]
VERSION = '1.0'
LEGACY_SPECIES = dict(Hsap='Homo_sapiens.faa', Cele='Caenorhabditis_elegans.fa',
               Shae='Schistosoma_haematobium.fa', Sjap='Schistosoma_japonicum.fa',
               Sman='Schistosoma_mansoni.fa')
SPECIES = dict(LEGACY_SPECIES, Dmel='Drosophila_melanogaster.faa',
               Mmus='Mus_musculus.faa', Xlae='Xenopus_laevis.faa')
SCHISTOSOMA = ('Shae', 'Sjap', 'Sman')
AA = set('ACDEFGHIKLMNPQRSTVWYBXZJUO')
CUSTOM_SPECIES = None


def configure_species(path):
    """Accept an explicit panel for extraction; never infer a panel from table columns."""
    global SPECIES, CUSTOM_SPECIES
    if path is None:
        return
    codes = json.loads(_layout_resource(Path(path)).read_text(encoding='utf-8'))
    if (not isinstance(codes, list) or len(codes) != len(set(codes))
            or len(codes) < 2
            or any(not re.fullmatch(r'[A-Za-z][A-Za-z0-9_]*', s) for s in codes)):
        raise ValueError('Species list must contain at least two unique safe codes')
    CUSTOM_SPECIES = codes
    SPECIES = {s: s+'.faa' for s in codes}


def profile_paths(profile):
    if profile == 'legacy5':
        return LEGACY_SPECIES, _layout_resource(ROOT/'03_qc'), _layout_resource(ROOT/'02_primary_proteomes/prepared')
    if profile == 'expanded8':
        return SPECIES, _layout_resource(ROOT/'03_qc/expanded8'), _layout_resource(ROOT/'02_primary_proteomes/prepared_8species')
    raise ValueError(f'Unknown species profile: {profile}')


def table_species(fields):
    active = [s for s in SPECIES if s in (fields or [])]
    allowed = (set(CUSTOM_SPECIES),) if CUSTOM_SPECIES is not None else (set(LEGACY_SPECIES), set(SPECIES))
    if set(active) not in allowed:
        raise ValueError('Expected a complete root table for the configured species panel')
    if CUSTOM_SPECIES is not None:
        extra = set(fields or []) - set(CUSTOM_SPECIES) - {'HOG','OG','Orthogroup','Gene Tree Parent Clade','Gene_Tree_Parent_Clade','Matched_candidates'}
        if extra:
            raise ValueError(f'Unexpected columns/species outside the configured panel: {sorted(extra)}')
    return active


def digest(path):
    with _layout_resource(path).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def records(path):
    with _layout_resource(path).open('rb') as f:
        compressed = f.read(2) == b'\x1f\x8b'
    opener = gzip.open if compressed else open
    with opener(path, 'rt', encoding='utf-8-sig') as f:
        header, seq = None, []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(seq)
                header, seq = line[1:], []
                if not header.strip():
                    raise ValueError(f'{path}: empty header')
            else:
                if header is None:
                    raise ValueError(f'{path}: sequence before header')
                seq.append(''.join(line.split()).upper())
        if header is not None:
            yield header, ''.join(seq)


def write_tsv(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with _layout_resource(path).open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fields, delimiter='\t', extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)


def read_tsv(path):
    with _layout_resource(path).open(encoding='utf-8-sig', newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        return reader.fieldnames, list(reader)


def candidates(path):
    ids = []
    for line in _layout_resource(path).read_text(encoding='utf-8-sig').splitlines():
        if not line.strip() or line.lstrip().startswith('#'):
            continue
        token = line.strip().removeprefix('>').split()[0]
        if token.lower() == 'protein_id':
            continue
        if not re.fullmatch(r'(Shae|Sjap|Sman)__[A-Za-z0-9_.-]+', token):
            raise ValueError(f'{path}: invalid candidate ID {token}')
        ids.append(token)
    return sorted(set(ids)), len(ids)


def prepare(args):
    species, qc, out = profile_paths(args.profile)
    qc.mkdir(parents=True, exist_ok=True)
    out.mkdir(parents=True, exist_ok=True)
    # A prior successful preparation cannot authorize a failed new preparation.
    manifest_path = _layout_resource(qc / 'preparation.json')
    manifest_path.write_text(json.dumps({'ready': False, 'status': 'preparing'}))
    mapping = {}
    if args.human_map:
        fields, rows = read_tsv(args.human_map)
        if not {'protein_id', 'gene_id'} <= set(fields or []):
            raise ValueError('Human mapping needs protein_id and gene_id TSV columns')
        for row in rows:
            p, g = row['protein_id'].strip(), row['gene_id'].strip()
            if not p or not g or (p in mapping and mapping[p] != g):
                raise ValueError(f'Invalid or conflicting human mapping: {p}')
            mapping[p] = g
    summary, idmap, missing = [], [], []
    manifest = {'pipeline_version':VERSION,'ready': False, 'created_utc': datetime.now(timezone.utc).isoformat(),
                'selection': 'longest sequence; ties by original protein ID', 'inputs': {}, 'outputs': {},
                'profile': args.profile, 'species': list(species), 'annotations': {}}
    if args.human_map:
        manifest['human_map'] = {'path': str(_layout_resource(args.human_map).resolve()), 'sha256': digest(args.human_map)}
    for sp, name in species.items():
        source = _layout_resource(ROOT / '01_raw_proteomes') / name
        manifest['inputs'][sp] = {'file': name, 'sha256': digest(source)}
        species_map = mapping if sp == 'Hsap' and args.human_map else None
        if sp in ('Hsap', 'Dmel', 'Mmus', 'Xlae') and species_map is None:
            from refseq_gene_map import parse_gff
            annotation = source.with_suffix('.gff')
            print(f'{sp}: reading {annotation.name}', flush=True)
            gene_sets, metadata = parse_gff(annotation)
            conflicts = {pid: sorted(gs) for pid, gs in gene_sets.items() if len(gs) != 1}
            if conflicts:
                raise ValueError(f'{sp}: conflicting GFF gene mappings: {list(conflicts.items())[:5]}')
            species_map = {pid: 'GeneID_'+next(iter(gs)) for pid, gs in gene_sets.items()}
            manifest['annotations'][sp] = {'file': annotation.name, 'sha256': digest(annotation),
                                          'metadata': metadata}
            write_tsv(_layout_resource(qc/f'{sp}_protein_gene.tsv'), ['protein_id','gene_id'],
                      [dict(protein_id=p, gene_id=g) for p,g in sorted(species_map.items())])
        groups, seen = defaultdict(list), set()
        count, stops, errors = 0, 0, []
        for header, seq in records(source):
            pid = header.split()[0]
            count += 1
            if pid in seen:
                errors.append(f'duplicate ID: {pid}')
            seen.add(pid)
            if seq.endswith('*'):
                seq = seq[:-1]
                stops += 1
            if not seq or set(seq) - AA:
                errors.append(f'empty/invalid sequence: {pid} ({"".join(sorted(set(seq)-AA))})')
            match = re.search(r'(?:^|\s)gene[=:](\S+)', header)
            gene = species_map.get(pid) if species_map is not None else (match.group(1) if match else None)
            if not gene:
                missing.append({'species': sp, 'protein_id': pid})
                continue
            if not re.fullmatch(r'[A-Za-z0-9_.-]+', gene):
                errors.append(f'unsafe gene ID: {gene}')
            groups[gene].append((pid, seq))
        if count == 0:
            errors.append('empty proteome')
        if errors:
            raise ValueError(f'{sp}: ' + '; '.join(errors[:20]))
        species_missing = sum(r['species'] == sp for r in missing)
        if not species_missing:
            target = _layout_resource(out / f'{sp}.faa')
            with _layout_resource(target).open('w', encoding='ascii', newline='\n') as f:
                for gene, variants in sorted(groups.items()):
                    chosen = sorted(variants, key=lambda x: (-len(x[1]), x[0]))[0]
                    stable = f'{sp}__{gene}'
                    f.write(f'>{stable}\n')
                    f.write('\n'.join(chosen[1][i:i+80] for i in range(0, len(chosen[1]), 80)) + '\n')
                    for pid, seq in variants:
                        idmap.append(dict(species=sp, original_id=pid, gene_id=gene, prepared_id=stable,
                                          representative_id=chosen[0], selected=int(pid == chosen[0]), length=len(seq)))
            manifest['outputs'][sp] = {'file': target.name, 'sha256': digest(target)}
        summary.append(dict(species=sp, raw_sequences=count, genes=len(groups), missing_gene_mapping=species_missing,
                            terminal_stops_removed=stops, status='blocked' if species_missing else 'prepared'))
    write_tsv(_layout_resource(qc/'proteome_summary.tsv'), list(summary[0]), summary)
    write_tsv(_layout_resource(qc/'missing_gene_mapping.tsv'), ['species', 'protein_id'], missing)
    write_tsv(_layout_resource(qc/'identifier_map.tsv'), ['species','original_id','gene_id','prepared_id','representative_id','selected','length'], idmap)
    aliases = defaultdict(set)
    for r in idmap:
        for alias in (r['prepared_id'], f"{r['species']}__{r['original_id']}"):
            aliases[alias].add(r['prepared_id'])
    audits, unresolved = [], 0
    for file in sorted((_layout_resource(ROOT/'candidate_ids')).glob('*.txt')):
        if file.stem.endswith('.clean'):
            continue
        ids, raw_count = candidates(file)
        if file.name == 'fibronectin_like_protein.txt' and len(ids) != 89:
            raise ValueError(f'Expected 89 fibronectin-like IDs, found {len(ids)}')
        for cid in ids:
            hits = aliases[cid]
            status = 'matched' if len(hits) == 1 else ('ambiguous' if hits else 'missing')
            unresolved += status != 'matched'
            audits.append(dict(candidate_set=file.stem, candidate_id=cid,
                               prepared_id=next(iter(hits)) if len(hits)==1 else '', status=status))
        print(f'{file.name}: {len(ids)} unique IDs ({raw_count-len(ids)} duplicates)')
    write_tsv(_layout_resource(qc/'candidate_audit.tsv'), ['candidate_set','candidate_id','prepared_id','status'], audits)
    manifest['ready'] = not missing and not unresolved
    manifest['unresolved_candidates'] = unresolved
    manifest['candidate_inputs'] = {p.name: digest(p) for p in sorted((_layout_resource(ROOT/'candidate_ids')).glob('*.txt'))}
    manifest_path.write_text(json.dumps(manifest, indent=2)+'\n')
    for row in summary:
        print(row)
    print(f'Inference ready: {manifest["ready"]}; unresolved candidates: {unresolved}')
    return 0 if manifest['ready'] else 2


def run(args):
    species, qc, folder = profile_paths(args.profile)
    manifest = json.loads((_layout_resource(qc/'preparation.json')).read_text())
    if not manifest['ready']:
        raise ValueError(f'Preparation is blocked; inspect {qc}.')
    if set(manifest['outputs']) != set(species) or set(manifest['inputs']) != set(species):
        raise ValueError('Preparation manifest does not match the requested species profile')
    if set(p.name for p in folder.iterdir()) != {f'{s}.faa' for s in species}:
        raise ValueError(f'Prepared directory must contain exactly {len(species)} expected FASTAs')
    for sp, record in manifest['outputs'].items():
        if digest(folder/record['file']) != record['sha256']:
            raise ValueError(f'{sp}: prepared FASTA changed; rerun preparation')
    for sp, record in manifest['inputs'].items():
        if digest(_layout_resource(ROOT/'01_raw_proteomes')/record['file']) != record['sha256']:
            raise ValueError(f'{sp}: raw input changed; rerun preparation')
    for sp, record in manifest.get('annotations', {}).items():
        if digest(_layout_resource(ROOT/'01_raw_proteomes')/record['file']) != record['sha256']:
            raise ValueError(f'{sp}: annotation changed; rerun preparation')
    if args.threads < 1 or args.analysis_threads < 1:
        raise ValueError('Thread counts must be positive')
    dest = _layout_resource(ROOT/'04_results')/args.name
    if _layout_resource(dest).exists() or _layout_resource(Path(args.name)).name != args.name:
        raise ValueError('Use a new, simple run name; existing runs cannot be overwritten')
    cmd = ['orthofinder','-f',str(folder),'-t',str(args.threads),'-a',str(args.analysis_threads),'-o',str(dest)]
    print('Command:', cmd)
    if args.dry_run:
        return 0
    help_result = subprocess.run(['orthofinder','-h'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    if help_result.returncode:
        raise ValueError('OrthoFinder help/version check failed: '+help_result.stdout)
    logdir = qc/args.name
    logdir.mkdir(exist_ok=False)
    (_layout_resource(logdir/'orthofinder_help.txt')).write_text(help_result.stdout)
    (_layout_resource(logdir/'run.json')).write_text(json.dumps({'command':cmd,'preparation':manifest},indent=2))
    with (_layout_resource(logdir/'orthofinder.log')).open('w') as log:
        result = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT)
    (_layout_resource(logdir/'exit_code.txt')).write_text(str(result.returncode)+'\n')
    return result.returncode


def members(value):
    return {v.strip() for v in value.split(',') if v.strip()}


def resolve_group_table(result, requested):
    if requested != 'auto':
        return _layout_resource(result/'Phylogenetic_Hierarchical_Orthogroups')/requested
    traditional = _layout_resource(result/'Phylogenetic_Hierarchical_Orthogroups/N0.tsv')
    if _layout_resource(traditional).exists():
        return traditional
    log = _layout_resource(result/'Log.txt')
    version = re.search(r'OrthoFinder version (\d+)\.(\d+)', _layout_resource(log).read_text() if _layout_resource(log).exists() else '')
    modern = _layout_resource(result/'Orthogroups/Orthogroups.tsv')
    if version and tuple(map(int, version.groups())) >= (3, 1) and _layout_resource(modern).exists():
        return modern
    native = _layout_resource(result/'Orthogroups/Orthogroups.txt')
    if version and tuple(map(int, version.groups())) >= (3, 1) and _layout_resource(native).exists():
        return native
    raise ValueError('No root HOG/orthogroup table found. Expected N0.tsv or an OrthoFinder >=3.1 orthogroup output.')


def read_native_orthogroups(path, species):
    """Adapt OrthoFinder's OG: member-list output to the TSV rows used by extraction."""
    observed = set()
    rows = []
    with _layout_resource(path).open(encoding='utf-8-sig') as stream:
        for line_number, line in enumerate(stream, 1):
            line = line.strip()
            if not line:
                continue
            group, separator, values = line.partition(':')
            members_in_group = values.split()
            if not separator or not group.strip().startswith('OG') or not members_in_group:
                raise ValueError(f'Malformed OrthoFinder orthogroup line {line_number}: {path}')
            row = {'OG': group.strip()}
            grouped = defaultdict(list)
            for protein in members_in_group:
                matches = [code for code in species if protein.startswith(code+'__')]
                if len(matches) != 1:
                    raise ValueError(f'Cannot map OrthoFinder member {protein!r} to one configured species in {path}')
                grouped[matches[0]].append(protein)
                observed.add(matches[0])
            row.update({code:', '.join(grouped.get(code,[])) for code in species})
            rows.append(row)
    if not rows:
        raise ValueError(f'No orthogroups found in {path}')
    return ['OG', *[code for code in species if code in observed]], rows


def pairwise_tables(file):
    fields, rows = read_tsv(file)
    cols = [s for s in SPECIES if s in (fields or [])]
    if len(cols) == 2:
        return [(cols, rows)]
    if len(cols) == 1 and {'Species', 'Orthologs'} <= set(fields or []):
        source = cols[0]
        grouped = defaultdict(list)
        for row in rows:
            target = row['Species']
            if target not in SPECIES or target == source:
                raise ValueError(f'Unexpected target species in {file}: {target}')
            grouped[target].append(dict(row, **{target: row['Orthologs']}))
        return [([source, target], group) for target, group in grouped.items()]
    raise ValueError(f'Unexpected orthologue schema: {file}')


def extract(args):
    result = _layout_resource(args.results).resolve()
    hogfile = resolve_group_table(result, args.hog)
    print(f'Using group table: {hogfile}')
    fields, hogs = (read_native_orthogroups(hogfile, list(CUSTOM_SPECIES or SPECIES))
                    if hogfile.name=='Orthogroups.txt' else read_tsv(hogfile))
    active_species = table_species(fields)
    targets_path = getattr(args, 'target_species', None)
    targets = json.loads(_layout_resource(targets_path).read_text(encoding='utf-8')) if targets_path else list(SCHISTOSOMA)
    if not targets or len(targets)!=len(set(targets)) or not set(targets) < set(active_species):
        raise ValueError('Target species must be a nonempty subset of the panel with at least one reference')
    references = [s for s in active_species if s not in targets]
    profile = 'custom' if CUSTOM_SPECIES is not None else 'expanded8' if len(active_species) == 8 else 'legacy5'
    audit_file = getattr(args, 'candidate_audit', None)
    if profile == 'custom' and audit_file is None:
        raise ValueError('Custom species panels require --candidate-audit')
    if audit_file is None:
        _, qc, _ = profile_paths(profile)
        audit_file = _layout_resource(qc/'candidate_audit.tsv')
    _, audit = read_tsv(audit_file)
    subset = [r for r in audit if r['candidate_set'] == args.candidate_set]
    if not subset:
        raise ValueError('Unknown/empty candidate set')
    wanted = {r['prepared_id'] for r in subset if r['status']=='matched'}
    if any(cid.split('__',1)[0] not in targets for cid in wanted):
        raise ValueError('Candidate mapping includes a species not declared as a target')
    index = defaultdict(list)
    selected = []
    for row in hogs:
        hits = set.union(*(members(row[s]) for s in active_species)) & wanted
        if hits:
            selected.append(dict(row, Matched_candidates=', '.join(sorted(hits))))
            for cid in hits:
                index[cid].append(row)
    pairs, links, seen = [], defaultdict(list), set()
    pairfiles = sorted((_layout_resource(result/'Orthologues')).rglob('*.tsv'))
    if not pairfiles:
        raise ValueError('No pairwise orthologue tables found; results may be incomplete')
    for file, cols, rows in ((file, cols, rows) for file in pairfiles for cols, rows in pairwise_tables(file)):
        for sp in targets:
            for ref in references:
                if set(cols) != {sp,ref}:
                    continue
                for row in rows:
                    source, target = members(row[sp]), members(row[ref])
                    if not target:
                        continue
                    relation = ('one' if len(source)==1 else 'many')+'-to-'+('one' if len(target)==1 else 'many')
                    for cid in source & wanted:
                        key = (cid,ref,tuple(sorted(source)),tuple(sorted(target)))
                        if key in seen:
                            continue
                        seen.add(key)
                        item = dict(candidate_id=cid, reference_species=ref, orthologues=', '.join(sorted(target)),
                                    relationship=relation, source_group_size=len(source), target_group_size=len(target),
                                    source_file=str(file.relative_to(result)))
                        pairs.append(item)
                        links[cid,ref].append(item)
    report = []
    for candidate in subset:
        cid = candidate['prepared_id']
        groups = index[cid]
        row = dict(candidate_id=candidate['candidate_id'], prepared_id=cid, mapping_status=candidate['status'],
                   HOG='; '.join(g.get('HOG',g.get('Orthogroup','')) for g in groups),
                   Orthogroup='; '.join(g.get('OG',g.get('Orthogroup','')) for g in groups))
        for ref in references:
            direct = links[cid,ref]
            row[ref+'_orthologues'] = ', '.join(sorted(set().union(*(members(p['orthologues']) for p in direct))))
            row[ref+'_relationship'] = '; '.join(sorted({p['relationship'] for p in direct})) if direct else (
                'candidate absent from prepared inputs' if not cid else
                'orthogroup member without direct pairwise orthologue' if any(members(g[ref]) for g in groups) else
                'no direct pairwise orthologue detected')
        row['HOG_status'] = ('not mapped to input' if not cid else 'not found in selected HOG table' if not groups else
                             'target-only in sampled HOG' if all(not members(g[ref]) for g in groups for ref in references) else 'shared HOG')
        for sp in active_species:
            row[sp+'_HOG_copy_counts'] = '; '.join(str(len(members(g[sp]))) for g in groups)
        report.append(row)
    out = _layout_resource(args.output).resolve()
    if _layout_resource(out).exists():
        raise ValueError('Extraction output exists; choose a new directory')
    out.mkdir(parents=True)
    write_tsv(_layout_resource(out/'candidate_HOGs.tsv'),fields+['Matched_candidates'],selected)
    write_tsv(_layout_resource(out/'direct_orthologues.tsv'), ['candidate_id','reference_species','orthologues','relationship','source_group_size','target_group_size','source_file'],pairs)
    write_tsv(_layout_resource(out/'candidate_report.tsv'),list(report[0]),report)
    dupfile = _layout_resource(result/'Gene_Duplication_Events/Duplications.tsv')
    if _layout_resource(dupfile).exists():
        df, rows = read_tsv(dupfile)
        gene_cols = [c for c in df if c.startswith('Genes')]
        write_tsv(_layout_resource(out/'candidate_duplications.tsv'), df, [r for r in rows if any(members(r[c]) & wanted for c in gene_cols)])
    (_layout_resource(out/'provenance.json')).write_text(json.dumps({'pipeline_version':VERSION,'results':str(result),'hog_table':str(hogfile),
        'hog_sha256':digest(hogfile),'candidate_set':args.candidate_set,
        'species': active_species, 'target_species': targets, 'profile': profile, 'candidate_audit_sha256':digest(audit_file),
        'pairwise_files':{str(p.relative_to(result)):digest(p) for p in pairfiles}},indent=2)+'\n')
    print(f'{len(report)} candidates; {len(selected)} HOG rows; {len(pairs)} direct orthology relationships -> {out}')
    return 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--version',action='version',version='adhesome-orthology '+VERSION)
    commands = parser.add_subparsers(dest='command', required=True)
    p = commands.add_parser('prepare')
    p.add_argument('--profile',choices=['expanded8','legacy5'],default='expanded8')
    p.add_argument('--human-map', type=Path)
    p.set_defaults(func=prepare)
    p = commands.add_parser('run')
    p.add_argument('--profile',choices=['expanded8','legacy5'],default='expanded8')
    p.add_argument('--name',required=True)
    p.add_argument('--threads',type=int,default=_default_thread_count(8),
                   help='CPU budget; defaults to the scheduler allocation or up to 8 local CPUs')
    p.add_argument('--analysis-threads',type=int,default=None,
                   help='OrthoFinder algorithm CPU budget; defaults to --threads')
    p.add_argument('--dry-run',action='store_true')
    p.set_defaults(func=run)
    p = commands.add_parser('extract')
    p.add_argument('--results',type=Path,required=True)
    p.add_argument('--hog',default='auto',help='auto detects the version-specific root table; or specify N0.tsv etc.')
    p.add_argument('--candidate-set',default='fibronectin_like_protein')
    p.add_argument('--candidate-audit', type=Path, help='Explicit mapping audit for a workflow run')
    p.add_argument('--species-list', type=Path, help='JSON list of all species codes in a custom panel')
    p.add_argument('--target-species', type=Path, help='JSON list of target species; other panel species are references')
    p.add_argument('--output',type=Path,required=True)
    p.set_defaults(func=extract)
    args = parser.parse_args()
    _warn_slow_wsl_mounts(value for value in vars(args).values() if isinstance(value,Path))
    if args.command=='run' and args.analysis_threads is None:
        args.analysis_threads=args.threads
    try:
        configure_species(getattr(args, 'species_list', None))
        return args.func(args)
    except (ValueError,OSError,KeyError) as exc:
        print(f'ERROR: {exc}',file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
