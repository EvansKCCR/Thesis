#!/usr/bin/env python3
"""Recursive, source-matched protein homologue discovery and reciprocal BLASTP.

Python >=3.9 and NCBI BLAST+ are required. No third-party Python dependencies.
Reference files define independent reference sets; target counts are unrestricted.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
from datetime import datetime, timezone
import math
import json
import re
import shutil
import subprocess
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

EXTS = ('.fa', '.faa', '.fasta', '.fas', '.fsa', '.pep', '.aa')
FIELDS = ('qseqid', 'sseqid', 'pident', 'length', 'qlen', 'slen',
          'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore')
OUTFMT = '6 ' + ' '.join(FIELDS)
VALID_AA = set('ACDEFGHIKLMNPQRSTVWYBXZJUO*')
RBH_STATUSES = {'RBH_candidate','RBH_candidate_augmented_reference','provisional_reference_only_RBH'}


@dataclass(frozen=True)
class Protein:
    original_id: str
    description: str
    seq: str
    internal_id: str


def open_text(path: Path):
    return gzip.open(path, 'rt', encoding='utf-8') if path.suffix.lower() == '.gz' else path.open(encoding='utf-8')


def read_fasta(path: Path, prefix: str, qc_rows=None, role="input") -> list[Protein]:
    records: list[Protein] = []
    used: set[str] = set()
    desc = None
    chunks = []
    stop_records = []

    def append_one(header, sequence_lines):
        if header is None:
            return
        if not header or not header.split():
            raise ValueError(f'Empty FASTA header in {path}')
        original_id = header.split()[0]
        seq = ''.join(sequence_lines).replace(' ', '').upper()
        if not original_id or not seq:
            raise ValueError(f'Empty header or sequence in {path}: {header!r}')
        if original_id in used:
            raise ValueError(f'Duplicate FASTA ID {original_id!r} in {path}; IDs must be unique.')
        invalid = set(seq) - VALID_AA
        if invalid:
            raise ValueError(f'Non-protein or aligned sequence for {original_id} in {path}: {invalid}. '
                             'Supply unaligned amino-acid FASTA.')
        if '*' in seq.rstrip('*'):
            stop_records.append(original_id)
            if qc_rows is not None:
                qc_rows.append(dict(role=role,input_file=str(path.resolve()),protein_id=original_id,
                    internal_stop_count=seq.rstrip('*').count('*'),action='retained_in_normalized_sequence'))
        if not seq.rstrip('*'):
            raise ValueError(f'Sequence {original_id} in {path} contains no amino-acid residues.')
        used.add(original_id)
        records.append(Protein(original_id, header, seq.rstrip('*'), f'{prefix}{len(records)+1:07d}'))

    with open_text(path) as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            if line.startswith('>'):
                append_one(desc, chunks)
                desc, chunks = line[1:].strip(), []
            elif desc is None:
                raise ValueError(f'{path} is not FASTA: sequence before first header.')
            else:
                chunks.append(line)
        append_one(desc, chunks)
    if not records:
        raise ValueError(f'No proteins found in {path}')
    if stop_records:
        examples=', '.join(stop_records[:5])
        remaining=f' (+{len(stop_records)-5} more)' if len(stop_records)>5 else ''
        print(f'WARNING: {path}: {len(stop_records)} proteins contain internal stop markers; '
              f'sequences retained. Examples: {examples}{remaining}. '
              'Successful runs list every affected protein in sequence_qc.tsv.',file=sys.stderr)
    return records


def fasta_files(path: Path) -> list[Path]:
    if path.is_file():
        if not any(path.name.lower().endswith(e) or path.name.lower().endswith(e+'.gz') for e in EXTS):
            raise ValueError(f'Unsupported protein FASTA format: {path}')
        return [path]
    if not path.is_dir():
        raise FileNotFoundError(path)
    found = []
    for p in sorted(path.rglob("*")):
        name = p.name.lower()
        if p.is_file() and any(name.endswith(e) or name.endswith(e + '.gz') for e in EXTS):
            found.append(p)
    if not found:
        raise ValueError(f'No FASTA files in {path}; expected .faa, .fa, .fasta, etc. (optionally .gz)')
    return found


def species_label(path: Path) -> str:
    name = path.name
    if name.lower().endswith('.gz'):
        name = name[:-3]
    for ext in EXTS:
        if name.lower().endswith(ext):
            name = name[:-len(ext)]
            break
    name = re.sub(r'[^A-Za-z0-9_.-]', '_', name)
    if not name:
        raise ValueError(f'Cannot construct species name from {path.name}')
    return name


def write_normalized(path: Path, records: list[Protein]):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', encoding='utf-8') as f:
        for r in records:
            f.write(f'>{r.internal_id}\n')
            for i in range(0, len(r.seq), 80):
                f.write(r.seq[i:i+80] + '\n')


def write_tsv(path: Path, names, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(names), delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def append_fasta(path: Path, items):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', encoding='utf-8') as f:
        for identifier, metadata, sequence in items:
            f.write(f'>{identifier} {metadata}\n')
            for i in range(0, len(sequence), 80):
                f.write(sequence[i:i+80] + '\n')


def run(cmd: list[str], log, stdout_path: Path | None = None):
    log.write(' '.join(str(x) for x in cmd) + '\n')
    log.flush()
    if stdout_path is None:
        result = subprocess.run(cmd, capture_output=True, text=True)
    else:
        with stdout_path.open('w', encoding='utf-8') as handle:
            result = subprocess.run(cmd, stdout=handle, stderr=subprocess.PIPE, text=True)
    if result.returncode:
        raise RuntimeError('Command failed:\n' + ' '.join(cmd) + '\n' + (result.stderr or ''))


def parse_blast(path: Path):
    """BLAST outfmt 6, one HSP per query-subject pair."""
    grouped = defaultdict(list)
    with path.open(encoding='utf-8') as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            cols = line.split('\t')
            if len(cols) != len(FIELDS):
                raise ValueError(f'Unexpected BLAST output ({len(cols)} columns) in {path}')
            x = dict(zip(FIELDS, cols))
            for key in ('pident', 'evalue', 'bitscore'):
                x[key] = float(x[key])
            for key in ('length', 'qlen', 'slen', 'qstart', 'qend', 'sstart', 'send'):
                x[key] = int(x[key])
            x['qcov'] = 100 * (abs(x['qend'] - x['qstart']) + 1) / x['qlen']
            x['scov'] = 100 * (abs(x['send'] - x['sstart']) + 1) / x['slen']
            grouped[x['qseqid']].append(x)
    for rows in grouped.values():
        rows.sort(key=lambda x: (-x['bitscore'], x['evalue'], -x['qcov'],
                                 -x['scov'], -x['pident'], x['sseqid']))
    return grouped


def passes(hit, args):
    return (hit['evalue'] <= args.evalue and hit['pident'] >= args.min_identity
            and hit['qcov'] >= args.min_qcov and hit['scov'] >= args.min_scov)


def bit_gap(hits):
    return None if len(hits) < 2 else hits[0]['bitscore'] - hits[1]['bitscore']


def ambiguous(hits, args):
    gap = bit_gap(hits)
    return gap is not None and gap <= args.min_bit_gap + 1e-8


def fmt(value):
    if value is None:
        return ''
    return f'{value:.4g}' if isinstance(value, float) else str(value)


def metrics(hit, prefix, id_map):
    if hit is None:
        return {prefix + k: '' for k in ('id', 'bitscore', 'evalue', 'pident', 'qcov', 'scov', 'length')}
    return {
        prefix + 'id': id_map[hit['sseqid']].original_id,
        prefix + 'bitscore': fmt(hit['bitscore']),
        prefix + 'evalue': fmt(hit['evalue']),
        prefix + 'pident': fmt(hit['pident']),
        prefix + 'qcov': fmt(hit['qcov']),
        prefix + 'scov': fmt(hit['scov']),
        prefix + 'length': hit['length'],
    }


def background_mapping(refs, background, augment_missing=False):
    """Map exact sequences; optionally add missing curated references explicitly.

    All existing competitors remain unchanged. Multiple exact matches still fail.
    Invalid mappings never leave a partially augmented background.
    """
    by_id={p.original_id:p for p in background}
    by_seq=defaultdict(list)
    for p in background:by_seq[p.seq].append(p)
    mapping,methods,additions,problems={},{},[],[]
    internal_ids={p.internal_id for p in background};original_ids=set(by_id)
    for index,r in enumerate(refs,1):
        same_id=by_id.get(r.original_id)
        alias_matches={token:by_id[token] for token in r.original_id.split('|')
                       if token in by_id and by_id[token].seq==r.seq}
        token_match=next(iter(alias_matches.values())) if len(alias_matches)==1 else None
        if same_id and same_id.seq==r.seq:
            match,method=same_id,'id_and_sequence'
        elif token_match is not None:
            match,method=token_match,'identifier_token_and_sequence'
        else:
            options=by_seq.get(r.seq,[])
            if len(options)==1:
                match,method=options[0],'identical_sequence'
            elif not options and augment_missing:
                base='curated_added|'+r.original_id;identifier=base;number=1
                while identifier in original_ids:
                    number+=1;identifier=f'{base}|{number}'
                internal=f'A{index:07d}';number=1
                while internal in internal_ids:
                    internal=f'A{index:07d}_{number}';number+=1
                match=Protein(identifier,r.description+' [curated sequence added to reverse background]',r.seq,internal)
                additions.append(match);original_ids.add(identifier);internal_ids.add(internal)
                method='curated_reference_added'
            else:
                problems.append(f'{r.original_id!r}: found {len(options)} exact background matches')
                continue
        mapping[r.internal_id]=match.internal_id;methods[r.internal_id]=method
    if problems:
        raise ValueError('Curated reference/background mapping failed:\n  - '+'\n  - '.join(problems)+
            '\nUse matching source-proteome releases/isoforms, or explicitly choose '
            '--augment-reference-background to add missing curated sequences while retaining every '
            'background competitor. Multiple exact matches still require reconciliation; '
            'no approximate sequence mapping is used.')
    background.extend(additions)
    return mapping,methods


ROOT = Path(__file__).resolve().parent
PROJECT_ROOT = ROOT.parent/'integrin_adhesome_analysis'
sys.path.insert(0, str(PROJECT_ROOT/'resource_library/scripts'))
from input_dependencies import InputRegistry
SOURCE_SPECIES_ALIASES = {'candidate':'homo_sapiens'}  # Explicit user-confirmed source label.
AUDIT_FIELDS = ['reference_set','species','reference_id','reference_description','reference_file',
    'reference_proteome','source_assignment','mapping_method','mapped_background_id','best_target_id',
    'forward_bitscore','forward_evalue','forward_pident','forward_qcov','forward_scov','forward_length',
    'forward_second_bit_gap','reverse_best_id','reverse_bitscore','reverse_evalue','reverse_pident',
    'reverse_qcov','reverse_scov','reverse_length','reverse_second_bit_gap','status','reverse_database_mode']
HIT_FIELDS = ['reference_set','species','reference_id','target_id','bitscore','evalue','pident',
    'qcov','scov','length','qlen','slen','qstart','qend','sstart','send','passes_filters','rank']
MAP_FIELDS = ['reference_set','reference_file','reference_id','reference_internal_id','reference_proteome',
    'reference_proteome_file','source_assignment','background_id','background_internal_id','mapping_method','added_to_background']
QC_FIELDS = ['role','input_file','protein_id','internal_stop_count','action']


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def dump(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + '\n', encoding='utf-8')


def default_directory(names):
    found = [ROOT/name for name in names if (ROOT/name).is_dir()]
    if len(found) > 1:
        raise ValueError('Multiple input-directory aliases exist; choose one explicitly: ' + str(found))
    return found[0] if found else ROOT/names[0]


def collection(path):
    """Recursive deterministic names retain relative directories, not just basenames."""
    path = Path(path).resolve()
    base = path.parent if path.is_file() else path
    files = fasta_files(path)
    names = {}
    for file in files:
        relative = file.relative_to(base)
        label = '__'.join([*relative.parts[:-1], species_label(file)])
        label = re.sub(r'[^A-Za-z0-9_.-]', '_', label)
        if label in ('', '.', '..') or label.casefold() in {name.casefold() for name in names}:
            raise ValueError(f'Input paths produce the same output label {label!r}: {names[label]} and {file}')
        names[label] = file
    return base, names


def name_key(value):
    return re.sub(r'[^a-z0-9]', '', str(value).casefold())


def source_rules(path):
    if path is None:
        return []
    with path.open(encoding='utf-8-sig', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        required = {'reference_file','reference_id','reference_proteome'}
        if not required <= set(reader.fieldnames or []):
            raise ValueError('Source map needs reference_file, reference_id and reference_proteome columns')
        rules = list(reader)
    seen = {}
    for row in rules:
        key = (row['reference_file'], row['reference_id'])
        if key in seen and seen[key] != row['reference_proteome']:
            raise ValueError('Conflicting source-map rows: ' + str(key))
        seen[key] = row['reference_proteome']
    return rules


def assign_source(query, backgrounds, reference_base, background_base, rules, sequence_sources, single_background):
    relative = query['file'].relative_to(reference_base).as_posix()
    applicable = [row for row in rules if row['reference_file'].replace('\\','/') == relative
                  and row['reference_id'] in ('*', query['protein'].original_id)]
    exact = [row for row in applicable if row['reference_id'] != '*']
    applicable = exact or applicable
    def aliases(label, bg):
        return {name_key(value) for value in (label, species_label(bg['file']),
                 *bg['file'].relative_to(background_base).parts[:-1])}
    if applicable:
        requested = {row['reference_proteome'].replace('\\','/') for row in applicable}
        requested_labels = {'__'.join([*Path(value).parts[:-1], species_label(Path(value))]) for value in requested}
        choices = [label for label, bg in backgrounds.items() if label in requested | requested_labels
                   or bg['file'].relative_to(background_base).as_posix() in requested]
        method = 'explicit_source_map'
    else:
        original_prefix = query['protein'].original_id.split('|', 1)[0]
        prefix = SOURCE_SPECIES_ALIASES.get(original_prefix.casefold(), original_prefix)
        choices = [label for label, bg in backgrounds.items() if name_key(prefix) in aliases(label, bg)]
        method = ('user_confirmed_header_alias:'+original_prefix+'='+prefix
                  if prefix != original_prefix else 'reference_header_species')
        if not choices and re.fullmatch(r'[A-Za-z]+_[A-Za-z]+', prefix):
            raise ValueError(f'{relative}: {query["protein"].original_id}: no matching proteome for '
                             f'explicit source species {prefix!r}; supply it or an explicit --source-map.')
        if not choices:
            parts = Path(relative).parts[:-1]
            stem = re.sub(r'(?:_references?|_refs)$', '', species_label(query['file']), flags=re.I)
            labels = {name_key(part) for part in [*parts, stem]}
            choices = [label for label, bg in backgrounds.items() if labels & aliases(label, bg)]
            method = 'reference_path_species'
        if not choices:
            choices = sorted(sequence_sources.get(query['protein'].seq, set()))
            method = 'unique_exact_sequence_source'
    if len(choices) != 1:
        raise ValueError(f'{relative}: {query["protein"].original_id}: expected one corresponding source '
                         f'proteome, found {choices}. Use species-prefixed headers, mirrored species folders, '
                         'or --source-map. No approximate source assignment is performed.')
    return choices[0], method


def prepare_inputs(args):
    reference_base, reference_files = collection(args.references)
    target_base, target_files = collection(args.proteomes)
    requested = set(getattr(args, 'reference_set', None) or [])
    if requested:
        if requested-set(reference_files):
            raise ValueError('Unidentified reference sets: '+str(sorted(requested-set(reference_files))))
        reference_files = {k:v for k,v in reference_files.items() if k in requested}
    qc, queries, groups = [], [], {}
    for label, file in reference_files.items():
        groups[label] = []
        for protein in read_fasta(file, 'R', qc, role='curated_reference'):
            protein = Protein(protein.original_id, protein.description, protein.seq, f'R{len(queries)+1:07d}')
            query = dict(reference_set=label, file=file, protein=protein)
            groups[label].append(query)
            queries.append(query)
    backgrounds, sequence_sources, problems = {}, defaultdict(set), []
    rules = source_rules(args.source_map)
    if args.reference_only:
        background_base = reference_base
        for label, group in groups.items():
            key = 'curated__' + label
            records = [Protein(q['protein'].original_id, q['protein'].description,
                               q['protein'].seq, 'B' + q['protein'].internal_id[1:]) for q in group]
            backgrounds[key] = dict(file=reference_files[label], proteins=records, reference_only=True)
            for query in group:
                query.update(background=key, source_assignment='explicit_curated_reference_only')
    else:
        background_base, background_files = collection(args.reference_background)
        selected = set(getattr(args, 'reference_proteome', None) or [])
        if selected:
            if selected-set(background_files):
                raise ValueError('Unidentified reference proteomes: '+str(sorted(selected-set(background_files))))
            background_files = {k:v for k,v in background_files.items() if k in selected}
        # Identify known source species using filenames/headers/maps before
        # reading proteomes. Only unlabelled queries need exact-sequence source
        # discovery across the supplied candidate collection.
        stubs = {label:dict(file=file) for label,file in background_files.items()}
        needed, unlabelled = set(), False
        for query in queries:
            try:
                label, method = assign_source(query, stubs, reference_base, background_base, rules, {}, False)
                needed.add(label)
            except ValueError:
                prefix = query['protein'].original_id.split('|',1)[0]
                prefix = SOURCE_SPECIES_ALIASES.get(prefix.casefold(), prefix)
                relative = query['file'].relative_to(reference_base).as_posix()
                mapped = any(r['reference_file'].replace('\\','/') == relative and
                             r['reference_id'] in ('*',query['protein'].original_id) for r in rules)
                if not mapped and not re.fullmatch(r'[A-Za-z]+_[A-Za-z]+', prefix):
                    unlabelled = True
        for i, (label, file) in enumerate(background_files.items(), 1):
            if not unlabelled and label not in needed:
                continue
            records = read_fasta(file, f'B{i:03d}_', qc, role='reference_background')
            backgrounds[label] = dict(file=file, proteins=records, reference_only=False)
            for protein in records:
                sequence_sources[protein.seq].add(label)
        for query in queries:
            try:
                label, method = assign_source(query, backgrounds, reference_base, background_base,
                                              rules, sequence_sources, len(backgrounds) == 1)
                query.update(background=label, source_assignment=method)
            except ValueError as exc:
                problems.append(str(exc))
    mappings = []
    for label, bg in backgrounds.items():
        assigned = [q for q in queries if q.get('background') == label]
        if not assigned:
            continue
        try:
            mapping, methods = background_mapping([q['protein'] for q in assigned], bg['proteins'],
                                                   args.augment_reference_background)
        except ValueError as exc:
            problems.append(label + ': ' + str(exc))
            continue
        by_id = {p.internal_id:p for p in bg['proteins']}
        for query in assigned:
            protein = query['protein']
            query.update(mapped=mapping[protein.internal_id], mapping_method=methods[protein.internal_id])
            mappings.append(dict(reference_set=query['reference_set'], reference_file=str(query['file']),
                reference_id=protein.original_id, reference_internal_id=protein.internal_id,
                reference_proteome=label, reference_proteome_file=str(bg['file']),
                source_assignment=query['source_assignment'], background_id=by_id[query['mapped']].original_id,
                background_internal_id=query['mapped'], mapping_method=query['mapping_method'],
                added_to_background=int(query['mapping_method']=='curated_reference_added')))
    targets = []
    for label, file in target_files.items():
        records = read_fasta(file, 'P', qc, role='target')
        targets.append(dict(species=label, file=str(file), proteins=len(records)))
    inputs = sorted(set(reference_files.values()) | set(target_files.values()) |
                    ({bg['file'] for bg in backgrounds.values()} if not args.reference_only else set()))
    if args.source_map is not None:
        inputs.append(args.source_map.resolve())
    return dict(queries=queries, groups=groups, backgrounds=backgrounds, target_files=target_files,
                targets=targets, qc=qc, mappings=mappings, problems=problems,
                input_hashes={str(file):digest(file) for file in inputs})


def input_report(prepared):
    return dict(status='inputs_invalid' if prepared['problems'] else 'inputs_valid',
        analysis_launched=False, files_written=False,
        reference_sets=[dict(reference_set=label, queries=len(group), file=str(group[0]['file']))
                        for label, group in prepared['groups'].items()],
        curated_references=len(prepared['queries']), target_proteomes=prepared['targets'],
        reference_proteomes=[dict(reference_proteome=label, file=str(bg['file']),
            proteins=len(bg['proteins'])) for label, bg in prepared['backgrounds'].items()],
        curated_sequences_added=sum(row['added_to_background'] for row in prepared['mappings']),
        reference_mappings=prepared['mappings'], proteins_with_internal_stop_markers=len(prepared['qc']),
        errors=prepared['problems'], input_decisions=[dict(action='unidentified_input',stage='homology.source_assignment',
            mandatory=True,decision='stop_stage',reason=problem,
            remedy='Provide the matching raw reference proteome or an explicit --source-map; no unrelated source is substituted.')
            for problem in prepared['problems']])


def cached_command(command, log, checkpoint, inputs, output=None, database=None, resume=False):
    signature = dict(command=[str(x) for x in command], executable_sha256=digest(shutil.which(command[0])),
                     input_sha256={str(path):digest(path) for path in inputs})
    if checkpoint.exists():
        saved = json.loads(checkpoint.read_text())
        if not resume or saved['signature'] != signature:
            raise ValueError('Cached command settings/inputs/tool differ: ' + str(checkpoint))
        if saved['status'] == 'complete':
            if any(not Path(path).is_file() or digest(path)!=value for path,value in saved['outputs'].items()):
                raise ValueError('Completed command outputs changed: ' + str(checkpoint))
            print('REUSE ' + checkpoint.stem, file=sys.stderr)
            log.write('REUSE '+checkpoint.stem+'\n')
            log.flush()
            return
    dump(checkpoint, dict(status='running', signature=signature))
    run(command, log, output)
    if database is not None:
        outputs = sorted(database.parent.glob(database.name + '.*'))
        if not all(any(path.name.endswith('.'+ext) for path in outputs) for ext in ('pin','phr','psq')):
            raise ValueError('Incomplete protein BLAST database: ' + str(database))
    else:
        outputs = [output]
        if output is None or not output.is_file():
            raise ValueError('Missing BLAST search output: ' + str(output))
    dump(checkpoint, dict(status='complete', signature=signature,
                          outputs={str(path):digest(path) for path in outputs}))


def export_candidates(path, supports, targets, queries, group=None):
    items = []
    for (species, tid), qids in sorted(supports.items()):
        qids = sorted(qid for qid in qids if group is None or queries[qid]['reference_set']==group)
        if not qids:
            continue
        protein = targets[(species, tid)]
        sets = sorted({queries[qid]['reference_set'] for qid in qids})
        items.append((species+'|'+protein.original_id,
                      'reference_sets='+','.join(sets)+' matching_reference_queries='+str(len(qids)), protein.seq))
    append_fasta(path, items)
    return len(items)


def execute(prepared, args, out):
    work = out/'work'
    work.mkdir(parents=True, exist_ok=True)
    query_by_id = {q['protein'].internal_id:q for q in prepared['queries']}
    curated = work/'curated_references.faa'
    write_normalized(curated, [q['protein'] for q in prepared['queries']])
    write_tsv(out/'reference_mapping.tsv', MAP_FIELDS, prepared['mappings'])
    used_backgrounds = sorted({q['background'] for q in prepared['queries']})
    background_ids = {label:{p.internal_id:p for p in bg['proteins']}
                      for label,bg in prepared['backgrounds'].items()}
    augmented_backgrounds = {q['background'] for q in prepared['queries']
                             if q['mapping_method']=='curated_reference_added'}
    audit, hits_table, summary = [], [], []
    all_supports, rbh_supports = defaultdict(set), defaultdict(set)
    # Retain sequences only for recovered candidates, not all target proteomes.
    candidate_proteins = {}
    log_mode = 'a' if args.resume else 'w'
    with (out/'commands.log').open(log_mode, encoding='utf-8') as log:
        databases = {}
        for label in used_backgrounds:
            print('[reference database] '+label, file=sys.stderr, flush=True)
            bg = prepared['backgrounds'][label]
            folder = work/'reference_databases'/label
            folder.mkdir(parents=True, exist_ok=True)
            normalized = folder/'proteins.faa'
            write_normalized(normalized, bg['proteins'])
            database = folder/'proteome'
            cached_command(['makeblastdb','-in',str(normalized),'-dbtype','prot','-parse_seqids','-out',str(database)],
                log, work/'checkpoints'/('reference_db__'+label+'.json'), [normalized], database=database,
                resume=args.resume)
            databases[label] = database
        for index, (species, file) in enumerate(prepared['target_files'].items(),1):
            print(f'[target {index}/{len(prepared["target_files"])}] {species}: forward/reverse searches',
                  file=sys.stderr, flush=True)
            proteins = read_fasta(file,'P')
            by_id = {p.internal_id:p for p in proteins}
            normalized = work/(species+'.faa')
            write_normalized(normalized, proteins)
            database = work/(species+'.db')
            cached_command(['makeblastdb','-in',str(normalized),'-dbtype','prot','-parse_seqids','-out',str(database)],
                log, work/'checkpoints'/('target_db__'+species+'.json'), [normalized], database=database,
                resume=args.resume)
            forward_file = work/(species+'.forward.tsv')
            search_options = ['-evalue',str(args.evalue),'-max_target_seqs',str(args.max_targets),
                              '-max_hsps','1','-num_threads',str(args.threads),'-outfmt',OUTFMT]
            cached_command(['blastp','-query',str(curated),'-db',str(database),*search_options],
                log, work/'checkpoints'/('forward__'+species+'.json'), [curated,normalized],
                output=forward_file, resume=args.resume)
            forward = parse_blast(forward_file)
            for qid, hits in forward.items():
                if qid not in query_by_id:
                    raise ValueError('Unknown reference query ID in BLAST output: '+qid)
                query = query_by_id[qid]
                for rank, hit in enumerate(hits,1):
                    tid = hit['sseqid']
                    protein = by_id[tid]
                    hits_table.append(dict(reference_set=query['reference_set'],species=species,
                        reference_id=query['protein'].original_id,target_id=protein.original_id,
                        **{k:fmt(hit[k]) for k in ('bitscore','evalue','pident','qcov','scov','length','qlen','slen',
                                                  'qstart','qend','sstart','send')},
                        passes_filters=int(passes(hit,args)),rank=rank))
                    if passes(hit,args):
                        all_supports[(species,tid)].add(qid)
                        candidate_proteins[(species,tid)] = protein
            reverse = {}
            for label in used_backgrounds:
                candidates = sorted({hits[0]['sseqid'] for qid,hits in forward.items()
                    if hits and query_by_id[qid]['background']==label and passes(hits[0],args)})
                reverse[label] = {}
                if not candidates:
                    continue
                print(f'[reverse] {species} -> {label}: {len(candidates)} candidate proteins',
                      file=sys.stderr, flush=True)
                folder = work/'reverse'/species/label
                folder.mkdir(parents=True, exist_ok=True)
                query_file, reverse_file = folder/'query.faa', folder/'hits.reverse.tsv'
                write_normalized(query_file,[by_id[tid] for tid in candidates])
                bg_input = work/'reference_databases'/label/'proteins.faa'
                cached_command(['blastp','-query',str(query_file),'-db',str(databases[label]),*search_options],
                    log, work/'checkpoints'/('reverse__'+species+'__'+label+'.json'),[query_file,bg_input],
                    output=reverse_file,resume=args.resume)
                reverse[label] = parse_blast(reverse_file)
            counts = defaultdict(int)
            for qid, query in query_by_id.items():
                hits = forward.get(qid,[])
                top = hits[0] if hits else None
                label = query['background']
                bg = prepared['backgrounds'][label]
                bg_ids = background_ids[label]
                rev = reverse[label].get(top['sseqid'],[]) if top else []
                best = rev[0] if rev else None
                if not top:status='no_forward_hit'
                elif not passes(top,args):status='forward_below_filters'
                elif ambiguous(hits,args):status='ambiguous_forward_best'
                elif not best:status='no_reverse_hit'
                elif not passes(best,args):status='reverse_below_filters'
                elif ambiguous(rev,args):status='ambiguous_reverse_best'
                elif best['sseqid']!=query['mapped']:status='reverse_matches_other_reference'
                elif args.reference_only:status='provisional_reference_only_RBH'
                elif query['mapping_method']=='curated_reference_added':status='RBH_candidate_augmented_reference'
                else:status='RBH_candidate'
                mode = ('curated_sequences_only' if args.reference_only else
                        'full_reference_proteome_plus_curated' if label in augmented_backgrounds
                        else 'full_reference_proteome')
                row = dict(reference_set=query['reference_set'],species=species,
                    reference_id=query['protein'].original_id,reference_description=query['protein'].description,
                    reference_file=str(query['file']),reference_proteome=label,
                    source_assignment=query['source_assignment'],mapping_method=query['mapping_method'],
                    mapped_background_id=bg_ids[query['mapped']].original_id,
                    best_target_id=by_id[top['sseqid']].original_id if top else '',
                    forward_second_bit_gap=fmt(bit_gap(hits)),reverse_second_bit_gap=fmt(bit_gap(rev)),
                    reverse_best_id=bg_ids[best['sseqid']].original_id if best else '',
                    status=status,reverse_database_mode=mode)
                row.update(metrics(top,'forward_',by_id))
                row.update(metrics(best,'reverse_',bg_ids))
                audit.append(row)
                counts[status] += 1
                if status in RBH_STATUSES:
                    rbh_supports[(species,top['sseqid'])].add(qid)
            summary.append(dict(species=species,proteins_in_proteome=len(proteins),reference_queries=len(query_by_id),
                references_with_forward_hits=len(forward),unique_forward_qualifying_proteins=sum(
                    key[0]==species for key in all_supports),unique_rbh_candidates=sum(key[0]==species for key in rbh_supports),
                accepted_reference_target_pairs=sum(counts[name] for name in RBH_STATUSES),
                ambiguous_forward=counts['ambiguous_forward_best'],ambiguous_reverse=counts['ambiguous_reverse_best']))
    write_tsv(out/'homology_audit.tsv',AUDIT_FIELDS,audit)
    write_tsv(out/'all_forward_hits.tsv',HIT_FIELDS,hits_table)
    write_tsv(out/'targets_summary.tsv',list(summary[0]),summary)
    export_candidates(out/'all_forward_qualifying.faa',all_supports,candidate_proteins,query_by_id)
    export_candidates(out/'rbh_candidates.faa',rbh_supports,candidate_proteins,query_by_id)
    for group in prepared['groups']:
        folder = out/'reference_sets'/group
        write_tsv(folder/'homology_audit.tsv',AUDIT_FIELDS,[row for row in audit if row['reference_set']==group])
        export_candidates(folder/'all_forward_qualifying.faa',all_supports,candidate_proteins,query_by_id,group)
        export_candidates(folder/'rbh_candidates.faa',rbh_supports,candidate_proteins,query_by_id,group)
    write_tsv(out/'sequence_qc.tsv',QC_FIELDS,prepared['qc'])
    return dict(reference_target_pairs=len(audit),unique_forward_qualifying_proteins=len(all_supports),
                unique_rbh_candidates=len(rbh_supports),reference_sets=len(prepared['groups']),
                target_proteomes=len(prepared['target_files']))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--reference-dir','--references',dest='references',type=Path,
                        help='Curated FASTA file or recursive directory; default reference_fasta (also refrence_fasta)')
    parser.add_argument('--target-dir','--proteomes',dest='proteomes',type=Path,
                        help='Target FASTA file or recursive directory; default target_proteomes; any target count')
    parser.add_argument('--background-dir','--reference-background',dest='reference_background',type=Path,
                        help='Complete source proteomes, recursively read; default reference_proteomes (also reference_proteome)')
    parser.add_argument('--source-map',type=Path,help='Optional TSV matching reference files/IDs to background files')
    parser.add_argument('--reference-set',action='append',help='Select an exact recursive reference-set label; repeat for several families')
    parser.add_argument('--reference-proteome',action='append',help='Restrict possible source-proteome labels; repeat for several')
    parser.add_argument('--dry-run',action='store_true',help='Dataset-free structural deployment validation')
    parser.add_argument('--clean-state',action='store_true',help='Treat biological inputs as absent during --dry-run')
    parser.add_argument('--report-dir',type=Path,help='Optional small deployment-validation report directory')
    parser.add_argument('--reference-only',action='store_true',help='Explicit provisional reverse searches against each curated reference set only')
    parser.add_argument('--augment-reference-background',action='store_true',
                        help='Explicitly add missing curated sequences to their matched source proteomes; keep all competitors and label additions')
    parser.add_argument('--check-inputs',action='store_true',help='Read-only FASTA/source validation; no BLAST or output writes')
    parser.add_argument('--resume',action='store_true',help='Reuse verified completed database/search checkpoints')
    parser.add_argument('--outdir',type=Path,default=ROOT/'results/homology_run')
    parser.add_argument('--threads',type=int,default=4)
    parser.add_argument('--evalue',type=float,default=1e-5)
    parser.add_argument('--min-identity',type=float,default=25)
    parser.add_argument('--min-qcov',type=float,default=50)
    parser.add_argument('--min-scov',type=float,default=50)
    parser.add_argument('--min-bit-gap',type=float,default=0)
    parser.add_argument('--max-targets',type=int,default=500)
    args = parser.parse_args(argv)
    if args.reference_only and args.augment_reference_background:
        parser.error('--augment-reference-background requires a full background, not --reference-only')
    if args.reference_only and args.reference_background is not None:
        parser.error('Choose --reference-only or a background directory, not both')
    if not all(math.isfinite(x) for x in (args.evalue,args.min_identity,args.min_qcov,args.min_scov,args.min_bit_gap)):
        parser.error('BLAST cutoffs must be finite numbers')
    if args.threads<1 or args.max_targets<2 or args.min_bit_gap<0 or args.evalue<=0:
        parser.error('Require threads>=1, max-targets>=2, bit-gap>=0 and evalue>0')
    if any(not 0<=x<=100 for x in (args.min_identity,args.min_qcov,args.min_scov)):
        parser.error('Identity/coverage percentages must be within [0,100]')
    if args.dry_run:
        from deployment_validation import main as validate_deployment
        flags=['--dry-run','--workflow','homology']
        if args.reference_only:flags.append('--reference-only')
        for name,flag in [('references','--curated-references'),('proteomes','--target-raw'),('reference_background','--reference-raw')]:
            if getattr(args,name):flags += [flag,str(getattr(args,name))]
        if args.clean_state: flags.append('--clean-state')
        if args.report_dir: flags += ['--report-dir',str(args.report_dir)]
        return validate_deployment(flags)
    args.references = (args.references or default_directory(['reference_fasta','refrence_fasta'])).resolve()
    registry = InputRegistry(PROJECT_ROOT)
    args.proteomes = registry.require('target_proteomes','homology.target_discovery',args.proteomes,allow_file=True)
    registry.require('curated_references','homology.reference_queries',args.references,allow_file=True)
    if not args.reference_only:
        args.reference_background = registry.require('reference_proteomes','homology.reference_sources',args.reference_background,allow_file=True)
    if args.source_map is None and (ROOT/'config/reference_sources.tsv').is_file():
        args.source_map = ROOT/'config/reference_sources.tsv'
    if args.source_map is not None:
        args.source_map = args.source_map.resolve()
    out = args.outdir.resolve()
    for source in (args.references,args.proteomes,args.reference_background):
        if source is not None and source.is_dir() and (out==source or source in out.parents):
            raise ValueError('--outdir must be outside recursively scanned input directories')
    prepared = prepare_inputs(args)
    report = input_report(prepared)
    if args.check_inputs:
        print(json.dumps(report,indent=2))
        return 1 if prepared['problems'] else 0
    if prepared['problems']:
        raise ValueError('Input validation failed:\n'+'\n'.join(prepared['problems']))
    for program in ('blastp','makeblastdb'):
        if shutil.which(program) is None:
            parser.error(program+' is not on PATH; activate a BLAST+ environment')
    signature = dict(script_sha256=digest(__file__),input_sha256=prepared['input_hashes'],
        reference_only=args.reference_only,augment_reference_background=args.augment_reference_background,
        threads=args.threads,evalue=args.evalue,min_identity=args.min_identity,min_qcov=args.min_qcov,
        min_scov=args.min_scov,min_bit_gap=args.min_bit_gap,max_targets=args.max_targets,
        executables={name:dict(path=shutil.which(name),sha256=digest(shutil.which(name)))
                     for name in ('blastp','makeblastdb')})
    config = out/'run_config.json'
    if args.resume:
        if not config.is_file() or json.loads(config.read_text())!=signature:
            raise ValueError('--resume requires the saved inputs/code/settings/tools; otherwise choose a new --outdir')
        complete = out/'run_manifest.json'
        if complete.exists():
            saved = json.loads(complete.read_text())
            if any(not (out/name).is_file() or digest(out/name)!=value for name,value in saved['output_sha256'].items()):
                raise ValueError('Completed outputs changed; refusing resume')
            print('REUSE completed homology run: '+str(out))
            return 0
    elif out.exists() and any(out.iterdir()):
        raise ValueError('Output is not empty; choose a new --outdir or use --resume for this pipeline')
    out.mkdir(parents=True,exist_ok=True)
    if not config.exists():
        dump(config,signature)
    dump(out/'run_state.json',dict(status='running',updated_utc=datetime.now(timezone.utc).isoformat()))
    try:
        result = execute(prepared,args,out)
        output_hashes={p.relative_to(out).as_posix():digest(p) for p in out.rglob('*')
            if p.is_file() and p.name not in {'run_manifest.json','run_state.json','commands.log'} and 'work' not in p.relative_to(out).parts}
        dump(out/'run_manifest.json',dict(status='complete',analysis_launched=True,summary=result,input_report=report,
            background_policy='per_reference_source_proteome',output_sha256=output_hashes,
            completed_utc=datetime.now(timezone.utc).isoformat()))
        dump(out/'run_state.json',dict(status='complete',summary=result))
    except (OSError,ValueError,RuntimeError,subprocess.CalledProcessError) as exc:
        dump(out/'run_state.json',dict(status='failed',error=str(exc)))
        raise
    print('Homologue discovery complete: '+str(out)+'\n'+json.dumps(result,indent=2))
    return 0


if __name__ == '__main__':
    try:
        sys.exit(main())
    except (OSError,ValueError,RuntimeError,subprocess.CalledProcessError) as exc:
        print('ERROR: '+str(exc),file=sys.stderr)
        sys.exit(1)


