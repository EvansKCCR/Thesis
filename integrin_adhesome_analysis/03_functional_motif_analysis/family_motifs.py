#!/usr/bin/env python3
"""Family-level protein motif inspection and exploratory discovery (Python >=3.10)."""
import argparse
import csv
import hashlib
import json
import math
import platform
import re
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

AA = set('ACDEFGHIKLMNPQRSTVWY')
HERE = Path(__file__).resolve().parent


def read_fasta(path):
    records, seen = [], set()
    header, chunks = None, []

    def flush():
        if header is None:
            return
        ident = header.split()[0]
        if ident in seen:
            raise ValueError(f'{path}: duplicate identifier {ident}')
        seen.add(ident)
        raw = ''.join(chunks).upper()
        seq = raw[:-1] if raw.endswith('*') else raw
        if not seq or set(seq) - (AA | set('XBZJUO')):
            raise ValueError(f'{path}: empty/invalid protein sequence for {ident}')
        records.append(dict(id=ident, header=header, sequence=seq,
                            terminal_stop=raw.endswith('*')))

    for line in Path(path).read_text(encoding='utf-8-sig').splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            flush()
            header, chunks = line[1:].strip(), []
            if not header:
                raise ValueError(f'{path}: empty FASTA header')
        else:
            if header is None:
                raise ValueError(f'{path}: sequence before header')
            chunks.append(''.join(line.split()))
    flush()
    if not records:
        raise ValueError(f'{path}: no sequences')
    return records


def write_tsv(path, rows, fields):
    with Path(path).open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def unique(records):
    return list({r['sequence']: r for r in reversed(records)}.values())


def occurrences(pattern, sequence):
    """Overlapping matches; ambiguous amino acids cannot establish a motif."""
    for match in re.finditer('(?=(' + pattern + '))', sequence):
        peptide = match.group(1)
        if peptide and set(peptide) <= AA:
            yield match.start() + 1, match.start() + len(peptide), peptide


def load_regions(path, families):
    regions = defaultdict(list)
    if not path:
        return regions
    lengths = {(f, r['id']): len(r['sequence']) for f, rs in families.items() for r in rs}
    with Path(path).open(encoding='utf-8-sig', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        required = {'family', 'sequence_id', 'region', 'start', 'end'}
        if not required <= set(reader.fieldnames or []):
            raise ValueError('Regions TSV requires family, sequence_id, region, start, end')
        for row in reader:
            key = row['family'], row['sequence_id']
            a, b = int(row['start']), int(row['end'])
            if key not in lengths or not 1 <= a <= b <= lengths[key]:
                raise ValueError(f'Unknown sequence or invalid region coordinates: {row}')
            regions[key].append((row['region'], a, b))
    return regions


def context_state(spec, start, end, length, regions):
    context = spec.get('context', 'whole_sequence')
    if context == 'whole_sequence':
        return 'sequence_only'
    names = spec.get('regions', [])
    supplied = [(a, b) for name, a, b in regions if name in names]
    if supplied:
        return 'context_supported' if any(a <= start <= end <= b for a, b in supplied) else 'outside_context'
    if context == 'c_terminal_window':
        return 'tail_window_only' if start > max(0, length - spec.get('window', 100)) else 'outside_context'
    return 'context_unverified'


def inspect(families, library, regions):
    hits, states, prevalence = [], [], []
    for family, records in families.items():
        config = library['families'][family]
        for record in records:
            specs = [library['motifs'][m] | {'motif': m} for m in config['motifs']]
            if not specs:
                states.append(dict(family=family, sequence_id=record['id'], motif='',
                                   state='domain_architecture_required', hits=0))
            for spec in specs:
                found = []
                for a, b, peptide in occurrences(spec['pattern'], record['sequence']):
                    state = context_state(spec, a, b, len(record['sequence']), regions[family, record['id']])
                    found.append(state)
                    hits.append(dict(family=family, sequence_id=record['id'], motif=spec['motif'],
                                     start=a, end=b, peptide=peptide, context_state=state,
                                     relative_start=round(a / len(record['sequence']), 5),
                                     flank=record['sequence'][max(0, a-11):b+10],
                                     interpretation=spec['interpretation'], source=spec['source']))
                state = ('detected_in_context' if 'context_supported' in found else
                         'detected_provisional' if any(s != 'outside_context' for s in found) else
                         'outside_context_only' if found else 'not_detected')
                states.append(dict(family=family, sequence_id=record['id'], motif=spec['motif'], state=state, hits=len(found)))
        target = unique(records)
        target_seqs = {r['sequence'] for r in target}
        background = unique([r for f, rs in families.items() if f != family for r in rs if r['sequence'] not in target_seqs])
        for motif in config['motifs']:
            spec = library['motifs'][motif]
            a = sum(bool(list(occurrences(spec['pattern'], r['sequence']))) for r in target)
            b = sum(bool(list(occurrences(spec['pattern'], r['sequence']))) for r in background)
            prevalence.append(dict(family=family, motif=motif, target_positive=a, target_total=len(target),
                                   other_positive=b, other_total=len(background),
                                   target_fraction=round(a/len(target), 5),
                                   other_fraction=round(b/len(background), 5) if background else '',
                                   comparison='whole_sequence_exact_deduplicated_descriptive'))
    return hits, states, prevalence


def entropy(word):
    return -sum((n/len(word))*math.log2(n/len(word)) for n in Counter(word).values())


def inspect_domains(families, library, regions):
    rows = []
    for family, records in families.items():
        alternatives = library['families'][family].get('domain_alternatives', [])
        for r in records:
            # Unique intervals prevent duplicate imported rows inflating repeat counts.
            counts = Counter(name for name, a, b in set(regions[family, r['id']]) if name.startswith('PF'))
            missing = [{name: n-counts[name] for name, n in alt.items() if counts[name] < n} for alt in alternatives]
            rows.append(dict(family=family, sequence_id=r['id'],
                             state='not_tested_missing_domains' if not counts else
                             'domain_count_rule_supported' if any(not m for m in missing) else 'domain_count_rule_not_met',
                             observed=json.dumps(dict(counts), sort_keys=True),
                             missing_alternatives=json.dumps(missing, sort_keys=True),
                             interpretation='Count screen only; domain order, topology, homology and upstream classification rules are not re-evaluated.'))
    return rows


def discover(families, min_k, max_k, min_fraction, top):
    """Exact words are leads, not statistical evidence or independent validation."""
    indexes = {}
    for family, records in families.items():
        index = defaultdict(set)
        for r in unique(records):
            seq = r['sequence']
            for k in range(min_k, max_k+1):
                for i in range(len(seq)-k+1):
                    word = seq[i:i+k]
                    if set(word) <= AA and entropy(word) >= 1.5:
                        index[word].add(seq)
        indexes[family] = index
    rows, sites = [], []
    for family, records in families.items():
        targets = {r['sequence'] for r in records}
        background = {r['sequence'] for f, rs in families.items() if f != family for r in rs} - targets
        if len(targets) < 2:
            continue
        candidates = []
        for word, support in indexes[family].items():
            if len(support) < max(2, math.ceil(min_fraction * len(targets))):
                continue
            negatives = set().union(*(idx.get(word, set()) for f, idx in indexes.items() if f != family)) & background
            p = len(support)/len(targets)
            q = len(negatives)/len(background) if background else 0
            candidates.append(dict(family=family, motif=word, target_positive=len(support), target_total=len(targets),
                                   other_positive=len(negatives), other_total=len(background), target_fraction=round(p, 5),
                                   other_fraction=round(q, 5) if background else '', prevalence_difference=round(p-q, 5),
                                   status='exploratory_exact_word_no_significance'))
        candidates.sort(key=lambda r: (-r['prevalence_difference'], -r['target_fraction'], -len(r['motif']), r['motif']))
        chosen = []
        for row in candidates:
            word = row['motif']
            if any(word in prev['motif'] and indexes[family][word] == indexes[family][prev['motif']] for prev in chosen):
                continue
            chosen.append(row)
            if len(chosen) >= top:
                break
        rows.extend(chosen)
        for row in chosen:
            for r in records:
                for a, b, peptide in occurrences(row['motif'], r['sequence']):
                    sites.append(dict(family=family, sequence_id=r['id'], motif=peptide, start=a, end=b))
    return rows, sites


def write_fasta(path, records):
    # MEME does not support U/O; preserve length with X, retain original in inspection.
    path.write_text(''.join('>'+r['id']+'\n'+re.sub('[BZJUO]', 'X', r['sequence'])+'\n' for r in records), encoding='utf-8')


def check_meme(path):
    text = Path(path).read_text(encoding='utf-8')
    if not re.search(r'^MEME version', text, re.M) or not re.search(r'^ALPHABET\s*=\s*ACDEFGHIKLMNPQRSTVWY\s*$', text, re.M):
        raise ValueError('Known database must be MEME format with standard protein ALPHABET= ACDEFGHIKLMNPQRSTVWY')
    if not re.search(r'^MOTIF\s+\S+', text, re.M):
        raise ValueError('Known protein database contains no motifs')


def external(families, out, args):
    stages, commands = [], []
    tools = {name: shutil.which(name) for name in ('streme', 'fimo', 'tomtom')}
    versions = {}
    for name, exe in tools.items():
        if exe and args.external == 'run':
            try:
                p = subprocess.run([exe, '--version' if name == 'streme' else '-version'], capture_output=True, text=True, timeout=15)
                versions[name] = (p.stdout+p.stderr).strip()[:2000]
            except (OSError, subprocess.TimeoutExpired) as exc:
                versions[name] = str(exc)

    def stage(family, name, command, available, enabled=True, reason=''):
        command = [str(x) for x in command]
        status = 'skipped'
        if enabled:
            commands.append(dict(family=family, stage=name, argv=command))
            if args.external == 'run' and available:
                with (out/family/(name+'.log')).open('w', encoding='utf-8') as log:
                    try:
                        proc = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, timeout=args.timeout)
                        status = 'completed' if proc.returncode == 0 else 'failed'
                        reason = f'exit={proc.returncode}'
                    except (OSError, subprocess.TimeoutExpired) as exc:
                        status, reason = 'failed', str(exc)
            else:
                status = 'planned' if args.external == 'plan' else 'skipped'
                reason = 'external execution not requested' if args.external != 'run' else 'executable unavailable'
        stages.append(dict(family=family, stage=name, status=status, reason=reason))
        return status

    for family, records in families.items():
        folder = out/family
        folder.mkdir()
        positive = unique(records)
        seqs = {r['sequence'] for r in positive}
        negative = unique([r for f, rs in families.items() if f != family for r in rs if r['sequence'] not in seqs])
        write_fasta(folder/'positive.fasta', positive)
        # Globally unique IDs even when the same ID occurs in distinct family files.
        write_fasta(folder/'negative.fasta', [r | {'id': f'background_{i:06d}'} for i, r in enumerate(negative)])
        base = [tools['streme'] or 'streme', '--protein', '--p', folder/'positive.fasta', '--oc', folder/'streme',
                '--minw', args.min_width, '--maxw', args.max_width, '--nmotifs', args.top, '--seed', args.seed]
        if negative:
            base += ['--n', folder/'negative.fasta']
        status = stage(family, 'streme', base, tools['streme'], len(positive) >= args.min_streme,
                       f'fewer than {args.min_streme} distinct sequences; small-sample discovery requires explicit override')
        discovered = folder/'streme'/'streme.txt'
        ready = status == 'planned' or (status == 'completed' and discovered.exists() and re.search(r'^MOTIF\s', discovered.read_text(), re.M))
        if status == 'completed' and not ready:
            stages.append(dict(family=family, stage='discovered_motifs', status='skipped', reason='no motifs reported'))
        stage(family, 'fimo_discovered', [tools['fimo'] or 'fimo', '--oc', folder/'fimo_discovered', '--qv-thresh', '--thresh', args.qvalue,
              discovered, folder/'positive.fasta'], tools['fimo'], bool(ready), 'no completed/planned discovery motifs')
        stage(family, 'fimo_known', [tools['fimo'] or 'fimo', '--oc', folder/'fimo_known', '--qv-thresh', '--thresh', args.qvalue,
              args.known_meme or '', folder/'positive.fasta'], tools['fimo'], bool(args.known_meme), 'no protein PWM database supplied')
        stage(family, 'tomtom', [tools['tomtom'] or 'tomtom', '-oc', folder/'tomtom', '-dist', 'ed', '-thresh', args.qvalue,
              discovered, args.known_meme or ''], tools['tomtom'], bool(ready and args.known_meme), 'requires discovered motifs and protein PWM database')
    (out/'commands.json').write_text(json.dumps(commands, indent=2), encoding='utf-8')
    return stages, tools, versions


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir', type=Path, default=HERE/'family_specific_candidate_library')
    parser.add_argument('--library', type=Path, default=HERE/'motif_library.json')
    parser.add_argument('--out', type=Path, required=True, help='New output directory; existing outputs are never overwritten')
    parser.add_argument('--regions', type=Path, help='Optional 1-based inclusive domain/topology TSV')
    parser.add_argument('--known-meme', type=Path, help='Optional real protein PWM database, not converted regexes')
    parser.add_argument('--elm-dir', type=Path, help='Directory of downloaded ELM tables; optional offline annotation')
    parser.add_argument('--elm-focus', type=Path, default=HERE/'elm_resources'/'elm_family_focus.json', help='ELM family review priorities')
    parser.add_argument('--external', choices=['skip', 'plan', 'run'], default='plan')
    parser.add_argument('--min-streme', type=int, default=50)
    parser.add_argument('--min-width', type=int, default=6)
    parser.add_argument('--max-width', type=int, default=15)
    parser.add_argument('--min-k', type=int, default=4)
    parser.add_argument('--max-k', type=int, default=6)
    parser.add_argument('--min-fraction', type=float, default=0.5)
    parser.add_argument('--top', type=int, default=20)
    parser.add_argument('--qvalue', type=float, default=0.05)
    parser.add_argument('--seed', type=int, default=1729)
    parser.add_argument('--timeout', type=int, default=3600)
    args = parser.parse_args(argv)
    if not (2 <= args.min_k <= args.max_k <= 15 and 3 <= args.min_width <= args.max_width <= 30 and
            0 < args.min_fraction <= 1 and 0 < args.qvalue <= 1 and args.top > 0 and args.min_streme >= 2 and args.timeout > 0):
        parser.error('Invalid discovery width, fraction, threshold, count or timeout')
    if args.out.exists():
        parser.error('Output directory already exists; choose a new --out')
    library = json.loads(args.library.read_text(encoding='utf-8'))
    for name, spec in library['motifs'].items():
        if re.compile(spec['pattern']).match(''):
            raise ValueError(f'Empty-matching motif: {name}')
    files = sorted(p for p in args.input_dir.iterdir() if p.suffix.lower() in ('.fa', '.faa', '.fasta'))
    if not files:
        parser.error('No FASTA files in input directory')
    families = {}
    for path in files:
        if path.stem not in library['families'] or path.stem in families:
            parser.error(f'Unknown or duplicate family: {path.stem}')
        families[path.stem] = read_fasta(path)
    regions = load_regions(args.regions, families)
    elm_resources, elm_focus, elm_summary = None, None, None
    if args.elm_dir:
        import elm_annotation
        elm_resources = elm_annotation.load_resources(args.elm_dir, read_fasta)
        elm_focus = elm_annotation.load_focus(args.elm_focus, elm_resources['classes'], families)
        # Validate reference coordinates before creating any output.
        elm_annotation.reference_audit(elm_resources)
    if args.known_meme:
        check_meme(args.known_meme)
        args.known_meme = args.known_meme.resolve()
    args.out = args.out.resolve()
    args.out.mkdir(parents=True)
    out = args.out
    (out/'motif_library.snapshot.json').write_text(json.dumps(library, indent=2), encoding='utf-8')
    qc, summary = [], []
    membership = defaultdict(set)
    for f, records in families.items():
        for r in records:
            membership[r['sequence']].add(f)
    for f, records in families.items():
        first = {}
        for r in records:
            seq = r['sequence']
            qc.append(dict(family=f, sequence_id=r['id'], length=len(seq),
                           duplicate_of=first.get(seq, ''), ambiguous=sum(x not in AA for x in seq),
                           terminal_stop_removed=r['terminal_stop'], shared_families=','.join(sorted(membership[seq]-{f})),
                           sha256=hashlib.sha256(seq.encode()).hexdigest()))
            first.setdefault(seq, r['id'])
        summary.append(dict(family=f, sequences=len(records), unique_sequences=len(first),
                            min_length=min(len(r['sequence']) for r in records), max_length=max(len(r['sequence']) for r in records),
                            configured_motifs=len(library['families'][f]['motifs']), note=library['families'][f]['note']))
    hits, states, prevalence = inspect(families, library, regions)
    words, sites = discover(families, args.min_k, args.max_k, args.min_fraction, args.top)
    domain_rows = inspect_domains(families, library, regions)
    for row in words:
        annotations = [name for name, spec in library['motifs'].items() if list(occurrences(spec['pattern'], row['motif']))]
        row['overlapping_known_patterns'] = ','.join(annotations)
        row['annotation'] = 'sequence_pattern_overlap_only' if annotations else 'unannotated'
    if elm_resources:
        print('Scanning supplied ELM classes against distinct protein sequences...', flush=True)
        elm_summary = elm_annotation.run(elm_resources, families, regions, elm_focus, sites, out, write_tsv)
    stages, tools, versions = external(families, out, args)
    tables = [
        ('sequence_qc', qc, 'family sequence_id length duplicate_of ambiguous terminal_stop_removed shared_families sha256'),
        ('family_summary', summary, 'family sequences unique_sequences min_length max_length configured_motifs note'),
        ('known_hits', hits, 'family sequence_id motif start end peptide context_state relative_start flank interpretation source'),
        ('inspection_states', states, 'family sequence_id motif state hits'),
        ('known_prevalence', prevalence, 'family motif target_positive target_total other_positive other_total target_fraction other_fraction comparison'),
        ('discovered_words', words, 'family motif target_positive target_total other_positive other_total target_fraction other_fraction prevalence_difference status overlapping_known_patterns annotation'),
        ('domain_inspection', domain_rows, 'family sequence_id state observed missing_alternatives interpretation'),
        ('discovered_sites', sites, 'family sequence_id motif start end'),
        ('external_stages', stages, 'family stage status reason')]
    for name, rows, fields in tables:
        write_tsv(out/(name+'.tsv'), rows, fields.split())
    manifest = dict(created_utc=datetime.now(timezone.utc).isoformat(), python=platform.python_version(),
                    parameters={k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()},
                    inputs={str(p.resolve()): digest(p) for p in files+[args.library, Path(__file__)] +
                            ([args.regions] if args.regions else []) + ([args.known_meme] if args.known_meme else []) +
                            (elm_resources['paths']+[args.elm_focus, HERE/'elm_annotation.py'] if elm_resources else [])},
                    elm=elm_summary,
                    tool_paths=tools, tool_versions=versions,
                    inspection_completed=True,
                    external_execution_complete=all(s['status']=='completed' for s in stages),
                    completed=not any(s['status']=='failed' for s in stages))
    (out/'manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    report = ['# Family motif inspection', '',
              f'{len(families)} families; {sum(len(rs) for rs in families.values())} input records; {len(hits)} known-pattern occurrences.', '',
              'All coordinates are 1-based inclusive. Family labels are supplied candidate assignments, not validated identities.',
              'Short patterns and exact-word discovery do not prove function, binding, phosphorylation, or family membership.',
              'Context-unverified hits require domain/topology evidence. Tail windows are not predicted cytoplasmic tails.',
              'Discovery compares exact-deduplicated families against other supplied candidate families, excluding shared sequences.',
              'These controls are not a matched proteome background. Length, composition, isoforms and shared ancestry can bias results.',
              'Exact words are exploratory, ranked on training data, with no p-values or claims of independent validation.',
              'FIMO scans of discovered motifs on their training sequences are localization only. External q-values are per tool run, not globally corrected.', '',
              '| Family | Records | Distinct sequences | Length range | Configured patterns |',
              '|---|---:|---:|---:|---:|']
    report += [f"| {r['family']} | {r['sequences']} | {r['unique_sequences']} | {r['min_length']}–{r['max_length']} | {r['configured_motifs']} |" for r in summary]
    report += ['', '## External stages', '', ', '.join(f'{k}: {v}' for k,v in Counter(s['status'] for s in stages).items()), '',
               'See external_stages.tsv for every completed, planned, failed or skipped step and commands.json for exact argument arrays.',
               'Protein PWM annotation requires --known-meme. No PWM database or domain coordinates are invented.', '']
    if elm_summary:
        report += ['', '## ELM annotation', '',
                   f"Scanned {elm_summary['classes']} ELM classes on {elm_summary['unique_sequences']} distinct sequences: "
                   f"{elm_summary['regex_hit_rows']} candidate occurrences across input records; "
                   f"{elm_summary['family_focus_hit_rows']} occurrences in the family review panel.", '',
                   f"Reference check: {elm_summary['exact_reference_spans_recovered']}/{elm_summary['reference_instances']} supplied instance spans recovered exactly.", '',
                   'Start with elm_family_focus_hits.tsv, elm_inspection_states.tsv and elm_family_prevalence.tsv. '
                   'elm_hits.tsv retains all regex candidates; family priority is a review aid, not confidence or classification.',
                   'ELM probabilities are class metadata, not hit p-values. Disorder, PTMs, cellular localization, taxonomic applicability and binding are not inferred.',
                   'Domain overlaps flag context for review; absence of a supplied domain does not imply disorder. Structured SLiMs are not automatically rejected.',
                   'Interaction domains identify potential binders, not domains required on motif-bearing proteins. No candidate interaction edges are inferred.',
                   'KEGG/Reactome links remain associations of ELM reference proteins; they are not pathway assignments to the candidates.',
                   'elm_discovery_annotation.tsv links discovered sites to ELM matches in the original full sequence, preserving terminal anchors and flanks.',
                   'See elm_resource_audit.tsv for optional missing files, download metadata and unresolved historical class references.', '']
    (out/'REPORT.md').write_text('\n'.join(report), encoding='utf-8')
    print(f'Wrote {out / "REPORT.md"}')
    return 1 if any(s['status']=='failed' for s in stages) else 0


if __name__ == '__main__':
    try:
        sys.exit(main())
    except (ValueError, OSError, KeyError) as exc:
        print(f'Error: {exc}', file=sys.stderr)
        sys.exit(2)
