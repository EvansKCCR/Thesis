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
VERSION = '1.0'


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
    def ordered_groups_present(intervals, groups):
        last = 0
        for choices in groups:
            matches = [(a, b) for name, a, b in intervals if name in choices and a > last]
            if not matches: return False
            last = min(matches)[0]
        return True
    rows = []
    for family, records in families.items():
        spec = library['families'][family]
        alternatives = spec.get('domain_alternatives', [])
        for r in records:
            # Unique intervals prevent duplicate imported rows inflating repeat counts.
            intervals = sorted((name,a,b) for name,a,b in set(regions[family, r['id']]) if name.startswith('PF'))
            counts = Counter(name for name, a, b in intervals)
            missing = [{name: n-counts[name] for name, n in alt.items() if counts[name] < n} for alt in alternatives]
            count_supported = any(not m for m in missing)
            order_supported = ordered_groups_present(intervals, spec['ordered_domain_groups']) if spec.get('support_requires_order') else None
            rows.append(dict(family=family, sequence_id=r['id'],
                             state=('not_configured_screening_hypothesis' if spec.get('classification_level')=='screening_hypothesis' else
                                    'not_configured_component_hypothesis') if not alternatives else 'not_tested_missing_domains' if not counts else
                             'domain_count_order_rule_supported' if count_supported and order_supported else
                             'domain_count_order_conflict' if count_supported and order_supported is False else
                             'domain_count_rule_supported' if count_supported else 'domain_count_rule_not_met',
                             order_status='supported' if order_supported else 'not_met' if order_supported is False else 'not_required',
                             observed=json.dumps(dict(counts), sort_keys=True),
                             missing_alternatives=json.dumps(missing, sort_keys=True),
                             interpretation='Pfam count/order screen only; homology and family identity remain unverified.'))
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
    import curated_annotation as curated
    import master_library
    import family_identity_review
    import import_sequence_context
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--version',action='version',version='adhesome-motif '+VERSION)
    parser.add_argument('--input-dir', type=Path, default=HERE/'family_specific_candidate_library')
    parser.add_argument('--library', type=Path, help='Current JSON/XLSX motif snapshot; defaults to adhesion_core.json')
    parser.add_argument('--motif-scope',choices=['adhesion','master'],default='adhesion',help='53 adhesion rules by default, or the complete 347-rule snapshot')
    parser.add_argument('--rules', type=Path, default=HERE/'Motif_library'/'context_rules.json')
    parser.add_argument('--out', type=Path, required=True, help='New output directory')
    parser.add_argument('--regions', type=Path, help='1-based inclusive domain/context TSV')
    parser.add_argument('--topology-dir', type=Path, help='Directory containing predicted_topologies.3line files recursively')
    parser.add_argument('--topology-file', type=Path, help='Single DeepTMHMM predicted_topologies.3line file')
    parser.add_argument('--topology-zip', type=Path, nargs='+', default=[], help='Additional DeepTMHMM ZIP archives, read without extraction')
    parser.add_argument('--context-dir', type=Path, help='Directory with DeepLoc, SignalP and NCBI CDD files; defaults to current adhesome_candidate data for the bundled FASTAs')
    parser.add_argument('--evidence', type=Path, help='Sequence-hash-verified, site-specific evidence TSV')
    parser.add_argument('--known-meme', type=Path, help='Optional actual protein PWM database')
    parser.add_argument('--external', choices=['skip','plan','run'], default='skip')
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
    args=parser.parse_args(argv)
    if not (2 <= args.min_k <= args.max_k <= 15 and 3 <= args.min_width <= args.max_width <= 30 and
            0 < args.min_fraction <= 1 and 0 < args.qvalue <= 1 and args.top > 0 and args.min_streme >= 2 and args.timeout > 0):
        parser.error('Invalid discovery parameters')
    if args.out.exists(): parser.error('Output directory already exists; choose a new --out')
    rules=json.loads(args.rules.read_text(encoding='utf-8'))
    domain_source=(args.rules.resolve().parent/rules['domain_rule_source']).resolve()
    current_rules=json.loads(domain_source.read_text(encoding='utf-8'))
    rules['families'].update(current_rules['families'])
    for alias,canonical in current_rules.get('family_aliases',{}).items():
        if alias not in rules['families'] and canonical in rules['families']:
            rules['families'][alias]=rules['families'][canonical]
    for alias,canonical in rules.get('family_aliases',{}).items():
        rules['families'][alias]=rules['families'][canonical]
    for spec in rules['families'].values():
        if spec.get('ordered_domains') and not spec.get('ordered_domain_groups'):
            spec['ordered_domain_groups']=[[p] for p in spec['ordered_domains']]
    args.library=args.library or HERE/'Motif_library'/('adhesion_core.json' if args.motif_scope=='adhesion' else 'motif_master.json')
    snapshot=master_library.load(args.library,args.motif_scope)
    files=sorted(p for p in args.input_dir.iterdir() if p.suffix.lower() in ('.fa','.faa','.fasta'))
    if not files: parser.error('No FASTA files')
    families={}
    for path in files:
        if path.stem not in rules['families'] or path.stem in families: parser.error('Unknown or duplicate family: '+path.stem)
        families[path.stem]=read_fasta(path)
    regions=load_regions(args.regions,families)
    if args.context_dir is None and args.input_dir.resolve()==(HERE/'family_specific_candidate_library').resolve():
        args.context_dir=HERE.parent/'02_topology_subcellular localization'/'adhesome_candidate'
    topology_paths=sorted(args.topology_dir.rglob('predicted_topologies.3line')) if args.topology_dir else []
    if args.topology_dir and not topology_paths: parser.error('No predicted_topologies.3line files found')
    if args.topology_file: topology_paths.append(args.topology_file)
    topology_paths += args.topology_zip
    if args.context_dir and not topology_paths:
        topology_paths.append(args.context_dir/'deeptmhmm_predicted_topologies.3line')
    topology_audit=curated.import_topology(topology_paths,families,regions)
    sequence_context,cdd_hits=import_sequence_context.load(args.context_dir,families,topology_audit,regions,include_all=True) if args.context_dir else ({},[])
    evidence=master_library.load_evidence(args.evidence,families,snapshot)
    if args.known_meme:
        check_meme(args.known_meme)
        args.known_meme=args.known_meme.resolve()
    args.out=args.out.resolve(); out=args.out
    out.mkdir(parents=True)
    write_tsv(out/'sequence_context.tsv',sequence_context.values(),
              'family sequence_id verification deeploc_localizations deeploc_cytoplasm deeploc_extracellular signalp_prediction signalp_cleavage_after cdd_specific_hits'.split())
    cdd_fields='family sequence_id accession short_name start end evalue bitscore hit_type model_identifier incomplete superfamily model_coverage specific_bitscore_cutoff specific_evalue_cutoff category_source interpretation'.split()
    write_tsv(out/'CDD_hits.tsv',cdd_hits,cdd_fields)
    write_tsv(out/'CDD_specific_hits.tsv',[h for h in cdd_hits if h['hit_type']=='specific' and float(h['evalue'])<=1e-3],cdd_fields)
    identity_rows=family_identity_review.review(families,regions)
    write_tsv(out/"PINCH_Parvin_identity_review.tsv",identity_rows,"family sequence_id length sequence_sha256 LIM_domains CH_domains actin_domains spectrin_domains hypothesis review_priority identity_status missing_evidence interpretation".split())
    domain_rows=inspect_domains(families,rules,regions)
    print('Scanning current evidence-based motif definitions and generating exploratory candidates...',flush=True)
    discovery_families={name:records for name,records in families.items()
                        if rules['families'][name].get('classification_level') not in ('component_hypothesis','screening_hypothesis')}
    words,sites=discover(discovery_families,args.min_k,args.max_k,args.min_fraction,args.top)
    master_result=master_library.scan(snapshot,families,regions,evidence,out,write_tsv,sequence_context)
    result=dict(annotation_rows=master_result['sequence_annotation_rows'],context_filters_implemented=False)
    stages,tool_paths,versions=external(discovery_families,out,args)
    stages.extend(dict(family=name,stage='family_discovery',status='skipped',
                       reason='screening_bucket_has_no_assigned_family' if rules['families'][name].get('classification_level')=='screening_hypothesis'
                       else 'component_bucket_has_no_assigned_family')
                  for name in families if name not in discovery_families)
    membership=defaultdict(set)
    for family,records in families.items():
        for r in records: membership[r['sequence']].add(family)
    qc=[]; summary=[]
    for family,records in families.items():
        first={}
        for r in records:
            seq=r['sequence']
            qc.append(dict(family=family,sequence_id=r['id'],length=len(seq),duplicate_of=first.get(seq,''),
                ambiguous=sum(x not in AA for x in seq),terminal_stop_removed=r['terminal_stop'],
                shared_families=';'.join(sorted(membership[seq]-{family})),sha256=hashlib.sha256(seq.encode()).hexdigest()))
            first.setdefault(seq,r['id'])
        summary.append(dict(family=family,sequences=len(records),unique_sequences=len(first),
                            min_length=min(len(r['sequence']) for r in records),max_length=max(len(r['sequence']) for r in records)))
    tables=[('sequence_qc',qc,'family sequence_id length duplicate_of ambiguous terminal_stop_removed shared_families sha256'),
            ('family_summary',summary,'family sequences unique_sequences min_length max_length'),
            ('domain_inspection',domain_rows,'family sequence_id state order_status observed missing_alternatives interpretation'),
            ('topology_audit',topology_audit,'family sequence_id status source'),
            ('discovered_words',words,'family motif target_positive target_total other_positive other_total target_fraction other_fraction prevalence_difference status'),
            ('discovered_sites',sites,'family sequence_id motif start end'),
            ('external_stages',stages,'family stage status reason')]
    for name,rows,fields in tables: write_tsv(out/(name+'.tsv'),rows,fields.split())
    write_tsv(out/'used_regions.tsv',[dict(family=f,sequence_id=i,region=n,start=a,end=b) for (f,i),rs in regions.items() for n,a,b in sorted(set(rs))],
              ['family','sequence_id','region','start','end'])
    inputs=files+[args.library,args.rules,domain_source,Path(__file__),HERE/'curated_annotation.py',HERE/'master_library.py',HERE/'family_identity_review.py',HERE/'import_sequence_context.py']+topology_paths
    if args.context_dir:
        inputs += [args.context_dir/name for name in ('DeepLoc_adhesome_candidate_list.csv','SignalP_adhesome_candidate_list.txt','NCBI_CDD_adhesome_hitdata.txt')]
    inputs += [x for x in (args.regions,args.evidence,args.known_meme) if x]
    failed=any(s['status']=='failed' for s in stages)
    manifest=dict(pipeline_version=VERSION,created_utc=datetime.now(timezone.utc).isoformat(),python=platform.python_version(),
        parameters={k:str(v) if isinstance(v,Path) else [str(x) for x in v] if isinstance(v,list) else v for k,v in vars(args).items()},
        inputs={str(p.resolve()):digest(p) for p in inputs},tool_paths=tool_paths,tool_versions=versions,
        inspection_completed=True,external_execution_complete=all(s['status']=='completed' for s in stages),completed=not failed,curated=result,master=master_result)
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    report=['# Family motif inspection', '',
        f"{len(families)} families; {sum(len(r) for r in families.values())} records; {master_result['master_records']} master library entries.",
        f"{master_result['sequence_annotation_rows']} sequence annotation rows, including {master_result['generic_hit_rows']} generic catalogue matches.",
        master_result['warning'],
        *master_result['summary_discrepancies'],
        'PINCH_Parvin_identity_review.tsv retains original labels and lists architecture triage and missing identity evidence.',
        f'Sequence context: {sum(r["verification"]=="exact_topology_sequence" for r in sequence_context.values())} family records have exact topology sequences; {len(cdd_hits)} CDD domain predictions are retained; the scoring export preserves the previous specific-hit selection.',
        'DeepLoc, SignalP and CDD are predictions. Their joins and coordinates are in sequence_context.tsv and CDD_specific_hits.tsv; they do not establish orthology or motif function.',
        'motif_annotations.tsv contains current definition matches; adhesion_motif_annotations.tsv is the adhesion subset. All hits remain sequence candidates; context requirements and external evidence are recorded without assigning functional confidence.',
        'Pathway links are class-reference associations, not candidate pathway assignments. Disabled legacy records are retained in legacy_review_catalogue.tsv and never scanned.',
        '', 'See [curated findings](CURATED_FINDINGS.md) for confidence tiers, evidence gaps and output interpretation.',
        '', 'The current XLSX is the sole motif source. No deleted ELM resource or previous motif panel is needed.',
        'Library rule confidence is not candidate confidence. Functional hypotheses are conditional, not validated functions.',
        'Domain counts are supporting screens, not full family classification. Exact-sequence topology imports are predictions.',
        'Exploratory words and cross-family prevalence are descriptive; they do not prove statistical enrichment or orthologous conservation.',
        '', '## Execution', '',
        'Topology: '+str(dict(Counter(r['status'] for r in topology_audit))),
        'External stages: '+str(dict(Counter(r['status'] for r in stages))),
        'STREME/FIMO/TOMTOM require installed tools and an actual protein PWM database where applicable. The XLSX is not a PWM database.', '']
    (out/'REPORT.md').write_text('\n'.join(report),encoding='utf-8')
    print(f'Wrote {out / "REPORT.md"}')
    return 1 if failed else 0


if __name__ == '__main__':
    try:
        sys.exit(main())
    except (ValueError, OSError, KeyError) as exc:
        print(f'Error: {exc}', file=sys.stderr)
        sys.exit(2)
