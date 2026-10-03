#!/usr/bin/env python3
"""Auditable preparation, domain scanning, homology recovery and evidence integration."""
import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
from adhesome_core import fasta, write_fasta, table, write_table, digest, save_json, pf, sequence_dict
from parse_adhesome_domtblout import parse, FIELDS

PROJECT_ROOT = Path(__file__).resolve().parent
CONFIG = PROJECT_ROOT / 'config'
REFERENCE_SET = PROJECT_ROOT / 'references' / 'schistosome_v5'

MODULES = {'extracellular_matrix', 'integrin_receptors', 'membrane_proximal_adaptors',
           'actin_linking', 'focal_adhesion_scaffolds', 'downstream_signalling'}

def prepare(config, output):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    settings = json.loads(Path(config).read_text(encoding='utf-8'))
    qc = settings['qc']
    if not (qc['min_length'] > 0 and 0 <= qc['max_ambiguous_fraction'] <= 1 and qc['max_ambiguous_run'] >= 0):
        raise ValueError('Invalid sequence QC thresholds')
    codes = [s['code'] for s in settings['proteomes']]
    if len(codes) != len(set(codes)) or any(not re.fullmatch(r'[A-Za-z0-9_]+', c) for c in codes):
        raise ValueError('Species codes must be unique safe identifiers')
    summaries = {}
    for species in settings['proteomes']:
        source = Path(config).resolve().parent / species['path']
        records, identifiers, groups = [], set(), defaultdict(list)
        for header, sequence in fasta(source):
            identifier = header.split()[0]
            if identifier in identifiers:
                raise ValueError(f'Duplicate protein ID in {source}: {identifier}')
            identifiers.add(identifier)
            match = re.search(species['gene_regex'], header)
            if not match:
                raise ValueError(f'No gene mapping for {identifier}; supply a validated gene_regex')
            gene = match.group(1)
            if not re.fullmatch(r'[A-Za-z0-9_.-]+', gene):
                raise ValueError(f'Unsafe gene ID: {gene}')
            sequence = sequence[:-1] if sequence.endswith('*') else sequence
            row = dict(original_id=identifier, gene_id=gene, original_header=header,
                       protein_id=species['code'] + '__' + gene, length=len(sequence), sequence=sequence)
            groups[gene].append(row)
            records.append(row)
        if not records:
            raise ValueError(f'Empty FASTA: {source}')
        for gene, isoforms in groups.items():
            # Collapse first, then QC. Do not silently replace a bad longest isoform.
            winner = sorted(isoforms, key=lambda r: (-r['length'], r['original_id']))[0]
            for row in isoforms:
                sequence = row['sequence']
                ambiguous = sum(c in 'XBZJ' for c in sequence)
                row['ambiguous_fraction'] = ambiguous / len(sequence) if sequence else 1.0
                runs = [len(m.group()) for m in re.finditer('[XBZJ]+', sequence)]
                row['longest_ambiguous_run'] = max(runs, default=0)
                reasons = []
                if row is not winner:
                    reasons.append('shorter_isoform_or_tie')
                if set(sequence) - set('ACDEFGHIKLMNPQRSTVWYXBZJUO'):
                    reasons.append('invalid_residue_or_internal_stop')
                if row['length'] < settings['qc']['min_length']:
                    reasons.append('below_global_min_length')
                if row['ambiguous_fraction'] > settings['qc']['max_ambiguous_fraction']:
                    reasons.append('ambiguous_fraction')
                if row['longest_ambiguous_run'] > settings['qc']['max_ambiguous_run']:
                    reasons.append('ambiguous_run')
                row['decision'] = ';'.join(reasons) if reasons else 'retained'
        retained = sorted((r for r in records if r['decision'] == 'retained'), key=lambda r: r['protein_id'])
        prefix = output / species['code']
        write_fasta(str(prefix) + '.faa', [(r['protein_id'], r['sequence']) for r in retained])
        write_table(str(prefix) + '.id_map.csv', records,
                    ['protein_id','original_id','gene_id','length','ambiguous_fraction','longest_ambiguous_run','decision','original_header'])
        summaries[species['code']] = dict(species=species['species'], assembly=species['assembly'],
            assembly_status='user_supplied_not_verified_from_FASTA', source_url=species.get('source_url'),
            input=str(source), input_sha256=digest(source), input_proteins=len(records), genes=len(groups),
            retained=len(retained), excluded=len(records)-len(retained),
            decisions=dict(Counter(r['decision'] for r in records)), output_sha256=digest(str(prefix) + '.faa'))
    save_json(output / 'preparation.json', dict(config=settings, config_sha256=digest(config), species=summaries))
    print(json.dumps(summaries, indent=2))

def run(command, log):
    if not shutil.which(command[0]):
        raise RuntimeError(f'{command[0]} is not installed/on PATH. See README_adhesome_pfam.md.')
    executable = shutil.which(command[0])
    flag = '-version' if command[0] == 'blastp' else '-h'
    version = subprocess.run([executable, flag], capture_output=True, text=True, check=False)
    save_json(str(log) + '.command.json', dict(command=[str(x) for x in command],
        executable=executable, executable_sha256=digest(executable),
        version_output=(version.stdout + version.stderr)[:8000]))
    with open(log, 'w', encoding='utf-8') as handle:
        subprocess.run([str(x) for x in command], stdout=handle, stderr=subprocess.STDOUT, check=True)

def scan(args):
    output = Path(args.out)
    output.mkdir(parents=True, exist_ok=True)
    metadata = json.loads(Path(args.library + '.metadata.json').read_text(encoding='utf-8'))
    if digest(args.library) != metadata['library_sha256']:
        raise ValueError('Library checksum mismatch')
    if digest(args.manifest) != metadata['manifest_sha256']:
        raise ValueError('Manifest changed since library extraction; rebuild library')
    mapping = {pf(r['Pfam accession']): r['Domain or family'] for r in table(args.manifest)}
    commands, parts = [], []
    models = metadata['models']
    if metadata.get('full_pfam'):
        mapping.update({key:value['name'] for key,value in models.items() if key not in mapping})
    for kind in ('ga', 'fallback'):
        library = args.library + '.' + kind + '.hmm'
        if digest(library) != metadata['split_sha256'][kind]:
            raise ValueError(f'{kind} library checksum mismatch')
        if not Path(library).stat().st_size:
            continue
        command = ['hmmpress', '-f', library]
        commands.append(command)
        run(command, output / (kind + '.hmmpress.log'))
        part = output / (kind + '.domtblout')
        command = ['hmmscan', '--cpu', str(args.cpu), '--noali', '-Z', str(len(models)),
                   '--domZ', str(len(models)), '--domtblout', str(part)]
        # Shared search space prevents split-database E-value inflation/deflation.
        command += ['--cut_ga'] if kind == 'ga' else ['-E', '100', '--domE', '100']
        command += [library, args.proteome]
        commands.append(command)
        run(command, output / (kind + '.hmmscan.log'))
        parts.append(part)
    combined = output / 'combined.domtblout'
    combined.write_text(''.join(p.read_text(encoding='utf-8') for p in parts), encoding='utf-8')
    rows = parse(combined, mapping, models, args.max_i, args.max_c, args.min_coverage)
    write_table(output / 'domains.csv', [r for r in rows if r['decision'] == 'retained'], FIELDS)
    write_table(output / 'domains.audit.csv', rows, FIELDS)
    save_json(output / 'scan.json', dict(commands=commands, proteome_sha256=digest(args.proteome),
        metadata_sha256=digest(args.library + '.metadata.json'), thresholds=dict(i=args.max_i,c=args.max_c,
        profile_coverage=args.min_coverage), search_space=len(models), python=sys.version))

def catalogue(path, sequence_path, rules=None):
    rows = table(path)
    required = ['reference_id','family','module','species','accession','evidence_source','experimental_evidence']
    seqs = sequence_dict(sequence_path)
    if not rows:
        raise ValueError('Reference catalogue is empty; populate experimentally supported references first')
    found = {}
    for row in rows:
        if any(not row.get(k, '').strip() for k in required):
            raise ValueError('Reference catalogue has missing required metadata')
        key = row['reference_id']
        if key not in seqs or key in found or row['module'] not in MODULES:
            raise ValueError(f'Invalid/duplicate reference or module: {key}')
        if rules is not None and row['family'] not in rules:
            raise ValueError(f'Family has no classification rule: {row["family"]}')
        if row.get('sequence_sha256'):
            import hashlib
            if hashlib.sha256(seqs[key].encode()).hexdigest() != row['sequence_sha256']:
                raise ValueError(f'Catalogue/sequence checksum mismatch: {key}')
        found[key] = row
    if set(seqs) != set(found):
        raise ValueError('Reference FASTA/catalogue IDs do not match exactly')
    return found

BLAST_FIELDS = 'qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen'

def blast_hits(path, evalue, coverage):
    hits = []
    with open(path, encoding='utf-8') as handle:
        for line in handle:
            x = line.split()
            if len(x) != 12:
                raise ValueError(f'Malformed BLAST row in {path}')
            qcov = (abs(int(x[5])-int(x[4]))+1)/int(x[10])
            scov = (abs(int(x[7])-int(x[6]))+1)/int(x[11])
            if float(x[8]) <= evalue and min(qcov, scov) >= coverage:
                hits.append(dict(query=x[0], subject=x[1], score=float(x[9]), evalue=float(x[8]), qcov=qcov, scov=scov,
                                 subject_start=int(x[6]),subject_end=int(x[7])))
    return hits

def homology(args):
    refs = catalogue(args.catalogue, args.references)
    reference_sequences = sequence_dict(args.references)
    broad_rows = table(args.metazoan_catalogue)
    broad = {r['reference_id']: r['family'] for r in broad_rows}
    broad_sequences = sequence_dict(args.metazoan)
    broad_ids = set(broad_sequences)
    if len(broad) != len(broad_rows) or set(broad) != broad_ids:
        raise ValueError('Broad metazoan catalogue/FASTA mismatch or duplicate IDs')
    if not any(not family or family == 'non_adhesome' for family in broad.values()):
        raise ValueError('Reciprocal database must include annotated non-adhesome competitors')
    output = Path(args.out)
    output.mkdir(parents=True, exist_ok=True)
    import hashlib
    for row in broad_rows:
        if row.get('sequence_sha256') and hashlib.sha256(broad_sequences[row['reference_id']].encode()).hexdigest() != row['sequence_sha256']:
            raise ValueError(f'Reciprocal catalogue/sequence checksum mismatch: {row["reference_id"]}')
    snapshot = output / 'reference_snapshot'
    snapshot.mkdir(exist_ok=True)
    paths = {'search_references.faa':args.references,'search_catalogue.csv':args.catalogue,
             'metazoan_references.faa':args.metazoan,'metazoan_catalogue.csv':args.metazoan_catalogue}
    lock_path = Path(args.references).resolve().parent / 'reference_lock.json'
    if lock_path.exists():
        lock = json.loads(lock_path.read_text(encoding='utf-8'))
        for name, source in paths.items():
            if lock['files'].get(Path(source).name) != digest(source):
                raise ValueError(f'Reference lock mismatch: {name}; rebuild or select a consistent reference set')
    # Every run keeps the exact search sequences, not just a pointer to mutable databases.
    for name, source in paths.items():
        destination = snapshot / name
        if destination.exists() and digest(destination) != digest(source):
            raise ValueError('Reference snapshot differs; choose a new output directory for this reference release')
        if Path(source).resolve() != destination.resolve():
            shutil.copy2(source,destination)
    for name in ('reference_lock.json','experimental_annotations.csv','resolution_audit.json'):
        source = Path(args.references).resolve().parent / name
        if source.exists():
            shutil.copy2(source,snapshot/name)
    commands = []
    for label, query, subject in [('forward', str(snapshot/'search_references.faa'), args.proteome),
                                  ('reverse', args.proteome, str(snapshot/'metazoan_references.faa'))]:
        command = ['blastp', '-query', query, '-subject', subject, '-evalue', str(args.evalue),
                   '-num_threads', str(args.cpu), '-max_target_seqs', '1000000', '-max_hsps', '1',
                   '-outfmt', '6 ' + BLAST_FIELDS, '-out', str(output / (label + '.tsv'))]
        commands.append(command)
        run(command, output / (label + '.log'))
    forward = blast_hits(output / 'forward.tsv', args.evalue, args.coverage)
    reverse = blast_hits(output / 'reverse.tsv', args.evalue, args.coverage)
    best = defaultdict(list)
    for hit in reverse:
        best[hit['query']].append(hit)
    rows = []
    for hit in forward:
        family = refs[hit['query']]['family']
        reverse_hits = best.get(hit['subject'], [])
        top_score = max((h['score'] for h in reverse_hits), default=0)
        top = [h for h in reverse_hits if h['score'] >= top_score * (1-args.tie_fraction)]
        families = {broad[h['subject']] for h in top}
        supported = families == {family}
        rows.append(dict(protein_id=hit['subject'], family=family, reference_id=hit['query'],
                         similarity='true', reciprocal='true' if supported else 'false',
                         evalue=hit['evalue'], bitscore=hit['score'], query_coverage=hit['qcov'],
                         subject_coverage=hit['scov'], reverse_families=';'.join(sorted(families)),
                         subject_start=hit['subject_start'],subject_end=hit['subject_end'],
                         reverse_reference_ids=';'.join(sorted({h['subject'] for h in top})),
                         reference_accession=refs[hit['query']]['accession'],reference_species=refs[hit['query']]['species'],
                         reference_evidence_source=refs[hit['query']]['evidence_source'],
                         reference_evidence_status=refs[hit['query']].get('evidence_status','legacy_unverified'),
                         reference_tier=refs[hit['query']].get('reference_tier','legacy_unverified'),
                         reference_phylogenetic_group=refs[hit['query']].get('phylogenetic_group',''),
                         reverse_reference_tiers=';'.join(sorted({r.get('reference_tier','legacy_unverified') for r in broad_rows if r['reference_id'] in {h['subject'] for h in top}})),
                         reference_sequence_sha256=hashlib.sha256(reference_sequences[hit['query']].encode()).hexdigest(),
                         source=str(snapshot/'search_catalogue.csv')))
    write_table(output / 'homology.csv', rows, ['protein_id','family','reference_id','similarity','reciprocal',
                'evalue','bitscore','query_coverage','subject_coverage','subject_start','subject_end','reverse_families','reverse_reference_ids',
                'reference_accession','reference_species','reference_evidence_source','reference_evidence_status','reference_tier',
                'reference_phylogenetic_group','reverse_reference_tiers','reference_sequence_sha256','source'])
    matched_ids = sorted({r['reference_id'] for r in rows})
    write_fasta(output/'matched_reference_sequences.faa',[(rid,reference_sequences[rid]) for rid in matched_ids])
    write_table(output/'matched_reference_catalogue.csv',[refs[rid] for rid in matched_ids],list(next(iter(refs.values()))))
    save_json(output / 'homology.json', dict(commands=commands, evalue=args.evalue, coverage=args.coverage,
        tie_fraction=args.tie_fraction, reference_snapshot=str(snapshot),snapshot_source_mapping={name:str(p) for name,p in paths.items()},
        snapshot_checksums={p.name:digest(p) for p in snapshot.iterdir() if p.is_file()},
        checksums={str(p): digest(p) for p in [args.proteome,args.references,
        args.catalogue,args.metazoan,args.metazoan_catalogue]}))

def architecture_matches(domains, rule):
    if rule.get('ordered_domain_groups'):
        from adhesome_classification import architecture_evidence
        evidence=architecture_evidence(domains,rule)
        return evidence['complete'] and evidence['order']=='consistent'
    counts = Counter(d['pfam_accession'] for d in domains)
    if not any(all(counts[p] >= n for p, n in group.items()) for group in rule['domain_alternatives']):
        return False
    order = rule.get('ordered_domains', [])
    position = 0
    for domain in domains:
        if position < len(order) and domain['pfam_accession'] == order[position]:
            position += 1
    return position == len(order)

def profiles(args):
    """Search externally curated, alignment-derived family HMMs (hmmsearch orientation)."""
    mapping = {r['model_name']: r['family'] for r in table(args.profile_map)}
    output = Path(args.out)
    output.mkdir(parents=True, exist_ok=True)
    command = ['hmmsearch', '--cpu', str(args.cpu), '--noali', '-E', '100', '--domE', '100',
               '--domtblout', str(output / 'profiles.domtblout'), args.profiles, args.proteome]
    run(command, output / 'hmmsearch.log')
    rows = []
    with open(output / 'profiles.domtblout', encoding='utf-8') as handle:
        for line in handle:
            if not line.strip() or line.startswith('#'):
                continue
            x = line.split(maxsplit=22)
            if len(x) < 22 or x[3] not in mapping:
                raise ValueError('Malformed profile result or unmapped model')
            coverage = (int(x[16])-int(x[15])+1)/int(x[5])
            if float(x[12]) <= args.max_i and float(x[11]) <= args.max_c and coverage >= args.min_coverage:
                rows.append(dict(protein_id=x[0], family=mapping[x[3]], profile='true',
                            source=f'hmmsearch:{x[3]}:{digest(args.profiles)}'))
    write_table(output / 'profile_evidence.csv', rows, ['protein_id','family','profile','source'])
    save_json(output / 'profiles.json', dict(command=command, i=args.max_i,c=args.max_c,
        coverage=args.min_coverage, checksums={str(p):digest(p) for p in [args.proteome,args.profiles,args.profile_map]}))

def classify_strict(args):
    config = json.loads(Path(args.rules).read_text(encoding='utf-8'))
    rules = config['families']
    from adhesome_motifs import load_library,scan_motifs,noncatalytic_alignment
    motif_path=(Path(args.rules).resolve().parent/config['motif_library']).resolve() if config.get('motif_library') else None
    motif_library=load_library(motif_path) if motif_path else {}
    sequences = sequence_dict(args.proteome)
    for family, rule in rules.items():
        if rule['module'] not in MODULES or rule['min_length'] < 1 or not rule['domain_alternatives']:
            raise ValueError(f'Invalid family rule: {family}')
    domains, evidence = defaultdict(list), defaultdict(list)
    for row in table(args.domains):
        if row['query_name'] not in sequences:
            raise ValueError(f'Unknown domain protein: {row["query_name"]}')
        domains[row['query_name']].append(row)
    for path in (args.homology, *(args.evidence or [])):
        if path:
            for row in table(path):
                if row['protein_id'] not in sequences or row['family'] not in rules:
                    raise ValueError('Unknown evidence protein/family')
                if not row.get('source'):
                    raise ValueError('Evidence must include a source')
                evidence[(row['protein_id'], row['family'])].append(row)
    results, architecture_rows = [], []
    for identifier, sequence in sequences.items():
        hits = sorted(domains[identifier], key=lambda r: (int(r['ali_from']), int(r['ali_to'])))
        architecture = ';'.join(f'{r["pfam_accession"]}:{r["ali_from"]}-{r["ali_to"]}' for r in hits)
        if hits:
            architecture_rows.append(dict(protein_id=identifier, length=len(sequence), architecture=architecture))
        for family, rule in rules.items():
            support = evidence[(identifier, family)]
            def supported(key):
                values = {r.get(key, '').lower() for r in support}
                return 'true' in values and 'false' not in values
            matches = architecture_matches(hits, rule)
            if rule.get('requires_family_evidence') and not any(supported(k) for k in ('similarity','profile','reciprocal','phylogeny')):
                continue
            if not matches and not supported('similarity') and not supported('profile'):
                continue
            length_ok = len(sequence) >= rule['min_length']
            motif_sequence = sequence[-rule['motif_c_terminal_window']:] if rule.get('motif_c_terminal_window') else sequence
            motifs = all(re.search(pattern, motif_sequence) for pattern in rule.get('motifs', []))
            observations=scan_motifs(sequence,hits,rule,motif_library)
            extended=noncatalytic_alignment(support,hits,rule['noncatalytic_alignment']) if rule.get('noncatalytic_alignment') else True
            topology = not rule.get('requires_topology', False) or supported('topology')
            extra = any(supported(key) for key in ('reciprocal','phylogeny','topology','motif'))
            required = length_ok and motifs and topology
            homology_support = [r for r in support if r.get('similarity','').lower() == 'true']
            anchor_only = bool(homology_support) and all(r.get('reference_tier') == 'computational_anchor' for r in homology_support)
            # Missing diagnostic domains can be rescued only by independently supplied phylogenetic support.
            accepted = not rule.get('screening_only',False) and required and supported('similarity') and supported('reciprocal') and (matches or supported('phylogeny')) and (extended or supported('phylogeny')) and (not anchor_only or supported('phylogeny') or supported('profile'))
            if not length_ok:
                status, reason = 'excluded', 'below_family_minimum_length'
            elif accepted:
                status, reason = 'supported_candidate', 'convergent_evidence'
            elif (matches and rule.get('diagnostic', False)) or extra:
                status, reason = 'provisional', 'additional_validation_required'
            else:
                status, reason = 'excluded', 'generic_domain_or_similarity_only'
            results.append(dict(protein_id=identifier, family=family, module=rule['module'], status=status,
                reason=reason, length=len(sequence), minimum_length=rule['min_length'], domain_match=matches,
                motif_check=bool(motifs), topology_check=topology, reciprocal=supported('reciprocal'),
                motif_details=';'.join(f"{r['motif']}:{r['state']}:{r['start']}-{r['end']}" for r in observations),
                architecture=architecture, evidence_sources=';'.join(sorted({r['source'] for r in support}))))
    # Conflicting supported family assignments remain provisional, never silently pick one.
    by_protein = defaultdict(list)
    for row in results:
        if row['status'] == 'supported_candidate':
            by_protein[row['protein_id']].append(row)
    for rows in by_protein.values():
        if len(rows) > 1:
            for row in rows:
                row.update(status='provisional', reason='conflicting_family_assignments')
    output = Path(args.out)
    output.mkdir(parents=True, exist_ok=True)
    write_table(output / 'candidates.csv', results, ['protein_id','family','module','status','reason','length',
        'minimum_length','domain_match','motif_check','motif_details','topology_check','reciprocal','architecture','evidence_sources'])
    write_table(output / 'architectures.csv', architecture_rows, ['protein_id','length','architecture'])
    selected = {r['protein_id'] for r in results if r['status'] != 'excluded'}
    write_fasta(output / 'candidates.faa', [(key, sequences[key]) for key in sorted(selected)])
    save_json(output / 'classification.json', dict(rule_status=config['status'],
        counts=dict(Counter(r['status'] for r in results)), checksums={str(p):digest(p) for p in
        [args.proteome,args.domains,args.rules,args.homology,*(args.evidence or []),motif_path] if p}))

def classify(args):
    if getattr(args, 'classification_mode', 'comparative') == 'strict':
        return classify_strict(args)
    from adhesome_classification import classify_comparative
    return classify_comparative(args)


def discover_proteomes(root, settings, excluded=None):
    """Locate one reference FASTA per configured species without merging releases."""
    root = Path(root).resolve()
    excluded = Path(excluded).resolve() if excluded else None
    if not root.is_dir():
        raise ValueError(f'Input directory does not exist: {root}')
    matches = defaultdict(list)
    ignored = []
    for path in sorted(root.rglob('*')):
        if not path.is_file() or path.suffix.lower() not in {'.fa', '.faa', '.fasta', '.fas', '.fsa'}:
            continue
        if excluded and (path.resolve() == excluded or excluded in path.resolve().parents):
            continue
        # Match filenames and folders relative to the input root, never unrelated ancestors.
        label = re.sub(r'[^a-z0-9]', '', str(path.relative_to(root)).lower())
        # Explicit filenames/codes take precedence over broad species-name hints.
        candidates = [species['code'] for species in settings['proteomes']
                      if path.name.lower() == Path(species['path']).name.lower()
                      or path.stem.lower() == species['code'].lower()]
        if not candidates:
            candidates = [species['code'] for species in settings['proteomes']
                          if species['species'].split()[-1].lower() in label]
        if len(candidates) > 1:
            raise ValueError(f'Ambiguous species assignment: {path}')
        if candidates:
            matches[candidates[0]].append(path.resolve())
        else:
            ignored.append(str(path))
    for species in settings['proteomes']:
        paths = matches[species['code']]
        if len(paths) != 1:
            raise ValueError(f'{species["species"]}: expected exactly one FASTA, found {len(paths)}: {paths}')
    return {code: paths[0] for code, paths in matches.items()}, ignored


def batch(args):
    settings = json.loads(Path(args.config).read_text(encoding='utf-8'))
    output = Path(args.out).resolve()
    inputs, ignored = discover_proteomes(args.input_dir, settings, output)
    if not args.prepare_only:
        if not args.library:
            raise ValueError('Batch scanning requires --library; use --prepare-only for preprocessing alone')
        for program in ('hmmpress', 'hmmscan'):
            if not shutil.which(program):
                raise RuntimeError(f'{program} is not installed/on PATH; no batch scans started')
        if not Path(args.library + '.metadata.json').is_file():
            raise ValueError('Build the HMM library and metadata before batch scanning')
    output.mkdir(parents=True, exist_ok=True)
    summary = dict(input_root=str(Path(args.input_dir).resolve()), ignored_fastas=ignored,
                   mode='prepare_only' if args.prepare_only else 'prepare_scan_classify', species={})
    for species in settings['proteomes']:
        code = species['code']
        folder = output / code
        folder.mkdir(parents=True, exist_ok=True)
        entry = dict(species)
        entry['path'] = str(inputs[code])
        config = folder / 'proteome.json'
        save_json(config, dict(qc=settings['qc'], proteomes=[entry]))
        status = dict(species=species['species'], input=str(inputs[code]), output=str(folder), status='running')
        summary['species'][code] = status
        save_json(output / 'batch.json', summary)
        print(f'Processing {species["species"]}: {inputs[code]}', flush=True)
        try:
            prepare(config, folder / 'prepared')
            if not args.prepare_only:
                scan_args = argparse.Namespace(**vars(args))
                scan_args.proteome = str(folder / 'prepared' / (code + '.faa'))
                scan_args.out = str(folder / 'domains')
                scan(scan_args)
                homology_path = None
                if args.reference_set:
                    references = Path(args.reference_set)
                    forward_name='search_references.faa' if (references/'search_references.faa').exists() else 'adhesome_references.faa'
                    catalogue_name='search_catalogue.csv' if forward_name=='search_references.faa' else 'reference_catalogue.csv'
                    homology(argparse.Namespace(proteome=scan_args.proteome,references=str(references/forward_name),
                        catalogue=str(references/catalogue_name),metazoan=str(references/'metazoan_references.faa'),
                        metazoan_catalogue=str(references/'metazoan_catalogue.csv'),out=str(folder/'homology'),cpu=args.cpu,
                        evalue=1e-5,coverage=.3,tie_fraction=.01))
                    homology_path = str(folder/'homology/homology.csv')
                classify(argparse.Namespace(proteome=scan_args.proteome,
                    domains=str(folder / 'domains/domains.csv'), rules=args.rules,
                    homology=homology_path, evidence=None, manifest=args.manifest,
                    classification_mode=args.classification_mode, out=str(folder / 'candidates')))
            status['status'] = ('prepared' if args.prepare_only else
                                'completed_with_homology' if args.reference_set else 'completed_domain_screen')
        except Exception as error:
            status.update(status='failed', error=str(error))
            save_json(output / 'batch.json', summary)
            raise
        save_json(output / 'batch.json', summary)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='stage', required=True)
    p = sub.add_parser('batch', help='Recursively discover and process all three species into separate folders')
    p.add_argument('--input-dir', default='proteome_fasta_files')
    p.add_argument('--config', default=str(CONFIG / 'proteomes.json'))
    p.add_argument('--out', default='results/species')
    p.add_argument('--library')
    p.add_argument('--prepare-only', action='store_true')
    p.add_argument('--reference-set', help='Paired reference directory; also run BLASTP and integrate its evidence')
    p.add_argument('--classification-mode', choices=['comparative','strict'], default='comparative')
    p.add_argument('--manifest', default=str(CONFIG / 'adhesome_pfam_manifest.csv'))
    p.add_argument('--rules', default=str(CONFIG / 'family_rules.json'))
    p.add_argument('--cpu', type=int, default=4)
    p.add_argument('--max-i', type=float, default=1e-3)
    p.add_argument('--max-c', type=float, default=1e-3)
    p.add_argument('--min-coverage', type=float, default=0.35)
    p = sub.add_parser('prepare')
    p.add_argument('--config', default=str(CONFIG / 'proteomes.json'))
    p.add_argument('--out', default='results/prepared')
    p = sub.add_parser('scan')
    p.add_argument('--proteome', required=True)
    p.add_argument('--library', required=True)
    p.add_argument('--manifest', default=str(CONFIG / 'adhesome_pfam_manifest.csv'))
    p.add_argument('--out', required=True)
    p.add_argument('--cpu', type=int, default=4)
    p.add_argument('--max-i', type=float, default=1e-3)
    p.add_argument('--max-c', type=float, default=1e-3)
    p.add_argument('--min-coverage', type=float, default=0.35)
    p = sub.add_parser('homology')
    for name in ('proteome','out'):
        p.add_argument('--' + name, required=True)
    for name, filename in [('references','search_references.faa'),('catalogue','search_catalogue.csv'),
                            ('metazoan','metazoan_references.faa'),('metazoan-catalogue','metazoan_catalogue.csv')]:
        p.add_argument('--'+name, default=str(REFERENCE_SET/filename))
    p.add_argument('--cpu', type=int, default=4)
    p.add_argument('--evalue', type=float, default=1e-5)
    p.add_argument('--coverage', type=float, default=0.3)
    p.add_argument('--tie-fraction', type=float, default=0.01)
    p = sub.add_parser('profiles')
    for name in ('proteome','profiles','profile-map','out'):
        p.add_argument('--' + name, required=True)
    p.add_argument('--cpu', type=int, default=4)
    p.add_argument('--max-i', type=float, default=1e-3)
    p.add_argument('--max-c', type=float, default=1e-3)
    p.add_argument('--min-coverage', type=float, default=0.35)
    p = sub.add_parser('classify')
    p.add_argument('--classification-mode', choices=['comparative','strict'], default='comparative')
    p.add_argument('--manifest', default=str(CONFIG / 'adhesome_pfam_manifest.csv'))
    p.add_argument('--proteome', required=True)
    p.add_argument('--domains', required=True)
    p.add_argument('--rules', default=str(CONFIG / 'family_rules.json'))
    p.add_argument('--homology')
    p.add_argument('--evidence', action='append', help='Repeat for profile/topology/phylogeny evidence files')
    p.add_argument('--out', required=True)
    args = parser.parse_args()
    for key in ('coverage','min_coverage','tie_fraction'):
        value = getattr(args, key, None)
        if value is not None and not 0 <= value < 1:
            parser.error(f'{key} must be >=0 and <1')
    for key in ('max_i','max_c','evalue','cpu'):
        value = getattr(args, key, None)
        if value is not None and not 0 < value < float('inf'):
            parser.error(f'{key} must be positive and finite')
    if args.stage == 'prepare':
        prepare(args.config, args.out)
    else:
        globals()[args.stage](args)

if __name__ == '__main__':
    main()
