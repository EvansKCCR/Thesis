#!/usr/bin/env python3
"""Auditable preparation, domain scanning, homology recovery and evidence integration."""

# Resolve source resources independently of the caller working directory.
import sys as _layout_sys
from pathlib import Path as _LayoutPath
_LAYOUT_BUNDLE_ROOT = next(p for p in _LayoutPath(__file__).resolve().parents
                          if (p/'config/layout_paths.json').is_file() and (p/'pipeline_paths.py').is_file())
if str(_LAYOUT_BUNDLE_ROOT) not in _layout_sys.path:
    _layout_sys.path.insert(0, str(_LAYOUT_BUNDLE_ROOT))
from pipeline_paths import activate as _layout_activate, resource_path as _layout_resource, default_thread_count as _default_thread_count
_layout_activate()
_ORTHOLOGY_SCRIPTS = _LAYOUT_BUNDLE_ROOT / 'resource_library' / 'orthology_phylogeny' / 'scripts'
if str(_ORTHOLOGY_SCRIPTS) not in _layout_sys.path:
    _layout_sys.path.insert(0, str(_ORTHOLOGY_SCRIPTS))

import argparse
from collections import Counter, defaultdict
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from adhesome_core import fasta, write_fasta, table, write_table, digest, save_json, pf, sequence_dict, VERSION, detection_evidence, validate_protein_input
from adhesome_components import annotate_components, load_crosswalk, CANDIDATE_FIELDS
from reference_identity import identity_evidence, IDENTITY_FIELDS
from parse_adhesome_domtblout import parse, FIELDS

PROJECT_ROOT = _layout_resource(Path(__file__)).resolve().parent
CONFIG = _layout_resource(PROJECT_ROOT / 'config')
REFERENCE_SET = _layout_resource(PROJECT_ROOT / 'references' / 'adhesome_v1.0')

MODULES = {'extracellular_matrix', 'integrin_receptors', 'membrane_proximal_adaptors',
           'actin_linking', 'focal_adhesion_scaffolds', 'downstream_signalling','other_LIM'}

from prepare_reference import prepare as prepare_gene_level, infer_gene
from pipeline import read_tsv

def prepare(config, output):
    output = _layout_resource(Path(output))
    settings = json.loads(_layout_resource(Path(config)).read_text(encoding='utf-8'))
    qc = settings['qc']
    if not (qc['min_length'] > 0 and 0 <= qc['max_ambiguous_fraction'] <= 1 and qc['max_ambiguous_run'] >= 0):
        raise ValueError('Invalid sequence QC thresholds')
    codes = [s['code'] for s in settings['proteomes']]
    if len(codes) != len(set(codes)) or any(not re.fullmatch(r'[A-Za-z0-9_]+', c) for c in codes):
        raise ValueError('Species codes must be unique safe identifiers')
    config_path=_layout_resource(Path(config)).resolve()
    for species in settings['proteomes']:
        validate_protein_input(config_path.parent/species['path'])
    output.mkdir(parents=True, exist_ok=True)
    summaries = {}
    for species in settings['proteomes']:
        source = config_path.parent / species['path']
        mapping=species.get('protein_gene_tsv')
        mapping_path=(config_path.parent/mapping).resolve() if mapping else None
        if mapping_path is not None and not mapping_path.is_file():
            raise ValueError(f'{species["code"]}: configured protein_gene.tsv not found: {mapping_path}')
        prep_dir=output/'_sequence_qc'/species['code']
        prepared=prepare_gene_level(species['code'],source,mapping_path,prep_dir,
            internal_stop_policy='mask',gene_regex=species.get('gene_regex'),qc=qc,
            cdhit=species.get('cdhit_executable','cd-hit'))
        prefix=output/species['code']
        shutil.copyfile(prepared,str(prefix)+'.faa')
        _, map_rows = read_tsv(prep_dir/'identifier_map.tsv')
        mapped=[]
        for row in map_rows:
            mapped.append(dict(protein_id=row['prepared_id'],original_id=row['original_id'],gene_id=row['gene_id'],
                length=row['length'],ambiguous_fraction=row['ambiguous_fraction'],
                longest_ambiguous_run=row['longest_ambiguous_run'],decision=row['decision'],original_header=row['original_header']))
        write_table(str(prefix)+'.id_map.csv',mapped,
                    ['protein_id','original_id','gene_id','length','ambiguous_fraction','longest_ambiguous_run','decision','original_header'])
        prep_meta=json.loads((prep_dir/'preparation.json').read_text(encoding='utf-8'))
        summaries[species['code']] = dict(species=species['species'], assembly=species.get('assembly','not_supplied'),
            assembly_status='user_supplied_not_verified_from_FASTA', source_url=species.get('source_url'),
            input=str(source), input_sha256=digest(source), input_proteins=prep_meta['input_protein_ids'],
            genes=prep_meta['genes'],retained=prep_meta['sequence_qc']['retained'],
            excluded=prep_meta['sequence_qc']['excluded'],sequence_qc=prep_meta['sequence_qc'],
            output_sha256=digest(str(prefix)+'.faa'))
    save_json(_layout_resource(output / 'preparation.json'), dict(pipeline_version=VERSION,config=settings, config_sha256=digest(config), species=summaries))
    print(json.dumps(summaries, indent=2))

def run(command, log):
    if not shutil.which(command[0]):
        raise RuntimeError(f'{command[0]} is not installed/on PATH. See docs/MANUAL.md in the bundle root.')
    executable = shutil.which(command[0])
    flag = '-version' if command[0] == 'blastp' else '-h'
    version = subprocess.run([executable, flag], capture_output=True, text=True, check=False)
    save_json(str(log) + '.command.json', dict(command=[str(x) for x in command],
        executable=executable, executable_sha256=digest(executable),
        version_output=(version.stdout + version.stderr)[:8000]))
    started=time.perf_counter()
    print(f'{datetime.now(timezone.utc).isoformat()} | TOOL START {command[0]}: {log}',flush=True)
    try:
        with open(_layout_resource(log), 'w', encoding='utf-8') as handle:
            subprocess.run([str(x) for x in command], stdout=handle, stderr=subprocess.STDOUT, check=True)
    except Exception as exc:
        print(f'{datetime.now(timezone.utc).isoformat()} | TOOL FAILED {command[0]}; elapsed {time.perf_counter()-started:.1f}s; {exc}',flush=True)
        raise
    print(f'{datetime.now(timezone.utc).isoformat()} | TOOL COMPLETE {command[0]}; elapsed {time.perf_counter()-started:.1f}s',flush=True)


def timed_stage(label, action):
    started=time.perf_counter();stamp=datetime.now(timezone.utc).isoformat()
    print(f'{stamp} | STAGE {label} START',flush=True)
    try:result=action()
    except Exception as exc:
        print(f'{datetime.now(timezone.utc).isoformat()} | STAGE {label} FAILED; elapsed {time.perf_counter()-started:.1f}s; {exc}',flush=True)
        raise
    print(f'{datetime.now(timezone.utc).isoformat()} | STAGE {label} COMPLETE; elapsed {time.perf_counter()-started:.1f}s',flush=True)
    return result

def scan(args):
    validate_protein_input(args.proteome,hmmer_limit=True)
    output = _layout_resource(Path(args.out))
    output.mkdir(parents=True, exist_ok=True)
    metadata = json.loads(_layout_resource(Path(args.library + '.metadata.json')).read_text(encoding='utf-8'))
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
        if not _layout_resource(Path(library)).stat().st_size:
            continue
        if os.environ.get('SCHISTO_ISOLATED_BRANCH'):
            # hmmpress writes sidecar indexes. Each parallel gate needs its own
            # writable runtime index, never concurrent writes to the frozen HMMs.
            local_library=output/(kind+'.library.hmm')
            if not _layout_resource(local_library).exists() or digest(local_library)!=digest(library):
                shutil.copyfile(library,local_library)
            library=str(local_library)
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
    combined = _layout_resource(output / 'combined.domtblout')
    combined.write_text(''.join(_layout_resource(p).read_text(encoding='utf-8') for p in parts), encoding='utf-8')
    rows = parse(combined, mapping, models, args.max_i, args.max_c, args.min_coverage)
    write_table(_layout_resource(output / 'domains.csv'), [r for r in rows if r['decision'] == 'retained'], FIELDS)
    write_table(_layout_resource(output / 'domains.audit.csv'), rows, FIELDS)
    save_json(_layout_resource(output / 'scan.json'), dict(pipeline_version=VERSION,commands=commands, proteome_sha256=digest(args.proteome),
        metadata_sha256=digest(args.library + '.metadata.json'), thresholds=dict(i=args.max_i,c=args.max_c,
        profile_coverage=args.min_coverage), search_space=len(models), python=sys.version,
        artifact_checksums={name:digest(output/name) for name in ('domains.csv','domains.audit.csv','combined.domtblout')}))

def catalogue(path, sequence_path, rules=None):
    rows = table(path)
    required = ['reference_id','family','module','species','accession','evidence_source']
    seqs = sequence_dict(sequence_path)
    if not rows:
        raise ValueError('Reference catalogue is empty; supply sequence-verified family references')
    found = {}
    for row in rows:
        if any(not row.get(k, '').strip() for k in required):
            raise ValueError('Reference catalogue has missing required metadata')
        key = row['reference_id']
        component_only = row['family']=='component_reference' and row.get('role')=='crosswalk_component' and bool(row.get('component_symbol'))
        if key not in seqs or key in found or (row['module'] not in MODULES and not component_only):
            raise ValueError(f'Invalid/duplicate reference or module: {key}')
        if rules is not None and row['family'] not in rules and not component_only:
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
    with open(_layout_resource(path), encoding='utf-8') as handle:
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


def reciprocal_agreement(reference, top_rows):
    """Compete by accession, allowing several hierarchical labels on one sequence.

    Component-only matches remain competitors for family queries unless that
    accession has an explicitly configured family label in this same panel.
    """
    by_accession = defaultdict(list)
    for row in top_rows:
        by_accession[row.get('accession') or row['reference_id']].append(row)
    if reference['family'] == 'component_reference':
        symbol = reference.get('component_symbol')
        return bool(by_accession) and all(symbol and any(row.get('component_symbol')==symbol for row in rows)
                                         for rows in by_accession.values())
    return bool(by_accession) and all(any(row['family']==reference['family'] for row in rows)
                                     for rows in by_accession.values())

def homology(args):
    refs = catalogue(args.catalogue, args.references)
    reference_sequences = sequence_dict(args.references)
    broad_rows = table(args.metazoan_catalogue)
    broad = {r['reference_id']: r['family'] for r in broad_rows}
    broad_metadata = {r['reference_id']: r for r in broad_rows}
    broad_sequences = sequence_dict(args.metazoan)
    broad_ids = set(broad_sequences)
    if len(broad) != len(broad_rows) or set(broad) != broad_ids:
        raise ValueError('Broad metazoan catalogue/FASTA mismatch or duplicate IDs')
    if not any(not family or family == 'non_adhesome' for family in broad.values()):
        raise ValueError('Reciprocal database must include annotated non-adhesome competitors')
    output = _layout_resource(Path(args.out))
    output.mkdir(parents=True, exist_ok=True)
    import hashlib
    for row in broad_rows:
        if row.get('sequence_sha256') and hashlib.sha256(broad_sequences[row['reference_id']].encode()).hexdigest() != row['sequence_sha256']:
            raise ValueError(f'Reciprocal catalogue/sequence checksum mismatch: {row["reference_id"]}')
    snapshot = _layout_resource(output / 'reference_snapshot')
    snapshot.mkdir(exist_ok=True)
    paths = {'search_references.faa':args.references,'search_catalogue.csv':args.catalogue,
             'metazoan_references.faa':args.metazoan,'metazoan_catalogue.csv':args.metazoan_catalogue}
    lock_path = _layout_resource(Path(args.references).resolve().parent / 'reference_lock.json')
    if _layout_resource(lock_path).exists():
        lock = json.loads(_layout_resource(lock_path).read_text(encoding='utf-8'))
        for name, source in paths.items():
            if lock['files'].get(_layout_resource(Path(source)).name) != digest(source):
                raise ValueError(f'Reference lock mismatch: {name}; rebuild or select a consistent reference set')
    # Every run keeps the exact search sequences, not just a pointer to mutable databases.
    for name, source in paths.items():
        destination = snapshot / name
        if _layout_resource(destination).exists() and digest(destination) != digest(source):
            raise ValueError('Reference snapshot differs; choose a new output directory for this reference release')
        if _layout_resource(Path(source)).resolve() != _layout_resource(destination).resolve():
            shutil.copy2(source,destination)
    for name in ('reference_lock.json','experimental_annotations.csv','resolution_audit.json',
                 'component_catalogue.csv','component_family_catalogue.csv','component_coverage.csv','component_expansion.json'):
        source = _layout_resource(Path(args.references)).resolve().parent / name
        if _layout_resource(source).exists():
            shutil.copy2(source,snapshot/name)
    commands = []
    for label, query, subject in [('forward', str(_layout_resource(snapshot/'search_references.faa')), args.proteome),
                                  ('reverse', args.proteome, str(_layout_resource(snapshot/'metazoan_references.faa')))]:
        command = ['blastp', '-query', query, '-subject', subject, '-evalue', str(args.evalue),
                   '-num_threads', str(args.cpu), '-max_target_seqs', '1000000', '-max_hsps', '1',
                   '-outfmt', '6 ' + BLAST_FIELDS, '-out', str(output / (label + '.tsv'))]
        commands.append(command)
        run(command, output / (label + '.log'))
    forward = blast_hits(_layout_resource(output / 'forward.tsv'), args.evalue, args.coverage)
    reverse = blast_hits(_layout_resource(output / 'reverse.tsv'), args.evalue, args.coverage)
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
        supported = reciprocal_agreement(refs[hit['query']], [broad_metadata[h['subject']] for h in top])
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
                         component_symbol=refs[hit['query']].get('component_symbol',''),
                         source=str(_layout_resource(snapshot/'search_catalogue.csv'))))
    write_table(_layout_resource(output / 'homology.csv'), rows, ['protein_id','family','reference_id','similarity','reciprocal',
                'evalue','bitscore','query_coverage','subject_coverage','subject_start','subject_end','reverse_families','reverse_reference_ids',
                'reference_accession','reference_species','reference_evidence_source','reference_evidence_status','reference_tier',
                'reference_phylogenetic_group','reverse_reference_tiers','reference_sequence_sha256','component_symbol','source'])
    matched_ids = sorted({r['reference_id'] for r in rows})
    write_fasta(_layout_resource(output/'matched_reference_sequences.faa'),[(rid,reference_sequences[rid]) for rid in matched_ids])
    write_table(_layout_resource(output/'matched_reference_catalogue.csv'),[refs[rid] for rid in matched_ids],list(next(iter(refs.values()))))
    save_json(_layout_resource(output / 'homology.json'), dict(pipeline_version=VERSION,commands=commands, evalue=args.evalue, coverage=args.coverage,
        artifact_checksums={name:digest(output/name) for name in ('forward.tsv','reverse.tsv','homology.csv')},
        tie_fraction=args.tie_fraction, reference_snapshot=str(snapshot),snapshot_source_mapping={name:str(p) for name,p in paths.items()},
        snapshot_checksums={p.name:digest(p) for p in snapshot.iterdir() if _layout_resource(p).is_file()},
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
    output = _layout_resource(Path(args.out))
    output.mkdir(parents=True, exist_ok=True)
    command = ['hmmsearch', '--cpu', str(args.cpu), '--noali', '-E', '100', '--domE', '100',
               '--domtblout', str(_layout_resource(output / 'profiles.domtblout')), args.profiles, args.proteome]
    run(command, _layout_resource(output / 'hmmsearch.log'))
    rows = []
    with open(_layout_resource(output / 'profiles.domtblout'), encoding='utf-8') as handle:
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
    write_table(_layout_resource(output / 'profile_evidence.csv'), rows, ['protein_id','family','profile','source'])
    save_json(_layout_resource(output / 'profiles.json'), dict(pipeline_version=VERSION,command=command, i=args.max_i,c=args.max_c,
        coverage=args.min_coverage, checksums={str(p):digest(p) for p in [args.proteome,args.profiles,args.profile_map]}))

def classify_strict(args):
    config = json.loads(_layout_resource(Path(args.rules)).read_text(encoding='utf-8'))
    rules = config['families']
    from adhesome_motifs import load_library,scan_motifs,noncatalytic_alignment
    motif_path=(_layout_resource(Path(args.rules).resolve().parent/config['motif_library'])).resolve() if config.get('motif_library') else None
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
                detection_evidence(row)
                component_only = row['family']=='component_reference' and bool(row.get('component_symbol'))
                if row['protein_id'] not in sequences or (row['family'] not in rules and not component_only):
                    raise ValueError('Unknown evidence protein/family')
                if not row.get('source'):
                    raise ValueError('Evidence must include a source')
                evidence[(row['protein_id'], row['family'])].append(row)
    results, architecture_rows = [], []
    classification_started=time.perf_counter();total_sequences=len(sequences)
    for sequence_index,(identifier, sequence) in enumerate(sequences.items(),1):
        if sequence_index==1 or sequence_index%100==0 or sequence_index==total_sequences:
            elapsed=time.perf_counter()-classification_started
            eta=(elapsed/sequence_index)*(total_sequences-sequence_index) if sequence_index else 0
            print(f'{datetime.now(timezone.utc).isoformat()} | CLASSIFY {sequence_index}/{total_sequences} ({100*sequence_index/max(1,total_sequences):.1f}%); elapsed {elapsed:.1f}s; ETA {eta:.1f}s; protein {identifier}',flush=True)
        hits = sorted(domains[identifier], key=lambda r: (int(r['ali_from']), int(r['ali_to'])))
        architecture = ';'.join(f'{r["pfam_accession"]}:{r["ali_from"]}-{r["ali_to"]}' for r in hits)
        if hits:
            architecture_rows.append(dict(protein_id=identifier, length=len(sequence), architecture=architecture))
        for family, rule in rules.items():
            support = evidence[(identifier, family)]
            resolved_identity,identity_details=identity_evidence(support,rule)
            def supported(key):
                values = {r.get(key, '').lower() for r in support}
                return 'true' in values and 'false' not in values
            matches = architecture_matches(hits, rule)
            if rule.get('requires_family_evidence') and not any(supported(k) for k in ('similarity','profile','reciprocal')):
                continue
            if not matches and not supported('similarity') and not supported('profile'):
                continue
            length_ok = len(sequence) >= rule['min_length']
            motif_sequence = sequence[-rule['motif_c_terminal_window']:] if rule.get('motif_c_terminal_window') else sequence
            motifs = all(re.search(pattern, motif_sequence) for pattern in rule.get('motifs', []))
            observations=scan_motifs(sequence,hits,rule,motif_library)
            extended=noncatalytic_alignment(support,hits,rule['noncatalytic_alignment']) if rule.get('noncatalytic_alignment') else True
            extra = any(supported(key) for key in ('reciprocal','motif'))
            required = length_ok and motifs
            homology_support = [r for r in support if r.get('similarity','').lower() == 'true']
            anchor_only = bool(homology_support) and all(r.get('reference_tier') == 'computational_anchor' for r in homology_support)
            screening=bool(rule.get('screening_only') or rule.get('classification_level')=='screen')
            accepted = not screening and (not rule.get('requires_identity_evidence') or resolved_identity) and required and supported('similarity') and supported('reciprocal') and matches and extended and (not anchor_only or supported('profile'))
            if not length_ok:
                status, reason = 'excluded', 'below_family_minimum_length'
            elif accepted:
                status, reason = 'supported_candidate', 'convergent_evidence'
            elif (matches and rule.get('diagnostic', False)) or extra:
                status, reason = 'provisional', 'additional_validation_required'
            else:
                status, reason = 'excluded', 'generic_domain_or_similarity_only'
            results.append(dict(protein_id=identifier, family=family, module=rule['module'], status=status,
                **identity_details,
                reason=reason, length=len(sequence), minimum_length=rule['min_length'], domain_match=matches,
                motif_check=bool(motifs), reciprocal=supported('reciprocal'),
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
    output = _layout_resource(Path(args.out))
    output.mkdir(parents=True, exist_ok=True)
    fields=['protein_id','family','module','status','reason','length',*IDENTITY_FIELDS,
        'minimum_length','domain_match','motif_check','motif_details','reciprocal','architecture','evidence_sources']
    crosswalk=getattr(args,'crosswalk',None)
    if crosswalk:
        annotate_components(crosswalk,getattr(args,'manifest',None),[r for hits in domains.values() for r in hits],evidence,results,output)
        fields+=CANDIDATE_FIELDS
    write_table(_layout_resource(output / 'candidates.csv'), results, fields)
    write_table(_layout_resource(output / 'architectures.csv'), architecture_rows, ['protein_id','length','architecture'])
    selected = {r['protein_id'] for r in results if r['status'] != 'excluded'}
    write_fasta(_layout_resource(output / 'candidates.faa'), [(key, sequences[key]) for key in sorted(selected)])
    save_json(_layout_resource(output / 'classification.json'), dict(pipeline_version=VERSION,stage='detection',deferred_stages=['topology','localization','phylogeny'],rule_status=config['status'],
        counts=dict(Counter(r['status'] for r in results)), checksums={str(p):digest(p) for p in
        [args.proteome,args.domains,args.rules,args.homology,*(args.evidence or []),motif_path,crosswalk] if p}))

def classify(args):
    if getattr(args, 'classification_mode', 'comparative') == 'strict':
        return classify_strict(args)
    from adhesome_classification import classify_comparative
    return classify_comparative(args)


def discover_proteomes(root, settings, excluded=None):
    """Locate one reference FASTA per configured species without merging releases."""
    root = _layout_resource(Path(root)).resolve()
    excluded = _layout_resource(Path(excluded)).resolve() if excluded else None
    if not _layout_resource(root).is_dir():
        alternate=_layout_resource(root.parent.parent/'input')/root.name if root.parent.name=='inputs' else None
        hint=f'; use --input-dir {alternate}' if alternate and _layout_resource(alternate).is_dir() else ''
        raise ValueError(f'Input directory does not exist: {root}'+hint)
    matches = defaultdict(list)
    ignored = []
    for path in sorted(_layout_resource(root).rglob('*')):
        if not _layout_resource(path).is_file() or path.suffix.lower() not in {'.fa', '.faa', '.fasta', '.fas', '.fsa'}:
            continue
        if excluded and (_layout_resource(path).resolve() == excluded or excluded in _layout_resource(path).resolve().parents):
            continue
        # Match filenames and folders relative to the input root, never unrelated ancestors.
        label = re.sub(r'[^a-z0-9]', '', str(path.relative_to(root)).lower())
        # Explicit filenames/codes take precedence over broad species-name hints.
        candidates = [species['code'] for species in settings['proteomes']
                      if path.name.lower() == _layout_resource(Path(species['path'])).name.lower()
                      or path.stem.lower() == species['code'].lower()]
        if not candidates:
            candidates = [species['code'] for species in settings['proteomes']
                          if species['species'].split()[-1].lower() in label]
        if len(candidates) > 1:
            raise ValueError(f'Ambiguous species assignment: {path}')
        if candidates:
            matches[candidates[0]].append(_layout_resource(path).resolve())
        else:
            ignored.append(str(path))
    for species in settings['proteomes']:
        paths = matches[species['code']]
        if len(paths) != 1:
            raise ValueError(f'{species["species"]}: expected exactly one FASTA, found {len(paths)}: {paths}')
    return {code: paths[0] for code, paths in matches.items()}, ignored


def check_inputs(root,settings,excluded=None):
    """Check all input mappings before output creation or external analysis."""
    paths,ignored=discover_proteomes(root,settings,excluded)
    summary={}
    for entry in settings['proteomes']:
        code=entry['code'];path=paths[code]
        gene_regex=entry.get('gene_regex')
        if gene_regex:
            try: pattern=re.compile(gene_regex)
            except re.error as exc:raise ValueError(f'{code}: invalid gene_regex: {exc}')
            if pattern.groups<1:raise ValueError(f'{code}: gene_regex must capture a gene identifier')
        identifiers=set();genes=set();count=0
        with _layout_resource(path).open(encoding='utf-8-sig') as handle:
            for line in handle:
                line=line.strip()
                if not line.startswith('>'):continue
                header=line[1:].strip();count+=1
                if not header:raise ValueError(f'{code}: empty FASTA header in {path}')
                identifier=header.split()[0]
                identifiers.add(identifier)
                inferred,_,_=infer_gene(header,identifier,entry.get('gene_regex'))
                genes.add(inferred or identifier)
        if not count:raise ValueError(f'{code}: no FASTA records in {path}')
        validate_protein_input(path)
        summary[code]=dict(fasta=str(path),proteins=count,genes=len(genes))
    return paths,ignored,summary


def batch(args):
    settings = json.loads(_layout_resource(Path(args.config)).read_text(encoding='utf-8'))
    output = _layout_resource(Path(args.out)).resolve()
    resuming=bool(getattr(args,'resume',False))
    excluded=set(getattr(args,'exclude_target',None) or [])
    if resuming and not (_layout_resource(output/'batch.json')).is_file():
        raise ValueError('No detection batch checkpoint found for --resume: '+str(_layout_resource(output/'batch.json')))
    if excluded:
        settings=dict(settings,proteomes=[entry for entry in settings['proteomes'] if entry['code'] not in excluded])
    if not settings['proteomes']:raise ValueError('No active targets remain after exclusions')
    inputs, ignored, input_summary = check_inputs(args.input_dir,settings,output)
    requested=set(getattr(args,'species',None) or [])
    configured={entry['code'] for entry in settings['proteomes']}
    if requested-configured:raise ValueError('Unknown configured targets: '+', '.join(sorted(requested-configured)))
    if getattr(args,'check_inputs',False):
        print(json.dumps(dict(input_root=str(_layout_resource(Path(args.input_dir)).resolve()),species=input_summary,ignored_fastas=ignored),indent=2))
        return
    crosswalk=getattr(args,'crosswalk',None)
    if crosswalk:
        load_crosswalk(crosswalk,args.manifest)
    if not args.prepare_only:
        if not args.library:
            raise ValueError('Batch scanning requires --library; use --prepare-only for preprocessing alone')
        for program in ('hmmpress', 'hmmscan'):
            if not shutil.which(program):
                raise RuntimeError(f'{program} is not installed/on PATH; no batch scans started')
        if not _layout_resource(Path(args.library + '.metadata.json')).is_file():
            raise ValueError('Build the HMM library and metadata before batch scanning')
    summary = dict(pipeline_version=VERSION,stage='detection',deferred_stages=['topology','localization','phylogeny'],input_root=str(_layout_resource(Path(args.input_dir)).resolve()), ignored_fastas=ignored,
                   mode='prepare_only' if args.prepare_only else 'prepare_scan_classify', species={})
    if (requested or resuming) and (_layout_resource(output/'batch.json')).exists():
        previous=json.loads((_layout_resource(output/'batch.json')).read_text(encoding='utf-8'))
        if previous.get('input_root')!=summary['input_root'] or previous.get('mode')!=summary['mode']:
            raise ValueError('Selected-target recovery requires the same input directory and analysis mode')
        removed=set(previous.get('species',{}))-configured
        complete_states={'completed_with_homology','completed_domain_screen','prepared'}
        if removed and (not resuming or any(previous['species'][code]['status'] in complete_states for code in removed)):
            raise ValueError('Current target config omits completed targets; retain them when resuming this run')
        if resuming:
            pending={code for code in configured if previous.get('species',{}).get(code,{}).get('status') not in complete_states}
            requested=(requested & pending) if requested else pending
        for entry in settings['proteomes']:
            code=entry['code']
            if code in requested or code not in previous.get('species',{}):continue
            if resuming and previous['species'][code]['status'] not in complete_states:continue
            stored=_layout_resource(output/code/'proteome.json')
            if not _layout_resource(stored).is_file():raise ValueError(f'{code}: prior configuration is missing; use a new run')
            old=json.loads(_layout_resource(stored).read_text(encoding='utf-8'))
            expected=dict(entry,path=str(inputs[code]))
            if old.get('qc')!=settings['qc'] or old.get('proteomes')!=[expected]:
                raise ValueError(f'{code}: unselected target settings changed; include it in --species or use a new run')
            provenance=_layout_resource(output/code/'prepared/preparation.json')
            record=json.loads(_layout_resource(provenance).read_text(encoding='utf-8'))['species'][code]
            if record['input_sha256']!=digest(inputs[code]) or record['output_sha256']!=digest(_layout_resource(output/code/'prepared'/f'{code}.faa')):
                raise ValueError(f'{code}: preserved input/prepared FASTA changed; include it in --species or use a new run')
            if not args.prepare_only:
                scan_meta=json.loads((_layout_resource(output/code/'domains/scan.json')).read_text(encoding='utf-8'))
                thresholds=dict(i=args.max_i,c=args.max_c,profile_coverage=args.min_coverage)
                if scan_meta['thresholds']!=thresholds or scan_meta['metadata_sha256']!=digest(args.library+'.metadata.json'):
                    raise ValueError(f'{code}: scan settings/library changed; use a new run')
                classification=json.loads((_layout_resource(output/code/'candidates/classification.json')).read_text(encoding='utf-8'))
                if classification.get('mode','comparative')!=args.classification_mode:
                    raise ValueError(f'{code}: classification mode changed; use a new run')
                for name,sha in classification['checksums'].items():
                    if digest(name)!=sha:raise ValueError(f'{code}: preserved classification input changed: {name}')
                for resource in (args.rules,args.manifest,crosswalk):
                    if resource and classification['checksums'].get(str(_layout_resource(Path(resource)).resolve()))!=digest(resource):
                        raise ValueError(f'{code}: classification resources changed; use a new run')
        summary['species']={code:record for code,record in previous['species'].items() if code in configured}
        if resuming:
            summary['excluded_targets']=dict(previous.get('excluded_targets',{}))
            for code in sorted(removed):
                summary['excluded_targets'][code]=dict(previous['species'][code],
                    exclusion_reason='user_requested_exclusion' if code in excluded else 'removed_from_target_configuration')
        summary['selected_target_attempts']=previous.get('selected_target_attempts',[])+[dict(targets=sorted(requested),previous_status={code:previous['species'].get(code,{}).get('status') for code in requested})]
    output.mkdir(parents=True, exist_ok=True)
    for target_index,species in enumerate(settings['proteomes'],1):
        code = species['code']
        if (requested or resuming) and code not in requested:
            if resuming:print(f'REUSE completed target {code}',flush=True)
            continue
        folder = output / code
        if not resuming and requested and summary['species'].get(code,{}).get('status')=='failed' and _layout_resource(folder).exists():
            archives=_layout_resource(output/'archived_attempts')
            if _layout_resource(folder).resolve().parent!=output or _layout_resource(archives).resolve().parent!=output:
                raise ValueError('Failed-attempt archive must remain inside the detection run')
            archives.mkdir(exist_ok=True)
            number=1
            while (_layout_resource(archives/f'{code}_{number:03}')).exists():number+=1
            archive=_layout_resource(archives/f'{code}_{number:03}')
            if _layout_resource(archive).resolve().parent!=_layout_resource(archives).resolve():raise ValueError('Unsafe failed-attempt archive path')
            shutil.move(str(folder),str(archive))
        folder.mkdir(parents=True, exist_ok=True)
        entry = dict(species)
        entry['path'] = str(inputs[code])
        config = _layout_resource(folder / 'proteome.json')
        target_settings=dict(qc=settings['qc'],proteomes=[entry])
        old_preparation=_layout_resource(folder/'prepared/preparation.json')
        reuse_preparation=False
        if resuming and _layout_resource(config).is_file() and _layout_resource(old_preparation).is_file():
            prior=json.loads(_layout_resource(old_preparation).read_text(encoding='utf-8')).get('species',{}).get(code,{})
            reuse_preparation=(json.loads(_layout_resource(config).read_text(encoding='utf-8'))==target_settings
                and prior.get('input_sha256')==digest(inputs[code])
                and (_layout_resource(folder/'prepared'/f'{code}.faa')).is_file()
                and prior.get('output_sha256')==digest(_layout_resource(folder/'prepared'/f'{code}.faa')))
        if not reuse_preparation:save_json(config,target_settings)
        status = dict(species=species['species'], input=str(inputs[code]), output=str(folder), status='running')
        summary['species'][code] = status
        save_json(_layout_resource(output / 'batch.json'), summary)
        print(f'{datetime.now(timezone.utc).isoformat()} | TARGET {target_index}/{len(settings["proteomes"])} START {species["species"]}: {inputs[code]}', flush=True)
        try:
            if reuse_preparation:print(f'REUSE {code} preparation',flush=True)
            else:timed_stage(f'{code} sequence preparation',lambda:prepare(config, _layout_resource(folder / 'prepared')))
            if not args.prepare_only:
                scan_args = argparse.Namespace(**vars(args))
                scan_args.proteome = str(_layout_resource(folder / 'prepared') / (code + '.faa'))
                scan_args.out = str(_layout_resource(folder / 'domains'))
                scan_receipt=_layout_resource(folder/'domains/scan.json')
                reuse_scan=False
                if resuming and _layout_resource(scan_receipt).is_file():
                    prior=json.loads(_layout_resource(scan_receipt).read_text(encoding='utf-8'))
                    artifacts=prior.get('artifact_checksums',{})
                    reuse_scan=(prior.get('proteome_sha256')==digest(scan_args.proteome)
                        and prior.get('metadata_sha256')==digest(args.library+'.metadata.json')
                        and prior.get('thresholds')==dict(i=args.max_i,c=args.max_c,profile_coverage=args.min_coverage)
                        and set(artifacts)=={'domains.csv','domains.audit.csv','combined.domtblout'}
                        and all((_layout_resource(folder/'domains'/name)).is_file() and digest(_layout_resource(folder/'domains')/name)==sha for name,sha in artifacts.items()))
                if reuse_scan:print(f'REUSE {code} completed domain scan',flush=True)
                else:timed_stage(f'{code} Pfam domain scan',lambda:scan(scan_args))
                homology_path = None
                if args.reference_set:
                    references = _layout_resource(Path(args.reference_set))
                    forward_name='search_references.faa' if (_layout_resource(references/'search_references.faa')).exists() else 'adhesome_references.faa'
                    catalogue_name='search_catalogue.csv' if forward_name=='search_references.faa' else 'reference_catalogue.csv'
                    timed_stage(f'{code} reciprocal homology search',lambda:homology(argparse.Namespace(proteome=scan_args.proteome,references=str(references/forward_name),
                        catalogue=str(references/catalogue_name),metazoan=str(_layout_resource(references/'metazoan_references.faa')),
                        metazoan_catalogue=str(_layout_resource(references/'metazoan_catalogue.csv')),out=str(_layout_resource(folder/'homology')),cpu=args.cpu,
                        evalue=1e-5,coverage=.3,tie_fraction=.01)))
                    homology_path = str(_layout_resource(folder/'homology/homology.csv'))
                timed_stage(f'{code} family classification',lambda:classify(argparse.Namespace(proteome=scan_args.proteome,
                    domains=str(_layout_resource(folder / 'domains/domains.csv')), rules=args.rules,
                    homology=homology_path, evidence=None, manifest=args.manifest,crosswalk=crosswalk,
                    classification_mode=args.classification_mode, out=str(_layout_resource(folder / 'candidates')))))
            status['status'] = ('prepared' if args.prepare_only else
                                'completed_with_homology' if args.reference_set else 'completed_domain_screen')
        except Exception as error:
            status.update(status='failed', error=str(error))
            save_json(_layout_resource(output / 'batch.json'), summary)
            raise
        save_json(_layout_resource(output / 'batch.json'), summary)
        print(f'{datetime.now(timezone.utc).isoformat()} | TARGET {target_index}/{len(settings["proteomes"])} COMPLETE {code}; state {status["status"]}',flush=True)

    save_json(_layout_resource(output / 'batch.json'), summary)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--version',action='version',version='adhesome-detection '+VERSION)
    sub = parser.add_subparsers(dest='stage', required=True)
    p = sub.add_parser('batch', help='Detect candidates in configured target species; external context belongs to resume')
    p.add_argument('--input-dir', default='proteome_fasta_files')
    p.add_argument('--config', default=str(_layout_resource(CONFIG / 'proteomes.json')))
    p.add_argument('--out', default='results/species')
    p.add_argument('--library')
    p.add_argument('--prepare-only', action='store_true')
    p.add_argument('--check-inputs',action='store_true',help='Validate all configured FASTA headers and gene mappings, then exit without creating outputs')
    p.add_argument('--species',action='append',help='Run only selected target codes while retaining unselected species and batch status')
    p.add_argument('--resume',action='store_true',help='Continue this batch in place; verify/reuse completed targets and stages without archiving or backups')
    p.add_argument('--exclude-target',action='append',help='Omit an unfinished target; keep its files untouched and record its exclusion')
    p.add_argument('--reference-set', help='Paired reference directory; also run BLASTP and integrate its evidence')
    p.add_argument('--classification-mode', choices=['comparative','strict'], default='comparative')
    p.add_argument('--manifest', default=str(_layout_resource(CONFIG / 'adhesome_pfam_manifest.csv')))
    p.add_argument('--crosswalk', default=str(_layout_resource(CONFIG / 'adhesome_component_crosswalk.csv')))
    p.add_argument('--rules', default=str(_layout_resource(CONFIG / 'family_rules.json')))
    p.add_argument('--cpu', type=int, default=_default_thread_count(4))
    p.add_argument('--max-i', type=float, default=1e-3)
    p.add_argument('--max-c', type=float, default=1e-3)
    p.add_argument('--min-coverage', type=float, default=0.35)
    p = sub.add_parser('prepare')
    p.add_argument('--config', default=str(_layout_resource(CONFIG / 'proteomes.json')))
    p.add_argument('--out', default='results/prepared')
    p = sub.add_parser('scan')
    p.add_argument('--proteome', required=True)
    p.add_argument('--library', required=True)
    p.add_argument('--manifest', default=str(_layout_resource(CONFIG / 'adhesome_pfam_manifest.csv')))
    p.add_argument('--out', required=True)
    p.add_argument('--cpu', type=int, default=_default_thread_count(4))
    p.add_argument('--max-i', type=float, default=1e-3)
    p.add_argument('--max-c', type=float, default=1e-3)
    p.add_argument('--min-coverage', type=float, default=0.35)
    p = sub.add_parser('homology')
    for name in ('proteome','out'):
        p.add_argument('--' + name, required=True)
    for name, filename in [('references','search_references.faa'),('catalogue','search_catalogue.csv'),
                            ('metazoan','metazoan_references.faa'),('metazoan-catalogue','metazoan_catalogue.csv')]:
        p.add_argument('--'+name, default=str(REFERENCE_SET/filename))
    p.add_argument('--cpu', type=int, default=_default_thread_count(4))
    p.add_argument('--evalue', type=float, default=1e-5)
    p.add_argument('--coverage', type=float, default=0.3)
    p.add_argument('--tie-fraction', type=float, default=0.01)
    p = sub.add_parser('profiles')
    for name in ('proteome','profiles','profile-map','out'):
        p.add_argument('--' + name, required=True)
    p.add_argument('--cpu', type=int, default=_default_thread_count(4))
    p.add_argument('--max-i', type=float, default=1e-3)
    p.add_argument('--max-c', type=float, default=1e-3)
    p.add_argument('--min-coverage', type=float, default=0.35)
    p = sub.add_parser('classify')
    p.add_argument('--classification-mode', choices=['comparative','strict'], default='comparative')
    p.add_argument('--manifest', default=str(_layout_resource(CONFIG / 'adhesome_pfam_manifest.csv')))
    p.add_argument('--crosswalk', default=str(_layout_resource(CONFIG / 'adhesome_component_crosswalk.csv')))
    p.add_argument('--proteome', required=True)
    p.add_argument('--domains', required=True)
    p.add_argument('--rules', default=str(_layout_resource(CONFIG / 'family_rules.json')))
    p.add_argument('--homology')
    p.add_argument('--evidence', action='append', help='Repeat for detection profile or motif evidence; return external context through resume')
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
