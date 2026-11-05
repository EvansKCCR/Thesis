#!/usr/bin/env python3
"""Resumable integration of detection, external context, motifs, orthology and audit.

Python >=3.11. Commands write new directories and never overwrite completed runs.
See docs/WORKFLOW.md for the input contracts and examples.
"""

# Resolve source resources independently of the caller working directory.
import sys as _layout_sys
from pathlib import Path as _LayoutPath
_LAYOUT_BUNDLE_ROOT = next(p for p in _LayoutPath(__file__).resolve().parents
                          if (p/'config/layout_paths.json').is_file() and (p/'pipeline_paths.py').is_file())
if str(_LAYOUT_BUNDLE_ROOT) not in _layout_sys.path:
    _layout_sys.path.insert(0, str(_LAYOUT_BUNDLE_ROOT))
from pipeline_paths import activate as _layout_activate, resource_path as _layout_resource
_layout_activate()


import argparse
import csv
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from identity_probability import resource_hashes as probability_resources
from evaluate_validation import after_audit as evaluate_controls,resource_hashes as validation_resources
from orthogroup_context import validate_weight,configuration as context_configuration
from stage2_evidence import default_profile
from input_dependencies import InputRegistry, raw_input

ROOT = _LAYOUT_BUNDLE_ROOT
VERSION = '1.0'
DEFAULT_FAMILY_COVERAGE = 'expanded'
DETECT = _layout_resource(ROOT/'detection')
MOTIF = _layout_resource(ROOT/'motif_annotation')
ORTH = _layout_resource(ROOT/'orthology_phylogeny/scripts')
AUDIT = _layout_resource(ROOT/'annotation/integrated_adhesome_pipeline.py')
AUDIT_POLICY = _layout_resource(AUDIT.with_name('audit_policy.json'))
AUDIT_EVIDENCE = _layout_resource(ROOT/'scripts/audit_join.py')
IDENTITY_CODE = _layout_resource(ROOT/'scripts/family_identity.py')
IDENTITY_LIBRARY = _layout_resource(AUDIT.with_name('family_identity_rules.json'))
AUDIT_REPORT = _layout_resource(ROOT/'scripts/audit_report.py')
FIBRONECTIN_FAMILIES = {'fibronectin_like', 'fibronectin_III_ligand', 'fibronectin_like_candidate'}
IDENTITY_FIELDS = ('reference_subfamily','subfamily_support_method','reciprocal_reference_subfamily',
                  'phylogeny_concordant','subfamily_profile_calibrated','top_profile_subfamily',
                  'profile_margin_pass','reference_accession','reference_dataset_version','profile_domain_gate_verified')


def digest(path):
    with _layout_resource(Path(path)).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def dump(path, value):
    _layout_resource(Path(path)).write_text(json.dumps(value, indent=2)+'\n', encoding='utf-8')


def audit_configuration(coverage=DEFAULT_FAMILY_COVERAGE):
    """Select family and identity rules as one pair; never alter recruitment."""
    if coverage not in ('original','expanded'):
        raise ValueError('Unknown audit family coverage: '+str(coverage))
    if coverage=='expanded':
        from family_catalogue import pair
        rules,identities=pair()
    else:
        rules=_layout_resource(DETECT/'config/family_rules.json');identities=IDENTITY_LIBRARY
    from evaluate_family_expansion import validate_pair
    validate_pair(json.loads(_layout_resource(rules).read_text()),json.loads(_layout_resource(identities).read_text()))
    return rules,identities


def read(path, delimiter='\t'):
    with _layout_resource(Path(path)).open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f, delimiter=delimiter))


def write(path, rows, fields=None):
    fields = fields or list(dict.fromkeys(k for row in rows for k in row))
    with _layout_resource(Path(path)).open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)


def fasta(path):
    result, ident, chunks = {}, None, []
    def save():
        if ident is not None:
            seq = ''.join(chunks).upper()
            if ident in result or not seq or set(seq)-set('ACDEFGHIKLMNPQRSTVWYBXZJUO'):
                raise ValueError(f'Duplicate ID, empty or invalid sequence: {path}: {ident}')
            result[ident] = seq
    for line in _layout_resource(Path(path)).read_text(encoding='utf-8-sig').splitlines():
        if line.startswith('>'):
            save()
            if not line[1:].strip():
                raise ValueError(f'Empty FASTA header: {path}')
            ident, chunks = line[1:].split()[0], []
        elif line.strip():
            if ident is None:
                raise ValueError(f'Sequence before header: {path}')
            chunks.append(''.join(line.split()))
    save()
    if not result:
        raise ValueError(f'Empty FASTA: {path}')
    return result


def write_fasta(path, seqs):
    _layout_resource(Path(path)).write_text(''.join(f'>{k}\n{v}\n' for k,v in sorted(seqs.items())), encoding='ascii')


def execute(command, logdir):
    logdir.mkdir(parents=True, exist_ok=True)
    n = len(list(_layout_resource(logdir).glob('command_*.json')))+1
    record = dict(command=[str(x) for x in command], cwd=str(ROOT), status='running')
    dest = _layout_resource(logdir/f'command_{n:03}.json')
    dump(dest, record)
    with (_layout_resource(logdir/f'command_{n:03}.log')).open('w', encoding='utf-8') as log:
        try:
            proc = subprocess.run(record['command'], cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
            record.update(status='complete' if proc.returncode == 0 else 'failed', exit_code=proc.returncode)
        except OSError as exc:
            record.update(status='failed', error=str(exc))
            dump(dest, record)
            raise
    dump(dest, record)
    if proc.returncode:
        raise ValueError(f'Command failed; inspect {logdir}/command_{n:03}.log')


def seal(folder, kind, **metadata):
    files = {str(p.relative_to(folder)): digest(p) for p in _layout_resource(folder).rglob('*') if _layout_resource(p).is_file() and p.name != 'workflow_manifest.json'}
    dump(_layout_resource(folder/'workflow_manifest.json'), dict(workflow_version=VERSION, kind=kind, files=files, **metadata))


def verify(folder, kind, required_files=None):
    data = json.loads((_layout_resource(folder/'workflow_manifest.json')).read_text(encoding='utf-8'))
    if data['kind'] != kind:
        raise ValueError(f'Expected {kind}: {folder}')
    files={name.replace('\\','/'):sha for name,sha in data['files'].items()}
    if required_files is not None:
        absent=set(required_files)-set(files)
        if absent:
            raise ValueError(f'Missing recorded audit inputs: {sorted(absent)}')
        files={name:files[name] for name in required_files}
    for name, value in files.items():
        if digest(folder/name) != value:
            raise ValueError(f'Changed workflow input: {folder/name}')
    return data


def start(args):
    if not args.detection and not args.input_dir:
        raise ValueError('Supply --input-dir with target proteomes, or --detection with existing detection output')
    config_path = getattr(args, 'targets_config', None)
    if not args.detection and config_path is None:config_path=raw_input('target_configuration')
    if config_path:InputRegistry().require('target_configuration','target.configuration',config_path)
    excluded=set(getattr(args,'exclude_target',None) or [])
    resuming=bool(getattr(args,'resume',False))
    expected = None
    if config_path:
        config = json.loads(_layout_resource(config_path).read_text(encoding='utf-8'))
        entries = [entry for entry in config['proteomes'] if entry['code'] not in excluded]
        config=dict(config,proteomes=entries)
        codes = [r['code'] for r in entries]
        if not codes or len(codes)!=len(set(codes)) or any(not re.fullmatch(r'[A-Za-z][A-Za-z0-9_]*', s) or '__' in s for s in codes):
            raise ValueError('Target config requires unique safe species codes without double underscores')
        for entry in entries:
            if not entry.get('species') or not entry.get('path') or not entry.get('gene_regex') or not entry.get('assembly'):
                raise ValueError('Each target requires code, species, assembly, path and gene_regex')
            if re.compile(entry['gene_regex']).groups < 1:
                raise ValueError('gene_regex must capture the explicit gene ID in its first group')
        expected = set(codes)
    if not args.detection and not config_path:
        raise ValueError('New detection requires --targets-config with the target species and gene-ID rules')
    out = _layout_resource(args.out).resolve()
    if not args.detection:
        InputRegistry().require('target_raw','detection',args.input_dir)
        sys.path.insert(0,str(DETECT))
        try:
            from adhesome_pipeline import check_inputs
            _,ignored,summary=check_inputs(args.input_dir,config,_layout_resource(out/'detection'))
        finally:
            sys.path.pop(0)
        if getattr(args,'check_inputs',False):
            print(json.dumps(dict(input_root=str(_layout_resource(args.input_dir).resolve()),species=summary,ignored_fastas=ignored),indent=2))
            return
    elif getattr(args,'check_inputs',False):
        raise ValueError('--check-inputs requires raw --input-dir and --targets-config')
    if resuming and (_layout_resource(out/'workflow_manifest.json')).exists():
        raise ValueError('This handoff is already complete; use its existing outputs')
    if resuming and not (_layout_resource(out/'detection/batch.json')).is_file() and not args.detection:
        raise ValueError('No failed detection batch found in '+str(_layout_resource(out/'detection')))
    out.mkdir(parents=True, exist_ok=resuming)
    detection = _layout_resource(args.detection).resolve() if args.detection else _layout_resource(out/'detection')
    if not args.detection:
        cmd = [sys.executable, _layout_resource(DETECT/'scripts/run_detection.py'), '--out', detection, '--cpu', args.threads]
        if args.input_dir:
            cmd += ['--input-dir', _layout_resource(args.input_dir).resolve()]
        cmd += ['--config', _layout_resource(config_path).resolve()]
        if resuming:cmd+=['--resume']
        for code in sorted(excluded):cmd+=['--exclude-target',code]
        execute(cmd, _layout_resource(out/'logs'))
    batch_path=_layout_resource(detection/'batch.json')
    if _layout_resource(batch_path).exists():
        statuses=json.loads(_layout_resource(batch_path).read_text(encoding='utf-8')).get('species',{})
        unfinished=[sp for sp in ((expected or set(statuses))-excluded) if statuses.get(sp,{}).get('status') not in ('completed_with_homology','completed_domain_screen')]
        if unfinished:raise ValueError('Detection is incomplete for '+', '.join(sorted(unfinished))+'. Finish those targets before exporting the handoff.')
    hypotheses, sequences, families, sources, species = [], {}, defaultdict(dict), {}, set()
    for table in sorted(_layout_resource(detection).glob('*/candidates/candidates.csv')):
        sp = table.parents[1].name
        if sp in excluded or (expected is not None and sp not in expected):continue
        if not re.fullmatch(r'[A-Za-z][A-Za-z0-9_]*', sp) or '__' in sp:
            raise ValueError(f'Unsafe target species code: {sp}')
        species.add(sp)
        source = _layout_resource(table.parents[1]/'prepared'/f'{sp}.faa')
        seqs = fasta(source)
        sources[str(table)] = digest(table)
        sources[str(source)] = digest(source)
        for row in read(table, ','):
            fibronectin = row['family'] in FIBRONECTIN_FAMILIES
            if row['status'] not in ('supported_candidate', 'provisional', 'conflicting') and not (fibronectin and row['status']=='unresolved'):
                continue
            pid, family = row['protein_id'], row['family']
            if not pid.startswith(sp+'__'):
                raise ValueError(f'Candidate ID does not match target species: {pid}')
            if not re.fullmatch(r'[A-Za-z0-9_-]+', family):
                raise ValueError(f'Unsafe family: {family}')
            seq = seqs[pid]
            if pid in sequences and sequences[pid] != seq:
                raise ValueError(f'Conflicting protein ID: {pid}')
            sequences[pid] = seq
            motif_family = 'fibronectin_like_candidate' if fibronectin else family
            families[motif_family][pid] = seq
            hypotheses.append(dict(row, Specie=sp, sequence_sha256=hashlib.sha256(seq.encode()).hexdigest(),
                                   motif_family=motif_family,
                                   export_basis='fibronectin_screening_followup_not_family_assignment' if fibronectin and row['status']=='unresolved' else 'detection_evidence_tier',
                                   hypothesis_tier='supported' if row['status']=='supported_candidate' else 'provisional'))
    if not species or not hypotheses or (expected is not None and species != expected):
        raise ValueError('Expected configured target species and nonempty supported/provisional hypotheses')
    # Preserve all initial hypotheses for each selected protein, including generic
    # or unresolved alternatives. Selection never adds a new protein at stage 2.
    seen_hypotheses={(r['Specie'],r['protein_id'],r['family']) for r in hypotheses}
    for table in sorted(_layout_resource(detection).glob('*/candidates/candidates.csv')):
        sp=table.parents[1].name
        if sp in excluded or (expected is not None and sp not in expected):continue
        for source_row in read(table,','):
            pid,family=source_row['protein_id'],source_row['family']
            if pid not in sequences or (sp,pid,family) in seen_hypotheses:continue
            if not pid.startswith(sp+'__') or not re.fullmatch(r'[A-Za-z0-9_-]+',family):
                raise ValueError(f'Invalid retained stage 1 hypothesis: {sp}, {pid}, {family}')
            seen_hypotheses.add((sp,pid,family))
            motif_family='fibronectin_like_candidate' if family in FIBRONECTIN_FAMILIES else family
            hypotheses.append(dict(source_row,Specie=sp,sequence_sha256=hashlib.sha256(sequences[pid].encode()).hexdigest(),
                motif_family=motif_family,export_basis='retained_stage1_alternative_for_selected_protein',hypothesis_tier='unresolved'))
            families[motif_family][pid]=sequences[pid]
    (_layout_resource(out/'families')).mkdir(exist_ok=resuming)
    for family, seqs in families.items():
        write_fasta(_layout_resource(out/'families'/f'{family}.faa'), seqs)
    write(_layout_resource(out/'hypotheses.tsv'), hypotheses)
    from classification_contract import stage_one_workbook
    stage_one_workbook(hypotheses,_layout_resource(out/'Stage_1_output.xlsx'))
    coverage_rows=[]
    for table in sorted(_layout_resource(detection).glob('*/candidates/reference_coverage.csv')):
        if table.parents[1].name not in species:continue
        coverage_rows.extend(dict(row, Specie=table.parents[1].name) for row in read(table, ','))
        sources[str(table)]=digest(table)
    if coverage_rows:
        write(_layout_resource(out/'reference_coverage.tsv'), coverage_rows)
    component_rows=[]
    component_fields=[]
    for table in sorted(_layout_resource(detection).glob('*/candidates/component_hypotheses.csv')):
        if table.parents[1].name not in species:continue
        rows=read(table,',')
        if rows:
            component_fields=list(rows[0])
            component_rows.extend(r for r in rows if r['protein_id'] in sequences)
        sources[str(table)]=digest(table)
    if component_fields:
        write(_layout_resource(out/'component_hypotheses.tsv'),component_rows,component_fields)
    fibronectin_rows = [r for r in hypotheses if r['family'] in FIBRONECTIN_FAMILIES]
    write(_layout_resource(out/'fibronectin_like_hypotheses.tsv'), fibronectin_rows, list(hypotheses[0]))
    write_fasta(_layout_resource(out/'fibronectin_like_candidates.faa'), {r['protein_id']:sequences[r['protein_id']] for r in fibronectin_rows})
    for tier in ('supported', 'provisional'):
        write(_layout_resource(out/f'{tier}_hypotheses.tsv'), [r for r in hypotheses if r['hypothesis_tier']==tier], list(hypotheses[0]))
    write_fasta(_layout_resource(out/'external_candidates.faa'), sequences)
    write(_layout_resource(out/'external_evidence_template.tsv'), [dict(protein_id=k, sequence_sha256=hashlib.sha256(v.encode()).hexdigest(),
          DeepLoc_2_1='', DeepTMHMM='', SignalP_v6_0='', NCBI_CDD_architecture='', source='',
          **{field:'' for field in IDENTITY_FIELDS}) for k,v in sorted(sequences.items())])
    domains = sorted(path for path in _layout_resource(detection).glob('*/domains/domains.csv') if path.parents[1].name in species)
    source_fastas = sorted(path for path in _layout_resource(detection).glob('*/prepared/*.faa') if path.parents[1].name in species)
    execute([sys.executable, _layout_resource(MOTIF/'import_domain_regions.py'), '--input-dir', _layout_resource(out/'families'),
             '--domains', *domains, '--source-fasta', *source_fastas, '--out', _layout_resource(out/'verified_regions.tsv')], _layout_resource(out/'logs'))
    provenance = json.loads((_layout_resource(out/'verified_regions.provenance.json')).read_text(encoding='utf-8'))
    if provenance['unverified_records']:
        raise ValueError('Unverified domain source sequences; inspect provenance')
    sources.update(provenance['inputs'])
    for sp in species:
        homology=_layout_resource(detection/sp/'homology')
        for path in [_layout_resource(homology/'homology.csv'),_layout_resource(homology/'homology.json'),*(_layout_resource(homology/'reference_snapshot')).glob('*')]:
            if _layout_resource(path).is_file():sources[str(_layout_resource(path).resolve())]=digest(path)
    for name in ('adhesome_pfam_manifest.csv','adhesome_component_crosswalk.csv','family_rules.json'):
        resource=_layout_resource(DETECT/'config')/name
        sources[str(resource)]=digest(resource)
    if config_path:
        sources[str(_layout_resource(config_path).resolve())] = digest(config_path)
    seal(out, 'handoff', detection=str(detection), sources=sources,
         query_unit='control_protein' if os.environ.get('SCHISTO_ISOLATED_BRANCH') in ('calibration','validation') else 'gene',
         completed_stages=['detection'],deferred_stages=['topology','localization','phylogeny'],
         excluded_target_species=sorted(set(json.loads(_layout_resource(batch_path).read_text()).get('excluded_targets',{}))|excluded) if _layout_resource(batch_path).exists() else sorted(excluded),
         target_species=sorted(species), state='awaiting_external_context', hypotheses=len(hypotheses), proteins=len(sequences))
    print(f'Handoff ready: {out}/external_candidates.faa ({len(sequences)} proteins; {len(hypotheses)} hypotheses)')


def context_rows(path, sequences):
    rows = read(path)
    seen = set()
    required = {'protein_id','sequence_sha256','DeepLoc_2_1','DeepTMHMM','SignalP_v6_0','NCBI_CDD_architecture','source'}
    if not rows or not required <= rows[0].keys():
        raise ValueError('External evidence must use the generated template columns')
    for row in rows:
        pid = row['protein_id']
        if pid in seen or pid not in sequences:
            raise ValueError(f'Duplicate/unknown external protein: {pid}')
        seen.add(pid)
        if row['sequence_sha256'] != hashlib.sha256(sequences[pid].encode()).hexdigest():
            raise ValueError(f'External evidence sequence hash mismatch: {pid}')
        if any(v for k,v in row.items() if k not in {'protein_id','sequence_sha256','source'}) and not row['source'].strip():
            raise ValueError(f'External evidence requires a source: {pid}')
    return rows


def resume(args):
    handoff = _layout_resource(args.handoff).resolve()
    verify(handoff, 'handoff')
    if not any((args.topology_file, args.context_dir, args.external_evidence, args.allow_missing_context)):
        raise ValueError('Supply external results or explicitly use --allow-missing-context')
    out = _layout_resource(args.out).resolve()
    seqs = fasta(_layout_resource(handoff/'external_candidates.faa'))
    external = context_rows(args.external_evidence, seqs) if args.external_evidence else []
    cmd = [sys.executable, _layout_resource(MOTIF/'family_motifs.py'), '--input-dir', _layout_resource(handoff/'families'),
           '--regions', _layout_resource(handoff/'verified_regions.tsv'), '--out', _layout_resource(out/'motifs'),
           '--motif-scope',getattr(args,'motif_scope','adhesion')]
    if getattr(args,'motif_library',None):cmd+=['--library',_layout_resource(args.motif_library).resolve()]
    inputs, copies = {}, []
    for flag in ('topology_file','context_dir','evidence'):
        value = getattr(args, flag)
        if value:
            paths = _layout_resource(value).rglob('*') if _layout_resource(value).is_dir() else [value]
            inputs.update({str(_layout_resource(p).resolve()):digest(p) for p in paths if _layout_resource(p).is_file()})
            if flag == 'context_dir':
                staged = _layout_resource(out/'context')
                aliases = {
                    'DeepLoc_adhesome_candidate_list.csv':'adhesome_DeepLoc_results.csv',
                    'SignalP_adhesome_candidate_list.txt':'adhesome_SiganP_prediction_results.txt',
                    'NCBI_CDD_adhesome_hitdata.txt':'adhesome_NCBI_CDD_hitdata.txt',
                    'deeptmhmm_predicted_topologies.3line':'adhesome_predicted_topologies.3line'}
                for canonical, alternative in aliases.items():
                    source = value/canonical if (_layout_resource(value/canonical)).exists() else value/alternative
                    if canonical.endswith('.3line') and args.topology_file:
                        continue
                    if not _layout_resource(source).is_file():
                        raise ValueError(f'Missing context input: {source}')
                    copies.append((source, staged/canonical))
                value = staged
            cmd += ['--'+flag.replace('_','-'), _layout_resource(value).resolve()]
    if args.external_evidence:
        inputs[str(_layout_resource(args.external_evidence).resolve())] = digest(args.external_evidence)
    if getattr(args,'motif_library',None):
        inputs[str(_layout_resource(args.motif_library).resolve())]=digest(args.motif_library)
    contract=dict(handoff_manifest=digest(_layout_resource(handoff/'workflow_manifest.json')),sources=inputs,
                  command=[str(x) for x in cmd[1:]],allow_missing_context=args.allow_missing_context)
    retrying=bool(getattr(args,'retry_failed',False))
    checkpoint=_layout_resource(out/'annotation_recovery.json')
    if _layout_resource(out).exists():
        if not retrying:
            raise ValueError('Annotation output already exists; use --retry-failed for a failed annotation before motif output, or recover for a failed downstream run')
        if (_layout_resource(out/'workflow_manifest.json')).exists():
            raise ValueError('Annotation is already complete; use recover --run with this annotation directory for downstream failures')
        if (_layout_resource(out/'motifs')).exists() or (_layout_resource(out/'downstream')).exists():
            raise ValueError('Retry requires an annotation failure before motif output was created; existing partial motif output must be reviewed separately')
        if _layout_resource(checkpoint).exists():
            previous=json.loads(_layout_resource(checkpoint).read_text(encoding='utf-8'))
            if previous.get('status')!='failed' or previous.get('contract')!=contract:
                raise ValueError('Failed annotation inputs/settings changed or attempt is not marked failed; retry refused')
        else:
            # Older releases did not checkpoint annotation. Adopt only the narrow
            # pre-output failure with an identical logged command and staged inputs.
            logs=sorted((_layout_resource(out/'logs')).glob('command_*.json'))
            previous=json.loads(_layout_resource(logs[-1]).read_text(encoding='utf-8')) if logs else {}
            if (previous.get('status')!='failed' or previous.get('command',[])[1:]!=contract['command']
                    or args.topology_file or args.evidence or args.external_evidence or not copies
                    or any(p.name not in ('context','logs') for p in out.iterdir())):
                raise ValueError('Cannot verify this older failed annotation for in-place retry')
            if any(not _layout_resource(dest).is_file() or digest(source)!=digest(dest) for source,dest in copies):
                raise ValueError('Returned context differs from the failed attempt; retry refused')
    elif retrying:
        raise ValueError('--retry-failed requires the existing failed annotation directory')
    out.mkdir(parents=True, exist_ok=retrying)
    state=dict(status='running',contract=contract)
    dump(checkpoint,state)
    try:
        for source,dest in copies:
            dest.parent.mkdir(parents=True,exist_ok=True)
            if _layout_resource(source).resolve()!=_layout_resource(dest).resolve():shutil.copyfile(source,dest)
        if args.external_evidence:write(_layout_resource(out/'external_evidence.tsv'),external)
        execute(cmd, _layout_resource(out/'logs'))
    except (ValueError,OSError,KeyError) as exc:
        state.update(status='failed',error=str(exc))
        dump(checkpoint,state)
        raise
    state['status']='complete'
    dump(checkpoint,state)
    seal(out, 'annotation', handoff=str(handoff), handoff_manifest=digest(_layout_resource(handoff/'workflow_manifest.json')), sources=inputs,
         context_policy='missing_allowed' if args.allow_missing_context else 'external_supplied')
    print(f'Annotation complete: {out}.')
    if not getattr(args, 'annotation_only', False):
        downstream(argparse.Namespace(annotation=out,out=_layout_resource(out/'downstream'),
            reference_raw=args.reference_raw,threads=args.threads,analysis_threads=args.analysis_threads,
            compare=args.compare,phylogeny=getattr(args,'phylogeny',False),
            benchmark_refogs=getattr(args,'benchmark_refogs',None),benchmark_input_dir=getattr(args,'benchmark_input_dir',None),
            benchmark_protocol=getattr(args,'benchmark_protocol','open_orthobench'),
            benchmark_scope=getattr(args,'benchmark_scope','all'),benchmark_selection=getattr(args,'benchmark_selection',None),
            benchmark_level=getattr(args,'benchmark_level','root'),benchmark_quality_policy=getattr(args,'benchmark_quality_policy',None),
            orthofinder_inflation=getattr(args,'orthofinder_inflation',None),orthofinder_search=getattr(args,'orthofinder_search',None),
            orthofinder_method=getattr(args,'orthofinder_method',None),
            iqtree_threads=args.iqtree_threads,allow_partial_phylogeny=args.allow_partial_phylogeny,
            phylogeny_grouping=args.phylogeny_grouping,
            family_coverage=getattr(args,'family_coverage',DEFAULT_FAMILY_COVERAGE),
            orthogroup_context_weight=getattr(args,'orthogroup_context_weight',None),
            phylogenetic_identity_resolution=getattr(args,'phylogenetic_identity_resolution',False),
            phylogenetic_identity_evidence=getattr(args,'phylogenetic_identity_evidence',None),
            stage2_profile=getattr(args,'stage2_profile',default_profile())))


def downstream(args):
    annotation=_layout_resource(args.annotation).resolve()
    meta=verify(annotation,'annotation')
    handoff=_layout_resource(Path(meta['handoff']))
    verify(handoff,'handoff')
    if digest(_layout_resource(handoff/'workflow_manifest.json'))!=meta['handoff_manifest']:
        raise ValueError('Annotation handoff changed')
    out=_layout_resource(args.out).resolve()
    recovering=getattr(args,'recovering',False)
    if recovering and (_layout_resource(out/'workflow_manifest.json')).exists():
        verify(out,'completed_workflow')
        prior=json.loads((_layout_resource(out/'recovery.json')).read_text())
        retry_trees=bool(getattr(args,'phylogeny',False) and prior['stages'].get('phylogeny_run',{}).get('status')!='complete')
        retry_benchmark=bool(getattr(args,'retry_benchmark',False) or prior['stages'].get('benchmark',{}).get('status')=='failed')
        if retry_trees or retry_benchmark:
            archive_attempt(_layout_resource(out/'workflow_manifest.json'))
        else:
            print(f'Run already completed: {out}')
            return
    out.mkdir(parents=True,exist_ok=recovering)
    checkpoint=_layout_resource(out/'recovery.json')
    if recovering:
        state=json.loads(_layout_resource(checkpoint).read_text())
        if state['annotation_sha256']!=digest(_layout_resource(annotation/'workflow_manifest.json')):
            raise ValueError('Annotation changed since failed run; create a new run')
        for key in ('threads','analysis_threads','iqtree_threads','phylogeny','phylogeny_grouping','phylogenetic_identity_resolution','phylogenetic_identity_evidence','benchmark_refogs','benchmark_input_dir','benchmark_protocol','benchmark_scope','benchmark_selection','benchmark_level','benchmark_quality_policy'):
            if hasattr(args,key):
                value=getattr(args,key)
                state['settings'][key]=str(value) if isinstance(value,Path) else value
        if getattr(args,'retry_benchmark',False):state['stages'].pop('benchmark',None)
        audit_meta=_layout_resource(out/'final/workflow_manifest.json')
        if state['stages'].get('audit',{}).get('status')=='complete' and _layout_resource(audit_meta).exists():
            saved=json.loads(_layout_resource(audit_meta).read_text())
            selected_rules,selected_identity=audit_configuration(state['settings'].get('family_coverage',DEFAULT_FAMILY_COVERAGE))
            if (saved.get('rules_sha256')!=digest(selected_rules) or saved.get('audit_script_sha256')!=digest(AUDIT)
                    or saved.get('audit_policy_sha256')!=digest(AUDIT_POLICY) or saved.get('audit_evidence_sha256')!=digest(AUDIT_EVIDENCE)
                    or saved.get('classification_contract_sha256')!=digest(_layout_resource(ROOT/'scripts/classification_contract.py'))
                    or saved.get('output_schema_sha256')!=digest(_layout_resource(ROOT/'templates/output_schema.json'))
                    or saved.get('identity_code_sha256')!=digest(IDENTITY_CODE)
                    or saved.get('stage2_profile','original')!=state['settings'].get('stage2_profile','original')
                    or saved.get('phylogenetic_identity_resolution',False)!=bool(getattr(args,'phylogenetic_identity_resolution',False))
                    or saved.get('phylogenetic_identity_evidence_sha256')!=(digest(args.phylogenetic_identity_evidence) if getattr(args,'phylogenetic_identity_resolution',False) and getattr(args,'phylogenetic_identity_evidence',None) and _layout_resource(args.phylogenetic_identity_evidence).is_file() else None)
                    or saved.get('stage2_evidence_code_sha256')!=digest(_layout_resource(ROOT/'scripts/stage2_evidence.py'))
                    or saved.get('stage2_evidence_rules_sha256')!=digest(_layout_resource(ROOT/'annotation/stage2_evidence_rules.json'))
                    or saved.get('identity_library_sha256')!=digest(selected_identity)
                    or saved.get('audit_report_sha256')!=digest(AUDIT_REPORT)
                    or saved.get('orthogroup_context_weight')!=getattr(args,'orthogroup_context_weight',None)
                    or saved.get('orthogroup_context_code_sha256')!=digest(_layout_resource(ROOT/'scripts/orthogroup_context.py'))
                    or saved.get('identity_probability_resources')!=probability_resources()
                    or ('validation_evaluation_resources' in saved and saved['validation_evaluation_resources']!=validation_resources())
                    or saved.get('motif_snapshot_sha256')!=digest(_layout_resource(MOTIF/'Motif_library/adhesion_core.json'))):
                state['stages'].pop('audit',None)
                if state['settings'].get('phylogeny_grouping')=='family':
                    state['stages'].pop('phylogeny_prepare',None);state['stages'].pop('phylogeny_run',None)
                dump(checkpoint,state)
                print('Current audit rules/code changed: refresh audit and dependent family-tree preparation; reuse orthology.',flush=True)
    else:
        state=dict(settings={k:str(v) if isinstance(v,Path) else v for k,v in vars(args).items() if k!='func'},
                   annotation_sha256=digest(_layout_resource(annotation/'workflow_manifest.json')),stages={})
        state['settings'].setdefault('stage2_profile',default_profile())
        for key in ('annotation','out','reference_raw','compare','benchmark_refogs','benchmark_input_dir','benchmark_selection','benchmark_quality_policy'):
            if state['settings'].get(key) is not None:
                state['settings'][key]=str(_layout_resource(Path(state['settings'][key])).resolve())
        dump(checkpoint,state)
    programs=[]
    if state['stages'].get('orthology',{}).get('status')!='complete' and not (_layout_resource(out/'orthology/inference_complete.json')).exists():
        programs.append('orthofinder')
    missing=[name for name in programs if not shutil.which(name)]
    if missing:
        raise ValueError('Missing executables for unfinished stages: '+', '.join(missing)+'. Install dependencies, then use recover --run '+str(out))
    def stage(name,folder,action,preserve=False):
        previous=state['stages'].get(name,{})
        if name=='benchmark' and previous.get('status')=='complete' and (_layout_resource(folder/'summary.json')).exists():
            data=json.loads((_layout_resource(folder/'summary.json')).read_text())
            configuration=dict(reference_raw=str(_layout_resource(args.reference_raw).resolve()),
                refogs=str(_layout_resource(args.benchmark_refogs).resolve()) if getattr(args,'benchmark_refogs',None) else None,
                input_dir=str(_layout_resource(args.benchmark_input_dir).resolve()) if getattr(args,'benchmark_input_dir',None) else None,
                scope=getattr(args,'benchmark_scope','all'),
                selection=str(_layout_resource(args.benchmark_selection).resolve()) if getattr(args,'benchmark_selection',None) else None,
                level=getattr(args,'benchmark_level','root'),
                quality_policy=str(_layout_resource(args.benchmark_quality_policy).resolve()) if getattr(args,'benchmark_quality_policy',None) else None)
            saved_configuration=dict(data.get('configuration',{}))
            saved_configuration.setdefault('scope','all');saved_configuration.setdefault('selection',None)
            saved_configuration.setdefault('level','root');saved_configuration.setdefault('quality_policy',None)
            changed=saved_configuration!=configuration
            if data.get('protocol')!=getattr(args,'benchmark_protocol','open_orthobench'):changed=True
            if data.get('status')=='unavailable' and _layout_resource(args.reference_raw).exists() and any(_layout_resource(args.reference_raw).glob('**/RefOGs')):changed=True
            if data.get('status')=='invalid_benchmark_inputs' and not data.get('inputs'):changed=True
            for filename,sha in data.get('inputs',{}).items():
                if not _layout_resource(Path(filename)).is_file() or digest(_layout_resource(Path(filename)))!=sha:changed=True;break
            if changed:
                previous={}
                state['stages'].pop(name,None)
        if previous.get('status')=='complete':
            for filename,sha in previous['files'].items():
                if digest(out/filename)!=sha:
                    raise ValueError(f'Completed checkpoint changed: {filename}')
            print(f'REUSE {name}',flush=True)
            return
        if _layout_resource(folder).exists() and not preserve:
            archive_attempt(folder)
        state['stages'][name]=dict(status='running')
        dump(checkpoint,state)
        try:
            action()
            files={str(p.relative_to(out)):digest(p) for p in _layout_resource(folder).rglob('*') if _layout_resource(p).is_file()}
            state['stages'][name]=dict(status='complete',files=files)
            dump(checkpoint,state)
        except Exception as exc:
            state['stages'][name]=dict(status='failed',error=str(exc))
            dump(checkpoint,state)
            raise
    reference_command=[sys.executable,_layout_resource(ORTH/'prepare_references.py'),'--raw',_layout_resource(args.reference_raw).resolve(),
                       '--handoff',handoff,'--out',_layout_resource(out/'references')]
    for flag in ('benchmark_refogs','benchmark_input_dir'):
        if getattr(args,flag,None):reference_command+=['--'+flag.replace('_','-'),getattr(args,flag)]
    for filename in getattr(args,'reference_file',None) or []:
        reference_command+=['--reference-file',filename]
    if getattr(args,'reference_panel',None):reference_command+=['--reference-panel',args.reference_panel]
    stage('references',_layout_resource(out/'references'),lambda:execute(reference_command,_layout_resource(out/'logs')))
    stage('panel',_layout_resource(out/'panel'),lambda:calibrate(argparse.Namespace(handoff=handoff,panel=_layout_resource(out/'references/panel.auto.json'),out=_layout_resource(out/'panel'))))
    stage('orthology',_layout_resource(out/'orthology'),lambda:orthology(argparse.Namespace(panel=_layout_resource(out/'panel'),out=_layout_resource(out/'orthology'),results=None,input_manifest=None,
        threads=args.threads,analysis_threads=args.analysis_threads,compare=args.compare,recovering=recovering,inside_workflow=True,
        orthofinder_inflation=getattr(args,'orthofinder_inflation',None),orthofinder_search=getattr(args,'orthofinder_search',None),
        orthofinder_method=getattr(args,'orthofinder_method',None))),preserve=True)
    optional_errors=[]
    benchmark_refreshed=state['stages'].get('benchmark',{}).get('status')!='complete'
    try:
        stage('benchmark',_layout_resource(out/'benchmark'),lambda:benchmark(argparse.Namespace(orthology=_layout_resource(out/'orthology'),out=_layout_resource(out/'benchmark'),
            reference_raw=args.reference_raw,refogs=getattr(args,'benchmark_refogs',None),input_dir=getattr(args,'benchmark_input_dir',None),
            protocol=getattr(args,'benchmark_protocol','open_orthobench'),
            scope=getattr(args,'benchmark_scope','all'),selection=getattr(args,'benchmark_selection',None),
            level=getattr(args,'benchmark_level','root'),quality_policy=getattr(args,'benchmark_quality_policy',None))))
    except Exception as exc:
        optional_errors.append('benchmark: '+str(exc))
        dump(_layout_resource(out/'benchmark_failure.json'),dict(status='failed',validated_weight=None,error=str(exc)))
    summary_path=_layout_resource(out/'benchmark/summary.json')
    if _layout_resource(summary_path).exists() and json.loads(_layout_resource(summary_path).read_text()).get('status')=='quality_rejected':
        optional_errors.append('benchmark_quality: measured orthology performance fails configured acceptance criteria')
    if state['stages'].get('audit',{}).get('status')=='complete':
        saved=json.loads((_layout_resource(out/'final/workflow_manifest.json')).read_text())
        if benchmark_refreshed or saved.get('benchmark_summary_sha256')!=(digest(summary_path) if _layout_resource(summary_path).exists() else None):
            state['stages'].pop('audit',None)
            if getattr(args,'phylogeny_grouping','orthogroup')=='family':
                state['stages'].pop('phylogeny_prepare',None);state['stages'].pop('phylogeny_run',None)
            dump(checkpoint,state)
    audit_refreshed=state['stages'].get('audit',{}).get('status')!='complete'
    stage('audit',_layout_resource(out/'final'),lambda:audit(argparse.Namespace(annotation=annotation,orthology=_layout_resource(out/'orthology'),out=_layout_resource(out/'final'),benchmark=_layout_resource(out/'benchmark'),
          family_coverage=getattr(args,'family_coverage',state['settings'].get('family_coverage',DEFAULT_FAMILY_COVERAGE)),
          orthogroup_context_weight=getattr(args,'orthogroup_context_weight',state['settings'].get('orthogroup_context_weight')),
          phylogenetic_identity_resolution=getattr(args,'phylogenetic_identity_resolution',False),
          phylogenetic_identity_evidence=getattr(args,'phylogenetic_identity_evidence',None),
          stage2_profile=getattr(args,'stage2_profile',state['settings'].get('stage2_profile','original')))))
    if audit_refreshed and getattr(args,'phylogeny_grouping','orthogroup')=='family':
        state['stages'].pop('phylogeny_prepare',None);state['stages'].pop('phylogeny_run',None)
        dump(checkpoint,state)
    print(f'Audit complete and available independently of phylogeny: {out}/final/Stage_2_output.xlsx',flush=True)
    candidates=_layout_resource(out/'orthology/candidates')
    if getattr(args,'phylogeny',False):
      try:
        dump(_layout_resource(out/'phylogeny_status.json'),dict(status='requested'))
        if read(_layout_resource(candidates/'candidate_HOGs.tsv')) or getattr(args,'phylogeny_grouping','orthogroup')=='family':
            if state['stages'].get('phylogeny_run',{}).get('status')!='complete':
                missing=[name for name in ('mafft','trimal','iqtree2') if not shutil.which(name)]
                if missing:
                    raise ValueError('Audit completed. Missing phylogeny executables: '+', '.join(missing)+'. Recover after installing them.')
            cmd=[sys.executable,_layout_resource(ORTH/'phylogeny.py'),'prepare','--hogs',_layout_resource(candidates/'candidate_HOGs.tsv'),
                 '--orthologues',_layout_resource(candidates/'direct_orthologues.tsv'),'--proteomes',_layout_resource(out/'panel/proteomes'),
                 '--species-list',_layout_resource(out/'panel/species.json'),'--output',_layout_resource(out/'phylogeny')]
            if getattr(args,'phylogeny_grouping','orthogroup')=='family':
                cmd+=['--audit-workbook',_layout_resource(out/'final/Stage_2_output.xlsx'),'--target-species',_layout_resource(out/'panel/targets.json'),
                      '--candidate-report',_layout_resource(candidates/'candidate_report.tsv')]
            dup=_layout_resource(candidates/'candidate_duplications.tsv')
            cmd+=['--duplications',dup] if _layout_resource(dup).exists() else ['--allow-missing-duplications']
            stage('phylogeny_prepare',_layout_resource(out/'phylogeny'),lambda:execute(cmd,_layout_resource(out/'logs')))
            tree_cmd=[sys.executable,_layout_resource(ORTH/'phylogeny.py'),'run','--project',_layout_resource(out/'phylogeny'),'--threads',args.threads,
                      '--iqtree-threads',getattr(args,'iqtree_threads',2)]
            if getattr(args,'allow_partial_phylogeny',False):
                tree_cmd+=['--allow-partial']
            if recovering:
                tree_cmd+=['--retry-failed']
            stage('phylogeny_run',_layout_resource(out/'phylogeny'),lambda:execute(tree_cmd,_layout_resource(out/'logs')),preserve=True)
            tree_summary=_layout_resource(out/'phylogeny/run_summary.json')
            dump(_layout_resource(out/'phylogeny_status.json'),dict(status=json.loads(_layout_resource(tree_summary).read_text()).get('status','complete') if _layout_resource(tree_summary).exists() else 'complete'))
        else:
            dump(_layout_resource(out/'phylogeny_status.json'),dict(status='no_candidate_groups',trees_inferred=0))
      except Exception as exc:
        optional_errors.append('phylogeny: '+str(exc))
        dump(_layout_resource(out/'phylogeny_status.json'),dict(status='failed',error=str(exc),audit_completed=True))
        print('Optional phylogeny failed; completed audit is retained: '+str(exc),flush=True)
    else:
        dump(_layout_resource(out/'phylogeny_status.json'),dict(status='not_requested'))
    summary=_layout_resource(out/'phylogeny/run_summary.json')
    status_file=_layout_resource(out/'phylogeny_status.json')
    tree_status=json.loads(_layout_resource(summary).read_text())['status'] if _layout_resource(summary).exists() else json.loads(_layout_resource(status_file).read_text()).get('status','not_run') if _layout_resource(status_file).exists() else 'not_run'
    if any(error.startswith('phylogeny:') for error in optional_errors):tree_status='failed'
    seal(out,'completed_workflow',annotation=str(annotation),phylogeny_status=tree_status,
         optional_errors=optional_errors,benchmark_status='failed' if any(error.startswith('benchmark:') for error in optional_errors) else json.loads(_layout_resource(summary_path).read_text()).get('status') if _layout_resource(summary_path).exists() else 'failed',
         completion_status='completed_with_optional_failures' if optional_errors or tree_status=='partial_failure' else 'complete')
    print(f'Workflow audit ready (phylogeny: {tree_status}): {out}/final/Stage_2_output.xlsx')


def archive_attempt(folder):
    """Preserve incomplete outputs beside their stage; never delete previous evidence."""
    if os.environ.get('SCHISTO_NO_ARCHIVES')=='1':
        raise ValueError('Partial stage retained in place (no-backup policy): '+str(folder)+
                         '. Resume the unfinished tool or use a separate recomputation run; no archive was created.')
    number=1
    while _layout_resource(folder.with_name(folder.name+f'_attempt_{number:03d}')).exists():
        number+=1
    folder.rename(_layout_resource(folder.with_name(folder.name+f'_attempt_{number:03d}')))


def recover(args):
    run=_layout_resource(args.run).resolve()
    if not (_layout_resource(run/'recovery.json')).exists() and (_layout_resource(run/'downstream/recovery.json')).exists():
        run=_layout_resource(run/'downstream')
    if not (_layout_resource(run/'recovery.json')).exists() and getattr(args,'adopt_existing',False):
        if (_layout_resource(run/'downstream/panel')).exists():
            run=_layout_resource(run/'downstream')
        annotation=(_layout_resource(args.annotation or run.parent)).resolve()
        ann=verify(annotation,'annotation')
        if (_layout_resource(run/'workflow_manifest.json')).exists():
            raise ValueError('Existing run is sealed as complete; do not modify it in place')
        panel=verify(_layout_resource(run/'panel'),'panel')
        if panel['handoff_manifest']!=ann['handoff_manifest']:
            raise ValueError('Existing panel and annotation use different handoffs')
        generated=json.loads((_layout_resource(run/'references/panel.auto.json')).read_text())
        if generated!=json.loads((_layout_resource(run/'panel/panel.snapshot.json')).read_text()):
            raise ValueError('Generated reference panel differs from staged panel')
        state=dict(annotation_sha256=digest(_layout_resource(annotation/'workflow_manifest.json')),stages={},settings=dict(
            annotation=str(annotation),out=str(run),reference_raw=str(raw_input('reference_raw')),
            threads=args.threads or 8,analysis_threads=args.analysis_threads or 8,iqtree_threads=args.iqtree_threads or 2,
            compare=None,phylogeny=False,allow_partial_phylogeny=False))
        def adopt(name,files):
            state['stages'][name]=dict(status='complete',files={str(p.relative_to(run)):digest(p) for p in files if _layout_resource(p).is_file()})
        adopt('references',(_layout_resource(run/'references')).rglob('*'))
        adopt('panel',(_layout_resource(run/'panel')).rglob('*'))
        if (_layout_resource(run/'orthology/workflow_manifest.json')).exists():
            orth=verify(_layout_resource(run/'orthology'),'orthology')
            if orth['panel_manifest']!=digest(_layout_resource(run/'panel/workflow_manifest.json')):
                raise ValueError('Existing orthology used a different panel')
            adopt('orthology',(_layout_resource(run/'orthology')).rglob('*'))
        phylo=_layout_resource(run/'phylogeny/manifest.json')
        if _layout_resource(phylo).exists() and 'orthology' in state['stages'] and getattr(args,'phylogeny',False):
            preparation=json.loads(_layout_resource(phylo).read_text())
            state['settings']['phylogeny_grouping']=preparation.get('grouping','orthogroup')
            prepared_files=[phylo,_layout_resource(run/'phylogeny/groups.tsv')]
            for path,sha in preparation['inputs'].items():
                if digest(path)!=sha:
                    raise ValueError(f'Phylogeny preparation input changed: {path}')
            for group in preparation['groups']:
                for name,sha in group['files'].items():
                    path=_layout_resource(run/'phylogeny')/group['group']/name
                    if digest(path)!=sha:
                        raise ValueError(f'Prepared group changed: {path}')
                    prepared_files.append(path)
            adopt('phylogeny_prepare',prepared_files)
        dump(_layout_resource(run/'recovery.json'),state)
    if not (_layout_resource(run/'recovery.json')).exists():
        raise ValueError('No recovery checkpoint found. For a run created before checkpoints were added, use --adopt-existing and --annotation with its completed annotation directory.')
    state=json.loads((_layout_resource(run/'recovery.json')).read_text())
    settings=state['settings']
    settings.setdefault('stage2_profile','original')  # Older checkpoints retain historical behavior.
    for key in ('annotation','out','reference_raw','compare','benchmark_refogs','benchmark_input_dir','benchmark_selection','benchmark_quality_policy'):
        if settings.get(key) is not None:
            settings[key]=_layout_resource(Path(settings[key]))
    settings['out']=run
    settings['phylogeny']=settings.get('phylogeny',False)
    settings.setdefault('phylogenetic_identity_resolution',False)
    if settings.get('phylogenetic_identity_evidence'):settings['phylogenetic_identity_evidence']=_layout_resource(Path(settings['phylogenetic_identity_evidence']))
    settings['retry_benchmark']=getattr(args,'retry_benchmark',False)
    for key in ('benchmark_scope','benchmark_protocol','benchmark_selection','benchmark_level','benchmark_quality_policy'):
        if getattr(args,key,None) is not None:
            value=getattr(args,key)
            settings[key]=_layout_resource(value).resolve() if isinstance(value,Path) else value
    if getattr(args,'benchmark_scope',None)=='all' and getattr(args,'benchmark_selection',None) is None:
        settings['benchmark_selection']=None
    if getattr(args,'phylogeny',None) is not None:settings['phylogeny']=args.phylogeny
    if getattr(args,'phylogeny_grouping',None) is not None:settings['phylogeny_grouping']=args.phylogeny_grouping
    for key in ('threads','analysis_threads','iqtree_threads'):
        if getattr(args,key) is not None:
            settings[key]=getattr(args,key)
    if getattr(args,'phylogenetic_identity_resolution',None) is not None:settings['phylogenetic_identity_resolution']=args.phylogenetic_identity_resolution
    if getattr(args,'phylogenetic_identity_evidence',None) is not None:settings['phylogenetic_identity_evidence']=args.phylogenetic_identity_evidence
    settings['recovering']=True
    downstream(argparse.Namespace(**settings))


def calibrate(args):
    handoff = _layout_resource(args.handoff).resolve()
    handoff_meta = verify(handoff, 'handoff')
    targets = set(handoff_meta.get('target_species') or [r['Specie'] for r in read(_layout_resource(handoff/'hypotheses.tsv'))])
    config = json.loads(_layout_resource(args.panel).read_text(encoding='utf-8'))
    entries = config['proteomes']
    codes = [e['species'] for e in entries]
    if len(codes)!=len(set(codes)) or not targets < set(codes) or any(not re.fullmatch(r'[A-Za-z][A-Za-z0-9_]*', s) or '__' in s for s in codes):
        raise ValueError('Panel requires unique safe species codes, every handoff target, and at least one reference')
    out = _layout_resource(args.out).resolve()
    out.mkdir(parents=True, exist_ok=False)
    (_layout_resource(out/'proteomes')).mkdir()
    allseqs, sources, counts = {}, {}, []
    for entry in entries:
        sp = entry['species']
        source = (_layout_resource(args.panel.resolve().parent/entry['fasta'])).resolve()
        benchmark_unit=entry.get('unit_level')=='benchmark_protein' and sp not in targets
        control_unit=entry.get('unit_level')=='control_protein' and sp in targets and handoff_meta.get('query_unit')=='control_protein'
        if not entry.get('source') or (entry.get('gene_level') is not True and not benchmark_unit and not control_unit):
            raise ValueError(f'{sp}: declare gene_level=true and document the proteome source/release')
        seqs = fasta(source)
        if any(not re.fullmatch(re.escape(sp)+r'__[A-Za-z0-9_.-]+', k) for k in seqs):
            raise ValueError(f'{sp}: expected gene-level IDs {sp}__gene_id; prepare/collapse raw isoforms first')
        allseqs.update(seqs)
        write_fasta(_layout_resource(out/'proteomes'/f'{sp}.faa'), seqs)
        sources[str(source)] = digest(source)
        counts.append(dict(species=sp, genes=len(seqs) if not (benchmark_unit or control_unit) else '',units=len(seqs),
                           unit_level='benchmark_protein' if benchmark_unit else 'control_protein' if control_unit else 'gene',source=entry['source']))
    candidates = fasta(_layout_resource(handoff/'external_candidates.faa'))
    mapping = []
    for pid, seq in candidates.items():
        if allseqs.get(pid) != seq:
            raise ValueError(f'Candidate missing or sequence changed in whole-proteome panel: {pid}')
        mapping.append(dict(candidate_set='workflow',candidate_id=pid,prepared_id=pid,status='matched'))
    write(_layout_resource(out/'candidate_audit.tsv'), mapping)
    write(_layout_resource(out/'proteome_summary.tsv'), counts)
    dump(_layout_resource(out/'species.json'), codes)
    dump(_layout_resource(out/'targets.json'), sorted(targets))
    dump(_layout_resource(out/'panel.snapshot.json'), config)
    seal(out, 'panel', handoff=str(handoff), handoff_manifest=digest(_layout_resource(handoff/'workflow_manifest.json')), sources=sources,
         interpretation='Reference-panel sensitivity analysis, not independent biological validation')
    print(f'Panel ready: {out}. Rerun orthology for this panel; prior runs remain unchanged.')


def orthofinder_configuration(args):
    """Explicit inference controls; never retune completed assignments during recovery."""
    settings={key:getattr(args,'orthofinder_'+key,None) for key in ('inflation','search','method')}
    inflation=settings['inflation']
    if inflation is not None and (isinstance(inflation,bool) or not math.isfinite(inflation) or not 1<=inflation<=10):
        raise ValueError('OrthoFinder inflation must be finite and between 1 and 10')
    if settings['search'] not in (None,'diamond','diamond_ultra_sens','blast','mmseqs'):
        raise ValueError('Unsupported OrthoFinder search configuration')
    if settings['method'] not in (None,'msa','dendroblast'):
        raise ValueError('Unsupported OrthoFinder inference method')
    options=[]
    for name,flag in [('inflation','-I'),('search','-S'),('method','-M')]:
        if settings[name] is not None:options.extend([flag,str(settings[name])])
    return settings,options


def orthology(args):
    inference_settings,inference_options=orthofinder_configuration(args)
    if args.results and inference_options:
        raise ValueError('Inference tuning options cannot change imported results; run new inference instead')
    panel = _layout_resource(args.panel).resolve()
    panel_meta = verify(panel, 'panel')
    out = _layout_resource(args.out).resolve()
    recovering=getattr(args,'recovering',False)
    out.mkdir(parents=True, exist_ok=recovering)
    if args.results:
        results = _layout_resource(args.results).resolve()
        # Bind external results to the exact input panel, not merely matching species labels.
        if not args.input_manifest or digest(args.input_manifest) != digest(_layout_resource(panel/'workflow_manifest.json')):
            raise ValueError('External results require --input-manifest from the panel used for inference')
    else:
        receipt=_layout_resource(out/'inference_complete.json')
        if recovering and _layout_resource(receipt).exists():
            saved=json.loads(_layout_resource(receipt).read_text())
            if any(value is not None for value in inference_settings.values()) and saved.get('inference_settings')!=inference_settings:
                raise ValueError('Completed inference settings cannot be recalibrated in place; use a new orthology output')
            if saved['panel_manifest']!=digest(_layout_resource(panel/'workflow_manifest.json')):
                raise ValueError('Panel changed since inference')
            for name,sha in saved.get('files',{}).items():
                if digest(_layout_resource(out/'inference')/name)!=sha:
                    raise ValueError(f'Completed inference output changed: {name}')
        else:
            if (_layout_resource(out/'inference')).exists():
                archive_attempt(_layout_resource(out/'inference'))
            inference_command=['orthofinder','-f',_layout_resource(panel/'proteomes'),'-t',args.threads,'-a',args.analysis_threads,'-o',_layout_resource(out/'inference'),*inference_options]
            execute(inference_command, _layout_resource(out/'logs'))
            dump(receipt,dict(panel_manifest=digest(_layout_resource(panel/'workflow_manifest.json')),inference_settings=inference_settings,
                inference_command=[str(item) for item in inference_command],
                files={str(p.relative_to(_layout_resource(out/'inference'))):digest(p) for p in (_layout_resource(out/'inference')).rglob('*') if _layout_resource(p).is_file()}))
        hits = list((_layout_resource(out/'inference')).glob('Results_*'))
        if len(hits)!=1:
            raise ValueError('Cannot identify one Results_* directory; use explicit --results on a new output')
        results = hits[0]
    if recovering and (_layout_resource(out/'candidates')).exists():
        archive_attempt(_layout_resource(out/'candidates'))
    execute([sys.executable, _layout_resource(ORTH/'pipeline.py'),'extract','--results',results,
             '--candidate-set','workflow','--candidate-audit',_layout_resource(panel/'candidate_audit.tsv'),
             '--species-list',_layout_resource(panel/'species.json'),'--target-species',_layout_resource(panel/'targets.json'),'--output',_layout_resource(out/'candidates')], _layout_resource(out/'logs'))
    report = read(_layout_resource(out/'candidates/candidate_report.tsv'))
    reference_columns = [k for k in report[0] if k.endswith('_orthologues')]
    review = []
    for row in report:
        support = [k.removesuffix('_orthologues') for k in reference_columns if row[k]]
        ambiguous = any('many' in v for k,v in row.items() if k.endswith('_relationship'))
        review.append(dict(candidate_id=row['candidate_id'], reference_count=len(support),
                           references=';'.join(support), calibration_review='no_direct_reference' if not support else 'duplication_or_paralogy_review' if ambiguous else 'direct_reference_detected'))
    write(_layout_resource(out/'calibration_review.tsv'), review)
    if args.compare:
        previous = verify(_layout_resource(args.compare).resolve(), 'orthology')
        previous_panel = verify(_layout_resource(Path(previous['panel'])), 'panel')
        if previous_panel['handoff_manifest'] != panel_meta['handoff_manifest']:
            raise ValueError('Calibration comparison requires the same sequence-verified handoff')
        old = {r['candidate_id']:r for r in read(_layout_resource(args.compare/'calibration_review.tsv'))}
        if set(old) != {r['candidate_id'] for r in review}:
            raise ValueError('Calibration comparison requires the same candidate set')
        write(_layout_resource(out/'calibration_comparison.tsv'), [dict(candidate_id=r['candidate_id'],
              previous_references=old[r['candidate_id']]['references'], current_references=r['references'],
              previous_review=old[r['candidate_id']]['calibration_review'], current_review=r['calibration_review']) for r in review])
    if not getattr(args,'inside_workflow',False):
        benchmark_command(results,panel,getattr(args,'reference_raw',raw_input('reference_raw')),_layout_resource(out/'benchmark'),_layout_resource(out/'logs'),
                          getattr(args,'benchmark_refogs',None),getattr(args,'benchmark_input_dir',None),getattr(args,'benchmark_protocol','open_orthobench'),
                          getattr(args,'benchmark_scope','all'),getattr(args,'benchmark_selection',None),
                          getattr(args,'benchmark_level','root'),getattr(args,'benchmark_quality_policy',None))
    seal(out, 'orthology', panel=str(panel), panel_manifest=digest(_layout_resource(panel/'workflow_manifest.json')), results=str(results))
    print(f'Orthology complete: {out}; inspect calibration_review.tsv before final audit.')


def benchmark_entrypoint(scope='all'):
    if scope not in ('all','adhesome'):
        raise ValueError('Unknown benchmark scope: '+str(scope))
    return ORTH/('adhesome_benchmark_subset.py' if scope=='adhesome' else 'benchmark.py')


def benchmark_command(results,panel,raw,out,logs,refogs=None,input_dir=None,protocol='open_orthobench',scope='all',selection=None,level='root',quality_policy=None):
    command=[sys.executable,benchmark_entrypoint(scope)]
    if scope=='all':command+=['--workflow']
    command+=['--results',results,'--panel',panel,
             '--reference-raw',raw,'--out',out,'--protocol',protocol,'--scope',scope,'--level',level]
    if refogs:command+=['--refogs',refogs]
    if input_dir:command+=['--input-dir',input_dir]
    if selection:command+=['--selection',selection]
    if quality_policy:command+=['--quality-policy',quality_policy]
    execute(command,logs)


def benchmark(args):
    orthology_dir=_layout_resource(args.orthology).resolve()
    meta=verify(orthology_dir,'orthology',required_files=['candidates/candidate_report.tsv','candidates/provenance.json'])
    panel=_layout_resource(Path(meta['panel']))
    verify(panel,'panel')
    if digest(_layout_resource(panel/'workflow_manifest.json'))!=meta['panel_manifest']:
        raise ValueError('Benchmark orthology panel changed')
    out=_layout_resource(args.out).resolve()
    benchmark_command(_layout_resource(Path(meta['results'])),panel,_layout_resource(args.reference_raw).resolve(),out,_layout_resource(out.parent/'logs'),
                      getattr(args,'refogs',None),getattr(args,'input_dir',None),getattr(args,'protocol','open_orthobench'),
                      getattr(args,'scope','all'),getattr(args,'selection',None),getattr(args,'level','root'),getattr(args,'quality_policy',None))
    seal(out,'evolutionary_benchmark',orthology=str(orthology_dir),orthology_manifest=digest(_layout_resource(orthology_dir/'workflow_manifest.json')),
         panel_manifest=digest(_layout_resource(panel/'workflow_manifest.json')),benchmark_script_sha256=digest(_layout_resource(ORTH/'benchmark.py')),
         benchmark_adapter_sha256=digest(_layout_resource(ORTH/'benchmark_workflow.py')),
         benchmark_entrypoint=benchmark_entrypoint(getattr(args,'scope','all')).name,
         benchmark_entrypoint_sha256=digest(benchmark_entrypoint(getattr(args,'scope','all'))))


def benchmark_evidence(folder,orthology_dir,orthology_meta,panel):
    """Only apply a score bound to unchanged ground truth, predictions and panel."""
    result=dict(Orthology_benchmark_status='unavailable',Orthology_benchmark_weight='',
                Orthology_benchmark_source='',Orthology_benchmark_reason='No evolutionary benchmark report',
                Orthology_benchmark_scope='reference_cohort_orthogroup_inference',Orthology_benchmark_f1='',
                Orthology_benchmark_precision='',Orthology_benchmark_recall='',Orthology_benchmark_refogs='',
                Orthology_benchmark_acceptance='not_assessed')
    if not folder or not (_layout_resource(folder/'summary.json')).exists():
        failure=_layout_resource(folder.parent/'benchmark_failure.json') if folder else None
        if failure and _layout_resource(failure).exists():
            detail=json.loads(_layout_resource(failure).read_text())
            result.update(Orthology_benchmark_status='failed',Orthology_benchmark_reason=detail.get('error','Benchmark execution failed'))
        return result,None
    path=_layout_resource(folder/'summary.json')
    sha=digest(path)
    try:
        if (_layout_resource(folder/'workflow_manifest.json')).exists():
            meta=verify(folder,'evolutionary_benchmark')
            if meta['orthology_manifest']!=digest(_layout_resource(orthology_dir/'workflow_manifest.json')):
                raise ValueError('Benchmark belongs to another orthology analysis')
        elif folder==_layout_resource(orthology_dir/'benchmark'):
            verify(orthology_dir,'orthology',required_files=['benchmark/summary.json'])
        else:
            raise ValueError('Benchmark report lacks workflow provenance')
        data=json.loads(_layout_resource(path).read_text())
        result.update(Orthology_benchmark_status=data['status'],Orthology_benchmark_source=str(path),
                      Orthology_benchmark_scope=data.get('evaluated_scope','reference_cohort_orthogroup_inference'),
                      Orthology_benchmark_acceptance=data.get('quality_acceptance',{}).get('status','not_assessed'),
                      Orthology_benchmark_reason='; '.join(data.get('diagnostics',[])))
        if data.get('metrics'):
            result.update(Orthology_benchmark_f1=data['metrics']['f1'],Orthology_benchmark_precision=data['metrics']['precision'],
                          Orthology_benchmark_recall=data['metrics']['recall'],Orthology_benchmark_refogs=data['metrics']['evaluated_refogs'])
        if data['status'] in ('validated','quality_rejected'):
            if data.get('panel_manifest_sha256')!=digest(_layout_resource(panel/'workflow_manifest.json')):
                raise ValueError('Benchmark panel differs from audit panel')
            prediction=_layout_resource(Path(data['prediction_file']))
            if _layout_resource(Path(orthology_meta['results'])).resolve() not in _layout_resource(prediction).resolve().parents:
                raise ValueError('Benchmark used predictions outside the audited orthology results')
            for name,expected in data['inputs'].items():
                if digest(_layout_resource(Path(name)))!=expected:raise ValueError('Benchmark input changed: '+name)
        if data['status']=='validated' and data.get('quality_acceptance',{}).get('status')!='accepted':
            result.update(Orthology_benchmark_status='quality_not_assessed',Orthology_benchmark_weight='',
                          Orthology_benchmark_reason='Legacy benchmark has no accepted performance assessment; rerun benchmarking with the current policy')
        elif data['status']=='validated':
            weight=float(data['validated_weight'])
            if not math.isfinite(weight) or not 0<=weight<=1:raise ValueError('Benchmark weight must be finite and between zero and one')
            if not math.isfinite(float(data['metrics']['f1'])) or abs(float(data['metrics']['f1'])-weight)>1e-12:
                raise ValueError('Benchmark weight differs from the measured reference-cohort F1')
            if data.get('scope','all')=='adhesome':
                weights={family:float(metrics['f1']) for family,metrics in data['family_metrics'].items()}
                if (not weights or not set(weights)<=set(data['selected_families'])
                        or any(not math.isfinite(v) or not 0<=v<=1 for v in weights.values())):
                    raise ValueError('Invalid family-scoped benchmark weights')
                result['Orthology_benchmark_family_weights']=json.dumps(weights,sort_keys=True)
                result['Orthology_benchmark_reason']=('; '.join(data.get('diagnostics',[]))+ '; ' if data.get('diagnostics') else '')+'Weights apply only to exact curated families; other families are unvalidated'
            else:
                result['Orthology_benchmark_weight']=weight
            result.update(Orthology_benchmark_f1=data['metrics']['f1'],Orthology_benchmark_precision=data['metrics']['precision'],
                          Orthology_benchmark_recall=data['metrics']['recall'],Orthology_benchmark_refogs=data['metrics']['evaluated_refogs'])
    except (ValueError,OSError,KeyError,TypeError) as exc:
        result.update(Orthology_benchmark_status='invalid_provenance',Orthology_benchmark_weight='',Orthology_benchmark_reason=str(exc))
    return result,sha


def audit(args):
    registry=InputRegistry()
    for name in ('negative_exclusion_fasta','negative_exclusion_labels'):registry.require(name,'audit.negative_exclusion')
    context_weight=validate_weight(getattr(args,'orthogroup_context_weight',None))
    coverage=getattr(args,'family_coverage',DEFAULT_FAMILY_COVERAGE)
    from stage2_evidence import CONFIG as stage2_config,configuration as stage2_configuration
    stage2_profile=getattr(args,'stage2_profile',default_profile())
    if stage2_profile not in stage2_configuration()['profiles']:raise ValueError('Unknown Stage 2 profile')
    selected_rules,selected_identity=audit_configuration(coverage)
    annotation,orthology_dir=_layout_resource(args.annotation).resolve(),_layout_resource(args.orthology).resolve()
    ann=verify(annotation,'annotation')
    ort=verify(orthology_dir,'orthology',required_files=['candidates/candidate_report.tsv','candidates/provenance.json'])
    panel=_layout_resource(Path(ort['panel']));pan=verify(panel,'panel');handoff=_layout_resource(Path(ann['handoff']));hand=verify(handoff,'handoff')
    if ann['handoff_manifest']!=pan['handoff_manifest'] or digest(_layout_resource(handoff/'workflow_manifest.json'))!=ann['handoff_manifest']:
        raise ValueError('Annotation and orthology must use the same unchanged stage 1 handoff')
    if digest(_layout_resource(panel/'workflow_manifest.json'))!=ort['panel_manifest']:raise ValueError('Orthology panel manifest changed')
    rows=read(_layout_resource(handoff/'hypotheses.tsv'));original_count=len(rows)
    allowed={(r.get('Specie') or r['protein_id'].split('__',1)[0],r['protein_id']) for r in rows}
    seen={(sp,pid,r['family']) for r in rows for sp,pid in [(r.get('Specie') or r['protein_id'].split('__',1)[0],r['protein_id'])]}
    sequence_hashes={r['protein_id']:r.get('sequence_sha256','') for r in rows}
    source_tables={}
    if hand.get('detection'):
        for sp in sorted({key[0] for key in allowed}):
            table=_layout_resource(Path(hand['detection'])/sp/'candidates/candidates.csv')
            expected=hand.get('sources',{}).get(str(_layout_resource(table).resolve()))
            if expected is None:continue
            if digest(table)!=expected:raise ValueError('Stage 1 source table changed: '+str(table))
            source_tables[str(_layout_resource(table).resolve())]=expected
            for index,source_row in enumerate(read(table,','),2):
                pid=source_row['protein_id'];key=(sp,pid,source_row['family'])
                if (sp,pid) not in allowed or key in seen:continue
                seen.add(key)
                rows.append(dict(source_row,Specie=sp,sequence_sha256=sequence_hashes[pid],
                    Stage1_source_table=str(table),Stage1_source_row=index,motif_family=source_row['family']))
    benchmark_folder=getattr(args,'benchmark',None)
    if benchmark_folder is None:
        benchmark_folder=next((p for p in (_layout_resource(orthology_dir/'benchmark'),_layout_resource(orthology_dir.parent/'benchmark')) if (_layout_resource(p/'summary.json')).exists()),_layout_resource(orthology_dir/'benchmark'))
    benchmark_fields,benchmark_sha=benchmark_evidence(benchmark_folder,orthology_dir,ort,panel)
    from audit_join import join
    rows,support=join(annotation,handoff,panel,orthology_dir,rows,benchmark_fields)
    if {(r.get('Specie') or r['protein_id'].split('__',1)[0],r['protein_id']) for r in rows}!=allowed:
        raise ValueError('Stage 2 evidence join changed the stage 1 protein pool')
    out=_layout_resource(args.out).resolve();out.mkdir(parents=True,exist_ok=False)
    configuration=_layout_resource(out/'configuration');configuration.mkdir()
    dump(_layout_resource(configuration/'orthogroup_context.json'),context_configuration(context_weight))
    shutil.copyfile(stage2_config,_layout_resource(configuration/'stage2_evidence_rules.json'))
    rules=_layout_resource(configuration/'family_rules.json');identities=_layout_resource(configuration/'family_identity_rules.json')
    shutil.copyfile(selected_rules,rules);shutil.copyfile(selected_identity,identities)
    print(f"Audit family coverage: {coverage} ({len(json.loads(_layout_resource(rules).read_text())['families'])} family configurations; "
          f"{len(json.loads(_layout_resource(identities).read_text())['rules'])} identity rules).",flush=True)
    write(_layout_resource(out/'joined_evidence.tsv'),rows)
    dump(_layout_resource(out/'supporting_data.json'),support)
    dump(_layout_resource(out/'stage1_join_provenance.json'),dict(handoff_hypotheses=original_count,retained_hypotheses=len(rows),
        selected_species_proteins=len(allowed),source_tables=source_tables,protein_pool_preserved=True))
    from openpyxl import Workbook
    book=Workbook();sheet=book.active;sheet.title='candidates'
    fields=list(dict.fromkeys(key for row in rows for key in row if key!='NCBI_CDD_hits'));sheet.append(fields)
    for row in rows:sheet.append([row.get(key) for key in fields])
    book.save(_layout_resource(out/'audit_input.xlsx'))
    audit_command=[sys.executable,AUDIT,_layout_resource(out/'audit_input.xlsx'),_layout_resource(out/'Stage_2_output.xlsx'),'--rules',rules,
             '--identity-rules',identities,'--supporting-data',_layout_resource(out/'supporting_data.json')]
    if context_weight is not None:audit_command+=['--orthogroup-context-weight',context_weight]
    audit_command+=['--stage2-profile',stage2_profile]
    if getattr(args,'phylogenetic_identity_resolution',False):
        audit_command+=['--phylogenetic-identity-resolution']
        if getattr(args,'phylogenetic_identity_evidence',None):audit_command+=['--phylogenetic-identity-evidence',args.phylogenetic_identity_evidence]
    execute(audit_command,_layout_resource(out/'logs'))
    validation=(dict(status='isolated_branch_evaluation',branch=os.environ['SCHISTO_ISOLATED_BRANCH'])
                if os.environ.get('SCHISTO_ISOLATED_BRANCH') else evaluate_controls(out,handoff))
    motif_snapshot=_layout_resource(MOTIF/'Motif_library/adhesion_core.json')
    seal(out,'audit',annotation=str(annotation),orthology=str(orthology_dir),benchmark_summary_sha256=benchmark_sha,
        phylogenetic_identity_resolution=bool(getattr(args,'phylogenetic_identity_resolution',False)),
        phylogenetic_identity_evidence_sha256=(digest(args.phylogenetic_identity_evidence) if getattr(args,'phylogenetic_identity_resolution',False) and getattr(args,'phylogenetic_identity_evidence',None) and _layout_resource(args.phylogenetic_identity_evidence).is_file() else None),
        stage2_profile=stage2_profile,stage2_evidence_code_sha256=digest(_layout_resource(ROOT/'scripts/stage2_evidence.py')),
        stage2_evidence_rules_sha256=digest(stage2_config),
        family_coverage=coverage,family_rules_source=str(selected_rules),identity_library_source=str(selected_identity),
        orthogroup_context_weight=context_weight,orthogroup_context_configuration=context_configuration(context_weight),
        orthogroup_context_code_sha256=digest(_layout_resource(ROOT/'scripts/orthogroup_context.py')),
        motif_snapshot_sha256=digest(motif_snapshot),rules_sha256=digest(rules),audit_script_sha256=digest(AUDIT),
        audit_policy_sha256=digest(AUDIT_POLICY),audit_evidence_sha256=digest(AUDIT_EVIDENCE),
        classification_contract_sha256=digest(_layout_resource(ROOT/'scripts/classification_contract.py')),output_schema_sha256=digest(_layout_resource(ROOT/'templates/output_schema.json')),stage1_source_tables=source_tables,
        identity_code_sha256=digest(IDENTITY_CODE),identity_library_sha256=digest(identities),audit_report_sha256=digest(AUDIT_REPORT),
        identity_probability_resources=probability_resources(),
        validation_evaluation=validation,validation_evaluation_resources=validation_resources(),
        annotation_manifest=digest(_layout_resource(annotation/'workflow_manifest.json')),orthology_manifest=digest(_layout_resource(orthology_dir/'workflow_manifest.json')))
    print(f'Stage 2 family refinement complete: {out}/Stage_2_output.xlsx')
    if validation['status']=='complete':
        print(f"Validation target performance: {validation['directory']} (F1={validation['F1']}, PPV={validation['PPV']})",flush=True)
    elif validation['status']=='failed':
        print('Validation evaluation failed; completed audit retained: '+validation['error'],flush=True)


def refine_unresolved(args):
    """Run the independent post hoc adapter against an explicit completed audit."""
    script=ROOT.parent/'schisto_pfam_cdd_reclassify/schisto_pfam_cdd_reclassify.py'
    if not script.is_file():
        raise ValueError('Missing sibling post hoc refinement bundle: '+str(script))
    command=[sys.executable,'-B',str(script),'--input',str(args.input),'--sheet',args.sheet]
    if args.out is not None:command+=['--outdir',str(args.out)]
    if args.check_inputs:command+=['--check-inputs']
    if args.compare:command+=['--compare',str(args.compare)]
    if args.strict:command+=['--strict']
    if subprocess.run(command).returncode:
        raise ValueError('Post hoc refinement failed; inspect the preceding input/review error')


def report_refinement(args):
    """Evaluate saved proposals against labels without rerunning the audit/review."""
    script=ROOT.parent/'schisto_pfam_cdd_reclassify/report_posthoc.py'
    if not script.is_file():raise ValueError('Missing post hoc reporting script: '+str(script))
    command=[sys.executable,'-B',str(script),'--refinement',str(args.refinement),
        '--labels',str(args.labels),'--sequence-manifest',str(args.sequence_manifest),
        '--kind',args.kind,'--label-sheet',args.label_sheet,'--out',str(args.out)]
    if args.aliases:command+=['--aliases',str(args.aliases)]
    if subprocess.run(command).returncode:
        raise ValueError('Post hoc performance reporting failed; inspect the preceding input/report error')


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--version', action='version', version='adhesome-workflow '+VERSION)
    sub = p.add_subparsers(dest='command', required=True)
    s = sub.add_parser('start', help='Run/reuse detection and export supported/provisional hypotheses')
    s.add_argument('--detection', type=Path, help='Existing batch detection output root')
    s.add_argument('--input-dir', type=Path, default=raw_input('target_raw'))
    s.add_argument('--targets-config', type=Path, help='JSON target proteome configuration; required for new detection')
    s.add_argument('--check-inputs',action='store_true',help='Validate all target file/header/gene mappings and exit without creating a run')
    s.add_argument('--resume',action='store_true',help='Continue a failed start in the same directory; reuse completed targets without backups')
    s.add_argument('--exclude-target',action='append',help='Omit an unfinished target from detection and handoff export')
    s.add_argument('--threads', type=int, default=8)
    s.set_defaults(func=start)
    s = sub.add_parser('resume', help='Import external context and run motif annotation')
    s.add_argument('--motif-scope',choices=['adhesion','master'],default='adhesion')
    s.add_argument('--motif-library',type=Path,help='Current JSON/XLSX snapshot override')
    s.add_argument('--handoff', type=Path, required=True)
    for flag in ('topology-file','context-dir','external-evidence','evidence'):
        s.add_argument('--'+flag, type=Path)
    s.add_argument('--allow-missing-context', action='store_true')
    s.add_argument('--annotation-only',action='store_true',help='Stop after motif annotation instead of continuing automatically')
    s.add_argument('--retry-failed',action='store_true',help='Retry a failed annotation before motif output in the same directory; verify unchanged inputs and preserve logs without backups')
    s.set_defaults(func=resume)
    s = sub.add_parser('continue',help='Run reference preparation, orthology, evolutionary benchmarking and audit; phylogeny is opt-in')
    s.add_argument('--annotation',type=Path,required=True)
    s.set_defaults(func=downstream)
    for name in ('resume','continue'):
        command=sub.choices[name]
        command.add_argument('--reference-raw',type=Path,default=raw_input('reference_raw'))
        command.add_argument('--reference-file',action='append',help='Select relative reference FASTA paths; repeat; unselected proteomes are not loaded')
        command.add_argument('--reference-panel',help='Named explicit file list in config/input_registry.json')
        command.add_argument('--benchmark-refogs',type=Path,help='Override automatic RefOGs discovery under reference_raw')
        command.add_argument('--benchmark-input-dir',type=Path,help='Override automatic matching benchmark Input discovery')
        command.add_argument('--benchmark-protocol',choices=['open_orthobench','custom'],default='open_orthobench')
        command.add_argument('--benchmark-scope',choices=['all','adhesome'],default='all')
        command.add_argument('--benchmark-selection',type=Path,help='Curated adhesome RefOG/family TSV; used with --benchmark-scope adhesome')
        command.add_argument('--benchmark-level',choices=['root','cohort-lca'],default='root')
        command.add_argument('--benchmark-quality-policy',type=Path,help='Override bundled precision >=0.80 / F1 >=0.70 acceptance criteria')
        command.add_argument('--threads',type=int,default=8)
        command.add_argument('--analysis-threads',type=int,default=8)
        command.add_argument('--compare',type=Path,help='Previous orthology output for reference-panel comparison')
        phylo=command.add_mutually_exclusive_group()
        phylo.add_argument('--phylogeny',action='store_true',help='Opt in to optional trees after the audit; tree failures do not invalidate the completed audit')
        phylo.add_argument('--skip-phylogeny',dest='phylogeny',action='store_false',help='Compatibility alias for the default: do not run phylogeny')
        command.set_defaults(phylogeny=False)
        command.add_argument('--phylogeny-grouping',choices=['orthogroup','family'],default='orthogroup',help='One tree per orthogroup or union of groups per candidate family')
        command.add_argument('--iqtree-threads',type=int,default=2)
        command.add_argument('--allow-partial-phylogeny',action='store_true',help='Retain the tree tool partial-failure summary; audit completion is always independent')
    s = sub.add_parser('calibrate', help='Prepare a new complete gene-level reference panel')
    s.add_argument('--handoff', type=Path, required=True)
    s.add_argument('--panel', type=Path, required=True)
    s.set_defaults(func=calibrate)
    s = sub.add_parser('orthology', help='Run OrthoFinder or import results, then extract and review support')
    s.add_argument('--panel', type=Path, required=True)
    s.add_argument('--results', type=Path)
    s.add_argument('--input-manifest', type=Path)
    s.add_argument('--compare', type=Path)
    s.add_argument('--threads', type=int, default=16)
    s.add_argument('--analysis-threads', type=int, default=8)
    s.add_argument('--reference-raw',type=Path,default=raw_input('reference_raw'))
    s.add_argument('--benchmark-refogs',type=Path)
    s.add_argument('--benchmark-input-dir',type=Path)
    s.add_argument('--benchmark-protocol',choices=['open_orthobench','custom'],default='open_orthobench')
    s.add_argument('--benchmark-scope',choices=['all','adhesome'],default='all')
    s.add_argument('--benchmark-selection',type=Path)
    s.add_argument('--benchmark-level',choices=['root','cohort-lca'],default='root')
    s.add_argument('--benchmark-quality-policy',type=Path)
    s.set_defaults(func=orthology)
    s=sub.add_parser('benchmark',help='Evaluate completed orthology against supplied reference-proteome RefOGs')
    s.add_argument('--orthology',type=Path,required=True)
    s.add_argument('--reference-raw',type=Path,default=raw_input('reference_raw'))
    s.add_argument('--refogs',type=Path)
    s.add_argument('--input-dir',type=Path)
    s.add_argument('--protocol',choices=['open_orthobench','custom'],default='open_orthobench')
    s.add_argument('--scope',choices=['all','adhesome'],default='all')
    s.add_argument('--selection',type=Path,help='Curated RefOG/family TSV; default: reference_raw/adhesome_refog_selection.tsv')
    s.add_argument('--level',choices=['root','cohort-lca'],default='root')
    s.add_argument('--quality-policy',type=Path)
    s.set_defaults(func=benchmark)
    s = sub.add_parser('audit', help='Join all evidence and run conservative family audit')
    s.add_argument('--annotation', type=Path, required=True)
    s.add_argument('--orthology', type=Path, required=True)
    s.add_argument('--benchmark',type=Path,help='Sealed benchmark output; standalone orthology uses its own benchmark subdirectory by default')
    s.set_defaults(func=audit)
    for name in ('resume','continue','orthology'):
        sub.choices[name].add_argument('--orthofinder-inflation',type=float,help='Explicit MCL inflation for new inference (1-10); unset preserves the installed OrthoFinder default')
        sub.choices[name].add_argument('--orthofinder-search',choices=['diamond','diamond_ultra_sens','blast','mmseqs'])
        sub.choices[name].add_argument('--orthofinder-method',choices=['msa','dendroblast'])
    for name in ('resume','continue','audit'):
        sub.choices[name].add_argument('--phylogenetic-identity-resolution',action='store_true',help='Optional review using existing tree placements; independent of --phylogeny inference')
        sub.choices[name].add_argument('--phylogenetic-identity-evidence',type=Path,help='Supported-clade TSV; ignored unless identity resolution is enabled')
        from stage2_evidence import configuration as stage2_configuration
        sub.choices[name].add_argument('--stage2-profile',choices=list(stage2_configuration()['profiles']),default=default_profile(),
            help='Stage 2 profile; default calibrated_v1; precision_v1 is development-only. Saved recovery checkpoints retain their recorded profile.')
        sub.choices[name].add_argument('--orthogroup-context-weight',type=float,default=None,
            help='Optional non-probabilistic orthogroup-context coefficient (0-1); omitted preserves the expanded baseline; 0 disables context; proposed heuristic: 0.541')
        sub.choices[name].add_argument('--family-coverage',choices=['original','expanded'],default=DEFAULT_FAMILY_COVERAGE,
            help='Default: current canonical family/identity pair; original uses the detection family catalogue')
    for s in sub.choices.values():
        s.add_argument('--out', type=Path, required=True, help='New output directory')
    s=sub.add_parser('refine-unresolved',help='Optional post hoc Pfam/CDD review of a completed Stage 2 audit')
    s.add_argument('--input',type=Path,required=True,help='Stage_2_output.xlsx or its audit TSV/CSV export')
    s.add_argument('--sheet',default='Audited_family_assignment')
    s.add_argument('--out',type=Path,help='New separate review directory; required unless checking inputs')
    s.add_argument('--check-inputs',action='store_true',help='Validate input schema and accounting; write nothing')
    s.add_argument('--compare',type=Path,help='Prior post hoc review TSV for comparison')
    s.add_argument('--strict',action='store_true',help='Historical frozen reproduction counts only; omit for new audits')
    s.set_defaults(func=refine_unresolved)
    s=sub.add_parser('report-refinement',help='Report saved post hoc proposals against calibration/validation labels')
    s.add_argument('--refinement',type=Path,required=True)
    s.add_argument('--labels',type=Path,required=True)
    s.add_argument('--sequence-manifest',type=Path,required=True)
    s.add_argument('--kind',choices=['calibration','validation'],required=True)
    s.add_argument('--label-sheet',default='control_identity')
    s.add_argument('--aliases',type=Path,help='Explicit reporting-only family/proposal label equivalences')
    s.add_argument('--out',type=Path,required=True,help='New separate report directory')
    s.set_defaults(func=report_refinement)
    s=sub.add_parser('recover',help='Resume a checkpointed failed downstream run in place')
    s.add_argument('--run',type=Path,required=True,help='Failed downstream directory or its parent annotation directory')
    identity_switch=s.add_mutually_exclusive_group()
    identity_switch.add_argument('--phylogenetic-identity-resolution',action='store_true',default=None,help='Optional existing-tree identity review; family audit remains independent')
    identity_switch.add_argument('--skip-phylogenetic-identity-resolution',dest='phylogenetic_identity_resolution',action='store_false',default=None)
    s.add_argument('--phylogenetic-identity-evidence',type=Path)
    s.add_argument('--adopt-existing',action='store_true',help='Validate and checkpoint an older failed run that predates recovery support')
    s.add_argument('--annotation',type=Path,help='Completed annotation directory for adoption of an older failed run')
    s.add_argument('--retry-benchmark',action='store_true',help='Explicitly recompute a benchmark and refresh dependent audit while reusing completed orthology')
    s.add_argument('--benchmark-scope',choices=['all','adhesome'],default=None,help='Override saved benchmark scope when retrying the benchmark')
    s.add_argument('--benchmark-selection',type=Path,help='Override the saved curated RefOG/family TSV')
    s.add_argument('--benchmark-protocol',choices=['open_orthobench','custom'],default=None)
    s.add_argument('--benchmark-level',choices=['root','cohort-lca'],default=None)
    s.add_argument('--benchmark-quality-policy',type=Path)
    for name in ('threads','analysis-threads','iqtree-threads'):
        s.add_argument('--'+name,type=int,default=None)
    phylo=s.add_mutually_exclusive_group()
    phylo.add_argument('--phylogeny',action='store_true',default=None,help='Opt in to or retry optional phylogeny without rerunning completed orthology/audit')
    phylo.add_argument('--skip-phylogeny',dest='phylogeny',action='store_false',default=None)
    s.add_argument('--phylogeny-grouping',choices=['orthogroup','family'],default=None)
    s.set_defaults(func=recover)
    args = p.parse_args(argv)
    try:
        validate_weight(getattr(args,'orthogroup_context_weight',None))
        if hasattr(args,'orthofinder_inflation'):orthofinder_configuration(args)
        for name in ('threads','analysis_threads','iqtree_threads'):
            if getattr(args,name,None) is not None and getattr(args,name)<1:
                raise ValueError('Thread counts must be positive')
        args.func(args)
        return 0
    except (ValueError, OSError, KeyError) as exc:
        print(f'ERROR: {exc}', file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
