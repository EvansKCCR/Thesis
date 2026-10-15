#!/usr/bin/env python3
"""Resumable integration of detection, external context, motifs, orthology and audit.

Python >=3.11. Commands write new directories and never overwrite completed runs.
See WORKFLOW.md for the input contracts and examples.
"""
import argparse
import csv
import hashlib
import json
import re
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent
VERSION = '1.0'
DETECT = ROOT/'01_detection'
MOTIF = ROOT/'03_functional_motif_annotation'
ORTH = ROOT/'04_orthology_phylogenetic_classification/scripts'
AUDIT = ROOT/'05_conservative_adhesome_family_assignment_audit/integrated_adhesome_pipeline.py'
FIBRONECTIN_FAMILIES = {'fibronectin_like', 'fibronectin_III_ligand', 'fibronectin_like_candidate'}
IDENTITY_FIELDS = ('reference_subfamily','subfamily_support_method','reciprocal_reference_subfamily',
                  'phylogeny_concordant','subfamily_profile_calibrated','top_profile_subfamily',
                  'profile_margin_pass','reference_accession','reference_dataset_version','profile_domain_gate_verified')


def digest(path):
    with Path(path).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def dump(path, value):
    Path(path).write_text(json.dumps(value, indent=2)+'\n', encoding='utf-8')


def read(path, delimiter='\t'):
    with Path(path).open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f, delimiter=delimiter))


def write(path, rows, fields=None):
    fields = fields or list(dict.fromkeys(k for row in rows for k in row))
    with Path(path).open('w', encoding='utf-8', newline='') as f:
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
    for line in Path(path).read_text(encoding='utf-8-sig').splitlines():
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
    Path(path).write_text(''.join(f'>{k}\n{v}\n' for k,v in sorted(seqs.items())), encoding='ascii')


def execute(command, logdir):
    logdir.mkdir(parents=True, exist_ok=True)
    n = len(list(logdir.glob('command_*.json')))+1
    record = dict(command=[str(x) for x in command], cwd=str(ROOT), status='running')
    dest = logdir/f'command_{n:03}.json'
    dump(dest, record)
    with (logdir/f'command_{n:03}.log').open('w', encoding='utf-8') as log:
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
    files = {str(p.relative_to(folder)): digest(p) for p in folder.rglob('*') if p.is_file() and p.name != 'workflow_manifest.json'}
    dump(folder/'workflow_manifest.json', dict(workflow_version=VERSION, kind=kind, files=files, **metadata))


def verify(folder, kind, required_files=None):
    data = json.loads((folder/'workflow_manifest.json').read_text(encoding='utf-8'))
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
    expected = None
    if config_path:
        config = json.loads(config_path.read_text(encoding='utf-8'))
        entries = config['proteomes']
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
    out = args.out.resolve()
    if not args.detection:
        sys.path.insert(0,str(DETECT))
        try:
            from adhesome_pipeline import check_inputs
            _,ignored,summary=check_inputs(args.input_dir,config,out/'detection')
        finally:
            sys.path.pop(0)
        if getattr(args,'check_inputs',False):
            print(json.dumps(dict(input_root=str(args.input_dir.resolve()),species=summary,ignored_fastas=ignored),indent=2))
            return
    elif getattr(args,'check_inputs',False):
        raise ValueError('--check-inputs requires raw --input-dir and --targets-config')
    out.mkdir(parents=True, exist_ok=False)
    detection = args.detection.resolve() if args.detection else out/'detection'
    if not args.detection:
        cmd = [sys.executable, DETECT/'scripts/run_detection.py', '--out', detection, '--cpu', args.threads]
        if args.input_dir:
            cmd += ['--input-dir', args.input_dir.resolve()]
        cmd += ['--config', config_path.resolve()]
        execute(cmd, out/'logs')
    hypotheses, sequences, families, sources, species = [], {}, defaultdict(dict), {}, set()
    for table in sorted(detection.glob('*/candidates/candidates.csv')):
        sp = table.parents[1].name
        if not re.fullmatch(r'[A-Za-z][A-Za-z0-9_]*', sp) or '__' in sp:
            raise ValueError(f'Unsafe target species code: {sp}')
        species.add(sp)
        source = table.parents[1]/'prepared'/f'{sp}.faa'
        seqs = fasta(source)
        sources[str(table)] = digest(table)
        sources[str(source)] = digest(source)
        for row in read(table, ','):
            fibronectin = row['family'] in FIBRONECTIN_FAMILIES
            if row['status'] not in ('supported_candidate', 'provisional') and not (fibronectin and row['status']=='unresolved'):
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
    (out/'families').mkdir()
    for family, seqs in families.items():
        write_fasta(out/'families'/f'{family}.faa', seqs)
    write(out/'hypotheses.tsv', hypotheses)
    component_rows=[]
    component_fields=[]
    for table in sorted(detection.glob('*/candidates/component_hypotheses.csv')):
        rows=read(table,',')
        if rows:
            component_fields=list(rows[0])
            component_rows.extend(r for r in rows if r['protein_id'] in sequences)
        sources[str(table)]=digest(table)
    if component_fields:
        write(out/'component_hypotheses.tsv',component_rows,component_fields)
    fibronectin_rows = [r for r in hypotheses if r['family'] in FIBRONECTIN_FAMILIES]
    write(out/'fibronectin_like_hypotheses.tsv', fibronectin_rows, list(hypotheses[0]))
    write_fasta(out/'fibronectin_like_candidates.faa', {r['protein_id']:sequences[r['protein_id']] for r in fibronectin_rows})
    for tier in ('supported', 'provisional'):
        write(out/f'{tier}_hypotheses.tsv', [r for r in hypotheses if r['hypothesis_tier']==tier], list(hypotheses[0]))
    write_fasta(out/'external_candidates.faa', sequences)
    write(out/'external_evidence_template.tsv', [dict(protein_id=k, sequence_sha256=hashlib.sha256(v.encode()).hexdigest(),
          DeepLoc_2_1='', DeepTMHMM='', SignalP_v6_0='', NCBI_CDD_architecture='', source='',
          **{field:'' for field in IDENTITY_FIELDS}) for k,v in sorted(sequences.items())])
    domains = sorted(detection.glob('*/domains/domains.csv'))
    source_fastas = sorted(detection.glob('*/prepared/*.faa'))
    execute([sys.executable, MOTIF/'import_domain_regions.py', '--input-dir', out/'families',
             '--domains', *domains, '--source-fasta', *source_fastas, '--out', out/'verified_regions.tsv'], out/'logs')
    provenance = json.loads((out/'verified_regions.provenance.json').read_text(encoding='utf-8'))
    if provenance['unverified_records']:
        raise ValueError('Unverified domain source sequences; inspect provenance')
    sources.update(provenance['inputs'])
    for name in ('adhesome_pfam_manifest.csv','adhesome_component_crosswalk.csv','family_rules.json'):
        resource=DETECT/'config'/name
        sources[str(resource)]=digest(resource)
    if config_path:
        sources[str(config_path.resolve())] = digest(config_path)
    seal(out, 'handoff', detection=str(detection), sources=sources,
         completed_stages=['detection'],deferred_stages=['topology','localization','phylogeny'],
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
    handoff = args.handoff.resolve()
    verify(handoff, 'handoff')
    if not any((args.topology_file, args.context_dir, args.external_evidence, args.allow_missing_context)):
        raise ValueError('Supply external results or explicitly use --allow-missing-context')
    out = args.out.resolve()
    seqs = fasta(handoff/'external_candidates.faa')
    external = context_rows(args.external_evidence, seqs) if args.external_evidence else []
    out.mkdir(parents=True, exist_ok=False)
    cmd = [sys.executable, MOTIF/'family_motifs.py', '--input-dir', handoff/'families',
           '--regions', handoff/'verified_regions.tsv', '--out', out/'motifs',
           '--motif-scope',getattr(args,'motif_scope','adhesion')]
    if getattr(args,'motif_library',None):cmd+=['--library',args.motif_library.resolve()]
    inputs = {}
    for flag in ('topology_file','context_dir','evidence'):
        value = getattr(args, flag)
        if value:
            paths = value.rglob('*') if value.is_dir() else [value]
            inputs.update({str(p.resolve()):digest(p) for p in paths if p.is_file()})
            if flag == 'context_dir':
                staged = out/'context'
                staged.mkdir()
                aliases = {
                    'DeepLoc_adhesome_candidate_list.csv':'adhesome_DeepLoc_results.csv',
                    'SignalP_adhesome_candidate_list.txt':'adhesome_SiganP_prediction_results.txt',
                    'NCBI_CDD_adhesome_hitdata.txt':'adhesome_NCBI_CDD_hitdata.txt',
                    'deeptmhmm_predicted_topologies.3line':'adhesome_predicted_topologies.3line'}
                for canonical, alternative in aliases.items():
                    source = value/canonical if (value/canonical).exists() else value/alternative
                    if canonical.endswith('.3line') and args.topology_file:
                        continue
                    shutil.copyfile(source, staged/canonical)
                value = staged
            cmd += ['--'+flag.replace('_','-'), value.resolve()]
    if args.external_evidence:
        inputs[str(args.external_evidence.resolve())] = digest(args.external_evidence)
        write(out/'external_evidence.tsv', external)
    execute(cmd, out/'logs')
    seal(out, 'annotation', handoff=str(handoff), handoff_manifest=digest(handoff/'workflow_manifest.json'), sources=inputs,
         context_policy='missing_allowed' if args.allow_missing_context else 'external_supplied')
    print(f'Annotation complete: {out}.')
    if not getattr(args, 'annotation_only', False):
        downstream(argparse.Namespace(annotation=out,out=out/'downstream',
            reference_raw=args.reference_raw,threads=args.threads,analysis_threads=args.analysis_threads,
            compare=args.compare,skip_phylogeny=args.skip_phylogeny,
            iqtree_threads=args.iqtree_threads,allow_partial_phylogeny=args.allow_partial_phylogeny,
            phylogeny_grouping=args.phylogeny_grouping))


def downstream(args):
    annotation=args.annotation.resolve()
    meta=verify(annotation,'annotation')
    handoff=Path(meta['handoff'])
    verify(handoff,'handoff')
    if digest(handoff/'workflow_manifest.json')!=meta['handoff_manifest']:
        raise ValueError('Annotation handoff changed')
    out=args.out.resolve()
    recovering=getattr(args,'recovering',False)
    if recovering and (out/'workflow_manifest.json').exists():
        verify(out,'completed_workflow')
        print(f'Run already completed: {out}')
        return
    out.mkdir(parents=True,exist_ok=recovering)
    checkpoint=out/'recovery.json'
    if recovering:
        state=json.loads(checkpoint.read_text())
        if state['annotation_sha256']!=digest(annotation/'workflow_manifest.json'):
            raise ValueError('Annotation changed since failed run; create a new run')
        audit_meta=out/'final/workflow_manifest.json'
        if state['stages'].get('audit',{}).get('status')=='complete' and audit_meta.exists():
            saved=json.loads(audit_meta.read_text())
            if (saved.get('rules_sha256')!=digest(DETECT/'config/family_rules.json') or saved.get('audit_script_sha256')!=digest(AUDIT)
                    or saved.get('motif_snapshot_sha256')!=digest(MOTIF/'Motif_library/adhesion_core.json')):
                state['stages'].pop('audit',None)
                if state['settings'].get('phylogeny_grouping')=='family':
                    state['stages'].pop('phylogeny_prepare',None);state['stages'].pop('phylogeny_run',None)
                dump(checkpoint,state)
                print('Current audit rules/code changed: refresh audit and dependent family-tree preparation; reuse orthology.',flush=True)
    else:
        state=dict(settings={k:str(v) if isinstance(v,Path) else v for k,v in vars(args).items() if k!='func'},
                   annotation_sha256=digest(annotation/'workflow_manifest.json'),stages={})
        for key in ('annotation','out','reference_raw','compare'):
            if state['settings'].get(key) is not None:
                state['settings'][key]=str(Path(state['settings'][key]).resolve())
        dump(checkpoint,state)
    programs=[]
    if state['stages'].get('orthology',{}).get('status')!='complete' and not (out/'orthology/inference_complete.json').exists():
        programs.append('orthofinder')
    missing=[name for name in programs if not shutil.which(name)]
    if missing:
        raise ValueError('Missing executables for unfinished stages: '+', '.join(missing)+'. Install dependencies, then use recover --run '+str(out))
    def stage(name,folder,action,preserve=False):
        previous=state['stages'].get(name,{})
        if previous.get('status')=='complete':
            for filename,sha in previous['files'].items():
                if digest(out/filename)!=sha:
                    raise ValueError(f'Completed checkpoint changed: {filename}')
            print(f'REUSE {name}',flush=True)
            return
        if folder.exists() and not preserve:
            archive_attempt(folder)
        state['stages'][name]=dict(status='running')
        dump(checkpoint,state)
        try:
            action()
            files={str(p.relative_to(out)):digest(p) for p in folder.rglob('*') if p.is_file()}
            state['stages'][name]=dict(status='complete',files=files)
            dump(checkpoint,state)
        except Exception as exc:
            state['stages'][name]=dict(status='failed',error=str(exc))
            dump(checkpoint,state)
            raise
    stage('references',out/'references',lambda:execute([sys.executable,ORTH/'prepare_references.py','--raw',args.reference_raw.resolve(),
             '--handoff',handoff,'--out',out/'references'],out/'logs'))
    stage('panel',out/'panel',lambda:calibrate(argparse.Namespace(handoff=handoff,panel=out/'references/panel.auto.json',out=out/'panel')))
    stage('orthology',out/'orthology',lambda:orthology(argparse.Namespace(panel=out/'panel',out=out/'orthology',results=None,input_manifest=None,
        threads=args.threads,analysis_threads=args.analysis_threads,compare=args.compare,recovering=recovering)),preserve=True)
    stage('audit',out/'final',lambda:audit(argparse.Namespace(annotation=annotation,orthology=out/'orthology',out=out/'final')))
    print(f'Audit complete and available independently of phylogeny: {out}/final/audited_candidates.xlsx',flush=True)
    candidates=out/'orthology/candidates'
    if not args.skip_phylogeny:
        if read(candidates/'candidate_HOGs.tsv') or getattr(args,'phylogeny_grouping','orthogroup')=='family':
            if state['stages'].get('phylogeny_run',{}).get('status')!='complete':
                missing=[name for name in ('mafft','trimal','iqtree2') if not shutil.which(name)]
                if missing:
                    raise ValueError('Audit completed. Missing phylogeny executables: '+', '.join(missing)+'. Recover after installing them.')
            cmd=[sys.executable,ORTH/'phylogeny.py','prepare','--hogs',candidates/'candidate_HOGs.tsv',
                 '--orthologues',candidates/'direct_orthologues.tsv','--proteomes',out/'panel/proteomes',
                 '--species-list',out/'panel/species.json','--output',out/'phylogeny']
            if getattr(args,'phylogeny_grouping','orthogroup')=='family':
                cmd+=['--audit-workbook',out/'final/audited_candidates.xlsx','--target-species',out/'panel/targets.json',
                      '--candidate-report',candidates/'candidate_report.tsv']
            dup=candidates/'candidate_duplications.tsv'
            cmd+=['--duplications',dup] if dup.exists() else ['--allow-missing-duplications']
            stage('phylogeny_prepare',out/'phylogeny',lambda:execute(cmd,out/'logs'))
            tree_cmd=[sys.executable,ORTH/'phylogeny.py','run','--project',out/'phylogeny','--threads',args.threads,
                      '--iqtree-threads',getattr(args,'iqtree_threads',2)]
            if getattr(args,'allow_partial_phylogeny',False):
                tree_cmd+=['--allow-partial']
            if recovering:
                tree_cmd+=['--retry-failed']
            stage('phylogeny_run',out/'phylogeny',lambda:execute(tree_cmd,out/'logs'),preserve=True)
        else:
            dump(out/'phylogeny_status.json',dict(status='no_candidate_groups',trees_inferred=0))
    else:
        dump(out/'phylogeny_status.json',dict(status='explicitly_skipped'))
    summary=out/'phylogeny/run_summary.json'
    tree_status=json.loads(summary.read_text())['status'] if summary.exists() else 'not_run'
    seal(out,'completed_workflow',annotation=str(annotation),phylogeny_status=tree_status,
         completion_status='completed_with_phylogeny_failures' if tree_status=='partial_failure' else 'complete')
    print(f'Workflow audit ready (phylogeny: {tree_status}): {out}/final/audited_candidates.xlsx')


def archive_attempt(folder):
    """Preserve incomplete outputs beside their stage; never delete previous evidence."""
    number=1
    while folder.with_name(folder.name+f'_attempt_{number:03d}').exists():
        number+=1
    folder.rename(folder.with_name(folder.name+f'_attempt_{number:03d}'))


def recover(args):
    run=args.run.resolve()
    if not (run/'recovery.json').exists() and (run/'downstream/recovery.json').exists():
        run=run/'downstream'
    if not (run/'recovery.json').exists() and getattr(args,'adopt_existing',False):
        if (run/'downstream/panel').exists():
            run=run/'downstream'
        annotation=(args.annotation or run.parent).resolve()
        ann=verify(annotation,'annotation')
        if (run/'workflow_manifest.json').exists():
            raise ValueError('Existing run is sealed as complete; do not modify it in place')
        panel=verify(run/'panel','panel')
        if panel['handoff_manifest']!=ann['handoff_manifest']:
            raise ValueError('Existing panel and annotation use different handoffs')
        generated=json.loads((run/'references/panel.auto.json').read_text())
        if generated!=json.loads((run/'panel/panel.snapshot.json').read_text()):
            raise ValueError('Generated reference panel differs from staged panel')
        state=dict(annotation_sha256=digest(annotation/'workflow_manifest.json'),stages={},settings=dict(
            annotation=str(annotation),out=str(run),reference_raw=str(ROOT/'input/reference_raw'),
            threads=args.threads or 8,analysis_threads=args.analysis_threads or 8,iqtree_threads=args.iqtree_threads or 2,
            compare=None,skip_phylogeny=False,allow_partial_phylogeny=False))
        def adopt(name,files):
            state['stages'][name]=dict(status='complete',files={str(p.relative_to(run)):digest(p) for p in files if p.is_file()})
        adopt('references',(run/'references').rglob('*'))
        adopt('panel',(run/'panel').rglob('*'))
        if (run/'orthology/workflow_manifest.json').exists():
            orth=verify(run/'orthology','orthology')
            if orth['panel_manifest']!=digest(run/'panel/workflow_manifest.json'):
                raise ValueError('Existing orthology used a different panel')
            adopt('orthology',(run/'orthology').rglob('*'))
        phylo=run/'phylogeny/manifest.json'
        if phylo.exists() and 'orthology' in state['stages']:
            preparation=json.loads(phylo.read_text())
            state['settings']['phylogeny_grouping']=preparation.get('grouping','orthogroup')
            prepared_files=[phylo,run/'phylogeny/groups.tsv']
            for path,sha in preparation['inputs'].items():
                if digest(path)!=sha:
                    raise ValueError(f'Phylogeny preparation input changed: {path}')
            for group in preparation['groups']:
                for name,sha in group['files'].items():
                    path=run/'phylogeny'/group['group']/name
                    if digest(path)!=sha:
                        raise ValueError(f'Prepared group changed: {path}')
                    prepared_files.append(path)
            adopt('phylogeny_prepare',prepared_files)
        dump(run/'recovery.json',state)
    if not (run/'recovery.json').exists():
        raise ValueError('No recovery checkpoint found. For a run created before checkpoints were added, use --adopt-existing and --annotation with its completed annotation directory.')
    state=json.loads((run/'recovery.json').read_text())
    settings=state['settings']
    for key in ('annotation','out','reference_raw','compare'):
        if settings.get(key) is not None:
            settings[key]=Path(settings[key])
    settings['out']=run
    for key in ('threads','analysis_threads','iqtree_threads'):
        if getattr(args,key) is not None:
            settings[key]=getattr(args,key)
    settings['recovering']=True
    downstream(argparse.Namespace(**settings))


def calibrate(args):
    handoff = args.handoff.resolve()
    handoff_meta = verify(handoff, 'handoff')
    targets = set(handoff_meta.get('target_species') or [r['Specie'] for r in read(handoff/'hypotheses.tsv')])
    config = json.loads(args.panel.read_text(encoding='utf-8'))
    entries = config['proteomes']
    codes = [e['species'] for e in entries]
    if len(codes)!=len(set(codes)) or not targets < set(codes) or any(not re.fullmatch(r'[A-Za-z][A-Za-z0-9_]*', s) or '__' in s for s in codes):
        raise ValueError('Panel requires unique safe species codes, every handoff target, and at least one reference')
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    (out/'proteomes').mkdir()
    allseqs, sources, counts = {}, {}, []
    for entry in entries:
        sp = entry['species']
        source = (args.panel.resolve().parent/entry['fasta']).resolve()
        if not entry.get('source') or entry.get('gene_level') is not True:
            raise ValueError(f'{sp}: declare gene_level=true and document the proteome source/release')
        seqs = fasta(source)
        if any(not re.fullmatch(re.escape(sp)+r'__[A-Za-z0-9_.-]+', k) for k in seqs):
            raise ValueError(f'{sp}: expected gene-level IDs {sp}__gene_id; prepare/collapse raw isoforms first')
        allseqs.update(seqs)
        write_fasta(out/'proteomes'/f'{sp}.faa', seqs)
        sources[str(source)] = digest(source)
        counts.append(dict(species=sp, genes=len(seqs), source=entry['source']))
    candidates = fasta(handoff/'external_candidates.faa')
    mapping = []
    for pid, seq in candidates.items():
        if allseqs.get(pid) != seq:
            raise ValueError(f'Candidate missing or sequence changed in whole-proteome panel: {pid}')
        mapping.append(dict(candidate_set='workflow',candidate_id=pid,prepared_id=pid,status='matched'))
    write(out/'candidate_audit.tsv', mapping)
    write(out/'proteome_summary.tsv', counts)
    dump(out/'species.json', codes)
    dump(out/'targets.json', sorted(targets))
    dump(out/'panel.snapshot.json', config)
    seal(out, 'panel', handoff=str(handoff), handoff_manifest=digest(handoff/'workflow_manifest.json'), sources=sources,
         interpretation='Reference-panel sensitivity analysis, not independent biological validation')
    print(f'Panel ready: {out}. Rerun orthology for this panel; prior runs remain unchanged.')


def orthology(args):
    panel = args.panel.resolve()
    panel_meta = verify(panel, 'panel')
    out = args.out.resolve()
    recovering=getattr(args,'recovering',False)
    out.mkdir(parents=True, exist_ok=recovering)
    if args.results:
        results = args.results.resolve()
        # Bind external results to the exact input panel, not merely matching species labels.
        if not args.input_manifest or digest(args.input_manifest) != digest(panel/'workflow_manifest.json'):
            raise ValueError('External results require --input-manifest from the panel used for inference')
    else:
        receipt=out/'inference_complete.json'
        if recovering and receipt.exists():
            saved=json.loads(receipt.read_text())
            if saved['panel_manifest']!=digest(panel/'workflow_manifest.json'):
                raise ValueError('Panel changed since inference')
            for name,sha in saved.get('files',{}).items():
                if digest(out/'inference'/name)!=sha:
                    raise ValueError(f'Completed inference output changed: {name}')
        else:
            if (out/'inference').exists():
                archive_attempt(out/'inference')
            execute(['orthofinder','-f',panel/'proteomes','-t',args.threads,'-a',args.analysis_threads,'-o',out/'inference'], out/'logs')
            dump(receipt,dict(panel_manifest=digest(panel/'workflow_manifest.json'),
                files={str(p.relative_to(out/'inference')):digest(p) for p in (out/'inference').rglob('*') if p.is_file()}))
        hits = list((out/'inference').glob('Results_*'))
        if len(hits)!=1:
            raise ValueError('Cannot identify one Results_* directory; use explicit --results on a new output')
        results = hits[0]
    if recovering and (out/'candidates').exists():
        archive_attempt(out/'candidates')
    execute([sys.executable, ORTH/'pipeline.py','extract','--results',results,
             '--candidate-set','workflow','--candidate-audit',panel/'candidate_audit.tsv',
             '--species-list',panel/'species.json','--target-species',panel/'targets.json','--output',out/'candidates'], out/'logs')
    report = read(out/'candidates/candidate_report.tsv')
    reference_columns = [k for k in report[0] if k.endswith('_orthologues')]
    review = []
    for row in report:
        support = [k.removesuffix('_orthologues') for k in reference_columns if row[k]]
        ambiguous = any('many' in v for k,v in row.items() if k.endswith('_relationship'))
        review.append(dict(candidate_id=row['candidate_id'], reference_count=len(support),
                           references=';'.join(support), calibration_review='no_direct_reference' if not support else 'duplication_or_paralogy_review' if ambiguous else 'direct_reference_detected'))
    write(out/'calibration_review.tsv', review)
    if args.compare:
        previous = verify(args.compare.resolve(), 'orthology')
        previous_panel = verify(Path(previous['panel']), 'panel')
        if previous_panel['handoff_manifest'] != panel_meta['handoff_manifest']:
            raise ValueError('Calibration comparison requires the same sequence-verified handoff')
        old = {r['candidate_id']:r for r in read(args.compare/'calibration_review.tsv')}
        if set(old) != {r['candidate_id'] for r in review}:
            raise ValueError('Calibration comparison requires the same candidate set')
        write(out/'calibration_comparison.tsv', [dict(candidate_id=r['candidate_id'],
              previous_references=old[r['candidate_id']]['references'], current_references=r['references'],
              previous_review=old[r['candidate_id']]['calibration_review'], current_review=r['calibration_review']) for r in review])
    seal(out, 'orthology', panel=str(panel), panel_manifest=digest(panel/'workflow_manifest.json'), results=str(results))
    print(f'Orthology complete: {out}; inspect calibration_review.tsv before final audit.')


def audit(args):
    annotation, orthology_dir = args.annotation.resolve(), args.orthology.resolve()
    ann = verify(annotation, 'annotation')
    ort = verify(orthology_dir, 'orthology', required_files=['candidates/candidate_report.tsv','candidates/provenance.json'])
    panel = Path(ort['panel'])
    pan = verify(panel, 'panel')
    handoff = Path(ann['handoff'])
    verify(handoff, 'handoff')
    if ann['handoff_manifest'] != pan['handoff_manifest'] or digest(handoff/'workflow_manifest.json') != ann['handoff_manifest']:
        raise ValueError('Annotation and orthology must use the same unchanged handoff')
    if digest(panel/'workflow_manifest.json') != ort['panel_manifest']:
        raise ValueError('Orthology panel manifest changed')
    rows = read(handoff/'hypotheses.tsv')
    reports = {r['candidate_id']:r for r in read(orthology_dir/'candidates/candidate_report.tsv')}
    sites, contexts, cdd = defaultdict(list), {}, defaultdict(list)
    motif_snapshot=MOTIF/'Motif_library/adhesion_core.json'
    current_motif_ids={r['rule_id'] for r in json.loads(motif_snapshot.read_text())['rules']}
    for r in read(annotation/'motifs/adhesion_motif_annotations.tsv'):
        sites[r['family'],r['sequence_id']].append(r)
    context_path = annotation/'motifs/sequence_context.tsv'
    if context_path.exists():
        contexts = {(r['family'],r['sequence_id']):r for r in read(context_path)}
    cdd_path = annotation/'motifs/CDD_specific_hits.tsv'
    if cdd_path.exists():
        for r in read(cdd_path):
            cdd[r['family'],r['sequence_id']].append(f"{r['short_name']}:{r['start']}-{r['end']}")
    external = {r['protein_id']:r for r in read(annotation/'external_evidence.tsv')} if (annotation/'external_evidence.tsv').exists() else {}
    topology = defaultdict(set)
    for r in read(annotation/'motifs/used_regions.tsv'):
        if r['region'] in ('transmembrane', 'signal_peptide'):
            topology[r['family'],r['sequence_id']].add('TM' if r['region']=='transmembrane' else 'SP')
    names = {'Hsap':'H. sapiens','Cele':'C. elegans','Dmel':'D. melanogaster','Mmus':'M. musculus','Xlae':'X. laevis'}
    for row in rows:
        pid, fam = row['protein_id'], row.get('motif_family') or row['family']
        report = reports[pid]
        refs = [k.removesuffix('_orthologues') for k,v in report.items() if k.endswith('_orthologues') and v]
        row['Orthology_reference_support'] = '; '.join(names.get(s,s) for s in refs)
        row['Orthology_direct_reference_count'] = len(refs)
        row['Orthology_HOG_status'] = report['HOG_status']
        row['Orthology_HOG'] = report['HOG']
        row['Orthology_interpretation'] = 'Direct pairwise evidence only; family-specific reference labels not inferred'
        old_sites=sites[fam,pid]
        ss=[dict(r,candidate_tier='context_conflict' if r.get('candidate_tier')=='context_conflict' else 'sequence_only_candidate')
            for r in old_sites if r.get('motif_id') in current_motif_ids]
        row['MotifPrediction_ignored_superseded_count']=len(old_sites)-len(ss)
        row['MotifPrediction_candidate_tiers'] = ';'.join(sorted({r['candidate_tier'] for r in ss}))
        row['MotifPrediction_context_supported_count'] = sum(r['candidate_tier'] in ('context_supported_candidate','high_confidence_hypothesis') for r in ss)
        row['MotifPrediction_priority_site_count'] = sum(r['candidate_tier']=='sequence_only_candidate' for r in ss)
        ctx, ext = contexts.get((fam,pid), {}), external.get(pid, {})
        for field in IDENTITY_FIELDS:
            if ext.get(field):row[field]=ext[field]
        row['DeepLoc_2.1'] = ext.get('DeepLoc_2_1') or ctx.get('deeploc_localizations','')
        row['SignalP.v6.0'] = ext.get('SignalP_v6_0') or ctx.get('signalp_prediction','')
        row['DeepTMHMM'] = ext.get('DeepTMHMM') or ';'.join(sorted(topology[fam,pid]))
        row['NCBI_CDD_architecture'] = ';'.join(filter(None, [ext.get('NCBI_CDD_architecture',''), *cdd[fam,pid]]))
        row['External_evidence_source'] = ext.get('source','')
        row['domain_match'] = str(row['domain_match']).lower() in ('true','1')
    from openpyxl import Workbook
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    write(out/'joined_evidence.tsv', rows)
    wb = Workbook()
    sh = wb.active
    sh.title = 'candidates'
    fields = list(dict.fromkeys(k for r in rows for k in r))
    sh.append(fields)
    for row in rows:
        sh.append([row.get(k) for k in fields])
    wb.save(out/'audit_input.xlsx')
    audit_rules=DETECT/'config/family_rules.json'
    execute([sys.executable,AUDIT,out/'audit_input.xlsx',out/'audited_candidates.xlsx','--rules',audit_rules,'--keep-all-audited'],out/'logs')
    seal(out, 'audit', annotation=str(annotation), orthology=str(orthology_dir),
         motif_snapshot_sha256=digest(motif_snapshot),
         rules_sha256=digest(audit_rules),audit_script_sha256=digest(AUDIT),
         annotation_manifest=digest(annotation/'workflow_manifest.json'), orthology_manifest=digest(orthology_dir/'workflow_manifest.json'))
    print(f'Workflow complete: {out}/audited_candidates.xlsx')


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--version', action='version', version='adhesome-workflow '+VERSION)
    sub = p.add_subparsers(dest='command', required=True)
    s = sub.add_parser('start', help='Run/reuse detection and export supported/provisional hypotheses')
    s.add_argument('--detection', type=Path, help='Existing batch detection output root')
    s.add_argument('--input-dir', type=Path)
    s.add_argument('--targets-config', type=Path, help='JSON target proteome configuration; required for new detection')
    s.add_argument('--check-inputs',action='store_true',help='Validate all target file/header/gene mappings and exit without creating a run')
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
    s.set_defaults(func=resume)
    s = sub.add_parser('continue',help='Run automatic reference preparation, orthology, phylogeny and audit from an existing annotation')
    s.add_argument('--annotation',type=Path,required=True)
    s.set_defaults(func=downstream)
    for name in ('resume','continue'):
        command=sub.choices[name]
        command.add_argument('--reference-raw',type=Path,default=ROOT/'input/reference_raw')
        command.add_argument('--threads',type=int,default=8)
        command.add_argument('--analysis-threads',type=int,default=8)
        command.add_argument('--compare',type=Path,help='Previous orthology output for reference-panel comparison')
        command.add_argument('--skip-phylogeny',action='store_true')
        command.add_argument('--phylogeny-grouping',choices=['orthogroup','family'],default='orthogroup',help='One tree per orthogroup or union of groups per candidate family')
        command.add_argument('--iqtree-threads',type=int,default=2)
        command.add_argument('--allow-partial-phylogeny',action='store_true',help='Continue to audit while retaining an explicit partial-failure phylogeny summary')
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
    s.set_defaults(func=orthology)
    s = sub.add_parser('audit', help='Join all evidence and run conservative family audit')
    s.add_argument('--annotation', type=Path, required=True)
    s.add_argument('--orthology', type=Path, required=True)
    s.set_defaults(func=audit)
    for s in sub.choices.values():
        s.add_argument('--out', type=Path, required=True, help='New output directory')
    s=sub.add_parser('recover',help='Resume a checkpointed failed downstream run in place')
    s.add_argument('--run',type=Path,required=True,help='Failed downstream directory or its parent annotation directory')
    s.add_argument('--adopt-existing',action='store_true',help='Validate and checkpoint an older failed run that predates recovery support')
    s.add_argument('--annotation',type=Path,help='Completed annotation directory for adoption of an older failed run')
    for name in ('threads','analysis-threads','iqtree-threads'):
        s.add_argument('--'+name,type=int,default=None)
    s.set_defaults(func=recover)
    args = p.parse_args(argv)
    try:
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
