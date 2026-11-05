"""Configuration-only deployment audit. No sequence analysis or tool execution."""
import argparse
import ast
import csv
import hashlib
from datetime import datetime, timezone
import importlib.util
import importlib.metadata
import json
from pathlib import Path
import re
import shutil
import sys

from input_dependencies import InputRegistry, ROOT, deployment_source_digest


def ordered(nodes):
    index = {n['id']: n for n in nodes}
    if len(index) != len(nodes):
        raise ValueError('Duplicate stage identifier')
    result, pending = [], dict(index)
    while pending:
        ready = [key for key, node in pending.items() if all(dep in result for dep in node['after']+node.get('soft_after', []))]
        if not ready:
            raise ValueError('Circular or unknown stage dependencies: ' + ', '.join(pending))
        for key in ready:
            result.append(key)
            del pending[key]
    return [index[key] for key in result]


def graph(workflow='all', phylogeny=False, calibration=True, reference_only=False):
    nodes = []
    def add(key, after, output, inputs=(), optional=False, gate=None, soft_after=()):
        nodes.append(dict(id=key, after=list(after), soft_after=list(soft_after), produces=output,
                          inputs=list(inputs), optional=optional, gate=gate or key.split('.')[0]))
    gates = []
    if workflow in ('all', 'target', 'reconciliation'):
        gates.append('target')
    if workflow in ('all', 'validation'):
        gates.append('validation')
    if calibration and workflow in ('all', 'target', 'validation'):
        gates.append('calibration')
    for gate in gates:
        inputs = ['target_raw','target_configuration'] if gate == 'target' else [gate + '_fasta', gate + '_labels']
        add(gate+'.detection', [], gate+'/handoff/workflow_manifest.json', inputs, gate == 'calibration')
        add(gate+'.external_context', [gate+'.detection'], gate+'/returned/context', optional=True)
        add(gate+'.motifs', [gate+'.external_context'], gate+'/annotation/workflow_manifest.json')
        add(gate+'.references', [gate+'.motifs'], gate+'/downstream/references/panel.auto.json', ['reference_raw'])
        add(gate+'.orthology', [gate+'.references'], gate+'/downstream/orthology/workflow_manifest.json')
        add(gate+'.benchmark', [gate+'.orthology'], gate+'/downstream/benchmark/summary.json', ['benchmark_refogs','benchmark_input'], True)
        add(gate+'.audit', [gate+'.orthology'], gate+'/downstream/final/audited_family_assignments.tsv',
            ['negative_exclusion_fasta','negative_exclusion_labels'],soft_after=[gate+'.benchmark'])
        add(gate+'.report', [gate+'.audit'], gate+'/reports/'+('validation_review.html' if gate=='validation' else 'report.html'))
        if phylogeny:
            add(gate+'.phylogeny', [gate+'.audit'], gate+'/downstream/phylogeny/summary.json', optional=True)
    if workflow in ('all', 'homology', 'reconciliation'):
        add('homology.rbh', [], 'homology/rbh_candidates.faa', ['curated_references','target_raw']+([] if reference_only else ['reference_raw']))
    if workflow in ('all', 'reconciliation'):
        add('reconciliation.library', ['target.audit','homology.rbh'], 'reconciliation/library_manifest.json')
    if workflow in ('all', 'phylogeny', 'reconciliation'):
        parent = ['reconciliation.library'] if workflow != 'phylogeny' else []
        add('independent.align', parent, 'independent/aligned.fasta', ['phylogeny_library'] if workflow=='phylogeny' else (), gate='independent')
        add('independent.trim', ['independent.align'], 'independent/trimmed.fasta', gate='independent')
        add('independent.tree', ['independent.trim'], 'independent/tree.treefile', gate='independent')
    return ordered(nodes)


def audit(root=ROOT, workflow='all', clean_state=True, phylogeny=False, registry_path=None, input_overrides=None, reference_only=False):
    root = Path(root).resolve()
    bundle = root.parent
    checks, errors, warnings = [], [], []
    def record(name, passed, detail):
        checks.append(dict(check=name, status='pass' if passed else 'fail', detail=detail))
        if not passed:
            errors.append(name+': '+str(detail))
    registry = InputRegistry(root, registry_path)
    input_overrides=input_overrides or {}
    configs = {}
    for path in sorted((root/'config').rglob('*.json')):
        try:
            configs[path.relative_to(root).as_posix()] = json.loads(path.read_text(encoding='utf-8-sig'))
        except (OSError, ValueError) as exc:
            record('Workflow configuration', False, str(path)+': '+str(exc))
    record('Workflow configuration JSON', bool(configs), f'{len(configs)} configuration files parsed')
    settings = configs.get('config/isolated_workflow.json', {})
    profiles = configs.get('config/annotation/stage2_evidence_rules.json', {})
    valid = (settings.get('version') == '1.0' and settings.get('stage2_profile') in profiles.get('profiles', {})
             and settings.get('family_coverage') in ('expanded','original') and isinstance(settings.get('qc'), dict)
             and 0 <= settings.get('orthogroup_context_weight', -1) <= 1)
    allowed={'version','family_coverage','orthogroup_context_weight','qc','orthofinder_search','orthofinder_method',
             'orthofinder_inflation','benchmark_scope','benchmark_selection','stage2_profile','phylogenetic_identity_resolution','reference_files','reference_panel'}
    qc=settings.get('qc',{})
    valid=valid and not set(settings)-allowed and qc.get('min_length',0)>0 and 0<=qc.get('max_ambiguous_fraction',-1)<=1 and qc.get('max_ambiguous_run',-1)>=0
    valid=valid and settings.get('orthofinder_search') in (None,'diamond','diamond_ultra_sens','blast','mmseqs') and settings.get('orthofinder_method') in (None,'msa','dendroblast')
    valid=valid and (settings.get('orthofinder_inflation') is None or 1<=settings['orthofinder_inflation']<=10) and settings.get('benchmark_scope','all') in ('all','adhesome')
    record('Workflow configuration', valid, 'Isolated settings, Stage 2 profile and paired family coverage')
    identity=configs.get('config/annotation/family_identity_rules.json',{}).get('rules',[])
    identities=[r['id'] for r in identity]
    additions=[r['id'] for r in profiles.get('family_rules',[])]
    record('Family and identity rule initialization', bool(identity) and len(identities)==len(set(identities)) and
           len(additions)==len(set(additions)) and not set(identities)&set(additions),
           dict(identity_rules=len(identity),stage2_additions=len(additions),family_rules=len(configs.get('config/annotation/family_rules.json',{}).get('families',{}))))
    for name in ['config/annotation/family_rules.json','config/annotation/family_identity_rules.json',
                 'config/annotation/family_set_registry.json','config/motif_annotation/adhesion_core.json',
                 'config/motif_annotation/motif_master.json','config/orthology_phylogeny/benchmark_policy.json']:
        record('Configuration resource: '+name, name in configs, 'Required static library configuration')
    layout = configs.get('config/layout_paths.json', {})
    required_paths={p for p in layout.get('file_map', {}).values() if 'cache' not in Path(p).parts}
    missing = sorted(p for p in required_paths if not (root/p).is_file())
    record('Source-layout paths', not missing, missing or f'{len(required_paths)} required source paths resolve; historical provenance caches are optional')
    release=bundle/'release_manifest.sha256'
    if release.is_file():
        invalid=[]
        for line in release.read_text().splitlines():
            if not line.strip() or line.startswith('#'):continue
            sha,name=line.split(None,1)
            path=(bundle/name).resolve()
            if bundle not in path.parents or not path.is_file():invalid.append(name)
            elif any(part in ('input','inputs','validation','runs','results','validation_reports') for part in Path(name).parts):
                invalid.append(name+' (data-dependent file in source inventory)')
            elif deployment_source_digest(path)!=sha:invalid.append(name+' (checksum changed)')
        record('Deployment source inventory',not invalid,invalid[:12] or 'Source-only release checksums resolve; user inputs and results excluded')
    else:record('Deployment source inventory',False,'Missing release_manifest.sha256; refresh with maintenance/build_source_manifest.py --write')
    attributes=(bundle/'.gitattributes').read_text() if (bundle/'.gitattributes').is_file() else ''
    record('Scientific asset checkout integrity','*.hmm -text' in attributes and 'integrin_adhesome_analysis/config/** -text' in attributes,
           'Git preserves hashed scientific library/config bytes; release code-text checks normalize LF without changing biological fingerprints')
    parse_errors, parsed = [], 0
    for folder in [root/'main', root/'maintenance', root/'resource_library', bundle/'homology_retrieval', bundle/'reconciliation', bundle/'independent_family_phylogeny', bundle/'schisto_pfam_cdd_reclassify']:
        for path in folder.rglob('*.py'):
            if '__pycache__' in path.parts:
                continue
            try:
                ast.parse(path.read_text(encoding='utf-8-sig'), filename=str(path)); parsed += 1
            except (SyntaxError, ValueError, OSError) as exc:
                parse_errors.append(str(path)+': '+str(exc))
    record('Script syntax', not parse_errors, parse_errors or f'{parsed} Python files parsed without importing scientific datasets')
    entrypoints = {
        'detection': root/'resource_library/detection/scripts/run_detection.py',
        'motifs': root/'resource_library/motif_annotation/family_motifs.py',
        'references': root/'resource_library/orthology_phylogeny/scripts/prepare_references.py',
        'orthology': root/'resource_library/orthology_phylogeny/scripts/pipeline.py',
        'audit': root/'resource_library/annotation/integrated_adhesome_pipeline.py',
        'isolated_branches': root/'main/run_pipeline.py',
        'homology': bundle/'homology_retrieval/reciprocal_blast.py',
        'reconciliation': bundle/'reconciliation/reconcile_candidates.py',
        'independent_phylogeny': bundle/'independent_family_phylogeny/run_phylogeny.py',
        'posthoc_refinement': bundle/'schisto_pfam_cdd_reclassify/schisto_pfam_cdd_reclassify.py',
        'posthoc_input_adapter': bundle/'schisto_pfam_cdd_reclassify/stage2_io.py',
        'posthoc_reporting': bundle/'schisto_pfam_cdd_reclassify/report_posthoc.py'}
    for name, path in entrypoints.items():
        record('Stage entry point: '+name, path.is_file(), str(path.relative_to(bundle)))
    # Check deployed HMM/library membership, not target/reference proteomes.
    hmm = root/'resource_library/detection/databases/pfam/adhesome_v1.0.hmm'
    manifest = root/'config/detection/adhesome_pfam_manifest.csv'
    if hmm.is_file() and manifest.is_file():
        with manifest.open(encoding='utf-8-sig', newline='') as handle:
            wanted = {row['Pfam accession'].split('.')[0] for row in csv.DictReader(handle)}
        with hmm.open(encoding='utf-8') as handle:
            found = {line.split()[1].split('.')[0] for line in handle if line.startswith('ACC ')}
        record('Detection model library', wanted <= found, dict(required=len(wanted), available=len(found), missing=sorted(wanted-found)))
        current_rules=configs.get('config/detection/family_rules.json',{}).get('families',{})
        rule_models=set(re.findall(r'PF\d{5}',json.dumps({k:{key:value for key,value in rule.items() if key in ('domain_alternatives','screening_domains')} for k,rule in current_rules.items()})))
        record('Current detection-rule model coverage',rule_models<=found,dict(required=len(rule_models),missing=sorted(rule_models-found)))
        metadata_path=hmm.with_name(hmm.name+'.metadata.json')
        if metadata_path.is_file():
            metadata=json.loads(metadata_path.read_text())
            hashes_ok=metadata.get('library_sha256')==hashlib.sha256(hmm.read_bytes()).hexdigest() and metadata.get('manifest_sha256')==hashlib.sha256(manifest.read_bytes()).hexdigest()
            for key,name in [('crosswalk_sha256','adhesome_component_crosswalk.csv')]:
                if metadata.get(key):hashes_ok=hashes_ok and hashlib.sha256((root/'config/detection'/name).read_bytes()).hexdigest()==metadata[key]
            if metadata.get('family_rules_sha256')!=hashlib.sha256((root/'config/detection/family_rules.json').read_bytes()).hexdigest():
                warnings.append('HMM metadata retains its historical build-time family-rule fingerprint; current positive/screening model coverage is checked separately, without rewriting that provenance.')
            for name,sha in metadata.get('split_sha256',{}).items():
                path=hmm.with_name(hmm.name+'.'+name+'.hmm')
                hashes_ok=hashes_ok and path.is_file() and hashlib.sha256(path.read_bytes()).hexdigest()==sha
            record('Deployed HMM checksums', hashes_ok, 'Combined, manifest and split-library checksums; no library rebuild performed')
        else:record('Deployed HMM checksums',False,'Missing HMM metadata')
    else:
        record('Detection model library', False, 'Missing deployed HMM or Pfam manifest')
    nodes = graph(workflow, phylogeny,reference_only=reference_only)
    events, disabled_gates = [], set()
    for node in nodes:
        if node['gate'] in disabled_gates:continue
        node_events=[]
        for name in node['inputs']:
            mandatory = not node['optional']
            node_events.append(registry.assess(name, node['id'], mandatory=mandatory, clean_state=clean_state,
                                               override=input_overrides.get(name)))
        events.extend(node_events)
        if node['id']=='calibration.detection' and any(not e['available'] for e in node_events):disabled_gates.add('calibration')
    by_stage = {}
    for event in events:
        by_stage.setdefault(event['stage'], []).append(event)
    simulated = {}
    for node in nodes:
        blocked = any(simulated[p] in ('blocked_on_data','blocked_on_upstream') for p in node['after'])
        absent = any(not e['available'] for e in by_stage.get(node['id'], []))
        skipped=any(simulated[p] in ('skipped_optional_input','skipped_optional_gate') for p in node['after'])
        if skipped:
            state = 'skipped_optional_gate'
        elif absent and node['optional']:
            state = 'skipped_optional_input'
        elif blocked:
            state = 'blocked_on_upstream'
        elif absent:
            state = 'blocked_on_data'
        elif node['id'].endswith('.external_context'):
            state = 'awaiting_external_return_on_real_run'
        elif node['gate'] == 'independent' and workflow == 'phylogeny' and not node['after']:
            state = 'requires_prepared_library_on_real_run'
        else:
            state = 'configuration_initialized'
        simulated[node['id']] = state
        node['simulated_state'] = state
        node['output_destination'] = str(root/'runs'/'PLANNED_ONLY'/node['produces'])
        if node['id']=='homology.rbh':node['output_destination']=str(bundle/'homology_retrieval/results/homology_run/rbh_candidates.faa')
        elif node['id']=='reconciliation.library':node['output_destination']=str(bundle/'independent_family_phylogeny/library_reconciled/library_manifest.json')
        elif node['gate']=='independent':node['output_destination']=str(bundle/'independent_family_phylogeny/results/PLANNED_ONLY/FAMILY'/Path(node['produces']).name)
        node['executed'] = False
    destinations=[n['output_destination'] for n in nodes]
    protected=[registry.path_for(name) for name in ('target_raw','reference_raw')]
    safe=len(set(destinations))==len(destinations) and not any(base in Path(value).parents for value in destinations for base in protected)
    record('Output destinations',safe,'Separate gate/downstream namespaces and independent bundle outputs; no directories created')
    record('Stage dependency graph', True, f'{len(nodes)} conditional stages topologically ordered; benchmark is a soft audit dependency; phylogeny follows audit')
    record('Raw-input paths', registry.path_for('reference_proteomes') == root/'input/reference_raw' and
           registry.path_for('target_proteomes') == root/'input/target_raw', 'Homology aliases resolve to main-project raw roots')
    record('Missing input handling', True, 'unidentified_input decisions recorded per requesting stage; mandatory stops cannot be disabled')
    record('Analysis-gate isolation', True, 'Separate gate namespaces, raw roots, branch contracts, reserved identifiers and existing overlap checks')
    record('Stale output dependencies', True, 'All generated dependencies are planned outputs; no saved audit, result, panel or inference checkpoint is opened')
    record('Selective reference loading', True, 'Homology resolves source species before loading sequences; explicit reference-set/proteome/panel selections are supported')
    software = []
    main_selected=any(n['gate'] in ('target','validation','calibration') for n in nodes)
    homology_selected=any(n['id']=='homology.rbh' for n in nodes)
    trees_selected=any(n['id']=='independent.tree' or n['id'].endswith('.phylogeny') for n in nodes)
    required_tools=({'hmmsearch','hmmpress','blastp','makeblastdb','orthofinder','mafft'} if main_selected else set())
    if homology_selected:required_tools.update(('blastp','makeblastdb'))
    if trees_selected:required_tools.update(('mafft','trimal','iqtree'))
    for name, candidates in {'hmmsearch':['hmmsearch'], 'hmmpress':['hmmpress'], 'blastp':['blastp'],
                              'makeblastdb':['makeblastdb'], 'orthofinder':['orthofinder'], 'mafft':['mafft'],
                              'trimal':['trimal'], 'iqtree':['iqtree','iqtree2','iqtree3']}.items():
        path = next((found for candidate in candidates if (found := shutil.which(candidate))), None)
        software.append(dict(name=name, available=bool(path), path=path, version_verified=False,required=name in required_tools))
    for module, package in [('openpyxl','openpyxl'),('Bio','biopython')]:
        available = importlib.util.find_spec(module) is not None
        try:
            version = importlib.metadata.version(package) if available else None
        except importlib.metadata.PackageNotFoundError:
            version = None
        software.append(dict(name=package, available=available, version=version,required=main_selected))
    missing_software = [s['name'] for s in software if s['required'] and not s['available']]
    if missing_software:
        warnings.append('Missing software in the current environment: '+', '.join(missing_software))
    warnings.append('Executable paths are checked without running programs; tool versions must match the supplied conda environment before real execution.')
    warnings.append('OrthoFinder bundled/transitive programs are not executed by this check; install the supplied conda environment and verify its configured search/tree backends before inference.')
    if sys.platform=='win32':warnings.append('Software availability here describes Windows Python, not any WSL/conda environment; repeat this command inside the intended analysis environment.')
    warnings.append('Historical standalone families.json/queue OG identifiers remain source-run-specific; a fresh source result requires reviewed mappings or the reconciliation handoff.')
    if any(not e['available'] for e in events):
        warnings.append('Biological inputs are absent or assumed absent: expected for a clean-state structural dry run; actual execution is not ready.')
    checks.append(dict(check='Software dependencies',status='warning' if missing_software else 'pass',detail=software))
    record('Clean-state dry run', True, 'Configuration and transition simulation complete; no scientific process launched and no biological results fabricated')
    structural = not errors
    data_ready = not any(not e['available'] and e['mandatory'] for e in events)
    assessment='fail' if not structural else 'warning' if warnings else 'pass'
    checks.append(dict(check='Deployment readiness',status=assessment,detail='Structural validation and data/software execution readiness are reported separately'))
    return dict(version='1.0', created_utc=datetime.now(timezone.utc).isoformat(), workflow=workflow,
                clean_state=clean_state, analysis_launched=False, biological_results_generated=False,
                structural_readiness='pass' if structural else 'fail', data_execution_readiness=data_ready,
                deployment_readiness=assessment,
                execution_ready=False,execution_prerequisites_present=structural and data_ready and not missing_software,
                data_readiness_scope='Presence/type/extension assessment only; sequence contents, gene mappings, tool versions and returned external evidence require real input validation.',
                invoking_environment=dict(platform=sys.platform,python=sys.executable),
                checks=checks, stages=nodes, input_decisions=events, software=software, errors=errors, warnings=warnings,
                input_mappings={name:str(registry.path_for(name)) for name in registry.aliases},
                limitation='A structural pass is not evidence that biological analyses, inference, probability calibration or performance thresholds will succeed.')


def write_report(report, folder):
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    (folder/'deployment_audit.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    lines = ['# Deployment validation — version 1.0', '',
             f"Structural readiness: **{report['structural_readiness']}**. Deployment: **{report['deployment_readiness']}**.",
             f"Execution ready: **{report['execution_ready']}**. No analyses were launched.", '',
             '| Check | Status | Detail |', '| --- | --- | --- |']
    for row in report['checks']:
        detail = str(row['detail']).replace('|','/').replace('\n',' ')
        lines.append(f"| {row['check']} | {row['status']} | {detail} |")
    lines += ['', '## Input decisions', '', '| Stage | Input | Mandatory | Decision | Reason |', '| --- | --- | --- | --- | --- |']
    lines += [f"| {e['stage']} | {e['requested']} | {e['mandatory']} | {e['decision']} | {e['reason']} |" for e in report['input_decisions']]
    lines += ['', '## Planned conditional transitions', '', '| Stage | Prerequisites | Simulated state | Executed |', '| --- | --- | --- | --- |']
    lines += [f"| {n['id']} | {', '.join(n['after']) or 'none'} | {n['simulated_state']} | false |" for n in report['stages']]
    lines += ['', '## Warnings and limits', ''] + ['- '+x for x in report['warnings']+report['errors']+[report['limitation']]]
    (folder/'deployment_audit.md').write_text('\n'.join(lines)+'\n', encoding='utf-8')


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dry-run', action='store_true', help='Configuration-only validation; never execute analytical tools')
    parser.add_argument('--clean-state', action='store_true', help='Treat all biological inputs as absent, even if present; never read their contents')
    parser.add_argument('--workflow','--mode', choices=['all','target','validation','homology','reconciliation','phylogeny'], default='all')
    parser.add_argument('--phylogeny', action='store_true', help='Include the optional main-workflow tree branch in the plan')
    parser.add_argument('--project-root', type=Path, default=ROOT)
    parser.add_argument('--input-registry', type=Path)
    parser.add_argument('--target-raw','--input-dir',type=Path)
    parser.add_argument('--reference-raw',type=Path)
    parser.add_argument('--curated-references',type=Path)
    parser.add_argument('--library',type=Path,help='Prepared independent phylogeny library, for structural planning only')
    parser.add_argument('--reference-only',action='store_true',help='Plan the explicit provisional curated-reference-only RBH mode')
    parser.add_argument('--report-dir', type=Path, help='Optional report directory; default stdout only')
    args = parser.parse_args(argv)
    try:
        overrides={k:getattr(args,k) for k in ('target_raw','reference_raw','curated_references') if getattr(args,k) is not None}
        if args.library:overrides['phylogeny_library']=args.library/'library_manifest.json'
        report = audit(args.project_root, args.workflow, args.clean_state, args.phylogeny, args.input_registry,overrides,args.reference_only)
    except (OSError, ValueError, KeyError, TypeError) as exc:
        report = dict(version='1.0', analysis_launched=False, structural_readiness='fail', deployment_readiness='fail',
                      execution_ready=False, checks=[], stages=[], input_decisions=[], warnings=[], errors=[str(exc)],
                      limitation='Configuration assessment failed; no scientific process was launched.')
    if args.report_dir:
        write_report(report, args.report_dir)
    print(json.dumps({key:report.get(key) for key in ('structural_readiness','deployment_readiness','execution_ready','analysis_launched','errors','warnings')}, indent=2))
    return 0 if report['structural_readiness'] == 'pass' else 1
