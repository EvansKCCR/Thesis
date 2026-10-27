#!/usr/bin/env python3
"""Read-only bundle checks: paths, Python syntax, defaults and CLI --help only."""
import argparse
import ast
import json
from pathlib import Path
import re
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[1]
ENTRYPOINTS=[
    'scripts/run_pipeline.py','adhesome_workflow.py','scripts/report_audits.py',
    'tools/evaluate_precision_identity.py','tools/prepare_phylogenetic_identity.py','annotation/family_set_review_gate.py',
    'detection/adhesome_pipeline.py','detection/scripts/run_detection.py',
    'motif_annotation/family_motifs.py',
    'orthology_phylogeny/scripts/prepare_references.py',
    'orthology_phylogeny/scripts/pipeline.py',
    'orthology_phylogeny/scripts/phylogeny.py',
    'orthology_phylogeny/scripts/benchmark.py',
    'orthology_phylogeny/scripts/adhesome_benchmark_subset.py',
    'annotation/integrated_adhesome_pipeline.py',
    'scripts/evaluate_validation.py','scripts/evaluate_family_expansion.py',
    'tools/audit_existing_workflow.py','tools/evaluate_orthogroup_ablation.py',
    'tools/recalibrate_orthogroup.py','tools/optimize_stage2.py',
    'tools/stage2_optimization_report.py','tools/validate_resources.py']
RESOURCES=[
    'pipeline_paths.py','validation/calibration_set.xlsx','validation/calibration.faa',
    'annotation/family_set_registry.json','annotation/family_set_rule_coverage_summary.tsv',
    'validation/validation_set.csv','validation/validation.faa',
    'configs/isolated_workflow.json','configs/panel.example.json',
    'detection/config/family_rules.json','detection/config/adhesome_pfam_manifest.csv',
    'detection/config/adhesome_component_crosswalk.csv',
    'detection/databases/pfam/adhesome_v1.0.hmm',
    'detection/databases/pfam/adhesome_v1.0.hmm.metadata.json',
    'detection/references/adhesome_v1.0/search_references.faa',
    'detection/references/adhesome_v1.0/metazoan_references.faa',
    'motif_annotation/Motif_library/adhesion_core.json',
    'annotation/family_rules.json',
    'annotation/family_identity_rules.json',
    'annotation/stage2_evidence_rules.json',
    'annotation/audit_policy.json',
    'templates/output_schema.json','docs/ISOLATED_WORKFLOWS.md','README.md']
CONFIG_PATHS={
    'detection/config/family_rules.json':'motif_library',
    'annotation/family_rules.json':'motif_library',
    'motif_annotation/Motif_library/context_rules.json':'domain_rule_source',
    'annotation/audit_policy.json':'identity_probability_model'}


def guide_links():
    guides={ROOT/'README.md',*(ROOT/'README').glob('*.md'),*(ROOT/'docs').glob('*.md')}
    guides.update(ROOT/name for name in ('tools/README.md','tests/README.md','detection/README.md',
        'annotation/README.md','orthology_phylogeny/README.md','orthology_phylogeny/PHYLOGENY.md'))
    count=0
    for guide in sorted(guides):
        for target in re.findall(r'\]\(([^)]+)\)',guide.read_text(encoding='utf-8')):
            target=target.split('#')[0].split(' "')[0].strip('<>')
            if not target or re.match(r'\w+://|mailto:',target):continue
            if not (guide.parent/target).exists():
                raise ValueError('Broken guide link: '+str(guide.relative_to(ROOT))+' -> '+target)
            count+=1
    return count


def check(help_checks=True):
    missing=[name for name in ENTRYPOINTS+RESOURCES if not (ROOT/name).is_file()]
    if missing:raise ValueError('Missing bundle paths: '+', '.join(missing))
    for name,key in CONFIG_PATHS.items():
        source=ROOT/name;target=json.loads(source.read_text())[key]
        if not (source.parent/target).is_file():raise ValueError('Broken config resource: '+name+' -> '+target)
    links=guide_links()
    syntax=0
    for path in ROOT.rglob('*.py'):
        if {'runs','inputs','returned','results','__pycache__','.venv'}&set(path.relative_to(ROOT).parts):continue
        ast.parse(path.read_text(encoding='utf-8-sig'),filename=str(path));syntax+=1
    config=json.loads((ROOT/'configs/isolated_workflow.json').read_text())
    rules=json.loads((ROOT/'annotation/stage2_evidence_rules.json').read_text())
    if config['stage2_profile']!=rules['default_profile'] or rules['default_profile']!='calibrated_v1':
        raise ValueError('Main workflow and Stage 2 defaults do not select calibrated_v1')
    checked=[]
    if help_checks:
        for name in ENTRYPOINTS:
            # A different working directory detects accidental reliance on cwd.
            result=subprocess.run([sys.executable,'-B',str(ROOT/name),'--help'],cwd=ROOT.parent,
                capture_output=True,text=True,timeout=30)
            if result.returncode:raise ValueError(name+' --help failed: '+result.stderr.strip())
            checked.append(name)
    return dict(status='passed',root=str(ROOT),default_stage2_profile=rules['default_profile'],
        resource_paths_checked=len(RESOURCES),config_resource_links_checked=len(CONFIG_PATHS),
        markdown_file_links_checked=links,python_files_syntax_checked=syntax,help_only_entrypoints=checked,
        caller_cwd_independent=True,analysis_launched=False,stored_outputs_modified=False)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--syntax-only',action='store_true',help='Skip subprocess --help checks')
    args=parser.parse_args()
    try:print(json.dumps(check(not args.syntax_only),indent=2));return 0
    except (OSError,ValueError,SyntaxError,subprocess.TimeoutExpired) as exc:
        print('ERROR: '+str(exc),file=sys.stderr);return 1


if __name__=='__main__':raise SystemExit(main())
