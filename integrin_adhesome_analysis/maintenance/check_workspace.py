#!/usr/bin/env python3
"""Read-only bundle checks: paths, Python syntax, defaults and CLI --help only."""
import argparse
import ast
import json
from pathlib import Path
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[1]
ENTRYPOINTS=[
    'run_pipeline.py','adhesome_workflow.py','report_audits.py',
    '01_detection/adhesome_pipeline.py','01_detection/scripts/run_detection.py',
    '03_functional_motif_annotation/family_motifs.py',
    '04_orthology_phylogenetic_classification/scripts/prepare_references.py',
    '04_orthology_phylogenetic_classification/scripts/pipeline.py',
    '04_orthology_phylogenetic_classification/scripts/phylogeny.py',
    '04_orthology_phylogenetic_classification/scripts/benchmark.py',
    '04_orthology_phylogenetic_classification/scripts/adhesome_benchmark_subset.py',
    '05_conservative_adhesome_family_assignment_audit/integrated_adhesome_pipeline.py',
    'evaluate_validation.py','evaluate_family_expansion.py',
    'tools/audit_existing_workflow.py','tools/evaluate_orthogroup_ablation.py',
    'tools/recalibrate_orthogroup.py','tools/optimize_stage2.py',
    'tools/stage2_optimization_report.py','tools/validate_resources.py']
RESOURCES=[
    'configs/isolated_workflow.json','configs/panel.example.json',
    '01_detection/config/family_rules.json','01_detection/config/adhesome_pfam_manifest.csv',
    '01_detection/config/adhesome_component_crosswalk.csv',
    '01_detection/databases/pfam/adhesome_v1.0.hmm',
    '01_detection/databases/pfam/adhesome_v1.0.hmm.metadata.json',
    '01_detection/references/adhesome_v1.0/search_references.faa',
    '01_detection/references/adhesome_v1.0/metazoan_references.faa',
    '03_functional_motif_annotation/Motif_library/adhesion_core.json',
    '05_conservative_adhesome_family_assignment_audit/family_rules_expanded.json',
    '05_conservative_adhesome_family_assignment_audit/family_identity_rules_expanded.json',
    '05_conservative_adhesome_family_assignment_audit/stage2_evidence_rules.json',
    '05_conservative_adhesome_family_assignment_audit/audit_policy.json',
    'templates/output_schema.json','docs/ISOLATED_WORKFLOWS.md','README.md']


def check(help_checks=True):
    missing=[name for name in ENTRYPOINTS+RESOURCES if not (ROOT/name).is_file()]
    if missing:raise ValueError('Missing bundle paths: '+', '.join(missing))
    syntax=0
    for path in ROOT.rglob('*.py'):
        if {'runs','inputs','returned','results','__pycache__','.venv'}&set(path.relative_to(ROOT).parts):continue
        ast.parse(path.read_text(encoding='utf-8-sig'),filename=str(path));syntax+=1
    config=json.loads((ROOT/'configs/isolated_workflow.json').read_text())
    rules=json.loads((ROOT/'05_conservative_adhesome_family_assignment_audit/stage2_evidence_rules.json').read_text())
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
        resource_paths_checked=len(RESOURCES),python_files_syntax_checked=syntax,help_only_entrypoints=checked,
        caller_cwd_independent=True,analysis_launched=False,stored_outputs_modified=False)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--syntax-only',action='store_true',help='Skip subprocess --help checks')
    args=parser.parse_args()
    try:print(json.dumps(check(not args.syntax_only),indent=2));return 0
    except (OSError,ValueError,SyntaxError,subprocess.TimeoutExpired) as exc:
        print('ERROR: '+str(exc),file=sys.stderr);return 1


if __name__=='__main__':raise SystemExit(main())
