# Optional tools

Run from the pipeline root. These tools are separate from routine analysis:

| Command | Purpose |
| --- | --- |
| `python -B tools/validate_resources.py --allow-stale-calibration` | Read-only resource/library checks |
| `python -B tools/audit_existing_workflow.py --help` | Saved-run reuse and provenance diagnostics |
| `python -B tools/evaluate_orthogroup_ablation.py --help` | Explicit fixed-coefficient replay experiments |
| `python -B tools/recalibrate_orthogroup.py --help` | Explicit calibration-only exploration |
| `python -B tools/optimize_stage2.py --help` | Calibration selection, freeze and subsequent validation experiments |
| `python -B tools/stage2_optimization_report.py --help` | Reports from completed Stage 2 experiments |

Calling `--help` never starts an experiment. Actual experiment commands require
their explicit inputs and a separate output location. Routine workflows do not
launch these development tools automatically. Their resource paths resolve to
the parent pipeline bundle, including when invoked from another directory.

The runtime modules `evaluate_validation.py` and `evaluate_family_expansion.py`
remain in the pipeline root because the maintained classifier, rule validator
and reports import their shared contracts.
