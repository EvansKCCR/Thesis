# Integrin–adhesome working bundle — version 1.0

Use these maintained entry points:

| Entry point | Purpose |
| --- | --- |
| `scripts/run_pipeline.py` | New target or validation analyses, each with a separate calibration self-test |
| `adhesome_workflow.py` | Staged detection, external evidence return, orthology, audit and recovery |
| `scripts/report_audits.py` | Label-based performance reports from completed standalone audits |

New runs and standalone audits use the expanded family/identity pair and the
improved `calibrated_v1` Stage 2 profile by default. Explicit
`--stage2-profile original` retains the historical comparison option. Older
recovery checkpoints keep their recorded profile; an older checkpoint without
that field retains `original`. Nothing is reclassified during installation.

Start with [the workflow guide](docs/ISOLATED_WORKFLOWS.md).

```text
integrin_adhesome_analysis/
  adhesome_workflow.py       staged workflow entry point
  pipeline_paths.py          shared bundle-root and import resolution
  scripts/                  runtime libraries, run_pipeline.py and report_audits.py
  README.md
  environment.yml / requirements.txt / VERSION
  configs/                  run settings and editable panel example
  docs/                     workflow, validation and benchmark guides
  tools/                    optional development/evaluation commands
  tests/
    workflow/               workflow, defaults, reporting and recovery tests
    audit/                  family classifier regression tests
  maintenance/              checkpoint adapters and read-only workspace checks
  schisto_adhesome/          isolated branch orchestration and reporting
  detection/                detection code, rules, reference panel and HMMs
  validation/               labelled calibration/validation inputs
  motif_annotation/
  orthology_phylogeny/
  annotation/
  templates/                workbook/export contracts
  inputs/ / returned/ / runs/   installation-specific data and saved outputs
```

Shared runtime libraries now live in `scripts/`: evidence joins/exports, family
rules, identity/ontology, exclusions, orthogroup context and terminal evaluation.
`pipeline_paths.py` exposes the bundle root and these libraries to each entry
point, independently of the caller's working directory. Resource lookups use
the renamed stage directories shown above.

Calibration fitting and the exclusion gate explicitly read `validation/calibration_set.xlsx`
and `validation/calibration.faa`; validation sequences and auxiliary workbooks
in the same directory cannot replace that pair. Label and sequence inputs,
historical model fingerprints and saved-run manifests are not rewritten by this
path update. [Path repairs](maintenance/path_repairs.json) records the active
configuration-path changes; family rule definitions and reference sequences are
unchanged. Numerical probability calibration remains disabled when its original
fingerprint is stale; no fitting is performed by path checks.

Run from this bundle root. For example, new validation analyses start with:

```bash
conda activate adhesome-workflow-v1
python -B scripts/run_pipeline.py --mode validation --select-labelled-validation \
  --out runs/validation_v1 --threads 4
```

`runs/` examples refer to outputs created under this bundle. Previously stored
outputs remain in the separate `pipeline_validation` installation. To audit
those from the reorganized source, use their absolute paths:

```bash
STORED_RUNS=/mnt/c/Users/eaasa/Desktop/Evans_PHD_Thesis/pipeline_validation/integrin_adhesome_analysis/runs
python -B adhesome_workflow.py audit \
  --annotation "$STORED_RUNS/annotation" \
  --orthology "$STORED_RUNS/annotation/downstream/orthology" \
  --benchmark "$STORED_RUNS/annotation/downstream/benchmark" \
  --out runs/target_audit_optimized \
  --stage2-profile calibrated_v1 --orthogroup-context-weight 1.0
```

Use a new output directory. Older frozen runs keep their recorded contracts;
this reorganization does not migrate them or rerun completed stages.

Optional scripts moved under `tools/`; their imports and resource roots resolve
from the script location. Their new command paths are listed in
[tools/README.md](tools/README.md). Tests and CI use the central `tests/` tree.
The move inventory is [maintenance/workspace_organization.json](maintenance/workspace_organization.json).

Check the organized bundle without analysis:

```bash
python -B maintenance/check_workspace.py
python -B tools/validate_resources.py --allow-stale-calibration
```

The first command checks paths, configuration targets, Markdown links, syntax and CLI help only. The second verifies
resource checksums without fitting a numerical identity model. Evidence-based
family confidence does not require experimental annotation or a probability fit.

For synthetic regression tests, see [tests/README.md](tests/README.md).
