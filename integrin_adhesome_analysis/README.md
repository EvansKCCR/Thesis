# Integrin–adhesome working bundle — version 1.0

Use these maintained entry points:

| Entry point | Purpose |
| --- | --- |
| `run_pipeline.py` | New target or validation analyses, each with a separate calibration self-test |
| `adhesome_workflow.py` | Staged detection, external evidence return, orthology, audit and recovery |
| `report_audits.py` | Label-based performance reports from completed standalone audits |

New runs and standalone audits use the expanded family/identity pair and the
improved `calibrated_v1` Stage 2 profile by default. Explicit
`--stage2-profile original` retains the historical comparison option. Older
recovery checkpoints keep their recorded profile; an older checkpoint without
that field retains `original`. Nothing is reclassified during installation.

Start with [the workflow guide](docs/ISOLATED_WORKFLOWS.md).

```text
integrin_adhesome_analysis/
  run_pipeline.py / adhesome_workflow.py / report_audits.py
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
  01_detection/             detection code, rules, reference panel and HMMs
  02_calibration_set/       labelled calibration/validation inputs
  03_functional_motif_annotation/
  04_orthology_phylogenetic_classification/
  05_conservative_adhesome_family_assignment_audit/
  templates/                workbook/export contracts
  inputs/ / returned/ / runs/   installation-specific data and saved outputs
```

The root Python modules beside the three entry points are shared runtime
libraries: evidence joins/exports, family rules, identity/ontology, exclusions,
orthogroup context and terminal evaluation. Stage scripts import them directly.
They remain at stable paths because provenance and calibrated-model fingerprints
record those locations. They are active dependencies, not abandoned scripts.

Optional scripts moved under `tools/`; their imports and resource roots resolve
from the script location. Their new command paths are listed in
[tools/README.md](tools/README.md). Tests and CI use the central `tests/` tree.
The move inventory is [maintenance/workspace_organization.json](maintenance/workspace_organization.json).

Check the organized bundle without analysis:

```bash
python -B maintenance/check_workspace.py
python -B tools/validate_resources.py --allow-stale-calibration
```

The first command checks paths, syntax and CLI help only. The second verifies
resource checksums without fitting a numerical identity model. Evidence-based
family confidence does not require experimental annotation or a probability fit.

For synthetic regression tests, see [tests/README.md](tests/README.md).
