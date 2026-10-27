# Integrinâ€“adhesome workflow â€” version 1.0

Run commands from this directory in the existing `adhesome-workflow-v1` environment.

| Entry point | Purpose |
| --- | --- |
| `scripts/run_pipeline.py` | Isolated target or validation analysis, with the parallel calibration self-test |
| `adhesome_workflow.py` | Detection handoff, external evidence return, orthology, audit and recovery |
| `scripts/report_audits.py` | Reference-label reports for completed calibration/validation audits |

Stage 2 uses the current canonical pair, `annotation/family_rules.json` and
`annotation/family_identity_rules.json`. These contain 222 family configurations
and 197 identity rules. The 159-set registry is
`annotation/family_set_registry.json`; a separately supplied reconciled registry
takes precedence if present. Coverage is not a claim of biological correctness.

The main algorithm remains **`calibrated_v1`**. The current change expands the
family and identity libraries; it does not enable stricter calling thresholds,
inspection gates or automatic-promotion restrictions. Those calibration-derived
flags have been removed from the supplied active libraries and are not recreated
from the historical registry. `src`, `Src` and `SRC` normalize to `Src`; `SCR` is
not an alias. `precision_v1` is an explicit development comparison only and is
never selected by routine runs.

Phylogenetic identity resolution is **optional and off by default**. It is
independent of optional tree inference (`--phylogeny`). A completed family audit
does not require trees, successful tree inference or resolved paralogy.

```bash
python -B adhesome_workflow.py audit \
  --annotation /absolute/path/to/stored/annotation \
  --orthology /absolute/path/to/stored/orthology \
  --out runs/identity_precision_audit --stage2-profile calibrated_v1
```

To review existing supported tree placements, add both:

```bash
--phylogenetic-identity-resolution \
--phylogenetic-identity-evidence /absolute/path/to/phylogenetic_identity_evidence.tsv
```

Missing or invalid optional evidence is reported and does not stop the family
audit. Tree concordance does not automatically transfer vertebrate paralog names.

See [the identity protocol](docs/IDENTITY_PRECISION.md) for evidence requirements,
new outputs, comparison commands and limitations, and
[the workflow guide](docs/ISOLATED_WORKFLOWS.md) for new isolated runs.

The detection HMM library now contains **214 models**, including the 38 missing
models added on 26 October 2026. All Pfam accessions referenced by the current
family rules are bundled. The original 176 model records and GA thresholds,
`calibrated_v1`, and all analytical code/configuration are unchanged. See the
[expansion receipt](maintenance/pfam_expansion_20261026.json) and
[detection README](detection/README.md). No scans, analyses or backups were created.
Existing component-reference label provenance still predates the current rules
and was not refreshed by this HMM-only update. The numerical probability model
remains stale pending a separate local refit; family confidence does not depend
on experimental annotation or that fit.

The source layout is `scripts/`, `detection/`, `validation/`, `motif_annotation/`,
`orthology_phylogeny/`, `annotation/`, `configs/`, `templates/`, `tools/` and
`tests/`. Shared imports resolve through `pipeline_paths.py`. Stored outputs in
the separate `pipeline_validation` installation remain at their original paths.

Checks without analysis:

```bash
python -B maintenance/check_workspace.py
python -B -m unittest discover -s tests -p 'test_*.py'
```
