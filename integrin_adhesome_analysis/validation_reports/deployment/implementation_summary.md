# Pre-deployment audit outcome — 3 November 2026

Structural readiness: **PASS**. Deployment readiness: **WARNING**. No biological
analysis was launched. No deleted dataset, prior result, backup or duplicate
workflow directory was recreated.

The complete conditional plan contains 32 stages. The deployed library has
214 Pfam models, covering all 172 current positive/screening-rule models. Family,
identity, motif, confidence, orthology and benchmark decision thresholds were
preserved. The source-only release inventory now contains 669 deployable files;
the previous inventory had 762 unresolved historical paths.

## Implemented corrections

- Main raw reference and target roots are `input/reference_raw/` and
  `input/target_raw/`. Homology aliases resolve there, with explicit input flags
  retained and no fallback to historical local proteome outputs.
- Source species/maps and selected families determine which RBH proteomes are
  loaded. Explicit orthology panels freeze/load selected FASTAs and their own
  GFF/TSV sidecars. Unlabelled source discovery remains exact, with explicit maps
  available to restrict it.
- `unidentified_input` records stage, role/path, mandatory status, reason, decision
  and remedy. Optional missing self-tests skip; mandatory data stops the stage.
- Child branch PYTHONPATH and modular helper imports now resolve inside this
  bundle. Calibration/validation/control data remain outside the target stream,
  including known calibration IDs when only their label file remains.
- Calibration self-tests are optional when missing; mandatory negative-exclusion
  libraries are preserved, not silently disabled. Their default files currently
  coincide with the curated calibration bank. Optional unavailable model/overlap
  comparisons remain unknown, rather than being called overlap-free.
- Reconciliation requires an explicit current audit; no dated deleted audit is
  selected. Each workflow supports dataset-free structural planning.
- Git checkout rules preserve hashed scientific assets. Release text checks are
  LF-normalized; all biological/model/run byte fingerprints remain strict.

See the [changed-file map and operating instructions](../../docs/MAINTENANCE.md#clean-state-deployment-validation)
and [structured deployment report](deployment_audit.md) for implementation details.

## Automated checks

460 tests passed; 9 tests skipped because Biopython is unavailable in
the invoking environment; 0 failures and 0 errors. Detailed group
results are in [test_summary.json](test_summary.json). Four independent dry-run
entry points initialized successfully. A real target input preflight with missing
mandatory data returned the expected failure without creating an analysis run.
The isolated branch CLI also initialized with the corrected package search path.

## Remaining blockers and recommendations

- Supply actual raw proteomes and a validated target configuration before real
  detection. Supply the mandatory curated exclusion library before the audit.
  Validation, optional self-test and benchmark operations require their own
  applicable inputs; missing truth cannot produce a validated score.
- The software check describes the invoking Windows Python, not WSL. Repeat it
  inside the intended conda environment. Bioinformatics executables and
  Biopython are absent here; executable versions and OrthoFinder transitive
  backends were not executed/probed. Use the supplied environment definition and
  run the skipped tests where Biopython is installed.
- The HMM metadata's family-rule hash records its historical build. Current
  positive/screening model coverage is verified without rewriting that history.
- Historical standalone OG configurations remain source-run-specific. Review
  mappings against new OrthoFinder results or use the current reconciliation
  handoff. Existing checkpoint fingerprints are not relaxed.
- Include the byte-preserving scientific-asset differences together with
  `.gitattributes` when publishing these local changes; setting checkout rules
  alone does not repair already-normalized historical Git blobs. No commit,
  push or remote deployment was performed here.

A structural pass validates configuration and dependency planning. It is not a
claim that biological analyses, independent performance evaluation, probability
calibration or scientific acceptance thresholds will succeed.
