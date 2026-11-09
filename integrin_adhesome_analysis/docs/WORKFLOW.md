# Integrin-adhesome workflow 1.0: operating workflow

Use [MANUAL.md](MANUAL.md) for evidence/resource definitions and
[MAINTENANCE.md](MAINTENANCE.md) for repairs and experiments. Every command below
is for the user to execute locally. Reading this guide does not launch analysis.
Choose new output directories for new scientific configurations and use actual
absolute stored-run paths when reusing the separate installation.

## Clean-state structural validation

```bash
python -B maintenance/validate_deployment.py --dry-run --clean-state --phylogeny \
  --report-dir validation_reports/deployment
```

This plans every conditional workflow without opening target/reference or
calibration/validation datasets. See [deployment maintenance](MAINTENANCE.md#clean-state-deployment-validation)
for input roles, selective references, readiness limits and the changed-file map.
`--preflight-only` remains a real-data preflight; it is not a clean-state dry run.

## 1. Choose the bundle root and environment

```bash
conda activate adhesome-workflow-v1
cd /mnt/c/Users/eaasa/Desktop/Evans_PHD_Thesis/Integrin_analysis/Thesis/integrin_adhesome_analysis
```

For a genuinely new environment only, first use `conda env create -f config/environment.yml`.
Existing inference checkpoints require their original environment/tool versions;
do not replace OrthoFinder 3.1.5 during recovery. All commands in this guide run
from the root, not from individual stage/script directories.

## 2. Prepare targets and raw references

Create your own `config/targets.json` from `config/targets.example.json` and put the
configured complete protein FASTAs under `input/target_raw/`. Each entry has an
explicit gene-ID mapping; codes and species are user-defined. Do not select
validation/calibration controls for an isolated production target analysis.

Put ordinary reference FASTAs and matching GFF or protein/gene TSV files in
`input/reference_raw/`. Keep original benchmark files under
`input/reference_raw/BENCHMARKS/Input/` and `BENCHMARKS/RefOGs/` as described in
the manual. Add ordinary reference proteomes here before creating a new panel.
The resumed workflow prepares them and writes panel configuration automatically.
The optional targeted human FASTA requires its `.reference.json` scope/hash
sidecar and is not substituted for a complete proteome without declaring scope.

## 3. Recommended isolated validation or target workflow

### Validation and calibration self-test

Start validation and the calibration self-test when its datasets are available:

```bash
python -B main/run_pipeline.py --mode validation \
  --select-labelled-validation --out runs/isolated_validation --threads 4
```

The first run stops at each active branch’s external-evidence handoff. Missing optional calibration inputs skip that self-test with a recorded reason. It selects exactly the
1,737 labelled validation IDs and records the 73 unlabelled FASTA records as
exclusions. Calibration has its separate 581-record population. A total budget
of four threads is shared across the two concurrent branches. Add
`--schedule serial` for serial scheduling, or use `--threads 1`.
`--preflight-only` validates composition/leakage before running stages.

```text
runs/isolated_validation/
  run_contract.json
  preflight.json
  run_summary.json
  calibration_validation_overlap.tsv
  calibration/
    input_manifest.json
    inputs/sequence_manifest.tsv
    handoff/external_candidates.faa
    handoff/hypotheses.tsv
  validation/
    input_manifest.json
    excluded_unlabelled.tsv
    inputs/sequence_manifest.tsv
    handoff/external_candidates.faa
    handoff/hypotheses.tsv
```

Submit each branch's own `external_candidates.faa` to external services, keeping
its exact exported IDs and sequences. Do not submit the original source FASTA
and join returned results by a gene name or human accession.

### Production targets and calibration self-test

```bash
python -B main/run_pipeline.py --mode target \
  --input-dir input/target_raw --targets-config config/targets.json \
  --out runs/isolated_targets --threads 4
```

Target mode runs the independent calibration gate alongside targets; it does
not run independent validation. Control IDs/codes are rejected from production
selection. Completed production-species detection may be reused explicitly with
`--reuse-target-detection /absolute/path/to/handoff/detection` if sequence,
configuration, domain-library and resource hashes still match. This never imports
old mixed control orthogroups; changed inference composition needs new orthology.

## 4. Return external evidence and resume isolated branches

Use separate directories such as `returned/isolated_calibration/` and
`returned/isolated_validation/`, with the existing context contract:

```text
DeepLoc_adhesome_candidate_list.csv
SignalP_adhesome_candidate_list.txt
NCBI_CDD_adhesome_hitdata.txt
deeptmhmm_predicted_topologies.3line
```

Supported filename alternatives remain accepted by the importer. Preserve exact
candidate IDs and full-sequence correspondence before transferring coordinates.

```bash
python -B main/run_pipeline.py --mode validation \
  --select-labelled-validation --out runs/isolated_validation --threads 4 \
  --resume --calibration-context returned/isolated_calibration \
  --validation-context returned/isolated_validation
```

For target mode repeat its original command with `--resume`,
`--target-context returned/isolated_targets` and the appropriate
`--calibration-context`. Completed resource_library/detection/calibration/orthology checkpoints
are verified and reused. Missing evidence is an external handoff state, not a
reason to silently invent topology. `--allow-missing-context` is an explicit
alternative study policy, not the chosen returned-evidence policy.

Exit status 2 means waiting for external context, 1 means a gate failed, and 0
means the requested gates completed. Check `branch.log`, `branch_status.json`,
`stages.json` and `run_summary.json`. A failed branch does not cancel its peer.
Each branch has its own annotation, panel, orthology, benchmark, final audit and
reports. Optional integrated tree inference is a separate operation; this mode
wrapper does not automatically launch it.

During a run, the terminal shows timestamped `BRANCH` and `STAGE` start,
completion, reuse and failure messages with elapsed time. Each branch's live
output is prefixed `[calibration]`, `[validation]` or `[target]` and is also
saved to that branch's `branch.log`. Long-running commands stream their output;
OrthoFinder's own search and MSA/tree progress counters therefore remain visible
while running. Command-specific stdout/stderr logs and stage receipts remain
available for later diagnosis and resume.

New validation runs/resumes automatically create
`validation/reports/validation_review.html` and `.md` after the audit, and refresh
calibration comparisons once both branches finish. To regenerate only a saved
isolated review:

```bash
python -B -m schisto_adhesome.validation_review --run runs/isolated_validation
```

## 5. Staged workflow for existing or modular runs

### Input preflight and recruitment

```bash
python -B main/adhesome_workflow.py start --targets-config config/targets.json \
  --input-dir input/target_raw --out runs/handoff --threads 4 --check-inputs
```

The preflight does not start analysis or create the handoff. After it passes,
remove `--check-inputs` to run detection. To export from a completed detection
batch, use `--detection /absolute/path/to/detection` instead of `--input-dir`.
An incomplete detection batch cannot be sealed as a complete handoff.

The handoff includes `supported_hypotheses.tsv`, `provisional_hypotheses.tsv`,
`hypotheses.tsv`, `external_candidates.faa`, component hypotheses, reference
coverage and provenance. Alternatives for the same protein are retained,
including fibronectin-like candidates. External evidence comes back against
these exported sequence IDs.

### Resume through annotation and downstream stages

```bash
python -B main/adhesome_workflow.py resume --handoff runs/handoff \
  --context-dir returned/context --reference-raw input/reference_raw \
  --out runs/annotation --threads 4
```

This performs context/motif annotation, automatic reference preparation/panel
configuration, orthology, benchmarking and the final audit. It does not pause
for hand-written reference JSON. `--annotation-only` deliberately stops after
annotation. `--topology-file returned/predicted_topologies.3line` is available
for a topology-only return. `--motif-scope master` explicitly selects the full
347-rule library instead of the adhesion subset.

The terminal reports each checkpointed stage with a UTC timestamp, elapsed
time, and `START`, `COMPLETE`, `REUSE` or `FAILED` status. Tool output is streamed
to the terminal and copied to `logs/command_NNN.log`; OrthoFinder's native
progress counters and time estimates are visible during inference rather than
appearing only after the command ends. The same progress messages are retained
in the parent run/branch log.

If annotation is already complete, use its actual sealed directory:

```bash
python -B main/adhesome_workflow.py continue --annotation runs/annotation \
  --reference-raw input/reference_raw --out runs/downstream_new --threads 4
```

`continue` creates a new downstream analysis. It is useful after reference
expansion; it is not the command to resume a partially inferred old partition.
For a failed existing downstream checkpoint use recovery as below.

```text
runs/annotation/downstream/
  references/panel.auto.json
  panel/
  orthology/
  benchmark/summary.json
  final/Stage_2_output.xlsx
  final/joined_evidence.tsv
  recovery.json
  workflow_manifest.json
```

### Expand a reference panel

Add appropriate FASTAs/GFF/TSV sidecars to the raw-reference directory, then use
a new `continue` output. The `calibrate` and `orthology` subcommands also expose
modular panel/inference controls; inspect their `--help` before choosing an
explicit manual route. `--compare /absolute/path/to/previous/orthology` records
panel comparisons where supported. Existing panels/results remain fixed.
The hand-written panel format in `config/panel.example.json` is an optional
advanced route, not a required step in the automatic workflow.

## 6. Audit replay without repeating searches or inference

Use current `calibrated_v1` and the canonical paired rules with already completed,
matching annotation and orthology:

```bash
python -B main/adhesome_workflow.py audit \
  --annotation /absolute/path/to/stored/annotation \
  --orthology /absolute/path/to/stored/orthology \
  --benchmark /absolute/path/to/stored/benchmark \
  --out runs/audit_current \
  --family-coverage expanded --stage2-profile calibrated_v1 \
  --orthogroup-context-weight 1.0
```

Omit `--benchmark` if no compatible report is available. The fixed coefficient
preserves the isolated routine baseline; omit the fixed override only for an
explicit comparison with the core's benchmark-dependent policy. The command
does not rerun detection, annotation, OrthoFinder or phylogeny. New HMM models
cannot appear in cached domain observations merely through an audit replay.
The annotation folder must contain its actual `workflow_manifest.json`; a
documentation example such as `runs/resumed` is not created by recovery.

### Calibration and validation reports from standalone audits

Both audit folders must be completed/sealed. Use the original isolated run and
its reference labels. For the saved installation:

```bash
INSTALL=/mnt/c/Users/eaasa/Desktop/Evans_PHD_Thesis/pipeline_validation/integrin_adhesome_analysis

python -B main/report_audits.py \
  --run "$INSTALL/runs/isolated_validation" \
  --calibration-audit runs/audits_current/calibration \
  --validation-audit runs/audits_current/validation \
  --calibration-labels "$INSTALL/02_calibration_set/calibration_set.xlsx" \
  --validation-labels "$INSTALL/02_calibration_set/validation_set.csv" \
  --out runs/reports_current
```

Replace the audit paths with the actual outputs. Labels default to the frozen
run contract if their flags are omitted. A repackaged calibration workbook may
be allowed **only** with `--allow-equivalent-calibration-label-bytes`, which
still requires every evaluation field, ID, sequence and cohort selection to
match the frozen manifest. It never relaxes validation-label hash checks.

The separate output provides `report.html`,
`calibration/calibration_review.html`, `validation/validation_review.html`,
metrics, family metrics, errors, predictions and reporting provenance. It does
not refit rules/calibration or rewrite either audit/original dashboard.

For real targets, report the completed target audit separately:

```bash
python -B -m schisto_adhesome.reports --audit runs/audits_current/target \
  --out runs/target_report_current --exclude-species validation
```

The exclusion is an exact stored species code and affects reporting only. It
handles historical mixed-target runs without modifying their sealed audits.

## 7. Benchmark existing inference or select a subset

The inference panel must already contain the matching original benchmark
proteomes. Adding them afterward needs new inference.

```bash
python -B main/adhesome_workflow.py benchmark \
  --orthology /absolute/path/to/stored/orthology \
  --reference-raw input/reference_raw --level root --out runs/benchmark_review
```

Explicit alternatives:

```text
--level cohort-lca
--scope adhesome --selection input/reference_raw/adhesome_refog_selection.tsv
--protocol custom
--quality-policy /absolute/path/to/quality_policy.json
```

Select only options appropriate to your declared evaluation. The workflow forms
are `--benchmark-level`, `--benchmark-scope`, `--benchmark-selection`,
`--benchmark-protocol`, `--benchmark-quality-policy`. Subset curation uses the
header-only template in `resource_library/templates/adhesome_refog_selection.tsv`. Keep full
predicted groups and preserve root/cohort-level distinctions. A failed full
benchmark does not trigger subset fallback. Use the resulting sealed report
explicitly in a new audit. Benchmark input status and precision/F1 acceptance
remain separate in `summary.json` and `acceptance_checks.tsv`.

## 8. Optional integrated phylogeny and identity review

Add `--phylogeny` to staged resume/continue to request inference **after** audit.
For family trees also add `--phylogeny-grouping family`. Grouping alone does not
launch trees. Missing maintenance, bad alignments or tree failures preserve the audit.
Integrated IQ-TREE threads default independently to two; use
`--iqtree-threads 2` explicitly when reproducing that choice.

To prepare audited A/B/C family trees from existing evidence without repeating
orthology:

```bash
python -B resource_library/orthology_phylogeny/scripts/phylogeny.py prepare \
  --hogs /absolute/path/to/orthology/candidates/candidate_HOGs.tsv \
  --orthologues /absolute/path/to/orthology/candidates/direct_orthologues.tsv \
  --duplications /absolute/path/to/orthology/candidates/candidate_duplications.tsv \
  --proteomes /absolute/path/to/panel/proteomes \
  --species-list /absolute/path/to/panel/species.json \
  --target-species /absolute/path/to/panel/targets.json \
  --audit-workbook runs/audit_current/Stage_2_output.xlsx \
  --candidate-report /absolute/path/to/orthology/candidates/candidate_report.tsv \
  --output runs/family_phylogeny
```

If duplication evidence is genuinely unavailable, explicitly use
`--allow-missing-duplications` instead of its path. Repeat `--family` for selected
audit family names. Inspect `groups.tsv`, `family_candidate_coverage.tsv`,
`family_orthogroups.tsv`, per-family `tips.tsv` and the preparation manifest
before inference. Target members without source orthogroups remain traceable.

```bash
python -B resource_library/orthology_phylogeny/scripts/phylogeny.py run \
  --project runs/family_phylogeny --threads 2 --iqtree-threads 2
```

For interruption recovery repeat with `--resume`. To retry reviewed failed
groups in place, add `--resume --retry-failed`; `--hog GROUP_ID` limits a group.
Do not confuse this integrated CLI or its filenames with the sibling standalone
family-phylogeny runner described in that bundle's own docs.

Existing-tree identity review is a separate optional operation:

```bash
python -B maintenance/prepare_phylogenetic_identity.py \
  --project /absolute/path/to/stored/phylogeny \
  --reference-labels /absolute/path/to/curated_reference_tip_labels.tsv \
  --out runs/identity_tree_review
```

Without reference labels this creates an inventory/template, not inferred
identities. Labels require `sequence_id`, `family`, `species`, `source` and exact
tip IDs. Add both `--phylogenetic-identity-resolution` and
`--phylogenetic-identity-evidence runs/identity_tree_review/phylogenetic_identity_evidence.tsv`
to an audit only when the evidence is supplied. Isolated workflows use the
corresponding branch-specific calibration/validation/target evidence flags.
See the manual for identity scope and rooting/reconciliation requirements.

## 9. Resume an existing failed run

For unfinished resource_library/detection/handoff, `start --resume` verifies and preserves
completed target evidence. A known invalid transcript target can be explicitly
excluded with `--exclude-target Sguin`; do not silently omit successful targets.
This is different from the `resume --handoff` external-return command.

```bash
python -B main/adhesome_workflow.py start --resume --exclude-target Sguin \
  --input-dir input/target_raw --targets-config config/targets.json \
  --out runs/handoff --threads 4
```

For a saved automatic downstream recovery checkpoint:

```bash
python -B main/adhesome_workflow.py recover \
  --run /absolute/path/to/existing/downstream --threads 2 --iqtree-threads 2
```

Use only actual saved paths and avoid a second writer to the same run.
`--adopt-existing --annotation ACTUAL_ANNOTATION_PATH` is an explicit legacy
adoption route requiring verified panel/reference provenance. Adding
`--phylogeny`, `--retry-benchmark`, or changing dependent-stage inputs can replace
stage outputs and invoke archival behavior in the legacy core. These are not
blanket no-backup operations. For no-archive benchmark/audit refresh, use new
standalone output directories instead. The isolated wrapper refuses legacy
recovery paths that would silently archive partial evidence.

See [MAINTENANCE.md](MAINTENANCE.md) for input/reference repairs, saved-version
constraints, probability refits and resource checks. No guide command changes
the scientific policies merely by being documented here.

## Optional independent RBH / main-audit reconciliation

The [reconciliation guide](../../reconciliation/docs/WORKFLOW.md) documents a
separate candidate handoff and `library_reconciled` family-tree library. The
three scientific workflows remain independently runnable. Existing audit
classifications, family libraries and analysis checkpoints are unchanged. Use
the new `rbh_audit_reconciled_v1` run name for this library; its source annotations
and short tree identifiers are traceable in the accompanying TSV files.

## Optional post hoc refinement of unresolved Stage 2 assignments

Post hoc refinement needs a saved Stage 2 workbook, but it does not require
orthology. You can create a preliminary Stage 2 audit immediately after motif
annotation, review its Pfam/CDD evidence, and later repeat the refinement against
the orthology-enriched audit. Keep each pass in a separate output directory.

First produce motif annotation without automatically entering the downstream
reference/orthology workflow:

```bash
python -B main/adhesome_workflow.py resume \
  --handoff runs/handoff/validation \
  --context-dir returned/validation \
  --out runs/annotation/validation \
  --annotation-only
```

Run a preliminary audit from annotation alone. Omitting `--orthology` marks
orthology and benchmark fields as not run; the family audit uses the available
Stage 1, external-context, Pfam/CDD and motif evidence:

```bash
python -B main/adhesome_workflow.py audit \
  --annotation runs/annotation/validation \
  --out runs/preliminary_audit/validation
```

The same command works on a production-target annotation. Validation branches
also generate the control performance evaluation when the validation target is
present. This preliminary audit is a separate snapshot; it does not occupy or
modify the later `annotation/downstream/final` result.

Run the independent post hoc refinement and report before orthology:

```bash
python -B main/adhesome_workflow.py refine-unresolved \
  --input runs/preliminary_audit/validation/Stage_2_output.xlsx \
  --out runs/posthoc_domain_review_v1/validation

python -B main/adhesome_workflow.py report-refinement \
  --refinement runs/posthoc_domain_review_v1/validation \
  --labels validation/validation_set.csv \
  --sequence-manifest runs/annotation/validation/inputs/sequence_manifest.tsv \
  --kind validation \
  --out runs/posthoc_performance_v1/validation
```

`--check-inputs` validates the input without writing outputs; `--out` is optional
in that mode. The audit TSV/CSV export is also accepted. Use a new output folder
outside the sealed source audit. This invokes the sibling
[`schisto_pfam_cdd_reclassify`](../../schisto_pfam_cdd_reclassify/README_local_Pfam_CDD_reclassification.md)
adapter and does not restart detection, orthology, audit, calibration or phylogeny.
The refinement can be run again later against the second audit; choose a new
`--out` directory for each pass.

The separate `Stage_2_posthoc_refinement.tsv` preserves every source column and
adds `posthoc_*` proposal fields. Supported and negative-excluded assignments are
preserved; unresolved, ambiguous and provisional assignments enter review.
`Schisto_unassigned_Pfam_CDD_classification.tsv` contains only that review subset.
Review tiers A–D are heuristic domain evidence tiers, distinct from audit grades
and calibrated probabilities. They do not change validation performance metrics
or promote original assignments automatically. The source workbook is unchanged.

To evaluate saved calibration or validation refinement outputs against their
reference labels, without repeating classification:

```bash
python -B main/adhesome_workflow.py report-refinement \
  --refinement /absolute/path/to/completed/posthoc/validation \
  --labels /absolute/path/to/validation_set.csv \
  --sequence-manifest /absolute/path/to/validation/inputs/sequence_manifest.tsv \
  --kind validation \
  --out runs/posthoc_performance_v1/validation
```

For calibration select `--kind calibration`, `calibration_set.xlsx` and the
calibration branch's own sequence manifest. Original supported calls, tier A
proposals and tier A+B proposals receive separate performance tables and confusion
counts. Full-cohort metrics retain missing audit positives as false negatives;
proposal-only metrics also report the eligible unresolved subset. Family-label
coverage gaps are explicit. No combined/promoted audit result is generated.

After orthology completes, refine the new enriched audit and rerun the report in
a new directory. Supply `--orthogroups` with that branch's saved `Orthogroups.tsv` and
`--panel-proteomes` with its `downstream/panel/proteomes` directory to add
positive-control reconstruction counts by species and match type. These inputs
must belong to the same saved branch. Recovery requires a full-length identical
sequence in the control's orthogroup from the same-species reference proteome.
ID agreement does not affect recovery. Available identical sequences missing
from the group and other same-species members are reported separately. Species
without a mapped reference are not assessable. See the sibling guide for
mapping overrides.

To run orthology and the orthology-enriched audit after the preliminary pass,
continue from the same annotation in a new downstream directory:

```bash
python -B main/adhesome_workflow.py continue \
  --annotation runs/annotation/validation \
  --out runs/annotation/validation/downstream \
  --threads 4 --analysis-threads 4
```

`continue` prepares references, runs orthology and the benchmark, then writes
the enriched audit at `runs/annotation/validation/downstream/final/Stage_2_output.xlsx`.
If a downstream stage fails, resume it with
`python -B main/adhesome_workflow.py recover --run runs/annotation/validation/downstream`.
Once the enriched audit and reports are
ready, run `recover --run runs/annotation/validation/downstream --phylogeny` to
start or resume the optional phylogeny without repeating completed orthology or
audit stages.
