# Parallel workflow gates — version 1.0

Run commands from the directory containing `run_pipeline.py`, in the existing
`adhesome-workflow-v1` Linux/WSL environment. The expanded family/identity pair,
motif library, exclusion gates and Atlas workbook remain the shared analytical
core. `adhesome_workflow.py` remains available for legacy staged runs.

## Start independent validation and its mandatory self-test

```bash
conda activate adhesome-workflow-v1
cd /mnt/c/Users/eaasa/Desktop/Evans_PHD_Thesis/pipeline_validation/integrin_adhesome_analysis
python run_pipeline.py --mode validation --select-labelled-validation \
  --out runs/isolated_validation --threads 4
```

Calibration contains exactly 581 controls: 234 positive, 119 hard negative and
228 background negative. Validation selects the 1,737 labelled IDs: 902 positive
orthologues and 835 negatives. The supplied FASTA has 1,810 records. The explicit
selection flag records the other 73 IDs in `validation/excluded_unlabelled.tsv`;
it never edits the scientific input. Count, identifier, class, species/taxid,
labelled length, residue and conflicting-label checks run before detection.
No extra short-protein filter has been introduced; the existing sequence QC is
recorded in the configuration and failed recruitment remains in the denominator.

Two processes run concurrently, with two threads per branch under the example
budget. `--schedule serial` or `--threads 1` schedules the independent gates
sequentially. A failure in either gate does not cancel the other gate.

The default first run stops at **both external-evidence handoffs**:

```text
runs/isolated_validation/
  run_contract.json                 frozen code/configuration/dataset/reference hashes
  preflight.json                    composition and leakage checks
  calibration_validation_overlap.tsv
  run_summary.json                 separate gate states; not completed metrics
  calibration/
    input_manifest.json
    inputs/sequence_manifest.tsv   original ID → branch/species protein ID
    handoff/external_candidates.faa
    handoff/hypotheses.tsv          supported/provisional hypotheses, including fibronectin-like
  validation/
    input_manifest.json
    excluded_unlabelled.tsv
    inputs/sequence_manifest.tsv
    handoff/external_candidates.faa
    handoff/hypotheses.tsv
```

Use **each branch's own candidate FASTA** for external tools. Preserve its exported
IDs. Do not submit the original FASTA and try to join on names or human accessions.
The sequence manifests retain original IDs and trusted labels for reporting;
analysis FASTA headers contain only opaque protein units and explicit unit IDs.
Controls remain protein units, not claims of complete gene-level proteomes.

## Return external evidence and resume

Place each branch's returned files in its own directory, for example
`returned/isolated_calibration/` and `returned/isolated_validation/`. Each directory
uses the existing context contract:

- `DeepLoc_adhesome_candidate_list.csv`
- `SignalP_adhesome_candidate_list.txt`
- `NCBI_CDD_adhesome_hitdata.txt`
- `deeptmhmm_predicted_topologies.3line`

The existing alternative filenames remain accepted by the core importer. Keep
sequence IDs and sequences unchanged. Then repeat the original command with
`--resume` and both context paths:

```bash
python run_pipeline.py --mode validation --select-labelled-validation \
  --out runs/isolated_validation --threads 4 --resume \
  --calibration-context returned/isolated_calibration \
  --validation-context returned/isolated_validation
```

Completed detection is reused. Context import/motifs, automatic raw-reference
preparation, isolated OrthoFinder inference, the original benchmark and the family
audit continue through the existing core. Each branch has its own `downstream/`
panel, inference, benchmark, final workbook and `reports/`. No query branch can
modify another branch's orthogroups. Phylogeny is not run by this mode wrapper;
it remains an optional separate operation through the existing phylogeny CLI.

`--allow-missing-context` is an explicit alternative for a study that accepts
missing external evidence. It is **not** the chosen policy for the current run.
Exit status 2 means external evidence is awaited; 1 means a gate failed; 0 means
both selected gates completed. See `branch.log`, `branch_status.json`, `stages.json`
and `run_summary.json`, not just the presence of an output directory.

## Production target mode

Use a target configuration containing only the intended production species.
Keep the existing general `proteomes` schema: `code`, `species`, `assembly`,
`path`, `gene_regex`; absolute FASTA paths are recommended. This is not restricted
to three named proteomes. The configured QC must match the frozen mode config.
The input folder may contain other files, but control files must never be selected.

```bash
python run_pipeline.py --mode target --input-dir inputs/targets \
  --targets-config targets.production.json --out runs/isolated_targets --threads 4
```

Resume using `--resume --target-context returned/isolated_targets
--calibration-context returned/isolated_calibration` on the same command.
Target mode does not execute independent validation. It checks validation IDs for
production leakage only. Production codes `validation`, `calibration`, `Cal_*`
and `Val_*`, or selected evaluation-control IDs, are rejected.

`--reuse-target-detection /absolute/path/to/runs/handoff/detection` imports only
the explicitly selected completed production species after sequence, input,
domain-library and classification-resource verification. An old mixed detection
batch can therefore supply its seven valid target species without importing its
control species. Reuse refuses stale resources. Existing raw returned target
context may be supplied again when its sequence IDs still match. Domain/motif
features remain reusable; **old mixed orthogroups must be recomputed**.

## Frozen configuration and calibration exploration

`--config configs/isolated_workflow.json` selects the default expanded configuration.
Its fixed context coefficient is **1.0**, preserving neutral context in the latest
quality-rejected-benchmark baseline. The original CLI's benchmark-dependent
baseline remains available for controlled comparison. An isolated routine run
requires an explicit fixed coefficient; benchmark scores are diagnostics and
cannot silently refit it. This coefficient is not a protein probability.

Separate explicit exploration, **after isolated calibration has completed**:

```bash
python recalibrate_orthogroup.py --run runs/isolated_validation \
  --out runs/calibration_weight_exploration
```

This reuses sealed calibration joined evidence and compares the frozen baseline,
the original benchmark policy, 0.541 and zero. Repeat `--additional-weight 0.25`
for justified extra arms. It checks that independent evidence is unchanged. It
never loads validation labels for fitting, automatically chooses weights, writes
an approved configuration, or interprets 0.541 as a probability. Explicitly freeze
a justified new configuration before final validation; do not select it by its
validation score. A changed analytical contract requires a distinct scientific run,
not relabelling old outputs.

## Reports and limits

`reports/metrics.tsv` preserves recruitment, supported and provisional tiers,
confusion counts, PPV/precision, recall, F1, specificity, FPR and balanced accuracy.
Separate strata cover species, negative classes, length bins, exact calibration,
reference and saved-training overlap. Family micro/macro scores have explicit
unmapped-label and undefined-denominator accounting. Negative labels remain
exclusion gates, not family classes. False positives, false negatives and wrong
family calls remain traceable in `errors.tsv` and `predictions.tsv`.

Ten validation positives exactly match calibration sequences. Shared orthology
reference overlap is also checked against the **whole raw reference panel**, not
just curated detection seeds. All records remain; overlap-free strata support
the appropriately qualified interpretation. Reference homologues and rule-development
history mean that an absent exact match alone cannot prove complete independence.
No experimental-evidence requirement has been added to family confidence.

`--compare-calibration /previous/calibration/reports` writes shared-control
assignment drift; without a previous report drift is explicitly unavailable.
Target reports retain the standard `Stage_2_output.xlsx` and compact assignment
TSV for Atlas ingestion, check candidate uniqueness, and summarize species/family,
confidence and unresolved calls.

Run the read-only historical audit separately:

```bash
python audit_existing_workflow.py --runs runs --out reports/existing_result_audit
```

It classifies feature reuse, reports lengths, historical invalid inputs and mixed
orthogroup membership. Add `--corrected-downstream runs/isolated_validation/validation/downstream`
only after corrected inference. Orthogroup pairwise F1 and family-classification F1
are different measures; the reports do not conflate them or assert that short
proteins caused the old F1 of 0.541. Full predicted benchmark groups remain intact.

No archives, `.bak` files or project copies are created. Completed sealed stages
are reused and hash-verified. Partial tool outputs are retained in place; if the
legacy recovery would archive them, this wrapper stops with an explicit diagnostic
instead. Tool-specific partial-stage recovery still needs review; the wrapper does
not silently delete evidence or pretend a partial inference is complete.
