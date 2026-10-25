# Benchmark configuration and performance acceptance — version 1.0

The requested acceptance criteria are **Open_Orthobench pairwise precision at
least 0.80 and pairwise F1 at least 0.70**, assessed on unrounded values. No minimum
recall or exact-reconstruction fraction has been imposed. Here precision is the
benchmark's pairwise PPV, not the ranking metric called average precision (area
under a precision-recall curve).

The bundled policy is
`04_orthology_phylogenetic_classification/config/benchmark_policy.json`.
To use a separately specified policy, pass `--quality-policy PATH` to `benchmark`
or `--benchmark-quality-policy PATH` to workflow commands. Policies and scoring
code are hashed into the benchmark report and checked by the audit.

## What the existing result shows

Read-only inspection of the saved report found precision 0.4003824, recall
0.8345527, F1 0.5411462 and 20 exactly reconstructed RefOGs out of 70. The F1 is
consistent with the reported precision and recall. Large numbers of false pairs
in merged groups are the main precision problem; the scoring formula is not a
parameter that can repair the inferred partitions.

The run used OrthoFinder **3.1.5**, with DIAMOND, MSA inference and MCL inflation
1.2, in a panel of 22 proteomes that included the 12 benchmark proteomes. The
benchmark is therefore an **embedded reference-cohort evaluation**, not a repeat
of inference on the original 12 proteomes alone. Extra taxa and the reference-like
validation target can affect the inferred groups. Reports now record this context.

OrthoFinder 3.1.5 deliberately moved root `N0.tsv` content into
`Orthogroups/Orthogroups.tsv`: the saved report's use of that filename is correct.
The updated selector records the version and selected node and prefers the renamed
table for 3.1.5+, even if a stale `N0.tsv` is also present. See the
[official output specification](https://github.com/OrthoFinder/OrthoFinder#output-files).

The original `benchmark.py` and its Open_Orthobench equations remain unchanged.
The primary score still uses the upstream `n_refog - 1` normalization and
RefOG-specific uncertain-assignment exclusions, and is checked against the
supplied scorer. Additional per-RefOG macro averages are labelled separately;
they cannot replace the primary metrics for acceptance. Full benchmark-proteome
membership is retained when scoring each predicted group.

## Refresh the report without rerunning inference

Run the following locally from the main validation installation. Use a new output
directory; this leaves existing reports intact and creates no backup directories:

```bash
cd /mnt/c/Users/eaasa/Desktop/Evans_PHD_Thesis/pipeline_validation/integrin_adhesome_analysis
conda activate adhesome-workflow-v1
python adhesome_workflow.py benchmark \
  --orthology runs/annotation/downstream/orthology \
  --reference-raw inputs/reference_raw \
  --level root \
  --out runs/benchmark_acceptance_review
```

For unchanged predictions and inputs, the primary scores should remain unchanged;
the new acceptance policy will reject that performance. It does not reclassify
the 40% precision result as an acceptable analysis.

New reports include:

- `summary.json`: input status, prediction tree level, inference cohort context,
  original primary metrics, exact-reconstruction count/fraction, merge/fragment
  counts and the separate `quality_acceptance` decision.
- `refog_scores.tsv`: per-RefOG precision, recall, F1, false pairs, missing members,
  foreign-member counts and error categories.
- `orthogroup_error_matrix.tsv`: counts of exact groups, merge errors, splits,
  losses and mixed errors. This is a group error matrix; a true-negative count
  is not invented for uncurated protein pairs.
- `acceptance_checks.tsv`: observed primary metrics, required minima and pass/fail.

`status=validated` now requires both valid inputs and accepted performance.
`status=quality_rejected` retains all measured scores but supplies no accepted
benchmark weight. The audit can still complete, and the workflow records
`completed_with_optional_failures` with `benchmark_status=quality_rejected`.
This is not acceptance of the analysis. Legacy reports without a performance
assessment cannot supply an accepted weight until benchmarking is refreshed.

## Compare the benchmark cohort's ancestral level explicitly

The original root evaluation remains the default. An optional alternative finds
the smallest node in the labelled species tree containing **all** mapped benchmark
species. It selects that node from species membership, without examining accuracy
scores or choosing the best-performing node. In the current saved tree the 12
benchmark species occur together below node N2.

```bash
python adhesome_workflow.py benchmark \
  --orthology runs/annotation/downstream/orthology \
  --reference-raw inputs/reference_raw \
  --level cohort-lca \
  --out runs/benchmark_cohort_lca_review
```

This reuses existing HOGs without rerunning inference. Missing tree labels, missing
benchmark species or a missing selected-node table cause an explicit error; the
code never substitutes a more favorable node. The report hashes the tree and
declares the different evaluation level. Keep the original root result for
comparison, and do not describe different levels as an identical benchmark.
Use `--benchmark-level cohort-lca` in a newly resumed workflow when that level is
the predetermined evaluation scope.

## Recalibrate the underlying assignments

Improving assignments requires new inference; it cannot be accomplished by
changing precision/F1 formulas or reusing a completed partition with new parameter
labels. The workflow now exposes `--orthofinder-inflation`,
`--orthofinder-search` and `--orthofinder-method` for new inference. Unspecified
values retain the installed OrthoFinder defaults. Requested settings and the exact
command are recorded in `orthology/inference_complete.json`. Tuning options are
rejected for imported results, and recovery cannot relabel completed inference.

A **candidate calibration arm**, not a validated improvement, is stronger search
sensitivity with moderate inflation 1.5:

```bash
python adhesome_workflow.py orthology \
  --panel runs/annotation/downstream/panel \
  --reference-raw inputs/reference_raw \
  --orthofinder-search diamond_ultra_sens \
  --orthofinder-method msa \
  --orthofinder-inflation 1.5 \
  --benchmark-level root \
  --compare runs/annotation/downstream/orthology \
  --threads 4 --analysis-threads 4 \
  --out runs/orthology_calibration_i15
```

This command runs expensive inference **only when you execute it locally**.
Higher inflation can split merged clusters but may reduce recall; stronger search
can also change merges. If warranted, evaluate a separately declared 2.0 arm in
another output, keeping the same panel, ground truth, scoring level and acceptance
criteria. Benchmarking runs automatically after inference. Compare all arms,
including failures, rather than dropping difficult RefOGs or selecting only
successes. Predetermined held-out groups or independent curated groups are needed
to assess transfer after parameter selection; tuning on these 70 RefOGs is not an
independent test of the chosen parameters.

The resulting audit can reuse the completed annotation:

```bash
python adhesome_workflow.py audit \
  --annotation runs/annotation \
  --orthology runs/orthology_calibration_i15 \
  --out runs/final_calibration_i15
```

No improved performance is claimed until a local run passes the stated criteria.
The optional `adhesome_benchmark_subset.py` remains explicitly selected via
`--scope adhesome` or `--benchmark-scope adhesome`; it still requires independent
curated RefOG/family links, including fibronectin-like links where available.
Failure of the original benchmark never triggers an automatic switch to a subset.
