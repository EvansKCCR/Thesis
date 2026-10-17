# Validation target performance — version 1.0

Every new audit automatically evaluates a linked `validation` target using the
trusted `control_identity` labels in `02_calibration_set/calibration_set.xlsx`
and their sequences in `02_calibration_set/calibration.faa`. No experimental
annotation is required. This report measures computational reconstruction and
family annotation; it never changes predictions, confidence or evidence weights.

The default validation target is discovered from its code, original FASTA basename
or `Pipeline validation` species label. Ordinary biological targets are not treated
as labelled controls. Existing audit commands and automatic resume/continue runs
need no new flags or hand-written JSON.

## Dual gates and reported policies

Gate 1 is recruitment into the exported stage 1 hypothesis pool. Gate 2 requires
a complete family architecture gate and an unambiguous audited family assignment,
with hard/background-negative exclusion applied. The primary `dual_gate_supported`
policy accepts grades A/B at this resolution. `dual_gate_strong` separately requires
strong grade A calls. Provisional grade C calls appear in an explicitly exploratory
coverage policy; they do not count as passing the complete family gate.

The report includes stage 1 recruitment, audit decisions **before** negative
exclusion, final dual-gate decisions and strong-call decisions. This shows the
effect of the configured negative controls. Negative proteins remain exclusion
references, not positive family classes or identity calibration samples.

## Outputs after audit

For the automatic workflow, open:

```text
runs/annotation/downstream/final/validation_evaluation/performance_report.html
```

For a standalone `audit --out runs/final`, use `runs/final/validation_evaluation`.
That directory contains:

- `summary.json`: cohort coverage, primary metrics, reference/calibration overlap,
  policy definitions, input and evaluation-resource hashes.
- `metrics.tsv`: TP, FP, TN, FN, PPV, recall, F1, specificity, NPV, accuracy and
  balanced accuracy, both pooled and contrasting positives with each negative class.
- `confusion_matrices.tsv`: actual-class rows and predicted-decision columns for
  each policy and control population.
- `family_metrics.tsv`: correctness of the reported family/lineage assignment.
- `family_confusion_matrix.tsv`: admissible trusted family labels versus predictions.
- `control_predictions.tsv`: every input control, including unscored controls,
  missed positives, preparation decisions, gate outcomes, names, exclusion reasons
  and exact reference/calibration overlaps.
- `performance_report.html`: readable tables and confusion matrices without extra
  software dependencies. The nine-sheet audit workbook remains unchanged.

All metric denominators start from the **whole labelled validation input**, not
just detected candidates. A labelled positive lost during preparation or stage 1
is a false negative. Rejected negatives are true negatives. Unlabelled proteins
are not assumed negative. Duplicate hypotheses collapse to one control decision.

Labels are joined by accession and checked against the trusted sequence where
available. A known sequence mismatch or unknown label is unscored. A labelled
accession missing from the duplicate calibration FASTA remains in the primary
cohort, with `accession_label_only` recorded; a separate sequence-verified cohort
is reported for comparison. The currently staged files have 595 input sequences,
583 labels, 12 unlabelled input sequences and two labelled accessions absent from
the trusted FASTA. Thus the primary cohort has 583 controls and the sequence-verified
comparison has 581, before any subsequently detected mismatch. These coverage facts
are not performance estimates. Counts are recalculated from the actual files.

## Membership and correct family names

Membership PPV is `TP / (TP + FP)`, recall is `TP / (TP + FN)`, and
F1 is `2TP / (2TP + FP + FN)`. Undefined denominators are null in JSON,
blank in TSV and `n.a.` in HTML. A missing result is never changed to zero.
These definitions follow the standard
[precision](https://scikit-learn.org/stable/modules/generated/sklearn.metrics.precision_score.html),
[F1](https://scikit-learn.org/stable/modules/generated/sklearn.metrics.f1_score.html) and
[confusion-matrix](https://scikit-learn.org/stable/modules/generated/sklearn.metrics.confusion_matrix.html) conventions.

A positive accepted under the wrong family is still a membership-positive call,
but it is a family-assignment error. Family micro metrics count a wrong family as
one FP and one FN. Trusted family labels use the existing gold-label ontology,
including acceptable broad/specific family scopes; labels without a mapped family
are excluded only from family correctness and counted explicitly. This does not
measure exact gene, isoform or paralogue naming accuracy.

`validation.faa` overlaps the adhesome reference panel and the calibration controls.
Its scores are a **reference-panel self-test**, not independent held-out or
novel-lineage performance. Post-exclusion specificity also reflects the configured
negative-control lookup/gates. The report exposes these overlaps and the
before-exclusion results rather than describing them as independent accuracy.
Open_Orthobench F1 continues to measure orthogroups and weight orthology only;
validation-target F1/PPV are separate and are not fed back into audit scoring.

## Evaluate an already completed audit

Use a separate report directory so the sealed audit remains unchanged:

```bash
python evaluate_validation.py --audit runs/annotation/downstream/final --out runs/validation_review
```

`--target`, `--labels` and `--control-fasta` allow explicit alternatives. A failed
evaluation records its error while preserving the completed audit. A workflow
without a validation target records `validation_target_not_present` and does not
invent scores. No searches, orthology or phylogeny are rerun by evaluation.
