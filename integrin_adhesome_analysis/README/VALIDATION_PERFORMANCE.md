# Validation target performance — version 1.0

Every new audit automatically evaluates a linked `validation` target using the
trusted `control_identity` labels in `validation/calibration_set.xlsx`
and their sequences in `validation/calibration.faa`. No experimental
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
python scripts/evaluate_validation.py --audit runs/annotation/downstream/final --out runs/validation_review
```

`--target`, `--labels` and `--control-fasta` allow explicit alternatives. A failed
evaluation records its error while preserving the completed audit. A workflow
without a validation target records `validation_target_not_present` and does not
invent scores. No searches, orthology or phylogeny are rerun by evaluation.

## Validation revision — 18 October 2026, version 1.0

The evaluator now emits the following additional reports after a future audit:

| Report | What it separates |
|---|---|
| `exclusion_ablation.tsv` | Complete family decisions before any exclusions, with only exact-sequence exclusions removed, and after all exclusions; pooled, hard-negative and background-negative populations |
| `exclusion_impact.tsv` | Counts excluded by each recorded gate and how many had passed the family gate |
| `confidence_tiers.tsv` | Mutually exclusive supported, provisional-only, ambiguous-family, recruited-unresolved, negative-excluded and not-recruited groups; recruitment metrics remain separate |
| `family_mapping_review.tsv` | Every scored positive, its trusted name/genes, expected families, label source and architecture-rule coverage |
| `family_label_coverage.tsv` | Labels covered by configured architecture gates, labels outside those gates, and unresolved mappings |
| `missed_positive_review.tsv` | Unsupported positives and incorrect accepted families, with preparation/recruitment, architecture, CDD, order and ambiguity diagnostics |
| `family_misclassifications.tsv` | Incorrect accepted family names without conflating them with membership false positives |
| `per_expected_family.tsv` | Recovery and correct/wrong supported calls for each admissible expected-family group |
| `sequence_group_conflicts.tsv` | Exact-sequence groups whose trusted labels or saved decisions disagree |

The exact-exclusion ablation retains diagnostic negative gates, including when a
sequence triggered both gates. New exclusions record every applicable gate in
`exclusion_gates`. An older exact-match record that lists only its first gate
cannot determine the exact-only counterfactual. It is explicitly unevaluable for
that policy, with its denominator reported; the before-all-exclusions comparison
remains available. No control is silently treated as an accepted or rejected
counterfactual when the provenance is missing.

`metrics.tsv` also includes calibration-overlapping/non-overlapping controls,
forward/reverse-reference overlap strata, calibration-by-forward-overlap groups,
controls without any exact reference/calibration overlap, and a comparison that
counts consistent exact-sequence groups once. Empty strata have undefined rates,
not perfect performance. Conflicting exact-sequence groups are excluded only from
the sequence-deduplicated comparison and remain in the primary control cohort.

`scripts/label_ontology.py` expands trusted-name family labels beyond the small original
architecture catalogue. Explicit `trusted_family` labels take precedence. A name
such as "microtubule-actin cross-linking factor" is not mapped to actin simply
because it contains that word. These mappings are evaluation labels and never
classifier inputs. Mapping a coronin, a transporter or an enzyme family does not
invent a corresponding architecture gate. The original workbook is unchanged;
unresolved or questionable trusted names remain visible for review.

Stage 2 now matches diagnostic CDD patterns against individual observed models,
so an anchored model pattern works even when that model is not the first listed
hit. Additional gated refinements cover receptor PTP alpha/epsilon, SHP, MACF1,
cortactin/HS1, LASP, Ena/VASP, ERM, BCAR/CAS, FHL, CSRP and coronin. They require
the configured domain types/copies, order where specified, an observed diagnostic
CDD and at least one compatible Pfam anchor. These refinements apply only within
their declared stage 1 parent/screening scopes and retain the same protein pool.
Generic PTP domains do not pass the PTP-PEST complete gate without PEST-family CDD.
Alpha-actinin gates require two CH domains, bounded spectrin copy/length ranges
and diagnostic CDD; long spectraplakin architectures cannot pass that complete
gate. Shared LIM repeats also require a diagnostic family label for a supported
paxillin, zyxin or PDLIM call. Experimental metadata is never a gate or penalty.

These changes were informed by the supplied control annotations and observed
errors. Their next control scores will still be development/self-test results.
Overlap stratification and saved-output ablation do not constitute independently
retrained holdouts. For independent transfer assessment, freeze the rules, hold
out whole family/lineage and exact-sequence groups from reference and calibration
resources, and rerun those groups locally. No such inference is performed during
code installation.

The old numerical probability model is deliberately left unchanged and becomes
stale after the scoring/ontology change. Strong/moderate/provisional evidence
confidence and family/protein names remain available. Numerical probabilities
stay unavailable until an explicitly requested local refit; calibration overlap
is not an acceptance requirement. Resource-only checks may use:

```bash
python tools/validate_resources.py --allow-stale-calibration
```

To test the revised audit later using already completed annotation and orthology,
write to a new audit output (the command is not executed by installation):

```bash
python adhesome_workflow.py audit --annotation runs/annotation --orthology runs/annotation/downstream/orthology --out runs/final_validation_revision
```

This reuses those upstream results and writes new audit/evaluation results. It
does not rerun detection, annotation, OrthoFinder or phylogeny. Existing sealed
audit reports retain their original scores.

## Stage 2 routing and annotation resolution — version 1.0

The reviewed coronin, cortactin/HS1, Ena/VASP, BCAR/CAS, ERM, spectraplakin and
PTP refinements can now dispatch from an unrelated preliminary hypothesis when
an observed diagnostic CDD model and a compatible coordinate-bearing Pfam anchor
are present. Each rule must still pass its full domain-copy, ordering, coordinate
and applicable length gates. Dispatch is explicitly enabled on reviewed rules;
it does not search the full catalogue to replace an arbitrary hypothesis. Target
names, accessions, expected labels and calibration membership never dispatch a
rule. The exported Stage 1 protein pool remains unchanged.

Observed `PF06297` is recognized as PET with its Pfam provenance. Testin now has a
separate PET + three ordered LIM rule; Prickle retains a separate five-LIM rule.
CDD/Pfam alternative models on the same footprint count once. A CDD summary such
as `LIM [3 regions]` supplies no coordinates and cannot manufacture three domains.
Observed bracketed CDD model/coordinate strings are accepted, and models outside
the full sequence are rejected.

New final fields are `family_support_level`, `annotation_resolution`,
`specific_identity_status`, `unresolved_specific_families` and a compact
`identity_gate_trace`. A supported broad-family call requires a complete explicit
broad core; it does not promote arbitrary grade C calls. Supported broad actin,
Src/Csk and PET-LIM cores can preserve useful annotation without inventing an Arp
subtype, a specific LIM scaffold, a gene or a paralogue. Full distinct-family
support is reported separately. Orthology-unresolved complete architecture
remains grade B; experimental evidence is never a gate or penalty.

After a future audit, these files explain the actual evaluated evidence:

- `stage2_rule_gate_trace.tsv`: applicable per-rule dispatch, diagnostic CDD,
  Pfam-anchor, domain-copy, order, coordinate and length outcomes.
- `stage2_candidate_gate_trace.tsv`: one resolved protein disposition, including
  global negative exclusion and the compact gate explanation.
- `validation_evaluation/stage2_gate_diagnostics.tsv`: control-level loss and
  limitation codes, including `NO_STAGE1_HYPOTHESIS`,
  `NO_IDENTITY_RULE_REACHED`, `MISSING_REQUIRED_DOMAIN`,
  `DOMAIN_ORDER_CONFLICT`, `ORTHOLOGY_UNRESOLVED`, `COMPETING_FAMILY` and
  `NEGATIVE_EXCLUSION`. Unrecognized aliases are reported only when established
  by a saved trace; missing diagnostic evidence is not called an alias failure.
- `family_macro_metrics.tsv` and `family_per_group_metrics.tsv`: equal-family
  comparisons alongside pooled micro metrics, with mapped/unmapped denominators.

Detailed rule traces stay outside Excel cells to avoid its character limit.
Existing legacy outputs without a trace remain explicitly unavailable for that
diagnostic; specific resolution is not guessed from a confidence adjective.

The trusted identity ontology no longer treats recruitment-catalogue family
labels as gold identity. GAB1 is distinct from GRB2, ICAP1 from integrin beta,
ARPC2 from ACTR2/ARP2, PPM1F from ILK, migfilin from filamin, and MACF1 from actin.
Testin and Prickle have distinct expected identities. Explicit trusted-family
labels take precedence; inconsistent own-gene/name labels require review. Valid
broad ancestors are reported separately and receive credit only in the broad
family policy. Mapped-family correctness is not labelled whole-cohort family
sensitivity. Candidate retention, broad support, distinct-family support,
provisional annotation and ambiguity are separate policies.

The current hard/background-control membership exclusion remains unchanged.
`negative_for` entries in the uncurated 92-row proposed benchmark are not imported
as global membership negatives. That template's missing accession/sequence/source
links must be curated before it is an operational benchmark. It supplies no
calibration weight or independent performance claim.

For a controlled local comparison, first reevaluate the unchanged old audit with
the corrected gold ontology, then test the revised audit using the same completed
annotation and orthology. These commands are provided for the user's local run;
installation does not execute them:

```bash
python scripts/evaluate_validation.py --audit runs/annotation/downstream/final --out runs/baseline_corrected_labels
python adhesome_workflow.py audit --annotation runs/annotation --orthology runs/annotation/downstream/orthology --out runs/stage2_recall_review
```

Compare the same-policy metrics and overlap strata in those report directories.
The original audit and its original scores remain untouched. New resolution
fields or traces may be unavailable for the baseline; they are not retroactively
inferred. Corrected-label scores should not be compared directly to an older
72-mapped-control headline without using the same label scope and denominator.
The annotation workbook's domain examples guide rule tests; trusted protein names
or gold family labels are never injected into the classifier as identity evidence.
