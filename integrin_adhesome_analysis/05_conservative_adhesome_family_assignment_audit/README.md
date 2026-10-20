# Stage 2 family refinement — version 1.0

`integrated_adhesome_pipeline.py` refines original stage 1 family hypotheses using
the selected family/identity pair and `audit_policy.json`. The integrated workflow
defaults to `family_rules_expanded.json` and `family_identity_rules_expanded.json`.
One representative is selected per species/protein. Original alternatives and
collapsed/recovery counts remain auditable. No new proteins or component calls
are created. Experimental annotation is not required.
Missing experimental evidence never lowers assignment or annotation confidence.
An unambiguous grade A call is strong for both family and protein identity at the
reported scope, using the passed family gate and computational orthology support.
CDD is mandatory only for families whose specificity rule requires it; optional
CDD is not a further confidence cap.

From the workflow root:

```bash
python adhesome_workflow.py audit --annotation runs/annotation --orthology runs/annotation/downstream/orthology --benchmark runs/annotation/downstream/benchmark --out runs/final_two_stage
```

The result is `runs/final_two_stage/Stage_2_output.xlsx`, with the nine sheets in
the supplied stage 2 template. Orthology columns use actual panel species.
Architecture/orthology/motif/localization maxima are 4/3/2/1. Complete architecture
plus orthology gives A; complete architecture without enough orthology gives B;
partial architecture plus orthology gives provisional C; insufficient architecture
gives unassigned D. Diagnostic CDD and ordered architectures constrain shared
families. Only complete PINCH/parvin hypotheses may be recovered in the same pool.

Current motif context support/conflicts are preserved. Superseded motif IDs are
excluded. Verified benchmark F1 weights orthology only, with raw/weighted scores
exported. Optional phylogeny runs after audit and cannot prevent completion.
See [the full protocol](../../README/COMPUTATIONAL_AUDIT.md) for gates, scoring,
outputs, evidence scope and cached stage 1 refresh instructions.

Run `python -B -m unittest discover -p 'test_*.py'` here for audit regression tests.


The final table now uses `putative_audited_family` and `protein_name`, with one
`assignment_confidence` field and no MULocDeep, empty or unused columns. A matching
`audited_family_assignments.tsv` is exported. Observed CDD footprints complement
Pfam architecture without double-counting overlapping domain copies. The current
scoped identity refinements use the selected identity rules; inferred example
features are not observations. Positive adhesome controls alone calibrate the
reported family/name correctness estimate, with uncertainty and control counts.
Hard/background negatives act only as exclusion gates; matching proteins are
dropped from final assignments and recorded in `negative_control_exclusions.tsv`.
Retained plus excluded proteins must account for every stage 1 input protein.
Shared domains cannot trigger exclusion by themselves. The bundled model uses
reference-overlapping controls and does not establish novel-lineage accuracy.
See the full protocol for recalibration and reference-panel requirements.

## Run the expanded pair through the main audit

From the installation containing your completed `runs/annotation` directory:

```bash
python adhesome_workflow.py audit \
  --annotation runs/annotation \
  --orthology runs/annotation/downstream/orthology \
  --family-coverage expanded \
  --out runs/final_expanded
```

This audits all existing Stage 1 candidates using both expanded JSON files.
It reuses annotation/orthology and writes a new audit rather than restarting
recruitment. Expect the message `Audit family coverage: expanded (190 family
configurations; 136 identity rules).` The result is
`runs/final_expanded/Stage_2_output.xlsx`; the validation target's metrics are in
`runs/final_expanded/validation_evaluation/` when that target is present.
Use an unused output directory. Expanded coverage is the default for `audit`,
`resume` and `continue`; the flag above is optional. Original rules remain
available with `--family-coverage original` for comparison only. Completed
annotation should be reused with `audit` above.

The selected pair is recorded under the audit's `configuration/` and verified
by both the audit and evaluation receipts. A numerical probability model fitted
to the original pair cannot be applied to expanded calls. Computational family
confidence remains governed by the observed family gates and orthology.
Stage 1 and its Pfam scan library are unchanged: this option cannot supply
missing domain observations or recruit Defensin-6.

## Expanded family coverage: isolated development comparison

From the `integrin_adhesome_analysis` directory, use the existing handoff and
completed audit's full `joined_evidence.tsv`. Replace `/path/to/existing/installation`
with the installation containing those runs; choose a new output directory:

```bash
python evaluate_family_expansion.py \
  --handoff /path/to/existing/installation/runs/handoff \
  --evidence /path/to/existing/installation/runs/annotation/downstream/final/joined_evidence.tsv \
  --out development/family_coverage_review
```

The script installs `family_rules_expanded.json` and
`family_identity_rules_expanded.json` together in the development output's
`installed_pair/`. It replays only Stage 2 for the existing validation target,
using unchanged observed Pfam/CDD, orthology, motif and topology evidence. The
baseline and expanded arms share the current classifier, frozen gold labels and
negative gates. Expected family labels never supply missing target domains.
Numerical calibration is disabled for this paired experiment; it is not refitted.

Read `report.md`, `regression_checks.tsv`, `new_identity_positives.tsv` and
`unresolved_identity_proposals.tsv`. Each arm also reports membership versus
family correctness, broad versus specific support, provisional calls,
family macro-F1, reference/calibration overlap and negative-exclusion ablations.
All 235 positive controls remain in membership-recall denominators, including
the unrecruited screening-only Defensin-6. The 115 identity proposals and 16
screening-only positives are evaluated separately.

The tool rejects unsupported identity-rule fields or parent routes and verifies
source checksums before/after replay. It neither installs expanded rules into
production nor runs Stage 1, external tools or phylogeny. Expanded JSON rules do
not add missing Pfam HMMs: inspect `missing_model_review.tsv` and the scan-library
coverage in `summary.json` before considering a production library update.

## Maintained validation tools

Use workflow-root `evaluate_family_expansion.py` for paired rule validation and
`evaluate_validation.py` for labelled control evaluation. Evolutionary benchmarking
uses the maintained scripts in `04_orthology_phylogenetic_classification/scripts/`.
The unused one-off expansion generator/checker and obsolete benchmark scorer with
its uncurated placeholder panel were removed on 2026-10-19. Expanded rules,
expansion provenance, labelled calibration data and active benchmark inputs remain.
