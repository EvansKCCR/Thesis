# Post hoc refinement of unresolved Stage 2 annotations — version 1.0

This independent step reads **Stage_2_output.xlsx produced by the family audit**,
using its `Audited_family_assignment` worksheet. It reinterprets existing Pfam
and NCBI CDD annotations for unresolved, ambiguous and provisional proteins.
Original audit assignments, names, support, grades, probabilities, motifs and
benchmark evidence remain in their original columns. Proposals are exported
in separate `posthoc_*` columns and a separate review table.

Python 3.10+ is sufficient for the standalone script, which uses only the
standard library. The main workflow entry point requires its usual Python 3.11+
environment. No domain search, BLAST, OrthoFinder, alignment or phylogeny is
launched. This optional step does not change audit acceptance or performance
metrics automatically.

## Run after a completed audit

From `Thesis/integrin_adhesome_analysis`:

```bash
python -B main/adhesome_workflow.py refine-unresolved \
  --input /absolute/path/to/audit/Stage_2_output.xlsx \
  --out runs/posthoc_domain_review_v1
```

The standalone entry point works from any directory:

```bash
python -B /absolute/path/to/Thesis/schisto_pfam_cdd_reclassify/schisto_pfam_cdd_reclassify.py \
  --input /absolute/path/to/audit/Stage_2_output.xlsx \
  --outdir /absolute/path/to/new/posthoc_domain_review_v1
```

The audit's `audited_family_assignments.tsv` or an equivalent CSV is also accepted.
`--sheet` overrides the XLSX worksheet and is ignored for TSV/CSV. Choose a new
output directory outside the sealed source audit; existing directories are
refused. No backups or automatic archives are created.

To validate schema, support labels, duplicate keys and record accounting without
creating proposals or writing files:

```bash
python -B /absolute/path/to/Thesis/schisto_pfam_cdd_reclassify/schisto_pfam_cdd_reclassify.py \
  --input /absolute/path/to/audit/Stage_2_output.xlsx --check-inputs
```

## Input contract and scope

- Required for nonempty inputs: `protein_id`, `Specie` (or `Species`/`species`),
  `family_support_level`, and at least one Pfam/CDD architecture column.
- `Pfam_architecture` also accepts `architecture` or `pfam_architecture`;
  `NCBI_CDD_architecture` accepts `cdd_architecture`. Family names accept
  `putative_audited_family` or `assigned_family`.
- Every source column is retained. Orthogroup, length, Stage 1 hypothesis,
  topology and localization fields are optional. The audit may omit entirely
  blank columns; this does not prevent review or invent missing evidence.
- Both underscore and older space-separated support labels are accepted.
  `specific_supported` and `broad_supported` are preserved and supply optional
  within-input orthogroup context. They are not reclassified.
- `unresolved`, `unassigned`, `provisional`, `ambiguous`, `conflicting`,
  `unsupported` and `not_supported` enter the review pool. Unknown/blank support
  labels fail validation rather than silently changing the review scope.
- Explicit source negative/exclusion flags receive no proposal. Proteins already
  removed by the audit are not imported again from screening or exclusion tables.
- The key is **species + protein ID**. The same ID in two species represents two
  records. Identical duplicates within a species are collapsed and logged;
  conflicting duplicates stop before output creation.
- A valid empty final audit produces header-only reports and zero counts.

## Rule interpretation

The supplied ordered `specific()` and `broad()` rule bodies are retained from
ruleset **2026-11-04-frozen-r1**. The pipeline-compatible adapter is version
**1.0**. Narrow rules run before broad domain/superfamily rules. Matching produces
a proposed family/class, rationale, functional module and review priority.

| Post hoc tier | Interpretation |
| --- | --- |
| A | Diagnostic domain heuristic matched; inspect the exact rationale. Some original rules use Pfam architecture alone. |
| B | Family-associated Pfam/domain architecture |
| C | Broad domain or catalytic superfamily; specific identity unresolved |
| D | No discriminating signature among the supplied rules |

These tiers are distinct from original audit grades and are not calibrated
probabilities or supported audit assignments. CDD model names alone do not
establish CDD specific-hit status or complete architecture. Modules are proposed
functional categories. Experimental evidence is not a prerequisite.

Supported proteins in the same single orthogroup provide descriptive context.
They cannot upgrade a C/D tier or overwrite the audited family. Multi-group
strings are not treated as one group; tied supported labels are displayed
explicitly. `src` and `SRC` compare case-insensitively; `SCR` is not an alias.
FN3-only candidates remain domain-level proposals rather than a specific
fibronectin identity.

## Outputs

| File | Content |
| --- | --- |
| `Stage_2_posthoc_refinement.tsv` | All unique source species/protein records; every original column plus `posthoc_*` fields |
| `Schisto_unassigned_Pfam_CDD_classification.tsv` | Reviewed subset, proposed family/class, rationale, Pfam/CDD and original audit labels |
| `scope_accounting.tsv` | Per-record review/preserve/exclude decision and reason |
| `duplicate_protein_ids.tsv` | Exact duplicate rows collapsed within species |
| `tier_summary.tsv` | Counts for all four tiers, including zeros |
| `family_summary.tsv` | Proposed family/class counts and species/orthogroup coverage |
| `module_summary.tsv` | Proposed module counts |
| `priority_review.tsv` | P1 subset; retains its header when empty |
| `run_manifest.json` | Input/code/ruleset provenance, column mapping, accounting and output SHA-256 values |

The source XLSX is never rewritten. TSV outputs open in Excel, R or Python.
Supported originals retain blank proposals and `supported_preserved` status;
explicit exclusions have `negative_excluded` status. No authoritative audit
family, protein identity, confidence or probability is silently substituted.

## Historical reproducibility and maintenance

The former script expected `adhesome_candidates_list.xlsx` and historical counts
(3,281 reviewed; A/B/C/D = 609/270/2,239/163). Those describe its historical run,
not expected results for current audits. `--strict` retains this optional count
check; **omit it for new Stage 2 outputs**. `--compare PRIOR_REVIEW.tsv` compares
fields by species/protein key. A requested strict/comparison failure returns
nonzero and is recorded in `run_manifest.json`.

For rule changes, update `RULESET_VERSION` and regression-test diagnostic examples.
Main audit family/identity rule changes do not automatically change this separate
rule set. Promotion into the main audit requires an explicit reviewed rule change
and a new audit. Post hoc proposals do not automatically become corrected
validation calls.

Synthetic checks from this bundle:

```bash
python -B -m unittest discover -s tests -v
```

See the main [workflow guide](../integrin_adhesome_analysis/docs/WORKFLOW.md).
The standalone entry point can review completed audits in a separate
`pipeline_validation` installation without changing that installation's code.

## Calibration and validation performance reports

Report **saved refinement outputs** against their reference labels with the
separate `report_posthoc.py` entry point. This command does not rerun refinement,
the audit, searches or inference. It does not fit or recalibrate a model.

```bash
python -B /absolute/path/to/Thesis/schisto_pfam_cdd_reclassify/report_posthoc.py \
  --refinement /absolute/path/to/completed/posthoc/calibration \
  --labels /absolute/path/to/calibration_set.xlsx \
  --sequence-manifest /absolute/path/to/calibration/inputs/sequence_manifest.tsv \
  --kind calibration \
  --out /absolute/path/to/new/posthoc_reports/calibration
```

For validation use `--kind validation`, `validation_set.csv`, and the matching
validation branch's `inputs/sequence_manifest.tsv`. The manifest maps internal
audit IDs back to exact labelled IDs. Validation joins use `fasta_id`; calibration
may use the explicit accession mapping. No fuzzy matching, gene-name joins,
identifier version stripping or prediction-derived truth is used. Every labelled
ID must have a unique manifest mapping. Unknown/unlabelled refinement records are
reported separately and excluded from metric denominators. Missing audit records
remain no-calls in the full labelled cohort, including false negatives for positives.

The matching main-workflow command is `main/adhesome_workflow.py report-refinement`
with the same options. Both entry points resolve the sibling main bundle's
existing trusted-name evaluation ontology and require Python 3.11+. Explicit
`trusted_family` labels take priority; other expected families resolve from trusted
names/genes in the label table. Prediction domains and reference-panel recruitment
labels do not determine gold identities. Exact protein/paralogue identity is not
scored by this family-level report.

Reports keep three policies separate:

- `original_supported`: original supported audit family calls, over all labels.
- `posthoc_A`: tier A proposals only, over all labels and the eligible review subset.
- `posthoc_AB`: tier A+B proposals only, over all labels and the eligible review subset.

There is **no combined or promoted audit policy**. C/D proposals remain descriptive.
Binary PPV, recall, F1, specificity and confusion counts measure issuing a call or
proposal for positive controls. Family PPV/recall/F1 additionally require a correct
mapped expected family. A wrong called family contributes both a family FP and FN.
Unmapped positive truth labels are retained in binary metrics and explicitly
excluded from family metrics. Unmapped A/B proposals count as incorrect when the
positive truth family is evaluable. Zero-denominator rates are not estimable.

The reporting-only `config/evaluation_labels.json` defines explicit equivalences
between descriptive proposal names and canonical family labels. It does not change
the classifier. Broad domain names or orthogroup suggestions are not automatically
equated to a protein identity. `src`/`SRC` are aliases; `SCR` is not. NF2/Merlin is
not merged into canonical ERM, and TRIO is not merged into spectrin. Inspect the
unmapped-label/proposal tables when expanding these evaluation equivalences.

Each new report directory contains:

| File | Content |
| --- | --- |
| `performance_report.html`, `performance_report.md` | Separate original-audit/proposal performance and record accounting |
| `binary_metrics.tsv` | Binary confusion counts, PPV, recall, F1, specificity and accuracy by policy/population |
| `family_metrics.tsv` | Family correctness metrics with evaluable/unmapped truth counts |
| `family_confusion.tsv` | Expected versus predicted family counts, including no-call/unmapped categories |
| `predictions.tsv` | One labelled record per row, preserving original calls alongside distinct proposals |
| `family_misclassifications.tsv` | Incorrect family calls/proposals and negative-control false positives |
| `unmapped_labels.tsv`, `unmapped_proposals.tsv` | Explicit label-coverage gaps |
| `unmatched_refinement_records.tsv` | Unlabelled source records excluded from metrics |
| `report_manifest.json` | Label/mapping/refinement/config/ontology hashes and reporting provenance |

Use a new report directory outside the refinement and frozen-input directories.
The original audit and refinement tables are never edited. Output checksums and
completed refinement status are verified first. Reports are descriptive evaluation
of the supplied panels; no independent generalization claim is made, and
reference/training overlap is not estimated by this reporting step.

### Positive-control orthogroup reconstruction

Add both `--orthogroups /path/to/Orthogroups.tsv` and
`--panel-proteomes /path/to/panel/proteomes` to either reporting entry point.
Use the orthogroups, FASTAs and sequence manifest from the **same saved branch**.
This counts all labelled positives regardless of their audit or proposal status.
It verifies positive-control sequences against the saved target FASTAs and reads
the selected same-species reference FASTAs. No orthology inference is launched.

The HTML and Markdown reports count recovery only when the control's saved
orthogroup contains a same-species reference protein with an identical full
amino acid sequence. Protein IDs locate records but never establish a match.
Sequence comparison ignores case and one terminal stop. No alignment, partial
sequence or similarity threshold is used. The species table distinguishes an
identical sequence recovered in the group, an identical sequence present in the
reference panel but not recovered in the group, and controls with no identical
reference sequence. Other same-species group members are shown separately.
Rates use two denominators: positives with a same-species reference, and
positives with an identical sequence available in that reference. A species
without a mapped reference is not assessable for this comparison.

Reference columns are identified by the exact species name with underscores,
or by `Bench_<species_name>_<assembly_or_release>`. The labelled target columns
are excluded from reference matching. For differently named or ambiguous inputs,
use `--reference-column 'Mus musculus=your_mouse_column'`. An empty value marks a
species not assessable. The mapping is recorded in the report manifest.

`orthogroup_reconstruction.tsv` contains the species totals and rates.
`orthogroup_reconstruction_controls.tsv` lists the assigned group, match type,
and matching reference IDs for each positive control. The manifest records the
input FASTA and orthogroup checksums. Positive labels select controls after
inference and do not enter the sequence comparison. The saved groups included
the controls as input sequences, so these counts measure recovery of known
sequences in the same group. They do not establish unseen-sequence recovery,
one-to-one orthology or exact reconstruction of a complete orthogroup.

### Blind, reference-only orthogroup reconstruction

For a blind test, use `blind_orthogroup_test.py` in two phases. The preparation
step takes reference FASTAs from the saved panel, puts only those proteomes in
`reference_proteomes/`, and writes positive controls to a separate FASTA outside
that inference directory. It checks every extracted control against the saved
sequence manifest. Run OrthoFinder only on `reference_proteomes/`. After it
finishes, evaluate the withheld controls against its `Orthogroups.tsv` by exact
full-length sequence identity. IDs are retained for locating FASTA records and
reporting matches; they are never a match criterion. This test reports only
controls with an identical reference sequence that map to one orthogroup as
reconstructed. Controls with no identical sequence in the available reference
proteomes are reported separately and are not counted as reconstruction
failures among identity-assessable controls.

Example from the validation installation in WSL (set `KIND=calibration` and
the corresponding branch when testing calibration instead):

```bash
cd /mnt/c/Users/eaasa/Desktop/Evans_PHD_Thesis/pipeline_validation/integrin_adhesome_analysis
SCRIPT=/mnt/c/Users/eaasa/Desktop/Evans_PHD_Thesis/Integrin_analysis/Thesis/schisto_pfam_cdd_reclassify/blind_orthogroup_test.py
KIND=validation
BRANCH="runs/isolated_validation/$KIND"
BLIND="runs/blind_reference_only/$KIND"

python -B "$SCRIPT" prepare \
  --panel-proteomes "$BRANCH/downstream/panel/proteomes" \
  --target-species "$BRANCH/downstream/panel/targets.json" \
  --controls "$BRANCH/inputs/sequence_manifest.tsv" \
  --out "$BLIND"

orthofinder -f "$BLIND/reference_proteomes" -t 4 -a 4 -o "$BLIND/orthofinder"

find "$BLIND/orthofinder" -type f -path '*/Orthogroups/Orthogroups.tsv' -print
ORTHOGROUPS=$(find "$BLIND/orthofinder" -type f -path '*/Orthogroups/Orthogroups.tsv' -print -quit)
test -n "$ORTHOGROUPS" || { echo 'Orthogroups.tsv not found'; exit 1; }
python -B "$SCRIPT" evaluate \
  --prepared "$BLIND" \
  --orthogroups "$ORTHOGROUPS" \
  --out "$BLIND/results"
```

Use a new `BLIND` directory each time; the script refuses an existing
destination rather than overwrite it. The OrthoFinder input folder is
`$BLIND/reference_proteomes` only. Do not point OrthoFinder at the parent
`$BLIND`, because that parent also contains the withheld controls. Results are
written to `blind_reconstruction_by_species.tsv`,
`blind_reconstruction_controls.tsv`, and `blind_test_results.json`. Exact
identity is deliberately conservative: it does not use approximate similarity,
alignment, or inferred paralogue assignment for controls lacking an identical
reference sequence.
