# Adhesome-subset orthology benchmark — version 1.0

Performance must pass the current policy (pairwise precision >=0.80 and F1 >=0.70)
before an accepted benchmark weight is exported. Measured scores remain visible
when performance is rejected. See [BENCHMARK_RECALIBRATION.md](BENCHMARK_RECALIBRATION.md)
for acceptance, tree-level comparison and inference tuning.

Use `--scope adhesome` on the standalone benchmark command, or
`--benchmark-scope adhesome` on `resume`, `continue` and `orthology`.
The original `benchmark.py` remains the default (`all` scope). Explicitly selecting
`adhesome` routes execution to the separate `orthology_phylogeny/scripts/adhesome_benchmark_subset.py`.
The subset script defaults to adhesome scope and shares the original scoring
equations through the workflow adapter. An original-benchmark failure does not
switch automatically to subset scoring; the user must select it. Scope selects
curated reference orthogroups for scoring; it does not reduce the inference
proteomes, alter family classification rules or repeat OrthoFinder inference.

## Prepare the curated selection

Keep the complete original Open_Orthobench `Input` and `RefOGs` bundle under
`inputs/reference_raw/BENCHMARKS/`. Original dataset integrity checks still apply
to the complete bundle even when only selected RefOGs are scored. Independently
curated custom bundles require explicit `--protocol custom` (or
`--benchmark-protocol custom` in a resumed workflow).

Copy the header-only template `templates/adhesome_refog_selection.tsv` to
`inputs/reference_raw/adhesome_refog_selection.tsv`. Populate one independently
curated RefOG/family association per row, using these tab-separated columns:

| Column | Required content |
| --- | --- |
| `refog` | Existing ground-truth filename stem, without `.txt`. |
| `family` | Exact canonical family ID from `annotation/family_rules_expanded.json`. |
| `curation_source` | Traceable publication, curated record or curation-file reference supporting this RefOG/family association independently of the evaluated predictions. |
| `curation_status` | `curated`. Proposed or provisional ground truth is not scored. |

Fibronectin scopes `fibronectin_like`, `fibronectin_III_ligand` and
`fibronectin_like_candidate` are accepted explicitly, including when a scope is
preserved as a hypothesis rather than a configured audit gate. Add an association
only if curated ground truth supports that scope. Shared FN3 domains alone are
insufficient. The template contains no invented RefOG-to-family links; the
original RefOG membership files do not supply these labels themselves.

Several RefOGs can map to one family. Their original partitions remain separate:
they are never merged into a single truth group. A RefOG can have more than one
curated family association if its documented scope justifies each; it is scored
once in the overall subset result. Avoid assigning a specific family label to
ground truth curated only at a broader family scope.

Use `--selection PATH` on the standalone command, or `--benchmark-selection PATH`
in the workflow, to override the default selection path. Missing, empty, duplicate,
unknown-family or unknown-RefOG selections cannot silently fall back to all groups
or produce an audit weight. A source citation is recorded for traceability; the
software cannot independently establish that the supplied curation is correct.

The subset entry point also accepts `--results`, `--panel`, `--reference-raw`,
`--selection`, `--protocol` and `--out` directly. For audit-ready sealed reports,
use `adhesome_workflow.py benchmark --scope adhesome` as shown below; it supplies
the verified result/panel paths and invokes `adhesome_benchmark_subset.py`.

## Reuse completed orthology

From the main validation installation in WSL:

```bash
cd /mnt/c/Users/eaasa/Desktop/Evans_PHD_Thesis/Integrin_analysis/Thesis/integrin_adhesome_analysis
conda activate adhesome-workflow-v1
python adhesome_workflow.py benchmark \
  --orthology runs/annotation/downstream/orthology \
  --reference-raw inputs/reference_raw \
  --scope adhesome \
  --out runs/adhesome_benchmark
python adhesome_workflow.py audit \
  --annotation runs/annotation \
  --orthology runs/annotation/downstream/orthology \
  --benchmark runs/adhesome_benchmark \
  --out runs/final_adhesome_benchmarked
```

Use new output directories. These commands reuse completed inference; the panel
must already contain the corresponding complete benchmark proteomes. Adding new
proteomes requires a new panel/inference run. The default expanded audit rules
remain selected. No phylogeny is required.

For a new resumed workflow, add the following to its normal command:

```text
--benchmark-scope adhesome --benchmark-selection inputs/reference_raw/adhesome_refog_selection.tsv
```

For a checkpointed run, `recover --run PATH --retry-benchmark` accepts the same
two flags and saves them for later recovery. Existing recovery behavior archives
replaced stage outputs. To retain existing runs without creating such archives,
use the standalone benchmark and new audit directories shown above.

## Scoring and outputs

The scorer reads the full root `N0.tsv` or `Orthogroups.tsv`. Complete predicted
membership within the benchmark-proteome universe is retained: unselected RefOG
members and proteins outside every curated RefOG still penalize incorrect merges
with selected groups. Groups and members are never filtered to adhesome candidates.
As in the original workflow, proteins belonging to other, unlabelled target or
reference proteomes are outside the ground-truth comparison universe; they are
not assumed to be false positives simply because the benchmark lacks their labels.
Missing selected truth members count against recall, and the supplied
`low_certainty_assignments` exclusions remain RefOG-specific.

Outputs include:

- `summary.json`: scoped precision/PPV, recall, F1, family metrics, coverage,
  provenance hashes and validation status.
- `refog_scores.tsv`: selected RefOG pair counts, fragmentation, missing members
  and exact-group correctness.
- `family_coverage.tsv`: expanded family catalogue plus the three fibronectin
  scopes, with `evaluated`, `too_few_certain_proteins` or
  `no_curated_ground_truth` status. Untested families have blank accuracy fields.
- `selection.tsv`: the selected associations, curation sources and RefOG hashes.
- `identifier_mapping.tsv` and `upstream_scores.log`: full cohort ID checks and
  the original Open_Orthobench scorer output for the selected groups.

Only validated per-family F1 values weight orthogroup-context evidence for that
exact audited family. A refined family cannot inherit another family's weight.
Untested families have a neutral baseline weight and are explicitly reported as
`family_not_benchmarked`; this is not a validated accuracy estimate. Zero measured
F1 is retained as zero. The subset aggregate is descriptive, not a global audit
weight or a protein identity probability. Domain, identity and negative-exclusion
gates remain unchanged, and benchmarking failure cannot block the final audit.

Changing the selection, ground truth, proteomes, predictions or relevant code
invalidates provenance until benchmarking is refreshed. Subset metrics do not
claim coverage of all adhesome families or independent accuracy on novel lineages.
