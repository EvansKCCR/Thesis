# Automatic evolutionary benchmarking — workflow 1.0

Run commands from the directory containing `adhesome_workflow.py`:

```bash
cd /mnt/c/Users/eaasa/Desktop/Evans_PHD_Thesis/Integrin_analysis/Thesis/integrin_adhesome_analysis
conda activate adhesome-workflow-v1
```

The supplied `04_orthology_phylogenetic_classification/scripts/benchmark.py`
retains the [Open_Orthobench](https://github.com/davidemms/Open_Orthobench)
pairwise scorer. Its new `--workflow` mode imports the workflow adapter, restores
reference identifiers and produces machine-readable reports. No additional Python
package is required. Cite Emms & Kelly, *Benchmarking Orthogroup Inference Accuracy:
Revisiting Orthobench*, [doi:10.1093/gbe/evaa211](https://doi.org/10.1093/gbe/evaa211),
and the original Trachana et al. benchmark, as requested by the upstream project.

## Prepare the original benchmark files

Extract the original benchmark bundle under **`inputs/reference_raw`** (plural
`inputs`). Retain its `Input` and `RefOGs` directories:

```text
inputs/reference_raw/
  BENCHMARKS/
    Input/
      <the 12 original reference protein FASTAs>
    RefOGs/
      RefOG001.txt
      ...
      RefOG070.txt
      low_certainty_assignments/
        <any supplied uncertain assignments>
  OtherReference.faa                  # optional additional reference proteome
  OtherReference.gff3                 # when needed for that ordinary proteome
```

Do not rename the protein IDs or replace the original proteomes with newer
releases. Gene symbols, human crosswalk components and predicted groups are not
ground truth. The original protocol checks the supplied benchmark's 12 proteomes,
70 RefOGs, 1,945 reference memberships and 251,378 input proteins. Incomplete or
replaced bundles cannot supply a validated audit weight. Files are supplied by
the user; the workflow does not silently download a different dataset.

The workflow discovers a unique `RefOGs` folder under `reference_raw` and the
adjacent `Input` folder. A flat arrangement with matching original FASTAs beside
`RefOGs` also works. Multiple benchmark bundles require explicit selection:
`--benchmark-refogs PATH --benchmark-input-dir PATH`. No JSON is written manually.

Benchmark proteomes are automatically included in the inference panel with
`Bench_` species codes. Original benchmark proteins remain separate units and
are not collapsed to guessed genes. Preparation records this protein-level scope,
original-to-prepared IDs, sequence hashes and RefOG provenance. Ordinary reference
proteomes still use explicit gene IDs, GFF or mapping sidecars as before. Exact
duplicate FASTA files are not inserted twice.

## Run automatically after orthology

After returning external topology/localization evidence:

```bash
python adhesome_workflow.py resume --handoff runs/handoff --context-dir returned/context --reference-raw inputs/reference_raw --out runs/resumed --threads 4
```

For completed motif annotation, reuse it without repeating detection or motifs:

```bash
python adhesome_workflow.py continue --annotation runs/annotation --reference-raw inputs/reference_raw --out runs/downstream_benchmarked --threads 4
```

The order is reference preparation → panel → orthology → evolutionary benchmark
→ final audit. **Phylogeny is disabled by default.** No pause or hand-written panel
configuration is needed between these stages.

The benchmark reads full root `N0.tsv`, falling back to full `Orthogroups.tsv`;
candidate-filtered HOG tables are never scored. Exact protein IDs and sequences
link benchmark inputs to selected prepared references. All benchmark reference
proteins, including background proteins outside the curated RefOGs, must be
represented by the inference panel. Target proteins and unrelated reference
species are excluded from the benchmark projection. They cannot be counted as
false positives against truth that does not cover them.

The adapter uses the supplied pairwise equations, including normalization by
RefOG size and exclusion of uncertain assignments. Zero true-positive pairs are
handled explicitly. Duplicate memberships, ambiguous mappings, changed sequences,
missing ground truth or insufficient root-group coverage cannot produce a
validated weight. The original minimum 50% clustered-input coverage is retained.

## Reports and audit weighting

For an automatic downstream run:

```text
downstream/benchmark/
  summary.json                 # status, precision, recall, F1, measured weight
  identifier_mapping.tsv       # expected original IDs and prepared matches
  refog_scores.tsv             # per-RefOG pair errors and fragmentation
  upstream_scores.log          # output from the supplied scorer
  workflow_manifest.json       # linkage to the exact orthology/panel
downstream/final/
  joined_evidence.tsv
  audited_candidates.xlsx
```

The optional detail files exist only when their processing step was reached.
`summary.json` always records the outcome. Its validated weight is the measured
pairwise F1 on the reference cohort, on a 0–1 scale. This is a run-level weighting
factor, not a probability that a target family assignment is correct.

Audit rows preserve benchmark status, source, F1, precision, recall and RefOG
count. The measured weight scales only the orthogroup-context score. Raw and
weighted scores and the weight used remain visible. Independent curated identity,
concordant RBH identity, calibrated profiles and candidate-specific phylogenetic
identity are not discounted by a benchmark of orthogroup partitions. Architecture,
identity and topology gates still determine conservative classification; a high
benchmark score cannot establish an unsupported family identity.

Weights are applied only when report hashes, original inputs, predictions and
panel linkage verify. Unavailable, invalid or failed benchmarking records a neutral
unvalidated baseline; it never fabricates an accuracy score or blocks the final
audit. Missing data should be corrected in a new appropriate panel/run. For a
deliberately custom independently curated benchmark, explicitly select
`--benchmark-protocol custom` (`--protocol custom` on the standalone benchmark
command). Such a report is labelled custom and does not claim the original
Open_Orthobench dataset's validation scope.

## Benchmark an already completed orthology run

This reuses inference results, provided their panel already included the matching
original reference proteomes:

```bash
python adhesome_workflow.py benchmark --orthology runs/annotation/downstream/orthology --reference-raw inputs/reference_raw --out runs/evolutionary_benchmark
python adhesome_workflow.py audit --annotation runs/annotation --orthology runs/annotation/downstream/orthology --benchmark runs/evolutionary_benchmark --out runs/final_benchmarked
```

Adding benchmark proteomes after inference requires a new panel and orthology run;
the script cannot retrospectively infer their orthogroups. A standalone `orthology`
command automatically writes its benchmark report inside its output directory,
and `audit` discovers that report without an extra flag.

## Optional trees and recovery

Add `--phylogeny` to `resume` or `continue` to run trees **after the completed
audit**. Use `--phylogeny --phylogeny-grouping family` for unified family trees.
Missing tools, alignment QC failures and tree failures are recorded as optional
failures while the audit and workflow completion record remain available.
`--skip-phylogeny` remains an alias for the default. `--allow-partial-phylogeny`
controls the tree runner's partial summary; it is not required to obtain an audit.

```bash
python adhesome_workflow.py recover --run runs/resumed/downstream --threads 2
python adhesome_workflow.py recover --run runs/resumed/downstream --phylogeny --iqtree-threads 2
python adhesome_workflow.py recover --run runs/resumed/downstream --retry-benchmark
```

The second command explicitly opts in to or retries optional trees. Verified
orthology, benchmark and audit stages are reused. A previous completion manifest
is archived before an explicit optional-tree retry; completed results are not
deleted. Older recovery files do not automatically opt in to phylogeny.

A failed benchmark execution is retried by recovery without repeating completed
orthology. `--retry-benchmark` explicitly requests recomputation even for a sealed
report and refreshes its dependent audit; prior completion/report outputs are
archived. It cannot add proteomes absent from the completed inference panel.
