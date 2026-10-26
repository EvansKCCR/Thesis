# Automatic orthology, evolutionary benchmarking and independent audit

Run every command from the package root containing `adhesome_workflow.py`:

```bash
cd /mnt/c/Users/eaasa/Desktop/Evans_PHD_Thesis/Integrin_analysis/Thesis/integrin_adhesome_analysis
conda activate adhesome-workflow-v1
```

The default flow is annotation → reference preparation → panel → orthology →
evolutionary benchmark → final audit. **Phylogeny is opt-in** and runs after the
audit. No manual panel JSON is required.

## Reference files and benchmark data

The default folder is **`inputs/reference_raw`** (plural `inputs`). Override it
with `--reference-raw PATH`. Ordinary reference FASTAs use safe species-code
filenames such as `RefA.faa`; `.fa`, `.faa`, `.fasta` and gzip variants are accepted.

For ordinary reference proteomes, automatic gene preparation uses a matching
`Code.protein_gene.tsv`/`Code.tsv`, matching NCBI-style `Code.gff`/`Code.gff3`
(CDS protein_id and Dbxref=GeneID), or explicit `gene=ID`/`gene:ID` headers.
Already prepared `Code__gene_id` headers also work. Every protein needs an explicit
mapping; the longest isoform per gene is selected and selection hashes are kept.
Gene relationships are never guessed from arbitrary identifiers.

For the original Open_Orthobench dataset, preserve its layout and protein IDs:

```text
inputs/reference_raw/
  BENCHMARKS/
    Input/                         # all 12 original protein FASTAs
    RefOGs/                        # RefOG001.txt through RefOG070.txt
      low_certainty_assignments/   # when supplied
  OtherReference.faa               # optional ordinary reference proteome
  OtherReference.gff3
```

The benchmark inputs are discovered automatically, prepared without collapsing
original benchmark protein units, and included in OrthoFinder. Their protein-level
scope is explicit in the generated panel. Do not supply invented gene mappings
or create ground truth from the inferred orthogroups. The original RefOGs are
required for validated scores. See [benchmark details](EVOLUTIONARY_BENCHMARKING.md)
for original-dataset integrity checks, file selection and audit weight policy.

## Resume or continue

After detection and external topology/localization prediction:

```bash
python adhesome_workflow.py resume --handoff runs/handoff --context-dir returned/context --reference-raw inputs/reference_raw --out runs/resumed --threads 4
```

For topology alone, use `--topology-file returned/predicted_topologies.3line`
instead of `--context-dir`. Use the exact handoff associated with the returned
sequences. `--annotation-only` deliberately stops after motif annotation.

If annotation already completed, reuse its actual directory:

```bash
python adhesome_workflow.py continue --annotation runs/annotation --reference-raw inputs/reference_raw --out runs/downstream_benchmarked --threads 4
```

Use new output directories. Supported/provisional and fibronectin-like hypotheses
follow the same path. Targets with no exported candidates remain in the inference
panel. Reference preparation, panel JSON generation, orthology, benchmarking and
audit proceed without an intervening manual preparation pause.

```text
runs/resumed/downstream/
  references/panel.auto.json       # generated automatically
  panel/                          # exact target/reference inference inputs
  orthology/                      # inference and extracted candidate evidence
  benchmark/summary.json          # status, reference-cohort metrics and weight
  benchmark/identifier_mapping.tsv
  benchmark/refog_scores.tsv
  final/Stage_2_output.xlsx    # final result, independent of optional trees
  final/joined_evidence.tsv
  recovery.json
  workflow_manifest.json
```

A valid benchmark's F1 weights orthogroup context only. Missing or incompatible
benchmark data receive an explicit status and no validated weight; the audit
still completes using its neutral unvalidated baseline. Architecture and resolving
identity gates remain mandatory. Per-candidate family identity is not inferred
from the benchmark's aggregate score.

## Optional phylogeny

Add `--phylogeny` to resume/continue to run MAFFT → trimAl → IQ-TREE after audit.
For unified family trees add both `--phylogeny --phylogeny-grouping family`.
Grouping alone does not enable trees. `--skip-phylogeny` remains a compatibility
alias for the default.

IQ-TREE defaults to two threads independently of OrthoFinder/MAFFT. Override with
`--iqtree-threads 2`. Missing tree tools, alignment QC failures and nonzero tree
commands are recorded in status/recovery files while the completed audit remains
available. The workflow returns a completion record with optional failures.
`--allow-partial-phylogeny` controls the tree runner's summary; it is not needed
for the audit to complete.

| Symptom | Interpretation |
| --- | --- |
| IQ-TREE nonzero exit, including 2 | Inspect stderr and gene_tree.log; the exit code alone does not establish a thread/resource problem. |
| Excessive gaps/no residues | Alignment QC failure; reducing threads does not repair alignment quality. |
| Changed alignment identifiers | Membership failure; inspect transport-ID mapping and membership logs. |
| Fewer than four sequences | Intentional tree skip. |

No taxa or alignment QC constraints are silently removed to make a tree succeed.
See [family-tree instructions](FAMILY_PHYLOGENY.md) for standalone preparation.

## Recovery and panel expansion

```bash
python adhesome_workflow.py recover --run runs/resumed/downstream --threads 2
python adhesome_workflow.py recover --run runs/resumed/downstream --phylogeny --iqtree-threads 2
python adhesome_workflow.py recover --run runs/resumed/downstream --retry-benchmark
```

Recovery validates and reuses completed checkpoints. The second command explicitly
opts in to or retries optional trees and archives a previous downstream completion
manifest before doing so. Completed orthology, benchmark and audit are reused;
failed tree groups retain their attempt logs. Older checkpoints do not implicitly
opt in to phylogeny. Do not recover while another process is writing to the run.

Failed benchmark execution is also retried while reusing completed orthology.
`--retry-benchmark` explicitly recomputes a report and refreshes the dependent
audit, archiving prior outputs. It does not infer new proteomes retrospectively.

For older failed runs without a checkpoint, use `recover --adopt-existing
--annotation ACTUAL_ANNOTATION_DIRECTORY` after verified reference/panel preparation.
For changed reference proteomes, use a new `continue` run and optionally
`--compare PREVIOUS_ORTHOLOGY_DIRECTORY`; earlier panels remain unchanged.
Adding benchmark proteomes after inference also requires a new panel/inference.

To benchmark and audit completed orthology whose panel already included the
original benchmark proteomes, without repeating inference:

```bash
python adhesome_workflow.py benchmark --orthology runs/annotation/downstream/orthology --reference-raw inputs/reference_raw --out runs/evolutionary_benchmark
python adhesome_workflow.py audit --annotation runs/annotation --orthology runs/annotation/downstream/orthology --benchmark runs/evolutionary_benchmark --out runs/final_benchmarked
```

Use directories that actually exist: `workflow_manifest.json` must be present in
the annotation directory and the corresponding orthology/panel provenance must
match. An example path such as `runs/resumed` is not created by recovery itself.
