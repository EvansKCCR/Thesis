# Automated orthology and phylogeny after resume

Run all commands from the package root containing `adhesome_workflow.py`.
On this machine, in WSL:

```bash
cd /mnt/c/Users/eaasa/Desktop/Evans_PHD_Thesis/Integrin_analysis/Thesis/integrin_adhesome_analysis
conda activate adhesome-workflow-v1
```

Install the environment from `environment.yml` first if needed. OrthoFinder,
MAFFT, trimAl, IQ-TREE and the Python dependencies must be available. No reference
panel JSON or manual gene-level FASTA preparation is required for the default run.

## Put the reference inputs in one folder

The default folder is **`input/reference_raw`** (singular `input`), relative to
the package root. Override it with `--reference-raw /path/to/references`.
Use a unique species code as each FASTA filename stem:

```text
input/reference_raw/
  RefA.faa
  RefA.gff3
  RefB.fasta
  RefB.protein_gene.tsv
  RefC.fa
```

Reference FASTAs must contain the complete reference proteome, not candidate
subsets. `.fa`, `.faa`, `.fasta` and their `.gz` variants are accepted. Files are
read directly from this folder, not recursively. Codes start with a letter and
contain letters, digits or single underscores; they cannot reuse a target code.

The workflow resolves gene identities automatically in this order:

1. Matching `Code.protein_gene.tsv` or `Code.tsv`, if supplied. Columns are
   `protein_id` and `gene_id`, separated by tabs; protein IDs exactly match FASTA
   header first tokens. Duplicate mapping files are rejected.
2. Matching `Code.gff` or `Code.gff3`. Supported annotations are NCBI-style CDS
   records carrying `protein_id` and `Dbxref=GeneID:...`. The annotation must match
   the FASTA assembly/release. Other GFF schemas need an explicit mapping TSV;
   the workflow does not guess gene relationships.
3. If neither is supplied, explicit `gene=ID` or `gene:ID` FASTA header fields,
   or already prepared `Code__gene_id` identifiers.

Every protein must have an unambiguous gene ID. The longest protein per gene is
selected (ties by protein ID), one terminal stop is removed, and prepared IDs
become `Code__gene_id`. Invalid/duplicate sequences, missing mappings and ambiguous
sidecars stop with an error. The raw files remain unchanged. Automatically
created mapping tables, source hashes and selection audits are retained.

## Resume once and let the workflow continue

After the detection handoff and external topology/localization analysis:

```bash
python adhesome_workflow.py resume --handoff runs/handoff --context-dir returned/context --out runs/resumed --threads 8
```

For topology alone, replace `--context-dir returned/context` with
`--topology-file returned/predicted_topologies.3line`. Existing return contracts
and exact-sequence checks remain active. No JSON needs to be written by the user.

The command executes, in order:

1. Import returned context and run motif annotation.
2. Discover and prepare all reference proteomes in `input/reference_raw`.
3. Reuse the handoff's complete target proteomes, verify their hashes, and write
   the panel configuration automatically.
4. Stage the panel and run whole-proteome OrthoFinder.
5. Extract candidate orthology, join evidence and write the conservative audited workbook.
6. Independently prepare per-HOG projects and run MAFFT → trimAl → IQ-TREE.

Supported/provisional and fibronectin-like hypotheses all follow this same route.
Targets with no exported candidates remain part of the whole-proteome panel.
The external predictors are the planned handoff; no further manual preparation
pause is required when valid inputs and dependencies are available.

## Find the outputs

```text
runs/resumed/
  motifs/                           # annotation output
  workflow_manifest.json            # annotation completion record
  downstream/
    references/                     # prepared references and source provenance
      panel.auto.json               # generated; do not edit
    panel/                          # exact inference inputs and target/reference roles
    orthology/                      # full inference + extracted candidate evidence
    phylogeny/                      # per-HOG alignment, QC and tree outputs
    final/audited_candidates.xlsx    # final workbook
    workflow_manifest.json          # full downstream completion record
```

Groups with fewer than four sequences are recorded as ineligible for supported
tree inference. If no candidate groups are found, an explicit no-groups status
is written. Missing OrthoFinder duplication data remain explicitly **unknown**
in phylogeny manifests/status files; they do not prevent sequence-based trees
and are not interpreted as evidence of no duplication. QC/tool failures return
nonzero and preserve diagnostic logs; they are not silently marked complete.

## Continue an existing annotation or expand references

If motif annotation already finished, reuse it directly:

```bash
python adhesome_workflow.py continue --annotation runs/resumed --out runs/downstream_retry --threads 8
```

This uses the same automatic reference preparation and downstream sequence.
Use recover for a checkpointed failure (see below). For a new downstream analysis, use a new output directory. Continue reuses annotation
but reruns downstream computation; it does not resume midway through an external
OrthoFinder/IQ-TREE process.

To expand the panel, add the new complete reference FASTA and matching annotation
(if needed) to `input/reference_raw`, then run:

```bash
python adhesome_workflow.py continue --annotation runs/resumed --out runs/expanded --compare runs/resumed/downstream/orthology --threads 8
```

The new run builds a fresh panel, reruns inference/trees/audit and writes
`runs/expanded/orthology/calibration_comparison.tsv`. Earlier runs remain intact.
No panel JSON editing or repeated motif annotation is needed.

Optional controls: `--annotation-only` on `resume` restores an intentional stop
after motif annotation; `--skip-phylogeny` runs reference preparation, orthology
and audit without custom gene-tree inference. The lower-level manual modules
remain available in `MANUAL_ORTHOLOGY_PHYLOGENY.md` for specialist use.
## Phylogeny failures and recovery

IQ-TREE now defaults to **two threads independently of MAFFT/OrthoFinder**.
Override using `--iqtree-threads` on `resume`, `continue`, or `phylogeny.py run`.
The workflow writes `phylogeny/run_summary.tsv` and `run_summary.json`, with
separate complete, failed and skipped groups.

| Symptom | Recorded category | Action |
|---|---|---|
| IQ-TREE nonzero exit, including 2 | tool_failure at iqtree | Inspect the stderr path and tail in run.json, command stderr log, and gene_tree.log if present. Exit 2 alone does not prove a resource problem. |
| Excessive gaps, no standard residues, insufficient informative sites | alignment_qc_failure | Inspect aligned.sequence_qc.tsv / trimmed.sequence_qc.tsv and alignment quality; lowering threads cannot repair the data. |
| Alignment/trimming identifier changes | identifier_mismatch | Inspect membership JSON and trimal_identifiers.tsv. Short reversible IDs protect trimAl input; genuine membership changes still fail validation. |
| Fewer than four input sequences | skipped / too_few_taxa | Intentional exclusion from supported tree inference; not an execution failure. |

There is no automatic deletion of problematic taxa, relaxation of QC thresholds,
or fabricated tree for a failed alignment. A family/HOG with heterogeneous domain
architecture may require a biologically justified alignment strategy or curated
homologous region selection before rerunning in a new prepared project.

For an interrupted run whose downstream completion manifest has **not** been
written, retry only failed groups after reviewing their errors:

```bash
python 04_orthology_phylogenetic_classification/scripts/phylogeny.py run --project runs/resumed/downstream/phylogeny --threads 2 --iqtree-threads 2 --retry-failed
```

Previous failures move to `inference_attempt_001`, etc.; completed groups are
preserved. A run marked running/incomplete is not automatically overwritten.
After successful recovery, the final audit can be invoked without repeating
OrthoFinder:

```bash
python adhesome_workflow.py audit --annotation runs/resumed --orthology runs/resumed/downstream/orthology --out runs/final_recovered
```

Strict failure remains the default. To deliberately finish the audit even if some
groups fail in a new run, pass `--allow-partial-phylogeny` to resume/continue.
The phylogeny summary retains `partial_failure`, and the downstream manifest
records `completed_with_phylogeny_failures`. Such an audit does not imply all trees
succeeded. Do not mutate already sealed completed downstream runs; create a new
run/project for subsequent revisions. Low-level `phylogeny.py run --allow-partial`
is the equivalent explicit option when using the standalone runner.
## Resume a failed run from its checkpoints

After fixing the reported cause, run from the package root:

```bash
python adhesome_workflow.py recover --run runs/resumed --threads 2 --iqtree-threads 2
```

You may instead pass `--run runs/resumed/downstream`. For a run created by
`continue --out runs/expanded`, use `--run runs/expanded`. The failed run's saved
annotation, references, panel and settings are reused; optional thread flags
override resource settings for the retried work. Do not launch recovery while
another process is still writing to that run.

`recovery.json` records running/failed/completed checkpoints and file hashes for:
reference preparation, panel construction, orthology, phylogeny preparation,
phylogeny inference and final audit. Each completed stage is verified and reused.
A changed/missing completed output is rejected rather than silently reused.
Recovery is for the same scientific inputs/settings; use a new `continue` run
for reference-panel expansion or changed biological inputs.

- If extraction failed after OrthoFinder succeeded, its verified inference
  outputs are reused and only extraction/downstream stages are retried.
- If tree inference failed, references, panel, orthology and phylogeny preparation
  are reused. Failed tree groups are retried; completed groups are retained.
- If final audit failed, completed computational stages are reused.
- If an external OrthoFinder invocation itself failed, that invocation is rerun.
  Recovery is not a promise to resume inside OrthoFinder's search algorithm.
- Incomplete non-resumable stage folders are archived as `*_attempt_001`, etc.
  Old outputs/logs are preserved. Tree attempts retain their per-group archives.
- Annotation must have completed before these downstream checkpoints exist. If
  annotation itself failed, rerun `resume` into a new annotation output folder;
  the completed detection handoff is still reused.

A successfully completed sealed run is checked and left unchanged. Recovery
rejects tree attempts still marked running/incomplete rather than assuming an
active process has stopped; inspect those before attempting a retry.

### Older failed runs without recovery.json

For a pre-checkpoint run where reference preparation and panel validation
completed, adopt its verified existing stages:

```bash
python adhesome_workflow.py recover --run runs/resumed/downstream --adopt-existing --annotation runs/resumed --threads 2 --iqtree-threads 2
```

Adoption validates annotation/panel linkage, generated panel consistency,
orthology manifests if present, and prepared phylogeny input hashes. It then
creates checkpoints and resumes. It will not adopt a run already sealed complete.
If an older run failed before a valid reference panel was produced, use `continue`
from its completed annotation to create a new downstream run instead.
## Audit is independent of phylogeny

The audit checkpoint now runs immediately after successful orthology, before phylogeny preparation, dependency checks or inference. Missing tree tools and failed trees therefore cannot prevent an available annotation/orthology audit. Tree failures still return nonzero and remain recoverable; the workbook is retained in downstream/final. For an existing run, invoke the standalone audit command to obtain the workbook without retrying any trees.
