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
5. Extract candidate orthology, prepare per-HOG projects and run MAFFT → trimAl → IQ-TREE.
6. Join evidence and write the conservative audited workbook.

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
Use a new output directory after correcting a failure. This reuses annotation
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