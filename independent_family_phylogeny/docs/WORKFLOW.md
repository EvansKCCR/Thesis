# Independent family phylogeny: operating workflow

Run commands from the bundle root, not docs/. See [MANUAL.md](MANUAL.md) for
family membership and evidence interpretation and [MAINTENANCE.md](MAINTENANCE.md)
for recovery and quality diagnostics. No guide command runs automatically.

## Local setup and planning (WSL/Linux)

Run commands in WSL from this directory:

```bash
cd /mnt/c/Users/eaasa/Desktop/Evans_PHD_Thesis/Integrin_analysis/Thesis/independent_family_phylogeny
```

Use an existing environment if Python, MAFFT, trimAl and IQ-TREE 2 or later are
already installed and on `PATH`. Otherwise create the separate environment:

```bash
conda env create -f environment.yml
conda activate independent-family-phylogeny
```

Inspect file integrity and the planned commands; these do not run analysis or
write results:

```bash
python -B check_library.py
python -B run_phylogeny.py
```

`run_phylogeny.py` only executes external tools when **`--execute`** is supplied.
The runner locates `mafft`, `trimal`, and `iqtree` / `iqtree2` / `iqtree3` on PATH.
All paths are resolved from this bundle, so the source installation is not
needed to run its prepared FASTAs. The preparation and inference runner use the
Python standard library. The optional branch-evidence annotator additionally
requires Biopython, now included in `environment.yml`.


## Run all prepared families locally

**Check whether a driver is active before launching analysis.** Use the queued
continuation below to add the six new families without starting a competing
analysis or restarting actin. The generic commands in this section are for a
fresh run or a run whose previous driver has exited.

```bash
python -B run_phylogeny.py --execute
```

Each family is processed separately and families run sequentially. The default
commands are equivalent to:

```bash
mafft --thread 4 --ep 0 --genafpair --maxiterate 1000 input.fasta > aligned.fasta
trimal -in aligned.fasta -out trimmed.fasta -gappyout
iqtree -s trimmed.fasta -m MFP -B 1000 -T AUTO
```

The runner additionally sets the IQ-TREE output prefix and reproducible seed
`20261028`. It does not pass `--redo` or `-mredo`: interrupted IQ-TREE runs can
reuse their native `tree.ckp.gz` and/or ModelFinder `tree.model.gz` checkpoints.
Completed stages are reused only if their inputs, executable hash,
command and output checksums still match. It records commands and separate
stdout/stderr logs. No alignment, trimming or inference from the adhesome
workflow is invoked or reused.

IQ-TREE uses automatic thread selection as requested. To cap local resource use:

```bash
python -B run_phylogeny.py --execute --run-name limited_threads \
  --mafft-threads 4 --iqtree-threads 4
```

Run a single family, or repeat `--family` to select multiple families:

```bash
python -B run_phylogeny.py --execute --family src_like
python -B run_phylogeny.py --execute --family rho_gtpase_like --family fes_fps_like
python -B run_phylogeny.py --execute --family hsp70_like
```

To omit the mixed validation/control panel while retaining every actual
proteome-panel record, choose a distinct run name:

```bash
python -B run_phylogeny.py --execute \
  --run-name proteome_panels_only --exclude-panel validation
```

This selection omits exactly 15 records across the 15 inputs. The original
prepared library remains intact, and each run records its selected identifiers
and excluded panels.


## Alternative alignment / trimming and manual review

L-INS-i with the explicit occupancy threshold requested:

```bash
python -B run_phylogeny.py --execute --run-name localpair_gt02 \
  --alignment localpair --trim gt02
```

These options are independent: `--alignment localpair` switches MAFFT to
`--localpair --maxiterate 1000`; `--trim gt02` switches trimAl to `-gt 0.2`.
Use a new `--run-name` for a different combination or different resource options.

Stop after alignment or trimming for inspection:

```bash
python -B run_phylogeny.py --execute --stop-after align
python -B run_phylogeny.py --execute --stop-after trim
```

Then run the same command with `--stop-after tree` (the default) to continue.
The runner rejects changed identifiers, missing/all-gap sequences, unequal
alignment lengths and MAFFT residue changes. It also stops tree inference if
fewer than four distinct residue sequences remain. These are integrity checks,
not a biological alignment-quality assessment. Profilin starts with only four
distinct sequences, so it has little redundancy for sequence exclusion.


## Add the six prepared families after the active run

The already-running Python process loaded its nine-family list at startup.
Changing `families.json` cannot change that in-memory list. The additions are
therefore prepared in the existing `library/` and recorded in
`pending_families.json`; this descriptor records their preparation, not their current execution status.

To activate a local continuation, keep the current driver running and use a
second WSL terminal in the same original environment:

```bash
conda activate adhesome-workflow-v1
cd /mnt/c/Users/eaasa/Desktop/Evans_PHD_Thesis/Integrin_analysis/Thesis/independent_family_phylogeny
python -B run_queued_families.py --execute --wait
```

The continuation waits for **all nine original tree completion records and the
original driver's final run status**. While waiting it launches no analytical
tools. Once the original driver finishes, it verifies/reuses its outputs and
runs the six additions sequentially under the same
`results/proteome_panels_only/` directory. It preserves `genafpair`, `gappyout`,
MAFFT threads 4, IQ-TREE threads AUTO, seed 20261028 and exclusion of `validation`.
The expanded final run status covers all 15 families. A lock prevents two copies
of this queued continuation from starting together.

Leave this continuation terminal open. Closing it leaves the active original
run unaffected; repeat the continuation command later. If the original batch
reports failure, the continuation stops without starting the additions so that
the failed original work can be recovered first. It is not launched by this
preparation step or registered as an OS/Codex background task.

To preview the queued commands without executing or waiting:

```bash
python -B run_queued_families.py
```

`prepare_library.py --append` supports future additive extensions. It requires
the exact original orthogroup table and source manifest, refuses changes to
existing family definitions, verifies all existing exported hashes, and writes
only new family input directories plus refreshed shared library metadata. The
old family inputs and active run directory are not rewritten. Its manifest
extension history records previous manifest hashes and preserved family-file
hashes as provenance; it does not create backups. The existing branch annotator
will automatically include all 15 families after their trees are available.


## Generate branch evidence from completed family trees

`annotate_branches.py` is a standalone post-processing command. It reads existing
Newick trees, their `run_config.json` selections, the prepared tip metadata and
the matching full OrthoFinder evidence. It does not run alignment, trimming,
tree inference, rooting, species-tree reconciliation or the adhesome workflow.
It creates a new annotation directory and refuses to overwrite existing output.

If this environment was created before Biopython was added to `environment.yml`,
install that dependency once:

```bash
conda activate independent-family-phylogeny
conda install -c conda-forge biopython
```

Run locally from this bundle:

```bash
python -B annotate_branches.py \
  --run-dir results/proteome_panels_only \
  --orthofinder-results /mnt/c/Users/eaasa/Desktop/Evans_PHD_Thesis/pipeline_validation/integrin_adhesome_analysis/runs/annotation/downstream/orthology/inference/Results_Oct17 \
  --out branch_annotations/proteome_panels_only \
  --skip-incomplete
```

`--skip-incomplete` explicitly skips missing family trees and records the skipped
families in `annotation_summary.tsv`. It does not skip identifier, provenance or
membership errors. Omit this flag to require all prepared family trees. To process selected
families only, repeat `--family`, for example `--family src_like --family actin`.
If you renamed `tree.treefile` to `gene_tree.treefile`, add
`--tree-name gene_tree.treefile`. Other run metadata and library names should
retain their existing names.


Only annotate finished family trees. `--skip-incomplete` skips missing tree files;
it does not establish that an existing intermediate tree has finished. Check
`tree.complete.json` for each selected family and use explicit `--family` choices
if an unfinished family already has a `tree.treefile`.

To defer actin locally, first stop any waiting queued continuation, then interrupt
the original foreground driver once with Ctrl+C and wait for its shell prompt.
Retain native checkpoints. Confirm no actin IQ-TREE process remains before starting
other families against the same run. Use the original environment and parameters:

```bash
conda activate adhesome-workflow-v1
python -B run_phylogeny.py --execute --resume \
  --run-name proteome_panels_only --exclude-panel validation \
  --alignment genafpair --trim gappyout \
  --mafft-threads 4 --iqtree-threads AUTO --seed 20261028 \
  --family shp_ptpn6_11_like --family insr_like --family talin_like \
  --family hsp70_like --family laminin_like --family collagen_like --family tetraspanin_like
```

Completed families are checked/reused. When that driver exits, the full resume
command in [MAINTENANCE.md](MAINTENANCE.md) can return to actin. The queued
continuation deliberately waits for all original nine families and is not the
route for bypassing deferred actin. Do not leave competing drivers active.

## Optional independent RBH / main-audit reconciliation

The [reconciliation guide](../../reconciliation/docs/WORKFLOW.md) documents a
separate candidate handoff and `library_reconciled` family-tree library. The
three scientific workflows remain independently runnable. Existing audit
classifications, family libraries and analysis checkpoints are unchanged. Use
the new `rbh_audit_reconciled_v1` run name for this library; its source annotations
and short tree identifiers are traceable in the accompanying TSV files.

## Dataset-free deployment planning

```bash
python -B run_phylogeny.py --dry-run --clean-state
```

This checks the workflow structure without a prepared library or original
orthology results. Ordinary planning requires a real prepared library. A real
run requires `--execute`; it cannot be combined with `--dry-run`. Existing
`families.json` and the queued extension retain historical OG identifiers for
their original source run. Review mappings against a new OrthoFinder result
before using them; matching OG numbers from different runs do not establish
family identity. Fresh RBH/main-audit handoffs use the reconciliation adapter.
