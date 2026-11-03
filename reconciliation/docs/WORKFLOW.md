# Independent candidate reconciliation and phylogeny

The three scientific workflows remain independent:

| Workflow | Entry point from `Thesis/` | Responsibility |
| --- | --- | --- |
| Homology retrieval | `homology_retrieval/reciprocal_blast.py` | Curated-query searches and reciprocal BLAST evidence |
| Main adhesome | `integrin_adhesome_analysis/adhesome_workflow.py` | Detection, returned context, orthology, family/identity audit |
| Independent family phylogeny | `independent_family_phylogeny/run_phylogeny.py` | MAFFT, trimAl and IQ-TREE |

`reconciliation/reconcile_candidates.py` reads completed outputs and prepares a
handoff library. It never launches any of those workflows. The pipeline's
`calibrated_v1` algorithm, family rules, benchmark and stored classifications are
unchanged. Trees are not automatically returned to the main audit.

## Fresh-deployment state

No historical RBH run, audit or prepared family library is required for
structural validation. This source bundle does not recreate deleted datasets or
previous results. An actual preparation requires an explicit completed `--audit`
and a completed current RBH run. Source decisions and exclusion gates are
preserved; preparation does not reclassify proteins.

```bash
python -B reconciliation/reconcile_candidates.py --dry-run --clean-state
```

Only queried families with eligible target candidates receive new tree inputs.
Other accepted main families remain in the reconciliation catalogue. Generated
libraries and tree results are separate from all source scientific workflows.

## Inspect and run the prepared phylogeny locally

Commands below are for WSL/Linux. Use the environment containing Python, MAFFT,
trimAl and IQ-TREE (the phylogeny bundle's `environment.yml` defines it).

```bash
cd /mnt/c/Users/eaasa/Desktop/Evans_PHD_Thesis/Integrin_analysis/Thesis

# Read-only integrity checks; source workflow folders are not needed.
python -B reconciliation/reconcile_candidates.py --check
python -B independent_family_phylogeny/check_library.py \
  --library independent_family_phylogeny/library_reconciled

# Read-only command plan; no analysis or results written.
python -B independent_family_phylogeny/run_phylogeny.py \
  --library independent_family_phylogeny/library_reconciled \
  --run-name rbh_audit_reconciled_v1 \
  --mafft-threads 4 --iqtree-threads 4

# Launch locally only when ready.
python -B independent_family_phylogeny/run_phylogeny.py \
  --library independent_family_phylogeny/library_reconciled \
  --run-name rbh_audit_reconciled_v1 \
  --mafft-threads 4 --iqtree-threads 4 --execute
```

This runs families serially with:

```bash
mafft --thread 4 --ep 0 --genafpair --maxiterate 1000 selected_input.fasta > aligned.fasta
trimal -in aligned.fasta -out trimmed.fasta -gappyout
iqtree -s trimmed.fasta -m MFP -B 1000 -T 4 -seed 20261028 \
  --prefix tree -st AA -keep-ident
```

Use `--alignment localpair` or `--trim gt02` with a **new** `--run-name` to select
the alternative strategies. `--iqtree-threads AUTO` is also supported; keep the
same setting when resuming. `-keep-ident` retains identical reference/taxon tips.
The new library uses short `Câ€¦` target and `Râ€¦` reference identifiers; the full
original IDs and source metadata are in each family's `tips.tsv`.

Select a single family with `--family parvin_like`. Repeat `--family` to select
several. A failed family is reported while other selected families continue.

## Resume the new run

```bash
python -B independent_family_phylogeny/run_phylogeny.py \
  --library independent_family_phylogeny/library_reconciled \
  --run-name rbh_audit_reconciled_v1 \
  --mafft-threads 4 --iqtree-threads 4 --execute --resume
```

Completed stages are reused after hash checks; an unfinished tree retains its
native IQ-TREE checkpoint. Use the same input library, tool environment and
options. No `--redo` is added. Continue earlier phylogeny runs from their original
bundle and output folder; this handoff does not migrate existing checkpoints.

New outputs are under
`independent_family_phylogeny/results/rbh_audit_reconciled_v1/<family>/`:
`aligned.fasta`, `trimmed.fasta`, `tree.treefile`, `tree.iqtree`, `tree.contree`,
stage logs and checkpoint receipts. No tree analysis was run during preparation.

## Rebuild the handoff after a later completed audit or RBH run

The default action is a read-only reconciliation plan:

```bash
python -B reconciliation/reconcile_candidates.py --audit integrin_adhesome_analysis/runs/YOUR_COMPLETED_AUDIT
```

For new scientific outputs, select a new library directory. This preparation
command writes FASTA/TSV handoffs only:

```bash
python -B reconciliation/reconcile_candidates.py --prepare \
  --homology-run homology_retrieval/results/homology_run \
  --homology-inputs homology_retrieval \
  --audit integrin_adhesome_analysis/runs/NEW_TARGET_AUDIT \
  --out independent_family_phylogeny/library_reconciled_v2
```

`--audit` expects the directory containing `workflow_manifest.json`,
`audited_family_assignments.tsv` and `negative_control_exclusions.tsv`.
An existing output directory is refused. Edit the explicit crosswalk when adding
new RBH reference sets; unknown sets are refused rather than guessed. No sequence
identity threshold, audit confidence, family rule or negative gate is changed.

The main handoff is read from the audit's annotation provenance. If that
installation has moved, supply its actual location explicitly:

```bash
  --handoff /mnt/c/Users/eaasa/Desktop/Evans_PHD_Thesis/pipeline_validation/integrin_adhesome_analysis/runs/handoff
```

This folder must contain `detection/<species>/prepared/<species>.id_map.csv`,
`<species>.faa` and `preparation.json`. Original WSL/Windows path spellings are
translated for the current operating system. Relocated homology inputs are read
from `--homology-inputs` and checked against their original saved hashes. Frozen
source manifests are never rewritten.

## Reconciliation rules and files

- The target cohort is the actual RBH target-proteome list, with no fixed species
  quota. Main validation/calibration panels are not target candidates.
- All main `family_supported` and `provisional_family` calls in that cohort are
  retained with their original tier. RBH-only candidates retain their recorded
  RBH status, including augmented-reference or reference-only provisional status.
- Explicit hard/background-negative exclusions apply to both candidate sources.
  A main architecture failure alone is not treated as a negative-control match.
- Original IDs are joined through the recorded gene/isoform maps and exact amino
  acid sequences. Version suffixes are not stripped to invent matches. Identical
  alternate isoforms of the same mapped gene share one tip; different isoform
  sequences remain separate. Different genes/taxa with identical sequences are
  retained separately. A disagreement in an accepted main protein's source
  sequence stops preparation rather than replacing that sequence.
- Family disagreements are visible in `family_conflict`. A target may occur in
  each explicitly claimed family tree for investigation; no winning assignment
  or confidence promotion is inferred by reconciliation.
- Original source bytes are untouched. FASTA parsing normalizes letter case,
  whitespace and terminal stop markers; internal stops in tree inputs are
  refused. Source-file and normalized-sequence hashes are recorded separately.

`library_reconciled/` contains:

| File | Purpose |
| --- | --- |
| `reconciled_candidates.faa`, `candidates.tsv` | Full accepted-main/RBH target union and unchanged source annotations |
| `rbh_main_reconciliation.tsv` | Every unique RBH input's ID match, main decision and exclusion status |
| `rbh_supporting_evidence.tsv` | Original accepted query/hit evidence rows |
| `family_summary.tsv` | Mapped family coverage, counts and preparation/skipping reason |
| `negative_exclusion_gates.tsv`, `excluded_candidates.tsv` | Applicable original negative rules and candidates removed by them |
| `sequence_conflicts.tsv` | Different isoforms or sequence-release disagreements |
| `<family>/input.fasta`, `<family>/tips.tsv` | Tree input and exact tip-to-source mapping |
| `library_manifest.json` | Source/output hashes, settings, counts and provenance |

The prepared library is sufficient to run the independent trees; the main and
RBH installations are needed only to rebuild it. Recorded orthogroup IDs remain
context from the main audit, not new orthology/duplication evidence. The existing
`annotate_branches.py` expects original OrthoFinder identifiers and complete OG
membership; do not apply it directly to the short-ID reconciled library. A
separate mapping/evidence validation is required before branch-level claims.

## Maintenance checks

Only synthetic tests are needed for adapter changes:

```bash
python -B -m unittest discover -s reconciliation/tests -v
python -B -m unittest discover -s independent_family_phylogeny/tests -v
```

No backups, BLAST searches, main audits, alignments or trees are created by these
tests. Keep earlier source snapshots and active run directories in place.
