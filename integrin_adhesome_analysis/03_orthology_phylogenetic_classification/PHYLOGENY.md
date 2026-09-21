# Per-family phylogenetic analysis

`scripts/phylogeny.py` builds a separate protein gene tree for each selected row
of candidate_HOGs.tsv. It includes **all** HOG members from all five prepared
proteomes, not just the matched candidates. Unrelated HOGs are never concatenated.
This is a reproducible inference workflow for publication work; alignment quality,
domain homology, rooting and biological interpretation still require review.

## Prepared fibronectin families

The real inputs have already been extracted into `06_phylogenies/fibronectin`:
23 groups, of which 19 contain at least four sequences. Four three-sequence groups
are retained but skipped by the supported-tree workflow. No external alignment or
tree inference has been run by Codex. The supplied candidate_duplications.tsv has
no event rows; this is recorded as absent supplied evidence, not absence of duplication.

To reproduce extraction into a **new** output directory, run from orthofinder_project:

```bash
python scripts/phylogeny.py prepare \
  --hogs 05_candidate_results/fibronectin/candidate_HOGs.tsv \
  --orthologues 05_candidate_results/fibronectin/direct_orthologues.tsv \
  --duplications 05_candidate_results/fibronectin/candidate_duplications.tsv \
  --proteomes 02_primary_proteomes/prepared \
  --output 06_phylogenies/fibronectin_new
```

Optionally repeat `--hog OG0000137 --hog OG0000211` to select specific rows. Every
sequence must match an exact prepared FASTA ID. Missing members fail extraction.
An existing output directory is never overwritten. Each group contains members.faa,
tips.tsv (species/candidate labels), orthology_edges.tsv and duplication_evidence.tsv.
groups.tsv summarizes all groups; manifest.json records input and extracted file hashes.

## Install and infer on WSL/Linux

```bash
conda env create -f environment-phylogeny.yml
conda activate integrin-phylogeny
conda list --explicit > 06_phylogenies/fibronectin/conda-linux-explicit.txt
python scripts/phylogeny.py run --project 06_phylogenies/fibronectin --threads 8
```

The separate environment pins MAFFT, trimAl, IQ-TREE and Biopython. Its installation
has not been tested on this host. If your executable is named iqtree or iqtree3,
pass `--iqtree iqtree` or `--iqtree iqtree3`; the actual version output is recorded.
To run a single family first, add `--hog OG0000137`.

Each group's inference directory contains:

- aligned.faa: MAFFT amino-acid L-INS-i (`--localpair --maxiterate 1000`). Iterative
  refinement uses `--threadit 0` to reduce multithreading variability.
- trimmed.faa: trimAl `-automated1`, plus retained-column and sequence QC reports.
- gene_tree.treefile, gene_tree.iqtree, gene_tree.log and support outputs: IQ-TREE
  ModelFinder (`-m MFP`), 1,000 UFBoot replicates, 1,000 SH-aLRT replicates and
  `-bnni`. Seed defaults to 20260920; identical sequences are retained (`-keep-ident`).
- branch_evidence.tsv: tip-set-defined branches, original support labels, counts
  of supplied direct orthology pairs within each clade/split and cautious interpretation.
- run.json and command logs: settings, executable versions, commands, input/output
  hashes and completion or failure status.

The script checks preservation of sequence IDs and residues, equal alignment lengths,
trimmed-column provenance, at least 50 retained sites, no sequence with >90% gaps,
at least four distinct aligned sequences and at least one parsimony-informative site.
These are screening thresholds, not guarantees of reliable inference. They are
configurable with `--min-sites` and `--max-gap`; inspect failures before changing them.
Failures are logged per group and other eligible groups continue. Existing inference
directories are protected; for a revised analysis prepare a new project and select
the relevant HOGs. No silent resume or overwriting is performed.

## Rooting and interpretation

The inferred tree is unrooted by default. Do not assume human or C. elegans is a valid
outgroup for every family. For a biologically justified outgroup, run a single HOG
with `--outgroup EXACT_PREPARED_ID` (repeat for multiple outgroup tips). Rooting is
performed after inference and the chosen set must form an exclusive clade. The
original unrooted tree is preserved, with rooted.treefile written separately.

Use the following evidence together:

| Evidence | Interpretation and limit |
|---|---|
| Direct orthology edge | An OrthoFinder-inferred candidate/reference relationship; membership of one HOG alone is insufficient. The supplied table covers candidate-to-human/C. elegans comparisons, not all tip pairs. |
| Supported mixed-species clade | Consistent with shared ancestry; requires rooting and reconciliation to distinguish orthologues from ancient paralogues. |
| Supported single-species clade in a rooted tree | A possible lineage expansion. Establish duplication placement relative to speciation and check gene models before describing recent paralogues. |
| Duplication event | Original OrthoFinder reconciliation evidence. Old node IDs are not transferred to the new topology. Exact matches between complete daughter gene sets and rooted child clades are recorded. Other cases need manual reconciliation. |
| Missing duplication row | No evidence in the supplied candidate-filtered file; not evidence that no duplication occurred. Consult the full Duplications.tsv for family-wide events. |

Branch IDs B0001 etc. refer to the row's explicit tip set and are not IQ-TREE or
OrthoFinder node IDs. Support labels are retained as reported by IQ-TREE (typically
SH-aLRT/UFBoot); verify their order in the .iqtree report. Approximate thresholds
such as SH-aLRT >=80 and UFBoot >=95 are screening conventions, not posterior
probabilities or proof of orthology. Identical sequences can have unresolved order.

For fibronectin-like proteins, examine domain architectures and homologous coverage
before interpreting full-length trees. Repeated domains and nonhomologous regions
can mislead an otherwise well-supported tree. `--aligner einsi` enables MAFFT
E-INS-i for large unalignable regions. Domain-restricted analyses require separately
validated homologous coordinates; this script does not guess them. Compare alignment
and trimming sensitivity in a new project and inspect IQ-TREE composition tests,
convergence and model selection in its report. Reconcile supported gene trees with
an independently justified species tree before making duplication-timing claims.

## Methods and references

Validation: four regression tests pass in the existing WSL of3_env, covering
FASTA/QC checks, rooting and duplication matching, and the full inference driver
with simulated external tool outputs. This does not constitute a real alignment
or IQ-TREE benchmark. Run `python scripts/test_phylogeny.py` in the phylogeny
environment to repeat all tests. The inspected of3_env contains MAFFT, IQ-TREE
and Biopython but lacks trimAl; use the separate environment above for inference.

Report the actual software versions, selected substitution model from each .iqtree
file, original/retained alignment lengths, sequence counts, seed, rooting rationale,
support methods and filtering decisions. Archive input hashes, alignments, trees,
logs and the explicit Conda environment export. Do not report inferred outcomes
before those jobs finish and pass review.

- MAFFT algorithms: https://mafft.cbrc.jp/alignment/software/manual/manual.html
- trimAl automated1: https://trimal.readthedocs.io/en/latest/_static/supplementary_material.pdf
- IQ-TREE commands: https://www.iqtree.org/doc/Command-Reference
- IQ-TREE support guidance: https://www.iqtree.org/doc/Frequently-Asked-Questions
