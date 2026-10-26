# Optional manual route

The default resume command now runs downstream stages automatically. Use resume --annotation-only only if you deliberately want this manual route. See ORTHOLOGY_PHYLOGENY_GUIDE.md for the recommended automatic run.

# Orthology and phylogeny: follow this sequence

All commands below are run from **the integrated package root**, the folder
containing `adhesome_workflow.py`. Do not change into a numbered stage folder.
Examples use bash in Linux/WSL and single-line commands. `TargetA` is an example
user-configured target; `RefA` is an example reference code. Replace them consistently.

## 1. Open the correct terminal and directory

On this Windows machine, open a WSL terminal and run:

```bash
cd /mnt/c/Users/eaasa/Desktop/Evans_PHD_Thesis/Integrin_analysis/Thesis/integrin_adhesome_analysis
conda activate adhesome-workflow-v1
pwd
ls adhesome_workflow.py
```

On another Linux/HPC machine, copy the whole package and use `cd /path/to/integrin_adhesome_analysis` instead.
If the environment does not exist, first run `conda env create -f environment.yml`.
Conda installation itself is a prerequisite. Full tool installation has not been
validated here. Check that programs are available:

```bash
command -v python orthofinder mafft trimal iqtree2
```

Stay in this directory for every command below. Paths in commands resolve from
this directory. FASTA paths inside a panel JSON resolve from the JSON's directory;
keep panel JSON files at this package root to make the two conventions agree.

## 2. Know which files you already have

Complete `start` first, using your `targets.json`:

```bash
python adhesome_workflow.py start --targets-config targets.json --input-dir inputs/targets --out runs/handoff
```

If `runs/handoff` already completed successfully, reuse it; do not rerun this command
into the same folder. Orthology uses the complete prepared target FASTAs:

```text
runs/handoff/detection/TargetA/prepared/TargetA.faa
```

For several targets, there is one such path per species code. If you used
`start --detection` to reuse another batch, the complete target FASTAs remain
under that supplied batch directory; use those paths in your panel JSON.

**Do not use** `external_candidates.faa`, `fibronectin_like_candidates.faa` or
`families/*.faa` as OrthoFinder proteomes. Those are candidate subsets, not
complete gene-level proteomes. The workflow extracts candidates after inference.

External topology/motif return can run separately. It is required for a final
combined audit but is not a prerequisite for reference preparation or orthology.

## 3. Prepare each reference species

You need these two matching inputs from the same assembly/annotation release:

| File | Contents |
|---|---|
| `inputs/reference_raw/RefA.faa` | Complete protein FASTA, including isoforms if supplied |
| `inputs/reference_raw/RefA.protein_gene.tsv` | Explicit protein-to-gene mapping, one row per FASTA protein |

Create folders for your inputs (these contain user data, not bundled resources):

```bash
mkdir -p inputs/reference_raw
```

Download/copy the reference FASTA and its matching mapping into that directory.
The TSV must use literal tabs and these column names, for example:

```text
protein_id	gene_id
protein001.1	gene001
protein001.2	gene001
protein002.1	gene002
```

Here the displayed columns are tab-separated. `protein_id` must exactly match the
first token after `>` in a FASTA header, including accession versions. `gene_id`
is the explicit biological gene identifier, not a guessed accession prefix.
Use letters, numbers, underscores, dots or hyphens in gene IDs. Every FASTA
protein must be mapped; missing mappings stop preparation. Do not type literal
backslash-t characters into the TSV.

If you have an NCBI-style GFF with CDS `protein_id` and `Dbxref=GeneID:...`, this
helper can create the TSV (it is not a universal GFF parser):

```bash
python orthology_phylogeny/scripts/refseq_gene_map.py --gff inputs/reference_raw/RefA.gff --fasta inputs/reference_raw/RefA.faa --output inputs/reference_raw/RefA.protein_gene.tsv
```

Otherwise obtain the mapping from the annotation provider. Then prepare the reference:

```bash
python orthology_phylogeny/scripts/prepare_reference.py --species RefA --fasta inputs/reference_raw/RefA.faa --map inputs/reference_raw/RefA.protein_gene.tsv --out inputs/reference_gene_level/RefA
```

This creates:

```text
inputs/reference_gene_level/RefA/
  RefA.faa             # one longest protein per gene; headers RefA__gene_id
  identifier_map.tsv   # all original proteins and chosen representatives
  preparation.json     # input/output hashes and selection policy
```

Length ties are resolved by original protein ID. One terminal stop is removed;
internal stops, duplicate protein IDs and invalid sequences are rejected.
Repeat with RefB, RefC, etc. Each species needs a unique code and a new output
folder. If you already have a complete validated gene-level reference FASTA with
`RefA__gene_id` headers, use it directly and retain its preparation provenance.

## 4. Write the panel JSON at the package root

Create `panel.initial.json` with your actual target/reference paths. A minimal
one-target/one-reference example is:

```json
{
  "proteomes": [
    {"species":"TargetA", "fasta":"runs/handoff/detection/TargetA/prepared/TargetA.faa", "gene_level":true, "source":"Target assembly/release and source"},
    {"species":"RefA", "fasta":"inputs/reference_gene_level/RefA/RefA.faa", "gene_level":true, "source":"Reference assembly/release and source; preparation.json alongside FASTA"}
  ]
}
```

Include every target declared in the handoff, even targets with no candidates,
and at least one additional reference species. Replace example provenance with
real assembly/release and source details. There is no fixed species count.

Build and validate the panel:

```bash
python adhesome_workflow.py calibrate --handoff runs/handoff --panel panel.initial.json --out runs/panel_initial
```

This creates `runs/panel_initial/proteomes/`, `species.json`, `targets.json`,
`candidate_audit.tsv` and the panel manifest. This step validates/stages inputs;
it does **not** run OrthoFinder. Candidate sequences must exactly match the
complete target proteomes. Do not edit the staged panel after validation.

## 5. Run orthology

```bash
python adhesome_workflow.py orthology --panel runs/panel_initial --out runs/orthology_initial --threads 16 --analysis-threads 8
```

Adjust threads to your machine. The wrapper runs OrthoFinder, locates its result
folder and extracts candidate evidence. You do not need to manually navigate to
`Results_DATE` for a wrapper-run analysis. Successful outputs include:

```text
runs/orthology_initial/
  inference/Results_.../           # complete OrthoFinder outputs
  candidates/candidate_HOGs.tsv
  candidates/direct_orthologues.tsv
  candidates/candidate_report.tsv
  candidates/candidate_duplications.tsv  # only when source duplication data exist
  calibration_review.tsv
  workflow_manifest.json
```

If inference was run elsewhere, follow the external-result import instructions
in WORKFLOW.md. Do not silently substitute a different panel or inference run.

## 6. Prepare and run phylogenies from that same panel/run

Use **all files from `orthology_initial` and `panel_initial` together**:

```bash
python orthology_phylogeny/scripts/phylogeny.py prepare --hogs runs/orthology_initial/candidates/candidate_HOGs.tsv --orthologues runs/orthology_initial/candidates/direct_orthologues.tsv --duplications runs/orthology_initial/candidates/candidate_duplications.tsv --proteomes runs/panel_initial/proteomes --species-list runs/panel_initial/species.json --output runs/phylogeny_initial
python orthology_phylogeny/scripts/phylogeny.py run --project runs/phylogeny_initial --threads 8
```

Preparation creates one folder per HOG with sequences and member metadata;
`run` executes MAFFT → trimAl → IQ-TREE. Only eligible groups proceed through
inference; review each group's QC/status. Trees and alignment outputs are under
`runs/phylogeny_initial/<HOG>/inference/`.

If `candidate_duplications.tsv` is absent, inspect the original OrthoFinder
`Gene_Duplication_Events/Duplications.tsv` and extraction provenance. The current
phylogeny preparation requires actual duplication-event input; do not manufacture
an empty table. A source table genuinely containing zero events is different
from missing evidence. Resolve the missing source/extraction before proceeding.

No outgroup is assumed. Specify one only with biological justification. Trees
are review evidence and are not automatically converted to family audit scores.

## 7. Add references without mixing runs

Prepare the added reference with step 3. Copy `panel.initial.json` to
`panel.expanded.json`, append its entry, then run:

```bash
python adhesome_workflow.py calibrate --handoff runs/handoff --panel panel.expanded.json --out runs/panel_expanded
python adhesome_workflow.py orthology --panel runs/panel_expanded --compare runs/orthology_initial --out runs/orthology_expanded --threads 16 --analysis-threads 8
```

Review `calibration_comparison.tsv`. To rebuild trees, repeat step 6 replacing
**every** `_initial` with `_expanded`. Do not combine old HOGs with new proteomes.
Use a new folder for each panel/run; already completed outputs are not overwritten.

## 8. Final evidence audit

After `resume` has produced `runs/annotation`, choose the orthology run you reviewed:

```bash
python adhesome_workflow.py audit --annotation runs/annotation --orthology runs/orthology_initial --out runs/final_initial
```

For expanded references, use `runs/orthology_expanded` and `runs/final_expanded`.
The annotation and orthology must derive from the same handoff. See WORKFLOW.md
for the external topology/localization return and motif evidence requirements.