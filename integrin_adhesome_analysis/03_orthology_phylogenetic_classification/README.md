# Five-species OrthoFinder pipeline

For downstream per-HOG alignment, trimming, supported gene trees and interpretation,
see [PHYLOGENY.md](PHYLOGENY.md) and `scripts/phylogeny.py`.

Run inference on the complete gene-level proteomes of human, C. elegans,
S. haematobium, S. japonicum, and S. mansoni. Filter candidate families only
after inference. The original inputs and suggested_workflow.txt are preserved.

## Current preparation

The local preparation on 2026-09-20 found:

| Species | Raw proteins | Representative genes | Status |
|---|---:|---:|---|
| Hsap | 136,807 | 20,076 | Prepared |
| Cele | 28,577 | 19,984 | Prepared |
| Shae | 14,700 | 9,431 | Prepared |
| Sjap | 16,936 | 10,089 | Prepared |
| Sman | 10,935 | 9,896 | Prepared |

All 477 IDs in all_schistosoma_candidates.txt and all 89 fibronectin-like IDs
map to prepared inputs. These are separate sets; do not assume one contains the other.
The C. elegans input is gzip-compressed despite its .fa suffix; the parser detects
compression from file contents. Existing files immediately under
02_primary_proteomes are raw copies and are not used by this pipeline.

The supplied Homo_sapiens.gff identifies GRCh38.p14 (GCF_000001405.40), RefSeq
annotation release RS_2025_08. Its CDS protein accessions match all 136,807 FASTA
accessions, with no missing/conflicting mappings and no extra GFF protein IDs.
The generated human_protein_gene.tsv maps exact protein accessions to NCBI GeneIDs;
human_protein_gene.provenance.json records source hashes and annotation metadata.
The human FASTA contains RefSeq protein accessions but no explicit gene mapping.
For future input updates, obtain the matching assembly/release annotation and provide a
tab-separated file with `protein_id` and `gene_id` columns. Use exact versioned
protein accessions (for example NP_000005.3); gene IDs must uniquely identify genes.
Every input protein must be mapped. Do not infer genes by trimming protein accession
versions or parsing protein descriptions. Missing proteins are listed in
03_qc/missing_gene_mapping.tsv.

## Prepare and validate

Use Python 3.11 or newer. From this project directory:

```bash
python scripts/refseq_gene_map.py --gff 01_raw_proteomes/Homo_sapiens.gff --fasta 01_raw_proteomes/Homo_sapiens.faa --output human_protein_gene.tsv
python scripts/pipeline.py prepare --human-map human_protein_gene.tsv
python scripts/test_pipeline.py
```

Without the mapping, `prepare` still prepares the other species and writes QC,
but returns a nonzero status and marks inference blocked. It retains the longest
protein for each explicit gene ID, breaking length ties by original protein ID.
Prepared headers are `species__gene_id`. The complete isoform-to-representative
mapping is recorded in 03_qc/identifier_map.tsv. A candidate referring to an
alternative isoform maps to the chosen representative of the same gene.

Only one terminal stop is removed. Internal stops, invalid residues, duplicate
protein IDs, empty sequences and unsafe gene IDs fail validation. Proteome counts
and checksums document the supplied inputs; they do not independently prove
assembly completeness. Record download URLs, assembly accessions, annotation
releases and download dates alongside your raw files when known.

Outputs:

- 02_primary_proteomes/prepared/: validated gene-level FASTAs only.
- 03_qc/proteome_summary.tsv: sequence and gene counts.
- 03_qc/candidate_audit.tsv: every candidate and its mapping status.
- 03_qc/preparation.json: readiness, input/output hashes and selection policy.

## Install and run on Linux, WSL2 or an HPC system

The environment file pins OrthoFinder 3.1.5. Installation and a complete inference
run have not been executed on this Windows host. Transfer the whole project to the
Linux execution host, or access it through WSL, then run from the project directory:

```bash
conda env create -f environment.yml
conda activate integrin-orthofinder
conda list --explicit > 03_qc/conda-linux-explicit.txt
python scripts/pipeline.py prepare --human-map human_protein_gene.tsv
python scripts/pipeline.py run --name combined_01 --threads 16 --analysis-threads 8 --dry-run
python scripts/pipeline.py run --name combined_01 --threads 16 --analysis-threads 8
```

Choose thread counts appropriate to the machine or scheduler allocation. The
explicit Conda export records exact platform-specific dependency builds for later
recreation with `conda create -n integrin-orthofinder-replay --file ...`.
The runner checks readiness and input hashes, uses exactly five FASTAs, and gives
OrthoFinder a new output path. It does not reuse the pre-existing 04_results directory
as the direct `-o` target. Logs, help/version output, command, manifest and exit status
are saved under 03_qc/combined_01. The actual OrthoFinder results are under
04_results/combined_01; inspect that directory to identify the generated Results_* folder.
Do not automatically choose a result by lexicographic date sorting.

## Extract candidate families

Supply the exact completed results directory containing Orthologues and
Phylogenetic_Hierarchical_Orthogroups. For example, replace the placeholder below:

```bash
python scripts/pipeline.py extract --results 04_results/combined_01/Results_DATE --candidate-set fibronectin_like_protein --output 05_candidate_results/fibronectin
python scripts/pipeline.py extract --results 04_results/combined_01/Results_DATE --candidate-set all_schistosoma_candidates --output 05_candidate_results/all_candidates
```

The extractor uses exact comma-separated identifiers and exact species columns.
It writes candidate_HOGs.tsv, direct_orthologues.tsv, candidate_report.tsv,
candidate_duplications.tsv when available, and provenance.json. Unmatched
candidates remain visible in the report. Additional sets can be added as .txt
files in candidate_ids (optional protein_id header); rerun preparation first.

The extractor automatically uses N0.tsv when present, or
Orthogroups/Orthogroups.tsv when Log.txt confirms OrthoFinder 3.1 or newer.
The completed 3.1.5 run uses the latter layout; N1/N2/N3 tables describe descendant
nodes and must not substitute for the five-species root table. Combined species
orthologue tables and separate pairwise tables are supported and duplicate
relationships are collapsed. For clade-specific HOG analysis,
inspect the labelled species tree and the corresponding node table separately.
Pairwise relationships are inferred from Orthologues tables, never from shared
HOG membership. One-to-many and many-to-one are oriented from the Schistosoma
species to human or C. elegans; group sizes reflect each original orthology row.
Multiple relationships are retained. HOG copy counts are per species per selected
HOG, not a count of the candidate list.

“Schistosoma-only in sampled HOG” describes these five sampled proteomes; it does
not establish evolutionary restriction. A candidate absent from N0 is reported as
not found there, rather than conclusively unassigned. Duplication rows retain the
original support values and require gene-tree review. Pfam architecture, sequence
identity, coverage and final functional interpretation require separate evidence
and are deliberately not fabricated from OrthoFinder output.

Before interpreting results, inspect Species_Tree, assignment statistics,
Resolved_Gene_Trees and duplication events. Domain sharing alone does not establish
whole-protein orthology for the fibronectin-like candidates.

## Upstream references

- Input preparation: https://orthofinder.github.io/OrthoFinder/tutorials/gettting-input-data/
- Results: https://orthofinder.github.io/OrthoFinder/tutorials/guide-to-results/
- Software: https://github.com/OrthoFinder/OrthoFinder
- Package versions: https://anaconda.org/bioconda/orthofinder

The original suggested workflow is retained as background; use the commands here
to avoid double prefixes, uncollapsed isoforms, duplicate input files and treating
the candidate-file header as a protein ID.
