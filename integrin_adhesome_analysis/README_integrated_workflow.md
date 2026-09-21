# Integrin Adhesome Comparative Workflow

## Overview

This repository implements an auditable comparative-proteomics workflow for identifying, inspecting, and phylogenetically analysing candidate integrin–adhesome proteins. The workflow links Pfam/domain evidence, curated reference homology, family-level motif inspection, orthology inference, and per-family phylogenetic reconstruction while preserving intermediate evidence and provenance at each stage.

```text
Prepared proteomes
  -> Pfam/HMMER scan
  -> filtered domain evidence
  -> homology / reciprocal-reference evidence
  -> comparative candidate classification
  -> family-specific candidate libraries
  -> motif / ELM inspection
  -> whole-proteome OrthoFinder analysis
  -> candidate HOG / orthologue extraction
  -> per-HOG sequence extraction
  -> MAFFT -> trimAl -> IQ-TREE
  -> annotated gene trees
```

## Workflow stages

### 1. Pfam library construction and domain scanning

A focused adhesome Pfam HMM library is generated from the configured Pfam manifest. Models with curated gathering thresholds are separated from fallback models, and metadata and file hashes are retained. Prepared proteomes are scanned with HMMER to produce retained domain calls, complete audit tables, raw `domtblout` outputs, and scan provenance.

Key outputs:
- `domains.csv`
- `domains.audit.csv`
- `ga.domtblout`
- `fallback.domtblout`
- `scan.json`

### 2. Reference homology and comparative classification

Candidate proteins are compared with curated adhesome references and a broader metazoan reciprocal-reference panel. Forward and reciprocal homology, domain architecture, family rules, motif evidence, and optional profile or phylogenetic evidence are integrated into comparative classification.

Key outputs:
- `homology.csv`
- `matched_reference_sequences.faa`
- `matched_reference_catalogue.csv`
- `candidates.csv`
- `protein_summary.csv`
- `architectures.csv`
- `motif_observations.csv`
- `candidates.faa`
- `prioritized_candidates.faa`
- `classification.json`

### 3. Family-level motif inspection

Family-specific candidate FASTAs are inspected against a versioned conservative motif library. Known matches are reported with sequence context and optional verified Pfam/topology regions. Exploratory exact-word discovery and optional ELM annotation are treated as supporting evidence rather than standalone classifiers.

Optional MEME Suite analysis can use:
- STREME for de novo protein motif discovery
- FIMO for motif localization
- TOMTOM for comparison with a genuine protein PWM database

Key outputs include:
- `known_hits.tsv`
- `inspection_states.tsv`
- `known_prevalence.tsv`
- `domain_inspection.tsv`
- `discovered_words.tsv`
- `discovered_sites.tsv`
- ELM annotation tables when enabled
- `manifest.json`
- `REPORT.md`

### 4. Whole-proteome orthology inference

Orthology is inferred using complete gene-level proteomes from *Homo sapiens*, *Caenorhabditis elegans*, *Schistosoma haematobium*, *S. japonicum*, and *S. mansoni*. Protein isoforms are collapsed to one representative per explicit gene identifier before OrthoFinder analysis. Candidate filtering is applied only after whole-proteome inference.

Key outputs:
- `candidate_HOGs.tsv`
- `direct_orthologues.tsv`
- `candidate_report.tsv`
- `candidate_duplications.tsv`
- `provenance.json`

### 5. Per-HOG phylogenetic reconstruction

Each candidate-associated HOG is analysed separately. All HOG members from the five prepared proteomes are extracted, aligned with MAFFT, trimmed with trimAl, quality checked, and used for maximum-likelihood inference with IQ-TREE.

Default settings:
- MAFFT L-INS-i: `--localpair --maxiterate 1000`
- trimAl: `-automated1`
- IQ-TREE ModelFinder: `-m MFP`
- 1,000 ultrafast bootstrap replicates
- 1,000 SH-aLRT replicates
- BNNI optimization
- reproducible fixed seed

Per-HOG outputs include:
- `members.faa`
- `tips.tsv`
- `orthology_edges.tsv`
- `duplication_evidence.tsv`
- `inference/aligned.faa`
- `inference/trimmed.faa`
- `inference/gene_tree.treefile`
- optional `inference/rooted.treefile`
- `inference/branch_evidence.tsv`
- `inference/alignment_qc.json`
- `inference/sequence_qc.tsv`
- `inference/run.json`

Trees are unrooted by default. Rooting should only be performed with a biologically justified outgroup.

## Core scripts

- `build_adhesome_pfam_library.py` — builds focused Pfam HMM libraries.
- `parse_adhesome_domtblout.py` — parses and filters HMMER domain output.
- `adhesome_pipeline.py` — orchestrates preparation, domain scanning, homology and classification.
- `adhesome_classification.py` — integrates comparative evidence and assigns candidate status.
- `adhesome_motifs.py` — shared motif inspection functions.
- `build_reference_sets.py` — builds provenance-aware curated reference sets.
- `family_motifs.py` — family-level motif inspection and exploratory discovery.
- `elm_annotation.py` — offline ELM annotation.
- `import_domain_regions.py` — imports verified upstream domain coordinates.
- `refseq_gene_map.py` — maps human RefSeq proteins to NCBI GeneIDs.
- `pipeline.py` — prepares five-species proteomes, runs OrthoFinder and extracts candidate orthology.
- `phylogeny.py` — prepares and executes per-HOG phylogenetic analyses.

## Reproducibility

The workflow retains SHA-256 hashes, exact commands, software paths and versions, stage-specific manifests, raw search output, filtered evidence tables, candidate mapping audits, alignment QC, tree QC, and per-run provenance. Existing result directories are generally protected from silent overwrite.

## Interpretation guidance

- Pfam/domain architecture supports structural compatibility.
- Sequence similarity and reciprocal evidence support homology.
- Short motifs provide corroborating context but are not sufficient for family assignment.
- Shared HOG membership does not by itself establish direct orthology.
- Pairwise OrthoFinder orthology and duplication outputs should be interpreted together with gene-tree topology.
- Phylogenetic interpretation should consider alignment quality, homologous coverage, rooting, branch support, and gene-tree/species-tree reconciliation.

For repeated-domain proteins or proteins containing substantial nonhomologous regions, inspect alignment coverage and consider domain-restricted analyses when independently validated coordinates are available.

## Suggested analysis sequence

```bash
# Build reference resources
python build_adhesome_pfam_library.py ...
python build_reference_sets.py ...

# Scan and classify
python adhesome_pipeline.py scan ...
python adhesome_pipeline.py homology ...
python adhesome_pipeline.py classify ...

# Motif inspection
python family_motifs.py --out results/motif_inspection ...

# Whole-proteome orthology analysis
python pipeline.py prepare --human-map human_protein_gene.tsv
python pipeline.py run --name combined_01 --threads 16 --analysis-threads 8

# Candidate orthology extraction
python pipeline.py extract   --results <completed_orthofinder_results>   --candidate-set <candidate_set>   --output <candidate_output>

# Per-HOG phylogeny preparation
python phylogeny.py prepare   --hogs <candidate_output>/candidate_HOGs.tsv   --orthologues <candidate_output>/direct_orthologues.tsv   --duplications <candidate_output>/candidate_duplications.tsv   --proteomes 02_primary_proteomes/prepared   --output 06_phylogenies/<project_name>

# Phylogenetic inference
python phylogeny.py run   --project 06_phylogenies/<project_name>   --threads 8
```

## Main files to report

Primary biological results:
- `protein_summary.csv`
- `candidates.csv`
- `candidate_report.tsv`
- `candidate_HOGs.tsv`
- `direct_orthologues.tsv`
- `gene_tree.treefile` or a justified `rooted.treefile`

Supplementary/reproducibility files:
- `trimmed.faa`
- `gene_tree.iqtree`
- `branch_evidence.tsv`
- stage-specific manifests and provenance JSON files
- software/environment exports

## Manual review

Manual review is recommended for provisional or unresolved candidates, incomplete or rearranged domain architectures, weak or conflicting homology evidence, repeated-domain proteins, short motif matches lacking structural context, many-to-many orthology relationships, duplication events, and trees with uncertain alignment quality, support, or rooting.

Final family assignment and evolutionary interpretation should be based on convergent evidence rather than any single automated criterion.
