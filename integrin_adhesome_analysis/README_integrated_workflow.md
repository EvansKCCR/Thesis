# Integrin Adhesome Comparative Workflow

## Executable integration of stages 01–05

Use [`adhesome_workflow.py`](adhesome_workflow.py) and the command guide in
[`WORKFLOW.md`](WORKFLOW.md) for the resumable workflow. It exports supported and
provisional detection hypotheses for external topology/localization, imports
returned evidence, runs motif annotation, supports configurable complete-proteome
orthology panels and calibration comparisons, and joins evidence for the final
conservative family audit. [`workflow_panel.example.json`](workflow_panel.example.json)
defines the current eight-species panel. The overview below describes the scientific
modules; its historical five-species examples are superseded by the command guide.

## Overview

This repository implements an auditable comparative-proteomics workflow for identifying, annotating, comparing, and phylogenetically analysing candidate integrin–adhesome proteins. The integrated workflow spans focused Pfam/domain screening, curated reference homology, comparative classification, workbook-driven functional motif annotation, whole-proteome orthology inference, and per-HOG phylogenetic reconstruction.

```text
Prepared proteomes
  -> focused Pfam/HMMER scan
  -> filtered domain evidence
  -> curated reference homology / reciprocal evidence
  -> comparative candidate classification
  -> family-specific candidate libraries
  -> curated motif / functional annotation
       + verified Pfam regions
       + optional DeepTMHMM topology
       + optional site-specific evidence
       + optional STREME/FIMO/TOMTOM
  -> whole-proteome OrthoFinder inference
  -> candidate HOG / direct-orthologue extraction
  -> per-HOG sequence extraction
  -> MAFFT -> trimAl -> IQ-TREE
  -> annotated gene trees
```

## 1. Pfam/domain screening

A focused adhesome Pfam HMM library is built from the active Pfam manifest. Models with curated gathering thresholds and fallback models are retained separately, with hashes and metadata recorded. Prepared proteomes are scanned with HMMER and parsed into retained domain calls and complete audit tables.

Principal outputs:
- `domains.csv`
- `domains.audit.csv`
- `ga.domtblout`
- `fallback.domtblout`
- `scan.json`

## 2. Curated homology and comparative classification

Prepared proteins are compared against curated adhesome references and a broader metazoan reciprocal-reference panel. Forward similarity, reciprocal-family evidence, domain architecture, family rules, and optional profile/phylogeny evidence are integrated into candidate classification.

Principal outputs:
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

`protein_summary.csv` is the preferred one-row-per-protein classification summary; `candidates.csv` retains the detailed evidence matrix.

## 3. Curated functional motif annotation

The current annotation workflow is driven by:
- `Motif_library/adhesion_protein_SLiM_library.xlsx` — sole curated motif-pattern and literature source.
- `Motif_library/context_rules.json` — operational routing, target families, contextual requirements, and evidence gates.

The workbook is read without modification. Required columns, unique motif IDs, evidence tiers, screening recommendations, citation status, confidence scores, and regex syntax are validated. Formula cells and Excel errors are rejected.

Candidate motif matches are assigned transparent non-probabilistic tiers:
- `high_confidence_hypothesis`
- `context_supported_candidate`
- `sequence_only_candidate`
- `review_only`
- `context_conflict`

A high-confidence hypothesis requires an eligible curated entry, supported family architecture, and all configured contextual/evidence requirements to be supported. Workbook confidence scores describe the library rule and are not candidate probabilities.

### Verified domain/context input

`import_domain_regions.py` imports retained upstream Pfam intervals only when the candidate identifier and full amino-acid sequence exactly match the upstream search FASTA.

### Optional topology input

DeepTMHMM `predicted_topologies.3line` files or ZIP archives can be imported. Sequence IDs and full sequences must match exactly. A cytoplasmic tail is derived only from a single predicted transmembrane helix immediately followed by terminal inside residues. Predicted inside/outside states are treated as provisional topology evidence.

### Optional independent site evidence

A TSV supplied with `--evidence` must contain:

```text
family sequence_id sequence_sha256 motif_id start end criterion status source
```

Evidence is accepted only when identifiers, sequence hashes, coordinates, motif IDs, and criteria match the current run. Contradictory evidence takes precedence.

### Current annotation outputs

- `motif_annotations.tsv`
- `high_confidence_hypotheses.tsv`
- `functional_hypotheses.tsv`
- `reader_hypotheses.tsv`
- `inspection_states.tsv`
- `library_audit.tsv`
- `library_coverage.tsv`
- `motif_prevalence.tsv`
- `domain_inspection.tsv`
- `topology_audit.tsv`
- `used_regions.tsv`
- `discovered_words.tsv`
- `discovered_sites.tsv`
- `discovery_annotation.tsv`
- `curated_library.snapshot.json`
- `context_rules.snapshot.json`
- `curated_summary.json`
- `CURATED_FINDINGS.md`
- `manifest.json`
- `REPORT.md`

Motif matches remain conditional functional hypotheses; short-word recurrence, generic domain support, topology prediction, or workbook confidence alone does not establish function.

## 4. Optional exploratory/PWM motif analysis

Exploratory exact-word discovery is retained for lead generation. By default, 4–6 residue words are considered, with at least two supporting sequences and at least 50% within-family prevalence.

Optional MEME Suite stages:
- STREME — de novo protein motif discovery
- FIMO — motif localization
- TOMTOM — comparison with a genuine protein PWM database

External stages default to `skip`. The XLSX motif library contains regex rules, not PWMs.

## 5. Whole-proteome orthology inference

Orthology is inferred across complete gene-level proteomes from *Homo sapiens*, *Caenorhabditis elegans*, *Schistosoma haematobium*, *S. japonicum*, and *S. mansoni*. Isoforms are collapsed to one representative per explicit gene identifier before OrthoFinder analysis. Candidate filtering is performed only after whole-proteome inference.

Principal outputs:
- `candidate_HOGs.tsv`
- `direct_orthologues.tsv`
- `candidate_report.tsv`
- `candidate_duplications.tsv`
- `provenance.json`

## 6. Per-HOG phylogenetic reconstruction

Each candidate-associated HOG is analysed independently. All HOG members are extracted, aligned with MAFFT, trimmed with trimAl, quality checked, and analysed with IQ-TREE.

Default settings:
- MAFFT L-INS-i: `--localpair --maxiterate 1000`
- trimAl: `-automated1`
- IQ-TREE ModelFinder: `-m MFP`
- 1,000 ultrafast bootstrap replicates
- 1,000 SH-aLRT replicates
- BNNI optimization
- fixed reproducible seed

Per-HOG outputs:
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

Trees are unrooted by default; rooting should only use a biologically justified outgroup.

## Core scripts

- `build_adhesome_pfam_library.py`
- `parse_adhesome_domtblout.py`
- `adhesome_pipeline.py`
- `adhesome_classification.py`
- `adhesome_core.py`
- `adhesome_motifs.py`
- `build_reference_sets.py`
- `family_motifs.py`
- `curated_annotation.py`
- `import_domain_regions.py`
- `refseq_gene_map.py`
- `pipeline.py`
- `phylogeny.py`

## Main reference/configuration resources

- active adhesome Pfam manifest
- family classification rules
- curated adhesome reference FASTAs/catalogues
- metazoan reciprocal-reference set
- phylogenetic anchors/reference locks
- `Motif_library/adhesion_protein_SLiM_library.xlsx`
- `Motif_library/context_rules.json`
- optional site-specific evidence TSV
- optional DeepTMHMM topology files/ZIP archives
- optional protein MEME/PWM database
- five prepared gene-level species proteomes

## Reproducibility

The workflow retains SHA-256 hashes, exact command arrays, software paths and versions, raw search outputs, filtered evidence tables, candidate mapping audits, workbook/rule snapshots, topology/evidence audits, alignment/tree QC, and stage-specific manifests. Existing result directories are generally protected from silent overwrite.

## Interpretation guidance

- Pfam/domain architecture supports structural compatibility.
- Sequence similarity and reciprocal evidence support homology.
- Curated motif matches support functional hypotheses only when required context is satisfied.
- Workbook confidence scores are properties of library entries, not candidate probabilities.
- Missing evidence remains unknown.
- Shared HOG membership does not by itself establish direct orthology.
- Pairwise orthology and duplication evidence should be considered together with gene-tree topology.
- Phylogenetic interpretation should consider homologous coverage, alignment quality, rooting, branch support, and gene-tree/species-tree reconciliation.

## Suggested analysis sequence

```bash
# Pfam/reference resources
python build_adhesome_pfam_library.py ...
python build_reference_sets.py ...

# Scan and classify
python adhesome_pipeline.py scan ...
python adhesome_pipeline.py homology ...
python adhesome_pipeline.py classify ...

# Import verified domain regions
python import_domain_regions.py   --input-dir family_specific_candidate_library   --domains <retained_domains.csv>   --source-fasta <exact_upstream_search.faa>   --out results/verified_domain_regions.tsv

# Curated annotation
python family_motifs.py   --out results/family_inspection_curated   --regions results/verified_domain_regions.tsv   --topology-dir <deeptmhmm_directory>

# Whole-proteome orthology
python pipeline.py prepare --human-map human_protein_gene.tsv
python pipeline.py run --name combined_01 --threads 16 --analysis-threads 8

# Candidate orthology extraction
python pipeline.py extract   --results <completed_orthofinder_results>   --candidate-set <candidate_set>   --output <candidate_output>

# Per-HOG phylogeny
python phylogeny.py prepare   --hogs <candidate_output>/candidate_HOGs.tsv   --orthologues <candidate_output>/direct_orthologues.tsv   --duplications <candidate_output>/candidate_duplications.tsv   --proteomes 02_primary_proteomes/prepared   --output 06_phylogenies/<project_name>

python phylogeny.py run   --project 06_phylogenies/<project_name>   --threads 8
```

## Main files to report

Primary biological outputs:
- `protein_summary.csv`
- `candidates.csv`
- `functional_hypotheses.tsv`
- `high_confidence_hypotheses.tsv`
- `candidate_report.tsv`
- `candidate_HOGs.tsv`
- `direct_orthologues.tsv`
- `gene_tree.treefile` or justified `rooted.treefile`

Recommended supplementary/reproducibility outputs:
- `motif_annotations.tsv`
- `library_audit.tsv`
- `topology_audit.tsv`
- `used_regions.tsv`
- `trimmed.faa`
- `gene_tree.iqtree`
- `branch_evidence.tsv`
- stage-specific manifests/provenance JSON files
- environment exports

## Manual review

Manual review is recommended for unresolved candidate classifications, incomplete or rearranged domain architectures, weak/conflicting homology evidence, motif hypotheses with missing context requirements, topology conflicts, repeated-domain proteins, many-to-many orthology relationships, duplication events, and trees with uncertain alignment quality, branch support, or rooting.

Final family assignment, functional interpretation, and evolutionary interpretation should rely on convergent evidence rather than any single automated result.
