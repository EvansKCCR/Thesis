# Integrated adhesome workflow — version 1.0.0

This is a source distribution of integrated workflow for integrin adhesome mapping from proteomes. It contains scripts, tests,
curated HMM/reference/motif resources, configuration templates and dependency
specifications. Historical versions in resource names identify their provenance; the integrated
workflow release is 1.0.0.

## Setup and inputs

Run from this directory with Python 3.11 or newer. Full execution requires Linux,
WSL or an HPC environment with the external programs available on PATH:

```bash
conda env create -f environment.yml
conda activate adhesome-workflow-v1
python adhesome_workflow.py --version
```

The environment specifies Python, HMMER, BLAST+, OrthoFinder, MAFFT, trimAl,
IQ-TREE, Biopython and openpyxl. The environment specification has not been
installed/solver-validated as part of packaging. For annotation/audit only,
`python -m pip install -r requirements.txt` supplies Python dependencies;
this does not install the external programs. MEME Suite is optional and disabled
by default. External topology/localization prediction services are not bundled.

Supply one or more complete target proteomes under an input directory. Copy
`targets.example.json` to `targets.json` and add one entry per target: a unique
`code`, descriptive `species`, `assembly`, source URL, FASTA `path` and
`gene_regex` whose first capture group extracts an explicit gene ID. Codes may
contain letters, digits and single underscores, must start with a letter, and
must not contain double underscores. There is no fixed species set or count.

For `start --input-dir`, the batch runner discovers configured FASTA filenames
within that directory (or matches the species code/epithet); exactly one file is
required per configured target. Prepared IDs become `Code__gene_id`. The example
uses `TargetA.faa` and headers containing `gene=...` or `gene:...`. Adjust the
regex to your annotation format; do not guess genes by removing accession
versions. Preparation collapses isoforms. Candidate-only FASTAs must not be used
as complete proteomes. No previous result is required.

The bundled detection reference/HMM/motif resources remain biologically curated
resources; accepting arbitrary target species does not establish their taxonomic
coverage or validate cross-species functional transfer.

## 1. Detection and external handoff

```bash
python adhesome_workflow.py start --targets-config targets.json --input-dir inputs/targets --out runs/handoff --threads 8
```

To reuse a separately generated detection batch, replace `--input-dir` with
`--detection /path/to/detection_batch`. The first completed module exports:

- `supported_hypotheses.tsv` and `provisional_hypotheses.tsv`: every qualifying
  family hypothesis, including alternatives for the same protein.
- `hypotheses.tsv`: combined evidence table retaining upstream detection fields.
- `external_candidates.faa`: one full sequence per candidate protein.
- `external_evidence_template.tsv`: sequence-hash-bound return template.
- `families/`, `verified_regions.tsv`: sequence-checked stage-03 inputs.

Supported/provisional are evidence tiers, not experimentally validated function.
Unresolved detection rows remain in the underlying detection CSVs. Submit the
external FASTA to the desired topology/localization tools, preserving exact IDs
and full sequences. Copy the return template elsewhere before filling it: files
inside a sealed run must not be edited.

## 2. Return topology/localization and annotate motifs

For a DeepTMHMM 3line file:

```bash
python adhesome_workflow.py resume --handoff runs/handoff --topology-file returned/predicted_topologies.3line --out runs/annotation
```

For a combined predictor folder, use `--context-dir returned/context` instead.
That folder must contain these files (either name in each pair is accepted):

| Predictor | Canonical filename | Alternate filename |
|---|---|---|
| DeepTMHMM | deeptmhmm_predicted_topologies.3line | adhesome_predicted_topologies.3line |
| DeepLoc | DeepLoc_adhesome_candidate_list.csv | adhesome_DeepLoc_results.csv |
| SignalP | SignalP_adhesome_candidate_list.txt | adhesome_SiganP_prediction_results.txt |
| CDD | NCBI_CDD_adhesome_hitdata.txt | adhesome_NCBI_CDD_hitdata.txt |

DeepLoc, SignalP and CDD evidence requires an exact full-sequence match in the
associated topology input. Partial coverage is allowed and recorded in
`motifs/topology_audit.tsv` and `motifs/sequence_context.tsv`. Unknown/mismatched
records never create positive evidence. `--allow-missing-context` explicitly
permits annotation without returned predictions.

Optional `--external-evidence returned/evidence.tsv` accepts the filled template
for protein-level final-audit summaries only; it does not transfer motif positions
or satisfy site-specific motif gates. Known IDs, exact sequence hashes and a
nonempty source for nonempty evidence are required. Native files provide motif
context. Optional `--evidence` accepts the motif module's site-specific TSV:
`family sequence_id sequence_sha256 motif_id start end criterion status source`.
Use tabs, one-based inclusive intervals and supported/contradicted statuses.
Predictions and source assertions remain distinct from experimental validation.

## 3. Whole-proteome orthology and calibration

Edit `panel.example.json` to point at complete gene-level reference proteomes and
record their assembly/release, source and preparation provenance. Paths are
relative to the JSON file. Target entries point at this run's full detection
proteomes. Reference FASTAs must have unique `SpeciesCode__gene_id` headers and
one representative per explicit gene ID. Prepare/collapse raw reference isoforms
using a gene mapping before this step; never infer genes by trimming accession
versions. The existing stage-04 preparation scripts support the established
five/eight-species profiles and explicit human/GFF gene maps.

```bash
python adhesome_workflow.py calibrate --handoff runs/handoff --panel panel.example.json --out runs/panel_initial
python adhesome_workflow.py orthology --panel runs/panel_initial --out runs/orthology_initial --threads 16 --analysis-threads 8
```

The panel validates full target candidate-sequence identity and records hashes.
`gene_level: true` and `source` are required declarations; they do not prove
biological proteome completeness. OrthoFinder runs on the whole proteomes and
only then extracts candidate HOGs and direct orthologues.

Review `calibration_review.tsv` for absent reference evidence and ambiguous
relationships. To add references, copy the panel JSON, append complete reference
proteomes with new unique species codes, and create a new panel and inference:

```bash
python adhesome_workflow.py calibrate --handoff runs/handoff --panel panel.expanded.json --out runs/panel_expanded
python adhesome_workflow.py orthology --panel runs/panel_expanded --compare runs/orthology_initial --out runs/orthology_expanded
```

`calibration_comparison.tsv` compares reference support by stable candidate ID.
HOG IDs are not transferable between inference runs. This is panel sensitivity
analysis, not confidence fitting or independent validation. Both comparisons must
use the same sequence-verified handoff. Existing runs are preserved.

For inference executed elsewhere, use `--results /path/to/Results_DATE` and
`--input-manifest /path/to/the/panel/workflow_manifest.json` on the orthology
command. The manifest must match the selected panel; its association with an
external results directory is a user provenance assertion. Wrapper-executed runs
record the exact invocation directly.

Per-HOG phylogenies remain an optional review module. Use stage 04's
`scripts/phylogeny.py prepare` with the extracted HOG, direct-orthologue and
actual duplication-event TSVs, panel `proteomes/`, and
`--species-list runs/panel_expanded/species.json`; then use `phylogeny.py run`.
See its `--help` for required options. Do not fabricate missing duplication
records. Trees are not automatically converted into family-support scores.

## 4. Final conservative audit

```bash
python adhesome_workflow.py audit --annotation runs/annotation --orthology runs/orthology_expanded --out runs/final
```

The final module joins original detection hypotheses, native/returned context,
motif annotations and direct orthology evidence into `joined_evidence.tsv` and
`audit_input.xlsx`, then creates `audited_candidates.xlsx`. The workbook retains
all hypotheses plus collapsed species/protein rows. Architecture gates remain
mandatory; HOG membership alone is not direct orthology, and reference gene/family
labels are not fabricated. New reference species count in direct support.

All commands require new output directories. Manifests record workflow version,
input/output hashes and module provenance; failed commands retain diagnostic logs
without a completion manifest. Resume from the last completed module into a new
output directory after correcting inputs. The final audit rejects mismatched or
changed handoffs/panels. Re-run only orthology and audit when expanding a panel.

## Tests

```bash
python -B -m unittest discover -s workflow_tests -v
cd 01_detection
python -B -m unittest discover -s tests -v
cd ../03_functional_motif_annotation
python -B -m unittest discover -s test -v
cd ../04_orthology_phylogenetic_classification/scripts
python -B -m unittest discover -p "test_*.py" -v
```

Tests use temporary fixtures and do not require retained analysis outputs. Full
biological inference additionally requires the stated tools and user inputs.
## Fibronectin-like inclusion

Fibronectin-like candidates are included in the combined hypothesis tables,
external FASTA/return template, motif annotation, orthology panel mapping and
reference-expansion comparisons. Additional subset exports are
`fibronectin_like_hypotheses.tsv` and `fibronectin_like_candidates.faa`.

Detection labels `fibronectin_like`, `fibronectin_III_ligand` and
`fibronectin_like_candidate` route to motif family `fibronectin_like_candidate`.
Unresolved fibronectin screening rows are exported as provisional follow-up
hypotheses with original `status` and `family` preserved and an explicit
`export_basis`. This is not promotion to supported family identity. Other
unresolved rows remain excluded. Generic FNIII screening alone cannot pass the
final family-assignment gate. Competing family hypotheses remain distinct while
sequences are deduplicated. Returned evidence uses the same exact sequence/hash
checks. No separate orthology reference panel is needed.
All handoff target species, including targets with no exported hypotheses, must be represented in each orthology panel, with at least one additional reference species. Target membership is recorded in the handoff manifest and panel targets.json. With --detection, targets are inferred from the supplied detection batch; --targets-config can additionally validate the expected species set.
