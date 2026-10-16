# Default: continuous execution after resume

Follow [the automated-run guide](ORTHOLOGY_PHYLOGENY_GUIDE.md). Put reference FASTAs and matching GFF/TSV files in inputs/reference_raw, with the original Open_Orthobench Input/RefOGs bundle for evolutionary benchmarking. Resume runs annotation, reference preparation, panel construction, OrthoFinder, benchmarking and final audit. Phylogeny is opt-in with --phylogeny and cannot block the audit. The modular commands remain available; --annotation-only deliberately stops after annotation.

# Integrated adhesome workflow — version 1.0

The bundled search panel now covers all 230 crosswalk components. See
[component reference expansion](COMPONENT_REFERENCE_PANEL.md) for panel contents,
BLAST-only refresh of cached detection runs, and optional orthology calibration.

This is a source distribution of integrated workflow for integrin adhesome mapping from proteomes. It contains scripts, tests,
curated HMM/reference/motif resources, configuration templates and dependency specifications. Historical versions in resource names identify their provenance; the integrated
workflow release is 1.0.

## Setup and inputs

Run from `Thesis/integrin_adhesome_analysis` with Python 3.11 or newer. Full execution requires Linux,
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
`../README/targets.example.json` to `targets.json` in the workflow folder and add one entry per target: a unique
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

The active version 1.0 detection library contains all 176 updated Pfam accessions
and reads the 230-row component crosswalk. Component possibilities and reference
similarity are reported separately; crosswalk overlap does not promote family
assignment. See the [detection guide](../integrin_adhesome_analysis/01_detection/README.md).
`start` executes detection and exports its handoff. Topology, localization and
phylogeny are processed after returning evidence through `resume`.

```bash
python adhesome_workflow.py start --targets-config targets.json --input-dir input/targets --out runs/handoff --threads 4
```

Use the actual input folder name. A project may use `inputs/targets`; pass that
path explicitly when it exists. Before starting, append `--check-inputs` to the
command to check every target's FASTA headers and gene-ID mapping without creating
the output directory or running analysis. Detection performs this check before
creating a new run. Each target may use its own `gene_regex`; gene-tagged proteomes
and pipe-delimited validation FASTAs need different mappings.

If an earlier attempt failed before processing and left only a log directory,
choose a new `--out` name for the corrected start command. Keep the failed logs.

To reuse a separately generated detection batch, replace `--input-dir` with
`--detection /path/to/detection_batch`. The first completed module exports:

- `supported_hypotheses.tsv` and `provisional_hypotheses.tsv`: every qualifying
  family hypothesis, including alternatives for the same protein.
- `hypotheses.tsv`: combined evidence table retaining upstream detection fields.
- `component_hypotheses.tsv`: crosswalk matches for exported candidates, when available.
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

## 3. Automatic orthology, evolutionary benchmark and reference calibration

After motif annotation, the default `resume` command automatically reads
`inputs/reference_raw`, prepares references (preserving original benchmark protein units), constructs a panel from
all handoff targets plus those references, runs OrthoFinder, benchmarks reference
orthogroups, and writes the final audit. Optional trees require `--phylogeny` and run after audit. You do not prepare references or write
panel JSON manually. See [the automatic guide](ORTHOLOGY_PHYLOGENY_GUIDE.md) for
accepted annotations, outputs, and commands to expand the panel using `continue`.

Use `--annotation-only` only if you deliberately want to stop after motifs. The
optional low-level `calibrate`, `orthology` and `audit` commands remain documented
in [the manual route](MANUAL_ORTHOLOGY_PHYLOGENY.md).

See [evolutionary benchmarking](EVOLUTIONARY_BENCHMARKING.md) for the original
Input/RefOGs layout, integrity checks and the measured audit-context weight.

## 4. Final conservative audit

The [two-stage family scheme](COMPUTATIONAL_AUDIT.md) refines the same stage 1
protein pool. Complete family architecture plus orthology supports grade A;
experimental annotation and manual identity labels are not prerequisites.

For `resume --out runs/resumed`, the final workbook is
`runs/resumed/downstream/final/Stage_2_output.xlsx`. Its nine sheets follow the
supplied template; `Screening_hypothesis` retains original/recovered hypotheses,
and `Audited_family_assignment` contains one ranked decision per species/protein.
Current motifs, benchmark metrics/weight and actual reference-panel columns are
included. `joined_evidence.tsv`, `audit_input.xlsx` and provenance accompany it.

The annotation manifest records completion of annotation only. The downstream
manifest is written only after all requested downstream stages and the final
audit complete successfully. Failed stages preserve logs and return nonzero.
Use `continue --annotation runs/resumed --out runs/downstream_retry` after
correcting dependencies or input errors; already successful annotation is reused.
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

## Checkpoint recovery

Use python adhesome_workflow.py recover --run runs/resumed --threads 2 --iqtree-threads 2 to reuse verified completed downstream stages and retry unfinished work. Older failed runs can be adopted with --adopt-existing --annotation runs/resumed. See the recovery section in [the guide](ORTHOLOGY_PHYLOGENY_GUIDE.md) for limits and exact paths. Continue creates a new downstream run; recover repairs the existing failed run.

Current final outputs include compact `putative_audited_family`, `protein_name`
and evidence confidence, plus a matching TSV. Eligible calls include positive-control
calibrated assignment probabilities and uncertainty. Hard/background-negative
matches are dropped from final assignments and retained in
`negative_control_exclusions.tsv`. Retained and excluded proteins together account
for the stage 1 pool. See [naming and output policy](COMPUTATIONAL_AUDIT.md) for
calibration commands, exclusions and reference-overlap limits.
