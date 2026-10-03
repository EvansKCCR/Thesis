# Resumable adhesome workflow

`adhesome_workflow.py` integrates stages 01–05. Use Python 3.11+ and openpyxl
for the final audit. Detection requires HMMER/BLAST; whole-proteome inference
requires OrthoFinder on a suitable Linux/WSL/HPC environment. External predictors
are a deliberate handoff, so an unattended raw-proteome-to-result run is not
guaranteed. All commands use new output directories and record input hashes,
command arrays, logs and manifests. Run commands from this repository root.

```text
01 detection -> supported/provisional hypotheses + external FASTA
  -> 02 external topology/localization -> resume -> 03 motif annotation
  -> 04 complete-proteome panel -> OrthoFinder -> calibration review
       ^ add references, use a new panel/run, compare evidence |
  -> 05 joined evidence -> conservative audited workbook
```

## 1. First output: hypotheses and sequences

Reuse current detection, or omit `--detection` to run detection from the original
proteomes (`--input-dir` can select another input folder):

```bash
python adhesome_workflow.py start --detection 01_detection/results/species_reference_v5 --out workflow_runs/handoff
```

`supported_hypotheses.tsv` and `provisional_hypotheses.tsv` retain every matching
family hypothesis, including competing families for one protein. These labels
are upstream detection evidence tiers, not validated function. `hypotheses.tsv`
retains the original evidence columns. Unresolved/non-family hits are not
promoted; they remain in the original detection CSVs.

Submit `external_candidates.faa` to topology/localization predictors. Its exact
identifiers and full sequences must be preserved. `families/` and
`verified_regions.tsv` are generated for stage 03 with exact sequence checking.
Do not edit files inside the sealed handoff; copy templates elsewhere to fill.

## 2. Return external results and annotate

For the current stage-02 folder:

```bash
python adhesome_workflow.py resume --handoff workflow_runs/handoff --context-dir "02_topology_subcellular localization/adhesome_candidate" --out workflow_runs/annotation
```

The adapter accepts the existing `adhesome_DeepLoc_results.csv`,
`adhesome_SiganP_prediction_results.txt`, `adhesome_NCBI_CDD_hitdata.txt`,
`adhesome_predicted_topologies.3line` filenames, or the canonical names documented
by stage 03. It copies them into the new run using canonical names; originals
are unchanged. A context directory requires all four files. For topology alone:

```bash
python adhesome_workflow.py resume --handoff workflow_runs/handoff --topology-file returned/predicted_topologies.3line --out workflow_runs/annotation_topology
```

Missing/mismatched predictions are reported by `motifs/topology_audit.tsv` and
`motifs/sequence_context.tsv`; they never become positive evidence. Partial
coverage is allowed. `--allow-missing-context` explicitly permits annotation
without returned results. `--evidence` accepts the existing stage-03 site-specific
evidence TSV (family, ID, hash, motif, interval, criterion, status and source).

An optional `--external-evidence returned/evidence.tsv` accepts a filled copy of
`external_evidence_template.tsv`. Every row is checked for unique known IDs,
matching full-sequence SHA-256 and a source for nonempty predictions. These
protein-level summaries feed the final audit only; they do not supply motif
coordinates or satisfy site-specific functional gates. Use native files for
motif context. Leave unavailable evidence blank. Sources are user assertions,
not independently verified experiments. MULocDeep has no native adapter here;
retain its interpretation/source in the returned summary only when appropriate.

## 3. Orthology panel and calibration

```bash
python adhesome_workflow.py calibrate --handoff workflow_runs/handoff --panel workflow_panel.example.json --out workflow_runs/panel8
python adhesome_workflow.py orthology --panel workflow_runs/panel8 --out workflow_runs/orthology8 --threads 16 --analysis-threads 8
```

The JSON panel lists **complete gene-level proteomes**, including all three
targets and one or more references. Paths resolve relative to the JSON file.
Each entry requires a unique species code, `fasta`, `gene_level: true`, and a
documented `source` (prefer assembly/release/download provenance). FASTA headers
must be `SpeciesCode__gene_id`. Collapse raw isoforms using explicit gene maps
before this step; the existing stage-04 preparation supports the bundled
five/eight-species inputs. The wrapper validates IDs and sequences, not biological
proteome completeness. Target candidates must exactly match their panel sequences.
The example deliberately uses the current detection's complete target proteomes.

Inspect `calibration_review.tsv`: absent direct reference evidence and
many-to-many/duplication-related cases identify follow-up candidates. To expand,
copy the JSON, add complete gene-level reference proteomes and run `calibrate`
into a **new** directory. Rerun whole-proteome inference; old HOG IDs cannot be
transferred between panels. Compare by stable candidate ID:

```bash
python adhesome_workflow.py calibrate --handoff workflow_runs/handoff --panel expanded_panel.json --out workflow_runs/panel_expanded
python adhesome_workflow.py orthology --panel workflow_runs/panel_expanded --compare workflow_runs/orthology8 --out workflow_runs/orthology_expanded
```

`calibration_comparison.tsv` compares direct reference presence and review states.
This is reference-panel sensitivity analysis, not fitted confidence or independent
validation. Add biologically justified references and inspect gene trees; more
species do not automatically establish family identity.

For inference performed externally, carry the panel manifest with the run:

```bash
python adhesome_workflow.py orthology --panel workflow_runs/panel8 --results /path/to/Results_DATE --input-manifest /path/to/run/workflow_manifest.json --out workflow_runs/orthology_import
```

The supplied manifest must match the panel. This is an explicit provenance
assertion: the wrapper cannot prove that an arbitrary external results directory
was generated from it. Use the wrapper-run command for directly recorded execution.

Optional per-HOG phylogenies use the existing stage-04 workflow, now with custom
species support. Supply the selected extraction's `candidate_HOGs.tsv`,
`direct_orthologues.tsv`, `candidate_duplications.tsv`, panel `proteomes/`, and
`--species-list workflow_runs/panel8/species.json` to `scripts/phylogeny.py prepare`.
Then use its `run` command (see stage-04 PHYLOGENY.md). If duplication events were
not supplied, resolve that missing input before phylogenetic preparation. Trees
are an explicit review module; they are not silently converted into family scores.

## 4. Join evidence and audit

```bash
python adhesome_workflow.py audit --annotation workflow_runs/annotation --orthology workflow_runs/orthology8 --out workflow_runs/final
```

Outputs: `joined_evidence.tsv`, `audit_input.xlsx`, `audited_candidates.xlsx`,
audit provenance and workflow manifest. The audited workbook includes all original
hypotheses and one collapsed row per species/protein. The existing conservative
family gates remain active. Direct orthologue species counts support arbitrary
added references; HOG membership alone is never counted as direct orthology.
Reference gene/family labels are not fabricated from IDs. Families requiring
specific CDD/homology evidence can therefore remain unassigned.

The final audit requires annotation and orthology to originate from the same
unchanged handoff. Re-run `audit` into a new directory after a calibrated orthology
run. Missing evidence stays unknown. Inspect competing hypotheses, evidence
conflicts and failed architecture gates before interpreting a final assignment.

## Verification

```bash
python -m unittest discover -s workflow_tests -v
python -m unittest discover -s 04_orthology_phylogenetic_classification/scripts -p "test_*.py" -v
```

Failed commands retain logs and lack a completion manifest. Choose a new output
directory after correcting inputs. The successful earlier modules remain reusable.
