# Release 1.0.0

Integrated source package, organized on 2026-10-03.

- One versioned entry point for detection, external-evidence return, annotation,
  orthology reference-panel calibration/comparison and conservative audit.
- One active motif implementation with its curated workbook and context rules.
- Retained required HMM, curated reference, seed/cache and motif resources.
- Removed previous analysis datasets, predictions, prepared proteomes, QC,
  populated audit workbooks, historical motif implementations and result folders
  from this thesis package. Parent-project data were not modified.
- Dependency specifications and input/panel templates included.

`release_manifest.sha256` inventories the distributed files except itself.
Fibronectin inclusion update: six integration tests passed, including fibronectin screening export/native topology return and fibronectin candidates through custom/expanded reference panels and final audit. Original unresolved detection statuses are preserved. Generic FNIII screening cannot establish family assignment. No analysis outputs were retained.

General target update: configurable one-or-more target proteomes using targets.example.json and --targets-config. Handoff species determine panel target/reference roles; custom orthology extraction no longer assumes Schistosoma. Seven integration tests pass, including two-target raw preparation and single non-schistosome target handoff/expanded-panel audit. Six orthology/phylogeny tests pass; two Biopython-dependent tests remain skipped.

General-target edge cases: eight integration tests and 65 detection tests pass. Explicit FASTA filenames/codes take precedence over species-name hints. Targets without exported candidates remain required in the orthology panel and are not reclassified as reference species.

Preparation documentation update: added ORTHOLOGY_PHYLOGENY_GUIDE.md with a single working directory and explicit raw-reference, panel, orthology and phylogeny paths. Added prepare_reference.py for validated gene-level FASTA preparation from protein-to-gene TSVs. Seven stage-04 tests passed; two Biopython-dependent tests skipped. No biological runs were generated.

Continuous-resume update: default resume now discovers reference FASTAs/annotations under input/reference_raw, prepares references, generates the panel and executes orthology, phylogeny and final audit. Continue reuses existing annotation and supports panel comparisons. Missing duplication evidence is explicitly unknown. Nine integration tests and seven stage-04 tests passed; two Biopython-dependent tests skipped. External OrthoFinder/tree execution was simulated in the orchestration test; real biological inference was not run.

Phylogeny failure handling: IQ-TREE defaults to two independent threads; reversible trimAl IDs prevent truncation; failed-alignment diagnostics preserve membership and per-sequence QC; tool failures record stage, exit code and stderr tail. Group summaries include intentional skips. --retry-failed archives prior failures and preserves completed groups. Partial completion is opt-in and explicitly labelled. Eighteen relevant tests passed; two Biopython-dependent tests skipped. Reported real failure logs were not present in the thesis package; real IQ-TREE exit-2 causes remain unconfirmed.

Checkpoint recovery: added recover command with verified stage reuse, thread overrides, completed-inference reuse after extraction failure, per-group tree retry, archived failed attempts, and validated adoption of older failed runs. Nine integration tests passed, including simulated extraction failure, subsequent tree failure and recovery/adoption without rerunning OrthoFinder. External inference remains simulated in these tests.

Audit independence fix: audit now executes/checkpoints immediately after orthology, before any phylogeny work or tree dependency check. Regression test confirms a completed audited workbook is retained after simulated tree failure. Nine integration tests passed. The user's command_003 log confirms a phylogeny run failure after orthology completion; the separate Test_pipeline_v1 script was patched with a backup.

Unified family phylogeny: optional union of candidate-associated orthogroups per family, with reference members, protein-ID deduplication, family ambiguity/source-group provenance and candidate-ID mapping. Available standalone through --family-map and integrated through --phylogeny-grouping family. Ten stage-04 tests and nine integration tests passed; two Biopython-dependent tests skipped. Only fixture tests ran; the user will test real data locally.
