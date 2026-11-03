# Deployment validation — version 1.0

Structural readiness: **pass**. Deployment: **warning**.
Execution ready: **False**. No analyses were launched.

| Check | Status | Detail |
| --- | --- | --- |
| Workflow configuration JSON | pass | 19 configuration files parsed |
| Workflow configuration | pass | Isolated settings, Stage 2 profile and paired family coverage |
| Family and identity rule initialization | pass | {'identity_rules': 197, 'stage2_additions': 11, 'family_rules': 222} |
| Configuration resource: config/annotation/family_rules.json | pass | Required static library configuration |
| Configuration resource: config/annotation/family_identity_rules.json | pass | Required static library configuration |
| Configuration resource: config/annotation/family_set_registry.json | pass | Required static library configuration |
| Configuration resource: config/motif_annotation/adhesion_core.json | pass | Required static library configuration |
| Configuration resource: config/motif_annotation/motif_master.json | pass | Required static library configuration |
| Configuration resource: config/orthology_phylogeny/benchmark_policy.json | pass | Required static library configuration |
| Source-layout paths | pass | 589 required source paths resolve; historical provenance caches are optional |
| Deployment source inventory | pass | Source-only release checksums resolve; user inputs and results excluded |
| Scientific asset checkout integrity | pass | Git preserves hashed scientific library/config bytes; release code-text checks normalize LF without changing biological fingerprints |
| Script syntax | pass | 110 Python files parsed without importing scientific datasets |
| Stage entry point: detection | pass | integrin_adhesome_analysis/resource_library/detection/scripts/run_detection.py |
| Stage entry point: motifs | pass | integrin_adhesome_analysis/resource_library/motif_annotation/family_motifs.py |
| Stage entry point: references | pass | integrin_adhesome_analysis/resource_library/orthology_phylogeny/scripts/prepare_references.py |
| Stage entry point: orthology | pass | integrin_adhesome_analysis/resource_library/orthology_phylogeny/scripts/pipeline.py |
| Stage entry point: audit | pass | integrin_adhesome_analysis/resource_library/annotation/integrated_adhesome_pipeline.py |
| Stage entry point: isolated_branches | pass | integrin_adhesome_analysis/main/run_pipeline.py |
| Stage entry point: homology | pass | homology_retrieval/reciprocal_blast.py |
| Stage entry point: reconciliation | pass | reconciliation/reconcile_candidates.py |
| Stage entry point: independent_phylogeny | pass | independent_family_phylogeny/run_phylogeny.py |
| Detection model library | pass | {'required': 214, 'available': 214, 'missing': []} |
| Current detection-rule model coverage | pass | {'required': 172, 'missing': []} |
| Deployed HMM checksums | pass | Combined, manifest and split-library checksums; no library rebuild performed |
| Output destinations | pass | Separate gate/downstream namespaces and independent bundle outputs; no directories created |
| Stage dependency graph | pass | 32 conditional stages topologically ordered; benchmark is a soft audit dependency; phylogeny follows audit |
| Raw-input paths | pass | Homology aliases resolve to main-project raw roots |
| Missing input handling | pass | unidentified_input decisions recorded per requesting stage; mandatory stops cannot be disabled |
| Analysis-gate isolation | pass | Separate gate namespaces, raw roots, branch contracts, reserved identifiers and existing overlap checks |
| Stale output dependencies | pass | All generated dependencies are planned outputs; no saved audit, result, panel or inference checkpoint is opened |
| Selective reference loading | pass | Homology resolves source species before loading sequences; explicit reference-set/proteome/panel selections are supported |
| Software dependencies | pass | [{'name': 'hmmsearch', 'available': True, 'path': '/home/evans_adu/miniconda3/envs/adhesome-workflow-v1/bin/hmmsearch', 'version_verified': False, 'required': True}, {'name': 'hmmpress', 'available': True, 'path': '/home/evans_adu/miniconda3/envs/adhesome-workflow-v1/bin/hmmpress', 'version_verified': False, 'required': True}, {'name': 'blastp', 'available': True, 'path': '/home/evans_adu/miniconda3/envs/adhesome-workflow-v1/bin/blastp', 'version_verified': False, 'required': True}, {'name': 'makeblastdb', 'available': True, 'path': '/home/evans_adu/miniconda3/envs/adhesome-workflow-v1/bin/makeblastdb', 'version_verified': False, 'required': True}, {'name': 'orthofinder', 'available': True, 'path': '/home/evans_adu/miniconda3/envs/adhesome-workflow-v1/bin/orthofinder', 'version_verified': False, 'required': True}, {'name': 'mafft', 'available': True, 'path': '/home/evans_adu/miniconda3/envs/adhesome-workflow-v1/bin/mafft', 'version_verified': False, 'required': True}, {'name': 'trimal', 'available': True, 'path': '/home/evans_adu/miniconda3/envs/adhesome-workflow-v1/bin/trimal', 'version_verified': False, 'required': True}, {'name': 'iqtree', 'available': True, 'path': '/home/evans_adu/miniconda3/envs/adhesome-workflow-v1/bin/iqtree', 'version_verified': False, 'required': True}, {'name': 'openpyxl', 'available': True, 'version': '3.1.5', 'required': True}, {'name': 'biopython', 'available': True, 'version': '1.85', 'required': True}] |
| Clean-state dry run | pass | Configuration and transition simulation complete; no scientific process launched and no biological results fabricated |
| Deployment readiness | warning | Structural validation and data/software execution readiness are reported separately |

## Input decisions

| Stage | Input | Mandatory | Decision | Reason |
| --- | --- | --- | --- | --- |
| target.detection | target_raw | True | stop_stage | expected_absent_in_clean_state |
| target.detection | target_configuration | True | stop_stage | expected_absent_in_clean_state |
| validation.detection | validation_fasta | True | stop_stage | expected_absent_in_clean_state |
| validation.detection | validation_labels | True | stop_stage | expected_absent_in_clean_state |
| calibration.detection | calibration_fasta | False | skip_operation | expected_absent_in_clean_state |
| calibration.detection | calibration_labels | False | skip_operation | expected_absent_in_clean_state |
| homology.rbh | curated_references | True | stop_stage | expected_absent_in_clean_state |
| homology.rbh | target_raw | True | stop_stage | expected_absent_in_clean_state |
| homology.rbh | reference_raw | True | stop_stage | expected_absent_in_clean_state |
| target.references | reference_raw | True | stop_stage | expected_absent_in_clean_state |
| validation.references | reference_raw | True | stop_stage | expected_absent_in_clean_state |
| target.benchmark | benchmark_refogs | False | skip_operation | expected_absent_in_clean_state |
| target.benchmark | benchmark_input | False | skip_operation | expected_absent_in_clean_state |
| validation.benchmark | benchmark_refogs | False | skip_operation | expected_absent_in_clean_state |
| validation.benchmark | benchmark_input | False | skip_operation | expected_absent_in_clean_state |
| target.audit | negative_exclusion_fasta | True | stop_stage | expected_absent_in_clean_state |
| target.audit | negative_exclusion_labels | True | stop_stage | expected_absent_in_clean_state |
| validation.audit | negative_exclusion_fasta | True | stop_stage | expected_absent_in_clean_state |
| validation.audit | negative_exclusion_labels | True | stop_stage | expected_absent_in_clean_state |

## Planned conditional transitions

| Stage | Prerequisites | Simulated state | Executed |
| --- | --- | --- | --- |
| target.detection | none | blocked_on_data | false |
| validation.detection | none | blocked_on_data | false |
| calibration.detection | none | skipped_optional_input | false |
| homology.rbh | none | blocked_on_data | false |
| target.external_context | target.detection | blocked_on_upstream | false |
| validation.external_context | validation.detection | blocked_on_upstream | false |
| calibration.external_context | calibration.detection | skipped_optional_gate | false |
| target.motifs | target.external_context | blocked_on_upstream | false |
| validation.motifs | validation.external_context | blocked_on_upstream | false |
| calibration.motifs | calibration.external_context | skipped_optional_gate | false |
| target.references | target.motifs | blocked_on_upstream | false |
| validation.references | validation.motifs | blocked_on_upstream | false |
| calibration.references | calibration.motifs | skipped_optional_gate | false |
| target.orthology | target.references | blocked_on_upstream | false |
| validation.orthology | validation.references | blocked_on_upstream | false |
| calibration.orthology | calibration.references | skipped_optional_gate | false |
| target.benchmark | target.orthology | skipped_optional_input | false |
| validation.benchmark | validation.orthology | skipped_optional_input | false |
| calibration.benchmark | calibration.orthology | skipped_optional_gate | false |
| target.audit | target.orthology | blocked_on_upstream | false |
| validation.audit | validation.orthology | blocked_on_upstream | false |
| calibration.audit | calibration.orthology | skipped_optional_gate | false |
| target.report | target.audit | blocked_on_upstream | false |
| target.phylogeny | target.audit | blocked_on_upstream | false |
| validation.report | validation.audit | blocked_on_upstream | false |
| validation.phylogeny | validation.audit | blocked_on_upstream | false |
| calibration.report | calibration.audit | skipped_optional_gate | false |
| calibration.phylogeny | calibration.audit | skipped_optional_gate | false |
| reconciliation.library | target.audit, homology.rbh | blocked_on_upstream | false |
| independent.align | reconciliation.library | blocked_on_upstream | false |
| independent.trim | independent.align | blocked_on_upstream | false |
| independent.tree | independent.trim | blocked_on_upstream | false |

## Warnings and limits

- HMM metadata retains its historical build-time family-rule fingerprint; current positive/screening model coverage is checked separately, without rewriting that provenance.
- Executable paths are checked without running programs; tool versions must match the supplied conda environment before real execution.
- OrthoFinder bundled/transitive programs are not executed by this check; install the supplied conda environment and verify its configured search/tree backends before inference.
- Historical standalone families.json/queue OG identifiers remain source-run-specific; a fresh source result requires reviewed mappings or the reconciliation handoff.
- Biological inputs are absent or assumed absent: expected for a clean-state structural dry run; actual execution is not ready.
- A structural pass is not evidence that biological analyses, inference, probability calibration or performance thresholds will succeed.
