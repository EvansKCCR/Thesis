
Version 1.0 library and stage-boundary update (2026-10-14): rebuilt the active library as adhesome_v1.0.hmm with all 174 updated manifest accessions downloaded from the documented EMBL-EBI InterPro API. Model versions, URLs and hashes are retained in source_v1.0/download_provenance.json. The 230-component crosswalk is validated and used to annotate domain overlap and reference-similarity hypotheses without promoting detection status. Selected candidate crosswalk matches are preserved in the external handoff. Detection accepts domain/homology/profile/motif evidence and defers topology, localization and phylogeny to resume. Active runners, resources, version flags and new pipeline metadata report 1.0. The former run_domain_v5 entry point redirects to version 1.0; historical resource/source identifiers remain provenance. The resume default now matches input/reference_raw. All 67 detection, 9 integration and 11 runnable stage-04 fixture tests passed; two Biopython-dependent tests skipped. The handoff crosswalk/return fixture passed after its additional assertions. Twelve CLI version checks passed. No target biological inference was executed.

Current-resource consolidation (2026-10-14): the latest 96-entry family catalogue supersedes complete entries from the 26/32-entry catalogues and is the single canonical config/family_rules.json used by detection and the new reference-gated audit. Removed superseded rule JSONs, their legacy builder, the duplicate audit implementation and the old detection motif library; source hashes remain in resource_resolution.json. Added the two required Pfam gates PF00046/PF03535 and rebuilt the 176-model library; all 96 old model accessions are covered, so obsolete HMM copies and targeted sources were removed. Updated the motif adapter for the supplied evidence-based JSON/XLSX schema, 53 adhesion rules by default and optional 347-rule master scope. Current definitions retain their exclusions and no automatic context confidence is assigned. Resume accepts all current family entries and independently returned identity evidence. Audit preserves Resolved_candidates and review states; stale motif IDs cannot contribute scores. Family preparation reads the new audit schema, preserving broad-only C candidates with explicitly unassigned family identity. Failed-run recovery refreshes changed audits and dependent family preparation while reusing orthology. Canonical pipeline version remains 1.0; resource schema versions remain separate. 112 fixture tests passed; two Biopython-dependent tests skipped. Resource hashes/model coverage/motif schema checks and active CLI version smoke checks passed. No target biological inference or empirical accuracy benchmark was run. Migration commands are in CURRENT_RESOURCES.md.

Input preflight fix: the supplied command_001 log stopped before sequence processing because inputs/targets did not exist. Added a pre-output input check covering discovery, every FASTA header, duplicate IDs and per-target gene mappings, plus --check-inputs on start and detection entry points. Corrected the main guide's target-folder example. A separate mapping is required for the provided pipe-delimited validation panel; targets.validation.corrected.json preserves the proteome mappings and uses accession IDs for all 494 validation controls. The actual validation folder was renamed to inputs/targets during diagnosis; a read-only check passed all four inputs (10,935/14,700/16,936/494 records). All 68 detection and 10 workflow fixture tests passed; a final whitespace-equivalence check passed. No biological inference ran and the external validation installation was not modified. INPUT_PREFLIGHT.md lists the corrected command and updated files.

Validation-target retry fix (2026-10-15): confirmed the deployed targets.json still contained the gene-tag mapping for pipe-delimited validation headers. Added --species to the detection runner/batch stage for explicit selected-target processing while retaining other targets and batch status. Recovery verifies unselected configuration, raw/prepared FASTA hashes, scan parameters/library and classification inputs/resources. Failed selected target directories are archived under archived_attempts. Handoff export rejects incomplete batch statuses. Sixty-nine detection and ten workflow fixture tests passed, including byte-preserved completed outputs and refusal to reuse changed inputs. The actual Sman/Shae/Sjap completed input/provenance hashes were checked read-only. After approval, backed up and corrected the validation installation's targets.json and installed the three updated scripts; its input-only check passed all 494 validation records and all other target headers. Biological inference was not launched. RETRY_VALIDATION.md gives the validation-only retry and combined handoff-export commands.

Evidence reporting fix (2026-10-15): distinguish exact-family reference gaps from completed searches without qualifying hits; retain reciprocal not-tested reasons when no forward hit qualified. Expose parent/legacy-family evidence as context only, with no subfamily promotion. Single-domain/repeat-only architectures report order not_applicable; missing multidomain constraints remain not_configured with a reason. Added per-rule reference coverage tables and combined handoff export. A cached-evidence refresh utility creates a new detection run while preserving the original and validating FASTA/domain/homology hashes. New searches checksum their BLAST artifacts. Seventy-five detection and ten workflow fixture tests passed. All four real target/search snapshots were checked read-only. Installed five scripts and the guide in the validation project with backup evidence_reporting_backup_20261015_053207; no real analysis was launched or results modified. EVIDENCE_STATUS.md explains interpretation and commands. Pipeline version remains 1.0.

Crosswalk reference expansion (2026-10-15): retrieved and validated reviewed human sequences for all 230 crosswalk components. Retained the original 71 forward/142 reciprocal entries and added 230 component references plus 53 explicitly labelled family records, yielding 354 forward/425 reciprocal records and 48 represented family labels. Gene/taxon checks resolve ITGA7 to Q13683 and RAC1 to P63000; reject mismatched alternatives in DOCK1/PRKCA rows while preserving original pins in coverage/provenance. Paired searches now compare tied records by accession to handle hierarchical labels without false conflicts. Component-only hits remain provisional follow-up hypotheses and pass handoff/external return without inventing family identity; heterogeneous component buckets are excluded from pooled family motif discovery. Added explicit BLAST-only refresh with cached prepared/domain files. Exported a 242-accession targeted human orthology FASTA with generated scope/hash sidecar; automatic preparation preserves that scope. Resource validation checks full component coverage and paired-panel hashes/membership. Eighty detection, ten workflow and nineteen motif fixture tests passed. No real target detection, orthology or phylogeny was run. Pipeline version remains 1.0. See COMPONENT_REFERENCE_PANEL.md.

The expanded panel and seventeen script/config/test files were installed in the validation project after approval, with prior scripts and paired panel preserved in component_panel_backup_20261015_060720. The deployed resource validator confirms 230 component references, 354 forward records and 425 reciprocal records. Existing raw inputs, detection outputs and handoff snapshots were not modified.

Evolutionary benchmark and optional phylogeny (2026-10-15): adapted the supplied Open_Orthobench scorer with a workflow mode and standard-library adapter. The default inputs/reference_raw location accepts the original BENCHMARKS/Input and RefOGs bundle, with automatic discovery and preparation of original protein units; no invented gene collapse or hand-written JSON is required. Original-cohort integrity checks require the supplied scorer's 12 proteomes, 70 RefOGs, 1,945 memberships and 251,378 input proteins; custom curated cohorts require an explicit protocol flag. Scoring restores exact original IDs with sequence checks, uses full root orthogroups, projects to the benchmark cohort and retains uncertain-assignment handling. The measured reference-cohort F1 weights orthogroup context in the audit, with raw/weighted scores and provenance retained; it cannot create resolving identity or override architecture gates. Missing/invalid benchmarking remains explicit and nonfatal. Phylogeny is now disabled by default, enabled only by --phylogeny, and runs after audit; tree failures preserve audit/workflow completion and can be explicitly retried through recovery without rerunning completed orthology/audit. Version remains 1.0. Twenty-four stage-04 fixtures (two skipped without Biopython), twelve workflow fixtures and eight audit fixtures passed: 42 passed, two skipped. No empirical Open_Orthobench run or target inference was executed because the original dataset is not staged in this source package. See EVOLUTIONARY_BENCHMARKING.md.

Benchmark execution failures can also be recovered after an audit completed on the neutral baseline; completed orthology is reused and the dependent audit is refreshed. --retry-benchmark provides an explicit recomputation option for sealed reports, preserving prior completion/report outputs. The focused end-to-end fixture verifies both failure and successful recovery without another OrthoFinder invocation.

Benchmark reference import fix (2026-10-16): the validation run stopped before OrthoFinder because original benchmark peptide ENSDARP00000141840 contains internal stop markers. Benchmark-only preparation now removes one terminal stop and encodes remaining stops as X without deleting proteins or internal positions. Per-protein sequence_qc.tsv records positions, lengths and sequence hashes; preparation/source metadata retain the policy and counts. Matching uses the same explicit normalization, rejects any other sequence difference and hashes the normalization implementation and QC reports. Benchmark reports identify affected RefOG members and state that metrics describe the processed cohort, not unmodified-sequence performance. Ordinary references remain strict. All 14 staged reference proteomes imported successfully in temporary outputs (284,844 input proteins, 280,599 prepared units); all 251,378 benchmark proteins matched exactly under the declared policy. The 864 internal markers in 173 proteins were recorded, including two RefOG members; all 95 checked raw FASTA/RefOG file hashes remained unchanged. Stage-04 regression tests: 25 passed, two skipped without Biopython. All 12 workflow tests passed, including recovery from partial reference preparation with completed annotation preserved. No real OrthoFinder, benchmark accuracy inference or phylogeny was launched. Version remains 1.0; recover the failed validation run with recover --run runs/annotation/downstream --threads 4. See EVOLUTIONARY_BENCHMARKING.md.

Computational audit revision (2026-10-16): connected the recorded detection similarity/reciprocal searches to final classification. The wrapper verifies target sequences, hit-table receipts, catalog labels and reference sequence hashes; no manual identity fields, experimental annotation or new inference is required. Version 1.0 policy explicitly distinguishes resolving full-sequence reciprocal support, weaker domain-limited family support, unique exact reference sequence identity and distinctive architecture with corroboration. Grade A is available at family level; unsupported paralog specificity remains unresolved or is retained at a supported common parent. Duplicate library labels do not multiply evidence. Same-accession non_adhesome membership annotations do not negate protein identity; genuinely competing accessions remain constraints. Components lacking a family rule receive a separate component identity rather than an invented family. Added explicit protein_call, assignment_level, decision_confidence, decision_basis, evidence_channels and review fields. Corrected PDLIM from calponin homology PF00307 to PDZ PF00595 before LIM PF00412, using already bundled models; reference labels/sequences remain unchanged and metadata records the rule revision. Benchmark weighting remains scoped to orthogroup context; the standalone audit also discovers and verifies a sibling report. All 25 audit, 12 workflow and 80 detection tests passed (117 total); resource validation passed. A separate replay of 314 recorded validation proteins produced A/B/C/D counts 230/10/35/39, with 174 component, 39 family and 17 subfamily high-confidence calls. Of the A calls, 222 are exact reference sequence identities; the replay is not an independent accuracy estimate. Original annotation, search, orthology, benchmark and audit outputs were preserved. No new BLAST, OrthoFinder or phylogeny was run. See COMPUTATIONAL_AUDIT.md for the scheme and audit-only rerun command.


Two-stage template correction (2026-10-16): supersedes the preceding computational
audit revision and its illustrative grade counts. Stage 1 now emits general
family hypotheses and the supplied candidates-sheet layout; motif checks are
restricted to the specified ECM/receptor scope. Stage 2 evaluates only preserved
initial hypotheses plus complete PINCH/parvin recoveries within the same protein
pool. Restored mandatory diagnostic architecture/CDD gates and the 4/3/2/1 score;
A/B/C/D follow the supplied family-level protocol. Removed catalogue-wide label
replacement, exact-reference exemptions and separate component-call routes.
Stage_2_output.xlsx follows all nine template tabs and adds current motif,
benchmark, panel-species and provenance fields. Missing orthology is distinguished
from observed panel restriction. Context-supported motif evidence is retained.
Optional phylogeny consumes the final family table and remains independent of
audit completion. Recovery checks shared schema and contract hashes. Original
templates and existing run outputs are preserved. Version remains 1.0.

Validation for this correction: 159 regression checks passed across detection,
motifs, orthology/phylogeny, audit and workflow integration; two additional tests
were skipped because Biopython is unavailable. Resource validation passed.
A separate read-only replay retained all 314 recorded species/protein pairs,
evaluated 2,198 original hypotheses and added two complete recovery hypotheses
within that pool. Grades A/B/C/D were 96/0/112/106. These are classifications of
the existing validation handoff, not an independent accuracy estimate. Current
stage 1 classification was also refreshed successfully from cached detection
evidence into a separate review directory; no new BLAST, OrthoFinder or tree
inference was launched.

Generic screening buckets cannot pass stage 2 as assigned families, even with HOG
support. A dedicated regression check enforces this boundary. Named general
family rules remain eligible; legacy screening hypotheses remain in the audit trail.


Putative identity and report cleanup (2026-10-16, version 1.0): removed MULocDeep,
entirely blank and unused report columns; consolidated semantic aliases and moved
repeated per-species detail out of the final view. The final table uses
putative_audited_family, protein_name and assignment_confidence, with a matching
TSV. Added 45 scoped, architecture/CDD-based identity refinements from the supplied
examples, keeping inferred regions and lineage assertions separate from observed
evidence. Observed CDD footprints may complement Pfam domains; overlapping copies
are merged and CDD architecture is displayed N to C. Large CDD lists are supplied
through the JSON sidecar rather than Excel cells. Naming/report code and library
hashes invalidate stale audit checkpoints without repeating orthology. Numerical
identity probabilities in that initial replay awaited labelled calibration data; orthogroup F1
continues to weight orthology only. Optional family trees accept the compact family
column and legacy workbooks. A cached replay retained 314 proteins, with 43 populated
final columns and A/B/C/D counts 100/0/108/106. Those counts are not an accuracy
estimate. Tests and original run data remain separate from the source package.

Positive-control calibration and negative exclusion (2026-10-16, version 1.0):
the supplied 583-label/581-sequence calibration library now provides positive-only
identity correctness calibration. Hard/background negatives are exclusion gates,
never identity training classes. Exact negative matches, or a complete related
architecture with an observed negative-specific curated CDD model, are dropped
from final assignments and preserved in a separate exclusion TSV and screening
trail. Shared domains alone cannot exclude a protein. This supersedes earlier
release entries describing full final-pool retention or negative-control protein
identity calls. Retained plus excluded proteins must equal the unchanged input pool.
The current cached replay retains 236 proteins and excludes 78 from 314 inputs.
The fitted isotonic model uses 130 named positive controls, removes a query's exact
sequence group before estimation, and exports probability, approximate Wilson
95% bounds, control count, method and scope. No probability is emitted for
unassigned calls, unsupported positive-family scope, absent sequence hashes,
out-of-range scores, stale models or changed reference panels. Reference-overlapping
controls do not establish independent novel-lineage accuracy. Recalibration is
automatic from scored evidence with calibrate_identity.py; no manual JSON is needed.

Computational-confidence policy (2026-10-16, version 1.0): an unambiguous grade A
family/protein assignment is strong after its required architecture gate and
computational orthology support pass. Removed the additional optional CDD confidence
cap and the experimental-metadata requirement in reference import/building.
Experimental flags, assertions and annotation counts do not affect grades, scores,
names or confidence; reference experimental metadata is provenance only. Required
family-specific CDD gates, computational conflicts and negative exclusions remain
evidence constraints. Positive-control probability calibration is refitted after
the confidence-stratum change; its empirical correctness is not an experimental
validation weight.

In-place detection recovery (2026-10-17, version 1.0): diagnosed the Sguin SIGABRT
as a nucleotide mRNA file supplied as protein input, including records exceeding
HMMER's 100,000-residue limit. Protein-input preflight now catches this before
analysis. `start --resume --exclude-target Sguin` continues the same failed handoff,
verifies/reuses completed targets, processes pending targets and records exclusion
without archives or backups. Failed targets removed from targets.json are also
supported; completed targets cannot be silently removed. Excluded target files
remain untouched and are omitted from all exported analysis inputs. Preparation
and completed domain scans can be reused using verified receipts/artifact hashes.
All 90 detection and 14 workflow fixture tests passed. No sequence translation or
DNA-as-protein scanning is used. See DETECTION_RECOVERY.md.

Post-audit control evaluation (2026-10-17, version 1.0): added evaluate_validation.py
and automatic evaluation after audit calls for the validation target. Reports
compare recruitment, complete audited family gates before/after negative exclusion,
strong grade A calls and exploratory provisional coverage. They include F1, PPV,
recall, specificity, accuracy, binary and family confusion matrices, full-cohort
control accounting and source hashes. Positive controls lost before audit remain
false negatives; unlabelled controls and known sequence mismatches are explicitly
unscored. Labelled accessions lacking duplicate trusted sequences remain in the
primary cohort and are flagged, with a sequence-verified comparison also exported.
Wrong-family calls are reported separately from membership decisions. The report
is a reference-overlapping control self-test; its scores do not change calls,
confidence, calibration or orthology weights. Evaluation errors preserve completed
audit results. See VALIDATION_PERFORMANCE.md.

Annotation recovery (2026-10-17, version 1.0): accept the retained
unassigned_domain_candidate screening alternative in motif annotation, preserve
its sequence/domain/context evidence, and exclude it from pooled family motif
discovery. Added resume --retry-failed for verified failures before motif output;
the same directory and failed logs are retained without backups. Changed inputs,
completed annotations, partial motif output and running attempts are refused.
Older failures can be adopted from an identical failed command and unchanged
staged context. Detection is reused and no biological analysis is launched during
installation.
