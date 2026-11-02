# Integrin-adhesome workflow 1.0: reference manual

This is the canonical description of the maintained pipeline. Use
[WORKFLOW.md](WORKFLOW.md) for commands and [MAINTENANCE.md](MAINTENANCE.md) for
resource updates, troubleshooting and checks. Earlier operational guides were
consolidated on 1 November 2026; their historical defaults are not current policy.

## Scope and locations

The workflow reconstructs putative adhesome-associated protein families from
protein sequences and computational evidence. It accepts a general set of target
proteomes and preserves broad, supported and provisional hypotheses for external
follow-up. It does not require experimental annotation for strong computational
family or scoped protein-name assignments. Computational support does not itself
demonstrate physical complex membership, interaction or biological function.

| Location | Purpose |
| --- | --- |
| `main/adhesome_workflow.py` | Staged start, resume, continue, reference calibration, orthology, benchmark, audit and recovery |
| `main/run_pipeline.py` | Isolated validation or production target branch alongside a calibration self-test |
| `main/report_audits.py` | Reports from completed standalone calibration/validation audits |
| `resource_library/scripts/`, `resource_library/schisto_adhesome/`, `main/pipeline_paths.py` | Shared runtime modules and path resolution |
| `resource_library/detection/` | Protein preparation, Pfam scans, paired reference search and recruitment |
| `resource_library/motif_annotation/` | Returned external context and motif annotation |
| `resource_library/orthology_phylogeny/` | Reference preparation, OrthoFinder, benchmark and optional integrated trees |
| `resource_library/annotation/` | Final audit implementation and curation evidence |
| `validation/` | Labelled calibration and validation inputs and numerical calibration resources |
| `config/`, `resource_library/templates/` | All editable settings/rules; evidence templates remain resources |
| `maintenance/`, `tests/` | Explicit experiments, maintenance utilities and synthetic tests |
| `inputs/`, `runs/` | User inputs and generated analysis products; not documentation cleanup targets |

The maintained source root is
`C:/Users/eaasa/Desktop/Evans_PHD_Thesis/Integrin_analysis/Thesis/integrin_adhesome_analysis`.
The separate installed analysis directory is
`C:/Users/eaasa/Desktop/Evans_PHD_Thesis/pipeline_validation/integrin_adhesome_analysis`.
Its saved results keep their existing paths. Commands run from the source use
the source's current rules; they can read installed results through absolute
paths. Documentation maintenance does not synchronize installations or migrate
old run contracts.

The independent 15-family research phylogeny is a separate sibling bundle,
[independent_family_phylogeny](../../../independent_family_phylogeny/README.md).
Its full-orthogroup selection is different from integrated audited-candidate
family preparation. Neither bundle launches the other automatically.

## Organized source layout (2 November 2026)

`main/` contains `adhesome_workflow.py`, `run_pipeline.py`, `report_audits.py`
and the canonical `pipeline_paths.py`. Maintenance and resource-update utilities
are in `maintenance/`. Component code, shared runtime packages, curated reference
FASTA/HMM banks and evidence templates are in `resource_library/`.

All editable configuration lives in `config/`: the target/panel examples and
isolated workflow settings at its top level; detection manifests/crosswalk/rules
in `config/detection/`; family identity, registry and audit profiles in
`config/annotation/`; motif rules in `config/motif_annotation/`; benchmark policy
and optional tree environment in `config/orthology_phylogeny/`; output schema in
`config/templates/`. The main environment and requirements are also in `config/`.

Scientific inputs, saved results, numerical calibration models and provenance
are not editable workflow configuration and retain their dataset/resource
locations. The standalone family phylogeny and external validation installation
were not reorganized. Run commands from the bundle root, without changing into
`main/` or a component directory.

The explicit `config/layout_paths.json` map resolves old source paths. Small root
loaders retain the historical `python adhesome_workflow.py` command and
`pipeline_paths` imports. No full duplicate source tree is retained. Saved run
contracts are not rewritten: exact old hashes may match only the corresponding
recorded layout-only transition in `maintenance/layout_migration.json`; an unknown
hash, altered rule or later code edit still fails the frozen check. This does not
repair previously stale contracts or refit an identity-probability model.

## Environment and protein inputs

Full analysis uses Linux/WSL/HPC and the `adhesome-workflow-v1` Conda environment.
The specification pins Python 3.12, HMMER 3.4, BLAST+ 2.16.0, OrthoFinder 3.1.5,
MAFFT 7.526, trimAl 1.5.0, IQ-TREE 2.4.0, Biopython 1.85 and openpyxl 3.1.5.
`config/requirements.txt` supplies Python packages, not the external executables.
MEME is optional; external topology/localization services are not bundled.
Resume saved inference with its original tool versions and environment.

Target configuration uses a `proteomes` list with `code`, `species`, `assembly`,
FASTA `path`, `source_url` and an explicit `gene_regex` whose first capture group
extracts the gene ID. See [targets.example.json](../config/targets.example.json).
There is no fixed species set or target count. Use complete protein proteomes,
unique safe species codes without `__`, and explicit gene mappings. Preparation
selects the longest isoform per gene and exports `Code__gene_id` identifiers.
Do not infer genes by stripping arbitrary accession suffixes or relabel mRNA as
protein. Input preflight checks headers, duplicate IDs, mappings and residue QC.

[targets.validation.corrected.json](../config/targets.validation.corrected.json)
is retained only as a historical staged mixed-target example. It contains a
validation species and must not be used as an isolated production target config.

## Current resources and stable defaults

| Resource | Current role |
| --- | --- |
| `config/detection/family_rules.json` | 222 discovery family entries |
| `config/annotation/family_rules.json` | Canonical 222-family Stage 2 configuration |
| `config/annotation/family_identity_rules.json` | Canonical 197-rule identity library |
| `config/annotation/family_set_registry.json` | 159-set coverage/QC registry; not a production confidence cap |
| `config/annotation/stage2_evidence_rules.json`, `config/annotation/audit_policy.json` | Executable evidence profiles and audit policy |
| `config/detection/adhesome_pfam_manifest.csv` | 214 unique Pfam models |
| `resource_library/detection/databases/pfam/adhesome_v1.0.hmm` | All 214 GA-thresholded models; fallback split empty |
| `config/detection/adhesome_component_crosswalk.csv` | 230 component records |
| `resource_library/detection/references/adhesome_v1.0/` | Paired curated search/competition panel, 354 forward and 425 reciprocal records |
| `config/motif_annotation/adhesion_core.json` | Default 53 adhesion-associated scan rules |
| Motif master JSON in `config/motif_annotation/`; source workbooks in `resource_library/motif_annotation/Motif_library/` | Optional 347-rule master scope |
| `config/isolated_workflow.json` | `expanded`, `calibrated_v1`, fixed orthogroup context coefficient 1.0, full benchmark, optional identity resolution disabled |

The main algorithm remains **`calibrated_v1`**. Library expansion did not enable
stricter global thresholds. Calibration-derived inspection/no-auto-promotion
restrictions have been removed from active libraries and are not recreated from
historical registry QC fields. `precision_v1` is an explicitly selected development
comparison. It is not a routine default. `src`, `Src` and `SRC` normalize to `Src`;
`SCR` is not an alias.

`expanded` selects the current canonical family/identity pair together. The
compatibility selector `original` uses current files; it does not restore deleted
historical libraries. Historical comparisons must use a saved configuration
snapshot. Pipeline release version 1.0 is separate from resource schema versions,
model accession versions and historical source filenames.

The 38-model HMM expansion preserved the original 176 model records, versions and
GA thresholds. It does not add observations to old scans. The receipt is
[pfam_expansion_20261026.json](../maintenance/pfam_expansion_20261026.json).
Component-reference label provenance predates the latest family catalogue and
needs a separate verified refresh. The numerical family/name probability model
is stale pending an explicit local refit; do not falsify its provenance to make
it current. Computational family support remains available without that fit.

## Stages and evidence interpretation

The staged flow is:

```text
protein preparation -> Pfam/reference detection -> supported/provisional handoff
  -> external predictor return -> context and motif annotation
  -> automatic reference preparation/panel -> orthology -> benchmark -> final audit
  -> optional phylogeny
```

Stage 1 is broad recruitment, not a definitive identity classifier. Its handoff
retains alternatives for the same protein, including fibronectin-like hypotheses.
Hypothesis-row counts can exceed unique candidate counts. Topology, localization
and phylogeny are deferred to the resumed workflow; detection does not call
external predictor services.

The final audit refines the same protein pool, preserves original alternatives
and collapse/recovery counts, and cannot invent proteins or unobserved domains.
Architecture/orthology/motif/localization maxima are 4/3/2/1. Full architecture
with sufficient orthology can yield grade A, full architecture without enough
orthology grade B, partial architecture with orthology provisional C, and
insufficient architecture unassigned D. Family-specific diagnostic and identity
gates still apply. Shared kinase, PH, LIM, FERM or FN3 domains alone do not resolve
every specific family or paralog. Missing experimental evidence is not a downgrade.
Optional CDD is not an extra cap unless the family rule requires it.

Positive adhesome controls contribute classification/calibration truth. Hard and
background negatives are exclusion gates, never positive family training classes.
Qualifying negative matches are dropped from final assignments and retained in
`negative_control_exclusions.tsv`; a generic shared domain is insufficient for
exclusion. Retained plus excluded records must account for the Stage 1 protein
pool. Stage 1 misses remain in evaluation denominators.

### Detection availability fields

| Field/value | Interpretation |
| --- | --- |
| `order_check=not_applicable` | Single-type/repeat rule; there is no distinct domain-type order to test |
| `order_check=not_configured` | A multidomain order constraint was not configured; not a passing check |
| `order_check=incomplete_or_reordered` | Configured order was not recovered |
| `similarity=supported` | A qualifying exact-family labelled reference hit exists |
| `similarity=no_qualifying_hit` | Verified search and exact-family references exist, but no hit qualified |
| `similarity=not_tested` | Use its reason to distinguish absent family references from undocumented search completion |
| `reciprocal=not_tested` | Often no qualifying forward hit exists for a reciprocal assessment |
| `reciprocal=not_supported` | Forward support exists but reciprocal family agreement was not recovered |

Reason/coverage fields, run-local paired-search manifests and reference hashes
separate missing evidence from negative evidence. Parent/related-family context
does not become independent subfamily support. `reference_coverage.csv` and the
handoff `reference_coverage.tsv` expose gaps. Adding downstream orthology proteomes
does not retroactively add detection BLAST references or hits.

Crosswalk domain overlap is a component hypothesis, not identity evidence. The
component panel has reviewed gene/taxon-verified records, tied-accession handling,
catalogues, coverage reports, `reference_lock.json` and expansion provenance.
Unresolved positive component hits remain provisional
`component_reference_candidate` records; they pass external return but are
excluded from heterogeneous pooled family motif discovery. The optional
242-accession `HsaAdhesomeReferencePanel.faa` plus scope/hash sidecar is a targeted
human panel, not a complete proteome or evidence of absence. Do not represent it
as an additional human species beside a complete proteome for that taxon.

### CDD and motifs

The importer retains CDD hit category, accession/model, coordinates, E-value,
bit score, incomplete-domain status and reported model cutoffs/coverage when
available. `CDD_hits.tsv` contains all categories. Main-profile
`CDD_specific_hits.tsv` keeps the historical specific-only, E-value <=1e-3 intake.
Overlapping Pfam footprints observed through Pfam and CDD are not independent
votes. The stricter development profile adds specificity/competition safeguards;
it does not silently replace main-profile behavior.

Motif evidence describes definitions and contextual hypotheses, not functional
validation of a new regex hit. The supplied 2026-10-14 motif snapshot has 347
deduplicated scan patterns (53 adhesion, 294 general regulatory), 341 retained
ELM classes and six primary-backed additions. Original source-sheet counts take
precedence over stale workbook summary counts. Canonical definitions replace
legacy generic patterns; excluded parvin forward/reverse rules, generic
polyproline/PEST/basic patterns and reader-only architecture records are not
reintroduced. Exact-regex aliases are scanned once; distinct contextual classes
remain distinguishable. PINCH LIM architecture is not a SLiM scan rule.

Python regex scans use the original uppercase full sequence and overlapping start
positions. `$` means the true protein C terminus. Hits spanning X/B/Z/J, gaps or
stops are rejected. O is not automatically converted to P or treated as annotated
hydroxyproline. Collagen proline hydroxylation/triple-helical context and host-to-
parasite interaction transfer cannot be established from sequence matching.
Site-specific, reader, topology, accessibility, PTM and conservation evidence stay
explicit. No uncalibrated hit-level confidence percentages are assigned.

## References, orthology and reference-panel expansion

Ordinary raw reference FASTAs and applicable GFF/gene maps belong in
`inputs/reference_raw/`. FASTA/gzip variants are accepted. Gene preparation uses
matching `Code.protein_gene.tsv`/`Code.tsv`, NCBI-style `Code.gff`/`Code.gff3`
protein-to-GeneID CDS annotations, explicit `gene=ID`/`gene:ID` headers or already
prepared `Code__gene_id` headers. Every protein needs an explicit mapping; the
longest isoform per gene is retained. Automatic preparation writes panel JSON;
manual JSON editing is not a required pause after external return.

Reference calibration means adding appropriate reference proteomes to a new
panel/inference and comparing that result. It does not rename parameters on an
already inferred partition or silently mutate an existing panel. Fibronectin-like
references/hypotheses follow the same export, return and expansion path.

## Benchmark scope, acceptance and weights

The default is the original Open_Orthobench scorer with `all` scope. Original
proteomes and ground truth are supplied under:

```text
inputs/reference_raw/BENCHMARKS/
  Input/                         # 12 original FASTAs
  RefOGs/                        # RefOG001.txt ... RefOG070.txt
    low_certainty_assignments/   # supplied uncertainty exclusions
```

The original protocol expects 12 proteomes, 70 RefOGs, 1,945 memberships and
251,378 input proteins. Do not replace them with modern proteome releases,
rename their IDs, invent gene mappings, or derive ground truth from predictions.
Custom curated data require explicit `custom` protocol selection. Original
benchmark protein units are not collapsed to guessed genes. Benchmark-only
normalization removes one terminal stop and encodes internal stops as X in
prepared copies, with positions/hashes in `sequence_qc.tsv`; raw inputs stay
unchanged and scoring verifies the same declared processed sequence scope.

The explicit alternative `adhesome` scope invokes
`resource_library/orthology_phylogeny/scripts/adhesome_benchmark_subset.py`. It needs a curated
selection TSV with `refog`, `family`, `curation_source`, `curation_status=curated`.
Use [adhesome_refog_selection.tsv](../resource_library/templates/adhesome_refog_selection.tsv).
Several reference groups may map to one family, but truth groups are never
merged. Fibronectin-like scopes are eligible where curated associations exist;
FN3 domains alone do not create ground truth. Missing, empty or invalid selection
does not fall back to full scoring. A full-benchmark failure never auto-selects
the subset.

Both scorers retain complete predicted membership within the benchmark-proteome
universe. Foreign members and members of unselected truth groups still penalize
incorrect merges. Proteins in unlabelled non-benchmark taxa are outside that
comparison universe, not automatically false positives. Original uncertainty
exclusions remain RefOG-specific. OrthoFinder 3.1.5 root output is selected from
`Orthogroups/Orthogroups.tsv`; older root/node table conventions remain version
aware. Optional `cohort-lca` chooses the smallest species-tree node containing
the full benchmark cohort by taxonomy, never by the most favorable score.

Acceptance requires unrounded **pairwise precision >=0.80 and F1 >=0.70** under
`config/orthology_phylogeny/benchmark_policy.json`. No extra recall or exact-group
threshold is implied. Here precision means PPV, not area under a precision-recall
curve. Input validity and quality acceptance are separate. A historical result
(precision 0.4003824, recall 0.8345527, F1 0.5411462, 20/70 exact groups) is a
quality-rejected result, not a defect repaired by changing scoring formulas.

`quality_rejected` reports retain measured scores and do not export accepted
benchmark weights. The audit can finish despite missing/rejected/failed optional
benchmarking. An accepted score may weight orthogroup context only; it cannot
override architecture, identity or negative gates. Subset family weights apply
only to that exact tested family, and zero remains zero. Untested families have
an explicit neutral/unvalidated baseline. Aggregate F1 is not a protein identity
probability. Routine isolated runs freeze the context coefficient at 1.0 and do
not silently refit it from benchmark or validation scores.

## Calibration, validation and reporting contracts

The frozen isolated calibration population has 581 records: 234 positives,
119 hard negatives and 228 background negatives. The labelled validation table
selects 1,737 IDs: 902 positives and 835 negatives from a 1,810-record FASTA.
The other 73 are explicitly reported as unlabelled exclusions. These populations
are different from older 235-positive mixed reference-panel controls; historical
metrics from different cohorts are not interchangeable.

Calibration and validation/target branches share a frozen algorithm/configuration
contract but have separate inputs, candidate IDs, returned evidence and orthology
inference. A branch cannot mutate another branch's orthogroups. Production target
config reject evaluation-control codes/IDs. Exact calibration/reference/training
overlap and prior rule-development inspection limit independence claims, even
when exact matches are absent. Preserve overlap-free strata and label coverage.

Report Stage 1 recruitment, supported assignments and provisional annotations
separately. PPV=TP/(TP+FP), recall=TP/(TP+FN), F1=2TP/(2TP+FP+FN). Undefined
denominators remain null/blank/`n.a.`. A wrong accepted family can be a correct
membership call but a family FP and FN. Unmapped expected-family labels are
excluded only from family correctness and explicitly counted. Exact protein,
paralog or isoform naming accuracy is not established by membership F1.

Reports retain exclusions/ablation, missed positives, wrong family calls, family
micro/macro metrics, species/length/negative-class strata and calibration/reference/
training overlap. Evaluation labels are joined after inference, not supplied as
prediction evidence. Known-negative exclusion performance is distinguished from
generalization to unseen homologues. The Open_Orthobench score measures groups,
not the independent correctness of a family name.

New isolated validation runs/resumes automatically generate
`validation/reports/validation_review.html` and `.md` after audit. If calibration
finishes later, its comparison is refreshed without rerunning either analysis.
Standalone audit replays require a separate reporting command; they do not rewrite
the original dashboard. Reporting failures leave analytical checkpoints reusable.

## Final exports and optional evolutionary interpretation

`Stage_1_output.xlsx` retains broad recruitment hypotheses. The final
`Stage_2_output.xlsx` follows the nine-sheet template and includes actual panel
species. `audited_family_assignments.tsv` is the compact equivalent. Final
`putative_audited_family`, `protein_name` and `assignment_confidence` avoid duplicate,
unused/empty and MULocDeep columns. Supporting outputs include
`joined_evidence.tsv`, `negative_control_exclusions.tsv`, rule/gate traces,
`annotation_decisions.tsv`, `cdd_specificity_audit.tsv`, `family_set_review.tsv` and
`family_configuration_reconciliation.json`. Detailed adjudication stays outside
Excel to avoid cell truncation. Configuration snapshots and workflow manifests
record the exact evidence/rules used.

Phylogenetic inference (`--phylogeny`) and identity resolution from existing trees
(`--phylogenetic-identity-resolution`) are separately optional and disabled by
default. Neither may prevent the family audit from completing. Integrated family
tree preparation selects audited A/B/C candidates by family, keeps candidates
without orthologues, and adds linked reference members without an arbitrary cap.
It excludes grade D and negative-excluded proteins. This differs from the sibling
standalone bundle's complete explicitly selected orthogroups.

Optional identity evidence requires exact protein/sequence/tree hashes, curated
reference labels and supported pure clades. The current adapter checks at least
two reference proteins from two taxa, SH-aLRT >=80 and UFBoot >=95; these are
screening criteria, not calibrated probabilities. Subfamily, ortholog/co-ortholog,
paralog and duplication interpretation additionally need declared lineage scope,
justified rooting and reconciliation provenance. An arbitrary display root and
family concordance cannot transfer a vertebrate paralog name. Outcomes preserve
family support separately from `identity_supported`,
`family_supported_subfamily_unresolved`, `conflicting_identity` and
`identity_unresolved` at the reported scope.

See [WORKFLOW.md](WORKFLOW.md) for the operational commands and
[MAINTENANCE.md](MAINTENANCE.md) for recovery limitations and library upkeep.
