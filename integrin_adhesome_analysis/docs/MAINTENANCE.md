# Integrin-adhesome workflow 1.0: maintenance and troubleshooting

Use [MANUAL.md](MANUAL.md) for current policy and [WORKFLOW.md](WORKFLOW.md) for
routine commands. Run maintenance commands from the maintained bundle root in
`adhesome-workflow-v1`. Resource changes, experiments and recovery are explicit
user operations; documentation cleanup does not execute them.

## Protect reusable analyses

Keep raw/prepared FASTAs, returned predictor files, GFF/gene maps, labels,
reference/model caches, provenance, saved rule snapshots, manifests, logs,
alignments, trees and native checkpoints at their recorded paths. Do not rewrite
sealed contracts, resource hashes or model training receipts to make a changed
analysis appear unchanged. Use the original tool versions for partial inference.
Do not start concurrent writers to the same run.

Audit-only replay can reuse completed annotation and orthology but cannot add
missing Stage 1 observations. Changing the reference panel requires new inference.
Changing model/rule resources can invalidate frozen contracts and numerical
calibration. New report/audit outputs belong in new output directories; do not
call these backups or overwrite the original analysis.

Legacy `recover --retry-benchmark` and some dependent-output replacement routes
can archive previous outputs. Integrated phylogeny `--retry-failed` without
`--resume` also retains historical attempt-archive behavior. For a no-archive
tree retry use the verified in-place `--resume --retry-failed` route. The isolated
wrapper stops if a legacy recovery path would archive partial evidence. Inspect
the actual CLI and checkpoint before invoking recovery; do not promise that
every recovery mode is archive-free.

## Read-only checks and synthetic tests

```bash
python -B maintenance/check_workspace.py --syntax-only
python -B maintenance/check_workspace.py
python -B maintenance/expand_detection_pfam_library.py --check-only
python -B maintenance/validate_resources.py --allow-stale-calibration --audit-resources-only
```

The first command checks resource/config/documentation paths, syntax and defaults
without CLI subprocesses. The full workspace checker adds `--help` checks from
outside the bundle root. Neither performs biological analysis. Resource auditing
reports stale reference-label/calibration provenance honestly; complete HMM rule
coverage alone does not imply every other resource is current.

Synthetic tests are separate from real runs:

```bash
python -B -m unittest discover -s tests -p 'test_*.py'
(cd resource_library/detection && python -B -m unittest discover -s tests -p 'test_*.py')
(cd resource_library/motif_annotation && python -B -m unittest discover -s test -p 'test_*.py')
(cd resource_library/orthology_phylogeny/scripts && python -B -m unittest discover -p 'test_*.py')
```

Central workflow/audit tests use temporary fixtures and simulated inference.
They cover evidence joins, recovery, isolation, reporting, current profile defaults,
historical compatibility, aliases, family gates and negative exclusions. Fixture
success is not a biological performance estimate. Dependency-related skips must
remain visible; do not describe skipped parser/inference checks as passed.

## Troubleshooting guide

| Symptom | Check and appropriate response |
| --- | --- |
| Missing `workflow_manifest.json` | Confirm the actual annotation directory and successful sealing; do not substitute an example `runs/resumed` path |
| Missing phylogeny `manifest.json` | Verify the original prepared tree project; a new folder does not recreate its membership |
| No gene mapping for pipe-delimited IDs | Use a validated per-target regex or explicit protein/gene map; do not reuse gene-tag assumptions for controls |
| HMMER 100,000-residue failure or overwhelmingly nucleotide FASTA | Check whether an mRNA/transcript file was supplied; use its matching protein FASTA or explicitly exclude the target |
| Many `not_tested` similarity/RBH values | Inspect exact-family reference availability, search-completion receipts and related-family context; do not treat absent testing as absence of homology |
| Saved OrthoFinder checkpoint asks for 3.1.5 | Restore that version; retain the saved inference inputs and parameterization |
| IQ-TREE nonzero exit, including status 2 | Inspect exact command, stderr and tree log; the exit code alone does not establish a thread or memory problem |
| Excessive gaps/no residues | Review alignment/trim QC; changing thread counts does not repair sequence quality |
| Changed alignment/tree identifiers | Inspect transport-ID mapping, aliases and membership logs; never invent a successful join |
| Fewer than four tree taxa | Intentional skip where reported, not a software failure |
| Optional tree/benchmark failure | Keep the completed audit; inspect status/recovery files rather than requiring successful trees for classification |
| Label-file checksum mismatch in report generation | Use the original frozen labels. The calibration-only semantic-equivalence flag is narrowly scoped and does not relax validation truth |

Completed targets can be retained during `start --resume`; the historical Sguin
transcript input was explicitly excluded rather than translated or relabelled.
Preflight now rejects invalid protein inputs before expensive scanning. An old
standalone evidence CSV may document observed hits, but it cannot establish
search coverage for families with no rows.

## Pfam and component-reference resources

The current 214-model library includes the 38 accessions previously absent from
the expanded family catalogue. The original 176 records/versions/GA thresholds
were preserved and the fallback split is empty. Inspect
[pfam_expansion_20261026.json](../maintenance/pfam_expansion_20261026.json),
`resource_library/detection/databases/pfam/adhesome_v1.0.hmm.metadata.json` and the
`source_v1.0/download_provenance.json` record. Model source/version labels identify
provenance and do not change pipeline release version 1.0.

For a future reviewed rule-driven model addition:

```bash
python -B maintenance/expand_detection_pfam_library.py --check-only
python -B maintenance/expand_detection_pfam_library.py
```

The second command updates resources; it is not a read-only check. It retrieves
missing models from the documented InterPro Pfam API, preserves existing cached
model versions and analytical configuration, validates accession/state/GA data,
and rebuilds manifest-linked source/GA/fallback hashes without launching scans.
Normal adjacent temporary publication preserves inherited Windows permissions;
there is no project backup. Detection later creates its HMMER indexes when a scan
is deliberately launched. To rebuild from an already supplied complete source:

```bash
python -B resource_library/detection/build_adhesome_pfam_library.py \
  --pfam-hmm resource_library/detection/databases/pfam/source_v1.0/combined_source.hmm
```

Component reference coverage is distinct from family coverage. The active paired
panel covers all 230 crosswalk components, with 354 forward and 425 reciprocal
records; multiple labels can share one sequence and are tied by accession.
Historical coverage of 48 explicit family labels was measured against an older
96-family catalogue and must not be presented as coverage of the current 222
families. Current provenance needs a verified label refresh, not an edited hash.

For a reviewed future crosswalk update, publish a separate panel:

```bash
python -B resource_library/detection/expand_component_references.py \
  --base resource_library/detection/references/adhesome_v1.0 \
  --crosswalk config/detection/adhesome_component_crosswalk.csv \
  --out reference_panel_updated --repair-identity --workers 4
```

Reviewed human gene/taxon identity must be uniquely resolved. Historical pin
repairs include ITGA7 Q13683 and RAC1 P63000; rejected alternatives remain in
`component_coverage.csv`/`component_expansion.json`. Shared Pfams do not choose
family labels. `--offline --cache PATH` uses an existing populated source cache.
The base-seed builder cannot overwrite an already expanded component panel.

### Refresh cached detection evidence explicitly

Reclassification from unchanged caches, with new outputs:

```bash
python -B resource_library/detection/scripts/refresh_detection_evidence.py \
  --detection /absolute/path/to/completed/detection \
  --out runs/detection_reclassified
```

This copies reusable prepared/domain/homology caches and reclassifies them; it
does not run new HMMER/BLAST and cannot invent hits to added references. It needs
space for the new result, and is not a no-copy operation. To deliberately refresh
forward/reciprocal BLAST while reusing domain scans, add:

```text
--reference-set resource_library/detection/references/adhesome_v1.0 --cpu 4
```

Export that new detection through `main/adhesome_workflow.py start --detection ...`
to a new handoff before external return. The original batch/handoff remains
fixed. Newly added Pfam models require actual new scans for new domain evidence;
BLAST-only or audit-only refresh does not provide it.

## Rules, motifs and probability model maintenance

Install family and identity rules as the canonical pair. Validate rule IDs,
parent routes, executable fields, coverage and specificity without replacing
the main `calibrated_v1` algorithm. Keep `precision_v1` explicit and development
only. Do not restore calibration-derived automatic-promotion restrictions or
experimental-evidence confidence caps to production.

Keep motif definition evidence, source aliases, exclusion records and contextual
requirements together. The 2026-10-14 snapshot rebuild is not a live refresh of all
ELM classes. Its source-sheet totals (691 master rows, 353 ELM classes, 39 curated
entries) supersede a stale 694-row summary. The erroneous SLIM019 citation was
removed; it must not reappear as RGD evidence. Preserve excluded parvin orientation
rules and the prohibition on inferring hydroxyproline from FASTA O. After a new
motif resource revision, reimport/reannotate into a new annotation directory;
superseded motif IDs do not score merely because they exist in old results.

Numerical family/protein-name probabilities require positive labelled calibration
proteins, matching current evidence/rules and an explicit validated fitting
procedure. Negatives are exclusion gates, not identity-class examples. The stale
bundled model must remain disabled where its checks fail; evidence confidence can
still be strong after family gates pass. Refit is a separate analysis, never a
side effect of documentation or resource cleanup. Reference overlap and finite
label coverage must remain visible in probability reports.

The optional calibration FASTA retrieval utility in
`validation/build_calibration_set/` retrieves public reviewed human records using
requests and explicit component inputs. Its historical command accepts
`--components`, `--out`, `--background-count`, `--seed`; outputs include candidate,
hard/background-negative and combined FASTAs, mappings, issues and source records.
Its background design is enzyme-enriched and excludes candidate/hard-negative
IDs, aliases, sequence duplicates and candidate Pfams. It is not a proteome-wide
unbiased negative sample, an independent schistosome validation or an automatic
replacement for the frozen calibration dataset. Keep labelled inputs frozen
unless intentionally defining a new calibration experiment.

## Benchmark and classifier experiments

Do not adjust scoring equations or remove difficult ground-truth members to make
a failed benchmark pass. Keep root/cohort-level distinctions and the accepted
precision/F1 policy fixed before comparison. The historical embedded 22-proteome
run (12 benchmark proteomes included) had pairwise precision about 0.4004, recall
0.8346, F1 0.5411 and 20/70 exact groups; it was quality rejected. It is not an
independent rerun on only the original 12 proteomes.

New OrthoFinder settings require new inference and a separate output. For example,
an explicitly exploratory arm may use `--orthofinder-search diamond_ultra_sens`,
`--orthofinder-method msa`, `--orthofinder-inflation 1.5`, unchanged panel/ground
truth and `--compare PREVIOUS_ORTHOLOGY`. Inflation can alter merges and splits;
no improvement is established until measured. Imported or completed inference
cannot be relabelled with new settings. Tuning on the 70 RefOGs is not independent
testing of the chosen parameters.

Optional maintenance, all separate from routine inference:

| Entry point | Purpose |
| --- | --- |
| `maintenance/evaluate_precision_identity.py` | Calibration-first library comparison; same `calibrated_v1` by default, explicit `--comparison-set development` for stricter experiments |
| `resource_library/scripts/evaluate_family_expansion.py` | Paired cached-evidence family-library comparison with regression accounting |
| `resource_library/scripts/evaluate_validation.py` | Labelled membership/family evaluation contracts |
| `maintenance/optimize_stage2.py` | Calibration selection, freeze, then fixed validation ablations |
| `maintenance/stage2_optimization_report.py` | Reporting from completed optimization outputs |
| `maintenance/evaluate_orthogroup_ablation.py` | Explicit fixed-coefficient replay comparisons |
| `maintenance/recalibrate_orthogroup.py` | Calibration-only orthogroup-context exploration |
| `maintenance/audit_existing_workflow.py` | Historical provenance and reuse diagnostics |
| `maintenance/prepare_phylogenetic_identity.py` | Existing-tree inventory/concordance using independently curated labels |

Use `--help` to inspect exact required inputs without launching an experiment.
An actual experiment/replay does perform classification or scoring and needs its
own authorized output. The default cached identity comparison uses calibration
first (`--phase calibration`) and fixed subsequent validation (`--phase validation`);
changed code/libraries invalidate the frozen comparison. Do not choose rules,
weights or reporting scopes by their validation score.

Historical saved-control comparisons reported membership recovery gains under
expanded rules but also accepted-family errors, overlap limitations and unmapped
truth. They are exploratory results, not proof of perfect independent precision.
Preserve separate supported/provisional outputs, negative-exclusion ablations,
family-label coverage and overlap strata. A wrong family is not reclassified as
a broad abstention merely to improve the report. PTPRF/LAR, NF2/Merlin versus
canonical ERM, and TRIO versus spectrin remain useful family-specific regressions;
unrecruited Defensin-6 needs a separate recruitment-stage change.

## Documentation and cleanup policy

The only operational guides are `docs/MANUAL.md`, `docs/WORKFLOW.md` and this file.
The root README is a navigation index. Stage READMEs and duplicate `README/` guides
were merged and removed. The two JSON target examples were moved unchanged into
`config/`, and their test/reference paths updated. `maintenance/check_workspace.py`
checks the new guide locations and recursively validates current Markdown links
while excluding historical inputs/results. The workspace-level README points to
the consolidated guides. Historical path-repair/organization receipts retain
their original recorded filenames as provenance, not as current commands.

`resource_library/annotation/family_coverage_expansion_report.md` remains a historical curation
report accompanying its manifest/resource bundle, not an alternate operating
manual. Its old 96/57/190/136 counts describe that proposal, not the active 222/197
pair. Generated run READMEs/reports and raw resource documentation captured as
historical evidence are not rewritten by cleanup.

Only superseded operational Markdown files and disposable Python bytecode caches
outside input/output directories are removed. No scripts are deleted merely
because routine runs do not import an optional command. Keep configuration,
source caches, dependencies, tests, MEGA projects, failed-run diagnostics, sealed
outputs and checkpoints. No recursive deletion of runs, library resources,
reference inputs or separate installations is part of this cleanup.

The consolidation receipt is
[documentation_consolidation.json](../maintenance/documentation_consolidation.json).
It records superseded source-document hashes, relocated configurations and
verification results without preserving duplicate document bodies or creating
backups. After future edits rerun the workspace checker, keep relative guide
links current and update the receipt only as a new documented maintenance event.

## Source-layout maintenance

Edit settings and rules in `config/`, entry points in `main/`, explicit maintenance
utilities in `maintenance/`, and component/shared code or library assets in
`resource_library/`. Do not put editable policy beside a component implementation.
Do not move `inputs/`, `validation/`, returned evidence or `runs/` when tidying
source code. A learned model or saved run manifest is data/provenance, not a
replacement for a central workflow configuration.

`config/layout_paths.json` is the source-path compatibility map.
`maintenance/layout_migration.json` records the old and new source hashes and the
unchanged data inventory for this layout migration. Do not edit a migration
receipt to accept subsequent algorithm or scientific-rule changes. Start a new
frozen run for those changes. Existing numerical calibration is not silently
re-certified by a path migration.

Check physical guide/config links, imports, syntax and CLI help without analysis:

```bash
python -B maintenance/check_workspace.py
```

The checker launches only `--help` commands. It does not run detection, annotation,
orthology, phylogeny, audits, model fitting or performance-report generation.

Historical control-panel construction scripts are retained in
`resource_library/build_calibration_set/build_calibration_set.py` and
`resource_library/build_validation_set/retrieve_orthologues.py`. Their original
input tables/caches remain under `validation/build_calibration_set/` and
`validation/build_validation_set/`. Supply explicit input/output paths when using
these separate resource-construction utilities; they are not launched by the
main workflow or the read-only workspace checker.
