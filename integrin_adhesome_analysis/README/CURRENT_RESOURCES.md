# Current consolidated resources — workflow version 1.0

Detection and audit now share `detection/config/family_rules.json`: 96 current
hierarchical screening/family/subfamily entries. Every family from the older JSONs
was already present in the newest file. The newest complete entry wins on name
collisions. Superseded rule JSONs, the legacy builder and the duplicate audit
implementation were removed after consolidation; source hashes are recorded in
`config/resource_resolution.json`.

The Pfam library now contains 176 models, including the two gates introduced by
the current rules: PF00046 and PF03535. The old 96-model library and duplicated
targeted model sources were removed; all their accessions are covered by the current
library and its versioned source/provenance files.

The current audit is `annotation/integrated_adhesome_pipeline.py`.
It follows the supplied nine-sheet Stage_2_output.xlsx template and refines the
same input species/protein pool under the architecture-gated 4/3/2/1 scheme, with
hard/background-negative matches excluded from final assignments. Original
hypotheses, competing families, collapse counts and complete PINCH/parvin
recoveries are retained. No catalogue-wide replacement or exact-reference
component-call route is used. See [the two-stage protocol](COMPUTATIONAL_AUDIT.md).

Stage 03 uses the supplied evidence-based motif snapshot. Its default is the
53-rule `Motif_library/adhesion_core.json`. `--motif-scope master` selects all 347
deduplicated current rules. Current JSON and XLSX forms are supported. Legacy
excluded motifs, reader annotations and old forward/reverse parvin rules are not
reintroduced. `context_rules.json` points to the canonical domain rules rather than
duplicating an older domain catalogue. Sequence hits retain motif definition
evidence, context requirements, source aliases, domain/topology annotations and
site-specific evidence; their default tier remains `sequence_only_candidate`.
Definition evidence is not hit-level functional confidence.

## Local test installation

Copy the updated workflow source folder and associated resources into your test
installation. Copying only scripts would omit the consolidated JSONs and updated
HMM source. Keep your existing inputs and runs separately. Superseded files are
absent from this source distribution; current entry points use canonical paths.

From the test workflow root in your activated environment:

```bash
python tools/validate_resources.py
python adhesome_workflow.py --version
```

The validator should report 96 families, 176 Pfam models, 230 components, 53 adhesion
motifs and 347 master motifs. It checks model coverage, hashes, motif integrity,
hierarchy links and crosswalk consistency. Pipeline versions remain 1.0;
`schema_version` identifies resource file formats independently.

Use a new annotation directory to process an existing handoff with the current motifs:

```bash
python adhesome_workflow.py resume --handoff runs/handoff --context-dir returned/context --out runs/annotation_current --annotation-only
```

If you deliberately want the full motif snapshot, append `--motif-scope master`.
Return site evidence using current `MOTIF_...` IDs. The external-evidence template
also accepts independently computed reference/subfamily evidence with a source and
the exact candidate sequence hash. Gold benchmark labels are not prediction evidence.

An existing matching orthology run can be reused for the new annotation:

```bash
python adhesome_workflow.py audit --annotation runs/annotation_current --orthology runs/annotation/downstream/orthology --out runs/final_current
```

Both runs must refer to the same unchanged handoff and panel. Old motif IDs absent
from the current definitions are ignored in audit scoring and counted in
`MotifPrediction_ignored_superseded_count`. To apply the current classifications
to family phylogeny, prepare a new project using `--audit-workbook runs/final_current/Stage_2_output.xlsx`;
see [FAMILY_PHYLOGENY.md](FAMILY_PHYLOGENY.md). Grade C candidates retain their provisional family assignment for grouping.

For an unfinished automatic run, recovery refreshes changed audits and dependent
family-tree preparation while preserving completed orthology. Completed runs are
kept as prior analyses; use new output directories to apply new resources.

Validation uses fixtures and simulated external inference; two Biopython-dependent
tests require Biopython in the test environment. The supplied rule proposals remain
architecture constraints; fixture tests do not establish biological accuracy.
The identity probability model is fitted only on labelled positive adhesome
controls and reports its reference-overlap limitations. Negative controls provide
exclusion gates, never identity training classes. The validator checks both the
model provenance and exclusion resources.

Current final outputs include compact `putative_audited_family`, `protein_name`
and evidence confidence, plus a matching TSV. See [naming and output policy](COMPUTATIONAL_AUDIT.md).
