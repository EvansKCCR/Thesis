# Two-stage family classification — workflow 1.0

Stage 1 proposes general adhesome family hypotheses. Stage 2 refines those same
species/protein pairs using architecture, orthology, current motifs and returned
topology/localization. Experimental annotation is not required. This protocol
supersedes the earlier exact-reference/component-call audit design.

## Stage 1: recruitment and handoff

Sequence similarity, reciprocal homology, profiles and family architecture
recruit complete/partial architecture hypotheses or independent family evidence.
Complete compatible architecture with similarity and reciprocal support may be
supported; weaker evidence is provisional, conflicts are flagged, and generic
evidence remains unresolved. Multiple independently supported families become
provisional for stage 2. Subfamily search labels fold into general families.

Stage 1 motif checks are limited to collagen-like, fibronectin-domain/like,
vitronectin-like and integrin receptors. Topology, localization and phylogeny are
deferred; their template columns contain a deferred marker.

`start` exports `Stage_1_output.xlsx` with the template's `candidates` sheet,
`hypotheses.tsv`, candidate FASTA and the external-evidence return template.
All original alternatives for selected proteins are preserved, including
unresolved alternatives. Fibronectin-like follow-up hypotheses remain included.
Targets are configurable. Blank template rows are formatting, not candidate quotas.

## Stage 2: constrained family refinement

| Evidence | Maximum | Interpretation |
|---|---:|---|
| Architecture | 4 | Family domain types/copies, order, coordinates, length and diagnostic CDD |
| Orthology | 3 | One/two/at least three reference species in the candidate HOG score 1/2/3 before weighting |
| Motifs | 2 | Context-supported sites score 2; sequence-only or priority sites score 1; contradicted-only sites score 0 |
| Localization/topology | 1 | Compatible returned predictions add support |

The complete architecture gate requires score 4. PINCH-like requires exactly
five tandem LIM domains, length 250–500 aa inclusive and PINCH-specific CDD.
Parvin-like requires at least two ordered CH domains and parvin-specific CDD.
Syk requires SH2–SH2–kinase plus Syk/ZAP70 CDD. Csk requires SH3–SH2–kinase plus
Csk/Chk CDD. A compatible core lacking diagnostic CDD is `Src_Csk_like`.
Nck, Grb2, Crk and Shc require their configured ordered cores and specific CDD.
Supported stage 1 baselines are rechecked against domain types and copies.

| Grade | Decision |
|---|---|
| A | Complete architecture plus positive weighted orthology; retained at family level |
| B | Complete architecture without enough orthology; retained at family level |
| C | Partial architecture plus orthology; provisional family assignment |
| D | No required architecture support, or incomplete architecture without orthology; unassigned |

Observed target-only B assignments are restricted to the sampled panel; unavailable
orthology is explicitly labelled unavailable. No species is hardcoded.
Missing topology loses its bonus but does not override grades; conflicts remain
visible. Exact matches and motifs cannot bypass the architecture gate.

Only complete PINCH-like and parvin-like recoveries may be added for existing
stage 1 proteins. Hypotheses rank by assignment status, family decision, grade,
domain support, architecture-specific recovery, orthology, motif, localization,
total and candidate status. Specific CDD helps resolve equivalent family decisions.
One representative is selected per species/protein, with competing initial
hypotheses, ties, original/recovery counts and collapsed-record counts preserved.
Hard/background-negative matches are then excluded from the final assignments.
The pipeline asserts that retained plus excluded species/protein pairs equal the
input set, without overlaps or added proteins.
It does not search the entire catalogue for replacement labels or create a
separate component-call dataset.

## Template outputs and current evidence

`Stage_2_output.xlsx` follows the supplied nine sheets, in order:

1. `Pipeline_summary`
2. `Screening_hypothesis`
3. `Audited_family_assignment`
4. `Motif_annotations`
5. `Orthology_candidate_report`
6. `Orthology_HOGs`
7. `Direct_orthologues`
8. `NCBI_CDD_selected_hits`
9. `Topology_predictions`

The supplied nine-sheet structure is retained. The final assignment table uses
the compact named view described below; supporting tables retain populated evidence
columns. Reference-species details reflect the actual panel.
Current motif IDs are scored; superseded IDs are excluded and counted. Supplied
site-specific support/conflicts retain their tiers. Definitions alone cannot
promote a match to context-supported. Unavailable source fields remain blank.

Provenance-verified Open_Orthobench reference-cohort F1 scales orthology only.
Both raw and weighted orthology scores, status and source are exported. This F1
is not a candidate probability or family-specific accuracy estimate. Missing or
invalid benchmarking is explicitly unvalidated and uses neutral weight 1.

## Reaudit without repeating inference

From the directory containing `adhesome_workflow.py`, in the active environment:

```bash
python tools/validate_resources.py
python adhesome_workflow.py audit --annotation runs/annotation --orthology runs/annotation/downstream/orthology --benchmark runs/annotation/downstream/benchmark --out runs/final_two_stage
```

Read `runs/final_two_stage/Stage_2_output.xlsx`. Choose a new output directory.
This reuses the unchanged handoff, annotation, orthology and benchmark. It restores
omitted original hypotheses from checksum-verified detection tables only for the
same selected proteins. Existing results remain unchanged.

To refresh general stage 1 hypotheses from a completed cached detection batch:

```bash
python 01_detection/scripts/refresh_detection_evidence.py --detection runs/handoff/detection --out runs/detection_two_stage
python adhesome_workflow.py start --detection runs/detection_two_stage --out runs/handoff_two_stage
```

Then follow the external-evidence return/resume steps in [WORKFLOW.md](WORKFLOW.md).
The audit-only route refines an existing handoff; it does not silently replace it
with new detection. Phylogeny remains optional and runs after audit. A tree failure
cannot prevent the audit. Recovery detects changed audit code, policy, schema and
rules while preserving completed orthology.

Version remains 1.0. Template hashes are in `templates/output_schema.json`; audit
manifests hash the inputs, rules, policy, join adapter and schema. Regression
checks verify implementation and complete protein accounting, not biological accuracy.


## Putative names and compact final report

The final `Audited_family_assignment` view now exposes `putative_audited_family`,
`protein_name`, `assignment_confidence` and `protein_name_scope`. This replaces
repeated family-name aliases. Its matching `audited_family_assignments.tsv` has
the same populated columns. Blank, unused and MULocDeep columns are removed;
known duplicate aliases are consolidated. Numeric zero and False evidence remain
visible. Detailed per-species HOG members and source search fields remain in the
supporting sheets and joined evidence, rather than being repeated in the final
view. A manifest records which columns were omitted and why. The nine sheets
remain; an unavailable evidence sheet states that no records were supplied.

Observed coordinate-bearing CDD footprints can complement missing Pfam domain
models. Overlapping annotations of the same domain count once, and inferred
regions do not supply domain copies. The displayed CDD architecture is ordered
from N to C. Large structured CDD lists are consumed from the JSON sidecar rather
than Excel cells, avoiding Excel's 32,767-character cell limit.

The version 1.0 `family_identity_rules.json` adds 45 scoped refinements based on
the supplied architecture examples. Rules test only related stage 1 hypothesis
scopes and require a compatible observed core, an observed Pfam anchor and
diagnostic CDD support. Examples include HIP1/HIP1R, actin-related lineages,
spectrin, PDZ-LIM, Src/Csk/Abl/Tec, collagen/laminin subtypes and Ig/FN3 receptor
families. A single generic FN3, LIM or kinase domain does not name these proteins.
The example table is preserved as resource provenance, not per-protein gold labels.
An inferred transmembrane segment, loop, subclade or restricted lineage in that
table is not treated as an observation. Copy counts alone do not establish a
lineage-restricted or species-restricted paralogue.

An unambiguous grade A call has strong family and protein-identity confidence at
the reported scope: the required family architecture gate and computational
orthology support are sufficient. Optional CDD corroboration is not an extra
confidence gate; family-specific CDD remains required where the architecture rule
needs it to distinguish related families. Experimental evidence is neither required
nor scored, and its absence causes no penalty to assignment, annotation confidence
or evidence weights. Reference import accepts computational sources without an
experimental-evidence field. Experimental metadata, when present, is provenance only.
Distinct isoforms, gene symbols or phylogenetic subclades are not inferred from
broad HOG membership. Equally ranked or competing identities remain explicitly
ambiguous. No experimental annotation is required. This interpretation follows
[NCBI's distinction between specific CDD hits and broader domain evidence](https://www.ncbi.nlm.nih.gov/Structure/cdd/docs/cdd_how_to_protein_function.html).

`assignment_confidence` describes evidence convergence. The supplied
`02_calibration_set/calibration_set.xlsx` and `calibration.faa` now support an
empirical isotonic correctness estimate at the reported family/name scope.
`assignment_probability`, approximate Wilson 95% bounds, calibration count, method
and evidence scope accompany eligible assignments. These are not probabilities of
an exact gene, isoform or paralogue identity. Scores and orthogroup F1 are not
converted directly into probabilities. Small strata fall back to pooled controls;
scores outside the observed range, unassigned calls, missing sequence hashes or
changed reference panels receive no numerical estimate.

Only positive adhesome controls train this model. The current cached validation
replay supplies 130 named positive controls; 103 positive controls with abstained
assignments do not train correctness of a reported name. Both correct and incorrect
named predictions enter calibration; gold labels never set predicted names.
The labels contain 583 proteins and the FASTA contains 581; missing sequences are
reported. Exact query sequences are removed in full before fitting their estimate.
Many controls overlap the detection reference panel, so this measures conditional
control-cohort performance, not independent accuracy for novel proteins or lineages.
The review TSV retains expected labels, predictions, exclusions and reference
overlap. Methods follow [isotonic calibration](https://scikit-learn.org/stable/modules/calibration.html)
and [Wilson intervals](https://www.itl.nist.gov/div898/handbook/prc/section2/prc241.htm).

Hard-negative and background-negative controls are exclusion gates only. An exact
negative-control sequence is removed even if the classifier abstained. A novel
sequence can be excluded when a complete related family architecture and observed
curated CDD model identify a negative subfamily; that model must be absent from all
positive controls. Shared LIM, CH, spectrin or other domains alone cannot exclude
a protein. Exact positive controls are protected from inferred subfamily gates;
contradictory positive/negative labels on an identical sequence stop the run.
Excluded proteins remain in `Screening_hypothesis`, with control class, accession,
gate and reason, and are listed once in `negative_control_exclusions.tsv`.
They do not enter the final assignment table, probability calibration or optional
audited-family phylogeny. The latest replay accounts for 314 inputs as 236 retained
and 78 excluded proteins; these are fixture-cohort counts, not a target quota.

To refit after new controls or reference-panel expansion, first complete an audit
of the controls against that panel, then run from the workflow root:

```bash
python calibrate_identity.py --features runs/calibration_audit/identity_calibration_features.tsv --report-dir runs/identity_calibration_review
```

The model is written automatically to `02_calibration_set/identity_probability_model.json`.
Reaudit to apply it. No hand-written JSON is required. Prediction code, control
files, model and reference-panel hashes prevent silent reuse after changes.

The optional family phylogeny reader accepts `putative_audited_family` and remains
compatible with older `assigned_family` workbooks. Recovery detects changes to
the naming library, identity code and report presentation, preserving orthology.
