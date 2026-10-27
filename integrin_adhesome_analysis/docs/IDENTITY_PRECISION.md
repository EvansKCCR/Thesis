# Family identity and optional phylogenetic resolution â€” version 1.0

The current canonical family/identity pair supersedes the removed `_expanded`
filenames. `expanded` remains a CLI compatibility selector for that pair.
`original` selects the detection family catalogue with the current identity
library; it does not restore deleted historical libraries. Controlled comparisons
use the immutable configuration snapshots inside a saved run.

## Family resources and decisions

The supplied resources contain 222 Stage 1 families, 222 Stage 2 configurations,
197 identity rules and 159 registry sets covering 234 calibration positives.
No separate `family_set_registry_reconciled.json` was found. The supplied registry
is authoritative until one is provided. Each discovery family has an explicit
handling policy: diagnostic refinement, architecture-only family support,
review-required proposal or screening-only fallback. Parent lists control rule
dispatch and do not define evolutionary ancestry. Duplicate rule IDs and unknown
executable gate fields are rejected. Calibration accessions are provenance and
coverage metadata; they never route a query to an identity.

**Routine runs keep `calibrated_v1`.** The active family/identity libraries no longer contain calibration-derived inspection or no-auto-promotion flags. Registry QC metadata does not restrict production calls.

The explicit, non-default development profile `precision_v1` respects `no_auto_promotion`, `user_inspection_required`,
`candidate_rescue_only`, curated-specific-hit requirements and named repeat-model
requirements. An unvalidated proposal cannot become supported merely because it
matches its own calibration reference. Existing alternative architectures remain
in the family rules, including FAK FERM alternatives. Explicitly incomplete
sequences retain missing-region limitations; no missing domains or fragments are
invented and a partial core cannot justify a specific paralog.

Additional provenance reports retain the previous family confidence and add four scoped outcomes: `identity_supported`,
`family_supported_subfamily_unresolved`, `conflicting_identity` and
`identity_unresolved`. `supported_identity_level` states the scope. A supported
family call is not a demonstrated one-to-one ortholog or exact gene identity.
`retained_stage1_family` preserves the discovery hypothesis even when refinement
abstains. Negative controls remain exclusion gates, never positive identity
training classes. Experimental metadata does not reduce confidence.

## CDD evidence

The context importer retains hit type, model/accession, coordinates, E-value,
bit score, incomplete-domain status and superfamily. Model coverage/cutoffs are
preserved when supplied and remain unknown when absent. `CDD_hits.tsv` contains
all imported categories; `CDD_specific_hits.tsv` retains the previous specific-only export and E-value <=1e-3 selection for the unchanged main algorithm.
In the optional development precision profile, non-specific, superfamily and name-only observations cannot independently justify
a diagnostic specific identity. Explicit failed model cutoffs override a reported
specific category. The exact legacy specific-only export marker remains traceable
for cached evidence. Pfam models returned by CDD share a footprint with overlapping
Pfam observations and do not count as independent votes.

This follows [NCBI's specific-hit interpretation](https://www.ncbi.nlm.nih.gov/Structure/cdd/cdd_help.shtml).
A model name alone is not a threshold test. A generic kinase, FERM, LIM or FN3
core cannot independently establish a precise family in the stricter profile.

All complete diagnostic alternatives are considered, including alternatives not
eligible for automatic redispatch. Incompatible alternatives trigger review.
Raw E-values or bit scores from different models are not treated as comparable
identity scores. No uncalibrated numerical decision margin is invented.
`src`, `Src` and `SRC` normalize to `Src`; the original input spelling is retained.

## Optional phylogenetic identity resolution

Default: disabled. `--phylogeny` requests optional inference after audit and is
separate from `--phylogenetic-identity-resolution`, which reviews existing evidence.
Neither stage can prevent the family audit from completing.

Prepare existing tree evidence without inference:

```bash
python -B tools/prepare_phylogenetic_identity.py \
  --project /absolute/path/to/stored/phylogeny \
  --reference-labels /absolute/path/to/curated_reference_tip_labels.tsv \
  --out runs/phylogenetic_identity_review
```

Reference labels require `sequence_id`, `family`, `species`, `source` and exact
tree-tip identifiers. Omitting `--reference-labels` generates an inventory and
an empty evidence template; it does not invent identities. The adapter checks
prepared inputs, complete tip membership and pure labelled clades with at least
two references from two taxa, SH-aLRT >=80 and UFBoot >=95. These are conservative
screening thresholds, not fitted identity probabilities. Rooting and duplication
reconciliation are not inferred from an arbitrary display root.

An evidence TSV requires `protein_id`, `sequence_sha256`, `family`, `clade_id`,
`tree_file`, `tree_sha256`, `sh_alrt`, `ufboot`, `reference_count`,
`reference_taxa_count`, `reference_label_source`, `placement_status` and
`duplication_status`. `placement_status` must be `pure_supported_clade`.
Contradictory placements trigger review and preserve the original family call.
Subfamily/ortholog/co-ortholog/paralog interpretation additionally requires explicit
target-lineage scope, outgroup-supported rooting, reconciled duplication evidence
and relationship provenance. Family concordance alone does not establish these.

For an audit, add:

```bash
--phylogenetic-identity-resolution \
--phylogenetic-identity-evidence runs/phylogenetic_identity_review/phylogenetic_identity_evidence.tsv
```

Isolated `scripts/run_pipeline.py` also accepts the enable flag and separate
`--calibration-phylogenetic-evidence`, `--validation-phylogenetic-evidence` and
`--target-phylogenetic-evidence`. Branch evidence is frozen independently.
Missing or changed optional evidence is recorded without cancelling the audit.

## Outputs and reporting

Audits keep the workbook/TSV contracts and add:

- `family_configuration_reconciliation.json`: family policies, aliases, conflicts,
  missing model resources and architectural alternatives.
- `cdd_specificity_audit.tsv`: reported categories, model statistics and eligibility.
- `family_set_review.tsv`: registry review flags separate from family confidence.
- `annotation_decisions.tsv`: full candidate adjudication, competing evidence,
  identity scope and optional evolutionary interpretation.
- Existing candidate/rule gate traces retain detailed rejection explanations.

Detailed adjudication stays outside Excel to avoid cell truncation. Compact
outcomes and the provenance filename remain in final exports for Atlas dossiers.
The aggregate evidence score is not a correctness probability or proof of physical
adhesome membership. Family identity, evolutionary relationships and functional
context remain separate fields.

Calibration-first cached comparisons, with no searches or tree inference:

```bash
python -B tools/evaluate_precision_identity.py \
  --run /absolute/path/to/runs/isolated_validation \
  --out runs/identity_comparison --phase calibration
python -B tools/evaluate_precision_identity.py \
  --run /absolute/path/to/runs/isolated_validation \
  --out runs/identity_comparison --phase validation
```

The default comparison holds the previous `calibrated_v1` algorithm fixed and changes only the library snapshot. The second command refuses changed analytical code/libraries. Explicit `--comparison-set development` runs the stricter exploratory ablations; it does not select the main algorithm. Controls are joined
only after classification. Reports include Stage 1 recruitment, supported and
provisional policies, confusion matrices, family metrics, overlap strata,
candidate changes, redispatch/competition ablations and an explicitly experimental
architecture ablation. Experimental assignments are not production outputs.
Exact-label metrics and correctness at the reported broad scope are both retained;
correct broad abstention is not described as a wrong specific identity.

No independent adjudicated phylogenetic benchmark or held-out clade labels were
provided. Phylogenetic comparison therefore measures evidence concordance only.
The categorical family gates do not provide a calibrated continuous family score
for a defensible probability-based precisionâ€“recall curve. Report coverage losses
and abstentions explicitly and do not claim an accuracy improvement without data.

The 38-model Stage 1 HMM resource gap was resolved on 26 October 2026: the
library now contains 214 models and covers every Pfam accession referenced by
current family rules. Existing 176 model records and GA cutoffs were preserved.
See the [library receipt](../maintenance/pfam_expansion_20261026.json). Stored
analysis outputs were not rerun and retain their original model coverage.
Component-reference label provenance is a separate outstanding resource refresh.
