# Stage 2 Family Coverage Expansion — calibration-reference proposal

## Scope and baseline

Source: `calibration_set.xlsx`, original `family_rules.json`, previous recalibrated `family_identity_rules_recalibrated.json`, and current mapping / gate diagnostics. The existing **96 family configurations** and **57 identity rules** are preserved. All proposed changes are **reference-panel grounded**, not independently validated classifiers or measured new Stage 2 performance.

## Additions

- **94 family configurations** proposed, bringing `family_rules.json` from **96 to 190** entries.
- **79 evidence-gated identity rules** proposed, bringing `family_identity_rules` from **57 to 136**; they cover **115** existing positive controls. Each meets the specified observed Pfam architecture and matches at least one named `cdNNNNN` CDD segment per target positive.
- **15 architecture-screen configurations**, covering **16** positive controls, are screening only. They lack an observed *named curated CDD model* corroborating their proposed diagnostic pattern in one or more positive controls; they are NOT promoted into Stage 2 supported identity rules.
- Across both tiers, **131** of the previous **152** out-of-configuration positive controls now have a proposed configuration. Together with the original 78, this yields **209/235 (88.9%) potential configuration coverage**, *not observed family-recall improvement*.
- **19 family groups** covering **21** positives remain outside the new configurations; a separate five controls need trusted-family mapping review.
- **130/131** newly covered controls were already recruited by Stage 1. `Q01524` (Defensin-6, screen-only tier) was not; adding only Stage 2 rules will not recruit it.

## Priority newly proposed identity families

| Family | Positive reference controls | Status |
|---|---:|---|
| `RhoGEF_like` | 8 | Identity proposal, rerun required |
| `RhoGAP_like` | 7 | Identity proposal, rerun required |
| `ASAP_AGAP_like` | 4 | Identity proposal, rerun required |
| `PKC_like` | 4 | Identity proposal, rerun required |
| `tubulin_like` | 4 | Identity proposal, rerun required |
| `SORBS_like` | 3 | Identity proposal, rerun required |
| `ABI_like` | 3 | Identity proposal, rerun required |
| `GIT_like` | 2 | Identity proposal, rerun required |
| `INPP5_SHIP_like` | 2 | Identity proposal, rerun required |
| `MAPK_like` | 2 | Identity proposal, rerun required |
| `RNF5_185_like` | 2 | Identity proposal, rerun required |
| `Rho_GTPase_like` | 2 | Identity proposal, rerun required |
| `calpain_like` | 2 | Identity proposal, rerun required |

## Unconfigured / deferred families

The following were not promoted solely from ambiguous or generic evidence. Counts are reference positives, not novel genomes.

| Family | Count | Disposition |
|---|---:|---|
| `CIB_like` | 2 | No defensible specificity template selected from available evidence |
| `PPM_like` | 2 | Existing speculative screen could match other positives or controls |
| `LIMK_like` | 1 | Existing speculative screen could match other positives or controls |
| `PEAK_like` | 1 | No defensible specificity template selected from available evidence |
| `PPP6_regulatory_like` | 1 | No defensible specificity template selected from available evidence |
| `RACK1_like` | 1 | Incomplete raw Pfam/CDD corroboration |
| `RAPGEF_like` | 1 | Existing speculative screen could match other positives or controls |
| `RAVER_like` | 1 | No defensible specificity template selected from available evidence |
| `SPTLC_like` | 1 | No defensible specificity template selected from available evidence |
| `SRCIN1_like` | 1 | No defensible specificity template selected from available evidence |
| `TESK_like` | 1 | No defensible specificity template selected from available evidence |
| `THY1_like` | 1 | No defensible specificity template selected from available evidence |
| `calreticulin_like` | 1 | No defensible specificity template selected from available evidence |
| `layilin_like` | 1 | No defensible specificity template selected from available evidence |
| `myomesin_like` | 1 | No defensible specificity template selected from available evidence |
| `nexilin_like` | 1 | No defensible specificity template selected from available evidence |
| `polycystin_like` | 1 | No defensible specificity template selected from available evidence |
| `supervillin_like` | 1 | No defensible specificity template selected from available evidence |
| `synemin_like` | 1 | No defensible specificity template selected from available evidence |

## Static validation and safeguards

1. Preserved every existing family definition and existing identity-rule object unchanged; proposed additions are new keys only. No gate threshold, negative-exclusion policy, probability calibration, or Stage 1 detection configuration was modified.
2. For each proposed added rule, every assigned reference positive in the target expected-family label set matches both the specified Pfam architecture and named CDD text pattern. The static screen found **zero other-labelled positive** and **zero negative-control** simultaneous matches among the 583 calibrated controls. This is a *within-panel consistency check*, not an independent specificity estimate.
3. A named curated `cdNNNNN` hit was required for proposals that may reach Stage 2 identity refinement. When absent, the proposed rule is screening only and has no corresponding new `family_identity_rules` entry.
4. The revised `family_identity_rules` preserves the prior LAR, ERM, kindlin, and PTPN1 refinements from the previous calibration. New Merlin/NF2 and RhoGEF identity rules may assist two earlier incorrect accepted calls (NF2→ERM and TRIO→spectrin), but whether they actually correct them is **unverified until a full Stage 2 rerun**.
5. The modules, generic Stage 1 parent dispatch, and `allow_diagnostic_dispatch` fields are proposed for integration. Their runtime reachability and behavior under real competing hypotheses have not yet been tested.

## Mandatory implementation regression

- Save backups of original files; install the expanded `family_rules_expanded.json` and `family_identity_rules_expanded.json` together only in a new run / branch. Test JSON schema and family-ID / parent routes against the parser.
- Rerun Stage 1 only if its configured family catalogue uses `family_rules.json` for recruitment. Rerun Stage 2 after confirming that original candidate tables and identity dispatch can reach new rules.
- Verify the original 77/78 configured-positive supported calls remain correct and check whether PTPRF gets correct LAR support. Track newly supported **correct family** assignments among the additional 115 identity proposal controls; 16 screen-only positives cannot be claimed as supported based on this library alone.
- Explicit regression tests: NF2 (Merlin must not remain strong ERM), TRIO (must not remain strong spectrin), PTPRF (LAR route), and all 348 known negatives, with and without exact exclusion.
- Stratify precision, supported family recall, provisional/ambiguous rates, and macro-F1 by family-rule coverage tier; retain exact-reference/calibration-overlap stratification. Do not describe reference-overlapping controls as independent validation.
- Manually curate family memberships for the five mapping-review controls before treating them as gold labels. Avoid using updated predictions as their gold truth.

## Interpretation

The files are **configuration proposals** and an *auditable static rule-matching expansion*. They provide **no newly observed recall figure**. Actual output depends on architecture/order handling, competition, orthology evidence, Stage 1 routing, and downstream exclusions. Since many additions are supported by only one known reference protein, expanding to divergent species remains a future validation objective.
