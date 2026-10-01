# Curated adhesion motif annotation

The pipeline uses `Motif_library/exhaustive_comprehensive_adhesion_motif_library_PINCH_Parvin_updated.xlsx` as its sole motif source. The workbook is read without modification. Python 3.10+ is sufficient; XLSX reading uses the standard library.

## Run

From `03_functional_motif_annotation`:

```bash
python3 family_motifs.py --out results/family_inspection_expanded_context --regions results/verified_domain_regions_expanded.tsv
python3 -B -m unittest discover -s . -p 'test_*.py' -v
```
`--topology-dir`, `--topology-zip`, `--context-dir` and `--evidence` are optional; missing evidence never becomes positive evidence. External MEME stages default to `--external skip`; use `plan` to record commands or `run` to execute installed tools. FIMO/TOMTOM require a real protein PWM database through `--known-meme`. The workbook contains regexes, not PWMs.

## Library interpretation

The master workbook contains 691 records: 37 curated adhesion entries, 352 ELM classes and 302 legacy candidates. Six adhesion entries are annotation-only properties, including PINCH architecture. All source sheets and row numbers are preserved in `curated_library.snapshot.json`; its name is retained for compatibility and it contains the complete master workbook snapshot. The parser validates required columns, unique IDs, evidence tiers, scores and regexes, and rejects formula cells or Excel errors rather than relying on stale cached values.

`Motif_library/context_rules.json` supplies explicit routing and operational requirements, not additional motif patterns. It also preserves the domain-count alternatives from the upstream detection pipeline. New motif IDs require a reviewed rule before execution. Review those rules when changing a motif's biological meaning or context.

- Talin/kindlin recognition entries are scanned on integrin-beta partners, not talin/kindlin sequences.
- Annotation-only actin, kindlin, Src, ILK and profilin entries produce separate reader hypotheses without scanning their regexes.
- Profilin polyproline and filamin-binding entries are exploratory partner screens. The latter is restricted to supplied integrin-beta candidates. Neither is treated as a motif encoded by its named reader.
- Fibronectin source categories are explicitly routed to `fibronectin_like_candidate`. Unsupplied families, such as vitronectin, are reported in library coverage.
- Collagen `[OP]` becomes encoded `P`, following the workbook's precursor instruction. FASTA `O` is not interpreted as hydroxyproline. Original and normalized regexes are both recorded.
- Exact workbook patterns and terminal anchors are retained. Y/F variants and broader old ELM patterns are not added silently.
- Identical motif-site observations from different workbook rows share a `site_group_id`. They are not independent evidence.
- The workbook's literature, citation status and notes are preserved as supplied. They are not a claim that every citation was independently bibliographically verified. Indirect/unresolved entries cannot yield high-confidence hypotheses. Known source issues are surfaced in `library_audit.tsv`.

## Confidence and functional hypotheses

Workbook `confidence_score` rates evidence and portability of a **library rule**. It is never copied into candidate confidence or interpreted as a probability. Candidate tiers follow transparent gates:

| Tier | Meaning |
|---|---|
| `high_confidence_hypothesis` | Eligible Include entry (Tier A/B, score at least 75, direct/supporting citation), supported domain-count architecture, and every configured context/evidence requirement supported. Still a hypothesis, not experimental validation. |
| `context_supported_candidate` | Eligible library entry with domain-count support and compatible supported/provisional localization, but additional required evidence is missing. |
| `sequence_only_candidate` | Eligible library entry matched; architecture/localization support is incomplete. |
| `review_only` | Low-confidence, Tier C, Review-only or indirect/unresolved entry. Cannot be promoted merely by a regex match. |
| `context_conflict` | A required criterion is contradicted, including incompatible topology/location. Retained for review rather than silently discarded. |

Every high-confidence hypothesis requires independent orthology, positional conservation, accessibility and partner-compatibility evidence in addition to motif-specific criteria. Some motifs also require phosphorylation, hydroxylation, structural presentation, force exposure, mature-chain location or membrane-proximal/distal assignment. The explicit requirements are in `context_rules.json` and every result lists supported and missing/conflicting criteria.

A high-scoring library does not guarantee high-confidence candidate results. An empty `high_confidence_hypotheses.tsv` is an honest result when required evidence is unavailable. `functional_hypotheses.tsv` contains context-supported and high-confidence hypotheses, with conditional biological interpretations and next-validation requirements. Reader hypotheses remain separate and conditional on family identity and partner context.

## Domain and topology evidence

Regions use one-based inclusive coordinates with columns `family`, `sequence_id`, `region`, `start`, `end`. Versionless Pfam accessions identify domain intervals. The retained verified domain TSV and provenance were imported with exact sequence checking; `import_domain_regions.py` can regenerate them from upstream retained domain CSVs and their exact search FASTAs.

DeepTMHMM `predicted_topologies.3line` imports require identical identifiers **and full sequences**. Lengths and topology symbols are validated. Mismatched or missing records are reported without transferring coordinates. A cytoplasmic tail is derived only from a single predicted TM helix immediately followed by terminal inside residues. `I`/`O` are side predictions, not experimentally established compartments; outside can include organelle lumen and is only provisional extracellular evidence. Signal-peptide or membrane-contained motifs conflict with the corresponding soluble-compartment hypothesis.

Domain counts support architecture but do not reproduce upstream homology, domain-order or classification rules. Missing domains are unknown rather than automatic family rejection. Intracellular motifs overlapping domains require resolution before high-confidence promotion, except configured domain-embedded/structural motifs such as VBS and extracellular ligands. Absence of domain overlap does not prove disorder or accessibility.

## Adding independent site evidence

Pass a TSV with `--evidence`. Required columns:

```text
family sequence_id sequence_sha256 motif_id start end criterion status source
```

Use tabs, exact sequence IDs, the SHA-256 of the normalized full protein sequence (available in `sequence_qc.tsv`), a workbook motif ID, one-based inclusive evidence intervals, and `supported` or `contradicted`. `criterion` must be a requirement listed for that motif, or `localization`, `family_architecture`, or `domain_overlap_resolved`. A nonempty source must identify the supporting analysis, alignment, reference or experiment. `Motif_library/site_evidence_template.tsv` provides the header without fabricated evidence.

Evidence applies only when its interval fully contains the motif hit. Incorrect hashes, unknown IDs, invalid coordinates and unsupported criteria fail validation. Contradictory evidence takes precedence. Do not mark conservation supported merely because a short word recurs in several candidate sequences. Do not label a predicted compartment, generic domain count or workbook score as independent orthology/functional evidence. Supplied sources remain user assertions; this pipeline validates their linkage and scope, not the underlying experiment.

## Outputs

- `REPORT.md`, `CURATED_FINDINGS.md`, `curated_summary.json`: run overview, confidence counts and evidence gaps.
- `motif_annotations.tsv`: unified compact sequence-hit table for enabled master records. `adhesion_motif_annotations.tsv`: detailed adhesion evidence, library scores, candidate tiers, conditional functions and next validation.
- `high_confidence_hypotheses.tsv`, `functional_hypotheses.tsv`: filtered evidence-based review views.
- `reader_hypotheses.tsv`: reader properties, never claimed as encoded motif matches.
- `library_audit.tsv`, `library_coverage.tsv`, `inspection_states.tsv`: normalization, routing, absent families and explicit detection states.
- `motif_prevalence.tsv`: deduplicated target versus other-candidate prevalence, excluding shared sequences. These are unmatched controls and no enrichment p-values are claimed.
- `sequence_qc.tsv`, `family_summary.tsv`, `domain_inspection.tsv`, `topology_audit.tsv`, `used_regions.tsv`: input and context checks.
- `discovered_words.tsv`, `discovered_sites.tsv`, `discovery_annotation.tsv`: exploratory exact-word discovery and original-position overlaps with curated motifs. Word recurrence is neither alignment-based conservation nor independent validation.
- Workbook/rule snapshots and `manifest.json`: source values, operational requirements, input/code hashes and parameters.
- `commands.json`, `external_stages.tsv`, family subfolders: optional MEME commands, statuses and deduplicated protein inputs. External failures return nonzero; no tool is silently installed.

Exploratory discovery uses 4â€“6 residue exact words, at least two distinct supporting sequences and 50% family prevalence by default. It filters low-entropy words, ranks descriptive prevalence differences and suppresses contained redundant words with identical support. Related proteins, isoforms, sequence length and composition remain confounders. STREME uses protein mode and skips small families by default. FIMO localization on training sequences is not independent validation, and external q-values are per run.

The confidence design follows the need for context, conservation and accessibility described by [SLiMSearch](https://pmc.ncbi.nlm.nih.gov/articles/PMC3125787/) and the [ELM resource paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC8728146/). These methodological references do not create an ELM-data dependency.
## Unified master workbook

`Master_motif_library` is authoritative. Supporting sheets are never scanned independently or used to resurrect missing master records. Only scan-enabled, valid, nonempty, non-reader, non-excluded patterns are scanned. Disabled legacy candidates remain in `legacy_review_catalogue.tsv`; no sequence hits are fabricated for them. Adhesion citation status is joined by exact motif ID from `Adhesion_curated`. Master scores, patterns and flags take precedence over supporting sheets.

`motif_annotations.tsv` is now the unified compact table for all enabled adhesion and ELM matches. Detailed adhesion evidence, sources and functional hypotheses are in `adhesion_motif_annotations.tsv`. Generic ELM matches are labeled `ELM_context_unverified`, irrespective of class curation or probability. The supplied generic context prose does not contain enough machine-readable class-specific taxonomy, PTM and compartment constraints to validate those sites. Domain overlap is reported for review; no generic ELM match is promoted to high confidence.

`master_library_audit.tsv` records every master ID and effective scan state. `generic_class_counts.tsv` includes zero-hit and disabled classes. `pathway_reference_associations.tsv` keeps KEGG/Reactome associations at class level, with no transfer to candidate proteins. `source_workbook_provenance.tsv` preserves the three supplied source checksums as recorded provenance, not independently verified source contents. The current workbook itself is hashed in `manifest.json`.

The Library_summary reports 41 excluded/flagged records. This is reported in `master_summary.json` and `REPORT.md`.

## PINCH and unresolved Parvin

SLIM037 annotates PINCH architecture without scanning a linear motif. SLIM038/039 screen paxillin-like partner sequences only and remain provisional; a match cannot establish Parvin identity. Proline-containing matches conflict with the helix requirement. `PINCH_Parvin_identity_review.tsv` reviews all PINCH, actin and alpha-actinin inputs without relabeling or creating Parvin FASTA records. Overlapping same-accession domain intervals are merged before counting. Five LIM domains and tandem CH domains are follow-up screens, not family identification. Independent family-specific homology/CDD, orthology and IPP context remain required.


## Expanded families and predicted sequence context

The default FASTA library contains 26 families, including Csk, Crk, Grb2, Nck, Shc and Syk. Their Pfam count and order screens are mirrored from `01_detection/config/family_rules.json`; they are architecture triage only and do not establish homology or family identity. The current master workbook has no curated adhesion SLiM rows for Csk, Crk, Grb2, Nck, Shc and Syk families. Generic ELM class matches remain `ELM_context_unverified`. No family-specific motif is invented.

When the bundled FASTAs are used, the pipeline imports `deeptmhmm_predicted_topologies.3line`, `DeepLoc_adhesome_candidate_list.csv`, `SignalP_adhesome_candidate_list.txt` and `NCBI_CDD_adhesome_hitdata.txt` from `../02_topology_subcellular localization/adhesome_candidate` by default. For another FASTA set, pass `--context-dir` explicitly. `--topology-file` selects a single 3line file. DeepLoc, SignalP and CDD records are applied only when the same ID has an exact full-sequence match in the topology file. CDD imports specific hits with E-value at most 1e-3 and de-duplicates identical accession/interval calls; non-specific and superfamily calls are not transferred.

`sequence_context.tsv` reports imported predictions and verification per family record. `CDD_specific_hits.tsv` gives specific domain locations and source scores. CDD domain overlaps and predicted topology/compartment information appear beside hits in `motif_annotations.tsv`; curated hits also include SignalP cleavage and DeepLoc information in `adhesion_motif_annotations.tsv`. Concordant DeepLoc and SignalP can make an otherwise unknown localization provisional. Predicted evidence alone never satisfies orthology, conservation, accessibility, partner compatibility or high-confidence gates. 
