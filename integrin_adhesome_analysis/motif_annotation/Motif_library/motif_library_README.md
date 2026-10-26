# Evidence-based motif master library

## Purpose
Evidence-supported, nonredundant motif definitions for candidate production scanning; not a validated functional-hit classifier.

## Default scan scope
Adhesion_core and adhesion_core.json. General ELM regulatory motifs are available separately in the full master. Reader motifs are not diagnostic family signatures.

## Inclusion policy
ELM classes with >=1 curated instance in supplied snapshot and valid nonempty regex; additional exact motifs only with inspected primary experimental support. Zero-instance exclusion is a conservative policy, not a claim that the class is invalid.

## Deduplication
One scan per exact regex; aliases and reader annotations merged. Abbreviated direct-mapped rules replaced by canonical ELM patterns. Distinct overlapping classes retained because their reader/context specificity differs.

## Excluded content
Legacy candidates, provisional parvin orientation rules, generic PxxP/polyproline/PEST/basic patterns, domain-architecture records and unsupported generalized regexes. Exclusion_log retains IDs and reasons, not active scan rules.

## Evidence versus prediction
Evidence confidence refers to the motif definition. All new regex hits remain candidates. PTMs, reader compatibility, topology, accessibility, taxonomy and conservation need external assessment. No hit-level confidence percentages assigned.

## Species transfer
No schistosome-specific validation established by this rebuild. A supported host motif class does not prove a conserved parasite interaction.

## PINCH and parvin
PINCH LIM architecture is not a SLiM scan rule. Parvin forward/reverse candidate rules excluded. Canonical paxillin LD class is not a parvin-specific classifier.

## Regex semantics
Python re; uppercase protein sequence; scan original full-length sequence to preserve anchors; detect overlapping start positions. $ means true protein C terminus, not fragment boundary. Do not reverse protein sequences.

## Unknown residues
Reject hits spanning X, B, Z, J, gaps or stop symbols. O is pyrrolysine in standard FASTA, not a hydroxyproline annotation. No automatic O-to-P conversion.

## Collagen
Canonical ELM precursor regex retained. Proline hydroxylation and triple-helical context cannot be inferred from sequence matching.

## Live evidence coverage
Primary sources for added rules and ELM FAT-LD, VBS and EVH1 class pages inspected. Remaining ELM definitions inherit the uploaded snapshot evidence and were not individually live-revalidated.

## Release limitation
Rebuilt from uploaded ELM snapshot, not a complete refresh of the current online ELM release. No sensitivity/specificity benchmark or empirical false-discovery calibration performed.

## Citation correction
SLIM019 PMID 2578479 resolves to a Campylobacter antibody paper, not an RGD study. Removed as evidence; RGD consolidated into supported ELM class.

## Source reconciliation
Source master contains 691 rows, curated sheet 39 and ELM sheet 353. Source summary claims 694 master rows. Rebuild uses source-sheet contents, not stale summary totals.

## Version
2026-10-14 snapshot-rebuild v1.0

## Source SHA256
6fbb145ecdedd69143f754438b31d7f10f107c400e938724f4944c2489f967c0

## Method reference
https://pmc.ncbi.nlm.nih.gov/articles/PMC10927331/

## Rebuild counts
- Original master records: 691
- ELM source classes: 353
- Curated source entries: 39
- Production unique scan rules: 347
- Adhesion-associated scan rules: 53
- General regulatory scan rules: 294
- Retained ELM classes: 341
- Primary-backed additional rules: 6
- Legacy candidate rows excluded: 302
- ELM zero-instance classes excluded: 12
