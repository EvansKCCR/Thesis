# Detection evidence status (workflow 1.0)

Detection retains domain-based hypotheses for external follow-up even when their
identity remains unresolved. Each row is a protein/family hypothesis: a protein
can have several rows. Counts of rows are not counts of distinct sequences.

## Reading the checks

| Field/value | Meaning |
| --- | --- |
| `order_check=not_applicable` | The rule describes one domain type or repeat counts; there is no distinct domain-type order to test. Domain counts are still checked. |
| `order_check=not_configured` | A multidomain rule has no configured order constraint. This is a configuration limitation, not a passing check. |
| `order_check=incomplete_or_reordered` | A configured domain order was not recovered; missing modules and rearrangements are not distinguished by this field. |
| `similarity=supported` | Evidence contains a qualifying hit for this exact rule family. It does not establish a specific paralog or adhesome membership. |
| `similarity=no_qualifying_hit` | The paired search completed for this proteome, references exist for this exact family, but no qualifying hit was retained. It does not prove biological absence. |
| `similarity=not_tested`, reason `no_family_reference` | The search panel lacks references labelled for this exact family/subfamily. |
| `similarity=not_tested`, reason `search_completion_not_documented` | No verified completion record accompanies the evidence; missing rows cannot be interpreted as a completed search. |
| `reciprocal=not_tested`, reason `no_qualifying_forward_hit_for_reciprocal_assessment` | A paired search ran, but no forward family hit qualified for the reciprocal test. |
| `reciprocal=not_supported` | A qualifying forward hit exists but reciprocal family agreement was not recovered. |

`similarity_reason`, `reciprocal_reason`, `order_check_reason`, `homology_check`
and `reference_count` explain availability separately from evidence support.
`related_reference_families`, `related_reference_count`, `related_similarity`,
`related_reciprocal` and `related_evidence_sources` expose parent/legacy-family
context. These fields do not contribute positive subfamily evidence or populate
independent identity fields. Similarity to Syk-family references, for example,
does not resolve SYK versus ZAP70.

The classifier checks the run-local paired-search manifest, target FASTA hash
and reference snapshot hashes before interpreting a missing hit. New search runs
also checksum the BLAST artifacts; older completed manifests remain readable.
Standalone evidence CSV files can supply positive/negative evidence but cannot
document coverage of families without hits.

Each target has `candidates/reference_coverage.csv` covering every configured
rule, including rules with no candidate rows. A new handoff combines these as
`reference_coverage.tsv`. Families lacking direct references are flagged
`add_family_specific_references`. Adding reference proteomes at the downstream
orthology calibration stage does not retroactively add detection BLAST references:
the detection reference set must also be expanded and searched if detection-level
similarity/RBH evidence is required. A cached refresh cannot create missing evidence.

## Refresh existing detection without rerunning searches

From the project root in WSL (the directory containing `adhesome_workflow.py`):

```bash
python detection/scripts/refresh_detection_evidence.py --detection runs/handoff_retry/detection --out runs/detection_evidence_review
python adhesome_workflow.py start --detection runs/detection_evidence_review --targets-config targets.json --out runs/handoff_evidence_review
```

Both output directories must be new. The first command validates cached FASTA,
domain and homology input hashes, copies prepared/domain/homology directories,
and reclassifies all completed targets with the current rules and libraries.
It requires enough space for a second copy of those cached files. It launches no
HMMER, BLAST, OrthoFinder or phylogeny jobs. Original detection and handoff outputs
remain unchanged. Use the new handoff for subsequent external-evidence return and
resume. The second command exports the new reports and runs only domain-region
import/verification. If inputs have changed or detection is incomplete, finish a
consistent detection run before attempting this refresh.

For installation, update `adhesome_workflow.py`, `detection/adhesome_pipeline.py`,
`detection/adhesome_classification.py`, and add `detection/evidence_coverage.py`
and `detection/scripts/refresh_detection_evidence.py` from this release.

## Observed validation run

Read-only inspection of `runs/handoff_ready/hypotheses.tsv` on 2026-10-15 found
6,735 exported hypothesis rows: 573 with supported similarity and 6,162 labelled
`not_tested` by the old reporting. Of those 6,162, 5,533 lacked exact-family
references and 629 had exact-family references but no qualifying hit. At that time
the search catalogue covered 26 family labels while the rules contained 96 entries.
These are reporting/coverage counts, not empirical sensitivity estimates. Real
analysis has not been rerun as part of this code update.

The later [crosswalk expansion](COMPONENT_REFERENCE_PANEL.md) covers all 230
components and raises explicit family-label coverage to 48. Use the documented
`--reference-set` refresh option to run BLAST against those new references while
reusing existing domain scans.
