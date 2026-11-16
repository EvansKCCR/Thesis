# CDD domain sequence library v1

Extracted from NCBI CD-Search results in the returned calibration and validation evidence, restricted to labelled positive adhesome controls. Coordinates are 1-based inclusive and fragments are cut from the exact prepared protein FASTAs used by each run.

`cdd_domain_sequences.faa` stores unique amino-acid fragments. Identical fragments share a FASTA ID. `cdd_domain_metadata.tsv` preserves each source hit, including source protein ID, panel origin, accession/PSSM ID, family annotation, CDD hit class, superfamily, coordinates, E-value and bit score. Thus one deduplicated sequence may have multiple provenance rows.

NCBI `specific` maps to `specific`, `non-specific` remains `non-specific`, and NCBI `superfamily` maps to the broader `broad` tier. These labels should remain separate features; they are not interchangeable evidence. `manifest.json` records extraction rules and counts. Its classifier reference-bank entry describes the calibration-only subset: 20,003 unique fragments across 41,798 traceable calibration hit records and 5,570 accessions. All 20,003 fragments are retained, even when the source positive control has no mapped trusted family; those fragments still contribute CDD class/accession/annotation evidence without being assigned an invented family label.

To rebuild from the saved validation runs and returned CDD tables, run from the main workflow bundle:

```bash
python -B resource_library/adhesome_ml/build_cdd_domain_library.py \
  --calibration-run /path/to/runs/isolated_validation/calibration \
  --validation-run /path/to/runs/isolated_validation/validation \
  --calibration-cdd /path/to/returned/isolated_calibration/NCBI_CDD_adhesome_hitdata.txt \
  --validation-cdd /path/to/returned/isolated_validation/NCBI_CDD_adhesome_hitdata.txt \
  --out resource_library/adhesome_ml/cdd_domain_library_v1
```

The command refuses an existing output directory. The v5.8 classifier verifies each selected calibration fragment against its recorded SHA-256 digest, then builds its frozen reference bank from calibration-positive records only. Validation-origin fragments never enter the bank. Validation and target proteins are searched against the same frozen calibration bank. Query and subject coverage use inclusive BLAST HSP coordinate spans divided by full query and subject lengths, capped at 1.0. Calibration records exclude same-protein references.

Family-resolution searches use the configurable E-value/coverage profiles in
`config/ml_classifier/classifier.json` and can retain CDD-fragment matches through
E-value `1e-2`. The historical `cdd:*` features used by the membership model remain
limited to E-value `1e-3`; the new `family_similarity:*` evidence is used only for
family assignment. Thus changing `--family-similarity-profile` does not alter the
membership probability or its cutoff. Optional per-rule `sequence_similarity_policy` entries
in `config/annotation/family_identity_rules.json` add family-specific E-value and
coverage acceptance routes without relying on exact sequence equality. These routes
do not replace the selected global profile. Prediction rows retain the supporting
fragment ID, CDD annotation/accession, identity, bit score, E-value, coverage and the
threshold source for each similarity match reported by the gate.

The `adhesome_family_rule_crosswalk.xlsx` workbook is an additional source for
family-discriminating classifier features. The loader pairs CDD accessions in reviewed
P1-P3 rows and configured family-pattern rows with the linked rules in
`config/annotation/family_identity_rules.json`. A feature activates from sequence
evidence only when the CDD accession has a calibration-bank HSP; companion features
measure how much of the linked rule's domain architecture is present and whether the
full rule constraints are satisfied. Proposed family/class text, same-protein family
context, calibration counts and the disabled alias-proposal sheet are not used as
query features. The crosswalk and family-rule hashes are saved with each v5.8 model.

Classifier v5.8 calibrates membership and family outputs separately using sigmoid
maps fit only on group-disjoint out-of-fold calibration predictions. CD-HIT groups
sequences at 98% identity and 90% coverage of both proteins for folds; exact-sequence
SHA-256 grouping is the recorded fallback. These grouping IDs never enter model features.
Raw forest votes remain diagnostic fields. Family probabilities are conditional on adhesome
association and are reported only for families meeting the configured minimum of 10
independent positive sequence-cluster groups and 5 negative groups; other family probabilities are
marked unavailable rather than replaced with raw votes. Validation is used only for
reporting Brier score, log loss, reliability bins/plots, and classification metrics.
The default membership, specific-family, margin, and broad-parent gates now operate on
the calibrated score scales. These probabilities reflect the labelled control
distribution and do not establish prevalence or universal biological certainty.
Reports preserve previous raw-vote versus calibrated-policy calls on the same validation
cohort. Parent-family probabilities sum only represented calibrated descendants and
record any unrepresented descendant classes per prediction.
