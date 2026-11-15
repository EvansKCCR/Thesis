# CDD domain sequence library v1

Extracted from NCBI CD-Search results in the returned calibration and validation evidence, restricted to labelled positive adhesome controls. Coordinates are 1-based inclusive and fragments are cut from the exact prepared protein FASTAs used by each run.

`cdd_domain_sequences.faa` stores unique amino-acid fragments. Identical fragments share a FASTA ID. `cdd_domain_metadata.tsv` preserves each source hit, including source protein ID, panel origin, accession/PSSM ID, family annotation, CDD hit class, superfamily, coordinates, E-value and bit score. Thus one deduplicated sequence may have multiple provenance rows.

NCBI `specific` maps to `specific`, `non-specific` remains `non-specific`, and NCBI `superfamily` maps to the broader `broad` tier. These labels should remain separate features; they are not interchangeable evidence. `manifest.json` records extraction rules and counts.

To rebuild from the saved validation runs and returned CDD tables, run from the main workflow bundle:

```bash
python -B resource_library/adhesome_ml/build_cdd_domain_library.py \
  --calibration-run /path/to/runs/isolated_validation/calibration \
  --validation-run /path/to/runs/isolated_validation/validation \
  --calibration-cdd /path/to/returned/isolated_calibration/NCBI_CDD_adhesome_hitdata.txt \
  --validation-cdd /path/to/returned/isolated_validation/NCBI_CDD_adhesome_hitdata.txt \
  --out resource_library/adhesome_ml/cdd_domain_library_v1
```

The command refuses an existing output directory. The v5 classifier uses calibration-origin fragments from this FASTA as its reference bank. Validation-origin fragments are never used as model references. It searches validation and target protein sequences against the frozen calibration bank and adds hit class, family, accession, similarity, query/subject coverage, family hit count, copy number and domain order features. Calibration records are scored with same-protein references excluded.
