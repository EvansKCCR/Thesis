# Canonical family resources — workflow version 1.0

`family_rules.json` is the current 222-entry discovery catalogue. The matching
Stage 2 resources are `annotation/family_rules.json` and
`annotation/family_identity_rules.json`, with 197 identity rules. The main
`calibrated_v1` calling algorithm and thresholds remain unchanged. Resource
schema numbers are separate from pipeline version 1.0.

`adhesome_pfam_manifest.csv` contains 214 unique Pfam accessions. The 38 accessions
missing from the updated family-rule coverage were appended on 26 October 2026;
the original manifest rows were retained. The active HMM library has 214 GA models
and an empty fallback split. Existing model versions and GA thresholds are
unchanged. Crosswalk content and its 230 component records were not modified.

[The library expansion receipt](../../maintenance/pfam_expansion_20261026.json)
records added model versions/URLs/hashes and the protected analytical-file hashes.
Download provenance is in `../databases/pfam/source_v1.0/download_provenance.json`.
Current motif definitions are in `motif_annotation/Motif_library/adhesion_core.json`;
full motif annotation occurs during resume.

Read-only checks from the workflow root:

```bash
python -B maintenance/expand_detection_pfam_library.py --check-only
python -B tools/validate_resources.py --allow-stale-calibration --audit-resources-only
```

The second command also reports the separate, older component-reference label
provenance and the stale numerical probability model. It does not relabel
references, fit models or run analysis. HMM model coverage is now complete for
the configured discovery rules; empirical recruitment performance requires a
future user-run scan. Existing scan/audit outputs are unchanged.
