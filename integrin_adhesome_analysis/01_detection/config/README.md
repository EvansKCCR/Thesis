# Canonical family rules — workflow version 1.0

`family_rules.json` is the single active 96-entry catalogue for detection and the
final reference-gated audit. It contains every family from the supplied older
catalogues; the newest complete entry takes precedence for duplicate names.
The resource retains schema 5.0 as a format identifier. Pipeline release/version
reports remain 1.0. `resource_resolution.json` records source filenames, hashes
and consolidation policy. Superseded JSONs and their legacy builder were removed.

Entries distinguish broad screens, families and subfamilies. Shared domains,
HOG membership, motifs and topology cannot establish a resolving paralog identity.
Specific assignments requiring identity evidence use one independently supported
route: curated reference accession, concordant phylogeny, concordant reciprocal
reference identity or a calibrated competing-profile winner. Without that evidence,
architecture-compatible candidates remain provisional detection hypotheses.

The updated Pfam library contains 176 models: the 174 manifest entries plus
PF00046 and PF03535 required by the current rules. Obsolete inherited motif
expressions were removed. Current motif definitions come from stage 03's
`Motif_library/adhesion_core.json`; full motif annotation occurs during resume.

Validate resource consistency from the workflow root:

```bash
python validate_resources.py
```

The audit uses `05_conservative_adhesome_family_assignment_audit/integrated_adhesome_pipeline.py`.
Its default rule path is this canonical file, and the integrated wrapper passes
that path explicitly. See the audit module README for output states and benchmarks.
The supplied rule thresholds remain uncalibrated proposals; fixture tests verify
software behavior, not empirical classification accuracy.
