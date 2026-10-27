# Detection module â€” version 1.0

Run from `Thesis/integrin_adhesome_analysis`:

```bash
python adhesome_workflow.py start --targets-config targets.json --input-dir inputs/targets --out runs/handoff --threads 4
```

Detection prepares gene-level target proteomes, scans Pfam domains, searches the
curated reference set, and exports supported/provisional family hypotheses.
Topology, subcellular localization and phylogeny belong to the resumed workflow.
Detection evidence files containing those results are rejected with a resume
instruction. External predictor results are imported through `resume`; predictor
services remain external.

The active runner is `scripts/run_detection.py`. The active library is
`databases/pfam/adhesome_v1.0.hmm`, containing all 214 accessions in the updated
`config/adhesome_pfam_manifest.csv`. The 230-row
`config/adhesome_component_crosswalk.csv` is validated against that manifest.
The canonical `config/family_rules.json` contains the 222 Stage 1 hierarchical
entries. Stage 2 defaults to the expanded family and identity rules in
`annotation/` (222 family configurations / 197 identity rules). The main
`calibrated_v1` algorithm, scoring gates and cutoffs are unchanged. Historical
library comparisons use saved run configurations. References are in
`references/adhesome_v1.0`, now covering all 230 crosswalk components with 354
forward and 425 reciprocal reference records. See
[the expanded panel guide](../README/COMPONENT_REFERENCE_PANEL.md) to refresh
BLAST evidence while reusing existing domain scans. The unused previous runner
alias has been removed. Historical seed/model filenames
identify source provenance rather than a different active pipeline version.

Each target's `candidates/` directory contains:

- `candidates.csv`: existing family evidence and crosswalk annotations.
- `component_hypotheses.csv`: possible component matches, reference accessions,
  functional categories, matched domains, signature completeness and similarity.
- `component_crosswalk.json`: crosswalk hash and coverage information.

`crosswalk_possible_components` records all domain-overlap possibilities;
`crosswalk_reference_components` records matches to reference accessions supported
by similarity. Shared domains and human gene symbols do not establish target
orthology, component identity or experimental adhesome membership. Candidate
status uses the current family rules, homology and required independent identity
evidence. New domains
without a family rule are retained as unresolved domain candidates for review;
adding the crosswalk does not invent 230 validated family rules. Components without
reference Pfams are recorded in the coverage report. Crosswalk matches for exported
candidate proteins are carried into `runs/handoff/component_hypotheses.tsv`.

Positive component-reference hits with no retained family hypothesis are exported
as provisional `component_reference_candidate` records. Their family identity stays
unassigned, and resume accepts their external evidence. These heterogeneous records
are excluded from pooled family motif discovery. `crosswalk_reciprocal_components`
and `reference_reciprocal` retain panel-relative component reciprocity separately.

The library expansion on 26 October 2026 adds the 38 models required by the
current family catalogue, taking the library from 176 to **214 models**. All
214 have GA thresholds; the fallback split is empty. Every original model record,
version and GA threshold is preserved. The manifest, combined source, GA/fallback
splits, checksums and per-model download provenance have been refreshed. This
resource update did not run `hmmpress`, detection, annotation, orthology or trees.

The [expansion receipt](../maintenance/pfam_expansion_20261026.json) lists all 38
added accessions, model versions, source URLs and hashes. Existing scans/results
retain their original coverage; new models are available to subsequent scans.
Component-reference label provenance remains a separate outstanding refresh and
was not rewritten by this HMM-only update.

Check HMM rule coverage without a scan:

```bash
python -B maintenance/expand_detection_pfam_library.py --check-only
```

To reproduce extraction from the supplied source:

```bash
python detection/build_adhesome_pfam_library.py --pfam-hmm detection/databases/pfam/source_v1.0/combined_source.hmm
```

To retrieve models after a future manual manifest expansion:

```bash
python detection/build_adhesome_pfam_library.py --download
```

Models come from the documented [InterPro Pfam HMM API](https://interpro-documentation.readthedocs.io/en/latest/api.html).
`source_v1.0/download_provenance.json` records model accession versions, URLs,
hashes and retrieval time. Library metadata verifies the manifest, crosswalk,
combined library and GA/fallback splits. Detection runs `hmmpress` on these splits
in the execution environment. HMMER and BLAST+ remain external dependencies.

See [WORKFLOW.md](../README/WORKFLOW.md) for setup and resume instructions.
The Conda environment name remains `adhesome-workflow-v1` for existing installations;
all pipeline release/version reports are 1.0. Run detection fixture tests from
this directory with `python -B -m unittest discover -s tests`.

For a future rule-driven expansion while preserving cached model versions:

```bash
python -B maintenance/expand_detection_pfam_library.py
```

This maintenance command prepares library resources only, preserves the existing
analytical files, and creates no backups. Detection creates any HMMER indexes
when the user later starts a scan.
