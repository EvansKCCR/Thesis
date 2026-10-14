# Detection module — version 1.0

Run from `Thesis/integrin_adhesome_analysis`:

```bash
python adhesome_workflow.py start --targets-config targets.json --input-dir input/targets --out runs/handoff --threads 4
```

Detection prepares gene-level target proteomes, scans Pfam domains, searches the
curated reference set, and exports supported/provisional family hypotheses.
Topology, subcellular localization and phylogeny belong to the resumed workflow.
Detection evidence files containing those results are rejected with a resume
instruction. External predictor results are imported through `resume`; predictor
services remain external.

The active runner is `scripts/run_detection.py`. The active library is
`databases/pfam/adhesome_v1.0.hmm`, containing all 176 accessions in the updated
`config/adhesome_pfam_manifest.csv`. The 230-row
`config/adhesome_component_crosswalk.csv` is validated against that manifest.
The canonical `config/family_rules.json` contains the current 96 hierarchical
entries and is shared with the reference-gated audit. References are in
`references/adhesome_v1.0`. The previous runner name remains
as a compatibility entry point to version 1.0. Historical seed/model filenames
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

The bundled library is ready for detection. To reproduce its extraction:

```bash
python 01_detection/build_adhesome_pfam_library.py --pfam-hmm 01_detection/databases/pfam/source_v1.0/combined_source.hmm
```

To retrieve models after a future manifest expansion:

```bash
python 01_detection/build_adhesome_pfam_library.py --download
```

Models come from the documented [InterPro Pfam HMM API](https://interpro-documentation.readthedocs.io/en/latest/api.html).
`source_v1.0/download_provenance.json` records model accession versions, URLs,
hashes and retrieval time. Library metadata verifies the manifest, crosswalk,
combined library and GA/fallback splits. Detection runs `hmmpress` on these splits
in the execution environment. HMMER and BLAST+ remain external dependencies.

See [WORKFLOW.md](../../README/WORKFLOW.md) for setup and resume instructions.
The Conda environment name remains `adhesome-workflow-v1` for existing installations;
all pipeline release/version reports are 1.0. Run detection fixture tests from
this directory with `python -B -m unittest discover -s tests`.
