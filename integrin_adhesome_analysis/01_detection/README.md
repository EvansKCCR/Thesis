# Comparative integrin-adhesome detection

The active v5 pipeline uses the 96-model Pfam library, motif-aware comparative family rules, 71 forward references and 142 reciprocal references. Syk, Csk, Nck, Grb2, Crk and Shc have separate family rules; Src remains included. See [v5 expansion and limitations](docs/SIGNALLING_EXPANSION_V5.md).Run commands from `01_detection` with Python 3.9+; full searches require HMMER 3 and NCBI BLAST+ (for example in WSL `structbio`).

## Run

```bash
python scripts/run_domain_v5.py --cpu 8
python -B -m unittest discover -s tests -v
```

The wrapper recursively discovers the three species in `proteome_fasta_files` and runs preparation, domain discovery, homology and classification. New results default to `results/species_reference_v5`; use `--out` to choose another directory. It uses the final rules, library and reference set without depending on older runs.

Rebuild resources offline into separate check directories:

```bash
python build_reference_sets.py --offline --strict-reviewed-curated --out results/reference_rebuild_check
python build_adhesome_pfam_library.py --pfam-hmm databases/pfam/targeted_v4/combined_source.hmm --out results/library_rebuild_check/adhesome_v4.hmm
```

For standalone parsing, supply `--model-metadata databases/pfam/adhesome_v4.hmm.metadata.json` and the active `--library-sheet config/adhesome_pfam_manifest.csv`. The pipeline handles HMM pressing automatically.

## Retained layout

- Root Python modules: final preparation, search, parsing, classification and reference/library builders.
- `config/`: active rules, manifest, proteome configuration and input templates.
- `proteome_fasta_files/`: original species proteomes.
- `references/schistosome_v5/`, `references/seeds/schistosome_v5/`: paired FASTAs/catalogues, pinned seeds and provenance.
- `references/cache/uniprot/`: raw records needed for offline rebuilding; `references/motifs/`: supplied motif resources.
- `databases/pfam/`: final library and indexes. `targeted_v3/` and `targeted_v4/` hold source records used by the final composite library, not alternative pipeline versions. Only the final combined source is retained.
- `results/species_reference_v4/`: latest species results, with prepared inputs inside each species directory.
- `results/reference_calibration_v4/`: latest reference-panel calibration.
- `scripts/`: full-run and reference-calibration utilities; `tests/`: current regression suite.
- `docs/`: final methods and cleanup inventory.

See [latest results](results/species_reference_v4/DOMAIN_RESULTS.md), [methods and limitations](docs/DOMAIN_RECALIBRATION_V4.md), and [reference provenance](references/README.md). Confidence labels are evidence tiers, not probabilities or proof of experimental adhesion function. Calibration uses the development panel, not independent validation.

Cleanup preserves latest result contents and downstream paths. Historical command logs retain original absolute paths; `docs/cleanup_manifest.json` records removals, checksums and prepared-input relocation. No older pipeline versions or result trees are required for a new full run.

V4 reference resources and `config/family_rules_v4.json` are retained for exact interpretation of existing results. V5 uses the unchanged v4 HMM library.
