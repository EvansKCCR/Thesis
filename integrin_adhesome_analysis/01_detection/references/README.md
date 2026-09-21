# Final reference resources

`schistosome_v4/` contains 55 forward-search sequences (26 curated seeds, 20 computational spiralian anchors and nine supplemental family references) and a targeted 126-sequence reciprocal panel. Active accession-pinned seeds are in `seeds/schistosome_v4/`. The original motif library is retained in `motifs/`.

Find a protein by `reference_id` in `schistosome_v4/search_catalogue.csv` or `metazoan_catalogue.csv`, then retrieve that identifier from its paired FASTA. Catalogues record accession, taxon, sequence version, sequence hash and evidence status. Individual FASTAs, experimental annotations and `reference_lock.json` retain exact provenance. Required raw UniProt responses remain in `cache/uniprot/` for offline replay.

```bash
python build_reference_sets.py --offline --strict-reviewed-curated --out results/reference_rebuild_check
```

Reviewed status and experimental annotations do not automatically verify adhesome membership or identify the exact isoform used in every experiment. Computational anchors are not experimentally demonstrated adhesome proteins. The reciprocal panel is targeted and human-weighted; `non_adhesome` denotes a competing family, not proof that a protein never participates in adhesion. The unresolved historical unc-60 seed is excluded from the final pinned inputs.

Each homology result retains its own reference snapshot, matched sequences, catalogues and checksums. Original paths in locks/logs are provenance from the time of execution; the project now resides under `01_detection`. See the [final methods](../docs/DOMAIN_RECALIBRATION_V4.md) and [cleanup manifest](../docs/cleanup_manifest.json).
