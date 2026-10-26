# Components FASTA preparation package
## Retrieval script
Requires Python 3 and requests, plus internet access.

```bash
python build_calibration_set.py --components components.csv --out output --background-count 232 --seed 20261014
```

Expected outputs: candidates.fasta, hard_negatives.fasta, background_negatives.fasta, components_with_negatives.fasta, mapping and sequence manifests, issues.csv, summary.json and retrieved source records.

## Explicit design assumptions
- uman reference sequences are used as the initial scope. This is not a schistosome-specific negative set.
- Hard-negative gene proposals are challenges for the specified target family, not claims that those proteins have no adhesion-related functions. Only uniquely resolved reviewed human entries are accepted.
- Background sampling uses the reviewed, independently enzyme-annotated subset of human reference proteome UP000005640. It excludes candidate and hard-negative accessions, exact sequence duplicates, candidate gene aliases, and entries sharing candidate Pfam identifiers. The seed is fixed. This is a conservative enzyme-enriched background, not an unbiased proteome-wide sample.
- No architecture analysis, homology search, train/test split, or experimental validation is performed by this script.