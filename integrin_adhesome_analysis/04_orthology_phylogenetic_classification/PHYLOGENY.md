# Optional per-HOG phylogenetic review

From the package root, after candidate orthology extraction:

```bash
python 04_orthology_phylogenetic_classification/scripts/phylogeny.py prepare --hogs runs/orthology/candidates/candidate_HOGs.tsv --orthologues runs/orthology/candidates/direct_orthologues.tsv --duplications runs/orthology/candidates/candidate_duplications.tsv --proteomes runs/panel/proteomes --species-list runs/panel/species.json --output runs/phylogeny
python 04_orthology_phylogenetic_classification/scripts/phylogeny.py run --project runs/phylogeny --threads 8
```

Use actual duplication-event output; if it is unavailable, resolve the missing
input before preparation. Each HOG is aligned independently with MAFFT, trimmed
with trimAl and inferred with IQ-TREE. Review alignment QC, homologous coverage,
branch support and duplication evidence. Root only with a justified outgroup.
Trees are not automatically converted into final family-assignment scores.
Use each subcommand's `--help` for options. No trees or prior projects are bundled.