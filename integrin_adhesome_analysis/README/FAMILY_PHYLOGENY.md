# Unified family phylogenies

The updated scripts can union all candidate-associated orthogroups for each
family into one alignment/tree project. The existing per-orthogroup mode remains
the default. This option does not rerun or merge the OrthoFinder inference itself.

## Test with your existing completed orthology

Copy the updated `04_orthology_phylogenetic_classification/scripts/phylogeny.py`
from this package to the corresponding location in `Test_pipeline_v1`.
For integrated resume/continue support, also copy the updated root
`adhesome_workflow.py`. No existing data or run directories need to be copied or
modified. These instructions use your actual run name: `runs/annotation`.

From the `Test_pipeline_v1` root in your activated environment:

```bash
python 04_orthology_phylogenetic_classification/scripts/phylogeny.py prepare --hogs runs/annotation/downstream/orthology/candidates/candidate_HOGs.tsv --orthologues runs/annotation/downstream/orthology/candidates/direct_orthologues.tsv --duplications runs/annotation/downstream/orthology/candidates/candidate_duplications.tsv --proteomes runs/annotation/downstream/panel/proteomes --species-list runs/annotation/downstream/panel/species.json --family-map runs/handoff/hypotheses.tsv --family-column motif_family --candidate-report runs/annotation/downstream/orthology/candidates/candidate_report.tsv --output runs/family_phylogeny
python 04_orthology_phylogenetic_classification/scripts/phylogeny.py run --project runs/family_phylogeny --threads 2 --iqtree-threads 2
```

If duplication data are unavailable, replace `--duplications PATH` with
`--allow-missing-duplications`; missing evidence remains explicitly unknown.
Use a new `--output` folder. To test one family first, append e.g.
`--family fibronectin_like_candidate` to the preparation command. Repeat
`--family` to select several. The run command then runs the prepared selections.

## What is grouped

The family-map TSV needs `protein_id` and the column chosen by `--family-column`.
The workflow's `motif_family` column normalizes fibronectin screening labels to
`fibronectin_like_candidate`. To keep original detection family labels separate,
use `--family-column family` instead. A curated two-column mapping TSV can also
be supplied when you want to control/review family assignments before inference.

Each candidate ID is translated through `candidate_report.tsv` to its prepared
gene ID. Every source group containing a candidate associated with the family
is included, with **all its target and reference members**. Protein IDs are
unioned: each ID appears only once in a family tree even when it occurs in
several source groups. Different IDs with identical sequences remain separate
biological tips. Sequences are aligned, not concatenated.

Candidate families are hypotheses. An ambiguous candidate/group can contribute
to multiple family projects; no winning family is silently selected. Family
membership of every included reference is not independently established by the
merge. Full-length alignment QC remains active: combining highly divergent or
multidomain groups can still require a reviewed homologous-region strategy.
The merge itself does not relax QC or prove whole-protein orthology.

## Outputs and provenance

```text
runs/family_phylogeny/
  family_orthogroups.tsv          # family -> source groups and triggering candidates
  family_candidate_coverage.tsv # included candidates and those absent from selected groups
  groups.tsv
  manifest.json                 # grouping mode and hashed family/mapping inputs
  family_<name>/
    members.faa
    tips.tsv                    # original ID, candidate flag, family and source groups
    orthology_edges.tsv
    duplication_evidence.tsv
    inference/                  # alignment, QC, tree, support and command logs
```

Shared orthogroups are visible in `family_orthogroups.tsv`, and per-tip source
memberships are retained in `tips.tsv`. A matched-candidate flag refers to that
family's candidates, not every candidate from an included group.

## Automatic workflows and recovery

For a new automatic run, add `--phylogeny-grouping family` to `resume` or
`continue`. Orthology and the independent audit still run as before; only the
phylogeny preparation grouping changes. Grouping is saved for checkpoint recovery.
`continue` performs a new downstream run; use the direct prepare/run commands
above to reuse existing orthology without repeating it.

To retry failed family trees in an unsealed standalone project:

```bash
python 04_orthology_phylogenetic_classification/scripts/phylogeny.py run --project runs/family_phylogeny --threads 2 --iqtree-threads 2 --retry-failed
```

Do not change grouping inside an existing prepared project; create a new output
folder. No biological analysis has been run for this update; tests used fixtures.