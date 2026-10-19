# Audited family phylogenies

## Resume the stopped validation project locally

The existing validation project is an orthogroup run. Resume its original
prepared membership in place; preparation and audit are not repeated:

```bash
cd /mnt/c/Users/eaasa/Desktop/Evans_PHD_Thesis/pipeline_validation/integrin_adhesome_analysis
conda activate adhesome-workflow-v1
python 04_orthology_phylogenetic_classification/scripts/phylogeny.py run \
  --project runs/annotation/downstream/phylogeny \
  --resume --threads 4 --iqtree-threads 2
```

`--resume` verifies and skips completed groups, reuses valid saved alignments,
and invokes IQ-TREE with its existing prefix/checkpoint. A partial MAFFT stage
must restart that stage because MAFFT has no checkpoint here. Existing failed
groups remain flagged and are skipped while interrupted/not-started groups
continue. Inspect `run_summary.tsv` afterwards: retained failures keep the
project's final status as partial failure. The alignment QC thresholds are
unchanged. Scientific settings and prepared membership must match the saved run.
Do not start a second process against the same project.

To retry a reviewed failed group without attempt backups, use the same command
with `--resume --retry-failed --hog GROUP_ID`. Without `--resume`, the historical
`--retry-failed` option retains its archived-attempt behavior. New resume runs
record stage hashes, numbered command logs and interruption state in `run.json`.
No inference directory is renamed or backed up by `--resume`.

The stopped project retains its original orthogroups. To use the current expanded
audit for newly grouped family trees, prepare a separate family project as below;
resuming does not change an existing project's sequence membership.

Integrated phylogeny is opt-in: supply both `--phylogeny --phylogeny-grouping family`
on resume/continue. The audit completes first and is retained even if trees fail.
Without `--phylogeny`, grouping alone does not run any trees. The standalone
preparation/run commands below remain explicit optional operations.

Family mode selects the current audit's `Audited_family_assignment` worksheet: grades
A, B and C, grouped by `putative_audited_family`. Grade D and other target proteins are
excluded. Selection is data driven, with no hard-coded species or count.
The current architecture-gated rules can change earlier classifications and counts.
Grade C candidates retain provisional family assignments for grouping. Older saved
audits using `Resolved_candidates` or `Collapsed_candidates` remain readable.

Selected targets remain included even without an orthologue or source group.
Reference proteins from linked orthogroups are retained and deduplicated by ID.
Total tree tips can therefore exceed the target count. No reference cap is
imposed. Reference membership is orthogroup evidence, not independent confirmation
of family assignment. Alignment QC remains mandatory.

Fibronectin-like candidates follow the same A/B/C rule. The broad fibronectin
hypothesis handoff and external-evidence return remain available; rejected audit
assignments do not enter these trees.

## Reuse completed orthology without running it again

In WSL, activate your workflow environment, then run from the test project root:

```bash
cd /mnt/c/Users/eaasa/Desktop/Evans_PHD_Thesis/pipeline_validation/integrin_adhesome_analysis
python 04_orthology_phylogenetic_classification/scripts/phylogeny.py prepare \
  --hogs runs/annotation/downstream/orthology/candidates/candidate_HOGs.tsv \
  --orthologues runs/annotation/downstream/orthology/candidates/direct_orthologues.tsv \
  --duplications runs/annotation/downstream/orthology/candidates/candidate_duplications.tsv \
  --proteomes runs/annotation/downstream/panel/proteomes \
  --species-list runs/annotation/downstream/panel/species.json \
  --target-species runs/annotation/downstream/panel/targets.json \
  --audit-workbook runs/final_expanded/Stage_2_output.xlsx \
  --candidate-report runs/annotation/downstream/orthology/candidates/candidate_report.tsv \
  --output runs/family_phylogeny_audited
```

If duplication evidence is unavailable, replace `--duplications PATH` with
`--allow-missing-duplications`. Use a new output folder. To prepare one family,
append `--family talin`; repeat for several families. Use audit `putative_audited_family`
names. Existing prepared trees are not rewritten.

Preparation does not launch alignment or tree inference. First inspect
`runs/family_phylogeny_audited/groups.tsv`: `candidates` and `reference_sequences`
report counts separately; `sequences` gives total tips per family.
`family_candidate_coverage.tsv` records original/prepared IDs, grades and targets
included without a source group. `family_orthogroups.tsv` records contributing
groups; its `member_count` is the original size before target filtering.
Per-family `tips.tsv` records actual membership. The manifest records selection
policy, grade counts and input hashes.

Then run:

```bash
python 04_orthology_phylogenetic_classification/scripts/phylogeny.py run \
  --project runs/family_phylogeny_audited --threads 2 --iqtree-threads 2
```

For interruption recovery, repeat this run command with `--resume`. To retry
reviewed failed families in place, use `--resume --retry-failed`. Completed
families are preserved. Fewer than four tips produces an intentional skip.

## Integrated execution

Add `--phylogeny-grouping family` to a new `resume` or `continue` run. The workflow
performs the independent audit after orthology, then passes that audit and the
panel target-species list to phylogeny preparation automatically.
`continue` starts a new downstream run; the commands above reuse existing results.
Recovery of an already prepared broad-family project retains its original
selection. Prepare a new project to apply the revised selection.

The previous `--family-map TSV --family-column COLUMN` mode remains available
explicitly for broad hypothesis exploration and includes entire source groups.
Do not use it for the restricted audited selection described here.

Current final outputs include compact `putative_audited_family`, `protein_name`
and evidence confidence, plus a matching TSV. See [naming and output policy](COMPUTATIONAL_AUDIT.md).
Hard/background-negative exclusions are absent from the audited assignment table
and therefore cannot enter its family phylogenies.
