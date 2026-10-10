# Audited family phylogenies

Family mode selects the final audit's `Collapsed_candidates` worksheet: grades
A, B and C, grouped by `assigned_family`. Grade D and other target proteins are
excluded. Selection is data driven, with no hard-coded species or count.
The inspected test audit contains 198 A + 29 B + 45 C = **272 unique targets**;
the 45 provisional candidates are already included in that total.

Selected targets remain included even without an orthologue or source group.
Reference proteins from linked orthogroups are retained and deduplicated by ID.
Total tree tips can therefore exceed the target count. No reference cap is
imposed. Reference membership is orthogroup evidence, not independent confirmation
of family assignment. Alignment QC remains mandatory.

Fibronectin-like candidates follow the same A/B/C rule. The broad fibronectin
hypothesis handoff and external-evidence return remain available; rejected audit
assignments do not enter these trees.

## Update the local test copy

Copy these updated files to the same locations in `Test_pipeline_v1`:

- `adhesome_workflow.py`
- `04_orthology_phylogenetic_classification/scripts/phylogeny.py`

## Reuse completed orthology without running it again

In WSL, activate your workflow environment, then run from the test project root:

```bash
cd /mnt/c/Users/eaasa/Desktop/Evans_PHD_Thesis/Test_pipeline_v1
python 04_orthology_phylogenetic_classification/scripts/phylogeny.py prepare \
  --hogs runs/annotation/downstream/orthology/candidates/candidate_HOGs.tsv \
  --orthologues runs/annotation/downstream/orthology/candidates/direct_orthologues.tsv \
  --duplications runs/annotation/downstream/orthology/candidates/candidate_duplications.tsv \
  --proteomes runs/annotation/downstream/panel/proteomes \
  --species-list runs/annotation/downstream/panel/species.json \
  --target-species runs/annotation/downstream/panel/targets.json \
  --audit-workbook runs/annotation/downstream/final/audited_candidates.xlsx \
  --candidate-report runs/annotation/downstream/orthology/candidates/candidate_report.tsv \
  --output runs/family_phylogeny_audited
```

If duplication evidence is unavailable, replace `--duplications PATH` with
`--allow-missing-duplications`. Use a new output folder. To prepare one family,
append `--family talin`; repeat for several families. Use audit `assigned_family`
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

For failure recovery, repeat this run command with `--retry-failed`. Completed
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
