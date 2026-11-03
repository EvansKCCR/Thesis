# Independent family phylogeny: maintenance and recovery

See [MANUAL.md](MANUAL.md) for data/evidence contracts and [WORKFLOW.md](WORKFLOW.md)
for local setup and normal runs.

## Outputs and retrying locally

Runtime outputs will be created only when you execute the runner:

```text
results/<run-name>/<family>/
  selected_input.fasta
  run_config.json
  aligned.fasta
  trimmed.fasta
  tree.treefile
  tree.contree
  tree.iqtree
  tree.log
  align.command.json / trim.command.json / tree.command.json
  align.stderr.log / trim.stderr.log / tree.stderr.log
  align.complete.json / trim.complete.json / tree.complete.json
```

IQ-TREE may also emit its normal checkpoint/support files. Overall requested
family status is in `results/<run-name>/run_status.json`. A family failure is
reported, the remaining families are attempted, and the runner returns nonzero
if any requested family failed.

To resume after an interruption or failure, **repeat your exact original
command**, optionally adding `--resume`, with the same options and run name.
Activate the same environment that created the saved checkpoints. Completed
stages are checked and reused; unfinished IQ-TREE work reuses its native
checkpoint files. If model selection stopped before `tree.ckp.gz` was created,
`tree.model.gz` retains completed model evaluations; the interrupted model and
some initial setup can repeat. Unfinished MAFFT/trimAl stages lack native
checkpoint recovery and would be retried. If you change input data, tool binaries
or settings after completing a stage, choose a new run name. There is no
automatic backup generation.

For the existing `proteome_panels_only` run, created in
`adhesome-workflow-v1`, run locally in WSL:

```bash
conda activate adhesome-workflow-v1
cd /mnt/c/Users/eaasa/Desktop/Evans_PHD_Thesis/Integrin_analysis/Thesis/independent_family_phylogeny
python -B run_phylogeny.py --execute --resume \
  --run-name proteome_panels_only --exclude-panel validation \
  --alignment genafpair --trim gappyout \
  --mafft-threads 4 --iqtree-threads AUTO --seed 20261028
```

This verifies and reuses completed family stages, resumes actin at its valid
native checkpoint, and runs other selected unfinished families. Read completion
receipts for current state; do not infer live progress from documentation.
The runner accepts historical completed tree fingerprints containing `--redo`
without rewriting those receipts or repeating their analyses. Incomplete
IQ-TREE stdout/stderr logs are appended and command attempts are recorded in
the existing command JSON. Inputs, model checkpoints and completed outputs are
not renamed, deleted or copied. The original `--redo` behaviour was corrected
on 29 October 2026.

See the IQ-TREE [checkpoint and model-selection documentation](https://iqtree.github.io/doc/Command-Reference)
for native `.ckp.gz` / `.model.gz` recovery. Do not add `--redo` or `-mredo` to a
resume command: those options invalidate the corresponding recovery behaviour.


## Re-preparation from a future orthogroup library

The current library is already prepared. This command is only for an explicitly
chosen future source; check the group IDs in `families.json` first:

```bash
python -B prepare_library.py \
  --results /path/to/OrthoFinder/Results_new \
  --out /path/to/new_family_library
python -B run_phylogeny.py --library /path/to/new_family_library
```

Preparation only reads the supplied results and exports selected inputs. It
refuses to overwrite an existing library. To verify against the original source
after moving between Windows and WSL, supply its local path explicitly:

```bash
python -B check_library.py --source-results \
  /mnt/c/Users/eaasa/Desktop/Evans_PHD_Thesis/pipeline_validation/integrin_adhesome_analysis/runs/annotation/downstream/orthology/inference/Results_Oct17
```


## Alignment failures and membership safeguards

A successful trimAl exit does not guarantee successful trimming QC. The runner
requires the same identifiers after trimming; it deliberately stops if trimAl
removes an all-gap row. Known saved cases are INSR `Scura__SCUR1_62910` and four
laminin members (the Tetraodon `ENSTNIP00000004565` record and Sinter
`SINT2_97670`, `SINT2_97680`, `SINT2_97690`). Their diagnostics and alignments are
retained. This is not repaired by changing thread count or editing a completion
receipt. Review homologous coverage and fragments, then document any intended
filter/alternative trimming in a separate reproducible configuration. The current
runner does not silently drop records or promote a partial tree to complete.

Non-converged UFBoot checks can extend tree search repeatedly. A running CPU
process and a newer native checkpoint distinguish slow computation from an old
log that simply has not flushed another iteration. Logged time-left estimates
are provisional and can change after a new best tree or convergence test. A
partial `.treefile` is not a final supported result. Preserve native checkpoints
and do not add `--redo` or `-mredo` when recovering original work.

## Integrity and synthetic checks

```bash
python -B check_library.py
python -B -m unittest discover -s tests -p 'test_*.py'
```

The library check reads hashes, identifiers, tip metadata and optional source
membership; it runs no alignment or inference. Tests use synthetic fixtures;
Biopython-dependent skips must be reported. `run_phylogeny.py` and
`run_queued_families.py` plan commands unless `--execute` is supplied.
`annotate_branches.py` is a real post-processing operation when invoked normally.

## Directory upkeep

Operational documentation is centralized in three files under docs/. The root
README is an index. Preserve all prepared `library/` files, `families.json`,
`pending_families.json`, lock files, source caches, stage receipts and all analysis
results, including failed-run diagnostics, native model/tree checkpoints and
user-managed `MEGA_analysis/` projects. An empty lock file is not evidence of
an active process, nor proof that a lock can be deleted safely.

The 1 November 2026 documentation cleanup does not rename scripts, change the
15-family/51-orthogroup selection, alter QC or restart a run. Stale prose saying
that no analyses had ever run was replaced by procedures that refer to current
completion receipts. No backups or copies of research runs were created.
