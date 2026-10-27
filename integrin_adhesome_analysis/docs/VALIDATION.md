# Main validation pipeline — version 1.0

Benchmark performance acceptance now requires **pairwise precision >=0.80 and
F1 >=0.70**. Valid input files alone do not imply accepted performance. See
[BENCHMARK_RECALIBRATION.md](BENCHMARK_RECALIBRATION.md) for the saved-result
diagnosis, report refresh, tree-level comparison and explicit inference tuning.

The main installation is:

```text
C:\Users\eaasa\Desktop\Evans_PHD_Thesis\pipeline_validation\integrin_adhesome_analysis
```

Expanded Stage 2 coverage is the default for `resume`, `continue` and `audit`:
222 family configurations and 197 identity rules. Both files are selected as
one pair. `--family-coverage original` is retained for explicit regression
comparisons; it is not the main validation configuration.

Adhesome-specific orthology benchmarking is available with `--scope adhesome`
on `benchmark` or `--benchmark-scope adhesome` in a resumed workflow. The original
`benchmark.py` remains the default; the explicitly selected alternative uses
`adhesome_benchmark_subset.py`. It requires
a curated RefOG/family selection TSV and retains full predicted benchmark-group
membership. See [ADHESOME_BENCHMARK.md](ADHESOME_BENCHMARK.md) for preparation,
completed-run reuse, fibronectin coverage and family-scoped audit weights.

## Current completed validation

The current expanded result is `runs/final_expanded/`:

- `Stage_2_output.xlsx`: all target audit assignments.
- `audited_family_assignments.tsv`: compact assignments.
- `validation_evaluation/`: labelled-control metrics and diagnostics.
- `configuration/`: exact paired rules used by that audit.
- `workflow_manifest.json`: evidence and configuration checksums.

Its upstream evidence is in `runs/handoff/` and `runs/annotation/`, with completed
orthology in `runs/annotation/downstream/orthology/`. Existing earlier audit
outputs are retained as historical results. The source distribution is maintained
in `Integrin_analysis/Thesis/integrin_adhesome_analysis`.

## Refresh Stage 2 using existing evidence

From the main installation, activate `adhesome-workflow-v1` and run:

```bash
python adhesome_workflow.py audit \
  --annotation runs/annotation \
  --orthology runs/annotation/downstream/orthology \
  --out runs/final_expanded_refresh
```

Choose an unused output path. Expect `Audit family coverage: expanded (222 family
configurations; 197 identity rules).` This command reuses completed upstream
evidence. It does not rerun Stage 1, OrthoFinder or optional phylogeny.

Stage 1 recruitment and its scan library are unchanged. Expanded family coverage
does not invent missing observed domains or recruit the missing Defensin-6
control. Reference-overlapping validation is not independent held-out validation.

## Resume the original stopped phylogeny

From this installation, with `adhesome-workflow-v1` active:

```bash
python orthology_phylogeny/scripts/phylogeny.py run \
  --project runs/annotation/downstream/phylogeny \
  --resume --threads 4 --iqtree-threads 2
```

This continues the existing orthogroup project with its original membership;
it does not prepare trees from the newer expanded audit. Valid completed trees,
alignments and IQ-TREE checkpoints are reused. Existing failed groups remain
flagged for review while unfinished groups continue. To retry them explicitly
in place after review, add `--retry-failed` alongside `--resume`. Resume creates
no attempt backups. See `../README/FAMILY_PHYLOGENY.md` for new family preparation.
