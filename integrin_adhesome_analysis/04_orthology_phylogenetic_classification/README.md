# Orthology, evolutionary benchmark and optional phylogeny — version 1.0

Run from the package root containing `adhesome_workflow.py`; do not change into
this stage directory. See [automatic preparation](../../README/ORTHOLOGY_PHYLOGENY_GUIDE.md)
and [evolutionary benchmarking](../../README/EVOLUTIONARY_BENCHMARKING.md).
The original benchmark Input/RefOGs bundle belongs in `inputs/reference_raw`.
Reference preparation, panel configuration, orthology, benchmarking and audit
run automatically. Benchmark protein IDs are preserved and explicitly scoped;
ordinary references retain explicit gene preparation. The supplied
`scripts/benchmark.py --workflow` uses `benchmark_workflow.py` for verified mapping
and structured results while retaining the supplied Open_Orthobench scorer.
Trees are disabled unless `--phylogeny` is requested, and their failure cannot
prevent a completed audit. All active pipeline/version metadata remain 1.0.
