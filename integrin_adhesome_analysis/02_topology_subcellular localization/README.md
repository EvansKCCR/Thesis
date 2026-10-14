# External topology/localization return — version 1.0

Detection (`start`) produces the supported/provisional hypothesis tables,
`external_candidates.faa` and a sequence-bound evidence template. Submit those
full sequences to external predictors after detection completes. Return the
results to `adhesome_workflow.py resume` with `--topology-file`, `--context-dir`
or `--external-evidence`. Resume verifies identifiers and sequences, imports
context, annotates motifs and continues through orthology, audit and phylogeny.

See [WORKFLOW.md](../../README/WORKFLOW.md). The predictors run externally;
detection does not execute or consume their analysis.
