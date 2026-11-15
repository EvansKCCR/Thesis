Thesis# Integrin-adhesome workflow 1.0

The maintained main profile is **calibrated_v1** with the current expanded family
and identity libraries. All operational documentation is in docs/:

- [Manual](docs/MANUAL.md): scientific scope, resources, evidence rules and outputs.
- [Workflow](docs/WORKFLOW.md): setup, handoffs, external return, audits and reports.
- [Maintenance](docs/MAINTENANCE.md): recovery, resource updates, diagnostics and tests.

Run commands from this directory. On Linux/WSL, use `make info` to inspect the
detected host and `make prepare` to provision/check the pinned environment and
run a no-analysis preflight. On module-managed HPC, load the site's environment
and use `make prepare USE_ACTIVE_ENV=1`. Start with the workflow guide; do not
change into a component directory for normal analysis. Inputs,
returned evidence and stored runs retain their recorded locations. The separate
[family phylogeny bundle](../../independent_family_phylogeny/README.md) has its own
documentation and is not automatically part of this pipeline.

The organized bundle has three code areas: `main/` for workflow and report entry
points, `maintenance/` for maintenance and explicit diagnostics, and
`resource_library/` for component implementations, shared modules and data.
Editable settings, family/identity rules, motif rules, benchmark policy and
installation specifications are centralized in `config/`. Labelled inputs and
learned calibration resources remain in `validation/`; run manifests and resource
provenance stay alongside their data.

Use `python -B main/adhesome_workflow.py --help` or
`python -B main/run_pipeline.py --help` from this bundle. The two small root
Python files are compatibility loaders, not duplicate implementations.

An optional [machine-learning classifier](docs/WORKFLOW.md#optional-machine-learning-classifier)
can train from saved calibration outputs and evaluate on a separate validation
branch. It can classify unknown FASTA files independently, or run after Stage 1
in target mode in 250-sequence batches when `--target-ml-model` is specified;
that shortcut skips Stage 2 and downstream audit/orthology analysis and produces
ML hypotheses rather than audited assignments. The v5.3 gate learns named families
from calibration reference labels and ranks all learned families for every candidate.
It verifies the saved Pfam scans against their FASTAs and rescans sequence motifs;
features include Pfam identity, copy number, score, E-value, coverage and domain order, plus
forward/reverse BLAST similarity and coverage evidence. Exact sequence lookup,
sequence hashes, k-mer neighbors and sequence composition are not classifier
features; records without informative evidence remain unresolved. Validation reports
separate membership performance from family recovery, including abstentions and
sequence-overlap strata. Calibration CDD domain fragments add class, family, accession,
similarity, coverage, order and copy-number features. Validation fragments stay out
of the reference bank. Stage 2 identity-rule features are evaluated from sequence data only: required domain groups, specific CDD patterns, group order, copy/maximum constraints and length ranges. The v5.3 feature set also parses boolean conditions and individual Pfam/CDD signatures from the post-hoc rules. It excludes post-hoc family proposals, rescue/final decisions, orthogroup context, topology, localization and signal peptide. Calibration fits the model; validation remains evaluation-only. Existing models must be retrained when changing feature schema.
