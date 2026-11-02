# Integrin-adhesome workflow 1.0

The maintained main profile is **calibrated_v1** with the current expanded family
and identity libraries. All operational documentation is in docs/:

- [Manual](docs/MANUAL.md): scientific scope, resources, evidence rules and outputs.
- [Workflow](docs/WORKFLOW.md): setup, handoffs, external return, audits and reports.
- [Maintenance](docs/MAINTENANCE.md): recovery, resource updates, diagnostics and tests.

Run commands from this directory in `adhesome-workflow-v1`. Start with the workflow
guide; do not change into a component directory for normal analysis. Inputs,
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
