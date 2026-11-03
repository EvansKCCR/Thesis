# Independent family phylogeny

Standalone analysis bundle for **15 families / 51 orthogroups**. It is separate
from the main adhesome workflow. All operational documentation is in docs/:

- [Manual](docs/MANUAL.md): family datasets, provenance, outputs and interpretation.
- [Workflow](docs/WORKFLOW.md): setup, local runs, queueing and branch annotation.
- [Maintenance](docs/MAINTENANCE.md): checkpoint recovery, QC failures and checks.

Run commands from this directory. Scripts plan unless their documented execution
option is supplied; branch annotation is an explicit post-processing command.
Check results/<run-name>/<family>/tree.complete.json for completion rather than
assuming an existing treefile is final. Inputs, results, native checkpoints and
MEGA_analysis projects are retained at their current paths.
