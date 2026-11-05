# Independent analysis workflows version 1.0

This workspace contains three independently runnable workflows. Each has independent execution, configuration, environment and results. Homology
reads the main project’s registered raw inputs.

| Bundle | Entry point | Instructions |
| --- | --- | --- |
| `integrin_adhesome_analysis/` | `adhesome_workflow.py` | [Main workflow](integrin_adhesome_analysis/docs/WORKFLOW.md) | [Manual](integrin_adhesome_analysis/docs/MANUAL.md) | [Maintenance](integrin_adhesome_analysis/docs/MAINTENANCE.md) |
| `homology_retrieval/` | `reciprocal_blast.py` | [Homology/RBH manual](homology_retrieval/RBH_README.md) |
| `independent_family_phylogeny/` | `run_phylogeny.py` | [Phylogeny workflow](independent_family_phylogeny/docs/WORKFLOW.md) | [Manual](independent_family_phylogeny/docs/MANUAL.md) | [Maintenance](independent_family_phylogeny/docs/MAINTENANCE.md) |

The main adhesome workflow uses the expanded family/identity library and
`calibrated_v1` Stage 2 profile. Its maintained scripts are in `main/`,
maintenance utilities in `maintenance/`, dependencies/rules in `config/`, and
runtime resources in `resource_library/`. The root `adhesome_workflow.py` remains
a compatible entry point. Each bundle's guide specifies its working directory.

[Candidate reconciliation and independent tree preparation](reconciliation/docs/WORKFLOW.md)
provides an explicit handoff between completed RBH results, the selected target
audit and a new family-phylogeny library. It retains source decisions and negative
exclusions. It launches no scientific workflow and does not feed tree results
back automatically. Existing libraries, stored audits and tree runs remain in
their original folders.
