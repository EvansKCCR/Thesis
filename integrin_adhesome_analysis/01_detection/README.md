# Detection module — integrated workflow v1

Use `python adhesome_workflow.py start --targets-config targets.json --input-dir inputs/targets --out runs/handoff`
from the package root. Supply one or more configured target proteomes with
explicit gene identifiers. No analysis input or output is bundled.

The active runner is `scripts/run_domain_v5.py`. The 96-model HMM resource retains
its historical `adhesome_v4` filename; current family rules and the v5 curated
reference set drive classification. Other retained library/seed versions support
resource provenance and regression fixtures, not alternate integrated workflows.
HMMER/BLAST are external dependencies. Builders and offline source/cache resources
are retained for reproducibility. Supported and provisional hypotheses remain
conditional on evidence and are exported by the root entry point.

See `../WORKFLOW.md` for the handoff and resumed analysis. Run the local regression
suite with `python -B -m unittest discover -s tests -v`.