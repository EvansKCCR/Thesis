# Functional motif annotation — integrated workflow v1

This is the single active motif implementation, consolidated from the previous
version1.3 source. Use the root `adhesome_workflow.py resume` command with the
current handoff and returned predictor files. Candidate family FASTAs are generated
by `start`; none are bundled here. The curated XLSX in `Motif_library/` and
`context_rules.json` remain the authoritative motif and routing resources.

Sequence matches are hypotheses. Library scores are not candidate probabilities.
High-confidence hypotheses require the configured independent context/evidence;
missing evidence remains unknown and contradictions are retained. Native topology
imports require exact IDs and full sequences before coordinates can transfer.
Generic ELM matches remain context-unverified. Optional MEME stages default to skip.

See `../WORKFLOW.md` for imports, templates and evidence contracts. Run
`python -B -m unittest discover -s test -v` from this module directory.