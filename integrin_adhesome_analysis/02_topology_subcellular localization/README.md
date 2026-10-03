# External evidence handoff

No analysis results are bundled. `adhesome_workflow.py start` writes supported
and provisional hypothesis tables, `external_candidates.faa`, and a hash-bound
external-evidence template. Submit that FASTA to external predictors without
changing identifiers or sequences. Keep returned results outside this source
directory, then use `resume` to import them. See `../WORKFLOW.md` for accepted
filenames, partial-coverage behavior, and native versus summary evidence.
