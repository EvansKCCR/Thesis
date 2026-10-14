# Reference-gated audit — version 1.0

The current implementation is `integrated_adhesome_pipeline.py`, consolidated from
the supplied updated classifier. Detection and audit share the single canonical
`../01_detection/config/family_rules.json`. The obsolete implementation and its
separate versioned-rule default were removed.

Run from the workflow root:

```bash
python adhesome_workflow.py audit --annotation runs/annotation --orthology runs/annotation/downstream/orthology --out runs/final_current
```

The wrapper joins sequence-verified context, current motif matches and orthology,
then invokes this audit. It does not invent reference/subfamily labels from HOG
membership. Independently computed identity evidence may be supplied in the
handoff's external-evidence template and returned through `resume`.

Outputs preserve source worksheets and add `Resolved_candidates`,
`All_hypotheses_audited`, `Manual_review_queue` and `Pipeline_summary`.
`classification_state` distinguishes `reference_supported`, `provisional_fragment`,
`broad_only`, `ambiguous` and `unresolved`. Broad-only grade C candidates carry
`broad_homology_group`; their `assigned_family` remains unassigned. Family phylogeny
preserves these as explicitly broad provisional groups. Grades/counts are data driven.

Old motif IDs absent from the current library are excluded from joined scoring,
with `MotifPrediction_ignored_superseded_count` recording the count. Current regex
matches remain sequence candidates; definition evidence does not establish function.
Recovery of an unfinished run refreshes changed audits and dependent family-tree
preparation while reusing completed orthology. Completed runs remain separate;
create a new audit output to apply new resources.

`adhesome_reference_benchmark_panel.tsv` is a prospective reference/negative-control
panel. The evaluator requires independent sequence/domain provenance and verified
exclusion of test references; it refuses an unpopulated benchmark. Expected gold
labels must not be copied into prediction evidence. No empirical accuracy metrics
are supplied by this integration. Run fixture tests from this directory with
`python -B -m unittest discover -p 'test_*.py'`.
