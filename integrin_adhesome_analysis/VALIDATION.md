# Main validation pipeline — version 1.0

The main installation is:

```text
C:\Users\eaasa\Desktop\Evans_PHD_Thesis\pipeline_validation\integrin_adhesome_analysis
```

Expanded Stage 2 coverage is the default for `resume`, `continue` and `audit`:
190 family configurations and 136 identity rules. Both files are selected as
one pair. `--family-coverage original` is retained for explicit regression
comparisons; it is not the main validation configuration.

## Current completed validation

The current expanded result is `runs/final_expanded/`:

- `Stage_2_output.xlsx`: all target audit assignments.
- `audited_family_assignments.tsv`: compact assignments.
- `validation_evaluation/`: labelled-control metrics and diagnostics.
- `configuration/`: exact paired rules used by that audit.
- `workflow_manifest.json`: evidence and configuration checksums.

Its upstream evidence is in `runs/handoff/` and `runs/annotation/`, with completed
orthology in `runs/annotation/downstream/orthology/`. Existing earlier audit
outputs are retained as historical results. The source distribution is maintained
in `Integrin_analysis/Thesis/integrin_adhesome_analysis`.

## Refresh Stage 2 using existing evidence

From the main installation, activate `adhesome-workflow-v1` and run:

```bash
python adhesome_workflow.py audit \
  --annotation runs/annotation \
  --orthology runs/annotation/downstream/orthology \
  --out runs/final_expanded_refresh
```

Choose an unused output path. Expect `Audit family coverage: expanded (190 family
configurations; 136 identity rules).` This command reuses completed upstream
evidence. It does not rerun Stage 1, OrthoFinder or optional phylogeny.

Stage 1 recruitment and its scan library are unchanged. Expanded family coverage
does not invent missing observed domains or recruit the missing Defensin-6
control. Reference-overlapping validation is not independent held-out validation.
