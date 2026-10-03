# Conservative family-assignment audit — integrated workflow v1

Use the root `adhesome_workflow.py audit` command to join detection, motif,
external-context and orthology evidence and create an audited workbook.
`integrated_adhesome_pipeline.py` remains callable for a correctly structured
input workbook. No result workbook or populated input workbook is bundled.

Family architecture gates are mandatory. Failed gates cannot be rescued by a
high summed score. Missing evidence stays unknown. All original hypotheses are
retained, with an additional collapsed row per species/protein. Arbitrary panel
reference species are supported via explicit direct-reference counts; gene or
family annotations are not invented from reference identifiers.