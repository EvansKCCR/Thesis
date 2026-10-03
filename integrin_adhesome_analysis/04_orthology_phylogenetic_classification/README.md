# Orthology and reference-panel calibration — integrated workflow v1

Use root `adhesome_workflow.py calibrate` and `orthology` commands. A configurable
JSON panel lists complete gene-level proteomes with `SpeciesCode__gene_id` headers.
Every target declared in the handoff and at least one additional reference are required. Add more
reference species by creating a new panel and rerunning full inference, then use
`--compare` to compare support by stable candidate ID. HOG IDs are run-specific.

No raw proteomes, prepared inputs, candidate lists, mapping tables, QC or inference
results are bundled. `scripts/pipeline.py prepare` remains a lower-level helper
for established five/eight-species layouts: supply its raw FASTAs, GFF/gene maps
and candidate lists first. The integrated panel interface consumes prepared
gene-level FASTAs and supports additional reference species without a fixed count.

Candidate HOG membership is separate from direct pairwise orthology. Reference
support and duplication ambiguity require review; calibration does not establish
function or independently validate the classifier. See `../WORKFLOW.md` and
`PHYLOGENY.md`. External tools are declared in the root environment file.