# Input-folder and validation-header fix

The reported command failed because it requested `inputs/targets` while the folder
was named `input/targets`. The folder was subsequently renamed to `inputs/targets`.
The original attempt reached no sequence preparation or scanning.

A second issue was present in `targets.json`: the validation panel uses headers
such as `candidate|P12814|ACTN1`, `hard_negative|...` and
`background_negative_provisional|...`, rather than `gene=` tags. Its mapping must
extract the accession between the first two pipes:

```json
"gene_regex": "^[^|]+\\|([^|]+)\\|"
```

`targets.validation.corrected.json` is a corrected copy of the supplied config.
It preserves the three proteome mappings, uses each actual filename, and assigns
this accession mapping to the validation panel. The input check verified all four
files: 10,935 Sman records, 14,700 Shae records, 16,936 Sjap records, and 494 validation
records. It did not run biological inference or write run outputs.

For the existing validation installation, copy that corrected JSON to the project
root as `targets.json`. From that root in WSL, start a new output directory:

```bash
python adhesome_workflow.py start --targets-config targets.json --input-dir inputs/targets --out runs/handoff_retry --threads 4
```

The source workflow additionally includes `--check-inputs`, which validates all
target files and gene mappings before creating output directories or invoking
external analysis tools. To use that new option in the validation installation,
copy these updated source files to their corresponding locations:

- `adhesome_workflow.py`
- `01_detection/adhesome_pipeline.py`
- `01_detection/scripts/run_detection.py`

Then append `--check-inputs` to the start command for an input-only check. Existing
`runs/handoff` logs are preserved; no downstream recovery is needed for the original
attempt because detection did not begin. The test installation was not modified
automatically.
