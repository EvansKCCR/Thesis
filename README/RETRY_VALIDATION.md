# Retry the validation target without repeating completed proteomes

The retry log shows the validation entry still used a `gene=` expression. Its
headers instead use `category|accession|gene`, so the accession mapping is:

```json
"gene_regex": "^[^|]+\\|([^|]+)\\|"
```

The current configuration copy is `targets.validation.corrected.json`. The Sman,
Shae and Sjap detection outputs in `runs/handoff_retry/detection` completed
successfully; their raw/prepared FASTA and classification-input hashes were verified.
The current detector supports explicit target selection and preserves other
species and their batch status. It refuses reuse when preserved inputs, settings
or classification resources have changed. The failed target directory is archived
under `detection/archived_attempts` before replacement.

Install the corrected `targets.json` and the updated `adhesome_workflow.py`,
`01_detection/adhesome_pipeline.py` and `01_detection/scripts/run_detection.py`.
Keep the full four-target config. From the validation project root in WSL:

```bash
python 01_detection/scripts/run_detection.py --config targets.json --input-dir inputs/targets --out runs/handoff_retry/detection --cpu 4 --species validation
```

After that command completes, export the combined handoff from all four existing
detection directories:

```bash
python adhesome_workflow.py start --detection runs/handoff_retry/detection --targets-config targets.json --out runs/handoff_ready
```

The second command does not rerun HMMER or BLAST. It refuses export if any configured
target still has an incomplete detection status. Existing handoff logs are retained;
choose another new output name if `runs/handoff_ready` already exists.

For a header-only check before running detection, append `--check-inputs` to the
first command. It creates no run outputs or external analysis jobs. The regression
tests simulate external tools; biological inference is left for the local test run.
