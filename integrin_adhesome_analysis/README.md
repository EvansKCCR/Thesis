# Integrin–adhesome integrated workflow v1.0.0

A modular source distribution integrating detection, external topology/
localization, motif annotation, whole-proteome orthology and conservative family
assignment. Start with [WORKFLOW.md](WORKFLOW.md) for setup and all commands.

```text
start -> supported/provisional hypotheses + external candidate FASTA
      -> external topology/localization analysis
resume -> sequence-verified context and functional motif annotation
calibrate -> complete target/reference proteome panel
orthology -> candidate orthology + reference-panel review/comparison
         -> optionally add references and rerun with a new panel
 audit -> joined evidence and conservative audited workbook
```

The first handoff is explicitly resumable. Reference-panel calibration can add
species without overwriting earlier analyses. Missing evidence remains unknown.

| Location | Purpose |
|---|---|
| adhesome_workflow.py | Single versioned command-line entry point |
| 01_detection | Detection code, family rules, HMM and curated reference resources |
| 02_topology_subcellular localization | External handoff documentation |
| 03_functional_motif_annotation | One active motif implementation and curated workbook/rules |
| 04_orthology_phylogenetic_classification | Orthology/preparation/phylogeny scripts |
| 05_conservative_adhesome_family_assignment_audit | Final conservative audit script |
| workflow_tests | Integration regression tests |
| environment.yml / requirements.txt | External-tool and Python dependency specifications |
| panel.example.json | Editable reference-panel template |

Fibronectin-like screening candidates are explicitly preserved as provisional follow-up hypotheses across export, external-evidence return and orthology calibration. See WORKFLOW.md for interpretation and subset outputs.

Target inputs are configurable: use targets.example.json for one or more species and pass --targets-config targets.json to start. There is no three-species restriction. panel.example.json illustrates one target and one reference; expand both templates as required.
