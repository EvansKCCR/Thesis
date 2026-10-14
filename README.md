# Integrin–adhesome integrated workflow v1.0

A modular source distribution integrating detection, external topology/
localization, motif annotation, whole-proteome orthology and conservative family
assignment. Start with [WORKFLOW.md](README/WORKFLOW.md) for setup and all commands.

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
| integrin_adhesome_analysis/adhesome_workflow.py | Single version 1.0 command-line entry point |
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

For step-by-step reference preparation and orthology/phylogeny commands, start with [ORTHOLOGY_PHYLOGENY_GUIDE.md](ORTHOLOGY_PHYLOGENY_GUIDE.md).

Detection uses the current 96-entry family catalogue, 176-model Pfam library and 230-component crosswalk. The reference-gated audit shares the same canonical rules. Resume scans the current 53-rule adhesion motif set, with the 347-rule master set available through `--motif-scope master`. `start` stops at the hypothesis/external-evidence handoff; topology, localization and phylogeny are part of the resumed workflow. See the [resource update guide](README/CURRENT_RESOURCES.md).

Default resume is continuous: import topology/localization -> motif annotation -> automatic reference preparation from input/reference_raw -> generated panel -> orthology -> independent audit -> phylogeny. See [automated-run instructions](README/ORTHOLOGY_PHYLOGENY_GUIDE.md). Use continue with an existing annotation for reference-panel expansion without repeating motif annotation.

To replicate the analysis from the integrin adhesome reconstruction [paper] (link to be inserted later), download the download the associated files [here](https://drive.google.com/drive/folders/1Y4YH8uxdom6Khsgdj0iUSA_Hts5eDgNU?usp=drive_link).
