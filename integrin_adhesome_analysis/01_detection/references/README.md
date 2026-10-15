# Reference resources — workflow 1.0

`adhesome_v1.0/` is the active paired reference panel. It includes all 230
components from `config/adhesome_component_crosswalk.csv`, existing family/ECM
references (including fibronectin-like proteins), spiralian anchors and
non-adhesome competitors. There are 354 forward and 425 reciprocal reference
records. Multiple explicit hierarchical labels can share one sequence; reciprocal
comparison groups tied records by accession.

`component_catalogue.csv` contains one reviewed, gene/taxon-verified human record
per component. `component_coverage.csv` records every original pin, selected
accession, rejected alternatives and explicit family labels. `component_expansion.json`
records provenance, and `reference_lock.json` hashes the published files.
All component references are available for component hypotheses; 48 current
family labels have direct reference coverage. The 96 family rules are independent
of the component inventory, so remaining subfamilies may still need calibration.

The expansion builder can read this panel, replace its component entries and
publish a new output without accumulating duplicates:

```bash
python 01_detection/expand_component_references.py --base 01_detection/references/adhesome_v1.0 --crosswalk 01_detection/config/adhesome_component_crosswalk.csv --out reference_panel_updated --repair-identity
```

Run that command from the workflow root. It retrieves public UniProt records and
requires unique reviewed human gene/taxon matches. It never infers family labels
from shared domains. `build_reference_sets.py` is the legacy/base-seed builder;
it cannot overwrite an already expanded component panel.

`HsaAdhesomeReferencePanel.faa` and its generated `.reference.json` sidecar form
an optional targeted 242-accession human orthology panel. They are not a complete
proteome. See [the crosswalk panel guide](../../../README/COMPONENT_REFERENCE_PANEL.md)
for BLAST-only refresh, output interpretation and calibration instructions.

Historical seed/panel names retain provenance; the active pipeline version is 1.0.
