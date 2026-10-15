# Crosswalk reference expansion — workflow 1.0

The active panel in `01_detection/references/adhesome_v1.0` covers all **230**
components in `01_detection/config/adhesome_component_crosswalk.csv` with one
reviewed human sequence per component. Original family references, spiralian
anchors, fibronectin-like references and non-adhesome competitors are retained.

| Resource | Records |
| --- | ---: |
| Crosswalk component references | 230 |
| Additional explicitly family-labelled component entries | 53 |
| Combined forward BLAST panel | 354 |
| Combined reciprocal competition panel | 425 |
| Configured family labels represented | 48 of 96 |
| Optional human orthology reference FASTA | 242 unique accessions |

The 354/425 counts are reference records: a sequence can have both a component
entry and explicitly configured hierarchical family labels. Reciprocal comparison
groups tied records by accession, so duplicate labels on the same sequence do not
create false family conflicts. Other components and non-adhesome entries remain
competitors. Component reciprocal agreement is panel-relative, not a one-to-one
orthology assignment.

## Identity and coverage

Sequences and metadata were retrieved from the [official UniProt API](https://www.uniprot.org/help/api).
Every selected component passed human taxon, gene-name/synonym or GeneID, reviewed
record and sequence validation. `component_coverage.csv` maps each original
crosswalk row to the selected accession. `component_expansion.json` retains source
URLs, retrieval/cache records, sequence/entry versions, hashes and rejected pins.
`reference_lock.json` hashes all published panel files. The input crosswalk remains
unchanged; its original assertions are preserved.

Two pins required unique human gene/taxon resolution: ITGA7 `Q4LE35` → `Q13683`
and RAC1 `P15154` → `P63000`. In multi-accession rows, DOCK1 selects `Q14185`
and rejects the ELMO1 accession `Q92556`; PRKCA selects `P17252` and rejects
accessions for other PKC proteins. Rejected alternatives remain in provenance.

`component_catalogue.csv` / `component_references.faa` contain the 230 component
entries. `component_family_catalogue.csv` / `component_family_references.faa`
contain extra family-labelled entries only where the current rules explicitly
name the verified reference gene, or the base panel already assigns that accession
to the family. Shared Pfam domains or functional categories never choose a family.
The expanded files used automatically by detection are `search_catalogue.csv`,
`search_references.faa`, `metazoan_catalogue.csv` and `metazoan_references.faa`.

Full component coverage does not cover every configured family/subfamily: many
rules describe paralogs absent from this crosswalk. These retain explicit reference
gaps in `reference_coverage.csv`. Added human references are sequence identity
references, not newly verified experimental adhesome evidence.

## Apply the expanded panel to existing domain scans

From the validation project root in WSL:

```bash
python 01_detection/scripts/refresh_detection_evidence.py --detection runs/handoff_retry/detection --reference-set 01_detection/references/adhesome_v1.0 --out runs/detection_expanded_panel --cpu 4
python adhesome_workflow.py start --detection runs/detection_expanded_panel --targets-config targets.json --out runs/handoff_expanded_panel
```

Use new output names. Supplying `--reference-set` explicitly runs forward/reverse
BLAST for every target against the expanded panel and reclassifies the results.
It reuses cached prepared FASTAs and HMM/domain scans; it does not rerun HMMER.
Without that flag, refresh only reclassifies existing homology evidence and cannot
create hits to new references. Original runs and reference snapshots remain
unchanged. `blastp` must be available for the first command.

The classifier annotates component matches in `component_hypotheses.csv`, including
`reference_similarity` and `reference_reciprocal`. Candidate and handoff tables
also include `crosswalk_reference_components` and `crosswalk_reciprocal_components`.
If a positive component hit has no retained family hypothesis, it is preserved as
a provisional `component_reference_candidate`, with family identity unassigned.
This follow-up bucket supports external topology/localization return and resume.
It is excluded from pooled family motif discovery; sequence motif scans remain
available, and the audit independently evaluates the canonical family rules.

## Optional orthology calibration

`HsaAdhesomeReferencePanel.faa` contains 242 unique human accessions: all 230
components plus existing human references outside the crosswalk, including
fibronectin, collagen, laminin and vitronectin references. Its explicit `gene=`
headers work with automatic reference preparation; no GFF or JSON is required.

To use this targeted panel, copy the FASTA and its generated
`HsaAdhesomeReferencePanel.reference.json` sidecar to the reference-raw directory used by resume
(default `input/reference_raw`, or the directory supplied through `--reference-raw`).
Choose a new downstream run/panel so prior orthology remains reproducible. Keep
adding appropriate complete reference proteomes for phylogenetic breadth and
orthology calibration. This FASTA is **not a complete proteome**: do not interpret
missing reference groups as biological absence. Avoid adding it as an additional
human species when a complete human proteome already represents the same taxon;
prefer the complete proteome in that case. It is not automatically inserted into
an existing OrthoFinder panel. Automatic preparation verifies the sidecar's FASTA
hash and preserves its targeted-panel scope in reference provenance.

## Rebuild after future crosswalk updates

Build a new panel without overwriting the active release:

```bash
python 01_detection/expand_component_references.py --base 01_detection/references/adhesome_v1.0 --crosswalk 01_detection/config/adhesome_component_crosswalk.csv --out reference_panel_updated --repair-identity --workers 4
```

Already expanded component entries are replaced rather than accumulated. Identity
ambiguities stop publication; sequence replacement requires a unique reviewed
human gene/taxon result and is recorded. `--offline --cache PATH` supports a
previously populated UniProt cache. `python validate_resources.py` verifies the
installed panel coverage, crosswalk/rule hashes, forward/reverse membership and
sequence hashes together with HMM/motif resources. Pipeline version remains 1.0.
