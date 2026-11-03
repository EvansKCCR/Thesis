# Independent family phylogeny preparation

The historical preparation covered 15 unified family analyses from the **51
orthogroups explicitly requested**: the original nine families plus six additions
prepared on 29 October 2026. It is separate from the adhesome workflow:
there are no hooks into its audit, benchmark, calibration or saved phylogeny.
HSP70-like proteins have their own family analysis under `associated_protein`.

This guide describes the prepared library and reusable methods. Live progress is recorded in results, not in documentation.

The source is the current stored target OrthoFinder result:

```text
C:/Users/eaasa/Desktop/Evans_PHD_Thesis/pipeline_validation/integrin_adhesome_analysis/runs/annotation/downstream/orthology/inference/Results_Oct17
```

The saved orthology workflow manifest identifies this result directory. Each
selected source FASTA was checked against its recorded workflow hash and against
the membership/species columns in `Orthogroups/Orthogroups.tsv`. Orthogroup numbers
are specific to this result; identically numbered groups from other runs must
not be substituted without checking membership.


## Prepared family library: original nine families

| Family / folder | Orthogroups | Records | Distinct sequences | Validation-panel records |
| --- | --- | ---: | ---: | ---: |
| Fes/Fps-like / `fes_fps_like` | OG0000625 | 59 | 56 | 0 |
| Rho GTPase-like / `rho_gtpase_like` | OG0000876, OG0008033, OG0009113 | 73 | 53 | 1 |
| Src-like / `src_like` | OG0000122, OG0001280 | 170 | 159 | 3 |
| Spectrin-like / `spectrin_like` | OG0001980 | 35 | 34 | 1 |
| Tubulin-like / `tubulin_like` | OG0000130, OG0009122, OG0012564 | 138 | 106 | 3 |
| ADF/cofilin / `adf_cofilin` | OG0000416, OG0002163, OG0002221 | 138 | 128 | 1 |
| Profilin-like / `profilin_like` | OG0013725 | 5 | 4 | 0 |
| Actin / `actin` | OG0008484, OG0000012, OG0006209, OG0006780, OG0006649, OG0007724, OG0000437, OG0002740, OG0014370 | 556 | 517 | 0 |
| HSP70-like / `hsp70_like` | OG0004053, OG0003628, OG0000082, OG0003875 | 232 | 200 | 1 |
| **Original subtotal** | **27** | **1,406** | **1,257** | **10** |

Six additional families are prepared in the same library:

| Family / folder | Orthogroups | Prepared records | Records in the current run after excluding `validation` |
| --- | --- | ---: | ---: |
| Talin-like / `talin_like` | OG0001037 | 48 | 47 |
| SHP PTPN6/11-like / `shp_ptpn6_11_like` | OG0001434 | 41 | 39 |
| INSR-like / `insr_like` | OG0000730 | 56 | 55 |
| Laminin-like / `laminin_like` | OG0000671, OG0000925, OG0001615, OG0002753 | 177 | 177 |
| Collagen-like / `collagen_like` | OG0000219, OG0000247, OG0000250, OG0009857 | 285 | 285 |
| Tetraspanin-like / `tetraspanin_like` | OG0000840, OG0000290, OG0009188, OG0000180, OG0001078, OG0000303, OG0001445, OG0010672, OG0000405, OG0006047, OG0000534, OG0006646, OG0011637 | 607 | 606 |
| **Added subtotal** | **24** | **1,214** | **1,209** |

The expanded library contains **2,620 records**, including 15 explicitly marked
mixed-control records. The `proteome_panels_only` selection contains 2,605
records across the 15 families. Preparation verified that all original family
FASTA/metadata files and 76 immutable run files remained byte-for-byte unchanged.
New inputs retain complete orthogroup membership; fragments and identical
sequences with different identifiers have not been removed.

For each `library/<family>/`:

- `input.fasta`: all records from the specified orthogroups combined into one
  family input, retaining original identifiers, full headers and amino-acid
  sequences. No sampling, deduplication, fragment removal or domain extraction.
- `tips.tsv`: sequence ID, source orthogroup, source panel, record role, length,
  original header and sequence SHA-256. Use this to label or colour tree tips.
- `panel_counts.tsv`: record counts by source panel.

Shared library files:

- `family_summary.tsv`: family/group counts, sequence counts and length ranges.
- `source_orthogroups.tsv`: all 51 source FASTA paths, counts and hashes.
- `library_manifest.json`: preparation provenance and checksums of exported files.
- `../families.json`: the requested family-to-orthogroup selection.

`validation` is a mixed reference/control panel, **not a biological species**.
Its 15 records are retained and explicitly marked `mixed_control_panel`. They
must not be counted as an extra lineage or as schistosome gene copies. An
explicit exclusion option is provided below. Identical sequences with different
IDs are preserved; `unique_sequences` reports their distinct sequence count.

Family labels describe the requested groupings. Combining orthogroups does not
independently establish that every member has the same specific identity. Full
lengths range from 29 to 4,337 amino acids across these sets; inspect the combined
alignments and retained homologous regions before interpreting unified trees.
No conserved-core coordinates or outgroups were guessed. Default outputs are
unrooted; rooting requires a justified outgroup or a separately documented choice.
The bundle is a separate phylogenetic analysis using existing orthology data,
not an independent validation dataset for the adhesome classifier.


## Output formats and branch evidence

For each processed family, the new output directory contains:

| File | Contents |
| --- | --- |
| `branch_evidence.tsv` | Internal-branch support labels, tip membership, source groups, orthology counts and duplication-event comparisons |
| `orthology_evidence.tsv` | Unique unordered protein pairs from native pair tables, with source file/row links and a flag indicating whether both tips are retained in this tree |
| `duplication_evidence.tsv` | Normalized original event child sets, original OrthoFinder support/type/node labels, source row and missing/excluded tips |
| `branch_evidence.provenance.json` | Input and output hashes, family selection, exclusions, rooting declaration and interpretation policies |

At the annotation-directory root, `evidence_sources.json` records checksums of
all evidence tables read and `annotation_summary.tsv` reports completed/skipped
families and evidence counts. Original trees, run checkpoints and alignments are
not changed. No backups are generated.

Evidence comes from **native full**
`Orthologues/Orthologues_*/*.tsv` and
`Gene_Duplication_Events/Duplications.tsv`. One direction is read per species
pair, and repeated protein pairs are counted once. The candidate-only
`candidates/candidate_duplications.tsv` can be empty even when the full
duplication table contains events, so it is not substituted. Edges are never
invented between different orthogroups simply because their groups were merged
into a family tree. Native evidence files are checked against the original saved
orthology manifest; the source membership table must match the library hash.

Identifiers are matched by exact original IDs, exact OrthoFinder species-added
prefix aliases, or explicit `WorkingDirectory/SequenceIDs.txt` aliases. An
unresolved evidence ID is an error. Tree tips must exactly match the selected
family IDs after the recorded panel exclusions; a changed original completed
tree is rejected when its saved checksum is available.

The first eight branch columns retain the previous table's names, followed by
explicit split/provenance fields. Counts now use full native orthology evidence,
so they can exceed the earlier candidate-only counts. `species` excludes the
mixed `validation` collection; those tips are separately identified.

By default, each nontrivial **unrooted bipartition** is represented once.
`tip_ids` is its smaller side (lexicographic tie-break), and
`complement_tip_ids` is the other side. Thus the displayed side does not imply
ancestry. Support labels are preserved without interpreting their numeric scale.
`direct_orthology_edges_within` counts unique edges wholly inside the displayed
side; `orthology_edges_across_split` and `orthology_edges_in_complement` provide
the other counts.

`duplication_rows_within` lists complete original events wholly within the
displayed side. `duplication_split_compatible_rows` lists events whose two
complete child sets fall on opposite sides. This is **topological compatibility**,
not proof of a duplication on that branch. Event identifiers such as `D00000042`
link to line 42 of the original full duplication TSV and to the sidecar table.
Events with excluded/missing genes remain in the sidecar and in the relevant
`incomplete_duplication_rows`; they are not pruned into artificial exact matches.
Zero matches do not establish that duplication did not occur.

`exact_duplication_rows` remains empty for unrooted trees. Only if you have
already supplied biologically rooted trees may you pass `--rooted`; this flag
does not root them. It enables comparison of complete original event child sets
to both immediate children of an existing bifurcating internal node. Original
OrthoFinder species-tree node labels refer to its original reconciliation, not
to a newly reconciled family tree. Root-node events remain in the event sidecar;
the tree root has no incoming branch and is not emitted as a branch row.


## Filename equivalents and analysis boundaries

| Previous integrated name | Standalone counterpart |
| --- | --- |
| `gene_tree.iqtree` | `tree.iqtree` |
| `gene_tree.treefile` | `tree.treefile` |
| `trimmed.faa` | `trimmed.fasta`, the same protein FASTA format |
| `run.json` | Information is split among `run_config.json`, stage command JSON and completion receipts; schemas are not interchangeable |
| `branch_evidence.tsv` | Generated only by the explicit post-processing annotator |

Do not rename tracked outputs when checkpoint reuse is needed. `tree.treefile` can
exist before analysis finishes; `tree.complete.json` is the runner's completion
receipt. The root `run_status.json` summarizes the requested invocation, not an
automatically inferred success state for every library family.

`MEGA_analysis/` contains user-managed MEGA inputs/projects (`.mas`, `.mdsx`). It is
retained independently of the Python runner and is not a disposable cache.

Use [WORKFLOW.md](WORKFLOW.md) for commands and [MAINTENANCE.md](MAINTENANCE.md)
for recovery, quality failures, checks and directory upkeep.

Fresh source installations do not contain the historical prepared library or
results. Configuration preserves its requested OG identifiers, which are specific
to that source inference. Use `run_phylogeny.py --dry-run --clean-state` for
structural checks; real execution requires a new verified library.
