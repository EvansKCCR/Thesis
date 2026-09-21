# Family-level motif inspection

Requires Python 3.10+; the inspection and exact-word discovery stages use only the standard library. Input FASTA filenames identify the 19 configured candidate families. Original FASTA files are never edited. Unknown families and invalid sequences fail explicitly. Outputs must go into a new directory.

From this directory:

```bash
python3 family_motifs.py --out results/my_inspection
python3 family_motifs.py --out results/with_domains --regions results/verified_domain_regions.tsv
python3 family_motifs.py --out results/with_meme --external run --known-meme /path/to/protein_motifs.meme
python3 -m unittest discover -s . -p 'test_family_motifs.py' -v
```

On Windows use a Python 3 executable explicitly if `python` resolves to an older installation. Run `family_motifs.py` directly on Windows, Linux or WSL.

## Inspection and annotation

`motif_library.json` is a versioned, editable panel with patterns, sources, intended families, context rules and interpretation limits. `build_library.py` regenerates the default panel and snapshots domain count alternatives from `01_detection/config/family_rules.json`; rebuilding overwrites panel customizations.

Family names are input hypotheses, not classification results. Some families have no sufficiently specific short linear motif in this panel. These receive explicit domain-architecture-required states. Domain count alternatives support inspection of all families but do not rerun upstream order/topology/homology checks or assign family identity. Human residue numbers, LD dipeptides, generic kinase motifs and structural binding surfaces are not turned into diagnostic regexes. The zyxin EVH1 screen uses a conservative aromatic-polyproline core; generic PxxP is separately labeled.

All coordinates are **1-based inclusive**. Overlapping matches are retained. Nested motifs and collagen repeat matches must not be counted as independent evidence. C-terminal windows are provisional; supply `cytoplasmic_tail` regions to establish positional context. Domain matches without supplied coordinates remain context-unverified, not validated sites. Absent motifs do not exclude a divergent/partial protein.

Optional region TSV columns:

```text
family	sequence_id	region	start	end
integrin_beta	example_id	PF00362	100	350
integrin_beta	example_id	cytoplasmic_tail	750	800
```

Use versionless Pfam accessions and exact FASTA identifiers. A domain region must fully contain the motif. Intervals are bounds-checked. A missing domain record means evidence unavailable, not a negative domain scan. `import_domain_regions.py` imports upstream retained `domains.csv` intervals only when an upstream FASTA confirms the exact identifier **and sequence**; it writes a provenance sidecar and lists unverified records. Its `--domains` and `--source-fasta` options accept multiple files. It does not guess isoform mappings or transfer coordinates between different sequences.

## Discovery

The built-in exploratory screen searches exact 4–6 residue words by default, requiring at least two distinct supporting sequences and 50% within-family prevalence. Words with Shannon entropy below 1.5 bits are omitted. Results rank by target-minus-background prevalence; redundant contained words with identical target support are suppressed. This is a lead-generation method, not a replacement for PWM discovery, alignment or statistical testing. Singletons cannot support discovery. Pattern overlap with the known panel is annotated as sequence overlap only.

Identical sequences are collapsed within each family for prevalence/discovery. Other supplied families form the comparison set; sequences shared with the target family are excluded, and the remainder are globally deduplicated. These are candidate-family contrasts, **not measured proteome specificity**. Length/composition differences, related proteins and near-identical isoforms remain confounders. No p-values are fabricated for words selected on the same dataset. Use independent, appropriately matched sequences and homology-aware splits for validation.

Optional MEME Suite execution uses STREME `--protein`, reproducible seed 1729 and explicit other-family negative sequences. With only one family, STREME uses its shuffled control default. The default minimum is 50 distinct target sequences; `--min-streme` can lower it for exploratory use, but small-set STREME scores must not be interpreted as significance. Negative holdout size also matters. See the [STREME tutorial](https://meme-suite.org/meme/doc/streme-tutorial.html).

`--external plan` (default) writes commands without execution. `run` executes installed binaries, logs their output, and records missing tools and skipped dependencies. Any executed stage failure returns nonzero. No tools are silently installed and no remote services receive sequences. STREME, FIMO and TOMTOM must be executable on PATH in the Python environment, such as inside WSL; a PowerShell alias named `fimo` is not the MEME executable.

A real protein PWM database with `ALPHABET= ACDEFGHIKLMNPQRSTVWY` is required for known FIMO/TOMTOM annotation. Regexes are not converted into fabricated PWMs. FIMO uses a q-value threshold (default 0.05); TOMTOM uses protein-compatible Euclidean distance. Tool q-values apply within each run, not across every family/run. FIMO localization on discovery training sequences is not independent validation. TOMTOM matches suggest similarity, not biological function. Raw tool outputs retain their native reports and statistics.

## Outputs

- `REPORT.md`, `family_summary.tsv`: dataset overview and limitations.
- `sequence_qc.tsv`: lengths, ambiguities, terminal-stop normalization, identical sequences and cross-family sharing.
- `known_hits.tsv`, `inspection_states.tsv`: coordinates, flanks, sources, contextual support and explicit absent/not-applicable states.
- `known_prevalence.tsv`: whole-sequence descriptive prevalence in target versus other candidates, independent of domain availability.
- `domain_inspection.tsv`: supplied domain count support and missing alternatives; no family reassignment.
- `discovered_words.tsv`, `discovered_sites.tsv`: exploratory leads, pattern-overlap annotations and original-sequence coordinates.
- `external_stages.tsv`, `commands.json`: stage status, reasons and exact argument arrays.
- `manifest.json`, `motif_library.snapshot.json`: SHA-256 input/code hashes, parameters, runtime and executed tool versions. `inspection_completed` is distinct from `external_execution_complete`.
- Family subfolders: deduplicated protein inputs and available external results/logs. Ambiguous B/Z/J/U/O become X only in external-tool FASTAs, preserving coordinates; original sequences are used for inspection.

Only the current ELM-enabled results are retained. See `results/family_inspection_v3_elm/ELM_FINDINGS.md` for interpreted findings.

## ELM resource integration

The ELM upgrade is enabled explicitly with `--elm-dir`. It runs offline, without additional packages or MEME executables. From this directory:

```bash
python3 family_motifs.py --out results/my_elm_run --elm-dir elm_resources --regions results/verified_domain_regions.tsv --external skip
python3 -m unittest discover -s . -p 'test_*.py' -v
```

The completed dataset run is `results/family_inspection_v3_elm/REPORT.md`. Original ELM downloads are retained in `elm_resources/`; superseded runs are removed. The manifest hashes all loaded ELM files, the ELM module and the family-focus configuration. Classes are required; ancillary exports are optional and their availability is recorded. Count metadata, table schemas, regex syntax, probability range, duplicate class identifiers and reference coordinates are checked before analysis. Files called `.csv` may contain tabs; the parser reads their actual delimiter and handles comment metadata and quoted headers.

All classes are scanned against the original full protein sequence, retaining overlapping start positions and respecting anchors/lookarounds. Matches spanning ambiguous residues are rejected. At a given start, regex alternation and greediness select one match length; the scanner does not enumerate every possible variable-length endpoint. Raw matches are candidates, not validated SLiMs. ELM `Probability` is preserved as class metadata and is never reported as a candidate p-value or posterior probability.

`elm_resources/elm_family_focus.json` prioritizes classes for review, including paxillin LD helices, talin/alpha-actinin vinculin-binding sites, zyxin/vinculin EVH1 ligands, integrin-beta PTB ligands, collagen ligand precursors and Src regulatory patterns. These mappings are explicit biological review hypotheses, not learned classifiers or universal family requirements. Families with empty focus lists still receive the complete scan. ELM now provides sourced LD/VBS regexes beyond the earlier conservative panel; the matches still require structural and accessibility assessment.

`elm_hits.tsv` contains all candidates. Start with the smaller `elm_family_focus_hits.tsv`, which adds context states and review notes. `elm_inspection_states.tsv` includes absent focused classes; `elm_family_prevalence.tsv` reports every family/class combination, including zeros, using exact-sequence deduplication and exclusion of shared target sequences from controls. Comparisons remain descriptive and have the same homology, length and composition limitations as the original workflow.

Domain overlap is recorded for review and does not automatically eliminate a site: some motifs are structured or become exposed through conformational change. No domain overlap is not evidence of disorder. Optional region rows named `disordered_region` can mark supplied disorder intervals; the pipeline does not run a disorder predictor. Taxonomic applicability, localization, phosphorylation, other modifications and binding are not inferred. Tail windows remain provisional without supplied cytoplasmic-tail regions.

ELM interaction domains describe a motif's **binding partner**. They are not required domains on the protein carrying the motif. `elm_partner_domain_candidates.tsv` reports compatible supplied domains but does not create protein-interaction edges or infer coexpression, colocalization or same-species pairing. `elm_class_annotations.tsv` retains descriptions, class probabilities, partner-domain links and reference literature identifiers. `elm_reference_pathways.tsv` keeps KEGG/Reactome associations tied to their original UniProt references; it never assigns those pathways to a candidate because of a regex match.

The supplied instance subset has five entries from human SRC/P12931, not a comprehensive set of experimentally tested ELM instances. `elm_reference_instances.tsv` records regex-span recovery, annotation logic, methods and literature separately. One supplied true-positive entry explicitly lists a sequence-prediction method, so the true-positive label must not be conflated with direct experimental evidence. A candidate receives an imported reference-instance link only if its entire sequence and hit coordinates match a supplied true-positive reference. No human residue numbering is transferred to schistosome proteins.

`elm_discovery_annotation.tsv` links exploratory discovered sites to ELM matches at the same original protein positions. It distinguishes containment from partial overlap. Isolated words are never rescanned as separate proteins, which would create false terminal matches and omit necessary flanks. This annotation is positional overlap, not functional validation or independent discovery validation.

`elm_resource_audit.tsv` reports unresolved historical class identifiers rather than guessing aliases. `elm_summary.json` and the main report give counts. ELM regex downloads are not a protein PWM database and do not replace `--known-meme` for FIMO/TOMTOM. See the [ELM resource paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC8728146/) for class/instance distinctions and the need for contextual filtering.

## Current directory layout

- `family_motifs.py`, `elm_annotation.py`: current analysis pipeline.
- `motif_library.json`, `build_library.py`: conservative panel and its optional regeneration utility.
- `import_domain_regions.py`: verified upstream domain-coordinate importer.
- `elm_resources/`: all ELM downloads and family-focus configuration.
- `family_specific_candidate_library/`: original family-specific protein inputs.
- `test_family_motifs.py`, `test_elm_annotation.py`: regression checks.
- `results/verified_domain_regions.tsv` and its provenance JSON: verified domain inputs.
- `results/family_inspection_v3_elm/`: the sole retained analysis run, including findings, tables, snapshots, provenance and family-specific tool inputs.

The retained run's `manifest.json` preserves its original execution paths and code hash. `layout_provenance.json` maps relocated inputs to their current paths, verifies their hashes, and records the later default-path-only code update. The results were not recomputed during directory cleanup. Per-family positive/negative FASTAs remain because they document the external-tool inputs, even though those tools were unavailable.
