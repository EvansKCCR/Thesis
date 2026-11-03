# General protein homologue discovery and reciprocal BLASTP

`reciprocal_blast.py` recursively discovers reference FASTAs, matches each
reference query to its corresponding complete source proteome, and searches
any number of target proteomes. It uses Python's standard library and NCBI
BLAST+ (`blastp`, `makeblastdb`).

## Directory layout

Run from the project root in WSL/Linux:

```text
homology_retrieval/
  reciprocal_blast.py
  RBH_README.md
  environment.yml
  reference_fasta/            one reference set/family per FASTA; subfolders allowed
  tests/                     synthetic regressions; BLAST is mocked
  results/                   generated only when an analysis is launched locally
```

Supported suffixes are `.fa`, `.faa`, `.fasta`, `.fas`, `.pep`, `.aa`, optionally
compressed with `.gz`. Every input directory is searched recursively. Each
reference FASTA defines an independent named set. There is no Parvin-specific
filter or restriction on species, reference count or target count. Identical
reference sequences with distinct IDs are retained; duplicate IDs within one
FASTA are rejected. IDs in separate files are internally namespaced.

Curated queries default to this bundle's `reference_fasta/` (the legacy
`refrence_fasta/` spelling remains supported). Complete source and target
proteomes default to the **main project**, independent of the caller's cwd:

```text
../integrin_adhesome_analysis/input/reference_raw/  # reference_proteomes
../integrin_adhesome_analysis/input/target_raw/     # target_proteomes
```

These logical names are registered in the main project's
`config/input_registry.json`. Local historical proteome/output directories are
never used as implicit fallbacks. Explicit input flags still accept directories
or individual protein FASTAs. Formats are `.fa`, `.faa`, `.fasta`, `.fas`, `.fsa`,
`.pep`, `.aa`, optionally gzip-compressed. GFF/TSV sidecars are not protein inputs.

## How references are matched to source proteomes

Each query is routed to **one complete source proteome**, in this order:

1. An explicit matching row in `--source-map`, if supplied.
2. The species prefix before `|` in its ID, matched to a background filename or
   species folder. Matching is case-insensitive for source labels only.
3. A species-labelled reference folder or reference filename matching a background.
4. An exact sequence found in one unique source-proteome file.
Unidentified queries stop with an explicit `unidentified_input` assessment. A sole
unrelated background is never substituted; supply `--source-map` when necessary.

The user-confirmed alias **`candidate = homo_sapiens` applies to every
`candidate|...` reference**, across all families. Original protein IDs remain
unchanged. Its trace is
`source_assignment=user_confirmed_header_alias:candidate=homo_sapiens`.
A species explicitly named in a reference header cannot silently route to a
different species if its background is absent.

Examples: `homo_sapiens|Q9NVD7|PARVA` maps to `Homo_sapiens.faa`, and
`mus_musculus|...` maps to `Mus_musculus/primary.faa`. Arbitrary non-source
prefixes such as `sp|...` can use unique exact-sequence matching or an explicit
map. Multiple releases matching one species are ambiguous: select one release
or map the required file explicitly. References are never assigned by
approximate sequence similarity.

After source assignment, the query must match an exact background sequence.
An exact matching protein-ID token can disambiguate a sequence shared by more
than one background record; reverse-score ties still prevent RBH acceptance.

For exceptional identifiers, an optional tab-delimited source map contains:

```text
reference_file	reference_id	reference_proteome
nested/family.faa	some_unlabelled_id	Homo_sapiens.faa
```

Use literal TAB separators. Paths are
relative to their input roots. `reference_id=*` maps an entire reference file;
an exact ID rule takes priority over `*`. Background paths or discovered labels
are accepted. If `config/reference_sources.tsv` exists, it is loaded by default;
otherwise no mapping file is needed. A supplied `--source-map` overrides it.

## Local environment

Reuse your BLAST+ environment, such as `parvin-rbh`, or create the optional
general-purpose environment once:

```bash
cd /mnt/c/Users/eaasa/Desktop/Evans_PHD_Thesis/Integrin_analysis/Thesis/homology_retrieval
conda env create -f environment.yml
conda activate homology-retrieval
```

## Check inputs, run and resume

Read-only input validation requires no BLAST executables and creates no output
files or directories:

```bash
python -B reciprocal_blast.py --augment-reference-background --check-inputs
```

Input counts depend on the supplied datasets and selected reference sets.
`--augment-reference-background` explicitly adds an unchanged copy
of each missing curated sequence to that source's reverse-search database.
Every native competitor is retained. Added records are labelled in the mapping,
audit and manifest. Source FASTAs are never edited. Multiple unresolved exact
matches still fail; the augmentation flag does not force a match or resolve ties.
Omit this flag to require strict exact source-proteome representation instead.

Start locally after the checks pass:

```bash
python -B reciprocal_blast.py \
  --augment-reference-background \
  --outdir results/homology_run \
  --threads 4
```

Explicit paths, including other directories or individual FASTA inputs:

```bash
python -B reciprocal_blast.py \
  --reference-dir reference_fasta \
  --background-dir ../integrin_adhesome_analysis/input/reference_raw \
  --target-dir ../integrin_adhesome_analysis/input/target_raw \
  --augment-reference-background \
  --outdir results/homology_run \
  --threads 4
```

The older flag names `--references`, `--reference-background` and `--proteomes`
remain aliases. No JSON configuration needs to be written manually.

After a failed or interrupted process has stopped, repeat the same settings
with `--resume`:

```bash
python -B reciprocal_blast.py \
  --augment-reference-background \
  --outdir results/homology_run \
  --threads 4 \
  --resume
```

Verified completed database builds and searches are reused; an unfinished
command reruns. A completed entire run is reused without launching tools.
Changed input bytes, source mapping, code, settings, executables or completed
outputs cause resume to fail. Use a new output directory for changed analyses.
The old Parvin workflow's output folders are not adopted as general-pipeline
checkpoints. Nonempty outputs are never overwritten as a fresh run. Place
outputs outside all recursively scanned input directories.

If complete source proteomes are intentionally unavailable, explicitly select
`--reference-only` instead of `--background-dir`. Reverse search then uses each
curated reference set alone and successful matches are labelled
`provisional_reference_only_RBH`. Missing or ambiguous full-background inputs
never silently switch to this mode.

## Search and evidence policy

Each target database is built once. All reference queries share one forward
BLASTP search per target, avoiding one proteome scan per family. Qualifying top
candidates are reverse-searched separately against their assigned source
proteomes. Background databases are also built once and reused across targets.

This **per-source reverse search** differs from the old script's concatenation
of all source species into one database. It implements the requested matching
of each reference to its corresponding proteome. Results should be interpreted
relative to the named source proteome, not a single combined reference database.

Original score and quality gates remain: rank all reported hits by bitscore,
require the top hit to pass thresholds, reject exact or configured near-ties,
and require the unique reverse best hit to map to the original curated protein.
All passing forward hits are exported for broader homologue discovery, including
non-RBH hits. No family identity is inferred from RBH alone.

Defaults: E-value `1e-5`, amino-acid identity `25%`, query coverage `50%`, subject
coverage `50%`, bit-score gap `0`, maximum targets `500`, one HSP per pair and
four threads. Coverage uses the best HSP's coordinate span divided by full
protein length. `--min-bit-gap 0` rejects exact ties. Do not use a maximum-target
cap of one, which removes competitor information.

Statuses include `RBH_candidate`, `RBH_candidate_augmented_reference`,
`provisional_reference_only_RBH`, `ambiguous_forward_best`,
`ambiguous_reverse_best`, `forward_below_filters`, `reverse_below_filters`,
`reverse_matches_other_reference`, `no_forward_hit` and `no_reverse_hit`.

Internal-stop markers remain in normalized sequences and are summarized once
per input read. `sequence_qc.tsv` lists affected input records. No competing
proteins are silently discarded or residues replaced.

## Outputs

| File | Content |
| --- | --- |
| `homology_audit.tsv` | Every reference-query Ãƒâ€” target-proteome pair, including failures/no hits; source assignment and both BLAST directions |
| `all_forward_hits.tsv` | All reported forward hits, ranks and quality flags |
| `all_forward_qualifying.faa` | Unique passing target proteins across all reference sets |
| `rbh_candidates.faa` | Unique target proteins with accepted RBH/provisional support |
| `reference_mapping.tsv` | Reference set/file/ID, assigned source proteome, exact mapping or labelled addition |
| `targets_summary.tsv` | Per-target sequence counts, unique candidate counts and accepted pair counts |
| `reference_sets/<set>/homology_audit.tsv` | Audit rows for that reference set |
| `reference_sets/<set>/*.faa` | Unique all-forward and accepted-RBH FASTAs for that set |
| `sequence_qc.tsv` | Internal-stop QC records |
| `run_config.json` | Frozen input, option, code and executable hashes |
| `run_state.json` | Running, failed or complete state |
| `run_manifest.json` | Completed-run summary, mappings and output hashes |
| `commands.log` | Actual commands and reuse diagnostics |
| `work/` | Normalized FASTAs, BLAST databases/search tables and command checkpoints |

Candidate FASTA IDs are `<target_label>|<original_protein_id>` and appear once
per output, even when multiple queries support the same protein. All supporting
pairs remain separate in the audit. Do not equate accepted pair counts with the
number of unique proteins. Relative subfolders are encoded in set/target labels
(e.g. `cladeA__Species`) to distinguish identical basenames. Colliding labels
are rejected rather than overwriting another input's outputs.

## Synthetic tests

```bash
python -B -m unittest discover -s tests -p 'test_*.py'
```

Tests cover recursive discovery, arbitrary target counts, source matching,
the `candidate` alias, augmentation, ambiguity rejection, unique exports,
zero-hit results, checkpoints and resume. External BLAST calls are mocked.

Reference: [NCBI BLAST+ documentation](https://www.ncbi.nlm.nih.gov/books/NBK279690/).

## Independent main-audit / family-tree handoff

Use the [reconciliation guide](../reconciliation/docs/WORKFLOW.md) to compare a
completed `results/homology_run/rbh_candidates.faa` with the selected completed
main adhesome audit and prepare independent MAFFTâ€“trimAlâ€“IQ-TREE inputs. The
adapter checks saved hashes and original target sequences, preserves RBH status
and main classification, and honors explicit negative exclusions. It does not
launch BLAST, re-audit candidates or start phylogeny. Existing RBH result
manifests retain their original source paths; relocated source inputs are read
from the current bundle only after hash verification.

## Clean-state validation and selective loading

```bash
python -B reciprocal_blast.py --dry-run --clean-state
```

This succeeds structurally without input FASTAs and never launches BLAST. Missing
mandatory data remains an execution blocker. `--check-inputs` instead validates
real sequence/source data and fails when mandatory inputs are absent.

Select families with repeated `--reference-set <exact-label>` and restrict source
panels with repeated `--reference-proteome <exact-label>`. Headers, approved
species aliases and explicit source maps identify needed proteomes before their
contents are read. Unrelated source proteomes are not opened. Unlabelled queries
may require exact-sequence source discovery across the supplied candidate panels;
provide an explicit map/selection to keep that discovery narrowly scoped.

The source-specific reverse searches, thresholds, tie rejection and output file
names are unchanged. Selected inputs and actual loaded source proteomes retain
normal run provenance and resume checks.
