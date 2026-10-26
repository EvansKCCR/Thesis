LOCAL ORTHOLOGUE RETRIEVAL

Python 3.9+ and requests required.

Install:
  python -m pip install requests
Run:
  python retrieve_orthologues.py --input candidates.tsv --out orthologue_results
Stricter confidence flag:
  python retrieve_orthologues.py --input candidates.tsv --out orthologue_high_confidence --require-high-confidence
Custom species:
  python retrieve_orthologues.py --input candidates.tsv --out mouse_rat --species mus_musculus rattus_norvegicus

Defaults: mouse, rat, zebrafish, chicken, Xenopus tropicalis; vertebrates Compara.
This is a defined panel, NOT all model organisms. Invertebrates and yeast are not
covered by the default panel/database. Changing species alone does not guarantee
coverage: use a server/compara dataset containing both human and target species.
The --server and --compara options allow explicit database selection.

Input: headered TSV with accession (or input_accession). Other columns optional.
All 235 source accessions from the supplied table are included in candidates.tsv.
No gene symbols are used as rescue mappings. UniProt accession mapping is required;
ambiguous mappings and failed validations are excluded, never guessed.

Acceptance: ortholog_one2one, unique target gene/protein, target species confirmed,
protein -> transcript -> target gene linkage checked, sequence ID and length checked.
The optional high-confidence flag requires is_high_confidence=1; absent flags fail.
One-to-one means gene-level orthology, not 100% identity, functional proof, or
isoform equivalence. FASTA contains the target protein reported by Compara.
Source UniProt accession redirects are recorded, not silently hidden.

Outputs:
  all_one2one_proteins.fasta : unique target proteins
  by_species/*.fasta : one file per requested species (may be empty)
  by_candidate/*.fasta : accepted nonhuman sequences per source accession
  source_mapping.tsv : accession-to-human-gene mapping audit
  pair_audit.tsv : every candidate x requested species, including errors/exclusions
  accepted_manifest.tsv : accepted rows linked to FASTA IDs and sequence SHA256
  raw_cache/*.json : raw successful API responses, URLs, headers, retrieval times
  run_metadata.json : input hash, API versions/releases if available, status counts

NO_ORTHOLOGUE_REPORTED means none reported by the selected dataset, not proof
of biological absence. Errors remain distinct from missing orthology.
Exit code 2 means at least one unresolved mapping/query/sequence validation error;
partial accepted output is still written. Inspect pair_audit.tsv.
Existing run directories are rejected to prevent stale FASTA contamination.
For a fresh retry choose a new --out; cached successful responses can be copied
into its raw_cache directory beforehand. --refresh ignores cache.

Validation: Python syntax and six offline control-flow checks passed. Synthetic
metadata was used only for offline tests. Live API execution could not be tested
in this chat environment. No biological sequences or orthology calls fabricated.

API documentation:
https://rest.ensembl.org/documentation/info/homology_species_gene_id
https://rest.ensembl.org/documentation/info/xref_external
https://rest.ensembl.org/documentation/info/lookup
https://rest.ensembl.org/documentation/info/sequence_id
