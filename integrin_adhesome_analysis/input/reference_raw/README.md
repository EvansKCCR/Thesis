# Reference inputs

Place complete reference protein FASTAs here: `Code.fa`, `Code.faa` or
`Code.fasta` (optionally gzip-compressed). Supply matching `Code.gff`/`Code.gff3`
with NCBI CDS protein_id/GeneID annotations, or `Code.protein_gene.tsv` containing
protein_id and gene_id columns. Explicit gene= fields or already prepared
Code__gene_id FASTA headers can be used without a sidecar.

Run resume from the package root. References, mappings and the OrthoFinder panel
are prepared automatically. See ../../ORTHOLOGY_PHYLOGENY_GUIDE.md.
