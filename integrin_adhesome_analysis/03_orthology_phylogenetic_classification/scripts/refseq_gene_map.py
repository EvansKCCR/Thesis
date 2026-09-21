#!/usr/bin/env python3
"""Extract exact RefSeq protein_id -> NCBI GeneID mappings from CDS features."""
import argparse
import gzip
import json
from collections import defaultdict
from pathlib import Path
from urllib.parse import unquote
from pipeline import digest, records, write_tsv


def parse_gff(path):
    with path.open('rb') as f:
        compressed = f.read(2) == b'\x1f\x8b'
    mappings = defaultdict(set)
    metadata = []
    with (gzip.open if compressed else open)(path, 'rt', encoding='utf-8-sig') as f:
        for number, line in enumerate(f, 1):
            if line.startswith('##FASTA'):
                break
            if line.startswith('#!'):
                metadata.append(line.strip())
            if line.startswith('#') or not line.strip():
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) != 9:
                raise ValueError(f'Invalid GFF line {number}')
            if fields[2] != 'CDS':
                continue
            attributes = dict(item.split('=', 1) for item in fields[8].split(';') if '=' in item)
            protein = unquote(attributes.get('protein_id', ''))
            if not protein:
                continue
            genes = {unquote(item).split(':', 1)[1] for item in attributes.get('Dbxref', '').split(',')
                     if unquote(item).startswith('GeneID:')}
            if not genes:
                raise ValueError(f'CDS protein {protein} has no explicit GeneID on line {number}')
            mappings[protein].update(genes)
    return mappings, metadata


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gff', type=Path, required=True)
    parser.add_argument('--fasta', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    mappings, metadata = parse_gff(args.gff)
    proteins = {h.split()[0] for h, _ in records(args.fasta)}
    problems = [dict(protein_id=p, status='missing' if not mappings.get(p) else 'conflicting',
                     gene_ids=','.join(sorted(mappings.get(p, []))))
                for p in sorted(proteins) if len(mappings.get(p, [])) != 1]
    write_tsv(args.output.with_suffix('.issues.tsv'), ['protein_id','status','gene_ids'], problems)
    if problems:
        raise ValueError(f'{len(problems)} missing/conflicting protein mappings; inspect issues TSV')
    rows = [dict(protein_id=p, gene_id='GeneID_'+next(iter(mappings[p]))) for p in sorted(proteins)]
    write_tsv(args.output, ['protein_id','gene_id'], rows)
    provenance = dict(gff=str(args.gff.resolve()), gff_sha256=digest(args.gff),
                      fasta=str(args.fasta.resolve()), fasta_sha256=digest(args.fasta),
                      annotation_metadata=metadata, mapped_proteins=len(rows),
                      genes=len({r['gene_id'] for r in rows}),
                      gff_proteins_absent_from_fasta=len(set(mappings)-proteins),
                      mapping_sha256=digest(args.output))
    args.output.with_suffix('.provenance.json').write_text(json.dumps(provenance, indent=2)+'\n')
    print(json.dumps(provenance, indent=2))


if __name__ == '__main__':
    main()
