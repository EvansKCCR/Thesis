"""Prepare one complete reference proteome using explicit protein-to-gene mapping."""
import argparse
import json
import re
from collections import defaultdict
from pathlib import Path
from pipeline import records, digest, read_tsv, write_tsv, AA


def prepare(species, source, mapping, out):
    if not re.fullmatch(r'[A-Za-z][A-Za-z0-9_]*', species) or '__' in species:
        raise ValueError('Use a species code starting with a letter and no double underscores')
    fields, rows = read_tsv(mapping)
    if not {'protein_id','gene_id'} <= set(fields or []):
        raise ValueError('Mapping TSV requires protein_id and gene_id columns')
    genes = {}
    for row in rows:
        pid, gene = row['protein_id'].strip(), row['gene_id'].strip()
        if not pid or pid in genes or not re.fullmatch(r'[A-Za-z0-9_.-]+', gene):
            raise ValueError(f'Duplicate/invalid mapping: {pid}')
        genes[pid] = gene
    groups, seen = defaultdict(list), set()
    for header, seq in records(source):
        pid = header.split()[0]
        if pid in seen or pid not in genes:
            raise ValueError(f'Duplicate protein or missing gene mapping: {pid}')
        seen.add(pid)
        if seq.endswith('*'):
            seq = seq[:-1]
        if not seq or set(seq)-AA:
            raise ValueError(f'Invalid sequence: {pid}')
        groups[genes[pid]].append((pid,seq))
    if not groups:
        raise ValueError('Empty reference proteome')
    out.mkdir(parents=True, exist_ok=False)
    target = out/f'{species}.faa'
    audit = []
    with target.open('w', encoding='ascii') as f:
        for gene, variants in sorted(groups.items()):
            chosen = sorted(variants, key=lambda x:(-len(x[1]),x[0]))[0]
            prepared = species+'__'+gene
            f.write(f'>{prepared}\n{chosen[1]}\n')
            for pid,seq in variants:
                audit.append(dict(protein_id=pid,gene_id=gene,prepared_id=prepared,
                                  representative_id=chosen[0],selected=int(pid==chosen[0]),length=len(seq)))
    write_tsv(out/'identifier_map.tsv', list(audit[0]), audit)
    (out/'preparation.json').write_text(json.dumps(dict(pipeline_version='1.0',species=species, proteins=len(seen),genes=len(groups),
        policy='longest isoform per explicit gene; ties by protein ID; one terminal stop removed',
        input_fasta=str(source.resolve()),input_sha256=digest(source),
        mapping=str(mapping.resolve()),mapping_sha256=digest(mapping),
        output_sha256=digest(target),unused_mapping_rows=len(set(genes)-seen)),indent=2)+'\n')
    return target


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--version',action='version',version='adhesome-gene-preparation 1.0')
    p.add_argument('--species',required=True)
    p.add_argument('--fasta',type=Path,required=True)
    p.add_argument('--map',dest='mapping',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True,help='New output directory')
    a=p.parse_args()
    try:
        print(prepare(a.species,a.fasta,a.mapping,a.out))
    except (ValueError,OSError,KeyError) as exc:
        p.exit(1,f'ERROR: {exc}\n')


if __name__=='__main__':
    main()
