"""Prepare one complete reference proteome using explicit protein-to-gene mapping."""
import argparse
import hashlib
import json
import re
from collections import defaultdict
from pathlib import Path
from pipeline import records, digest, read_tsv, write_tsv, AA


BENCHMARK_SEQUENCE_POLICY = 'remove one terminal stop; mask remaining stop markers as X; retain protein IDs and positions'


def normalize_sequence(seq, protein_id, *, internal_stop_policy='reject'):
    """Return the inference sequence; stop masking is explicit and benchmark-scoped."""
    if internal_stop_policy not in ('reject', 'mask'):
        raise ValueError('Unknown internal stop policy: '+internal_stop_policy)
    seq = seq.upper().removesuffix('*')
    if internal_stop_policy == 'mask':
        seq = seq.replace('*', 'X')
    invalid = set(seq)-AA
    if not seq or invalid:
        detail = 'empty after terminal-stop removal' if not seq else 'unsupported characters '+repr(''.join(sorted(invalid)))
        raise ValueError(f'Invalid sequence: {protein_id} ({detail})')
    return seq


def prepare(species, source, mapping, out, *, internal_stop_policy='reject'):
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
    groups, seen, qc = defaultdict(list), set(), []
    for header, seq in records(source):
        pid = header.split()[0]
        if pid in seen or pid not in genes:
            raise ValueError(f'Duplicate protein or missing gene mapping: {pid}')
        seen.add(pid)
        original = seq
        seq = normalize_sequence(seq, pid, internal_stop_policy=internal_stop_policy)
        if original != seq:
            stops = [str(i+1) for i,c in enumerate(original.removesuffix('*')) if c=='*']
            qc.append(dict(protein_id=pid,original_length=len(original),prepared_length=len(seq),
                terminal_stop_removed=int(original.endswith('*')),internal_stops_masked=len(stops),
                masked_positions_1based=';'.join(stops),
                original_sequence_sha256=hashlib.sha256(original.encode('ascii')).hexdigest(),
                prepared_sequence_sha256=hashlib.sha256(seq.encode('ascii')).hexdigest()))
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
    write_tsv(out/'sequence_qc.tsv', ['protein_id','original_length','prepared_length','terminal_stop_removed',
        'internal_stops_masked','masked_positions_1based','original_sequence_sha256','prepared_sequence_sha256'], qc)
    (out/'preparation.json').write_text(json.dumps(dict(pipeline_version='1.0',species=species, proteins=len(seen),genes=len(groups),
        policy='longest isoform per explicit gene; ties by protein ID; one terminal stop removed',
        internal_stop_policy=internal_stop_policy,
        sequence_policy=BENCHMARK_SEQUENCE_POLICY if internal_stop_policy=='mask' else 'remove one terminal stop; reject other non-amino-acid characters',
        sequence_qc=dict(changed_proteins=len(qc),terminal_stops_removed=sum(r['terminal_stop_removed'] for r in qc),
            internal_stop_proteins=sum(r['internal_stops_masked']>0 for r in qc),
            internal_stop_sites=sum(r['internal_stops_masked'] for r in qc),report_sha256=digest(out/'sequence_qc.tsv')),
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
