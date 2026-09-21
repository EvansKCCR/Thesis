#!/usr/bin/env python3
"""Extract versioned models and separate GA models from fallback models."""
import argparse
import re
from pathlib import Path
from adhesome_core import table, hmm_records, pf, digest, save_json

def build(sheet, source, output, full_pfam=False):
    wanted = {a for row in table(sheet) for a in re.findall(r'PF\d{5}', row['Pfam accession'])}
    if not wanted:
        raise ValueError('Empty Pfam manifest')
    found, blocks, ga, fallback = {}, [], [], []
    for meta, block in hmm_records(source):
        accession = pf(meta['ACC'][0])
        if not full_pfam and accession not in wanted:
            continue
        if accession in found:
            raise ValueError(f'Duplicate model: {accession}')
        thresholds = [float(x.rstrip(';')) for x in meta['GA']] if 'GA' in meta else None
        if thresholds is not None and len(thresholds) != 2:
            raise ValueError(f'Invalid GA: {accession}')
        found[accession] = dict(name=meta['NAME'][0], accession=meta['ACC'][0],
                                length=int(meta['LENG'][0]), ga=thresholds)
        blocks.append(block)
        (ga if thresholds is not None else fallback).append(block)
    missing = wanted - found.keys()
    if missing:
        raise ValueError('Missing Pfam accessions (check release/retired models): ' + ', '.join(sorted(missing)))
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(''.join(blocks), encoding='utf-8')
    Path(str(output) + '.ga.hmm').write_text(''.join(ga), encoding='utf-8')
    Path(str(output) + '.fallback.hmm').write_text(''.join(fallback), encoding='utf-8')
    save_json(str(output) + '.metadata.json', dict(source=str(source), source_sha256=digest(source),
              manifest_sha256=digest(sheet), library_sha256=digest(output), full_pfam=full_pfam,
              split_sha256={kind:digest(str(output)+'.'+kind+'.hmm') for kind in ('ga','fallback')}, models=found))
    return found

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sheet', default=str(Path(__file__).resolve().parent / 'config/adhesome_pfam_manifest.csv'))
    parser.add_argument('--pfam-hmm', required=True)
    parser.add_argument('--out', default='databases/pfam/adhesome_v4.hmm')
    parser.add_argument('--full-pfam', action='store_true', help='Preserve all Pfam models for complete observed architectures')
    args = parser.parse_args()
    models = build(args.sheet, args.pfam_hmm, args.out, args.full_pfam)
    print(f'Extracted {len(models)} models; scan stage will run hmmpress.')

if __name__ == '__main__':
    main()
