#!/usr/bin/env python3
"""Import retained upstream domain intervals only with exact sequence verification."""

# Resolve source resources independently of the caller working directory.
import sys as _layout_sys
from pathlib import Path as _LayoutPath
_LAYOUT_BUNDLE_ROOT = next(p for p in _LayoutPath(__file__).resolve().parents
                          if (p/'config/layout_paths.json').is_file() and (p/'pipeline_paths.py').is_file())
if str(_LAYOUT_BUNDLE_ROOT) not in _layout_sys.path:
    _layout_sys.path.insert(0, str(_LAYOUT_BUNDLE_ROOT))
from pipeline_paths import activate as _layout_activate, resource_path as _layout_resource
_layout_activate()

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path
from family_motifs import read_fasta, write_tsv, digest


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input-dir', type=Path, required=True)
    p.add_argument('--domains', type=Path, nargs='+', required=True)
    p.add_argument('--source-fasta', type=Path, nargs='+', required=True,
                   help='Exact upstream search FASTA(s), whose IDs match query_name')
    p.add_argument('--out', type=Path, required=True)
    args = p.parse_args()
    if _layout_resource(args.out).exists() or _layout_resource(args.out.with_suffix('.provenance.json')).exists():
        p.error('Output already exists')
    source = {}
    for path in args.source_fasta:
        for r in read_fasta(path):
            if r['id'] in source and source[r['id']] != r['sequence']:
                p.error('Conflicting source identifier: '+r['id'])
            source[r['id']] = r['sequence']
    domains = defaultdict(set)
    for path in args.domains:
        with _layout_resource(path).open(encoding='utf-8-sig', newline='') as handle:
            reader = csv.DictReader(handle)
            required = {'query_name', 'pfam_accession', 'ali_from', 'ali_to', 'decision'}
            if not required <= set(reader.fieldnames or []):
                p.error('Unexpected upstream domain schema: '+str(path))
            for row in reader:
                if row['decision'] == 'retained':
                    domains[row['query_name']].add((row['pfam_accession'].split('.')[0], int(row['ali_from']), int(row['ali_to'])))
    rows, unverified = [], []
    files = sorted(x for x in args.input_dir.iterdir() if x.suffix.lower() in ('.fa', '.faa', '.fasta'))
    for path in files:
        for r in read_fasta(path):
            if source.get(r['id']) != r['sequence']:
                unverified.append(path.stem+':'+r['id'])
                continue
            for name, a, b in sorted(domains[r['id']]):
                if not 1 <= a <= b <= len(r['sequence']):
                    p.error('Invalid upstream interval for '+r['id'])
                rows.append(dict(family=path.stem, sequence_id=r['id'], region=name, start=a, end=b))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    write_tsv(args.out, rows, ['family', 'sequence_id', 'region', 'start', 'end'])
    args.out.with_suffix('.provenance.json').write_text(json.dumps(dict(
        inputs={str(_layout_resource(x).resolve()):digest(x) for x in files+args.source_fasta+args.domains},
        imported_intervals=len(rows), unverified_records=unverified), indent=2), encoding='utf-8')
    print(f'Imported {len(rows)} intervals; {len(unverified)} records lacked exact sequence verification.')


if __name__ == '__main__':
    main()
