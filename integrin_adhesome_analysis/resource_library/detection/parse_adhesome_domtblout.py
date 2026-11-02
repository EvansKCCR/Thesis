#!/usr/bin/env python3
"""Filter hmmscan hits using GA/fallback thresholds and resolve same-family overlaps."""

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
import json
import math
from adhesome_core import pf, table, write_table

COLS = ['target_name','target_accession','target_length','query_name','query_accession','query_length',
        'full_evalue','full_score','full_bias','domain_number','domain_count','c_evalue','i_evalue',
        'domain_score','domain_bias','hmm_from','hmm_to','ali_from','ali_to','env_from','env_to','accuracy','description']
FIELDS = COLS + ['pfam_accession','profile_coverage','threshold_method','curated_domain_family','decision']

def parse(path, mapping, models, max_i=1e-3, max_c=1e-3, min_coverage=0.35, min_score=0):
    if not (0 < max_i < float('inf') and 0 < max_c < float('inf') and 0 <= min_coverage <= 1):
        raise ValueError('Invalid domain thresholds')
    rows = []
    with open(_layout_resource(path), encoding='utf-8') as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.lstrip().startswith('#'):
                continue
            values = line.rstrip().split(maxsplit=22)
            if len(values) < 22:
                raise ValueError(f'{path}:{number}: malformed domtblout')
            if len(values) == 22:
                values.append('')
            row = dict(zip(COLS, values))
            accession = pf(row['target_accession']) or pf(row['target_name'])
            if accession not in mapping:
                continue
            if accession not in models:
                raise ValueError(f'No threshold metadata for {accession}')
            meta = models[accession]
            if row['target_accession'] != meta['accession']:
                raise ValueError(f'Pfam version mismatch: {row["target_accession"]} vs {meta["accession"]}')
            for key in ('target_length','query_length','hmm_from','hmm_to','ali_from','ali_to','env_from','env_to'):
                row[key] = int(row[key])
            for key in ('full_score','domain_score','c_evalue','i_evalue'):
                row[key] = float(row[key])
                if not math.isfinite(row[key]):
                    raise ValueError(f'{path}:{number}: non-finite score/E-value')
            if 'length' in meta and row['target_length'] != meta['length']:
                raise ValueError('HMM model length mismatch')
            if not (1 <= row['hmm_from'] <= row['hmm_to'] <= row['target_length'] and
                    1 <= row['ali_from'] <= row['ali_to'] <= row['query_length']):
                raise ValueError(f'{path}:{number}: invalid coordinates')
            coverage = (row['hmm_to'] - row['hmm_from'] + 1) / row['target_length']
            ga = meta['ga']
            passed = (row['full_score'] >= ga[0] and row['domain_score'] >= ga[1]) if ga else (
                row['i_evalue'] <= max_i and row['c_evalue'] <= max_c and row['domain_score'] >= min_score)
            row.update(pfam_accession=accession, profile_coverage=coverage,
                       threshold_method='GA' if ga else 'Evalue', curated_domain_family=mapping[accession],
                       decision='retained' if passed and coverage >= min_coverage else 'threshold_rejected')
            rows.append(row)
    kept = {}
    for row in sorted(rows, key=lambda r: (-r['domain_score'], -r['profile_coverage'], r['ali_from'])):
        if row['decision'] != 'retained':
            continue
        group = kept.setdefault((row['query_name'], row['pfam_accession']), [])
        if any(max(row['ali_from'], other['ali_from']) <= min(row['ali_to'], other['ali_to']) for other in group):
            row['decision'] = 'same_family_overlap'
        else:
            group.append(row)
    return sorted(rows, key=lambda r: (r['query_name'], r['ali_from'], r['pfam_accession']))

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--domtblout', required=True)
    parser.add_argument('--library-sheet', required=True)
    parser.add_argument('--model-metadata', required=True)
    parser.add_argument('--out', required=True)
    parser.add_argument('--max-i-evalue', type=float, default=1e-3)
    parser.add_argument('--max-c-evalue', type=float, default=1e-3)
    parser.add_argument('--min-profile-coverage', type=float, default=0.35)
    parser.add_argument('--min-domain-score', type=float, default=0)
    args = parser.parse_args()
    mapping = {pf(r['Pfam accession']): r['Domain or family'] for r in table(args.library_sheet)}
    models = json.load(open(_layout_resource(args.model_metadata), encoding='utf-8'))['models']
    rows = parse(args.domtblout, mapping, models, args.max_i_evalue, args.max_c_evalue,
                 args.min_profile_coverage, args.min_domain_score)
    write_table(args.out, [r for r in rows if r['decision'] == 'retained'], FIELDS)
    write_table(args.out + '.audit.csv', rows, FIELDS)
    print(f'Retained {sum(r["decision"] == "retained" for r in rows)} / {len(rows)} curated hits')

if __name__ == '__main__':
    main()
