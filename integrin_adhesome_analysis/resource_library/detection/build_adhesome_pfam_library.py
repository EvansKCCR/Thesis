#!/usr/bin/env python3
"""Extract versioned models and separate GA models from fallback models."""

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
import re
import gzip
from datetime import datetime, timezone
from concurrent.futures import ThreadPoolExecutor
from urllib.request import urlopen, Request
from pathlib import Path
from adhesome_core import table, hmm_records, pf, digest, save_json, VERSION
from adhesome_components import load_crosswalk

def build(sheet, source, output, full_pfam=False, crosswalk=None):
    wanted = {a for row in table(sheet) for a in re.findall(r'PF\d{5}', row['Pfam accession'])}
    if not wanted:
        raise ValueError('Empty Pfam manifest')
    if crosswalk:
        load_crosswalk(crosswalk,sheet)
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
    output = _layout_resource(Path(output))
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(''.join(blocks), encoding='utf-8')
    _layout_resource(Path(str(output) + '.ga.hmm')).write_text(''.join(ga), encoding='utf-8')
    _layout_resource(Path(str(output) + '.fallback.hmm')).write_text(''.join(fallback), encoding='utf-8')
    save_json(str(output) + '.metadata.json', dict(pipeline_version=VERSION,source=str(source), source_sha256=digest(source),
              manifest_sha256=digest(sheet), library_sha256=digest(output), full_pfam=full_pfam,
              crosswalk_sha256=digest(crosswalk) if crosswalk else None,
              split_sha256={kind:digest(str(output)+'.'+kind+'.hmm') for kind in ('ga','fallback')}, models=found))
    return found


def download_models(sheet, folder):
    """Fetch all requested models together; retain version/hash/URL provenance."""
    wanted=sorted({a for row in table(sheet) for a in re.findall(r'PF\d{5}',row['Pfam accession'])})
    folder=_layout_resource(Path(folder))
    folder.mkdir(parents=True,exist_ok=True)
    def fetch(accession):
        url=f'https://www.ebi.ac.uk/interpro/api/entry/pfam/{accession}/?annotation=hmm'
        path=folder/(accession+'.hmm')
        if not _layout_resource(path).exists():
            with urlopen(Request(url,headers={'User-Agent':'adhesome-workflow/1.0'}),timeout=60) as response:
                data=response.read()
            if data[:2]==b'\x1f\x8b':
                data=gzip.decompress(data)
            text=data.decode('utf-8')
            if not text.startswith('HMMER3/') or not text.rstrip().endswith('//'):
                raise ValueError('Invalid downloaded HMM: '+accession)
            path.write_text(text,encoding='utf-8')
        records=list(hmm_records(path))
        if len(records)!=1 or pf(records[0][0]['ACC'][0])!=accession:
            raise ValueError('Downloaded accession mismatch: '+accession)
        return dict(accession=accession,model_version=records[0][0]['ACC'][0],url=url,sha256=digest(path))
    with ThreadPoolExecutor(max_workers=4) as pool:
        audit=list(pool.map(fetch,wanted))
    source=_layout_resource(folder/'combined_source.hmm')
    source.write_text(''.join((_layout_resource(folder/(a+'.hmm'))).read_text(encoding='utf-8') for a in wanted),encoding='utf-8')
    save_json(_layout_resource(folder/'download_provenance.json'),dict(pipeline_version=VERSION,
        retrieved_utc=datetime.now(timezone.utc).isoformat(),models=audit,manifest_sha256=digest(sheet)))
    return source

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--version',action='version',version='adhesome-library '+VERSION)
    parser.add_argument('--sheet', default=str(_layout_resource(Path(__file__).resolve().parent / 'config/adhesome_pfam_manifest.csv')))
    source=parser.add_mutually_exclusive_group(required=True)
    source.add_argument('--pfam-hmm')
    source.add_argument('--download',action='store_true',help='Fetch requested models from EMBL-EBI InterPro')
    parser.add_argument('--download-dir',type=Path,default=_layout_resource(Path(__file__).resolve().parent/'databases/pfam/source_v1.0'))
    parser.add_argument('--crosswalk',default=str(_layout_resource(Path(__file__).resolve().parent/'config/adhesome_component_crosswalk.csv')))
    parser.add_argument('--out', default=str(_layout_resource(Path(__file__).resolve().parent/'databases/pfam/adhesome_v1.0.hmm')))
    parser.add_argument('--full-pfam', action='store_true', help='Preserve all Pfam models for complete observed architectures')
    args = parser.parse_args()
    if args.download and args.full_pfam:
        parser.error('--full-pfam requires a complete local --pfam-hmm source')
    models = build(args.sheet, download_models(args.sheet,args.download_dir) if args.download else args.pfam_hmm,
                   args.out,args.full_pfam,args.crosswalk)
    print(f'Extracted {len(models)} models; scan stage will run hmmpress.')

if __name__ == '__main__':
    main()
