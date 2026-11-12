"""Version 1.0: reuse cached detection; optionally rerun only BLAST against a new panel."""

# Resolve source resources independently of the caller working directory.
import sys as _layout_sys
from pathlib import Path as _LayoutPath
_LAYOUT_BUNDLE_ROOT = next(p for p in _LayoutPath(__file__).resolve().parents
                          if (p/'config/layout_paths.json').is_file() and (p/'pipeline_paths.py').is_file())
if str(_LAYOUT_BUNDLE_ROOT) not in _layout_sys.path:
    _layout_sys.path.insert(0, str(_LAYOUT_BUNDLE_ROOT))
from pipeline_paths import activate as _layout_activate, resource_path as _layout_resource, default_thread_count as _default_thread_count
_layout_activate()

import argparse
import json
from pathlib import Path
import shutil
import sys

ROOT = _layout_resource(Path(__file__)).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from adhesome_core import digest, save_json
from adhesome_classification import classify_comparative
from evidence_coverage import search_coverage


def refresh(args):
    source, output = _layout_resource(args.detection).resolve(), _layout_resource(args.out).resolve()
    if _layout_resource(output).exists() or source == output or source in output.parents:
        raise ValueError('Choose a new output directory outside the source detection directory')
    reference_set = getattr(args, 'reference_set', None)
    if reference_set:
        reference_set = _layout_resource(reference_set).resolve()
        lock = json.loads((_layout_resource(reference_set/'reference_lock.json')).read_text(encoding='utf-8'))
        for name, expected in lock['files'].items():
            if digest(reference_set/name) != expected:
                raise ValueError('New reference panel checksum mismatch: '+name)
    batch_path = _layout_resource(source/'batch.json')
    batch = json.loads(_layout_resource(batch_path).read_text(encoding='utf-8'))
    targets = batch.get('species', {})
    if not targets or any(row.get('status') not in ('completed_with_homology', 'completed_domain_screen') for row in targets.values()):
        raise ValueError('Refresh requires completed detection for every target in batch.json')
    inputs = []
    for code, record in sorted(targets.items()):
        if not code or _layout_resource(Path(code)).name != code or '/' in code or '\\' in code or code in ('.', '..'):
            raise ValueError('Unsafe target code in batch.json')
        folder = source/code
        proteome = _layout_resource(folder/'prepared'/f'{code}.faa')
        domains = _layout_resource(folder/'domains/domains.csv')
        meta = json.loads((_layout_resource(folder/'candidates/classification.json')).read_text(encoding='utf-8'))
        if meta.get('mode', 'comparative') != 'comparative':
            raise ValueError('Refresh requires comparative classification: '+code)
        homology = _layout_resource(folder/'homology/homology.csv') if record['status'] == 'completed_with_homology' else None
        for path in [proteome, domains, *([homology] if homology else [])]:
            expected = {sha for name, sha in meta['checksums'].items() if name.replace('\\','/').split('/')[-1] == path.name}
            if expected != {digest(path)}:
                raise ValueError('Cached classification input missing or changed: '+str(path))
        search_coverage(homology, proteome)
        inputs.append((code, folder, homology is not None))
    output.mkdir(parents=True)
    batch['refreshed_from'] = str(source)
    batch['refresh_complete'] = False
    for record in targets.values():
        record.update(status='pending_refresh')
    save_json(_layout_resource(output/'batch.json'), batch)
    for code, folder, has_homology in inputs:
        destination = output/code
        destination.mkdir()
        for name in ('prepared', 'domains', *(['homology'] if has_homology and not reference_set else [])):
            shutil.copytree(folder/name, destination/name)
        if (_layout_resource(folder/'proteome.json')).exists():
            shutil.copy2(_layout_resource(folder/'proteome.json'), _layout_resource(destination/'proteome.json'))
        targets[code].update(output=str(destination), status='refreshing')
        save_json(_layout_resource(output/'batch.json'), batch)
        try:
            if reference_set:
                from adhesome_pipeline import homology
                homology(argparse.Namespace(proteome=str(_layout_resource(destination/'prepared'/f'{code}.faa')),
                    references=str(_layout_resource(reference_set/'search_references.faa')),catalogue=str(_layout_resource(reference_set/'search_catalogue.csv')),
                    metazoan=str(_layout_resource(reference_set/'metazoan_references.faa')),metazoan_catalogue=str(_layout_resource(reference_set/'metazoan_catalogue.csv')),
                    out=_layout_resource(destination/'homology'),cpu=getattr(args,'cpu',4),evalue=1e-5,coverage=.3,tie_fraction=.01))
            classify_comparative(argparse.Namespace(proteome=_layout_resource(destination/'prepared'/f'{code}.faa'),
                domains=_layout_resource(destination/'domains/domains.csv'), homology=_layout_resource(destination/'homology/homology.csv') if has_homology or reference_set else None,
                evidence=[], rules=args.rules, manifest=args.manifest, crosswalk=args.crosswalk, out=_layout_resource(destination/'candidates')))
        except Exception as exc:
            targets[code].update(status='failed', error=str(exc))
            save_json(_layout_resource(output/'batch.json'), batch)
            raise
        targets[code].update(status='completed_with_homology' if has_homology or reference_set else 'completed_domain_screen')
        print('Refreshed cached evidence: '+code, flush=True)
        save_json(_layout_resource(output/'batch.json'), batch)
    batch['refresh_complete'] = True
    save_json(_layout_resource(output/'batch.json'), batch)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--detection', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--reference-set',type=Path,help='Explicitly rerun forward/reverse BLAST with this paired panel, reusing existing prepared FASTA and domain scans')
    parser.add_argument('--cpu',type=int,default=_default_thread_count(4),
                        help='CPU budget; defaults to the scheduler allocation or up to 4 local CPUs')
    for name, filename in [('rules','family_rules.json'), ('manifest','adhesome_pfam_manifest.csv'),
                           ('crosswalk','adhesome_component_crosswalk.csv')]:
        parser.add_argument('--'+name, type=Path, default=_layout_resource(ROOT/'config')/filename)
    args = parser.parse_args()
    if args.cpu<1:parser.error('--cpu must be positive')
    try:
        refresh(args)
    except (ValueError, OSError) as exc:
        parser.exit(1, 'ERROR: '+str(exc)+'\n')


if __name__ == '__main__':
    main()
