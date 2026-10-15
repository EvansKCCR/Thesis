"""Version 1.0: reuse cached detection; optionally rerun only BLAST against a new panel."""
import argparse
import json
from pathlib import Path
import shutil
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from adhesome_core import digest, save_json
from adhesome_classification import classify_comparative
from evidence_coverage import search_coverage


def refresh(args):
    source, output = args.detection.resolve(), args.out.resolve()
    if output.exists() or source == output or source in output.parents:
        raise ValueError('Choose a new output directory outside the source detection directory')
    reference_set = getattr(args, 'reference_set', None)
    if reference_set:
        reference_set = reference_set.resolve()
        lock = json.loads((reference_set/'reference_lock.json').read_text(encoding='utf-8'))
        for name, expected in lock['files'].items():
            if digest(reference_set/name) != expected:
                raise ValueError('New reference panel checksum mismatch: '+name)
    batch_path = source/'batch.json'
    batch = json.loads(batch_path.read_text(encoding='utf-8'))
    targets = batch.get('species', {})
    if not targets or any(row.get('status') not in ('completed_with_homology', 'completed_domain_screen') for row in targets.values()):
        raise ValueError('Refresh requires completed detection for every target in batch.json')
    inputs = []
    for code, record in sorted(targets.items()):
        if not code or Path(code).name != code or '/' in code or '\\' in code or code in ('.', '..'):
            raise ValueError('Unsafe target code in batch.json')
        folder = source/code
        proteome = folder/'prepared'/f'{code}.faa'
        domains = folder/'domains/domains.csv'
        meta = json.loads((folder/'candidates/classification.json').read_text(encoding='utf-8'))
        if meta.get('mode', 'comparative') != 'comparative':
            raise ValueError('Refresh requires comparative classification: '+code)
        homology = folder/'homology/homology.csv' if record['status'] == 'completed_with_homology' else None
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
    save_json(output/'batch.json', batch)
    for code, folder, has_homology in inputs:
        destination = output/code
        destination.mkdir()
        for name in ('prepared', 'domains', *(['homology'] if has_homology and not reference_set else [])):
            shutil.copytree(folder/name, destination/name)
        if (folder/'proteome.json').exists():
            shutil.copy2(folder/'proteome.json', destination/'proteome.json')
        targets[code].update(output=str(destination), status='refreshing')
        save_json(output/'batch.json', batch)
        try:
            if reference_set:
                from adhesome_pipeline import homology
                homology(argparse.Namespace(proteome=str(destination/'prepared'/f'{code}.faa'),
                    references=str(reference_set/'search_references.faa'),catalogue=str(reference_set/'search_catalogue.csv'),
                    metazoan=str(reference_set/'metazoan_references.faa'),metazoan_catalogue=str(reference_set/'metazoan_catalogue.csv'),
                    out=destination/'homology',cpu=getattr(args,'cpu',4),evalue=1e-5,coverage=.3,tie_fraction=.01))
            classify_comparative(argparse.Namespace(proteome=destination/'prepared'/f'{code}.faa',
                domains=destination/'domains/domains.csv', homology=destination/'homology/homology.csv' if has_homology or reference_set else None,
                evidence=[], rules=args.rules, manifest=args.manifest, crosswalk=args.crosswalk, out=destination/'candidates'))
        except Exception as exc:
            targets[code].update(status='failed', error=str(exc))
            save_json(output/'batch.json', batch)
            raise
        targets[code].update(status='completed_with_homology' if has_homology or reference_set else 'completed_domain_screen')
        print('Refreshed cached evidence: '+code, flush=True)
        save_json(output/'batch.json', batch)
    batch['refresh_complete'] = True
    save_json(output/'batch.json', batch)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--detection', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--reference-set',type=Path,help='Explicitly rerun forward/reverse BLAST with this paired panel, reusing existing prepared FASTA and domain scans')
    parser.add_argument('--cpu',type=int,default=4)
    for name, filename in [('rules','family_rules.json'), ('manifest','adhesome_pfam_manifest.csv'),
                           ('crosswalk','adhesome_component_crosswalk.csv')]:
        parser.add_argument('--'+name, type=Path, default=ROOT/'config'/filename)
    args = parser.parse_args()
    if args.cpu<1:parser.error('--cpu must be positive')
    try:
        refresh(args)
    except (ValueError, OSError) as exc:
        parser.exit(1, 'ERROR: '+str(exc)+'\n')


if __name__ == '__main__':
    main()
