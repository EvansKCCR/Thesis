#!/usr/bin/env python3
"""Inventory deployable source/libraries, excluding user datasets and run products."""
import argparse
import hashlib
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'integrin_adhesome_analysis/resource_library/scripts'))
from input_dependencies import deployment_source_digest
EXCLUDED={'__pycache__','.git','cache','results','runs','input','inputs','validation','validation_reports'}


def source_files(root=ROOT):
    root=Path(root).resolve()
    main=root/'integrin_adhesome_analysis'
    selected=set()
    for folder in ['config','main','maintenance','resource_library','tests','docs']:
        selected.update(p for p in (main/folder).rglob('*') if p.is_file())
    for folder in ['homology_retrieval','independent_family_phylogeny','reconciliation','.github']:
        selected.update(p for p in (root/folder).rglob('*') if p.is_file() and p.suffix.lower() in ('.py','.md','.json','.yml','.yaml'))
    for folder in (root,main):
        selected.update(p for p in folder.iterdir() if p.is_file() and p.name in
                        ('README.md','LICENSE','VERSION','.gitattributes','.gitignore','adhesome_workflow.py','pipeline_paths.py'))
    return sorted(p for p in selected if not EXCLUDED.intersection(p.relative_to(root).parts) and
                  p.suffix not in ('.pyc','.lock') and not p.name.endswith('.complete.json') and 'reference_fasta' not in p.relative_to(root).parts)


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write',action='store_true',help='Refresh source-only release_manifest.sha256; no backups or analysis')
    args=parser.parse_args(argv)
    paths=source_files()
    if args.write:
        header='# Source inventory: code/document JSON/text use canonical LF; scientific assets retain byte hashes. Validate with maintenance/validate_deployment.py.\n'
        (ROOT/'release_manifest.sha256').write_text(header+''.join(deployment_source_digest(p)+'  '+p.relative_to(ROOT).as_posix()+'\n' for p in paths),encoding='utf-8')
    print(f'{len(paths)} source/library files; '+('source release manifest refreshed' if args.write else 'plan only; no files written'))
    return 0


if __name__=='__main__':raise SystemExit(main())
