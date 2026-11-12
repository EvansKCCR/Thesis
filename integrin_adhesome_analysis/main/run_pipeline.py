#!/usr/bin/env python3
"""Version 1.0: selected analysis and isolated self-test when its inputs exist."""

# Resolve source resources independently of the caller working directory.
import sys as _layout_sys
from pathlib import Path as _LayoutPath
_LAYOUT_BUNDLE_ROOT = next(p for p in _LayoutPath(__file__).resolve().parents
                          if (p/'config/layout_paths.json').is_file() and (p/'pipeline_paths.py').is_file())
if str(_LAYOUT_BUNDLE_ROOT) not in _layout_sys.path:
    _layout_sys.path.insert(0, str(_LAYOUT_BUNDLE_ROOT))
from pipeline_paths import (activate as _layout_activate, resource_path as _layout_resource,
                            default_thread_count as _default_thread_count,
                            warn_slow_wsl_mounts as _warn_slow_wsl_mounts)
_layout_activate()


import argparse
from pathlib import Path
import sys
from schisto_adhesome.datasets import ROOT


def main(argv=None):
    argv=list(sys.argv[1:] if argv is None else argv)
    if '--dry-run' in argv:
        from deployment_validation import main as validate_deployment
        return validate_deployment(argv)
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--version',action='version',version='schisto-adhesome 1.0')
    parser.add_argument('--mode',choices=['validation','target'],required=True)
    parser.add_argument('--config',type=Path,default=_layout_resource(ROOT/'configs/isolated_workflow.json'))
    parser.add_argument('--out',type=Path,required=True)
    from input_dependencies import raw_input
    parser.add_argument('--reference-raw',type=Path,default=raw_input('reference_raw'))
    parser.add_argument('--reference-file',dest='reference_files',action='append',help='Select relative raw reference FASTA paths; repeat')
    parser.add_argument('--reference-panel',help='Named reference file list in config/input_registry.json')
    parser.add_argument('--input-dir',type=Path,default=raw_input('target_raw'))
    parser.add_argument('--targets-config',type=Path,default=raw_input('target_configuration'))
    parser.add_argument('--reuse-target-detection',type=Path,help='Reuse only sequence/config-verified completed production species; controls are omitted')
    parser.add_argument('--compare-calibration',type=Path,help='Previous isolated calibration reports directory for assignment/metric drift diagnostics')
    for kind in ('calibration','validation'):
        parser.add_argument('--'+kind+'-fasta',type=Path,default=_layout_resource(ROOT/f'validation/{kind}.faa'))
        parser.add_argument('--'+kind+'-labels',type=Path,default=ROOT/('validation/calibration_set.xlsx' if kind=='calibration' else 'validation/validation_set.csv'))
    for kind in ('calibration','validation','target'):
        parser.add_argument('--'+kind+'-context',type=Path,help='External evidence directory for this branch only')
    parser.add_argument('--phylogenetic-identity-resolution',action='store_true',help='Optional existing-tree concordance review; family audits complete without it')
    for kind in ('calibration','validation','target'):
        parser.add_argument('--'+kind+'-phylogenetic-evidence',type=Path,help='Branch-specific supported-clade TSV; never shared between evaluation branches')
    parser.add_argument('--select-labelled-validation',action='store_true',help='Select label-table IDs and report extra FASTA records as excluded; never alter source FASTA')
    parser.add_argument('--allow-missing-context',action='store_true',help='Explicitly continue with missing external context, recording that limitation')
    parser.add_argument('--threads',type=int,default=_default_thread_count(4),
                        help='Total CPU budget across both branches; defaults to the scheduler allocation or up to 4 local CPUs')
    parser.add_argument('--phylogeny',action='store_true',help='Target mode only: run optional candidate phylogenies after orthology and audit')
    parser.add_argument('--phylogeny-grouping',choices=['orthogroup','family'],default='orthogroup',help='Target phylogeny grouping')
    parser.add_argument('--schedule',choices=['parallel','serial'],default='parallel')
    parser.add_argument('--preflight-only',action='store_true')
    parser.add_argument('--resume',action='store_true')
    args=parser.parse_args(argv)
    _warn_slow_wsl_mounts(value for value in vars(args).values() if isinstance(value,Path))
    if args.threads<1:parser.error('--threads must be positive')
    if args.mode=='validation' and args.phylogeny:
        parser.error('--phylogeny is a target-only option; validation/calibration skip mainstream orthology and phylogeny')
    from schisto_adhesome.workflows import run
    try:return run(args)
    except (OSError,ValueError,KeyError) as exc:
        print('ERROR: '+str(exc),file=sys.stderr);return 1


if __name__=='__main__':
    sys.exit(main())
