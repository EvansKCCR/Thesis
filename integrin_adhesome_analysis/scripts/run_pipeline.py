#!/usr/bin/env python3
"""Version 1.0: selected analysis plus mandatory isolated calibration self-test."""

# Resolve the organized bundle independently of the caller working directory.
import sys as _bundle_sys
from pathlib import Path as _BundlePath
_BUNDLE_ROOT=next(p for p in _BundlePath(__file__).resolve().parents if (p/'pipeline_paths.py').is_file())
if str(_BUNDLE_ROOT) not in _bundle_sys.path:_bundle_sys.path.insert(0,str(_BUNDLE_ROOT))
from pipeline_paths import activate as _activate_bundle_paths
_activate_bundle_paths()
import argparse
from pathlib import Path
import sys
from schisto_adhesome.datasets import ROOT
from schisto_adhesome.workflows import run


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--version',action='version',version='schisto-adhesome 1.0')
    parser.add_argument('--mode',choices=['validation','target'],required=True)
    parser.add_argument('--config',type=Path,default=ROOT/'configs/isolated_workflow.json')
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--reference-raw',type=Path,default=ROOT/'inputs/reference_raw')
    parser.add_argument('--input-dir',type=Path)
    parser.add_argument('--targets-config',type=Path)
    parser.add_argument('--reuse-target-detection',type=Path,help='Reuse only sequence/config-verified completed production species; controls are omitted')
    parser.add_argument('--compare-calibration',type=Path,help='Previous isolated calibration reports directory for assignment/metric drift diagnostics')
    for kind in ('calibration','validation'):
        parser.add_argument('--'+kind+'-fasta',type=Path,default=ROOT/f'validation/{kind}.faa')
        parser.add_argument('--'+kind+'-labels',type=Path,default=ROOT/('validation/calibration_set.xlsx' if kind=='calibration' else 'validation/validation_set.csv'))
    for kind in ('calibration','validation','target'):
        parser.add_argument('--'+kind+'-context',type=Path,help='External evidence directory for this branch only')
    parser.add_argument('--select-labelled-validation',action='store_true',help='Select label-table IDs and report extra FASTA records as excluded; never alter source FASTA')
    parser.add_argument('--allow-missing-context',action='store_true',help='Explicitly continue with missing external context, recording that limitation')
    parser.add_argument('--threads',type=int,default=4,help='Total budget across both branches')
    parser.add_argument('--schedule',choices=['parallel','serial'],default='parallel')
    parser.add_argument('--preflight-only',action='store_true')
    parser.add_argument('--resume',action='store_true')
    args=parser.parse_args(argv)
    if args.threads<1:parser.error('--threads must be positive')
    try:return run(args)
    except (OSError,ValueError,KeyError) as exc:
        print('ERROR: '+str(exc),file=sys.stderr);return 1


if __name__=='__main__':
    sys.exit(main())
