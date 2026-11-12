"""Version 1.0 detection: prepare, scan, homology and candidate hypotheses."""

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
import subprocess
import sys
from input_dependencies import raw_input


def main():
    root=_layout_resource(Path(__file__)).resolve().parents[1]
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--version',action='version',version='adhesome-detection 1.0')
    parser.add_argument('--cpu',type=int,default=_default_thread_count(8),
                        help='CPU budget; defaults to the scheduler allocation or up to 8 local CPUs')
    parser.add_argument('--config',type=Path,required=True,help='Target species and explicit gene-ID extraction rules')
    parser.add_argument('--input-dir',type=Path,default=raw_input('target_raw'))
    parser.add_argument('--out',type=Path,default=_layout_resource(root/'results/detection_v1.0'))
    parser.add_argument('--check-inputs',action='store_true',help='Validate target inputs without running detection')
    parser.add_argument('--species',action='append',help='Process only this configured target; repeat to select several and preserve other target results')
    parser.add_argument('--resume',action='store_true',help='Continue in place without backups; verify and retain completed targets')
    parser.add_argument('--exclude-target',action='append',help='Exclude an unfinished target from the active batch')
    for name,path in [('library','databases/pfam/adhesome_v1.0.hmm'),
                      ('manifest','config/adhesome_pfam_manifest.csv'),
                      ('crosswalk','config/adhesome_component_crosswalk.csv'),
                      ('reference-set','references/adhesome_v1.0'),('rules','config/family_rules.json')]:
        parser.add_argument('--'+name,type=Path,default=root/path)
    args=parser.parse_args()
    _warn_slow_wsl_mounts(value for value in vars(args).values() if isinstance(value,Path))
    if args.cpu<1:
        parser.error('--cpu must be positive')
    command=[sys.executable,str(_layout_resource(root/'adhesome_pipeline.py')),'batch',
             '--input-dir',str(_layout_resource(args.input_dir).resolve()),'--config',str(_layout_resource(args.config).resolve()),
             '--out',str(_layout_resource(args.out).resolve()),'--cpu',str(args.cpu)]
    for flag in ('library','manifest','crosswalk','reference_set','rules'):
        command+=['--'+flag.replace('_','-'),str(_layout_resource(getattr(args,flag)).resolve())]
    if args.check_inputs:command+=['--check-inputs']
    if args.resume:command+=['--resume']
    for code in args.exclude_target or []:command+=['--exclude-target',code]
    for code in args.species or []:command+=['--species',code]
    subprocess.run(command,check=True,cwd=root)


if __name__=='__main__':
    main()
