"""Version 1.0 detection: prepare, scan, homology and candidate hypotheses."""
import argparse
from pathlib import Path
import subprocess
import sys


def main():
    root=Path(__file__).resolve().parents[1]
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--version',action='version',version='adhesome-detection 1.0')
    parser.add_argument('--cpu',type=int,default=8)
    parser.add_argument('--config',type=Path,required=True,help='Target species and explicit gene-ID extraction rules')
    parser.add_argument('--input-dir',type=Path,default=root/'proteome_fasta_files')
    parser.add_argument('--out',type=Path,default=root/'results/detection_v1.0')
    parser.add_argument('--check-inputs',action='store_true',help='Validate target inputs without running detection')
    for name,path in [('library','databases/pfam/adhesome_v1.0.hmm'),
                      ('manifest','config/adhesome_pfam_manifest.csv'),
                      ('crosswalk','config/adhesome_component_crosswalk.csv'),
                      ('reference-set','references/adhesome_v1.0'),('rules','config/family_rules.json')]:
        parser.add_argument('--'+name,type=Path,default=root/path)
    args=parser.parse_args()
    if args.cpu<1:
        parser.error('--cpu must be positive')
    command=[sys.executable,str(root/'adhesome_pipeline.py'),'batch',
             '--input-dir',str(args.input_dir.resolve()),'--config',str(args.config.resolve()),
             '--out',str(args.out.resolve()),'--cpu',str(args.cpu)]
    for flag in ('library','manifest','crosswalk','reference_set','rules'):
        command+=['--'+flag.replace('_','-'),str(getattr(args,flag).resolve())]
    if args.check_inputs:command+=['--check-inputs']
    subprocess.run(command,check=True,cwd=root)


if __name__=='__main__':
    main()
