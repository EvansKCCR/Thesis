"""Run the final pipeline from original proteomes without historical outputs."""
import argparse
from pathlib import Path
import subprocess
import sys


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cpu', type=int, default=8)
    parser.add_argument('--config', type=Path, required=True, help='Target species and explicit gene-ID extraction rules')
    parser.add_argument('--input-dir', type=Path, default=root / 'proteome_fasta_files')
    parser.add_argument('--out', type=Path, default=root / 'results/species_reference_v5')
    args = parser.parse_args()
    if args.cpu < 1:
        parser.error('--cpu must be positive')
    subprocess.run([
        sys.executable, str(root / 'adhesome_pipeline.py'), 'batch',
        '--input-dir', str(args.input_dir.resolve()),
        '--config', str(args.config.resolve()),
        '--library', str(root / 'databases/pfam/adhesome_v4.hmm'),
        '--reference-set', str(root / 'references/schistosome_v5'),
        '--rules', str(root / 'config/family_rules.json'),
        '--out', str(args.out.resolve()), '--cpu', str(args.cpu),
    ], check=True, cwd=root)


if __name__ == '__main__':
    main()
