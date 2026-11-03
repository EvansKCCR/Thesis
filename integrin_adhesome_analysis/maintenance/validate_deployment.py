#!/usr/bin/env python3
"""Run the dataset-free deployment audit from any working directory."""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'resource_library/scripts'))
from deployment_validation import main

if __name__ == '__main__':
    raise SystemExit(main())
