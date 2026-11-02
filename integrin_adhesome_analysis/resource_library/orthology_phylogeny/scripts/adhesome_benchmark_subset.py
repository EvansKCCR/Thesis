#!/usr/bin/env python3
"""Explicit opt-in entry point for the curated adhesome-subset benchmark (1.0).

Shares the original Open_Orthobench scoring equations through benchmark_workflow.
This entry point defaults to adhesome scope and cannot select the all-RefOG mode.
"""

# Resolve source resources independently of the caller working directory.
import sys as _layout_sys
from pathlib import Path as _LayoutPath
_LAYOUT_BUNDLE_ROOT = next(p for p in _LayoutPath(__file__).resolve().parents
                          if (p/'config/layout_paths.json').is_file() and (p/'pipeline_paths.py').is_file())
if str(_LAYOUT_BUNDLE_ROOT) not in _layout_sys.path:
    _layout_sys.path.insert(0, str(_LAYOUT_BUNDLE_ROOT))
from pipeline_paths import activate as _layout_activate, resource_path as _layout_resource
_layout_activate()

import sys

from benchmark_workflow import main as workflow_main


def main(argv=None):
    return workflow_main(argv, subset_only=True)


if __name__ == '__main__':
    sys.exit(main())
