#!/usr/bin/env python3
"""Explicit opt-in entry point for the curated adhesome-subset benchmark (1.0).

Shares the original Open_Orthobench scoring equations through benchmark_workflow.
This entry point defaults to adhesome scope and cannot select the all-RefOG mode.
"""
import sys

from benchmark_workflow import main as workflow_main


def main(argv=None):
    return workflow_main(argv, subset_only=True)


if __name__ == '__main__':
    sys.exit(main())
