#!/usr/bin/env python3
"""Compatibility launcher; use main/adhesome_workflow.py for the organized CLI."""
from pathlib import Path as _CompatPath
_compat_source = _CompatPath(__file__).resolve().parent / 'main/adhesome_workflow.py'
__file__ = str(_compat_source)
exec(compile(_compat_source.read_bytes(), str(_compat_source), 'exec'))
