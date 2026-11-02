"""Compatibility import; canonical implementation is main/pipeline_paths.py."""
from pathlib import Path as _CompatPath
_compat_source = _CompatPath(__file__).resolve().parent / 'main/pipeline_paths.py'
__file__ = str(_compat_source)
exec(compile(_compat_source.read_bytes(), str(_compat_source), 'exec'))
