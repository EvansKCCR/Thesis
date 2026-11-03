"""Canonical bundle paths and explicit source-layout compatibility (version 1.0).

Only paths recorded in config/layout_paths.json are relocated. User datasets,
saved run paths, tool commands and external installations are not rebased.
"""
from pathlib import Path, PureWindowsPath
import hashlib
import json
import sys

ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / 'config'
MAIN = ROOT / 'main'
MAINTENANCE = ROOT / 'maintenance'
RESOURCE_LIBRARY = ROOT / 'resource_library'
SCRIPTS = RESOURCE_LIBRARY / 'scripts'
DETECTION = RESOURCE_LIBRARY / 'detection'
VALIDATION = ROOT / 'validation'
MOTIFS = RESOURCE_LIBRARY / 'motif_annotation'
ORTHOLOGY = RESOURCE_LIBRARY / 'orthology_phylogeny'
ANNOTATION = RESOURCE_LIBRARY / 'annotation'
LEGACY_FOLDERS = {
    '01_detection': 'detection', '02_calibration_set': 'validation',
    '03_functional_motif_annotation': 'motif_annotation',
    '04_orthology_phylogenetic_classification': 'orthology_phylogeny',
    '05_conservative_adhesome_family_assignment_audit': 'annotation'}
_LAYOUT = json.loads((CONFIG / 'layout_paths.json').read_text(encoding='utf-8'))
_RECEIPT_CACHE = None


def activate():
    """Make canonical CLI, shared modules and maintenance packages importable."""
    for folder in (ROOT, MAIN, SCRIPTS, MAINTENANCE, RESOURCE_LIBRARY, DETECTION, MOTIFS, ANNOTATION):
        if str(folder) not in sys.path:
            sys.path.insert(0, str(folder))
    return ROOT


def _map_key(key):
    key = str(key).replace('\\', '/').strip('/')
    original = key
    pieces = key.split('/')
    if pieces and pieces[0] in LEGACY_FOLDERS:
        key = '/'.join([LEGACY_FOLDERS[pieces[0]], *pieces[1:]])
    if key in _LAYOUT['file_map']:
        return _LAYOUT['file_map'][key]
    for source, destination in sorted(_LAYOUT['directory_map'].items(), key=lambda pair: -len(pair[0])):
        if key == source or key.startswith(source + '/'):
            return destination + key[len(source):]
    return key if key != original else None


def resource_path(value):
    """Resolve only known old source paths; preserve ordinary relative user paths."""
    path = Path(value)
    if path.is_absolute():
        normalized = path.resolve()
        try:
            key = normalized.relative_to(ROOT).as_posix()
        except ValueError:
            return path
    else:
        if PureWindowsPath(str(value)).is_absolute():
            return path
        # Normalize parent traversals only when the resulting path is in this bundle.
        normalized = path.resolve()
        try:
            absolute_key = normalized.relative_to(ROOT).as_posix()
        except ValueError:
            absolute_key = None
        key = absolute_key if absolute_key and _map_key(absolute_key) is not None else path.as_posix()
    relocated = _map_key(key)
    return ROOT / relocated if relocated is not None else path


def recorded_resource(name):
    """Resolve historical model/code keys without editing saved model receipts."""
    path = Path(str(name).replace('\\', '/'))
    relocated = resource_path(path)
    if relocated != path or path.is_absolute() or PureWindowsPath(str(name)).is_absolute():
        return relocated
    if len(path.parts) == 1 and (SCRIPTS / path).is_file():
        return SCRIPTS / path
    return ROOT / path


def recorded_digest_matches(name, expected):
    """Accept exact bytes or one specifically recorded layout-only transition.

    Unknown digests and subsequent source edits still fail. This is not a rule,
    threshold or parameter migration and never rewrites the saved contract.
    """
    global _RECEIPT_CACHE
    path = recorded_resource(name)
    if not path.is_file():
        return False
    actual = hashlib.sha256(path.read_bytes()).hexdigest()
    if actual == expected:
        return True
    receipt = MAINTENANCE / 'layout_migration.json'
    if not receipt.is_file():
        return False
    stamp = receipt.stat().st_mtime_ns
    if _RECEIPT_CACHE is None or _RECEIPT_CACHE[0] != stamp:
        data = json.loads(receipt.read_text())
        _RECEIPT_CACHE = (stamp, {entry['old_path']: entry for entry in data['files']})
    key = str(name).replace('\\', '/')
    try:
        if Path(name).is_absolute():
            key = Path(name).resolve().relative_to(ROOT).as_posix()
    except ValueError:
        return False
    entry = _RECEIPT_CACHE[1].get(key)
    return bool(entry and entry['old_sha256'] == expected and entry.get('new_sha256') == actual
                and path == ROOT / entry['new_path'])
