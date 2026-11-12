"""Canonical bundle paths and explicit source-layout compatibility (version 1.0).

Only paths recorded in config/layout_paths.json are relocated. User datasets,
saved run paths, tool commands and external installations are not rebased.
"""
from pathlib import Path, PureWindowsPath
import hashlib
import json
import os
import re
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


def allocated_cpu_count():
    """Return CPUs allocated to this process/job, then affinity, then host count."""
    for name in ('SLURM_CPUS_PER_TASK','PBS_NP','PBS_NUM_PPN','NCPUS','NSLOTS',
                 'LSB_DJOB_NUMPROC','SLURM_CPUS_ON_NODE'):
        match=re.match(r'\s*(\d+)',os.environ.get(name,''))
        if match and int(match.group(1))>0:
            allocated=int(match.group(1))
            try:
                affinity=len(os.sched_getaffinity(0))
            except (AttributeError,OSError):
                affinity=allocated
            return max(1,min(allocated,affinity))
    try:
        return max(1,len(os.sched_getaffinity(0)))
    except (AttributeError,OSError):
        return max(1,os.cpu_count() or 1)


def default_thread_count(local_max=8):
    """Use scheduler allocation on HPC; cap implicit desktop defaults locally."""
    has_allocation=any(os.environ.get(name) for name in
        ('SLURM_CPUS_PER_TASK','PBS_NP','PBS_NUM_PPN','NCPUS','NSLOTS',
         'LSB_DJOB_NUMPROC','SLURM_CPUS_ON_NODE'))
    count=allocated_cpu_count()
    return count if has_allocation else min(count,max(1,int(local_max)))


def configure_runtime_environment():
    """Use scheduler scratch and prevent nested BLAS/OpenMP oversubscription."""
    scheduler_tmp=os.environ.get('SLURM_TMPDIR')
    if not os.environ.get('TMPDIR') and scheduler_tmp and Path(scheduler_tmp).is_dir():
        os.environ['TMPDIR']=scheduler_tmp
    for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS',
                 'BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
        os.environ.setdefault(name,'1')


def warn_slow_wsl_mounts(paths):
    """Warn when large inputs/outputs are under a Windows-mounted WSL drive."""
    release=Path('/proc/sys/kernel/osrelease')
    try:
        is_wsl=bool(os.environ.get('WSL_DISTRO_NAME')) or 'microsoft' in release.read_text().casefold()
    except OSError:
        is_wsl=bool(os.environ.get('WSL_DISTRO_NAME'))
    if not is_wsl:
        return
    mounted=[]
    for value in paths:
        if value is None:
            continue
        try:
            path=Path(value).expanduser().resolve()
        except (OSError,RuntimeError):
            path=Path(value).expanduser()
        if path.is_absolute() and re.match(r'^/mnt/[a-z](?:/|$)',path.as_posix(),re.I):
            mounted.append(path)
    if mounted:
        print('WSL storage note: sequence data/output under /mnt/<drive> can be slow for OrthoFinder and many-file workloads. '
              'For faster I/O, keep the active bundle, inputs and run directory in the WSL Linux filesystem (for example under ~/work) '
              'or use fast persistent HPC scratch. Paths were not moved.',file=sys.stderr,flush=True)


def activate():
    """Make canonical CLI, shared modules and maintenance packages importable."""
    configure_runtime_environment()
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
