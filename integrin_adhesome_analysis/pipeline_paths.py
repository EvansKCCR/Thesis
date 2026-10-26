"""Canonical paths for the organized version 1.0 source bundle."""
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parent
SCRIPTS=ROOT/'scripts'
DETECTION=ROOT/'detection'
VALIDATION=ROOT/'validation'
MOTIFS=ROOT/'motif_annotation'
ORTHOLOGY=ROOT/'orthology_phylogeny'
ANNOTATION=ROOT/'annotation'
LEGACY_FOLDERS={
    '01_detection':'detection','02_calibration_set':'validation',
    '03_functional_motif_annotation':'motif_annotation',
    '04_orthology_phylogenetic_classification':'orthology_phylogeny',
    '05_conservative_adhesome_family_assignment_audit':'annotation'}


def activate():
    """Expose shared modules for root, component, tool and scripts/ entry points."""
    for folder in (ROOT,SCRIPTS):
        if str(folder) not in sys.path:sys.path.insert(0,str(folder))
    return ROOT


def recorded_resource(name):
    """Resolve old relative model-input keys without rewriting historical models.

    Existing absolute paths stay exact. This is only a source-library relocation
    adapter, never a rebase or repair of a saved biological run.
    """
    name=str(name).replace('\\','/');path=Path(name)
    if path.is_absolute():return path
    direct=ROOT/path
    if direct.exists():return direct
    parts=path.parts
    if parts and parts[0] in LEGACY_FOLDERS:
        return ROOT/LEGACY_FOLDERS[parts[0]]/Path(*parts[1:])
    if len(parts)==1 and (SCRIPTS/path).is_file():return SCRIPTS/path
    return direct
