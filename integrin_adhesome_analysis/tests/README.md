# Synthetic regression tests

From the maintained pipeline root, in its existing Python environment:

```bash
python -B -m unittest discover -s tests -p 'test_*.py'
```

`workflow/` covers orchestrator defaults, recovery, evidence joins, reporting,
benchmark integration and data isolation. `audit/` covers family architecture,
CDD/Pfam observations, confidence, exclusions and identity naming. These tests
use small temporary fixtures and simulate external inference; they do not run
stored target/calibration/validation analyses.

Component-owned tests remain with detection and orthology components when they
test internal component APIs. CI runs the central `tests/` collection. The
default-profile regressions verify both improved new-run defaults and explicit
historical/recovery compatibility.

Component suites use their own working directories. From the bundle root:

```bash
(cd detection && python -B -m unittest discover -s tests -p 'test_*.py')
(cd motif_annotation && python -B -m unittest discover -s test -p 'test_*.py')
(cd orthology_phylogeny/scripts && python -B -m unittest discover -p 'test_*.py')
```

The layout regressions also check mixed-panel input discovery, historical model
resource keys and coverage of relocated modules in the frozen-code inventory.
