"""Version 1.0 read-only audit of reusable features and historical mixed inference."""

# Resolve shared pipeline imports from the installed script, not the caller cwd.
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import argparse
from pathlib import Path
from schisto_adhesome.diagnostics import audit_existing

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--runs',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--corrected-downstream',type=Path)
    a=p.parse_args();print(audit_existing(a.runs,a.out,a.corrected_downstream))
