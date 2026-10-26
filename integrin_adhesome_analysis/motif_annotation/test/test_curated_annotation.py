
# Resolve the organized bundle independently of the caller working directory.
import sys as _bundle_sys
from pathlib import Path as _BundlePath
_BUNDLE_ROOT=next(p for p in _BundlePath(__file__).resolve().parents if (p/'pipeline_paths.py').is_file())
if str(_BUNDLE_ROOT) not in _bundle_sys.path:_bundle_sys.path.insert(0,str(_BUNDLE_ROOT))
from pipeline_paths import activate as _activate_bundle_paths
_activate_bundle_paths()
import unittest,tempfile,zipfile,re
from pathlib import Path
from collections import defaultdict
import curated_annotation as ca


class CuratedTests(unittest.TestCase):
    def test_topology_exact_sequences_and_tail(self):
            with tempfile.TemporaryDirectory() as d:
                p=Path(d)/'predicted_topologies.3line'
                p.write_text('>x | TM\nACDEFGHI\nOOMMMIII\n>y\nACDE\nIIII\n')
                fam={'a':[dict(id='x',sequence='ACDEFGHI'),dict(id='y',sequence='ACDF')]}
                rs=defaultdict(list)
                audit=ca.import_topology([p],fam,rs)
                self.assertEqual([r['status'] for r in audit],['exact_sequence_verified','sequence_mismatch'])
                self.assertIn(('cytoplasmic_tail',6,8),rs['a','x'])
                self.assertEqual(rs['a','y'],[])

    def test_zip_topology_without_extraction(self):
            with tempfile.TemporaryDirectory() as d:
                p=Path(d)/'top.zip'
                with zipfile.ZipFile(p,'w') as archive:
                    archive.writestr('nested/predicted_topologies.3line','>x\nACDEFGHI\nOOMMMIII\n')
                rs=defaultdict(list)
                rows=ca.import_topology([p],{'a':[dict(id='x',sequence='ACDEFGHI')]},rs)
                self.assertEqual(rows[0]['status'],'exact_sequence_verified')
                self.assertIn('!nested/predicted_topologies.3line',rows[0]['source'])
                self.assertFalse((Path(d)/'nested').exists())

    def test_full_sequence_terminal_anchor_and_unknown_residues(self):
        self.assertEqual(list(ca.occurrences(re.compile("P..P"),"PPPPP")),[(1,4,"PPPP"),(2,5,"PPPP")])
        self.assertEqual(list(ca.occurrences(re.compile("NP.Y"),"NPXY")),[])
        self.assertEqual(list(ca.occurrences(re.compile("RGD$"),"RGDA")),[])
