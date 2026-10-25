"""Improved defaults for new runs; saved checkpoints keep their own profile."""
import argparse
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
import adhesome_workflow as core
from stage2_evidence import default_profile,profile


class MainDefaultTests(unittest.TestCase):
    def test_all_new_cli_routes_select_improved_profile_without_analysis(self):
        for name,flags,action in [
            ('resume',['--handoff','fixture'],'resume'),
            ('continue',['--annotation','fixture'],'downstream'),
            ('audit',['--annotation','fixture','--orthology','fixture'],'audit')]:
            with patch.object(core,action) as called:
                self.assertEqual(core.main([name,*flags,'--out','unused_output']),0)
                self.assertEqual(called.call_args.args[0].stage2_profile,'calibrated_v1')
        self.assertEqual(default_profile(),'calibrated_v1')
        self.assertEqual(profile({})[0],'calibrated_v1')

    def test_explicit_original_remains_available(self):
        with patch.object(core,'audit') as called:
            core.main(['audit','--annotation','fixture','--orthology','fixture','--out','unused',
                       '--stage2-profile','original'])
            self.assertEqual(called.call_args.args[0].stage2_profile,'original')

    def test_recovery_keeps_explicit_and_historical_profiles(self):
        for saved,expected in [(None,'original'),('original','original'),('calibrated_v1','calibrated_v1')]:
            with tempfile.TemporaryDirectory() as directory:
                settings={'annotation':str(Path(directory)/'annotation')}
                if saved:settings['stage2_profile']=saved
                Path(directory,'recovery.json').write_text(json.dumps(dict(settings=settings,stages={})))
                with patch.object(core,'downstream') as called:
                    core.recover(argparse.Namespace(run=Path(directory),threads=None,analysis_threads=None,iqtree_threads=None))
                self.assertEqual(called.call_args.args[0].stage2_profile,expected)


if __name__=='__main__':unittest.main()
