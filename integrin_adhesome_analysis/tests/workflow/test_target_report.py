"""Saved target summaries exclude requested controls without changing audits."""

# Resolve source resources independently of the caller working directory.
import sys as _layout_sys
from pathlib import Path as _LayoutPath
_LAYOUT_BUNDLE_ROOT = next(p for p in _LayoutPath(__file__).resolve().parents
                          if (p/'config/layout_paths.json').is_file() and (p/'pipeline_paths.py').is_file())
if str(_LAYOUT_BUNDLE_ROOT) not in _layout_sys.path:
    _layout_sys.path.insert(0, str(_LAYOUT_BUNDLE_ROOT))
from pipeline_paths import activate as _layout_activate, resource_path as _layout_resource
_layout_activate()


import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT=_LAYOUT_BUNDLE_ROOT;sys.path.insert(0,str(ROOT))
from schisto_adhesome import reports
from schisto_adhesome.datasets import table
import adhesome_workflow as core


class TargetReportTests(unittest.TestCase):
    def fixture(self,directory):
        audit=_layout_resource(Path(directory)/'audit');audit.mkdir()
        table(_layout_resource(audit/'audited_family_assignments.tsv'),[
            dict(protein_id=code+'__p',Specie=code,putative_audited_family='profilin',grade='A',assignment_confidence='strong')
            for code in ('NovelTarget','validation')])
        (_layout_resource(audit/'Stage_2_output.xlsx')).write_bytes(b'saved-workbook-fixture')
        core.seal(audit,'audit',stage2_profile='calibrated_v1')
        return audit,_layout_resource(Path(directory)/'report')

    def test_report_only_excludes_controls_and_preserves_all_audit_bytes(self):
        with tempfile.TemporaryDirectory() as directory:
            audit,out=self.fixture(directory);before={p:_layout_resource(p).read_bytes() for p in audit.iterdir()}
            self.assertEqual(reports.main(['--audit',str(audit),'--out',str(out),'--exclude-species','validation']),0)
            summary=json.loads((_layout_resource(out/'summary.json')).read_text())
            self.assertEqual(summary['candidates'],1);self.assertEqual(summary['excluded_candidate_rows'],1)
            self.assertEqual(reports.read(_layout_resource(out/'candidate_counts.tsv'))[0]['species'],'NovelTarget')
            for path,value in before.items():self.assertEqual(_layout_resource(path).read_bytes(),value)

    def test_default_api_retains_every_stored_species(self):
        with tempfile.TemporaryDirectory() as directory:
            audit,out=self.fixture(directory)
            self.assertEqual(reports.target_report(audit,out)['candidates'],2)

    def test_refuses_existing_or_sealed_output_and_altered_audit(self):
        with tempfile.TemporaryDirectory() as directory:
            audit,out=self.fixture(directory);out.mkdir()
            self.assertEqual(reports.main(['--audit',str(audit),'--out',str(out)]),1)
            self.assertEqual(reports.main(['--audit',str(audit),'--out',str(_layout_resource(audit/'reports'))]),1)
            (_layout_resource(audit/'audited_family_assignments.tsv')).write_text('changed')
            self.assertEqual(reports.main(['--audit',str(audit),'--out',str(_layout_resource(Path(directory)/'new'))]),1)


if __name__=='__main__':unittest.main()
