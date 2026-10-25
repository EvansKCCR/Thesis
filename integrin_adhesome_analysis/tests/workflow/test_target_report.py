"""Saved target summaries exclude requested controls without changing audits."""
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from schisto_adhesome import reports
from schisto_adhesome.datasets import table
import adhesome_workflow as core


class TargetReportTests(unittest.TestCase):
    def fixture(self,directory):
        audit=Path(directory)/'audit';audit.mkdir()
        table(audit/'audited_family_assignments.tsv',[
            dict(protein_id=code+'__p',Specie=code,putative_audited_family='profilin',grade='A',assignment_confidence='strong')
            for code in ('NovelTarget','validation')])
        (audit/'Stage_2_output.xlsx').write_bytes(b'saved-workbook-fixture')
        core.seal(audit,'audit',stage2_profile='calibrated_v1')
        return audit,Path(directory)/'report'

    def test_report_only_excludes_controls_and_preserves_all_audit_bytes(self):
        with tempfile.TemporaryDirectory() as directory:
            audit,out=self.fixture(directory);before={p:p.read_bytes() for p in audit.iterdir()}
            self.assertEqual(reports.main(['--audit',str(audit),'--out',str(out),'--exclude-species','validation']),0)
            summary=json.loads((out/'summary.json').read_text())
            self.assertEqual(summary['candidates'],1);self.assertEqual(summary['excluded_candidate_rows'],1)
            self.assertEqual(reports.read(out/'candidate_counts.tsv')[0]['species'],'NovelTarget')
            for path,value in before.items():self.assertEqual(path.read_bytes(),value)

    def test_default_api_retains_every_stored_species(self):
        with tempfile.TemporaryDirectory() as directory:
            audit,out=self.fixture(directory)
            self.assertEqual(reports.target_report(audit,out)['candidates'],2)

    def test_refuses_existing_or_sealed_output_and_altered_audit(self):
        with tempfile.TemporaryDirectory() as directory:
            audit,out=self.fixture(directory);out.mkdir()
            self.assertEqual(reports.main(['--audit',str(audit),'--out',str(out)]),1)
            self.assertEqual(reports.main(['--audit',str(audit),'--out',str(audit/'reports')]),1)
            (audit/'audited_family_assignments.tsv').write_text('changed')
            self.assertEqual(reports.main(['--audit',str(audit),'--out',str(Path(directory)/'new')]),1)


if __name__=='__main__':unittest.main()
