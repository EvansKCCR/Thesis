"""Standalone audit reports use saved calls and their own labelled populations."""

# Resolve source resources independently of the caller working directory.
import sys as _layout_sys
from pathlib import Path as _LayoutPath
_LAYOUT_BUNDLE_ROOT = next(p for p in _LayoutPath(__file__).resolve().parents
                          if (p/'config/layout_paths.json').is_file() and (p/'pipeline_paths.py').is_file())
if str(_LAYOUT_BUNDLE_ROOT) not in _layout_sys.path:
    _layout_sys.path.insert(0, str(_LAYOUT_BUNDLE_ROOT))
from pipeline_paths import activate as _layout_activate, resource_path as _layout_resource
_layout_activate()


import argparse
import csv
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT=_LAYOUT_BUNDLE_ROOT;sys.path.insert(0,str(ROOT))
import report_audits as reporter
from schisto_adhesome import datasets as data
from schisto_adhesome.reports import read
import adhesome_workflow as core


class AuditReportTests(unittest.TestCase):
    def setUp(self):
        self.expected=patch.dict(data.EXPECTED,calibration=dict(positive=1,hard_negative=1,background_negative=1),
            validation=dict(positive=3,negative=1));self.expected.start();self.addCleanup(self.expected.stop)

    def fixture(self,base):
        from openpyxl import Workbook
        run=_layout_resource(base/'original');run.mkdir();contract={'datasets':{}};audits={}
        for kind in ('calibration','validation'):
            branch=run/kind;branch.mkdir()
            definitions=([('candidate|P1|ONE','Positive control','A'),('hard_negative|N1|NEG','Hard negative control','C'),
                          ('background_negative|N2|NEG','Background negative control','D')] if kind=='calibration' else
                         [('mus_musculus|G1|ONE','Positive control','A'),('mus_musculus|G2|MISS','Positive control','E'),
                          ('mus_musculus|G3|WRONG','Positive control','F'),('negative|N3|NEG','Negative control','D')])
            fasta=base/(kind+'.faa');labels=base/(kind+('.xlsx' if kind=='calibration' else '.csv'))
            fasta.write_text(''.join('>'+pid+' taxid=9606\nM'+aa*19+'\n' for pid,cls,aa in definitions))
            label_rows=[]
            for index,(pid,cls,aa) in enumerate(definitions):
                label=dict(**{'class':cls},genes='TEST'+str(index),trusted_family='test_family' if cls=='Positive control' else '')
                label.update(dict(accession=pid.split('|')[1],protein_name='Synthetic '+str(index),length=20) if kind=='calibration' else
                    dict(fasta_id=pid,mapped_accession='REF'+str(index),Protein_name='Synthetic '+str(index),sequence_length=20))
                label_rows.append(label)
            if kind=='calibration':
                book=Workbook();sheet=book.active;sheet.title='control_identity';fields=list(label_rows[0]);sheet.append(fields)
                for row in label_rows:sheet.append([row[k] for k in fields])
                book.save(labels)
            else:
                with _layout_resource(labels).open('w',newline='') as handle:
                    writer=csv.DictWriter(handle,list(label_rows[0]));writer.writeheader();writer.writerows(label_rows)
            trusted,summary=data.load_dataset(kind,fasta,labels,select_labelled=True);contract['datasets'][kind]=summary
            manifest=data.materialize(kind,trusted,_layout_resource(branch/'inputs'),{})
            handoff=_layout_resource(branch/'handoff');handoff.mkdir();included=[r for r in manifest if '|G2|' not in r['original_id']]
            data.table(_layout_resource(handoff/'hypotheses.tsv'),[dict(r,family='test_family') for r in included]);core.seal(handoff,'handoff')
            annotation=_layout_resource(branch/'annotation');annotation.mkdir()
            core.seal(annotation,'annotation',handoff=str(handoff),handoff_manifest=data.sha(_layout_resource(handoff/'workflow_manifest.json')))
            data.table(_layout_resource(branch/'reference_overlap.tsv'),[])
            audit=_layout_resource(base/'audits')/kind;audit.mkdir(parents=True);audits[kind]=audit;calls=[]
            for row in included:
                accepted=row['control_class']=='positive';family='wrong_family' if '|G3|' in row['original_id'] else 'test_family' if accepted else 'unassigned'
                calls.append(dict(protein_id=row['protein_id'],Specie=row['Specie'],family='test_family',preaudit_family='test_family',
                    audit_family=family,assigned_family=family,putative_audited_family=family,grade='A' if accepted else 'D',
                    domain_score_0_4=4 if accepted else 0,assignment_status='assigned' if accepted else 'unassigned',
                    family_decision='family_supported' if accepted else 'architecture_gate_failed',family_confidence='strong',
                    orthology_score_0_3=3,motif_score_0_2=0,localization_score_0_1=0,total_score_0_10=7 if accepted else 3))
            book=Workbook();sheet=book.active;sheet.title='Screening_hypothesis';fields=list(calls[0]);sheet.append(fields)
            for row in calls:sheet.append([row[k] for k in fields])
            book.save(_layout_resource(audit/'Stage_2_output.xlsx'));data.table(_layout_resource(audit/'audited_family_assignments.tsv'),calls)
            data.table(_layout_resource(audit/'negative_control_exclusions.tsv'),[])
            core.seal(audit,'audit',annotation=str(annotation),annotation_manifest=data.sha(_layout_resource(annotation/'workflow_manifest.json')),stage2_profile='calibrated_v1')
        data.dump(_layout_resource(run/'run_contract.json'),contract)
        resources=_layout_resource(base/'resources');folder=_layout_resource(resources/'detection/references/adhesome_v1.0');folder.mkdir(parents=True)
        for filename in ('search_catalogue.csv','metazoan_catalogue.csv'):(folder/filename).write_text('sequence_sha256\n')
        folder=_layout_resource(resources/'validation');folder.mkdir();data.dump(_layout_resource(folder/'identity_probability_model.json'),dict(samples=[]))
        for name in ('scripts/calibrate_identity.py','scripts/label_ontology.py','schisto_adhesome/reports.py'):
            path=resources/name;path.parent.mkdir(exist_ok=True);path.write_text('# fixture\n')
        args=argparse.Namespace(run=run,out=_layout_resource(base/'reports'),calibration_audit=audits['calibration'],validation_audit=audits['validation'],
                                calibration_labels=None,validation_labels=None)
        return args,resources

    def test_reports_new_audits_keeps_stage1_misses_and_family_errors(self):
        with tempfile.TemporaryDirectory() as directory:
            args,resources=self.fixture(_layout_resource(Path(directory)))
            before={p:_layout_resource(p).read_bytes() for p in _layout_resource(Path(directory)).rglob('*') if _layout_resource(p).is_file()}
            with patch.object(reporter,'ROOT',resources),patch.object(core,'execute',side_effect=AssertionError('No inference')):
                result=reporter.generate(args)
            self.assertTrue(_layout_resource(Path(result['validation_report'])).is_file())
            cal=json.loads((_layout_resource(args.out/'calibration/summary.json')).read_text());val=json.loads((_layout_resource(args.out/'validation/summary.json')).read_text())
            self.assertEqual((cal['supported']['TP'],cal['supported']['TN']),(1,2))
            self.assertEqual((val['supported']['TP'],val['supported']['FP'],val['supported']['FN']),(2,0,1))
            self.assertEqual((val['family']['TP'],val['family']['FP'],val['family']['FN']),(1,1,2))
            self.assertTrue((_layout_resource(args.out/'calibration/calibration_review.html')).is_file())
            provenance=json.loads((_layout_resource(args.out/'validation/validation_review_provenance.json')).read_text())
            self.assertEqual(_layout_resource(Path(provenance['source_files']['audit']['path'])).parent,args.validation_audit)
            self.assertEqual(_layout_resource(Path(provenance['source_files']['calibration_predictions']['path'])).parent,_layout_resource(args.out/'calibration'))
            for path,content in before.items():self.assertEqual(_layout_resource(path).read_bytes(),content)
            self.assertFalse(json.loads((_layout_resource(args.out/'reporting_provenance.json')).read_text())['analytical_pipeline_rerun'])

    def test_changed_labels_fail_before_output_is_created(self):
        with tempfile.TemporaryDirectory() as directory:
            args,resources=self.fixture(_layout_resource(Path(directory)))
            labels=_layout_resource(Path(json.loads((args.run/'run_contract.json').read_text())['datasets']['validation']['labels_file']))
            labels.write_text(_layout_resource(labels).read_text()+'\n')
            with self.assertRaisesRegex(ValueError,'reference labels differ'):reporter.generate(args)
            self.assertFalse(_layout_resource(args.out).exists())

    def test_equivalent_calibration_workbook_needs_explicit_option(self):
        from openpyxl import load_workbook
        with tempfile.TemporaryDirectory() as directory:
            args,resources=self.fixture(_layout_resource(Path(directory)));contract=json.loads((_layout_resource(args.run/'run_contract.json')).read_text())
            labels=_layout_resource(Path(contract['datasets']['calibration']['labels_file']))
            book=load_workbook(labels);book.properties.title='Workbook packaging updated';book.save(labels);book.close()
            with self.assertRaisesRegex(ValueError,'reference labels differ'):
                reporter.checked_cohort(args.run,'calibration',contract)
            rows,*_=reporter.checked_cohort(args.run,'calibration',contract,allow_equivalent_calibration_label_bytes=True)
            self.assertEqual(len(rows),3);self.assertFalse(_layout_resource(args.out).exists())

    def test_equivalence_option_cannot_accept_changed_calibration_truth(self):
        from openpyxl import load_workbook
        with tempfile.TemporaryDirectory() as directory:
            args,resources=self.fixture(_layout_resource(Path(directory)));contract=json.loads((_layout_resource(args.run/'run_contract.json')).read_text())
            labels=_layout_resource(Path(contract['datasets']['calibration']['labels_file']))
            book=load_workbook(labels);sheet=book['control_identity'];headings=[c.value for c in sheet[1]]
            sheet.cell(2,headings.index('trusted_family')+1,'altered_truth');book.save(labels);book.close()
            with self.assertRaisesRegex(ValueError,'metadata differ'):
                reporter.checked_cohort(args.run,'calibration',contract,allow_equivalent_calibration_label_bytes=True)
            self.assertFalse(_layout_resource(args.out).exists())

    def test_calibration_equivalence_option_never_relaxes_validation_hash(self):
        with tempfile.TemporaryDirectory() as directory:
            args,resources=self.fixture(_layout_resource(Path(directory)));contract=json.loads((_layout_resource(args.run/'run_contract.json')).read_text())
            labels=_layout_resource(Path(contract['datasets']['validation']['labels_file']));labels.write_text(_layout_resource(labels).read_text()+'\n')
            with self.assertRaisesRegex(ValueError,'reference labels differ'):
                reporter.checked_cohort(args.run,'validation',contract,allow_equivalent_calibration_label_bytes=True)
            self.assertFalse(_layout_resource(args.out).exists())

    def test_wrong_branch_audit_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            args,resources=self.fixture(_layout_resource(Path(directory)));args.validation_audit=args.calibration_audit
            with self.assertRaisesRegex(ValueError,'different annotation'):reporter.generate(args)
            self.assertFalse(_layout_resource(args.out).exists())

    def test_changed_reference_sequence_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            args,resources=self.fixture(_layout_resource(Path(directory)))
            path=_layout_resource(Path(json.loads((args.run/'run_contract.json').read_text())['datasets']['calibration']['protein_file']))
            path.write_text(_layout_resource(path).read_text().replace('MAA','MVA'))
            with self.assertRaisesRegex(ValueError,'reference FASTA differs'):reporter.generate(args)

    def test_changed_manifest_and_audit_seal_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            args,resources=self.fixture(_layout_resource(Path(directory)));path=_layout_resource(args.run/'calibration/inputs/sequence_manifest.tsv')
            text=_layout_resource(path).read_text();path.write_text(text.replace('Synthetic 0','Untrusted name'))
            with self.assertRaisesRegex(ValueError,'metadata differ'):reporter.generate(args)
            path.write_text(text)
            (_layout_resource(args.calibration_audit/'audited_family_assignments.tsv')).write_text('changed')
            with self.assertRaisesRegex(ValueError,'Changed workflow input'):reporter.generate(args)

    def test_existing_and_nested_outputs_are_refused(self):
        with tempfile.TemporaryDirectory() as directory:
            args,resources=self.fixture(_layout_resource(Path(directory)))
            for path in (_layout_resource(args.run/'reports'),_layout_resource(args.calibration_audit/'reports'),args.validation_audit,_layout_resource(Path(directory))):
                with self.assertRaises(ValueError):reporter.check_output(path,args.run,[args.calibration_audit,args.validation_audit])
            args.out.mkdir()
            with self.assertRaisesRegex(ValueError,'new report directory'):reporter.generate(args)


if __name__=='__main__':unittest.main()
