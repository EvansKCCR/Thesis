"""Standalone audit reports use saved calls and their own labelled populations."""
import argparse
import csv
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
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
        run=base/'original';run.mkdir();contract={'datasets':{}};audits={}
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
                with labels.open('w',newline='') as handle:
                    writer=csv.DictWriter(handle,list(label_rows[0]));writer.writeheader();writer.writerows(label_rows)
            trusted,summary=data.load_dataset(kind,fasta,labels,select_labelled=True);contract['datasets'][kind]=summary
            manifest=data.materialize(kind,trusted,branch/'inputs',{})
            handoff=branch/'handoff';handoff.mkdir();included=[r for r in manifest if '|G2|' not in r['original_id']]
            data.table(handoff/'hypotheses.tsv',[dict(r,family='test_family') for r in included]);core.seal(handoff,'handoff')
            annotation=branch/'annotation';annotation.mkdir()
            core.seal(annotation,'annotation',handoff=str(handoff),handoff_manifest=data.sha(handoff/'workflow_manifest.json'))
            data.table(branch/'reference_overlap.tsv',[])
            audit=base/'audits'/kind;audit.mkdir(parents=True);audits[kind]=audit;calls=[]
            for row in included:
                accepted=row['control_class']=='positive';family='wrong_family' if '|G3|' in row['original_id'] else 'test_family' if accepted else 'unassigned'
                calls.append(dict(protein_id=row['protein_id'],Specie=row['Specie'],family='test_family',preaudit_family='test_family',
                    audit_family=family,assigned_family=family,putative_audited_family=family,grade='A' if accepted else 'D',
                    domain_score_0_4=4 if accepted else 0,assignment_status='assigned' if accepted else 'unassigned',
                    family_decision='family_supported' if accepted else 'architecture_gate_failed',family_confidence='strong',
                    orthology_score_0_3=3,motif_score_0_2=0,localization_score_0_1=0,total_score_0_10=7 if accepted else 3))
            book=Workbook();sheet=book.active;sheet.title='Screening_hypothesis';fields=list(calls[0]);sheet.append(fields)
            for row in calls:sheet.append([row[k] for k in fields])
            book.save(audit/'Stage_2_output.xlsx');data.table(audit/'audited_family_assignments.tsv',calls)
            data.table(audit/'negative_control_exclusions.tsv',[])
            core.seal(audit,'audit',annotation=str(annotation),annotation_manifest=data.sha(annotation/'workflow_manifest.json'),stage2_profile='calibrated_v1')
        data.dump(run/'run_contract.json',contract)
        resources=base/'resources';folder=resources/'01_detection/references/adhesome_v1.0';folder.mkdir(parents=True)
        for filename in ('search_catalogue.csv','metazoan_catalogue.csv'):(folder/filename).write_text('sequence_sha256\n')
        folder=resources/'02_calibration_set';folder.mkdir();data.dump(folder/'identity_probability_model.json',dict(samples=[]))
        for name in ('calibrate_identity.py','label_ontology.py','schisto_adhesome/reports.py'):
            path=resources/name;path.parent.mkdir(exist_ok=True);path.write_text('# fixture\n')
        args=argparse.Namespace(run=run,out=base/'reports',calibration_audit=audits['calibration'],validation_audit=audits['validation'],
                                calibration_labels=None,validation_labels=None)
        return args,resources

    def test_reports_new_audits_keeps_stage1_misses_and_family_errors(self):
        with tempfile.TemporaryDirectory() as directory:
            args,resources=self.fixture(Path(directory))
            before={p:p.read_bytes() for p in Path(directory).rglob('*') if p.is_file()}
            with patch.object(reporter,'ROOT',resources),patch.object(core,'execute',side_effect=AssertionError('No inference')):
                result=reporter.generate(args)
            self.assertTrue(Path(result['validation_report']).is_file())
            cal=json.loads((args.out/'calibration/summary.json').read_text());val=json.loads((args.out/'validation/summary.json').read_text())
            self.assertEqual((cal['supported']['TP'],cal['supported']['TN']),(1,2))
            self.assertEqual((val['supported']['TP'],val['supported']['FP'],val['supported']['FN']),(2,0,1))
            self.assertEqual((val['family']['TP'],val['family']['FP'],val['family']['FN']),(1,1,2))
            self.assertTrue((args.out/'calibration/calibration_review.html').is_file())
            provenance=json.loads((args.out/'validation/validation_review_provenance.json').read_text())
            self.assertEqual(Path(provenance['source_files']['audit']['path']).parent,args.validation_audit)
            self.assertEqual(Path(provenance['source_files']['calibration_predictions']['path']).parent,args.out/'calibration')
            for path,content in before.items():self.assertEqual(path.read_bytes(),content)
            self.assertFalse(json.loads((args.out/'reporting_provenance.json').read_text())['analytical_pipeline_rerun'])

    def test_changed_labels_fail_before_output_is_created(self):
        with tempfile.TemporaryDirectory() as directory:
            args,resources=self.fixture(Path(directory))
            labels=Path(json.loads((args.run/'run_contract.json').read_text())['datasets']['validation']['labels_file'])
            labels.write_text(labels.read_text()+'\n')
            with self.assertRaisesRegex(ValueError,'reference labels differ'):reporter.generate(args)
            self.assertFalse(args.out.exists())

    def test_wrong_branch_audit_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            args,resources=self.fixture(Path(directory));args.validation_audit=args.calibration_audit
            with self.assertRaisesRegex(ValueError,'different annotation'):reporter.generate(args)
            self.assertFalse(args.out.exists())

    def test_changed_reference_sequence_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            args,resources=self.fixture(Path(directory))
            path=Path(json.loads((args.run/'run_contract.json').read_text())['datasets']['calibration']['protein_file'])
            path.write_text(path.read_text().replace('MAA','MVA'))
            with self.assertRaisesRegex(ValueError,'reference FASTA differs'):reporter.generate(args)

    def test_changed_manifest_and_audit_seal_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            args,resources=self.fixture(Path(directory));path=args.run/'calibration/inputs/sequence_manifest.tsv'
            text=path.read_text();path.write_text(text.replace('Synthetic 0','Untrusted name'))
            with self.assertRaisesRegex(ValueError,'metadata differ'):reporter.generate(args)
            path.write_text(text)
            (args.calibration_audit/'audited_family_assignments.tsv').write_text('changed')
            with self.assertRaisesRegex(ValueError,'Changed workflow input'):reporter.generate(args)

    def test_existing_and_nested_outputs_are_refused(self):
        with tempfile.TemporaryDirectory() as directory:
            args,resources=self.fixture(Path(directory))
            for path in (args.run/'reports',args.calibration_audit/'reports',args.validation_audit,Path(directory)):
                with self.assertRaises(ValueError):reporter.check_output(path,args.run,[args.calibration_audit,args.validation_audit])
            args.out.mkdir()
            with self.assertRaisesRegex(ValueError,'new report directory'):reporter.generate(args)


if __name__=='__main__':unittest.main()
