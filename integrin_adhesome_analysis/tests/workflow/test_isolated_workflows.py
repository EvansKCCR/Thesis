
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
from concurrent.futures import ThreadPoolExecutor
import csv
import json
from pathlib import Path
import sys
import tempfile
import threading
import unittest
from unittest.mock import patch

sys.path.insert(0,str(_layout_resource(Path(__file__)).resolve().parents[2]))
from schisto_adhesome import datasets as d
from schisto_adhesome import workflows as w
from schisto_adhesome.reports import summarize_predictions,metrics
import adhesome_workflow as core
import run_pipeline
from schisto_adhesome import branch as b


class IsolationTests(unittest.TestCase):
    def test_mode_is_required(self):
        with self.assertRaises(SystemExit):run_pipeline.main(['--out','unused'])

    def test_every_mode_schedules_calibration_and_selected_gate_in_parallel(self):
        for selected in ('target','validation'):
            barrier=threading.Barrier(2)
            def worker(kind):
                barrier.wait(timeout=3)
                return dict(status='complete',kind=kind,thread=threading.get_ident())
            result=w.launch_branches(['calibration',selected],worker,2)
            self.assertEqual(set(result),{'calibration',selected})
            self.assertTrue(all(v['status']=='complete' for v in result.values()))
            self.assertEqual(len({v['thread'] for v in result.values()}),2)

    def test_failure_does_not_cancel_self_test_and_serial_is_independent(self):
        for workers in (1,2):
            seen=[]
            def worker(kind):
                seen.append(kind)
                if kind=='target':raise ValueError('fixture failure')
                return dict(status='complete')
            result=w.launch_branches(['calibration','target'],worker,workers)
            self.assertEqual(set(seen),{'calibration','target'})
            self.assertEqual(result['calibration']['status'],'complete')
            self.assertEqual(result['target']['status'],'failed')

    def test_namespace_and_headers_do_not_expose_gold_labels(self):
        with tempfile.TemporaryDirectory() as td:
            root=_layout_resource(Path(td))
            rows=[dict(original_id='a|b',species='Homo sapiens',control_class='positive',
                sequence='MAAA',sequence_sha256=d.sequence_hash('MAAA'),length=4)]
            a=d.materialize('calibration',rows,_layout_resource(root/'cal'),{})
            b=d.materialize('validation',rows,_layout_resource(root/'val'),{})
            self.assertNotEqual(a[0]['protein_id'],b[0]['protein_id'])
            for p in _layout_resource(root).rglob('*.faa'):
                text=_layout_resource(p).read_text();self.assertNotIn('positive',text);self.assertNotIn('a|b',text)
            self.assertEqual(set(d.fasta(next((_layout_resource(root/'cal')).glob('*.faa'))).values()).pop()[1],'MAAA')

    def validation_fixture(self,root):
        proteins=_layout_resource(root/'validation.faa');labels=_layout_resource(root/'validation.csv')
        proteins.write_text('>mus_musculus|G|P\nMAAA\n>N|X taxid=9606\nMKKK\n>extra\nMSSS\n')
        with _layout_resource(labels).open('w',newline='') as h:
            writer=csv.DictWriter(h,['class','fasta_id','sequence_length']);writer.writeheader()
            writer.writerows([dict(zip(writer.fieldnames,['Positive control','mus_musculus|G|P',4])),
                             dict(zip(writer.fieldnames,['Negative control','N|X',4]))])
        return proteins,labels

    def test_labelled_selection_is_explicit_and_preserves_source_and_counts(self):
        with tempfile.TemporaryDirectory() as td,patch.dict(d.EXPECTED,validation={'positive':1,'negative':1}):
            root=_layout_resource(Path(td));proteins,labels=self.validation_fixture(root);original=_layout_resource(proteins).read_bytes()
            with self.assertRaisesRegex(ValueError,'unlabelled'):d.load_dataset('validation',proteins,labels)
            rows,summary=d.load_dataset('validation',proteins,labels,True)
            self.assertEqual(len(rows),2);self.assertEqual(summary['excluded_unlabelled_ids'],['extra'])
            self.assertEqual(_layout_resource(proteins).read_bytes(),original)
            proteins.write_text(_layout_resource(proteins).read_text().replace('MAAA','MA'))
            with self.assertRaisesRegex(ValueError,'length differs'):d.load_dataset('validation',proteins,labels,True)

    def test_invalid_residues_duplicates_and_nucleotides_fail(self):
        with tempfile.TemporaryDirectory() as td:
            p=_layout_resource(Path(td)/'x.faa')
            for text in ('>a\nMAAA\n>a\nMAAA\n','>a\nMA-\n','>a\n'+'ACGT'*30+'\n'):
                p.write_text(text)
                with self.assertRaises(ValueError):d.fasta(p)

    def test_frozen_resources_cannot_drift(self):
        with tempfile.TemporaryDirectory() as td:
            root=_layout_resource(Path(td));(_layout_resource(root/'core.py')).write_text('baseline');raw=_layout_resource(root/'ref');raw.mkdir();(_layout_resource(raw/'A.faa')).write_text('>a\nMAAA\n')
            contract=dict(algorithm_files={'core.py':d.sha(_layout_resource(root/'core.py'))},references=d.reference_inventory(raw))
            with patch.object(w,'ROOT',root):
                w.verify_frozen(contract)
                (_layout_resource(raw/'A.faa')).write_text('>a\nMKKK\n')
                with self.assertRaisesRegex(ValueError,'reference'):w.verify_frozen(contract)

    def test_no_backup_policy_preserves_partial_stage(self):
        with tempfile.TemporaryDirectory() as td,patch.dict('os.environ',SCHISTO_NO_ARCHIVES='1'):
            p=_layout_resource(Path(td)/'inference');p.mkdir();(_layout_resource(p/'partial')).write_text('evidence')
            with self.assertRaisesRegex(ValueError,'no-backup'):core.archive_attempt(p)
            self.assertEqual((_layout_resource(p/'partial')).read_text(),'evidence');self.assertEqual(list(p.parent.iterdir()),[p])

    def test_metrics_keep_missed_positives_and_undefined_denominators(self):
        rows=[dict(control_class='positive',supported=False),dict(control_class='negative',supported=False)]
        result=metrics(rows,'supported')
        self.assertEqual(result['FN'],1);self.assertEqual(result['TN'],1)
        self.assertEqual(result['F1'],0);self.assertIsNone(result['PPV'])
        self.assertEqual(result['balanced_accuracy'],.5)

    def test_post_classification_report_preserves_distinct_tiers_and_family_errors(self):
        def row(identifier,positive,called,family,expected):
            return dict(original_id=identifier,control_class='positive' if positive else 'negative',
                species='test',length_bin='<100',calibration_overlap=False,reference_overlap=False,training_overlap=False,
                grade='A' if called else 'C',stage1_recruited=True,supported=called,provisional=not called,
                including_provisional=True,before_exclusion=called,expected_families=expected,predicted_family=family)
        rows=[row('p',True,True,'wrong','right'),row('n',False,False,'unassigned','')]
        with tempfile.TemporaryDirectory() as td:
            result=summarize_predictions(rows,_layout_resource(Path(td)),'validation')
            self.assertEqual(result['supported']['TP'],1)
            self.assertEqual(result['family']['FP'],1);self.assertEqual(result['family']['FN'],1)
            self.assertEqual(result['family']['micro_F1'],0)
            self.assertIn('not independent',result['interpretation'])

    def test_target_gate_rejects_control_input_before_creating_run(self):
        with tempfile.TemporaryDirectory() as td:
            root=_layout_resource(Path(td));inputs=_layout_resource(root/'input');inputs.mkdir()
            (_layout_resource(inputs/'validation.faa')).write_text('>ctrl gene=ctrl\n'+'MAAA'*30+'\n')
            config=_layout_resource(root/'targets.json');settings=json.loads((_layout_resource(d.ROOT/'configs/isolated_workflow.json')).read_text())
            d.dump(config,dict(qc=settings['qc'],proteomes=[dict(code='validation',species='Homo sapiens',
                assembly='fixture',path='validation.faa',gene_regex=r'gene=(\S+)')]))
            runtime=_layout_resource(root/'settings.json');d.dump(runtime,settings)
            args=argparse.Namespace(mode='target',config=runtime,calibration_fasta=_layout_resource(root/'cal.faa'),
                calibration_labels=_layout_resource(root/'cal.xlsx'),validation_labels=_layout_resource(root/'missing.csv'),
                targets_config=config,input_dir=inputs,reference_raw=_layout_resource(root/'ref'),out=_layout_resource(root/'out'),resume=False)
            with patch.object(w,'load_dataset',return_value=([dict(original_id='ctrl')],{})):
                with self.assertRaisesRegex(ValueError,'control input'):w.prepare(args)
            self.assertFalse(_layout_resource(args.out).exists())

    def test_complete_branch_reuses_sealed_stages_and_retains_atlas_output(self):
        with tempfile.TemporaryDirectory() as td,patch.dict('os.environ',{},clear=False):
            root=_layout_resource(Path(td));target=_layout_resource(root/'target');target.mkdir()
            (_layout_resource(target/'targets.json')).write_text('{"proteomes":[]}')
            d.dump(_layout_resource(root/'run_contract.json'),dict(frozen_settings=dict(family_coverage='expanded',orthogroup_context_weight=1.0),
                references=dict(root=str(_layout_resource(root/'ref')),version='frozen'),datasets=dict(target=dict(input_root=str(root)))))
            commands=[]
            def fake(command):
                commands.append(command)
                out=_layout_resource(Path(command[command.index('--out')+1]));out.mkdir()
                kind={'start':'handoff','resume':'annotation','continue':'completed_workflow'}[command[0]]
                if command[0]=='continue':
                    final=_layout_resource(out/'final');final.mkdir()
                    (_layout_resource(final/'Stage_2_output.xlsx')).write_text('synthetic fixture workbook')
                    (_layout_resource(final/'audited_family_assignments.tsv')).write_text('protein_id\tSpecie\tgrade\tassigned_family\nTarget__p\tTarget\tA\tprofilin\n')
                core.seal(out,kind)
                return 0
            args=argparse.Namespace(run=root,kind='target',threads=1,reuse_detection=None,context_dir=None,
                                    allow_missing_context=True,compare_report=None)
            with patch.object(b,'verify_frozen'),patch.object(core,'main',side_effect=fake):
                self.assertEqual(b.execute(args),0);self.assertEqual(len(commands),3)
                original=(_layout_resource(target/'downstream/final/Stage_2_output.xlsx')).read_bytes()
                self.assertEqual(b.execute(args),0);self.assertEqual(len(commands),3)
                self.assertEqual((_layout_resource(target/'downstream/final/Stage_2_output.xlsx')).read_bytes(),original)
            state=json.loads((_layout_resource(target/'stages.json')).read_text())
            self.assertTrue(all(s['execution']=='reused' for s in state.values()))
            self.assertTrue((_layout_resource(target/'reports/report.html')).exists())
            self.assertFalse(list(_layout_resource(root).rglob('*attempt*')))

    def test_branch_exports_handoff_without_claiming_completed_evaluation(self):
        with tempfile.TemporaryDirectory() as td,patch.dict('os.environ',{},clear=False):
            root=_layout_resource(Path(td));branch=_layout_resource(root/'target');branch.mkdir()
            (_layout_resource(branch/'targets.json')).write_text('{"proteomes":[]}')
            d.dump(_layout_resource(root/'run_contract.json'),dict(frozen_settings={},datasets=dict(target=dict(input_root=str(root)))))
            def fake(command):
                out=_layout_resource(branch/'handoff');out.mkdir();core.seal(out,'handoff');return 0
            args=argparse.Namespace(run=root,kind='target',threads=1,reuse_detection=None,context_dir=None,
                                    allow_missing_context=False,compare_report=None)
            with patch.object(b,'verify_frozen'),patch.object(core,'main',side_effect=fake):
                self.assertEqual(b.execute(args),2)
            self.assertEqual(json.loads((_layout_resource(branch/'branch_status.json')).read_text())['status'],'awaiting_external_evidence')
            self.assertFalse((_layout_resource(branch/'reports')).exists())


if __name__=='__main__':unittest.main()
