"""Reporting-only fixtures: no biological pipeline/external tools are launched."""

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

sys.path.insert(0,str(_layout_resource(Path(__file__)).resolve().parents[2]))
from schisto_adhesome import validation_review as review
from schisto_adhesome import workflows, branch
from schisto_adhesome.datasets import dump, sha, table
from schisto_adhesome.reports import summarize_predictions
import adhesome_workflow as core


class ValidationReviewTests(unittest.TestCase):
    def fixture(self,base,excluded=False):
        run=_layout_resource(base/'different_output_name');folder=_layout_resource(run/'validation/reports');folder.mkdir(parents=True)
        labels=_layout_resource(base/'custom_reference.csv')
        source=[dict(fasta_id=key,**{'class':label},sequence_length=20,mapped_accession=key,
                     genes=gene,Protein_name=name) for key,label,gene,name in
                [('one','Positive control','GENE1','Protein one'),
                 ('two','Positive control','GENE2','Protein two'),
                 ('negative','Negative control','OTHER','Negative protein')]]
        with _layout_resource(labels).open('w',newline='') as handle:
            writer=csv.DictWriter(handle,list(source[0]));writer.writeheader();writer.writerows(source)
        rows=[]
        for index,reference in enumerate(source):
            is_excluded=excluded and index==2
            rows.append(dict(original_id=reference['fasta_id'],protein_id='Val__'+reference['fasta_id'],
                species='Example' if index<2 else 'Negative species',control_class='positive' if index<2 else 'negative',
                length=20,accession=reference['mapped_accession'],genes=reference['genes'],protein_name=reference['Protein_name'],
                expected_families='right' if index==0 else 'expected_two' if index==1 else '',
                predicted_family='right' if index==0 else 'unassigned' if is_excluded else 'wrong',
                stage1_recruited=True,supported=index==0,provisional=index!=0 and not is_excluded,
                including_provisional=not is_excluded,before_exclusion=index==0 or is_excluded,
                grade='A' if index==0 else 'not_accepted' if is_excluded else 'C',
                exclusion_gate='background_negative' if is_excluded else '',length_bin='<100',
                sequence_sha256=f'{index:064x}',
                calibration_overlap=False,reference_overlap=index==0,training_overlap=False))
        summarize_predictions(rows,folder,'validation')
        final=_layout_resource(run/'validation/downstream/final');final.mkdir(parents=True)
        audited=[dict(protein_id=r['protein_id'],putative_audited_family=r['predicted_family'],grade=r['grade'],
                      domain_score_0_4=4 if r['supported'] else 0,assignment_confidence='strong')
                 for r in rows if not (excluded and r['control_class']=='negative')]
        table(_layout_resource(final/'audited_family_assignments.tsv'),audited)
        table(_layout_resource(final/'negative_control_exclusions.tsv'),
              [dict(protein_id='Val__negative',exclusion_gate='background_negative')] if excluded else [])
        dump(_layout_resource(run/'run_contract.json'),dict(datasets=dict(validation=dict(labels_file=str(labels),labels_sha256=sha(labels),
             fasta_records=4,excluded_unlabelled_ids=['unlabelled']))))
        return run,rows,labels

    def add_calibration(self,run,rows):
        folder=_layout_resource(run/'calibration/reports');folder.mkdir(parents=True)
        summarize_predictions(rows,folder,'calibration')
        dump(_layout_resource(run/'calibration/branch_status.json'),dict(status='complete'))

    def test_arbitrary_cohort_and_filename_separates_membership_and_family_errors(self):
        with tempfile.TemporaryDirectory() as td:
            run,rows,labels=self.fixture(_layout_resource(Path(td)))
            inputs={p:_layout_resource(p).read_bytes() for p in _layout_resource(run).rglob('*') if _layout_resource(p).is_file()}
            result=review.generate(run)
            p=json.loads((_layout_resource(run/'validation/reports/validation_review_provenance.json')).read_text())
            self.assertEqual(p['labelled_records'],3)
            self.assertEqual(p['supported']['TP'],1)
            self.assertEqual(p['family_by_policy']['including_provisional']['FP'],2)
            self.assertAlmostEqual(p['family_by_policy']['including_provisional']['F1'],.4)
            text=_layout_resource(Path(result['report'])).read_text(encoding='utf-8')
            self.assertIn('custom_reference.csv',text);self.assertIn('3 labelled proteins',text)
            self.assertEqual(result['calibration_comparison'],'pending_or_unavailable')
            for path,value in inputs.items():self.assertEqual(_layout_resource(path).read_bytes(),value)
            self.assertFalse(p['analytical_pipeline_rerun']);self.assertFalse(p['backups_created'])

    def test_excluded_candidates_are_accounted_without_reappearing_in_final_audit(self):
        with tempfile.TemporaryDirectory() as td:
            run,rows,labels=self.fixture(_layout_resource(Path(td)),excluded=True)
            review.generate(run)
            records=review.read(_layout_resource(run/'validation/reports/validation_provisional_family_review.tsv'))
            self.assertEqual(len(records),1);self.assertEqual(records[0]['original_id'],'two')
            self.assertNotIn('Val__negative',[r['protein_id'] for r in review.read(_layout_resource(run/'validation/downstream/final/audited_family_assignments.tsv'))])

    def test_changed_reference_is_rejected_before_overwriting_report(self):
        with tempfile.TemporaryDirectory() as td:
            run,rows,labels=self.fixture(_layout_resource(Path(td)));review.generate(run)
            report=_layout_resource(run/'validation/reports/validation_review.html');before=_layout_resource(report).read_bytes()
            labels.write_text(_layout_resource(labels).read_text()+'\n')
            with self.assertRaisesRegex(ValueError,'reference changed'):review.generate(run)
            self.assertEqual(_layout_resource(report).read_bytes(),before)

    def test_duplicate_or_missing_truth_ids_are_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            run,rows,labels=self.fixture(_layout_resource(Path(td)))
            raw=review.label_rows(labels)
            for invalid in (raw+[raw[0]],raw[:-1]):
                with self.assertRaisesRegex(ValueError,'unique, exact join'):review.verify_reference(rows,invalid)

    def test_frozen_cohort_counts_cannot_silently_drift(self):
        with tempfile.TemporaryDirectory() as td:
            run,rows,labels=self.fixture(_layout_resource(Path(td)))
            path=_layout_resource(run/'run_contract.json');contract=json.loads(_layout_resource(path).read_text())
            contract['datasets']['validation']['labelled_records']=99;dump(path,contract)
            with self.assertRaisesRegex(ValueError,'cohort counts'):review.generate(run)

    def test_rejected_benchmark_is_reported_without_blocking_report_or_faking_weight(self):
        with tempfile.TemporaryDirectory() as td:
            run,rows,labels=self.fixture(_layout_resource(Path(td)));bench=_layout_resource(run/'validation/downstream/benchmark');bench.mkdir()
            dump(_layout_resource(bench/'summary.json'),dict(status='quality_rejected',validated_weight=None,scope='all',protocol='open_orthobench',
                metrics=dict(precision=.5,recall=.8,f1=.615),quality_acceptance=dict(checks={
                    'precision':dict(minimum=.8,passed=False),'f1':dict(minimum=.7,passed=False)})))
            review.generate(run)
            text=(_layout_resource(run/'validation/reports/validation_review.html')).read_text(encoding='utf-8')
            self.assertIn('quality_rejected',text);self.assertIn('validated weight: None',text)
            self.assertIn('50.0%',text)

    def test_pending_parallel_comparison_refreshes_when_calibration_finishes(self):
        with tempfile.TemporaryDirectory() as td:
            run,rows,labels=self.fixture(_layout_resource(Path(td)))
            self.assertEqual(review.generate(run)['calibration_comparison'],'pending_or_unavailable')
            self.add_calibration(run,rows)
            self.assertEqual(review.generate(run)['calibration_comparison'],'available')
            text=(_layout_resource(run/'validation/reports/validation_review.html')).read_text(encoding='utf-8')
            self.assertIn('Calibration contains 3 controls',text)
            self.assertNotIn('comparison is pending',text)

    def test_zero_provisional_errors_exports_header_only(self):
        with tempfile.TemporaryDirectory() as td:
            run,rows,labels=self.fixture(_layout_resource(Path(td)))
            for r in rows[1:]:
                r.update(provisional=False,including_provisional=False,predicted_family='unassigned',grade='D')
            summarize_predictions(rows,_layout_resource(run/'validation/reports'),'validation')
            audit=review.read(_layout_resource(run/'validation/downstream/final/audited_family_assignments.tsv'))
            for r in audit[1:]:r.update(putative_audited_family='unassigned',grade='D')
            table(_layout_resource(run/'validation/downstream/final/audited_family_assignments.tsv'),audit)
            review.generate(run)
            path=_layout_resource(run/'validation/reports/validation_provisional_family_review.tsv')
            self.assertIn('original_id',_layout_resource(path).read_text());self.assertEqual(review.read(path),[])

    def test_completed_validation_branch_auto_reports_without_launching_any_analysis(self):
        with tempfile.TemporaryDirectory() as td:
            base=_layout_resource(Path(td));run,rows,labels=self.fixture(base)
            contract=json.loads((_layout_resource(run/'run_contract.json')).read_text())
            contract.update(frozen_settings=dict(family_coverage='expanded'),references=dict(root='unused',version='fixture'))
            dump(_layout_resource(run/'run_contract.json'),contract)
            for folder,kind in [('handoff','handoff'),('annotation','annotation'),('downstream','completed_workflow')]:
                p=_layout_resource(run/'validation')/folder;p.mkdir(exist_ok=True);core.seal(p,kind)
            for kind in ('calibration','validation'):
                p=_layout_resource(run/kind/'inputs');p.mkdir(parents=True,exist_ok=True);table(_layout_resource(p/'sequence_manifest.tsv'),rows)
            table(_layout_resource(run/'validation/reference_overlap.tsv'),[])
            resource_root=_layout_resource(base/'resources');catalogue=_layout_resource(resource_root/'detection/references/adhesome_v1.0');catalogue.mkdir(parents=True)
            for name in ('search_catalogue.csv','metazoan_catalogue.csv'):(catalogue/name).write_text('sequence_sha256\n')
            model=_layout_resource(resource_root/'validation');model.mkdir();dump(_layout_resource(model/'identity_probability_model.json'),dict(samples=[]))
            args=argparse.Namespace(kind='validation',run=run,threads=1,context_dir=None,allow_missing_context=False,
                                    reuse_detection=None,compare_report=None)
            with patch.object(branch,'ROOT',resource_root),patch.object(branch,'verify_frozen'),\
                 patch.object(branch,'evaluate_branch',return_value=rows),\
                 patch.object(core,'main',side_effect=AssertionError('Analysis must not launch')):
                self.assertEqual(branch.execute(args),0)
            status=json.loads((_layout_resource(run/'validation/branch_status.json')).read_text())
            self.assertTrue(status['report'].endswith('validation_review.html'))
            self.assertTrue(all(s['execution']=='reused' for s in status['stages'].values()))

    def test_scheduler_refreshes_after_both_workers_without_rerunning_branches(self):
        with tempfile.TemporaryDirectory() as td:
            run,rows,labels=self.fixture(_layout_resource(Path(td)))
            contract=dict(contract_sha256='fixture',frozen_settings=dict(orthogroup_context_weight=1.0),datasets={'calibration':{},'validation':{}})
            args=argparse.Namespace(mode='validation',out=run,threads=2,schedule='parallel',preflight_only=False)
            def completed(kinds,worker,workers):
                self.add_calibration(run,rows)
                return {kind:dict(status='complete') for kind in kinds}
            with patch.object(workflows,'prepare',return_value=contract),patch.object(workflows,'verify_frozen'),\
                 patch.object(workflows,'launch_branches',side_effect=completed),\
                 patch.object(workflows.subprocess,'run',side_effect=AssertionError('No analysis subprocess')):
                self.assertEqual(workflows.run(args),0)
            summary=json.loads((_layout_resource(run/'run_summary.json')).read_text())
            self.assertEqual(summary['gates']['validation']['validation_review']['calibration_comparison'],'available')
            self.assertIn('validation/reports/validation_review.html',(_layout_resource(run/'README.md')).read_text())

    def test_report_failure_is_distinguished_without_restarting_completed_analysis(self):
        with tempfile.TemporaryDirectory() as td:
            run,rows,labels=self.fixture(_layout_resource(Path(td)));labels.write_text(_layout_resource(labels).read_text()+'\n')
            contract=dict(contract_sha256='fixture',frozen_settings=dict(orthogroup_context_weight=1.0),datasets={'calibration':{},'validation':{}})
            args=argparse.Namespace(mode='validation',out=run,threads=2,schedule='parallel',preflight_only=False)
            gates={kind:dict(status='complete') for kind in ('calibration','validation')}
            with patch.object(workflows,'prepare',return_value=contract),patch.object(workflows,'verify_frozen'),\
                 patch.object(workflows,'launch_branches',return_value=gates),\
                 patch.object(workflows.subprocess,'run',side_effect=AssertionError('No analysis subprocess')):
                self.assertEqual(workflows.run(args),1)
            state=json.loads((_layout_resource(run/'validation/branch_status.json')).read_text())
            self.assertIn('Validation reporting failed',state['error'])
            self.assertTrue((_layout_resource(run/'validation/downstream/final/audited_family_assignments.tsv')).is_file())


if __name__=='__main__':unittest.main()
