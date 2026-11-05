"""Synthetic saved refinements and independent labels; never score a stored run."""
import contextlib
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
import subprocess

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import schisto_pfam_cdd_reclassify as refine
import report_posthoc as report


class ReportingTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory(prefix='posthoc_reporting_synthetic_')
        self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
        self.audit=self.root/'audit';self.audit.mkdir()
        self.refinement=self.root/'refinement';self.refinement.mkdir()
        self.labels=self.root/'labels.csv'
        self.mapping=self.root/'inputs'/'sequence_manifest.tsv';self.mapping.parent.mkdir()
        self.out=self.root/'reports'
        # Two correct proposals, one incorrect positive and a false positive negative.
        # Missing label P4 is a full-cohort FN; N2 has an explicit exclusion.
        self.gold=[dict(fasta_id='P0',**{'class':'Positive control'},trusted_family='Src'),
            dict(fasta_id='P1',**{'class':'Positive control'},trusted_family='parvin_like'),
            dict(fasta_id='P2',**{'class':'Positive control'},trusted_family='cofilin'),
            dict(fasta_id='P3',**{'class':'Positive control'},trusted_family='NF2_merlin_like'),
            dict(fasta_id='P4',**{'class':'Positive control'},trusted_family='talin'),
            dict(fasta_id='N1',**{'class':'Hard negative control'},trusted_family=''),
            dict(fasta_id='N2',**{'class':'Background negative control'},trusted_family='')]
        self.entries=[self.entry('P0','specific_supported','SRC','supported_preserved'),
            self.entry('P1'),self.entry('P2',proposal='ADF/cofilin-domain protein',tier='B'),
            self.entry('P3',proposal='Ezrin/radixin/moesin-like FERM protein'),
            self.entry('N1'),self.entry('N2',status='negative_excluded')]
        self.save()

    def entry(self,name,support='unresolved',family='unassigned',status='reviewed',
              proposal='Parvin-like',tier='A'):
        return dict(Specie='Synthetic',protein_id='internal_'+name,family_support_level=support,
            putative_audited_family=family,protein_name='Original unchanged name',grade='D',length='20',
            posthoc_status=status,posthoc_proposed_family_or_class=proposal if status=='reviewed' else '',
            posthoc_evidence_tier=tier+': synthetic' if status=='reviewed' else '')

    def save(self):
        with self.labels.open('w',encoding='utf-8',newline='') as f:
            import csv
            writer=csv.DictWriter(f,fieldnames=list(self.gold[0]));writer.writeheader();writer.writerows(self.gold)
        mappings=[dict(original_id=r['fasta_id'],Specie='Synthetic',protein_id='internal_'+r['fasta_id'],
            control_class=report.category(r['class']),length='20') for r in self.gold]
        refine.write_tsv(self.mapping,mappings)
        path=self.refinement/'Stage_2_posthoc_refinement.tsv'
        refine.write_tsv(path,self.entries)
        receipt=dict(status='complete',source_record_accounting_passed=True,
            unique_species_protein_pairs=len(self.entries),input_sha256='synthetic_source_hash',
            input=str(self.audit/'Stage_2_output.xlsx'),column_mapping={
                k:k for k in ('Specie','protein_id','family_support_level','putative_audited_family','length')},
            output_sha256={path.name:refine.sha256(path)})
        (self.refinement/'run_manifest.json').write_text(json.dumps(receipt),encoding='utf-8')

    def run_report(self,*extra):
        with contextlib.redirect_stdout(io.StringIO()),contextlib.redirect_stderr(io.StringIO()):
            return report.main(['--refinement',str(self.refinement),'--labels',str(self.labels),
                '--sequence-manifest',str(self.mapping),'--kind','validation','--out',str(self.out),*extra])

    def table(self,name):return report.read_table(self.out/name)[0]

    def test_separate_binary_policies_full_cohort_and_unresolved_denominators(self):
        before={p:refine.sha256(p) for p in [self.labels,self.mapping,*self.refinement.iterdir()]}
        self.assertEqual(self.run_report(),0)
        rows=self.table('binary_metrics.tsv')
        original=next(r for r in rows if r['policy']=='original_supported')
        self.assertEqual([original[k] for k in ('TP','FP','TN','FN')],['1','0','2','4'])
        proposed=next(r for r in rows if r['policy']=='posthoc_AB' and r['population']=='full_labelled_cohort')
        self.assertEqual([proposed[k] for k in ('TP','FP','TN','FN')],['3','1','1','2'])
        self.assertAlmostEqual(float(proposed['PPV']),.75)
        self.assertAlmostEqual(float(proposed['F1']),6/9)
        eligible=next(r for r in rows if r['policy']=='posthoc_AB' and r['population']=='posthoc_eligible')
        self.assertEqual(eligible['n'],'4');self.assertEqual(eligible['recall'],'1.0')
        self.assertTrue(all('combined' not in r['policy'] for r in rows))
        self.assertEqual(before,{p:refine.sha256(p) for p in before})

    def test_wrong_family_counts_both_fp_and_fn_merlin_not_erm(self):
        self.assertEqual(self.run_report(),0)
        family=next(r for r in self.table('family_metrics.tsv') if r['policy']=='posthoc_AB' and r['population']=='full_labelled_cohort')
        self.assertEqual([family[k] for k in ('TP','FP','FN')],['2','2','3'])
        predictions=self.table('predictions.tsv')
        merlin=next(r for r in predictions if r['label_id']=='P3')
        self.assertEqual(merlin['posthoc_AB_family_correct'],'False')
        self.assertEqual(merlin['expected_families'],'NF2_merlin_like')
        self.assertEqual(next(r for r in predictions if r['label_id']=='P0')['original_supported_family_correct'],'True')

    def test_c_and_d_not_promoted_even_with_supported_orthogroup_hint(self):
        self.entries[1]['posthoc_proposed_family_or_class']='FN3 domain-containing protein [OG suggests fibronectin_like]'
        self.entries[1]['posthoc_evidence_tier']='C: domain superfamily';self.save()
        self.assertEqual(self.run_report(),0)
        row=next(r for r in self.table('predictions.tsv') if r['label_id']=='P1')
        self.assertEqual(row['posthoc_AB'],'False');self.assertEqual(row['posthoc_evaluation_family'],'')

    def test_unmapped_truth_explicitly_excluded_from_family_only(self):
        self.gold[1]['trusted_family']='';self.save()
        self.assertEqual(self.run_report(),0)
        unmapped=self.table('unmapped_labels.tsv')
        self.assertEqual([r['label_id'] for r in unmapped],['P1'])
        family=next(r for r in self.table('family_metrics.tsv') if r['policy']=='posthoc_AB' and r['population']=='full_labelled_cohort')
        self.assertEqual(family['unmapped_positive_labels'],'1')
        binary=next(r for r in self.table('binary_metrics.tsv') if r['policy']=='posthoc_AB' and r['population']=='full_labelled_cohort')
        self.assertEqual(binary['TP'],'3')

    def test_unmapped_proposal_is_error_not_removed_from_precision(self):
        self.entries[1]['posthoc_proposed_family_or_class']='Unmapped synthetic identity';self.save()
        self.assertEqual(self.run_report(),0)
        self.assertEqual([r['label_id'] for r in self.table('unmapped_proposals.tsv')],['P1'])
        family=next(r for r in self.table('family_metrics.tsv') if r['policy']=='posthoc_AB' and r['population']=='full_labelled_cohort')
        self.assertEqual([family[k] for k in ('TP','FP','FN')],['1','3','4'])

    def test_mutated_output_and_failed_run_refused_before_reports(self):
        path=self.refinement/'Stage_2_posthoc_refinement.tsv'
        path.write_text(path.read_text(encoding='utf-8')+'\n',encoding='utf-8')
        self.assertEqual(self.run_report(),1);self.assertFalse(self.out.exists())
        self.save();path=self.refinement/'run_manifest.json'
        data=json.loads(path.read_text());data['status']='validation_failed';path.write_text(json.dumps(data))
        self.assertEqual(self.run_report(),1);self.assertFalse(self.out.exists())

    def test_duplicate_truth_or_incomplete_mapping_refused(self):
        self.gold.append(dict(self.gold[0]));self.save()
        self.assertEqual(self.run_report(),1);self.assertFalse(self.out.exists())
        self.gold.pop();self.save()
        rows,fields=report.read_table(self.mapping);refine.write_tsv(self.mapping,rows[:-1],fields)
        self.assertEqual(self.run_report(),1);self.assertFalse(self.out.exists())

    def test_zero_denominators_are_not_perfect_scores(self):
        binary,family=report.metrics([], 'posthoc_AB')
        self.assertIsNone(binary['PPV']);self.assertIsNone(binary['F1'])
        self.assertIsNone(family['recall'])

    def test_no_matching_panel_and_audit_length_mismatch_refused(self):
        self.entries=[dict(r,protein_id='other_'+r['protein_id']) for r in self.entries];self.save()
        self.assertEqual(self.run_report(),1);self.assertFalse(self.out.exists())
        self.entries=[self.entry('P1')];self.entries[0]['length']='21';self.save()
        self.gold[1]['sequence_length']='20'
        for r in self.gold:r.setdefault('sequence_length','20')
        self.save()
        self.assertEqual(self.run_report(),1);self.assertFalse(self.out.exists())

    def test_src_alias_does_not_include_scr(self):
        config=json.loads((ROOT/'config/evaluation_labels.json').read_text(encoding='utf-8'))
        self.assertEqual(report.normalized('src',config),'Src')
        self.assertEqual(report.normalized('SRC',config),'Src')
        self.assertEqual(report.normalized('SCR',config),'SCR')

    def test_main_reporting_command_from_another_cwd(self):
        script=ROOT.parent/'integrin_adhesome_analysis/main/adhesome_workflow.py'
        result=subprocess.run([sys.executable,'-B',str(script),'report-refinement',
            '--refinement',str(self.refinement),'--labels',str(self.labels),
            '--sequence-manifest',str(self.mapping),'--kind','validation','--out',str(self.out)],
            cwd=self.root,capture_output=True,text=True)
        self.assertEqual(result.returncode,0,result.stderr)
        self.assertTrue((self.out/'performance_report.html').is_file())

    def test_html_escapes_source_labels_and_no_classification_called(self):
        self.entries[1]['posthoc_proposed_family_or_class']='<script>bad</script>';self.save()
        with patch.object(refine,'classify',side_effect=AssertionError('Classification rerun')):
            self.assertEqual(self.run_report(),0)
        html=(self.out/'performance_report.html').read_text(encoding='utf-8')
        self.assertIn('Original supported audit calls and post hoc proposals are scored separately',html)
        self.assertNotIn('<script>',html)
        manifest=json.loads((self.out/'report_manifest.json').read_text())
        self.assertFalse(manifest['classification_rerun']);self.assertFalse(manifest['combined_policy_generated'])

    def test_existing_report_and_nested_source_output_refused(self):
        self.out.mkdir();(self.out/'keep').write_text('keep')
        self.assertEqual(self.run_report(),1);self.assertEqual((self.out/'keep').read_text(),'keep')
        self.out=self.refinement/'reports'
        self.assertEqual(self.run_report(),1);self.assertFalse(self.out.exists())

    def test_calibration_accession_join_uses_explicit_manifest(self):
        from openpyxl import Workbook
        book=Workbook();sheet=book.active;sheet.title='control_identity'
        sheet.append(['class','accession','protein_name','genes','length'])
        sheet.append(['Positive control','P00001','Alpha-parvin','PARVA',20])
        self.labels=self.root/'calibration_set.xlsx';book.save(self.labels)
        self.entries=[self.entry('P00001')];self.save_calibration_mapping()
        with contextlib.redirect_stdout(io.StringIO()),contextlib.redirect_stderr(io.StringIO()):
            code=report.main(['--refinement',str(self.refinement),'--labels',str(self.labels),
                '--sequence-manifest',str(self.mapping),'--kind','calibration','--out',str(self.out)])
        self.assertEqual(code,0)
        self.assertEqual(self.table('predictions.tsv')[0]['posthoc_AB_family_correct'],'True')

    def save_calibration_mapping(self):
        refine.write_tsv(self.mapping,[dict(original_id='candidate|P00001|PARVA',accession='P00001',
            Specie='Synthetic',protein_id='internal_P00001',control_class='positive',length='20')])
        path=self.refinement/'Stage_2_posthoc_refinement.tsv';refine.write_tsv(path,self.entries)
        receipt=dict(status='complete',source_record_accounting_passed=True,unique_species_protein_pairs=1,
            input_sha256='fixture',input=str(self.audit/'Stage_2_output.xlsx'),column_mapping={
                k:k for k in ('Specie','protein_id','family_support_level','putative_audited_family','length')},
            output_sha256={path.name:refine.sha256(path)})
        (self.refinement/'run_manifest.json').write_text(json.dumps(receipt),encoding='utf-8')


if __name__=='__main__':unittest.main()
