"""Synthetic Stage 2 fixtures; no stored audit or biological search is executed."""
import contextlib
import csv
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import schisto_pfam_cdd_reclassify as review


class PosthocTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='posthoc_synthetic_')
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.source = self.root/'Stage_2.tsv'
        self.out = self.root/'review'

    def row(self, protein='synthetic', species='SyntheticA', support='unresolved', **extra):
        return dict(protein_id=protein,Specie=species,family_support_level=support,
            putative_audited_family='unassigned',protein_name='Unresolved protein',
            grade='D',Pfam_architecture='PF00307:1-20',
            NCBI_CDD_architecture='CH_PARV_rpt1:1-20',length='20',**extra)

    def write(self, rows, fields=None):
        fields = fields or list(dict.fromkeys(k for row in rows for k in row))
        review.write_tsv(self.source,rows,fields)

    def invoke(self,*extra):
        with contextlib.redirect_stdout(io.StringIO()),contextlib.redirect_stderr(io.StringIO()):
            return review.main(['--input',str(self.source),'--outdir',str(self.out),*extra])

    def table(self,name):
        with (self.out/name).open(encoding='utf-8',newline='') as f:
            reader=csv.DictReader(f,delimiter='\t')
            return list(reader),reader.fieldnames

    def test_current_schema_preserves_supported_and_all_original_fields(self):
        rows=[self.row('resolved',support='specific_supported'),self.row('broad',support='broad_supported'),
              self.row('pending',support='provisional'),self.row('ambiguous',support='ambiguous')]
        self.write(rows);before=self.source.read_bytes()
        self.assertEqual(self.invoke(),0)
        reviewed,_=self.table('Schisto_unassigned_Pfam_CDD_classification.tsv')
        self.assertEqual({r['Protein ID'] for r in reviewed},{'pending','ambiguous'})
        self.assertTrue(all(r['Suggested family / class']=='Parvin-like' for r in reviewed))
        reconciled,_=self.table('Stage_2_posthoc_refinement.tsv')
        for original,result in zip(rows,reconciled):
            self.assertTrue(all(result[k]==v for k,v in original.items()))
        self.assertEqual(reconciled[0]['posthoc_proposed_family_or_class'],'')
        self.assertEqual(self.source.read_bytes(),before)

    def test_optional_evidence_columns_missing_or_blank_are_supported(self):
        row=self.row();del row['NCBI_CDD_architecture'];del row['length']
        self.write([row]);self.assertEqual(self.invoke(),0)
        result,_=self.table('Schisto_unassigned_Pfam_CDD_classification.tsv')
        self.assertTrue(result[0]['Evidence tier'].startswith('C'))
        self.assertEqual(result[0]['DeepTMHMM'],'')
        self.assertEqual(result[0]['Protein length (aa)'],'')

    def test_pfam_only_versioned_accession_and_cdd_model_without_coordinates(self):
        row=self.row();row['Pfam_architecture']='PF00307.31';row['NCBI_CDD_architecture']='CH_PARV_rpt1'
        self.write([row]);self.assertEqual(self.invoke(),0)
        result,_=self.table('Schisto_unassigned_Pfam_CDD_classification.tsv')
        self.assertEqual(result[0]['Suggested family / class'],'Parvin-like')

    def test_cross_species_ids_retained_exact_duplicates_accounted(self):
        first=self.row();other=self.row(species='SyntheticB')
        self.write([first,first,other]);self.assertEqual(self.invoke(),0)
        results,_=self.table('Stage_2_posthoc_refinement.tsv')
        self.assertEqual(len(results),2)
        manifest=json.loads((self.out/'run_manifest.json').read_text(encoding='utf-8'))
        self.assertEqual(manifest['duplicate_rows'],1)
        self.assertEqual(manifest['unique_species_protein_pairs'],2)
        self.assertTrue(manifest['source_record_accounting_passed'])

    def test_conflicting_duplicate_fails_before_writing(self):
        first=self.row();second=dict(first,grade='C')
        self.write([first,second]);self.assertEqual(self.invoke(),1)
        self.assertFalse(self.out.exists())

    def test_legacy_space_labels_are_supported_preserved(self):
        self.write([self.row(support='Specific Supported'),self.row('two',support='broad supported')])
        self.assertEqual(self.invoke(),0)
        result,fields=self.table('Schisto_unassigned_Pfam_CDD_classification.tsv')
        self.assertEqual(result,[]);self.assertIn('Evidence tier',fields)
        priorities,fields=self.table('priority_review.tsv')
        self.assertEqual(priorities,[]);self.assertIn('Protein ID',fields)

    def test_negative_exclusions_cannot_be_reintroduced(self):
        self.write([self.row('negative',stage2_excluded='true'),
                    self.row('background',control_class='background_negative'),self.row('retained')])
        self.assertEqual(self.invoke(),0)
        result,_=self.table('Schisto_unassigned_Pfam_CDD_classification.tsv')
        self.assertEqual([r['Protein ID'] for r in result],['retained'])
        all_rows,_=self.table('Stage_2_posthoc_refinement.tsv')
        self.assertEqual([r['posthoc_status'] for r in all_rows],['negative_excluded','negative_excluded','reviewed'])

    def test_orthogroup_context_does_not_upgrade_broad_domains(self):
        rows=[self.row(str(i),support='specific_supported',Orthology_Orthogroup='OG0000001') for i in range(3)]
        for r in rows:r['putative_audited_family']='Parvin-like'
        target=self.row('pending',Orthology_Orthogroup='OG0000001');target['NCBI_CDD_architecture']='CH:1-20'
        self.write([*rows,target]);self.assertEqual(self.invoke(),0)
        result,_=self.table('Schisto_unassigned_Pfam_CDD_classification.tsv')
        self.assertTrue(result[0]['Evidence tier'].startswith('C'))
        self.assertIn('OG suggests Parvin-like',result[0]['Suggested family / class'])
        self.assertEqual(result[0]['Original audited family'],'unassigned')

    def test_check_inputs_and_invalid_schema_write_nothing(self):
        self.write([self.row()]);self.assertEqual(self.invoke('--check-inputs'),0)
        self.assertFalse(self.out.exists())
        row=self.row();del row['family_support_level'];self.write([row])
        self.assertEqual(self.invoke(),1);self.assertFalse(self.out.exists())

    def test_existing_directory_and_sealed_audit_nesting_refused(self):
        self.write([self.row()]);self.out.mkdir();sentinel=self.out/'retained.txt';sentinel.write_text('retained')
        self.assertEqual(self.invoke(),1);self.assertEqual(sentinel.read_text(),'retained')
        self.out=self.root/'nested'
        (self.root/'workflow_manifest.json').write_text('{}')
        self.assertEqual(self.invoke(),1);self.assertFalse(self.out.exists())

    def test_unknown_support_and_duplicate_headers_rejected(self):
        self.write([self.row(support='mistyped_supported')]);self.assertEqual(self.invoke(),1)
        self.source.write_text('Specie\tprotein_id\tprotein_id\n',encoding='utf-8')
        self.assertEqual(self.invoke(),1);self.assertFalse(self.out.exists())

    def test_all_negative_empty_audit_has_complete_header_only_reports(self):
        self.write([],['Specie','protein_id','putative_audited_family','protein_name','assignment_confidence','grade'])
        self.assertEqual(self.invoke(),0)
        results,fields=self.table('Stage_2_posthoc_refinement.tsv')
        self.assertEqual(results,[]);self.assertIn('posthoc_status',fields)
        summary,_=self.table('tier_summary.tsv')
        self.assertEqual([r['Candidates'] for r in summary],['0']*4)

    def test_comparison_uses_species_and_rejects_duplicate_keys(self):
        self.write([self.row(),self.row(species='SyntheticB')]);self.assertEqual(self.invoke(),0)
        reviewed,_=self.table('Schisto_unassigned_Pfam_CDD_classification.tsv')
        prior=self.out/'Schisto_unassigned_Pfam_CDD_classification.tsv'
        self.assertTrue(review.compare_prior(reviewed,prior)['all_rows_match'])
        with self.assertRaisesRegex(ValueError,'Ambiguous duplicate'):
            review.compare_prior(reviewed+[reviewed[0]],prior)

    @unittest.skipUnless(importlib.util.find_spec('openpyxl'),'openpyxl fixture writer unavailable')
    def test_real_audit_workbook_export_contract(self):
        main=ROOT.parent/'integrin_adhesome_analysis'
        sys.path.insert(0,str(main))
        from pipeline_paths import activate
        activate()
        from audit_report import prepare_sheet
        from classification_contract import write_sheet
        from openpyxl import Workbook
        raw=self.row('fixture')
        raw.update(assigned_family='unassigned',family_confidence='unresolved',preaudit_family='parvin')
        records,fields,_=prepare_sheet('Audited_family_assignment',[raw])
        book=Workbook();book.remove(book.active)
        write_sheet(book,2,'Audited_family_assignment',records,fields=fields)
        self.source=self.root/'Stage_2_output.xlsx';book.save(self.source)
        source_hash=hashlib.sha256(self.source.read_bytes()).hexdigest()
        self.assertEqual(self.invoke(),0)
        result,_=self.table('Schisto_unassigned_Pfam_CDD_classification.tsv')
        self.assertEqual(result[0]['Suggested family / class'],'Parvin-like')
        self.assertEqual(hashlib.sha256(self.source.read_bytes()).hexdigest(),source_hash)

    def test_entrypoint_from_another_working_directory(self):
        self.write([self.row()])
        result=subprocess.run([sys.executable,'-B',str(ROOT/'schisto_pfam_cdd_reclassify.py'),
            '--input',str(self.source),'--check-inputs'],cwd=self.root,capture_output=True,text=True)
        self.assertEqual(result.returncode,0,result.stderr)
        self.assertFalse(json.loads(result.stdout)['files_written'])

    def test_main_workflow_optional_command_dispatch(self):
        self.write([self.row()])
        main=ROOT.parent/'integrin_adhesome_analysis/main/adhesome_workflow.py'
        result=subprocess.run([sys.executable,'-B',str(main),'refine-unresolved',
            '--input',str(self.source),'--check-inputs'],cwd=self.root,capture_output=True,text=True)
        self.assertEqual(result.returncode,0,result.stderr)
        self.assertFalse(json.loads(result.stdout)['files_written'])


if __name__=='__main__':
    unittest.main()
