"""Synthetic tests of full broad gates, specific uncertainty and trace export."""

# Resolve the organized bundle independently of the caller working directory.
import sys as _bundle_sys
from pathlib import Path as _BundlePath
_BUNDLE_ROOT=next(p for p in _BundlePath(__file__).resolve().parents if (p/'pipeline_paths.py').is_file())
if str(_BUNDLE_ROOT) not in _bundle_sys.path:_bundle_sys.path.insert(0,str(_BUNDLE_ROOT))
from pipeline_paths import activate as _activate_bundle_paths
_activate_bundle_paths()

import sys as _sys
from pathlib import Path as _Path
_PIPELINE_ROOT=_Path(__file__).resolve().parents[2]
_sys.path.insert(0,str(_PIPELINE_ROOT))
_sys.path.insert(0,str(_PIPELINE_ROOT/"annotation"))
import csv
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from integrated_adhesome_pipeline import audit_record,run_xlsx

ROOT=Path(__file__).resolve().parents[2]


class Stage2ResolutionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.rules=json.loads((ROOT/'detection/config/family_rules.json').read_text());cls.rules['_stage2_profile']='calibrated_v1'

    def row(self,**values):
        row=dict(protein_id='Novel__g',Specie='Novel',family='actin',length=400,
                 architecture='PF00022:10-350',NCBI_CDD_architecture='',Orthology_HOG_reference_species_count=2)
        row.update(values);return row

    def test_broad_actin_support_does_not_claim_a_specific_arp(self):
        out=audit_record(self.row(NCBI_CDD_architecture='actin_Arp4_Arp6:10-350'),self.rules)
        self.assertEqual(out['assigned_family'],'actin');self.assertEqual(out['grade'],'A')
        self.assertEqual(out['family_support_level'],'broad_supported')
        self.assertEqual(out['family_confidence'],'strong');self.assertEqual(out['annotation_resolution'],'broad_family')
        self.assertEqual(out['unresolved_specific_families'],'arp4_like;arp6_like')
        self.assertNotIn('Arp4',out['protein_name'])

    def test_broad_pet_lim_has_its_own_full_observed_ordered_core(self):
        row=self.row(family='LIM_scaffold_screen',length=500,
            architecture='PF06297:10-100;PF00412:160-200;PF00412:250-300;PF00412:350-400',
            NCBI_CDD_architecture='PET:10-100')
        out=audit_record(row,self.rules)
        self.assertEqual(out['assigned_family'],'PET_LIM_like');self.assertEqual(out['family_support_level'],'broad_supported')
        self.assertNotIn('Testin',out['protein_name'])
        row['architecture']='PF06297:10-100;PF00412:160-200;PF00412:250-300'
        self.assertNotEqual(audit_record(row,self.rules)['family_support_level'],'broad_supported')
        row['architecture']='PF00412:10-50;PF00412:60-100;PF00412:120-160'
        self.assertEqual(audit_record(row,self.rules)['assigned_family'],'unassigned')

    def test_specific_complete_family_and_missing_orthology_are_separate_decisions(self):
        row=self.row(family='generic_kinase_screen',length=500,architecture='PF08953:10-60;PF00400:100-150',
            NCBI_CDD_architecture='DUF1899:10-60',Orthology_HOG_reference_species_count=0)
        out=audit_record(row,self.rules)
        self.assertEqual(out['assigned_family'],'coronin_like');self.assertEqual(out['grade'],'B')
        self.assertEqual(out['family_support_level'],'specific_supported')
        self.assertIn('ORTHOLOGY_UNRESOLVED',json.loads(out['identity_gate_trace'])['reason_codes'])

    def test_no_generic_ph_or_partial_family_promotion(self):
        out=audit_record(self.row(family='kindlin',architecture='PF00169:10-150',NCBI_CDD_architecture='PH:10-150'),self.rules)
        self.assertEqual(out['grade'],'D');self.assertEqual(out['family_support_level'],'unresolved')
        trace=json.loads(out['identity_gate_trace'])
        self.assertIn('MISSING_REQUIRED_DOMAIN',trace['reason_codes'])

    def test_export_keeps_rule_trace_outside_excel_and_emits_compact_candidate_trace(self):
        from openpyxl import Workbook,load_workbook
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);source=root/'input.xlsx';dest=root/'Stage_2_output.xlsx'
            row=self.row(family='generic_kinase_screen',length=500,architecture='PF08953:10-60;PF00400:100-150',
                         NCBI_CDD_architecture='DUF1899:10-60')
            book=Workbook();sheet=book.active;sheet.title='candidates';sheet.append(list(row));sheet.append(list(row.values()));book.save(source)
            run_xlsx(source,dest,'candidates',self.rules)
            book=load_workbook(dest,read_only=True);headers=next(book['Screening_hypothesis'].values)
            self.assertNotIn('_rule_gate_trace',headers);self.assertIn('identity_gate_trace',headers);book.close()
            with (root/'stage2_rule_gate_trace.tsv').open() as handle:trace=list(csv.DictReader(handle,delimiter='\t'))
            self.assertTrue(any(r['rule_id']=='coronin' and r['matched']=='True' for r in trace))
            with (root/'stage2_candidate_gate_trace.tsv').open() as handle:candidates=list(csv.DictReader(handle,delimiter='\t'))
            self.assertEqual(len(candidates),1);self.assertEqual(candidates[0]['assigned_family'],'coronin_like')
            self.assertLess(len(candidates[0]['identity_gate_trace']),32767)

    def test_synthetic_probability_fit_respects_declared_broad_name_scope(self):
        from openpyxl import Workbook
        from calibrate_identity import fit
        with tempfile.TemporaryDirectory(dir=ROOT) as td:
            root=Path(td);controls=root/'controls';controls.mkdir()
            book=Workbook();sheet=book.active;sheet.title='control_identity'
            sheet.append(['class','accession','genes','protein_name'])
            for key,subtype in [('A',4),('B',6)]:sheet.append(['Positive control',key,'ACTR'+str(subtype),'Actin-related protein '+str(subtype)])
            book.save(controls/'labels.xlsx');(controls/'controls.faa').write_text('>A\nMAAA\n>B\nMCCC\n')
            features=root/'features.tsv'
            with features.open('w',newline='') as handle:
                fields=['sequence_sha256','assigned_family','total_score_0_10','family_confidence','protein_name_scope','orthology_panel_sha256']
                writer=csv.DictWriter(handle,fieldnames=fields,delimiter='\t');writer.writeheader()
                for sequence in ('MAAA','MCCC'):writer.writerow(dict(sequence_sha256=hashlib.sha256(sequence.encode()).hexdigest(),
                    assigned_family='actin',total_score_0_10=7,family_confidence='strong',protein_name_scope='broad_family',orthology_panel_sha256='fixture'))
            fit(controls,features,root/'model.json',root/'report',minimum=2)
            model=json.loads((root/'model.json').read_text())
            self.assertEqual([r['correct'] for r in model['samples']],[1,1])
            self.assertIn('actin',model['context']['positive_family_scope'])


if __name__=='__main__':unittest.main()
