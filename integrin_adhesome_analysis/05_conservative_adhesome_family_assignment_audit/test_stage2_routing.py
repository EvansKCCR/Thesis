"""Synthetic routing, observed-domain and coordinate boundary regression tests."""
import json
from pathlib import Path
import unittest
from unittest.mock import patch

from integrated_adhesome_pipeline import audit_record,architecture,domains
from family_identity import (combined_domains,distinct_footprints,refinement_trace,
                             scoped_refinement,diagnostic_cdd,cdd_hits,library)

ROOT=Path(__file__).resolve().parents[1]


class Stage2RoutingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.rules=json.loads((ROOT/'01_detection/config/family_rules.json').read_text())

    def row(self,**values):
        row=dict(protein_id='Novel__sequence',family='generic_kinase_screen',length=1000,
                 architecture='',NCBI_CDD_architecture='',status='provisional',
                 Orthology_HOG_reference_species_count=2)
        row.update(values);return row

    def trace(self,row,rule_id):
        return next(t for t in refinement_trace(row,row['family'],combined_domains(row,domains(row))) if t['rule_id']==rule_id)

    def test_reviewed_components_dispatch_from_unrelated_preliminary_hypothesis(self):
        cases=[
            ('coronin','coronin_like','PF08953:5-70;PF00400:80-120','DUF1899:5-70'),
            ('cortactin','cortactin_HS1_like','PF02218:20-50;PF02218:80-110;PF14604:600-650','SH3_Cortactin:600-650'),
            ('ena_vasp','ENA_VASP_like','PF00568:5-110;PF08776:700-750','EVH1_Ena_VASP-like:5-110'),
            ('bcar_cas','BCAR_CAS_like','PF00018:4-60;PF12026:700-850','SH3_BCAR1:4-60'),
            ('erm','ERM_like','PF00373:50-200;PF00769:700-800','FERM_C_ERM:220-300'),
        ]
        for rule_id,family,pfam,cdd in cases:
            with self.subTest(rule=rule_id):
                row=self.row(architecture=pfam,NCBI_CDD_architecture=cdd)
                out=audit_record(row,self.rules)
                self.assertEqual(out['assigned_family'],family)
                self.assertEqual(out['grade'],'A')
                trace=self.trace(row,rule_id)
                self.assertTrue(trace['matched']);self.assertEqual(trace['routing_scope'],'diagnostic_dispatch')

    def test_dispatch_requires_diagnostic_and_observed_pfam_anchor(self):
        row=self.row(architecture='PF00568:5-110;PF08776:700-750')
        trace=self.trace(row,'ena_vasp')
        self.assertFalse(trace['reached']);self.assertIn('DIAGNOSTIC_CDD_UNCONFIRMED',trace['reason_codes'])
        row=self.row(NCBI_CDD_architecture='EVH1_Ena_VASP-like:5-110;VASP_tetra:700-750')
        trace=self.trace(row,'ena_vasp')
        self.assertFalse(trace['matched']);self.assertIn('PFAM_ANCHOR_MISSING',trace['reason_codes'])

    def test_dispatch_never_uses_target_gene_name_or_calibration_label(self):
        row=self.row(protein_id='VASP',genes='VASP',trusted_family='ENA_VASP_like',protein_name='VASP',
                     architecture='PF00568:5-110;PF08776:700-750')
        self.assertFalse(self.trace(row,'ena_vasp')['matched'])

    def test_unreviewed_legacy_rule_stays_within_its_parent_scope(self):
        row=self.row(family='profilin',length=2200,architecture='PF00373:10-200;PF01608:1800-2000',
                     NCBI_CDD_architecture='Talin_middle:500-900')
        trace=self.trace(row,'talin')
        self.assertTrue(trace['diagnostic_cdd_matched']);self.assertFalse(trace['reached'])
        self.assertFalse(trace['matched']);self.assertIn('NO_IDENTITY_RULE_REACHED',trace['reason_codes'])

    def test_testin_uses_pfam_pet_and_three_ordered_distinct_lim_footprints(self):
        row=self.row(family='paxillin_like',length=450,
            architecture='PF06297:109-193;PF00412:236-293;PF00412:299-354',
            NCBI_CDD_architecture='LIM1_Testin:236-293;LIM2_Testin:299-354;LIM3_Testin:361-419')
        hits=combined_domains(row,domains(row))
        self.assertTrue(any(h['model']=='PET' and 'Pfam' in h['sources'] for h in hits))
        out=audit_record(row,self.rules)
        self.assertEqual(out['assigned_family'],'Testin_like');self.assertEqual(out['grade'],'A')
        trace=self.trace(row,'testin')
        self.assertEqual(trace['domain_groups'][1]['observed'],3)
        self.assertFalse(self.trace(row,'prickle')['matched'])
        row['architecture']='PF00412:20-75;PF00412:80-135;PF06297:200-280'
        row['NCBI_CDD_architecture']='LIM1_Testin:20-75;LIM2_Testin:80-135;LIM3_Testin:140-195'
        self.assertIn('DOMAIN_ORDER_CONFLICT',self.trace(row,'testin')['reason_codes'])

    def test_prickle_keeps_its_separate_five_lim_rule(self):
        row=self.row(family='PET_LIM_like',architecture='PF06297:10-100;'+
            ';'.join(f'PF00412:{200+i*70}-{250+i*70}' for i in range(5)),
            NCBI_CDD_architecture='PET_prickle:10-100')
        out=audit_record(row,self.rules)
        self.assertEqual(out['assigned_family'],'Prickle_like')
        row['architecture']='PF06297:10-100;'+ ';'.join(f'PF00412:{200+i*70}-{250+i*70}' for i in range(3))
        self.assertFalse(self.trace(row,'prickle')['matched'])
        self.assertIn('MISSING_REQUIRED_DOMAIN',self.trace(row,'prickle')['reason_codes'])

    def test_receptor_cdd_alias_completes_second_ptp_without_pest_assignment(self):
        row=self.row(family='PTP_PEST',architecture='PF00102:210-500',
            NCBI_CDD_architecture='R-PTPc-A-1:214-509;R-PTPc-A-2:563-790')
        out=audit_record(row,self.rules)
        self.assertEqual(out['assigned_family'],'RPTP_alpha_epsilon');self.assertEqual(out['grade'],'A')
        row=self.row(family='PTP_PEST',architecture='PF00102:40-276',NCBI_CDD_architecture='PTPc-N1:16-292')
        self.assertEqual(audit_record(row,self.rules)['assigned_family'],'PTPN1_PTP1B_like')

    def test_observed_workbook_model_layouts_and_anchored_aliases(self):
        row=self.row(NCBI_CDD_architecture='unrelated [1-4]; cd14608 PTPc-N1 [16-292]')
        self.assertTrue(diagnostic_cdd(row,r'^PTPc[-_]N1$'))
        self.assertEqual(cdd_hits(row)[1]['short_name'],'PTPc-N1')
        self.assertEqual(combined_domains(row,[])[0]['model'],'PF00102')
        row['NCBI_CDD_architecture']='smart00132 LIM [3 regions]'
        self.assertEqual(cdd_hits(row),[])
        self.assertEqual(combined_domains(row,[]),[])

    def test_overlapping_alternative_models_and_cdd_do_not_multiply_copy_counts(self):
        row=self.row(architecture='PF02218:10-80;PF02218:12-78;PF14604:600-650',
            NCBI_CDD_architecture='HS1_rep:11-79;SH3_Cortactin:600-650;SH3_9:600-650')
        trace=self.trace(row,'cortactin')
        self.assertEqual(trace['domain_groups'][0]['observed'],1)
        self.assertFalse(trace['matched']);self.assertIn('MISSING_REQUIRED_DOMAIN',trace['reason_codes'])
        hits=combined_domains(row,domains(row))
        self.assertEqual(len(distinct_footprints(hits,{'PF00018','PF14604'})),1)

    def test_malformed_or_outside_sequence_cdd_cannot_supply_a_diagnostic(self):
        row=self.row(length=400,architecture='PF00568:5-110;PF08776:340-380',
                     NCBI_CDD_architecture='EVH1_Ena_VASP-like:500-600')
        self.assertFalse(diagnostic_cdd(row,'EVH1_Ena_VASP'))
        self.assertFalse(self.trace(row,'ena_vasp')['matched'])

    def test_trace_distinguishes_missing_domains_order_and_coordinates(self):
        row=self.row(architecture='PF00568:5-110',NCBI_CDD_architecture='EVH1_Ena_VASP-like:5-110')
        trace=self.trace(row,'ena_vasp')
        self.assertTrue(trace['reached']);self.assertIn('MISSING_REQUIRED_DOMAIN',trace['reason_codes'])
        row['architecture']='PF08776:5-45;PF00568:100-200'
        trace=self.trace(row,'ena_vasp')
        self.assertFalse(trace['order_passed']);self.assertIn('DOMAIN_ORDER_CONFLICT',trace['reason_codes'])
        row['architecture']='PF00568;PF08776:340-380'
        trace=self.trace(row,'ena_vasp')
        self.assertFalse(trace['coordinates_passed'])
        self.assertIn('DOMAIN_COORDINATES_UNCONFIRMED',trace['reason_codes'])
        json.dumps(trace)

    def test_multiple_complete_diagnostic_dispatches_preserve_competing_families(self):
        row=self.row(architecture='PF08953:5-50;PF00400:60-100;PF00568:200-300;PF08776:700-750',
                     NCBI_CDD_architecture='DUF1899:5-50;EVH1_Ena_VASP-like:200-300')
        result=scoped_refinement(row,row['family'],combined_domains(row,domains(row)))
        self.assertIsNone(result['rule']);self.assertEqual(result['competing_families'],['ENA_VASP_like','coronin_like'])


if __name__=='__main__':unittest.main()
