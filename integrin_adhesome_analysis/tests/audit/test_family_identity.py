
import sys as _sys
from pathlib import Path as _Path
_PIPELINE_ROOT=_Path(__file__).resolve().parents[2]
_sys.path.insert(0,str(_PIPELINE_ROOT))
_sys.path.insert(0,str(_PIPELINE_ROOT/"05_conservative_adhesome_family_assignment_audit"))
import json
from pathlib import Path
import unittest
import tempfile
from integrated_adhesome_pipeline import audit_record,classify_rows,architecture
from family_identity import combined_domains
from audit_report import prepare_sheet

ROOT=Path(__file__).resolve().parents[2]


class FamilyIdentityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.rules=json.loads((ROOT/'01_detection/config/family_rules.json').read_text())

    def row(self,**extra):
        row=dict(protein_id='Any__gene',Specie='Any',family='talin',length=2200,status='provisional',
            architecture='PF00373:10-200',NCBI_CDD_architecture='B41:10-200;I_LWEQ:2000-2150;Talin_middle:600-900',
            Orthology_HOG_reference_species_count=3,Orthology_analysis_completed=True)
        row.update(extra);return row

    def test_observed_cdd_completes_talin_architecture_and_supports_putative_identity(self):
        out=audit_record(self.row(),self.rules)
        self.assertEqual(out['assigned_family'],'talin');self.assertEqual(out['grade'],'A')
        self.assertEqual(out['family_confidence'],'strong');self.assertEqual(out['protein_name'],'Putative Talin')
        self.assertEqual(out['protein_name_scope'],'family_lineage')

    def test_an_inferred_region_does_not_complete_the_architecture(self):
        out=audit_record(self.row(NCBI_CDD_architecture='B41:10-200;[Inferred: I_LWEQ:2000-2150]'),self.rules)
        self.assertEqual(out['grade'],'D');self.assertNotEqual(out['family_confidence'],'strong')

    def test_complete_family_gate_and_orthology_are_strong_without_optional_cdd_or_experiments(self):
        row=self.row(family='actin',length=400,architecture='PF00022:10-390',NCBI_CDD_architecture='')
        out=audit_record(row,self.rules)
        self.assertEqual(out['grade'],'A');self.assertTrue(out['domain_gate_passed'])
        self.assertEqual(out['assigned_family'],'actin');self.assertEqual(out['family_confidence'],'strong')
        self.assertEqual(out['protein_name'],'Putative Actin-family protein')

    def test_experimental_flags_cannot_change_assignment_scores_name_or_confidence(self):
        row=self.row(family='actin',length=400,architecture='PF00022:10-390',NCBI_CDD_architecture='')
        variants=[{},dict(experimental_evidence='',adhesome_experiment_verified=False,experimental_annotation_count=0),
            dict(experimental_evidence='experimental annotation',adhesome_experiment_verified=True,experimental_annotation_count=10)]
        results=[audit_record(dict(row,**extra),self.rules) for extra in variants]
        for field in ('assigned_family','protein_name','family_confidence','protein_identity_confidence','grade',
                'domain_score_0_4','orthology_score_0_3','total_score_0_10'):
            self.assertEqual(len({result[field] for result in results}),1,field)

    def test_experimental_flags_cannot_bypass_a_required_family_gate(self):
        row=self.row(family='Syk',length=600,architecture='PF00017:30-100;PF07714:250-500',
            NCBI_CDD_architecture='',experimental_evidence='experiment',adhesome_experiment_verified=True)
        out=audit_record(row,self.rules)
        self.assertEqual(out['grade'],'D');self.assertFalse(out['domain_gate_passed'])
        self.assertNotEqual(out['family_confidence'],'strong')

    def test_cdd_label_alone_does_not_replace_missing_domain_core(self):
        out=audit_record(self.row(NCBI_CDD_architecture='Talin_middle:600-900'),self.rules)
        self.assertEqual(out['grade'],'D')

    def test_overlapping_pfam_cdd_footprints_do_not_multiply_lim_copies(self):
        row=self.row(length=450,architecture=';'.join(f'PF00412:{i*60+1}-{i*60+45}' for i in range(5)),
            NCBI_CDD_architecture=';'.join(f'LIM{i+1}_PINCH:{i*60+2}-{i*60+44}' for i in range(5)))
        from integrated_adhesome_pipeline import domains
        hits=combined_domains(row,domains(row))
        self.assertEqual(sum(h['model']=='PF00412' for h in hits),5)
        self.assertEqual(architecture(row,'PINCH_like',self.rules)['score'],4)

    def test_hip1_related_identity_requires_enth_and_ilweq(self):
        out=audit_record(self.row(architecture='PF01608:700-850',length=900,
            NCBI_CDD_architecture='ENTH:20-150;HIP1:20-150;I_LWEQ:700-850'),self.rules)
        self.assertEqual(out['assigned_family'],'HIP1_HIP1R_like');self.assertEqual(out['grade'],'A')
        self.assertIn('HIP1/HIP1R',out['protein_name'])

    def test_arp2_name_is_based_on_diagnostic_cdd_not_target_identifier(self):
        out=audit_record(self.row(family='actin',length=400,architecture='PF00022:10-390',
            NCBI_CDD_architecture='ASKHA_NBD_actin_Arp2:10-390'),self.rules)
        self.assertEqual(out['assigned_family'],'arp2_like');self.assertIn('Arp2-like',out['protein_name'])
        out=audit_record(self.row(protein_id='Any__ARP2',family='actin',length=400,architecture='PF00022:10-390',NCBI_CDD_architecture=''),self.rules)
        self.assertEqual(out['assigned_family'],'actin');self.assertNotIn('Arp2',out['protein_name'])

    def test_competing_arp_signatures_remain_ambiguous(self):
        out=audit_record(self.row(family='actin',length=400,architecture='PF00022:10-390',
            NCBI_CDD_architecture='actin_Arp4_Arp6:10-390'),self.rules)
        self.assertEqual(out['assigned_family'],'actin');self.assertEqual(out['family_support_level'],'broad_supported')
        self.assertEqual(out['specific_identity_status'],'specific_lineage_or_paralogue_unresolved')
        self.assertEqual(out['unresolved_specific_families'],'arp4_like;arp6_like')

    def test_src_specific_cdd_refines_core_to_src_instead_of_src_csk(self):
        out=audit_record(self.row(family='Src',length=550,architecture='PF00018:20-80;PF00017:100-190;PF07714:240-500',
            NCBI_CDD_architecture='SH3_Src:20-80;SH2_Src_Src:100-190;PTKc_Src:240-500'),self.rules)
        self.assertEqual(out['assigned_family'],'Src');self.assertEqual(out['family_confidence'],'strong')

    def test_single_fn3_cannot_name_robo_or_eph_receptor(self):
        out=audit_record(self.row(family='fibronectin_III_ligand',length=900,architecture='PF00041:100-170',
            NCBI_CDD_architecture='FN3:100-170;ROBO:100-170'),self.rules)
        self.assertEqual(out['assigned_family'],'unassigned')
        self.assertIn('fibronectin type III domain-containing',out['protein_name'])

    def test_scoped_refinement_does_not_search_unrelated_families_or_add_proteins(self):
        row=self.row(family='profilin')
        audited,resolved=classify_rows([row],self.rules)
        self.assertEqual(len(resolved),1);self.assertEqual(resolved[0]['protein_id'],row['protein_id'])
        self.assertNotEqual(resolved[0]['assigned_family'],'talin')

    def test_clean_final_view_removes_unused_duplicate_and_blank_columns(self):
        row=audit_record(self.row(),self.rules)
        row.update(MULocDeep_Subcellular_predicted_scores='0.5',unused_blank='',Orthology_direct_reference_count=0,
            collapsed_row_count=1,alternative_family_hypotheses='')
        records,fields,removed=prepare_sheet('Audited_family_assignment',[row])
        self.assertIn('protein_name',fields);self.assertIn('putative_audited_family',fields)
        self.assertNotIn('assigned_family',fields);self.assertNotIn('audit_family',fields)
        self.assertNotIn('MULocDeep_Subcellular_predicted_scores',fields)
        self.assertNotIn('alternative_family_hypotheses',fields)
        self.assertIn('Orthology_direct_reference_count',fields)
        self.assertEqual(records[0]['putative_audited_family'],'talin')
        self.assertNotIn('family_probability',fields)

    def test_supporting_views_keep_zero_false_and_remove_only_duplicate_aliases(self):
        rows=[dict(motif_id='M1',library_id='M1',reported_candidate_tier='sequence_only',candidate_tier='sequence_only',
            score=0,passed=False,blank=' ',MULocDeep_other='unused')]
        _,fields,_=prepare_sheet('Motif_annotations',rows)
        self.assertIn('score',fields);self.assertIn('passed',fields)
        self.assertNotIn('library_id',fields);self.assertNotIn('reported_candidate_tier',fields)
        self.assertNotIn('blank',fields);self.assertNotIn('MULocDeep_other',fields)

    def test_cdd_sidecar_replaces_excel_truncated_json_without_losing_architecture(self):
        from openpyxl import Workbook,load_workbook
        from integrated_adhesome_pipeline import run_xlsx
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);source=root/'input.xlsx';dest=root/'result.xlsx'
            row=self.row(NCBI_CDD_hits='[{"truncated":')
            book=Workbook();sheet=book.active;sheet.title='candidates'
            sheet.append(list(row));sheet.append(list(row.values()));book.save(source)
            support={'tables':{'NCBI_CDD_selected_hits':[
                dict(sequence_id=row['protein_id'],short_name='B41',accession='smart00295',start=10,end=200),
                dict(sequence_id=row['protein_id'],short_name='I_LWEQ',accession='cd00001',start=2000,end=2150)]}}
            run_xlsx(source,dest,'candidates',self.rules,support)
            book=load_workbook(dest,read_only=True);sheet=book['Audited_family_assignment']
            fields=next(sheet.values);record=dict(zip(fields,next(sheet.iter_rows(min_row=2,max_row=2,values_only=True))))
            self.assertEqual(record['putative_audited_family'],'talin');book.close()
            self.assertTrue((root/'audited_family_assignments.tsv').exists())
            with self.assertRaisesRegex(ValueError,'supporting-data'):run_xlsx(source,root/'no_sidecar.xlsx','candidates',self.rules)


if __name__=='__main__':unittest.main()
