
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
import json
from pathlib import Path
import tempfile
import unittest
from integrated_adhesome_pipeline import audit_record,classify_rows,run_xlsx,checksum,architecture,motif_score

ROOT=Path(__file__).resolve().parents[2]


class TwoStageAuditTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.rules=json.loads((ROOT/'detection/config/family_rules.json').read_text())

    def row(self,**extra):
        row=dict(protein_id='AnyTarget__g',Specie='AnyTarget',family='profilin',length=150,status='supported_candidate',
            domain_match=True,domain_type_fraction=1,domain_copy_fraction=1,architecture='PF00235:10-130',
            Orthology_HOG_status='shared HOG',Orthology_HOG_reference_species_count=3,Orthology_analysis_completed=True)
        row.update(extra);return row

    def arch(self,*models):return ';'.join(f'{model}:{10+i*60}-{50+i*60}' for i,model in enumerate(models))

    def test_full_architecture_and_orthology_grade_a_with_ten_point_scale(self):
        result=audit_record(self.row(MotifPrediction_context_supported_count=2,**{'DeepLoc_2.1':'Cytoplasm'}),self.rules)
        self.assertEqual([result[k] for k in ('domain_score_0_4','orthology_score_0_3','motif_score_0_2','localization_score_0_1','total_score_0_10')],[4,3,2,1,10])
        self.assertEqual(result['grade'],'A');self.assertEqual(result['assigned_family'],'profilin')

    def test_tested_target_only_complete_architecture_is_b_without_hardcoded_species(self):
        result=audit_record(self.row(Orthology_HOG_reference_species_count=0,Orthology_HOG_status='target-only in sampled HOG'),self.rules)
        self.assertEqual(result['grade'],'B');self.assertEqual(result['family_decision'],'target_panel_restricted')
        self.assertNotIn('Schistosoma',result['rationale'])

    def test_missing_orthology_is_not_called_observed_panel_restriction(self):
        result=audit_record(self.row(Orthology_HOG_reference_species_count=0,Orthology_HOG_status='',Orthology_analysis_completed=False),self.rules)
        self.assertEqual(result['grade'],'B');self.assertEqual(result['family_decision'],'orthology_unavailable')

    def test_unconfirmed_partial_stays_unresolved_and_legacy_profile_is_explicit(self):
        row=self.row(family='alpha_actinin',length=900,architecture='PF00307:10-110')
        self.assertEqual(audit_record(row,self.rules)['grade'],'D')
        self.assertEqual(audit_record(row,dict(self.rules,_stage2_profile='original'))['grade'],'C')
        row['Orthology_HOG_reference_species_count']=0
        self.assertEqual(audit_record(row,self.rules)['grade'],'D')

    def test_no_architecture_cannot_be_rescued_by_exact_identity_motifs_or_many_references(self):
        result=audit_record(self.row(architecture='',reference_subfamily='profilin',subfamily_support_method='curated_reference',
            reference_accession='EXACT',Computational_reference_evidence='[{"exact_sequence_match":true}]',MotifPrediction_context_supported_count=50),self.rules)
        self.assertEqual(result['grade'],'D');self.assertEqual(result['assigned_family'],'unassigned')

    def test_catalogue_is_not_searched_to_replace_the_supplied_family(self):
        result=audit_record(self.row(family='Nck'),self.rules)
        self.assertEqual(result['preaudit_family'],'Nck');self.assertEqual(result['grade'],'D')
        self.assertNotEqual(result['assigned_family'],'profilin')

    def test_generic_screening_bucket_cannot_become_an_assigned_family(self):
        result=audit_record(self.row(family='LIM_scaffold_screen',architecture=self.arch(*(['PF00412']*5)),length=450),self.rules)
        self.assertEqual(result['grade'],'D');self.assertEqual(result['assigned_family'],'unassigned')
        self.assertIn('screening bucket',result['architecture_evidence'])

    def test_specific_family_cdd_and_order_gates(self):
        cases=dict(PINCH_like=(['PF00412']*5,'LIM1_PINCH'),parvin_like=(['PF00307']*2,'CH_PARVIN'),
            Syk=(['PF00017','PF00017','PF07714'],'SH2_SYK'),Csk=(['PF00018','PF00017','PF07714'],'SH3_CSK'),
            Nck=(['PF00018','PF00018','PF00017'],'SH3_NCK2'),Grb2=(['PF00018','PF00017','PF00018'],'SH2_GRB2'),
            Crk=(['PF00017','PF00018'],'SH2_CRK'),Shc=(['PF00640','PF00017'],'PTB_SHC'))
        for family,(models,cdd) in cases.items():
            with self.subTest(family=family):
                row=self.row(family=family,length=450,architecture=self.arch(*models),NCBI_CDD_architecture=cdd)
                self.assertEqual(audit_record(row,self.rules)['domain_score_0_4'],4)
                self.assertEqual(audit_record(row,self.rules)['grade'],'A')
                row['NCBI_CDD_architecture']='generic domain'
                self.assertNotEqual(audit_record(row,self.rules)['grade'],'A')

    def test_nck_reordered_domains_do_not_pass_the_complete_gate(self):
        row=self.row(family='Nck',length=450,architecture=self.arch('PF00017','PF00018','PF00018'),NCBI_CDD_architecture='SH3_NCK')
        self.assertLess(audit_record(row,self.rules)['domain_score_0_4'],4)

    def test_pinch_length_and_copy_number_boundaries(self):
        for length,full in ((249,False),(250,True),(500,True),(501,False)):
            # Fit the five distinct domains within 250 aa.
            arch=';'.join(f'PF00412:{i*45+1}-{i*45+35}' for i in range(5))
            with self.subTest(length=length):self.assertEqual(architecture(self.row(length=length,architecture=arch,NCBI_CDD_architecture='PINCH'),'PINCH_like',self.rules)['score']==4,full)
        row=self.row(length=450,architecture=self.arch(*(['PF00412']*6)),NCBI_CDD_architecture='PINCH')
        self.assertNotEqual(architecture(row,'PINCH_like',self.rules)['score'],4)

    def test_src_csk_core_does_not_name_csk_without_specific_cdd(self):
        row=self.row(family='Src',length=500,architecture=self.arch('PF00018','PF00017','PF07714'))
        self.assertEqual(audit_record(row,self.rules)['audit_family'],'Src_Csk_like')
        row['NCBI_CDD_architecture']='SH3_CHK'
        self.assertEqual(audit_record(row,self.rules)['audit_family'],'Csk')

    def test_only_complete_pinch_and_parvin_recoveries_can_be_added_to_existing_pool(self):
        rows=[self.row(family='paxillin_like',length=450,architecture=self.arch(*(['PF00412']*5)),NCBI_CDD_architecture='LIM_PINCH'),
              self.row(protein_id='AnyTarget__h',family='alpha_actinin',length=450,architecture=self.arch('PF00307','PF00307'),NCBI_CDD_architecture='CH_PARVIN')]
        audited,resolved=classify_rows(rows,self.rules)
        self.assertEqual(len(audited),4);self.assertEqual(len(resolved),2)
        self.assertEqual({r['protein_id'] for r in resolved},{r['protein_id'] for r in rows})
        self.assertEqual({r['audit_family'] for r in resolved},{'PINCH_like','parvin_like'})
        self.assertTrue(all(r['recovered_hypothesis_count']==1 for r in resolved))
        rows[0]['NCBI_CDD_architecture']='generic LIM';rows[1]['architecture']=self.arch('PF00307')
        audited,resolved=classify_rows(rows,self.rules);self.assertEqual(len(audited),2)

    def test_original_alternatives_and_species_keys_are_preserved_without_mutating_inputs(self):
        rows=[self.row(family='profilin'),self.row(family='Nck'),self.row(Specie='AnotherTarget')]
        before=json.dumps(rows,sort_keys=True);audited,resolved=classify_rows(rows,self.rules)
        self.assertEqual(json.dumps(rows,sort_keys=True),before);self.assertEqual(len(audited),3);self.assertEqual(len(resolved),2)
        first=next(r for r in resolved if r['Specie']=='AnyTarget')
        self.assertEqual(first['collapsed_row_count'],2);self.assertEqual(first['alternative_family_hypotheses'],'Nck')

    def test_benchmark_weights_orthology_without_bypassing_architecture(self):
        result=audit_record(self.row(Orthology_benchmark_status='validated',Orthology_benchmark_weight=.5),self.rules)
        self.assertEqual(result['orthology_raw_score_0_3'],3);self.assertEqual(result['orthology_score_0_3'],1.5)
        result=audit_record(self.row(Orthology_benchmark_status='validated',Orthology_benchmark_weight=0),self.rules)
        self.assertEqual(result['grade'],'B');self.assertNotEqual(result['family_decision'],'target_panel_restricted')

    def test_adhesome_weight_is_measured_per_exact_audited_family(self):
        fields=dict(Orthology_benchmark_status='validated',
            Orthology_benchmark_scope='adhesome_subset_reference_orthogroup_inference',
            Orthology_benchmark_weight=.01,
            Orthology_benchmark_family_weights=json.dumps({'profilin':.5,'fibronectin_like':0.0}))
        result=audit_record(self.row(**fields),self.rules)
        self.assertEqual(result['orthology_score_0_3'],1.5)
        self.assertEqual(result['Orthology_benchmark_weight'],.5)
        self.assertEqual(result['orthology_benchmark_weight_status'],'validated_adhesome_family_weight')
        fields['Orthology_benchmark_family_weights']='{"fibronectin_like":0.2}'
        result=audit_record(self.row(**fields),self.rules)
        self.assertEqual(result['orthology_score_0_3'],3)
        self.assertEqual(result['Orthology_benchmark_weight'],'')
        self.assertEqual(result['orthology_benchmark_weight_status'],'family_not_benchmarked_neutral_weight')

    def test_subset_zero_weight_and_invalid_weights_do_not_leak_to_untested_families(self):
        from integrated_adhesome_pipeline import benchmark_weight
        row=dict(Orthology_benchmark_status='validated',
            Orthology_benchmark_scope='adhesome_subset_reference_orthogroup_inference',
            Orthology_benchmark_family_weights='{"fibronectin_like":0}')
        self.assertEqual(benchmark_weight(row,'fibronectin_like'),(0,'validated_adhesome_family_weight'))
        self.assertEqual(benchmark_weight(row,'fibronectin_III_ligand'),(1,'family_not_benchmarked_neutral_weight'))
        for invalid in ('[]','bad JSON','{"profilin":-1}','{"profilin":"nan"}'):
            row['Orthology_benchmark_family_weights']=invalid
            self.assertEqual(benchmark_weight(row,'profilin'),(1,'invalid_weight_ignored'))

    def test_refined_family_does_not_inherit_baseline_subset_weight(self):
        row=self.row(family='Src',length=500,architecture=self.arch('PF00018','PF00017','PF07714'),NCBI_CDD_architecture='SH3_CHK',
            Orthology_benchmark_status='validated',Orthology_benchmark_scope='adhesome_subset_reference_orthogroup_inference',
            Orthology_benchmark_family_weights='{"Src":0}')
        result=audit_record(row,self.rules)
        self.assertEqual(result['audit_family'],'Csk')
        self.assertEqual(result['orthology_score_0_3'],3)
        self.assertEqual(result['orthology_benchmark_weight_status'],'family_not_benchmarked_neutral_weight')

    def test_motif_context_sequence_and_conflict_scores(self):
        for row,expected in (({'MotifPrediction_context_supported_count':2},2),({'MotifPrediction_sequence_only_count':5},1),({'MotifPrediction_candidate_tiers':'context_conflict'},0)):
            self.assertEqual(motif_score(row),expected)

    def test_workbook_matches_user_sheet_contract_and_dynamic_reference_columns(self):
        from openpyxl import Workbook,load_workbook
        with tempfile.TemporaryDirectory() as td:
            source=Path(td)/'input.xlsx';dest=Path(td)/'Stage_2_output.xlsx'
            book=Workbook();sheet=book.active;sheet.title='candidates';row=self.row()
            sheet.append(list(row));sheet.append(list(row.values()));book.save(source);before=checksum(source)
            support=dict(reference_species=['NewRef'],target_species=['AnyTarget'],benchmark={'Orthology_benchmark_status':'validated'},tables={})
            meta=run_xlsx(source,dest,'candidates',self.rules,support)
            book=load_workbook(dest,read_only=True)
            self.assertEqual(book.sheetnames,['Pipeline_summary','Screening_hypothesis','Audited_family_assignment','Motif_annotations','Orthology_candidate_report','Orthology_HOGs','Direct_orthologues','NCBI_CDD_selected_hits','Topology_predictions'])
            fields=next(book['Audited_family_assignment'].values)
            self.assertNotIn('Orthology_Hsap_HOG_copy_count',fields)
            self.assertIn('putative_audited_family',fields);self.assertIn('protein_name',fields)
            self.assertIn('collapsed_row_count',fields)
            self.assertNotIn('audit_family',fields);self.assertNotIn('assigned_family',fields)
            book.close();self.assertEqual(checksum(source),before);self.assertTrue(meta['protein_pool_preserved'])


if __name__=='__main__':unittest.main()
