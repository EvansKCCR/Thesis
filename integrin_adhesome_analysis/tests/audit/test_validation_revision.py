"""Control-independent synthetic examples of diagnostic family boundaries."""

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
import unittest
from unittest.mock import patch
from integrated_adhesome_pipeline import audit_record,classify_rows
from calibrate_identity import expected_families
from control_exclusion import apply

ROOT=Path(__file__).resolve().parents[2]


class ValidationRevisionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.rules=json.loads((ROOT/'detection/config/family_rules.json').read_text())

    def row(self,**values):
        row=dict(protein_id='Novel__g',Specie='Novel',family='PTP_PEST',length=800,status='provisional',
                 architecture='PF00102:300-550',NCBI_CDD_architecture='',Orthology_HOG_reference_species_count=2)
        row.update(values);return row

    def test_ptpn1_and_receptor_ptp_diagnostic_models_do_not_receive_pest_names(self):
        ptpn=audit_record(self.row(NCBI_CDD_architecture='PTPc:290-550;PTPc-N1:300-550'),self.rules)
        self.assertEqual(ptpn['assigned_family'],'PTPN1_PTP1B_like');self.assertEqual(ptpn['grade'],'A')
        receptor=audit_record(self.row(architecture='PF00102:150-350;PF00102:450-680',
            NCBI_CDD_architecture='R-PTPc-A-1:150-350;R-PTPc-A-2:450-680'),self.rules)
        self.assertEqual(receptor['assigned_family'],'RPTP_alpha_epsilon');self.assertEqual(receptor['grade'],'A')
        self.assertNotEqual(audit_record(self.row(),self.rules)['grade'],'A')

    def test_macf_refinement_requires_distinct_ch_spectrin_and_terminal_gas2(self):
        row=self.row(family='alpha_actinin',length=4000,
            architecture='PF00307:10-100;PF00307:120-210;PF00435:500-600;PF00435:700-800;PF00435:900-1000;PF02187:3500-3600',
            NCBI_CDD_architecture='CH_MACF1_rpt2:120-210;GAS2:3500-3600')
        out=audit_record(row,self.rules)
        self.assertEqual(out['assigned_family'],'MACF1_spectraplakin_like');self.assertEqual(out['grade'],'A')
        row['NCBI_CDD_architecture']='CH:120-210'
        self.assertNotEqual(audit_record(row,self.rules)['grade'],'A')

    def test_generic_bucket_recovers_ena_vasp_only_with_complete_ordered_core(self):
        row=self.row(family='component_reference_candidate',length=500,
            architecture='PF00568:10-100;PF08776:400-450',NCBI_CDD_architecture='EVH1_Ena_VASP-like:10-100')
        audited,resolved=classify_rows([row],self.rules)
        self.assertEqual([r['protein_id'] for r in resolved],['Novel__g'])
        self.assertEqual(resolved[0]['assigned_family'],'ENA_VASP_like');self.assertEqual(resolved[0]['grade'],'A')
        row['architecture']='PF00568:10-100'
        self.assertEqual(audit_record(row,self.rules)['assigned_family'],'unassigned')
        row['architecture']='PF08776:10-50;PF00568:100-200'
        self.assertEqual(audit_record(row,self.rules)['assigned_family'],'unassigned')

    def test_unseen_erm_and_fhl_can_be_called_without_experimental_or_training_metadata(self):
        erm=self.row(family='talin',length=650,architecture='PF00373:20-200;PF00769:500-600',
                     NCBI_CDD_architecture='FERM_C_ERM:210-300')
        self.assertEqual(audit_record(erm,self.rules)['assigned_family'],'ERM_like')
        fhl=self.row(family='FHL_screen',length=300,
            architecture=';'.join(f'PF00412:{i*60+10}-{i*60+55}' for i in range(4)),NCBI_CDD_architecture='LIM1_FHL2:10-55')
        self.assertEqual(audit_record(fhl,self.rules)['assigned_family'],'FHL_like')
        fhl['NCBI_CDD_architecture']='generic LIM'
        self.assertEqual(audit_record(fhl,self.rules)['assigned_family'],'unassigned')
        with patch('identity_probability.annotate',return_value={'status':'no_labelled_probability_model'}):
            _,out=classify_rows([erm],self.rules)
        self.assertEqual(out[0]['grade'],'A');self.assertEqual(out[0]['family_confidence'],'strong')

    def test_anchored_cdd_model_pattern_matches_a_model_later_in_the_list(self):
        row=self.row(family='actin',length=400,architecture='PF00022:10-350',
                     NCBI_CDD_architecture='unrelated:1-5;Actin_cytoplasmic:10-350')
        self.assertEqual(audit_record(row,self.rules)['assigned_family'],'actin_conventional')

    def test_truth_does_not_treat_actin_associated_names_as_actin(self):
        self.assertEqual(expected_families({'protein_name':'Microtubule-actin cross-linking factor 1'},{}),{'MACF1_spectraplakin_like'})
        self.assertEqual(expected_families({'protein_name':'Coronin-1B'},{}),{'coronin_like'})
        self.assertNotIn('actin',expected_families({'protein_name':'Actin-binding protein without resolved family'},{}))
        self.assertEqual(expected_families({'protein_name':'Talin','trusted_family':'curated_override'},{}),{'curated_override'})

    def test_exact_and_diagnostic_negative_gates_are_both_preserved(self):
        negative={'class':'hard_negative control','accession':'NEG','negative_specific_models':{'specific_NEG'}}
        bank=dict(negative_by_hash={'hash':[negative]},negative_families={'PTP_PEST':[negative]},positive_hashes=set())
        row=self.row(sequence_sha256='hash',assigned_family='PTP_PEST',protein_name='Putative PTP',
                     domain_score_0_4=4,family_specific_CDD_support=True,NCBI_CDD_architecture='specific_NEG:300-550')
        with patch('control_exclusion.controls',return_value=bank):kept,excluded=apply([row])
        self.assertEqual(kept,[])
        self.assertEqual(excluded[0]['exclusion_gates'].split(';'),
            ['exact_negative_control_sequence','diagnostic_negative_subfamily_and_complete_architecture'])


if __name__=='__main__':unittest.main()
