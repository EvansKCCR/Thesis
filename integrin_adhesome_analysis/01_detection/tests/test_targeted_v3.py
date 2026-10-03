import json
from pathlib import Path
import unittest
from adhesome_classification import architecture_evidence
from adhesome_motifs import load_library,scan_motifs,noncatalytic_alignment

ROOT=Path(__file__).resolve().parents[1]

class TargetedTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rules=json.loads((ROOT/'config/family_rules.json').read_text())['families']
        cls.library=load_library(ROOT/'references/motifs/adhesome_motif.json')

    def hits(self,*accessions):
        return [dict(pfam_accession=p,ali_from=1+i*100,ali_to=90+i*100) for i,p in enumerate(accessions)]

    def test_alpha_single_repeat_is_not_diagnostic(self):
        evidence=architecture_evidence(self.hits('PF01839'),self.rules['integrin_alpha'])
        self.assertTrue(evidence['partial'])
        self.assertFalse(evidence['diagnostic'])
        self.assertFalse(evidence['specific'])

    def test_divergent_alpha_stalk_without_two_repeats(self):
        evidence=architecture_evidence(self.hits('PF01839','PF20805'),self.rules['integrin_alpha'])
        self.assertTrue(evidence['complete'])

    def test_fak_uses_focal_at_not_pikk_fat(self):
        rule=self.rules['FAK']
        self.assertFalse(architecture_evidence(self.hits('PF00373','PF07714','PF02259'),rule)['complete'])
        self.assertTrue(architecture_evidence(self.hits('PF00373','PF07714','PF03623'),rule)['complete'])
        self.assertEqual(architecture_evidence(self.hits('PF03623','PF07714','PF00373'),rule)['order'],'incomplete_or_reordered')

    def test_kindlin_does_not_require_detectable_ph(self):
        self.assertTrue(architecture_evidence(self.hits('PF18124','PF00373'),self.rules['kindlin'])['complete'])
        self.assertFalse(architecture_evidence(self.hits('PF00373','PF00169'),self.rules['kindlin'])['complete'])

    def test_fak_specific_ordered_core_supports_incomplete_architecture(self):
        arch=architecture_evidence(self.hits('PF21477','PF07714'),self.rules['FAK'])
        self.assertFalse(arch['complete'])
        self.assertTrue(arch['specific'])
        self.assertTrue(arch['support_order'])
        self.assertFalse(architecture_evidence(self.hits('PF00373','PF07714'),self.rules['FAK'])['support_order'])
        self.assertFalse(architecture_evidence(self.hits('PF07714','PF21477'),self.rules['FAK'])['support_order'])

    def test_alpha_tail_context_and_coordinates(self):
        seq='GFFKR'+'A'*200+'GFFKR'
        obs=scan_motifs(seq,[],self.rules['integrin_alpha'],self.library)
        self.assertEqual([(r['start'],r['end']) for r in obs],[(206,210)])

    def test_alpha_tail_variants_are_reported_exactly(self):
        obs=scan_motifs('A'*200+'GFFHR',[],self.rules['integrin_alpha'],self.library)
        self.assertEqual(obs[0]['sequence'],'GFFHR')
        self.assertIn('variants',obs[0]['motif'])

    def test_midas_requires_domain_context(self):
        seq='DASAS'+'A'*200
        obs=scan_motifs(seq,[],self.rules['integrin_beta'],self.library)
        self.assertEqual(next(r for r in obs if r['motif']=='MIDAS')['state'],'not_tested_missing_domain')

    def test_hydroxyproline_is_not_encoded_as_letter_o(self):
        obs=scan_motifs('GFPGER'+'A'*100,self.hits('PF01391'),self.rules['collagen_like'],self.library)
        self.assertEqual(obs[0]['state'],'detected')
        self.assertIn('not_established',obs[0]['interpretation'])

    def test_ptp_catalytic_only_alignment_not_sufficient(self):
        spec=self.rules['PTP_PEST']['noncatalytic_alignment']
        support=[dict(similarity='true',reciprocal='true',subject_start='1',subject_end='90')]
        self.assertFalse(noncatalytic_alignment(support,self.hits('PF00102'),spec))
        support[0]['subject_end']='170'
        self.assertTrue(noncatalytic_alignment(support,self.hits('PF00102'),spec))
        self.assertFalse(noncatalytic_alignment(support,[],spec))

    def test_unknown_motif_mapping_fails(self):
        with self.assertRaises(ValueError):
            scan_motifs('AAAA',[],dict(motif_checks=[dict(entity='x',name='x',pattern='A')]),self.library)

if __name__=='__main__':
    unittest.main()
