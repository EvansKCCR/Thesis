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
        cls.library=load_library(ROOT.parent/'03_functional_motif_annotation/Motif_library/adhesion_core.json')

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
        mid=next(k for k,r in self.library['current_rules'].items() if r['regex_python']=='GFFKR')
        obs=scan_motifs(seq,[],dict(motif_rule_ids=[mid]),self.library)
        self.assertEqual([(r['start'],r['end']) for r in obs],[(1,5),(206,210)])
        self.assertTrue(all(r['interpretation']=='sequence_candidate_only' for r in obs))

    def test_alpha_tail_variants_are_reported_exactly(self):
        mid=next(k for k,r in self.library['current_rules'].items() if r['regex_python']=='GFFKR')
        self.assertEqual(scan_motifs('A'*200+'GFFHR',[],dict(motif_rule_ids=[mid]),self.library),[])

    def test_midas_requires_domain_context(self):
        seq='DASAS'+'A'*200
        obs=scan_motifs(seq,[],self.rules['integrin_beta'],self.library)
        self.assertEqual(obs,[])

    def test_hydroxyproline_is_not_encoded_as_letter_o(self):
        mid=next(k for k,r in self.library['current_rules'].items() if 'SLIM012' in r['source_row_ids'])
        obs=scan_motifs('GFOGER',self.hits('PF01391'),dict(motif_rule_ids=[mid]),self.library)
        self.assertEqual(obs,[])

    def test_ptp_catalytic_only_alignment_not_sufficient(self):
        spec=self.rules['PTP_PEST']['noncatalytic_alignment']
        support=[dict(similarity='true',reciprocal='true',subject_start='1',subject_end='90')]
        self.assertFalse(noncatalytic_alignment(support,self.hits('PF00102'),spec))
        support[0]['subject_end']='170'
        self.assertTrue(noncatalytic_alignment(support,self.hits('PF00102'),spec))
        self.assertFalse(noncatalytic_alignment(support,[],spec))

    def test_unknown_motif_mapping_fails(self):
        with self.assertRaises(ValueError):
            scan_motifs('AAAA',[],dict(motif_rule_ids=['missing']),self.library)

if __name__=='__main__':
    unittest.main()
