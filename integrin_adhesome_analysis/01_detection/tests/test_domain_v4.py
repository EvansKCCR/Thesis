import argparse,json,tempfile,unittest
from pathlib import Path
from adhesome_classification import architecture_evidence,classify_comparative
from adhesome_core import table,write_table

ROOT=Path(__file__).resolve().parents[1]

class DomainV4Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rules=json.loads((ROOT/'config/family_rules.json').read_text())['families']

    def hits(self,*domains):
        return [dict(pfam_accession=p,ali_from=10+i*100,ali_to=80+i*100) for i,p in enumerate(domains)]

    def test_all_requested_models_are_in_active_library(self):
        models=json.loads((ROOT/'databases/pfam/adhesome_v1.0.hmm.metadata.json').read_text())['models']
        requested={'PF00039','PF00040','PF01033','PF00045','PF01839','PF08441','PF20805','PF20806',
                   'PF21477','PF00373','PF18038','PF03623','PF07714','PF12796'}
        self.assertTrue(requested<=models.keys())
        self.assertTrue(all(models[p]['ga'] for p in requested))

    def test_ilk_ank2_is_an_alternative_repeat_model(self):
        arch=architecture_evidence(self.hits('PF12796','PF07714'),self.rules['ILK'])
        self.assertTrue(arch['complete'])
        self.assertTrue(arch['support_order'])
        self.assertFalse(architecture_evidence(self.hits('PF12796'),self.rules['ILK'])['complete'])
        self.assertFalse(architecture_evidence(self.hits('PF07714','PF12796'),self.rules['ILK'])['support_order'])

    def test_ilk_legacy_alternative_is_preserved(self):
        self.assertTrue(architecture_evidence(self.hits('PF00023','PF00023','PF00069'),self.rules['ILK'])['complete'])

    def test_fibronectin_requires_a_combination_not_single_repeat(self):
        rule=self.rules['fibronectin_like']
        for p in ('PF00039','PF00040','PF00041'):
            arch=architecture_evidence(self.hits(p),rule)
            self.assertTrue(arch['partial'])
            self.assertFalse(arch['complete'])
            self.assertFalse(arch['specific'])
        self.assertTrue(architecture_evidence(self.hits('PF00039','PF00040','PF00041'),rule)['complete'])

    def test_vitronectin_requires_ordered_smb_hemopexin(self):
        rule=self.rules['vitronectin_like']
        for p in ('PF01033','PF00045'):
            self.assertFalse(architecture_evidence(self.hits(p),rule)['specific'])
        self.assertTrue(architecture_evidence(self.hits('PF01033','PF00045'),rule)['support_order'])
        self.assertFalse(architecture_evidence(self.hits('PF00045','PF01033'),rule)['support_order'])

    def test_fak_ferm_n_can_complete_ordered_architecture(self):
        rule=self.rules['FAK']
        self.assertTrue(architecture_evidence(self.hits('PF18038','PF07714','PF03623'),rule)['complete'])
        self.assertFalse(architecture_evidence(self.hits('PF18038'),rule)['specific'])

    def test_generic_ecm_screen_never_assigns_family(self):
        self.assertTrue(self.rules['generic_ECM_domain_screen']['screening_only'])
        self.assertTrue(self.rules['fibronectin_III_ligand']['screening_only'])

    def test_ilk_assignment_needs_reciprocal_and_ordered_domains(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder)
            (root/'p.fa').write_text('>p\n'+'A'*500+'\n')
            (root/'rules.json').write_text(json.dumps(dict(status='test',families={'ILK':self.rules['ILK']})))
            evidence=dict(protein_id='p',family='ILK',similarity='true',reciprocal='true',source='reference_panel',
                reference_subfamily='ILK',subfamily_support_method='curated_reference',reference_accession='FIXTURE')
            write_table(root/'e.csv',[evidence],list(evidence))
            args=argparse.Namespace(proteome=root/'p.fa',rules=root/'rules.json',domains=root/'d.csv',
                homology=root/'e.csv',evidence=[],manifest=None,out=root/'out')
            for ordered,status in [(('PF12796','PF07714'),'supported_candidate'),(('PF07714','PF12796'),'provisional')]:
                rows=[dict(query_name='p',**r) for r in self.hits(*ordered)]
                write_table(root/'d.csv',rows,list(rows[0]))
                classify_comparative(args)
                self.assertEqual(table(root/'out/candidates.csv')[0]['status'],status)

    def test_ecm_single_domains_not_promoted_even_with_similarity(self):
        cases=[('vitronectin_like',('PF01033',),'provisional'),
               ('vitronectin_like',('PF01033','PF00045'),'supported_candidate'),
               ('vitronectin_like',('PF00045','PF01033'),'provisional'),
               ('fibronectin_like',('PF00041',),'provisional'),
               ('fibronectin_like',('PF00039','PF00040','PF00041'),'supported_candidate')]
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder)
            (root/'p.fa').write_text('>p\n'+'A'*500+'\n')
            for family,domains,status in cases:
                config=dict(status='test',families={family:self.rules[family]},
                            motif_library=str(ROOT.parent/'03_functional_motif_annotation/Motif_library/adhesion_core.json'))
                (root/'rules.json').write_text(json.dumps(config))
                rows=[dict(query_name='p',**r) for r in self.hits(*domains)]
                write_table(root/'d.csv',rows,list(rows[0]))
                evidence=dict(protein_id='p',family=family,similarity='true',reciprocal='true',source='panel',
                    reference_subfamily=family,subfamily_support_method='curated_reference',reference_accession='FIXTURE')
                write_table(root/'e.csv',[evidence],list(evidence))
                classify_comparative(argparse.Namespace(proteome=root/'p.fa',rules=root/'rules.json',domains=root/'d.csv',
                    homology=root/'e.csv',evidence=[],manifest=None,out=root/'out'))
                self.assertEqual(table(root/'out/candidates.csv')[0]['status'],status)

if __name__=='__main__':
    unittest.main()
