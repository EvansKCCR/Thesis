"""Synthetic trusted-label regression cases; no calibration fit or analysis."""
from pathlib import Path
import sys
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from calibrate_identity import expected_families,expected_family_resolution,expected_broad_families
from label_ontology import named_family


class TrustedIdentityOntologyTests(unittest.TestCase):
    def test_binding_partners_do_not_become_gold_sequence_families(self):
        cases=[
            ('GAB1','GRB2-associated-binding protein 1','GAB_like',{'Grb2'}),
            ('ICAP1;ITGB1BP1','Integrin beta-1-binding protein 1','ICAP1_like',{'integrin_beta'}),
            ('ARC34;ARPC2','Actin-related protein 2/3 complex subunit 2','ARPC2_like',{'actin','arp2_like'}),
            ('PPM1F;POPX2','Protein phosphatase 1F','PPM_like',{'ILK'}),
            ('FBLIM1;FBLP1','Filamin-binding LIM protein 1','migfilin_like',{'filamin','zyxin_like'}),
            ('MACF1;ACF7','Microtubule-actin cross-linking factor 1','MACF1_spectraplakin_like',{'actin','actin_conventional','alpha_actinin'}),
        ]
        for genes,name,family,wrong in cases:
            with self.subTest(genes=genes):
                record=dict(genes=genes,protein_name=name)
                poisoned={gene:wrong for gene in genes.split(';')}
                resolution=expected_family_resolution(record,poisoned)
                self.assertEqual(resolution['families'],{family})
                self.assertFalse(resolution['family_label_review_required'])
                self.assertEqual(set(resolution['legacy_catalogue_families_ignored'].split(';')),wrong)

    def test_positive_identity_mapping_works_without_classifier_fields_or_catalogue(self):
        for genes,family in [('GAB1','GAB_like'),('ITGB1BP1','ICAP1_like'),('ARPC2','ARPC2_like'),
                ('PPM1F','PPM_like'),('FBLIM1','migfilin_like'),('MACF1','MACF1_spectraplakin_like')]:
            with self.subTest(genes=genes):
                result=expected_family_resolution(dict(genes=genes,protein_name=''))
                self.assertEqual(result['families'],{family})
                self.assertEqual(result['family_label_source'],'trusted_gene_identity_ontology')

    def test_actr2_is_not_arpc2_and_ancestors_have_separate_scope(self):
        arp=expected_families(dict(genes='ACTR2',protein_name='Actin-related protein 2'))
        arc=expected_families(dict(genes='ARPC2',protein_name='Actin-related protein 2/3 complex subunit 2'))
        self.assertEqual(arp,{'arp2_like'})
        self.assertEqual(expected_broad_families({},arp),{'actin'})
        self.assertEqual(arc,{'ARPC2_like'})
        self.assertEqual(expected_broad_families({},arc),set())
        # Other complex subunits are not silently treated as actin or ARP2.
        self.assertEqual(expected_families(dict(protein_name='Actin-related protein 2/3 complex subunit 3')),set())

    def test_testin_and_prickle_are_distinct_gold_identities(self):
        for gene,name,family in [('TES','Testin','Testin_like'),('PRICKLE1','Prickle-like protein 1','Prickle_like')]:
            with self.subTest(gene=gene):
                record=dict(genes=gene,protein_name=name)
                self.assertEqual(expected_families(record),{family})
                self.assertEqual(expected_broad_families(record),{'PET_LIM_like'})

    def test_explicit_trusted_family_always_takes_precedence(self):
        record=dict(genes='GAB1',protein_name='GRB2-associated-binding protein 1',trusted_family='curator_family; curator_scope')
        result=expected_family_resolution(record,{'GAB1':{'Grb2'}})
        self.assertEqual(result['families'],{'curator_family','curator_scope'})
        self.assertEqual(result['family_label_source'],'explicit_trusted_family')
        self.assertFalse(result['family_label_review_required'])

    def test_uncertain_gold_names_require_review_instead_of_catalogue_predictions(self):
        record=dict(genes='NOVEL1',protein_name='Uncharacterized adhesome-associated protein',
            assigned_family='talin',putative_audited_family='talin',family='talin',pfam='PF00373',NCBI_CDD='Talin')
        result=expected_family_resolution(record,{'NOVEL1':{'talin'}})
        self.assertEqual(result['families'],set())
        self.assertTrue(result['family_label_review_required'])
        self.assertEqual(result['family_label_source'],'requires_trusted_family_review')
        self.assertIn('not gold',result['family_label_review_reason'])
        self.assertEqual(result['legacy_catalogue_families_ignored'],'talin')

    def test_conflicting_trusted_own_gene_and_name_requires_curation(self):
        result=expected_family_resolution(dict(genes='GRB2',protein_name='GRB2-associated-binding protein 1'))
        self.assertEqual(result['families'],set())
        self.assertEqual(result['family_label_source'],'trusted_identity_conflict')
        self.assertTrue(result['family_label_review_required'])

    def test_current_own_identity_name_rules_remain_available(self):
        cases=[('Talin-1','talin'),('Integrin beta-1','integrin_beta'),('Integrin alpha-5','integrin_alpha'),
            ('Filamin-A','filamin'),('Integrin-linked protein kinase','ILK'),
            ('Growth factor receptor-bound protein 2','Grb2'),('Coronin-1B','coronin_like'),
            ('Cortactin','cortactin_HS1_like'),('Vasodilator-stimulated phosphoprotein','ENA_VASP_like'),
            ('Breast cancer anti-estrogen resistance protein 1','BCAR_CAS_like'),('Ezrin','ERM_like')]
        for name,family in cases:
            with self.subTest(name=name):self.assertIn(family,expected_families(dict(protein_name=name)))

    def test_unsupported_binding_descriptions_do_not_create_partner_identity(self):
        for name in ('Filamin-binding protein 9','GRB2-associated-binding protein X',
                'Integrin beta-1-binding protein 9','Integrin-linked kinase-associated phosphatase',
                'Fibronectin type III domain-containing protein 9'):
            with self.subTest(name=name):self.assertEqual(expected_families(dict(protein_name=name)),set())

    def test_broad_ancestors_are_explicit_and_do_not_match_arbitrary_substrings(self):
        self.assertEqual(expected_broad_families({}, {'Src','Csk','Testin_like','arp2_like','actin_conventional'}),
            {'Src_Csk_like','PET_LIM_like','actin'})
        self.assertEqual(expected_broad_families({}, {'ARPC2_like','GAB_like','other_named_identity'}),set())
        self.assertEqual(named_family(dict(protein_name='Testin')),'Testin_like')


if __name__=='__main__':unittest.main()
