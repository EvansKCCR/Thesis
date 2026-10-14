import argparse
import json
from pathlib import Path
import tempfile
import unittest
from adhesome_core import table, write_table, fasta
from adhesome_pipeline import classify
from adhesome_classification import architecture_evidence

class ComparativeTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.rule = dict(module='actin_linking',min_length=500,diagnostic=False,
            domain_alternatives=[{'PF00001':1,'PF00002':2}],ordered_domains=['PF00001','PF00002'])

    def tearDown(self):
        self.temp.cleanup()

    def classify(self, domains, length=200, evidence=None, rules=None):
        (self.root/'p.fa').write_text('>p\n'+'A'*length+'\n')
        (self.root/'r.json').write_text(json.dumps(dict(status='test',families=rules or {'test':self.rule})))
        rows = [dict(query_name='p',pfam_accession=domain,ali_from=1+i*20,ali_to=10+i*20,decision='retained') for i,domain in enumerate(domains)]
        write_table(self.root/'d.csv',rows,['query_name','pfam_accession','ali_from','ali_to','decision'])
        evidence_path = None
        if evidence:
            evidence_path = self.root/'e.csv'
            write_table(evidence_path,evidence,['protein_id','family','similarity','reciprocal','profile','phylogeny','topology','motif','reference_tier','source'])
        args = argparse.Namespace(proteome=self.root/'p.fa',rules=self.root/'r.json',domains=self.root/'d.csv',
            evidence=[evidence_path] if evidence_path else [],homology=None,out=self.root/'out',manifest=None)
        classify(args)
        return table(self.root/'out/candidates.csv')

    def test_short_incomplete_protein_retained(self):
        self.rule['screening_domains']=['PF00002']
        row = self.classify(['PF00002'])[0]
        self.assertEqual(row['status'],'provisional')
        self.assertEqual(row['assigned_family'],'')
        self.assertIn('below_expected_length',row['warnings'])
        self.assertIn('PF00001:1',row['missing_domains'])
        self.assertEqual(len(list(fasta(self.root/'out/candidates.faa'))),1)

    def test_reordered_domains_are_flagged_not_excluded(self):
        row = self.classify(['PF00002','PF00002','PF00001'])[0]
        self.assertEqual(row['status'],'provisional')
        self.assertEqual(row['order_check'],'incomplete_or_reordered')

    def test_generic_domain_retained_but_unassigned(self):
        self.rule['domain_alternatives']=[{'PF00001':1}]
        row = self.classify(['PF00001'])[0]
        self.assertEqual(row['status'],'unresolved')
        self.assertFalse(row['assigned_family'])
        self.assertEqual(len(list(fasta(self.root/'out/candidates.faa'))),1)
        self.assertEqual(list(fasta(self.root/'out/prioritized_candidates.faa')),[])

    def test_profile_only_recovery_without_domains(self):
        evidence=[dict(protein_id='p',family='test',profile='true',source='family_HMM')]
        row = self.classify([],evidence=evidence)[0]
        self.assertEqual(row['status'],'provisional')
        self.assertEqual(row['profile'],'supported')
        self.assertNotIn('topology_check',row)
        self.assertNotIn('phylogeny',row)

    def test_topology_is_deferred_to_resume(self):
        self.rule['requires_topology']=True
        self.rule['motifs']=['NP.Y']
        evidence=[dict(protein_id='p',family='test',similarity='true',reciprocal='true',source='reference_search')]
        row = self.classify(['PF00001','PF00002','PF00002'],evidence=evidence)[0]
        self.assertEqual(row['status'],'supported_candidate')
        self.assertEqual(row['assigned_family'],'test')
        self.assertEqual(row['motif_check'],'not_detected')
        evidence.append(dict(protein_id='p',family='test',topology='false',source='topology_prediction'))
        with self.assertRaisesRegex(ValueError,'resume'):
            self.classify(['PF00001','PF00002','PF00002'],evidence=evidence)

    def test_conflicting_family_assignments_are_not_resolved_arbitrarily(self):
        evidence=[dict(protein_id='p',family=family,similarity='true',reciprocal='true',source='search') for family in ['test','other']]
        rows = self.classify(['PF00001','PF00002','PF00002'],evidence=evidence,rules={'test':self.rule,'other':self.rule})
        self.assertTrue(all(r['assigned_family']=='' and r['status']=='provisional' for r in rows))

    def test_anchor_outside_expected_architecture_recovers_fragment(self):
        self.rule['screening_domains']=['PF00003']
        row = self.classify(['PF00003'])[0]
        self.assertEqual(row['status'],'provisional')
        self.assertEqual(row['matched_domains'],'PF00003')

    def test_invalid_evidence_fails(self):
        with self.assertRaisesRegex(ValueError,'true, false'):
            self.classify([],evidence=[dict(protein_id='p',family='test',profile='yes',source='test')])

    def test_generic_screening_anchor_does_not_confer_family_specificity(self):
        self.rule['screening_domains']=['PF00002']
        self.rule['screening_specific']=False
        row = self.classify(['PF00002'])[0]
        self.assertEqual(row['status'],'unresolved')
        self.assertEqual(row['assigned_family'],'')

    def test_unmodelled_curated_domain_is_preserved(self):
        self.classify(['PF00009'])
        manifest=self.root/'m.csv'
        manifest.write_text('Pfam accession,Domain or family\nPF00009,Generic domain\n')
        args=argparse.Namespace(proteome=self.root/'p.fa',rules=self.root/'r.json',domains=self.root/'d.csv',
            evidence=[],homology=None,out=self.root/'out',manifest=manifest)
        classify(args)
        rows=table(self.root/'out/candidates.csv')
        self.assertEqual(rows[0]['family'],'unassigned_domain_candidate')
        self.assertEqual(rows[0]['status'],'unresolved')

    def test_strict_mode_remains_available(self):
        self.classify(['PF00001','PF00002','PF00002'])
        args=argparse.Namespace(proteome=self.root/'p.fa',rules=self.root/'r.json',domains=self.root/'d.csv',
            evidence=[],homology=None,out=self.root/'strict',classification_mode='strict')
        classify(args)
        self.assertEqual(table(self.root/'strict/candidates.csv')[0]['status'],'excluded')

    def test_generic_bucket_never_promoted_by_homology(self):
        self.rule['screening_only']=True
        evidence=[dict(protein_id='p',family='test',similarity='true',reciprocal='true',source='search')]
        row=self.classify(['PF00001','PF00002','PF00002'],evidence=evidence)[0]
        self.assertEqual(row['status'],'unresolved')
        self.assertEqual(row['assigned_family'],'')

    def test_ptp_subfamily_requires_family_evidence(self):
        self.rule['requires_family_evidence']=True
        rows=self.classify(['PF00001','PF00002','PF00002'])
        self.assertFalse(any(r['family']=='test' for r in rows))
        evidence=[dict(protein_id='p',family='test',similarity='true',source='search')]
        rows=self.classify(['PF00001'],evidence=evidence)
        self.assertEqual(rows[0]['family'],'test')
        self.assertEqual(rows[0]['status'],'provisional')

    def test_computational_anchor_only_is_not_confirmatory(self):
        evidence=[dict(protein_id='p',family='test',similarity='true',reciprocal='true',reference_tier='computational_anchor',source='anchor')]
        row=self.classify(['PF00001','PF00002','PF00002'],evidence=evidence)[0]
        self.assertEqual(row['status'],'provisional')
        self.assertIn('computational_anchor_only',row['warnings'])
        evidence.append(dict(protein_id='p',family='test',phylogeny='true',source='independent_tree'))
        with self.assertRaisesRegex(ValueError,'resume'):
            self.classify(['PF00001','PF00002','PF00002'],evidence=evidence)

    def test_targeted_domain_order_gate_keeps_reordered_hit_provisional(self):
        self.rule['ordered_domain_groups']=[['PF00001'],['PF00002']]
        self.rule['support_requires_order']=True
        evidence=[dict(protein_id='p',family='test',similarity='true',reciprocal='true',source='search')]
        row=self.classify(['PF00002','PF00002','PF00001'],evidence=evidence)[0]
        self.assertEqual(row['status'],'provisional')
        row=self.classify(['PF00001','PF00002','PF00002'],evidence=evidence)[0]
        self.assertEqual(row['status'],'supported_candidate')

    def test_ptp_old_homology_without_coordinates_stays_provisional(self):
        self.rule['noncatalytic_alignment']=dict(domains=['PF00001'],minimum_residues=60)
        evidence=[dict(protein_id='p',family='test',similarity='true',reciprocal='true',source='old_search')]
        row=self.classify(['PF00001','PF00002','PF00002'],evidence=evidence)[0]
        self.assertEqual(row['status'],'provisional')
        self.assertIn('noncatalytic_family_evidence_missing',row['warnings'])

    def test_motif_alone_never_recruits_a_protein(self):
        self.rule['motifs']=['AAA']
        self.assertEqual(self.classify([]),[])

    def test_diagnostic_ordered_core_supports_family_with_missing_terminal_domain(self):
        self.rule.update(domain_alternatives=[{'PF00001':1,'PF00002':1,'PF00003':1}],
            diagnostic_domains=['PF00001'],support_requires_order=True,
            ordered_domain_groups=[['PF00001'],['PF00002'],['PF00003']],
            support_order_alternatives=[[['PF00001'],['PF00002']]])
        evidence=[dict(protein_id='p',family='test',similarity='true',reciprocal='true',source='search')]
        row=self.classify(['PF00001','PF00002'],evidence=evidence)[0]
        self.assertEqual(row['status'],'supported_candidate')
        self.assertEqual(row['support_order_check'],'consistent')
        self.assertIn('incomplete_expected_architecture',row['warnings'])

if __name__=='__main__':
    unittest.main()
