"""Regression tests for sequence feature symmetry and named-family resolution."""
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest

ML = Path(__file__).resolve().parents[2] / 'resource_library/adhesome_ml'
sys.path.insert(0, str(ML))
try:
    import sklearn
except ImportError:
    raise unittest.SkipTest('ML optional dependencies are not installed')
import classifier as c
import family_gate as g


class FamilyGateTests(unittest.TestCase):
    def setUp(self):
        self.row = {'protein_id':'query','ml_seq_len_aa':'100',
            'ml_domain_hits_json':json.dumps([{'source':'pfam','id':'PF00001','name':'PF00001',
                'start':1,'end':20,'score':42,'profile_coverage':.8,'evalue':1e-12}])}
        self.policy = dict(membership_probability_cutoff=.5,specific_family_probability_cutoff=.5,
            specific_family_probability_margin=.15,broad_parent_probability_cutoff=.6,
            minimum_family_examples=2)

    def test_hypotheses_and_audit_labels_cannot_change_features(self):
        cached = ([self.row],{},[])
        first = [{'protein_id':'query','stage1_family_hypotheses':'wrong','expected_families':'A'}]
        second = [{'protein_id':'query','stage1_family_hypotheses':'different','expected_families':'B'}]
        a,_=g.extract(first,{'query':'ACDACD'},None,[('m','ACD')],set(),c,cached)
        b,_=g.extract(second,{'query':'ACDACD'},None,[('m','ACD')],set(),c,cached)
        self.assertEqual(a,b)
        self.assertEqual(len(json.loads(a[0][0]['ml_motif_hits_json'])),2)
        self.assertTrue(any('copy_number' in key for key in a[0][1]))

    def test_domain_quality_survives_deduplication(self):
        vector=c.feature_dict(self.row,set())
        self.assertIn(42, [v for k,v in vector.items() if 'best_score' in k])
        self.assertIn(.8, [v for k,v in vector.items() if 'best_profile_coverage' in k])

    def test_competing_families_specific_and_broad(self):
        profiles={'A':['PF00001'],'B':['PF00001']}
        counts={'A':3,'B':3}
        args=(self.row,profiles,counts,{'A':'parent','B':'parent'},self.policy,c)
        self.assertEqual(g.decide(.9,[('B',.8),('A',.2)],*args)[:2],('specific','B'))
        self.assertEqual(g.decide(.9,[('A',.51),('B',.49)],*args)[:2],('broad','parent'))
        self.assertEqual(g.decide(.1,[('A',.9),('B',.1)],*args)[0],'non_adhesome')

    def test_missing_sequence_is_not_dropped(self):
        with self.assertRaisesRegex(ValueError,'No verified sequence'):
            g.extract([{'protein_id':'missing'}],{},None,[],set(),c,([],{},[]))

    def test_no_features_and_single_example_abstain(self):
        self.assertEqual(g.decide(.9,[('A',1)],{'protein_id':'x'},{'A':[]},{'A':3},{},self.policy,c)[0],'unresolved')
        self.assertEqual(g.decide(.9,[('A',1)],self.row,{'A':['PF00001']},{'A':1},{},self.policy,c)[0],'domain_only')

    def test_reference_labels_only(self):
        self.assertEqual(g.trusted_family({'predicted_family':'A','stage1_family_hypotheses':'A'},{},{}),'')
        self.assertEqual(g.trusted_family({'expected_families':'A;B'},{},{'A':'parent','B':'parent'}),'parent')

    def test_blast_custom_coordinates_and_directions(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory=Path(tmp)/'detection'; folder=directory/'species'/'homology'; folder.mkdir(parents=True)
            (folder/'forward.tsv').write_text('ref\tquery\t60\t80\t1\t80\t11\t90\t1e-20\t150\t100\t200\n')
            (folder/'reverse.tsv').write_text('query\tref\t60\t80\t11\t90\t1\t80\t1e-20\t150\t200\t100\n')
            rows,_,_=c.read_stage1_evidence(directory)
            self.assertEqual(len(rows),1)
            hits=json.loads(rows[0]['ml_homology_hits_json'])
            self.assertEqual(hits[0]['evalue'],1e-20)
            self.assertEqual(hits[0]['bitscore'],150)
            self.assertEqual(hits[0]['query_coverage'],.8)
            self.assertEqual(hits[1]['query_coverage'],.4)
            self.assertEqual(rows[0]['ml_reciprocal_pair_count'],'1')


if __name__=='__main__':
    unittest.main()
