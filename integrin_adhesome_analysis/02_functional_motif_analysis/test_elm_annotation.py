import csv
import json
import re
import tempfile
import unittest
from collections import defaultdict
from pathlib import Path
import elm_annotation as elm
import family_motifs as fm


class ElmTests(unittest.TestCase):
    def test_full_sequence_anchors_and_overlap(self):
        self.assertEqual(list(elm.matches(re.compile('^M?G..'), 'MGACMGAC')),[(1,4,'MGAC')])
        self.assertEqual(list(elm.matches(re.compile('P..P'), 'PPPPP')),[(1,4,'PPPP'),(2,5,'PPPP')])
        self.assertEqual(list(elm.matches(re.compile('RGD$'), 'ARGDA')),[])
        self.assertEqual(list(elm.matches(re.compile('RGD$'), 'ARGD')),[(2,4,'RGD')])
        self.assertEqual(list(elm.matches(re.compile('[^P]RGD'), 'XRGD')),[])

    def test_actual_resource_integrity_and_src_instances(self):
        resources=elm.load_resources(fm.HERE/'elm_resources',fm.read_fasta)
        self.assertEqual(len(resources['classes']),353)
        rows, exact=elm.reference_audit(resources)
        self.assertEqual(len(rows),5)
        self.assertTrue(all(r['state']=='exact_regex_span_recovered' for r in rows))
        self.assertEqual(len(exact),5)
        self.assertTrue(any(r['status']=='unresolved_class_references' for r in resources['audit']))
        self.assertEqual(len(resources['tables']['elm_kegg.csv']),15252)

    def test_comments_quoted_tsv_with_csv_suffix(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'test.csv'
            p.write_text('#version: 1\n"elm"\t"name"\n"LIG_X"\t"a, b"\n')
            rows,metadata=elm.read_table(p,{'elm','name'})
            self.assertEqual(rows,[{'elm':'LIG_X','name':'a, b'}])
            self.assertEqual(metadata,['#version: 1'])
            with self.assertRaises(ValueError):
                elm.read_table(p,{'missing'})

    def test_reference_no_coordinate_transfer(self):
        resources=dict(classes={'x':{'compiled':re.compile('RGD')}},references={'P1':'ARGDA'},
            tables={'elm_instances.tsv':[dict(Accession='I1', ELMIdentifier='x',Primary_Acc='P1',Start='2',End='4',InstanceLogic='true positive')]})
        _,exact=elm.reference_audit(resources)
        self.assertIn(('ARGDA','x',2,4),exact)
        self.assertNotIn(('CRGDC','x',2,4),exact)
        resources['tables']['elm_instances.tsv'][0]['InstanceLogic']='false positive'
        self.assertFalse(elm.reference_audit(resources)[1])

    def test_context_not_partner_domain(self):
        rule=dict(context='c_terminal_window',window=10,regions=['cytoplasmic_tail'])
        self.assertEqual(elm.focus_context(rule,95,99,100,[('PF00017',1,100)]),'terminal_window_only')
        self.assertEqual(elm.focus_context(rule,95,99,100,[('cytoplasmic_tail',96,100)]),'outside_expected_region')

    def test_end_to_end_original_site_annotation_and_background(self):
        resources=dict(classes={}, references={},audit=[],tables={name:[] for name in elm.RESOURCE_SCHEMAS})
        for ident,pattern in [('RGD','RGD'),('END','RGD$')]:
            resources['classes'][ident]=dict(Accession=ident,Regex=pattern,compiled=re.compile(pattern),Probability='0.01',Description='test')
        resources['tables']['elm_interaction_domains.tsv']=[{'ELM identifier':'RGD','Interaction Domain Id':'PF1'}]
        fam={'a':[dict(id='x',sequence='ARGDA'),dict(id='dup',sequence='ARGDA')],
             'b':[dict(id='shared',sequence='ARGDA'),dict(id='b',sequence='ACDEF')]}
        sites=[dict(family='a',sequence_id='x',motif='RGD',start=2,end=4)]
        focus={'families':{'a':[dict(elm_class='RGD',note='test')]}}
        regions=defaultdict(list,{('b','b'):[('PF1',1,5)]})
        with tempfile.TemporaryDirectory() as d:
            out=Path(d)
            summary=elm.run(resources,fam,regions,focus,sites,out,fm.write_tsv)
            def read(name):
                with (out/name).open(encoding='utf-8',newline='') as handle:
                    return list(csv.DictReader(handle,delimiter='\t'))
            self.assertEqual(summary['regex_hit_rows'],3)
            self.assertTrue(all(r['elm_class']=='RGD' for r in read('elm_discovery_annotation.tsv')))
            comparison=next(r for r in read('elm_family_prevalence.tsv') if r['family']=='a' and r['elm_class']=='RGD')
            self.assertEqual((comparison['target_total'],comparison['other_total'],comparison['other_positive']),('1','1','0'))
            self.assertEqual(read('elm_partner_domain_candidates.tsv')[0]['sequence_id'],'b')
            self.assertTrue(all(r['evidence']=='regex_candidate_only' for r in read('elm_hits.tsv')))


if __name__=='__main__':
    unittest.main()
