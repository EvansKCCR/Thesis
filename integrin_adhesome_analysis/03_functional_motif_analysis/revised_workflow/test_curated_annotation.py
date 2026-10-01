import csv
import hashlib
import json
import re
import tempfile
import unittest
import zipfile
from collections import defaultdict
from pathlib import Path
import curated_annotation as ca
import family_motifs as fm


class CuratedTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rules=json.loads((fm.HERE/'Motif_library/context_rules.json').read_text())
        cls.motifs,cls.sheets,cls.audit=ca.load_library(fm.HERE/'Motif_library/adhesion_protein_SLiM_library.xlsx',cls.rules)

    def test_workbook_and_routing(self):
        self.assertEqual(len(self.motifs),36)
        self.assertEqual(sum(m['reader_only'] for m in self.motifs),5)
        self.assertEqual(self.rules['motifs']['SLIM001']['target_families'],['integrin_beta'])
        self.assertNotIn('profilin',self.rules['motifs']['SLIM020']['target_families'])
        self.assertEqual(self.rules['motifs']['SLIM035']['target_families'],['integrin_beta'])
        self.assertTrue(all(m['compiled'] is None for m in self.motifs if m['reader_only']))

    def test_precursor_and_terminal_anchor(self):
        motif=next(m for m in self.motifs if m['motif_id']=='SLIM012')
        self.assertEqual(motif['scan_pattern'],'GFPGER')
        self.assertEqual(list(ca.occurrences(motif['compiled'],'GFOGER')),[])
        self.assertEqual(list(ca.occurrences(re.compile('YQPG$'),'YQPGE')),[])
        self.assertEqual(list(ca.occurrences(re.compile('P..P'),'PPPPP')),[(1,4,'PPPP'),(2,5,'PPPP')])
        self.assertEqual(list(ca.occurrences(re.compile('NP.Y'),'NPXY')),[])

    def test_score_is_not_candidate_confidence(self):
        m=next(m for m in self.motifs if m['motif_id']=='SLIM017')
        r=self.rules['motifs'][m['motif_id']]
        record=dict(id='x',sequence='ACDENPAY')
        result=ca.assess(m,r,'integrin_beta',record,5,8,[('cytoplasmic_tail',5,8)],'domain_count_rule_supported',{})
        self.assertEqual(m['confidence_score'],95)
        self.assertEqual(result[0],'context_supported_candidate')
        self.assertIn('positional_conservation',result[2])

    def test_high_gate_and_review_cap(self):
        m=next(m for m in self.motifs if m['motif_id']=='SLIM017')
        rule=self.rules['motifs'][m['motif_id']]
        record=dict(id='x',sequence='ACDENPAY')
        evidence={(('integrin_beta','x'),m['motif_id']):[dict(criterion=k,start=5,end=8,status='supported',source='test_reference') for k in rule['requirements']]}
        args=(rule,'integrin_beta',record,5,8,[('cytoplasmic_tail',5,8)],'domain_count_rule_supported',evidence)
        self.assertEqual(ca.assess(m,*args)[0],'high_confidence_hypothesis')
        self.assertEqual(ca.assess(m|{'screening_recommendation':'Review only'},*args)[0],'review_only')
        evidence[(('integrin_beta','x'),m['motif_id'])].append(dict(criterion='orthology',start=5,end=8,status='contradicted',source='test_conflict'))
        self.assertEqual(ca.assess(m,*args)[0],'context_conflict')

    def test_outside_not_proven_extracellular(self):
        self.assertEqual(ca.localization({'location':'extracellular'},2,4,10,[('predicted_outside',1,10)]),'provisional')
        self.assertEqual(ca.localization({'location':'cytoplasmic_tail'},2,4,10,[('cytoplasmic_tail',7,10)]),'contradicted')

    def test_topology_exact_sequences_and_tail(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'predicted_topologies.3line'
            p.write_text('>x | TM\nACDEFGHI\nOOMMMIII\n>y\nACDE\nIIII\n')
            fam={'a':[dict(id='x',sequence='ACDEFGHI'),dict(id='y',sequence='ACDF')]}
            rs=defaultdict(list)
            audit=ca.import_topology([p],fam,rs)
            self.assertEqual([r['status'] for r in audit],['exact_sequence_verified','sequence_mismatch'])
            self.assertIn(('cytoplasmic_tail',6,8),rs['a','x'])
            self.assertEqual(rs['a','y'],[])

    def test_stale_evidence_fails(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'e.tsv'
            p.write_text('family\tsequence_id\tsequence_sha256\tmotif_id\tstart\tend\tcriterion\tstatus\tsource\na\tx\twrong\tSLIM017\t1\t4\torthology\tsupported\tdoc\n')
            with self.assertRaises(ValueError): ca.load_evidence(p,{'a':[dict(id='x',sequence='NPAY')]},self.rules)

    def test_zip_topology_without_extraction(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'top.zip'
            with zipfile.ZipFile(p,'w') as archive:
                archive.writestr('nested/predicted_topologies.3line','>x\nACDEFGHI\nOOMMMIII\n')
            rs=defaultdict(list)
            rows=ca.import_topology([p],{'a':[dict(id='x',sequence='ACDEFGHI')]},rs)
            self.assertEqual(rows[0]['status'],'exact_sequence_verified')
            self.assertIn('!nested/predicted_topologies.3line',rows[0]['source'])
            self.assertFalse((Path(d)/'nested').exists())

    def test_reader_no_scan_duplicate_site_group(self):
        motifs=[m for m in self.motifs if m['motif_id'] in ('SLIM001','SLIM017','SLIM023')]
        fam={'integrin_beta':[dict(id='b',sequence='NPAY')],'kindlin':[dict(id='k',sequence='NPAY')]}
        domains=[dict(family=f,sequence_id=r['id'],state='domain_count_rule_supported') for f,rs in fam.items() for r in rs]
        with tempfile.TemporaryDirectory() as d:
            out=Path(d)
            result=ca.run(motifs,self.sheets,self.audit,self.rules,fam,defaultdict(list),domains,{},[],out,fm.write_tsv)
            self.assertEqual(result['annotation_rows'],2)
            self.assertEqual(result['distinct_sites'],1)
            self.assertEqual(result['high_confidence_rows'],0)
            with (out/'reader_hypotheses.tsv').open() as h:
                rows=list(csv.DictReader(h,delimiter='\t'))
            self.assertEqual(len(rows),1)
            self.assertEqual(rows[0]['sequence_id'],'k')


if __name__=='__main__': unittest.main()
