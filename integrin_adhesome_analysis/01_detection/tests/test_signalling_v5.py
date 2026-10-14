import argparse
import json
from pathlib import Path
import tempfile
import unittest
from adhesome_core import table, write_table, sequence_dict, digest
from adhesome_classification import architecture_evidence, classify_comparative

ROOT = Path(__file__).resolve().parents[1]
CASES = {
    'Syk': ['PF00017', 'PF00017', 'PF07714'],
    'Csk': ['PF00018', 'PF00017', 'PF07714'],
    'Nck': ['PF00018', 'PF00018', 'PF00017'],
    'Grb2': ['PF00018', 'PF00017', 'PF00018'],
    'Crk': ['PF00017', 'PF00018'],
    'Shc': ['PF00640', 'PF00017'],
}

class SignallingV5Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rules = json.loads((ROOT/'config/family_rules.json').read_text())['families']

    def hits(self, domains):
        return [dict(query_name='p', pfam_accession=p, ali_from=10+i*90, ali_to=75+i*90) for i,p in enumerate(domains)]

    def test_models_available_and_order_distinguishes_adaptors(self):
        models=json.loads((ROOT/'databases/pfam/adhesome_v1.0.hmm.metadata.json').read_text())['models']
        for family, domains in CASES.items():
            with self.subTest(family=family):
                self.assertTrue(set(domains)<=models.keys())
                a=architecture_evidence(self.hits(domains), self.rules[family])
                self.assertTrue(a['complete'] and a['support_order'])
        self.assertFalse(architecture_evidence(self.hits(CASES['Grb2']), self.rules['Nck'])['support_order'])
        self.assertFalse(architecture_evidence(self.hits(CASES['Nck']), self.rules['Grb2'])['support_order'])
        self.assertFalse(architecture_evidence(self.hits(['PF00017','PF07714']),self.rules['Syk'])['complete'])

    def test_new_families_need_reciprocal_evidence_and_order(self):
        for family, domains in CASES.items():
            for mode, reciprocal, expected in [('ordered','true','supported_candidate'),('ordered','false','provisional'),('incomplete','true','provisional'),('wrong_order','true','provisional')]:
                # Grb2's architecture is palindromic; move both SH3s before SH2 instead.
                observed=domains if mode=='ordered' else domains[:1] if mode=='incomplete' else CASES['Nck'] if family=='Grb2' else list(reversed(domains))
                with self.subTest(family=family,mode=mode,reciprocal=reciprocal), tempfile.TemporaryDirectory() as d:
                    p=Path(d)
                    (p/'p.fa').write_text('>p\n'+'A'*3000+'\n')
                    (p/'rules.json').write_text(json.dumps(dict(status='test',families={family:self.rules[family]})))
                    rows=self.hits(observed); write_table(p/'d.csv', rows, list(rows[0]))
                    evidence=dict(protein_id='p',family=family,similarity='true',reciprocal=reciprocal,source='test')
                    write_table(p/'e.csv',[evidence],list(evidence))
                    classify_comparative(argparse.Namespace(proteome=p/'p.fa',rules=p/'rules.json',domains=p/'d.csv',homology=p/'e.csv',evidence=[],manifest=None,out=p/'out'))
                    self.assertEqual(table(p/'out/candidates.csv')[0]['status'],expected)

    def test_domains_alone_do_not_assign_new_families(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d); (p/'p.fa').write_text('>p\n'+'A'*3000+'\n')
            (p/'rules.json').write_text(json.dumps(dict(status='test',families={f:self.rules[f] for f in CASES})))
            rows=self.hits(sum(CASES.values(),[]));write_table(p/'d.csv',rows,list(rows[0]))
            classify_comparative(argparse.Namespace(proteome=p/'p.fa',rules=p/'rules.json',domains=p/'d.csv',homology=None,evidence=[],manifest=None,out=p/'out'))
            self.assertFalse(any(r['assigned_family'] in CASES for r in table(p/'out/candidates.csv')))

    def test_paired_references_and_csk_labels(self):
        refs=ROOT/'references/adhesome_v1.0'
        forward={r['reference_id']:r for r in table(refs/'search_catalogue.csv')}
        reverse={r['reference_id']:r for r in table(refs/'metazoan_catalogue.csv')}
        self.assertEqual(len(forward),71); self.assertEqual(len(reverse),142)
        for prefix, rows in [('search',forward),('metazoan',reverse)]:
            self.assertEqual(set(sequence_dict(refs/(prefix+'_references.faa'))),set(rows))
        for rid,row in forward.items():
            self.assertEqual((row['family'],row['sequence_sha256']),(reverse[rid]['family'],reverse[rid]['sequence_sha256']))
        self.assertEqual(reverse['Hsa_CSK']['family'],'Csk')
        self.assertEqual(reverse['Dme_Shark']['family'],'non_adhesome')
        for family in CASES:
            self.assertTrue(any(r['family']==family for r in forward.values()))
        for name, expected in json.loads((refs/'reference_lock.json').read_text())['files'].items():
            self.assertEqual(digest(refs/name),expected)

if __name__=='__main__':
    unittest.main()
