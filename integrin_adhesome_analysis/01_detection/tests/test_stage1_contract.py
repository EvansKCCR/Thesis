import argparse,json,tempfile,unittest
from pathlib import Path
from adhesome_classification import classify_comparative
from adhesome_core import write_table,table
from classification_contract import general_family,general_rules

ROOT=Path(__file__).resolve().parents[1]


class StageOneContractTests(unittest.TestCase):
    def fixture(self,root,families,evidence,observed):
        (root/'p.fa').write_text('>TargetA__g\n'+'A'*600+'\n')
        config=dict(status='fixture',stage1_general_hypotheses=True,families=families)
        (root/'rules.json').write_text(json.dumps(config))
        hits=[dict(query_name='TargetA__g',pfam_accession=pf,ali_from=10+i*100,ali_to=80+i*100) for i,pf in enumerate(observed)]
        write_table(root/'domains.csv',hits,list(hits[0]))
        write_table(root/'e.csv',evidence,list(evidence[0]))
        classify_comparative(argparse.Namespace(proteome=root/'p.fa',rules=root/'rules.json',domains=root/'domains.csv',
            homology=root/'e.csv',evidence=[],manifest=None,out=root/'out'))
        return table(root/'out/candidates.csv')

    def rule(self,**extra):
        result=dict(module='actin_linking',classification_level='family',domain_alternatives=[{'PF00235':1}],min_length=100,requires_identity_evidence=True)
        result.update(extra);return result

    def evidence(self,family='profilin',**extra):
        result=dict(protein_id='TargetA__g',family=family,similarity='true',reciprocal='true',source='fixture')
        result.update(extra);return result

    def test_stage1_supports_general_family_without_manual_identity_fields(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);rows=self.fixture(root,{'profilin':self.rule()},[self.evidence()],['PF00235'])
            self.assertEqual(rows[0]['status'],'supported_candidate')
            self.assertEqual(rows[0]['topology_check'],'deferred_to_stage_2')
            self.assertTrue((root/'out/Stage_1_output.xlsx').exists())

    def test_partial_and_conflicting_family_evidence_are_distinguished(self):
        with tempfile.TemporaryDirectory() as td:
            rows=self.fixture(Path(td),{'profilin':self.rule(domain_alternatives=[{'PF00235':1,'PF00022':1}])},[self.evidence()],['PF00235'])
            self.assertEqual(rows[0]['status'],'provisional')
        with tempfile.TemporaryDirectory() as td:
            rows=self.fixture(Path(td),{'profilin':self.rule()},[self.evidence(reciprocal='false')],['PF00235'])
            self.assertEqual(rows[0]['status'],'conflicting')

    def test_independently_supported_competing_families_become_provisional(self):
        with tempfile.TemporaryDirectory() as td:
            rows=self.fixture(Path(td),{'family_a':self.rule(),'family_b':self.rule()},[self.evidence('family_a'),self.evidence('family_b')],['PF00235'])
            self.assertTrue(all(r['status']=='provisional' for r in rows))
            self.assertTrue(all(r['reason']=='multiple_supported_family_hypotheses' for r in rows))

    def test_gene_specific_library_entries_fold_to_general_stage1_hypotheses(self):
        config=json.loads((ROOT/'config/family_rules.json').read_text());families=config['families']
        for specific,general in [('PXN','paxillin_like'),('LPXN','paxillin_like'),('SRC','Src'),('SYK','Syk'),('ARP6','actin')]:
            self.assertEqual(general_family(specific,families),general)
        selected=general_rules(config)
        self.assertNotIn('PXN',selected);self.assertIn('PINCH_like',selected)
        self.assertEqual(set(config['stage1_motif_rule_ids']),{'integrin_alpha','integrin_beta','collagen_like','fibronectin_like','fibronectin_III_ligand','vitronectin_like'})


if __name__=='__main__':unittest.main()
