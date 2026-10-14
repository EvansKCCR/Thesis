import argparse
import json
from pathlib import Path
import re
import tempfile
import unittest
from adhesome_core import table, write_table, digest, hmm_records
from adhesome_pipeline import classify
from adhesome_components import load_crosswalk

ROOT=Path(__file__).resolve().parents[1]


class ComponentTests(unittest.TestCase):
    def test_active_library_matches_updated_manifest_and_crosswalk(self):
        manifest=ROOT/'config/adhesome_pfam_manifest.csv'
        crosswalk=ROOT/'config/adhesome_component_crosswalk.csv'
        library=ROOT/'databases/pfam/adhesome_v1.0.hmm'
        meta=json.loads(Path(str(library)+'.metadata.json').read_text())
        wanted={a for r in table(manifest) for a in re.findall(r'PF\d{5}',r['Pfam accession'])}
        self.assertEqual(set(meta['models']),wanted)
        self.assertEqual({r[0]['ACC'][0].split('.')[0] for r in hmm_records(library)},wanted)
        self.assertEqual(meta['manifest_sha256'],digest(manifest))
        self.assertEqual(meta['crosswalk_sha256'],digest(crosswalk))
        self.assertEqual(meta['library_sha256'],digest(library))
        self.assertEqual(meta['pipeline_version'],'1.0')
        for kind in ('ga','fallback'):
            self.assertEqual(meta['split_sha256'][kind],digest(str(library)+'.'+kind+'.hmm'))
        self.assertEqual(len(load_crosswalk(crosswalk,manifest)),230)

    def test_component_overlap_preserves_ambiguity_and_detection_status(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            (root/'p.faa').write_text('>p\n'+'A'*100+'\n')
            (root/'rules.json').write_text(json.dumps(dict(status='fixture',families=dict(
                test=dict(module='actin_linking',min_length=50,domain_alternatives=[{'PF00001':1}])))))
            write_table(root/'domains.csv',[dict(query_name='p',pfam_accession='PF00001',ali_from=1,ali_to=60,decision='retained')],
                        ['query_name','pfam_accession','ali_from','ali_to','decision'])
            write_table(root/'manifest.csv',[{'Pfam accession':'PF00001'}],['Pfam accession'])
            fields=['Official Symbol','Gene ID','Protein name','Swiss-Prot ID','Functional Category','FA','reference_accession','reference_pfams']
            rows=[dict(zip(fields,[s,'1',s,'P001','Actin regulation','Intrinsic Proteins','P001','PF00001'])) for s in ('COMP1','COMP2')]
            write_table(root/'crosswalk.csv',rows,fields)
            args=argparse.Namespace(proteome=root/'p.faa',domains=root/'domains.csv',rules=root/'rules.json',
                manifest=root/'manifest.csv',crosswalk=root/'crosswalk.csv',homology=None,evidence=[],out=root/'out')
            classify(args)
            row=table(root/'out/candidates.csv')[0]
            self.assertEqual(row['status'],'unresolved')
            self.assertEqual(row['assigned_family'],'')
            self.assertEqual(row['crosswalk_possible_components'],'COMP1;COMP2')
            self.assertEqual(row['crosswalk_reference_components'],'')
            self.assertTrue(all(r['reference_similarity']=='False' for r in table(root/'out/component_hypotheses.csv')))
            write_table(root/'evidence.csv',[dict(protein_id='p',family='test',similarity='true',reciprocal='true',reference_accession='P001',source='fixture')],
                        ['protein_id','family','similarity','reciprocal','reference_accession','source'])
            args.homology=root/'evidence.csv'
            classify(args)
            self.assertEqual(table(root/'out/candidates.csv')[0]['crosswalk_reference_components'],'COMP1;COMP2')
            self.assertTrue(all(r['reference_similarity']=='True' for r in table(root/'out/component_hypotheses.csv')))
            rows[0]['reference_pfams']='PF00002'
            write_table(root/'crosswalk.csv',rows,fields)
            with self.assertRaisesRegex(ValueError,'absent from manifest'):
                classify(args)
