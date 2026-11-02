
# Resolve source resources independently of the caller working directory.
import sys as _layout_sys
from pathlib import Path as _LayoutPath
_LAYOUT_BUNDLE_ROOT = next(p for p in _LayoutPath(__file__).resolve().parents
                          if (p/'config/layout_paths.json').is_file() and (p/'pipeline_paths.py').is_file())
if str(_LAYOUT_BUNDLE_ROOT) not in _layout_sys.path:
    _layout_sys.path.insert(0, str(_LAYOUT_BUNDLE_ROOT))
from pipeline_paths import activate as _layout_activate, resource_path as _layout_resource
_layout_activate()


import argparse
import json
from pathlib import Path
import re
import tempfile
import unittest
from adhesome_core import table, write_table, digest, hmm_records
from adhesome_pipeline import classify
from adhesome_components import load_crosswalk

ROOT=_layout_resource(Path(__file__)).resolve().parents[1]


class ComponentTests(unittest.TestCase):
    def test_active_library_matches_updated_manifest_and_crosswalk(self):
        manifest=_layout_resource(ROOT/'config/adhesome_pfam_manifest.csv')
        crosswalk=_layout_resource(ROOT/'config/adhesome_component_crosswalk.csv')
        library=_layout_resource(ROOT/'databases/pfam/adhesome_v1.0.hmm')
        meta=json.loads(_layout_resource(Path(str(library)+'.metadata.json')).read_text())
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
            root=_layout_resource(Path(td))
            (_layout_resource(root/'p.faa')).write_text('>p\n'+'A'*100+'\n')
            (_layout_resource(root/'rules.json')).write_text(json.dumps(dict(status='fixture',families=dict(
                test=dict(module='actin_linking',min_length=50,domain_alternatives=[{'PF00001':1}])))))
            write_table(_layout_resource(root/'domains.csv'),[dict(query_name='p',pfam_accession='PF00001',ali_from=1,ali_to=60,decision='retained')],
                        ['query_name','pfam_accession','ali_from','ali_to','decision'])
            write_table(_layout_resource(root/'manifest.csv'),[{'Pfam accession':'PF00001'}],['Pfam accession'])
            fields=['Official Symbol','Gene ID','Protein name','Swiss-Prot ID','Functional Category','FA','reference_accession','reference_pfams']
            rows=[dict(zip(fields,[s,'1',s,'P001','Actin regulation','Intrinsic Proteins','P001','PF00001'])) for s in ('COMP1','COMP2')]
            write_table(_layout_resource(root/'crosswalk.csv'),rows,fields)
            args=argparse.Namespace(proteome=_layout_resource(root/'p.faa'),domains=_layout_resource(root/'domains.csv'),rules=_layout_resource(root/'rules.json'),
                manifest=_layout_resource(root/'manifest.csv'),crosswalk=_layout_resource(root/'crosswalk.csv'),homology=None,evidence=[],out=_layout_resource(root/'out'))
            classify(args)
            row=table(_layout_resource(root/'out/candidates.csv'))[0]
            self.assertEqual(row['status'],'unresolved')
            self.assertEqual(row['assigned_family'],'')
            self.assertEqual(row['crosswalk_possible_components'],'COMP1;COMP2')
            self.assertEqual(row['crosswalk_reference_components'],'')
            self.assertTrue(all(r['reference_similarity']=='False' for r in table(_layout_resource(root/'out/component_hypotheses.csv'))))
            write_table(_layout_resource(root/'evidence.csv'),[dict(protein_id='p',family='test',similarity='true',reciprocal='true',reference_accession='P001',source='fixture')],
                        ['protein_id','family','similarity','reciprocal','reference_accession','source'])
            args.homology=_layout_resource(root/'evidence.csv')
            classify(args)
            self.assertEqual(table(_layout_resource(root/'out/candidates.csv'))[0]['crosswalk_reference_components'],'COMP1;COMP2')
            self.assertTrue(all(r['reference_similarity']=='True' for r in table(_layout_resource(root/'out/component_hypotheses.csv'))))
            # A hit to a component without a family rule must still reach handoff.
            write_table(_layout_resource(root/'domains.csv'),[],['query_name','pfam_accession','ali_from','ali_to','decision'])
            write_table(_layout_resource(root/'evidence.csv'),[dict(protein_id='p',family='component_reference',component_symbol='COMP1',
                similarity='true',reciprocal='false',reference_accession='P001',source='fixture')],
                ['protein_id','family','component_symbol','similarity','reciprocal','reference_accession','source'])
            classify(args)
            fallback=table(_layout_resource(root/'out/candidates.csv'))[0]
            self.assertEqual(fallback['family'],'component_reference_candidate')
            self.assertEqual(fallback['status'],'provisional')
            self.assertEqual(fallback['assigned_family'],'')
            self.assertIn('COMP1',fallback['crosswalk_reference_components'])
            rows[0]['reference_pfams']='PF00002'
            write_table(_layout_resource(root/'crosswalk.csv'),rows,fields)
            with self.assertRaisesRegex(ValueError,'absent from manifest'):
                classify(args)
