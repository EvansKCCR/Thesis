
# Resolve the organized bundle independently of the caller working directory.
import sys as _bundle_sys
from pathlib import Path as _BundlePath
_BUNDLE_ROOT=next(p for p in _BundlePath(__file__).resolve().parents if (p/'pipeline_paths.py').is_file())
if str(_BUNDLE_ROOT) not in _bundle_sys.path:_bundle_sys.path.insert(0,str(_BUNDLE_ROOT))
from pipeline_paths import activate as _activate_bundle_paths
_activate_bundle_paths()
import argparse
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from adhesome_core import write_table,table,sequence_dict
from adhesome_pipeline import catalogue
from build_reference_sets import resolve,valid_identity,experimental_annotations,build,UniProtClient

def entry(accession='P1',gene='GENE',taxon=9606,length=10,reviewed=True):
    return dict(primaryAccession=accession,genes=[dict(geneName=dict(value=gene))],organism=dict(taxonId=taxon,scientificName='Homo sapiens'),
        entryType='UniProtKB reviewed (Swiss-Prot)' if reviewed else 'UniProtKB unreviewed (TrEMBL)',
        proteinExistence='1: Evidence at protein level',sequence=dict(value='A'*length,length=length),
        entryAudit=dict(sequenceVersion=1,entryVersion=2))

class ReferenceTests(unittest.TestCase):
    def test_identity_requires_gene_and_taxon(self):
        self.assertTrue(valid_identity(entry(),dict(gene='GENE',taxon_id='9606')))
        self.assertFalse(valid_identity(entry(),dict(gene='OTHER',taxon_id='9606')))
        self.assertFalse(valid_identity(entry(),dict(gene='GENE',taxon_id='7227')))

    def test_pinned_anchor_identity_does_not_require_gene_symbol(self):
        row=dict(gene='',taxon_id='9606',accession='P1',role='phylogenetic_anchor',identity_mode='accession_taxon')
        self.assertTrue(valid_identity(entry(),row))
        self.assertFalse(valid_identity(entry(taxon=7227),row))
        self.assertFalse(valid_identity(entry(accession='P2'),row))
        row['role']='adhesome_seed'
        self.assertFalse(valid_identity(entry(),row))

    def test_mismatch_rejected_not_silently_used(self):
        client=unittest.mock.Mock()
        client.get.return_value=entry(gene='WRONG')
        with self.assertRaisesRegex(ValueError,'Identity mismatch'):
            resolve(dict(accession='P1',gene='GENE',taxon_id='9606'),client)

    def test_ambiguity_requires_policy(self):
        client=unittest.mock.Mock()
        client.get.side_effect=[dict(results=[entry(),entry('P2',length=20)])]
        with self.assertRaisesRegex(ValueError,'Ambiguous'):
            resolve(dict(accession='',gene='GENE',taxon_id='9606'),client)
        client.get.side_effect=[dict(results=[entry(),entry('P2',length=20)]),entry('P2',length=20)]
        result,method,choices=resolve(dict(accession='',gene='GENE',taxon_id='9606'),client,ambiguity_policy='prefer-complete')
        self.assertEqual(result['primaryAccession'],'P2')
        self.assertEqual(set(choices),{'P1','P2'})
        self.assertEqual(method,'ranked_gene_taxon_match')

    def test_experimental_assertions_not_inferred_from_reviewed_status(self):
        value=entry()
        self.assertEqual(experimental_annotations(value),[])
        value['comments']=[dict(texts=[dict(value='Binds actin',evidences=[dict(evidenceCode='ECO:0000269',source='PubMed',id='123')])])]
        self.assertEqual(experimental_annotations(value)[0]['id'],'123')

    def test_offline_does_not_use_network(self):
        with tempfile.TemporaryDirectory() as folder:
            with patch('urllib.request.urlopen') as request:
                with self.assertRaisesRegex(ValueError,'Offline cache missing'):
                    UniProtClient(folder,True).get('https://example.org')
                request.assert_not_called()

    def test_paired_build_and_sequence_identity_validation(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder)
            curated=[dict(reference_id='ref',family='actin',module='actin_linking',species='Homo sapiens',taxon_id='9606',gene='GENE',accession='P1',evidence_source='seed')]
            anchor=dict(curated[0],reference_id='anchor',accession='P3',gene='',role='phylogenetic_anchor',identity_mode='accession_taxon')
            broad=[curated[0],dict(curated[0],reference_id='competitor',family='non_adhesome',accession='P2'),anchor]
            write_table(root/'c.csv',curated,list(curated[0]))
            fields=list(curated[0])+['role','identity_mode']
            write_table(root/'b.csv',broad,fields)
            write_table(root/'anchors.csv',[anchor],fields)
            (root/'rules.json').write_text(json.dumps(dict(families={'actin':{}})))
            args=argparse.Namespace(curated_seed=root/'c.csv',metazoan_seed=root/'b.csv',anchor_seed=root/'anchors.csv',rules=root/'rules.json',out=root/'out',cache=root/'cache',
                offline=True,repair_identity=False,ambiguity_policy='error',strict_reviewed_curated=False)
            with patch.object(UniProtClient,'get',side_effect=lambda url:entry(url.rsplit('/',1)[-1].split('.')[0])):
                build(args)
            self.assertEqual(set(sequence_dict(root/'out/adhesome_references.faa')),{'ref'})
            self.assertEqual(set(sequence_dict(root/'out/search_references.faa')),{'ref','anchor'})
            self.assertEqual(table(root/'out/anchor_catalogue.csv')[0]['reference_tier'],'computational_anchor')
            self.assertEqual(table(root/'out/reference_catalogue.csv')[0]['adhesome_experiment_verified'],'false')
            self.assertEqual(table(root/'out/reference_catalogue.csv')[0]['evidence_status'],'reference_identity_verified')
            catalogue(root/'out/reference_catalogue.csv',root/'out/adhesome_references.faa')
            (root/'out/adhesome_references.faa').write_text('>ref\nCCCC\n')
            with self.assertRaisesRegex(ValueError,'checksum mismatch'):
                catalogue(root/'out/reference_catalogue.csv',root/'out/adhesome_references.faa')

if __name__=='__main__':
    unittest.main()
