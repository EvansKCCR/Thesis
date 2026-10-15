import argparse
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import shutil

from adhesome_classification import classify_comparative
from adhesome_core import digest, table, write_table
from evidence_coverage import search_coverage
from adhesome_classification import architecture_evidence


class EvidenceCoverageTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.folder = self.root/'source/TargetA'
        (self.folder/'prepared').mkdir(parents=True)
        self.proteome = self.folder/'prepared/TargetA.faa'
        self.proteome.write_text('>TargetA__one\n'+'A'*100+'\n>TargetA__two\n'+'A'*100+'\n')
        (self.folder/'domains').mkdir()
        self.domains = self.folder/'domains/domains.csv'
        write_table(self.domains, [dict(query_name='TargetA__'+name, pfam_accession='PF00001',
            ali_from=1, ali_to=50, decision='retained') for name in ('one', 'two')],
            ['query_name','pfam_accession','ali_from','ali_to','decision'])
        self.homology = self.folder/'homology/homology.csv'
        snapshot = self.homology.parent/'reference_snapshot'
        snapshot.mkdir(parents=True)
        write_table(snapshot/'search_catalogue.csv', [dict(reference_id='ref',family='broad')], ['reference_id','family'])
        for name in ('search_references.faa', 'metazoan_references.faa'):
            (snapshot/name).write_text('>ref\nAAAA\n')
        (snapshot/'metazoan_catalogue.csv').write_text('reference_id,family\nref,broad\n')
        write_table(self.homology, [dict(protein_id='TargetA__one',family='broad',similarity='true',reciprocal='true',source='fixture')],
                    ['protein_id','family','similarity','reciprocal','source'])
        self.manifest = self.homology.parent/'homology.json'
        self.manifest.write_text(json.dumps(dict(commands=[['blastp','-subject','/old/location/target.faa'],
            ['blastp','-query','/old/location/target.faa']], checksums={'/old/location/target.faa':digest(self.proteome)},
            snapshot_checksums={p.name:digest(p) for p in snapshot.iterdir()})))
        rule = dict(module='actin_linking',min_length=1,domain_alternatives=[{'PF00001':1}],diagnostic=True)
        self.rules = self.root/'rules.json'
        self.rules.write_text(json.dumps(dict(status='fixture',families=dict(broad=rule,
            specific=dict(rule,parent_family='broad',requires_identity_evidence=True)),family_aliases={})))
        self.args = argparse.Namespace(proteome=self.proteome,domains=self.domains,homology=self.homology,
            rules=self.rules,evidence=[],manifest=None,crosswalk=None,out=self.folder/'candidates')

    def tearDown(self):
        self.temp.cleanup()

    def test_missing_hit_distinguished_from_missing_subfamily_reference(self):
        classify_comparative(self.args)
        rows = {(r['protein_id'],r['family']):r for r in table(self.args.out/'candidates.csv')}
        absent = rows['TargetA__two','broad']
        self.assertEqual(absent['similarity'],'no_qualifying_hit')
        self.assertEqual(absent['similarity_reason'],'completed_search_no_qualifying_hit')
        self.assertEqual(absent['reciprocal'],'not_tested')
        self.assertEqual(absent['reciprocal_reason'],'no_qualifying_forward_hit_for_reciprocal_assessment')
        specific = rows['TargetA__one','specific']
        self.assertEqual(specific['similarity'],'not_tested')
        self.assertEqual(specific['similarity_reason'],'no_family_reference')
        self.assertEqual(specific['related_similarity'],'supported')
        self.assertEqual(specific['related_reciprocal'],'supported')
        self.assertEqual(specific['related_reference_families'],'broad')
        self.assertEqual(specific['assigned_family'],'')
        self.assertEqual(specific['reference_subfamily'],'')
        self.assertEqual(specific['status'],'provisional')
        self.assertEqual(specific['order_check'],'not_applicable')
        coverage = {r['family']:r for r in table(self.args.out/'reference_coverage.csv')}
        self.assertEqual(coverage['specific']['calibration_action'],'add_family_specific_references')
        self.assertEqual(coverage['broad']['no_qualifying_hit'],'1')

    def test_standalone_evidence_does_not_imply_completed_search(self):
        self.manifest.unlink()
        classify_comparative(self.args)
        row = next(r for r in table(self.args.out/'candidates.csv') if r['protein_id']=='TargetA__two' and r['family']=='broad')
        self.assertEqual(row['similarity'],'not_tested')
        self.assertEqual(row['reference_count'],'')
        self.assertEqual(row['similarity_reason'],'search_completion_not_documented')

    def test_multidomain_rule_without_order_is_still_not_configured(self):
        arch = architecture_evidence([], dict(domain_alternatives=[{'PF00001':1,'PF00002':1}]))
        self.assertEqual(arch['order'],'not_configured')
        self.assertEqual(arch['order_reason'],'no_domain_order_rule')

    def test_legacy_alias_is_context_only_for_child_family(self):
        config = json.loads(self.rules.read_text())
        config['families']['canonical_screen'] = dict(config['families']['broad'],screening_only=True)
        config['families']['specific']['parent_family'] = 'canonical_screen'
        config['family_aliases'] = dict(broad='canonical_screen')
        self.rules.write_text(json.dumps(config))
        classify_comparative(self.args)
        row = next(r for r in table(self.args.out/'candidates.csv') if r['family']=='specific' and r['protein_id']=='TargetA__one')
        self.assertEqual(row['related_similarity'],'supported')
        self.assertEqual(row['similarity'],'not_tested')
        self.assertEqual(row['assigned_family'],'')

    def test_changed_target_or_reference_snapshot_refused(self):
        self.proteome.write_text('>changed\nAAA\n')
        with self.assertRaisesRegex(ValueError,'target checksum'):
            search_coverage(self.homology,self.proteome)
        meta = json.loads(self.manifest.read_text())
        meta['checksums']['/old/location/target.faa'] = digest(self.proteome)
        self.manifest.write_text(json.dumps(meta))
        (self.homology.parent/'reference_snapshot/search_catalogue.csv').write_text('changed')
        with self.assertRaisesRegex(ValueError,'snapshot missing or changed'):
            search_coverage(self.homology,self.proteome)

    def test_refresh_reuses_cached_evidence_and_preserves_source(self):
        classify_comparative(self.args)
        source = self.folder.parent
        (source/'batch.json').write_text(json.dumps(dict(species=dict(TargetA=dict(status='completed_with_homology')))))
        original = {str(p.relative_to(source)):digest(p) for p in source.rglob('*') if p.is_file()}
        spec = importlib.util.spec_from_file_location('refresh_detection_evidence', Path(__file__).resolve().parents[1]/'scripts/refresh_detection_evidence.py')
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        args = argparse.Namespace(detection=source,out=self.root/'review',rules=self.rules,manifest=None,crosswalk=None)
        module.refresh(args)
        self.assertEqual(original,{str(p.relative_to(source)):digest(p) for p in source.rglob('*') if p.is_file()})
        self.assertTrue((args.out/'TargetA/candidates/reference_coverage.csv').is_file())
        batch = json.loads((args.out/'batch.json').read_text())
        self.assertTrue(batch['refresh_complete'])
        self.assertEqual(batch['species']['TargetA']['status'],'completed_with_homology')
        with self.assertRaisesRegex(ValueError,'new output directory'):
            module.refresh(args)
        # Explicit panel changes rerun only the paired search, preserving domains.
        panel=self.root/'panel';panel.mkdir()
        (panel/'reference_lock.json').write_text(json.dumps(dict(files={})))
        args.reference_set=panel;args.cpu=2;args.out=self.root/'researched'
        def fake_homology(search):
            self.assertEqual(search.cpu,2)
            self.assertEqual(Path(search.references),panel/'search_references.faa')
            shutil.copytree(self.homology.parent,search.out)
        with patch('adhesome_pipeline.homology',side_effect=fake_homology) as search:
            module.refresh(args)
            search.assert_called_once()
        self.assertEqual(digest(self.domains),digest(args.out/'TargetA/domains/domains.csv'))
        self.assertEqual(original,{str(p.relative_to(source)):digest(p) for p in source.rglob('*') if p.is_file()})
        self.domains.write_text('changed')
        args.out = self.root/'changed'
        with self.assertRaisesRegex(ValueError,'Cached classification input'):
            module.refresh(args)
        self.assertFalse(args.out.exists())


if __name__ == '__main__':
    unittest.main()
