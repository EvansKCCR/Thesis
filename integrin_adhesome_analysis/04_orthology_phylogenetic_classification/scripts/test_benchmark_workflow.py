import argparse
import json
from pathlib import Path
import tempfile
import unittest

from benchmark_workflow import run, measure
from pipeline import write_tsv


class BenchmarkWorkflowTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.root=Path(self.temp.name)
        self.raw=self.root/'reference_raw'
        self.inputs=self.raw/'BENCHMARKS/Input'
        self.refogs=self.raw/'BENCHMARKS/RefOGs'
        self.inputs.mkdir(parents=True);self.refogs.mkdir()
        (self.refogs/'RefOG001.txt').write_text('a1\nb1\n')
        (self.refogs/'RefOG002.txt').write_text('a2\nb2\n')
        self.panel=self.root/'panel';self.panel.mkdir()
        (self.panel/'targets.json').write_text(json.dumps(['TargetA']))
        entries=[]
        for species,genes in [('A',['a1','a2','a0']),('B',['b1','b2','b0'])]:
            source=self.inputs/(species+'.fa')
            source.write_text(''.join('>'+gene+'\nMAAA\n' for gene in genes))
            folder=self.root/'references'/species/'prepared';folder.mkdir(parents=True)
            prepared=folder/(species+'.faa')
            prepared.write_text(''.join('>'+species+'__'+gene+'\nMAAA\n' for gene in genes))
            write_tsv(folder/'identifier_map.tsv',['protein_id','prepared_id','selected'],[
                dict(protein_id=gene,prepared_id=species+'__'+gene,selected=1) for gene in genes])
            entries.append(dict(species=species,fasta=str(prepared),gene_level=False,unit_level='benchmark_protein',source='fixture'))
        (self.panel/'panel.snapshot.json').write_text(json.dumps(dict(proteomes=entries)))
        (self.panel/'workflow_manifest.json').write_text('{}')
        self.results=self.root/'results'
        (self.results/'Phylogenetic_Hierarchical_Orthogroups').mkdir(parents=True)
        self.prediction=self.results/'Phylogenetic_Hierarchical_Orthogroups/N0.tsv'
        self.predict([('A__a1','B__b1'),('A__a2','B__b2'),('A__a0','B__b0')])

    def tearDown(self):
        self.temp.cleanup()

    def predict(self,groups):
        rows=[]
        for i,group in enumerate(groups):
            rows.append(dict(HOG='HOG'+str(i),OG='OG'+str(i),**{'Gene Tree Parent Clade':'N0'},
                A=', '.join(g for g in group if g.startswith('A__')),
                B=', '.join(g for g in group if g.startswith('B__')),
                TargetA='TargetA__candidate'))
        write_tsv(self.prediction,['HOG','OG','Gene Tree Parent Clade','A','B','TargetA'],rows)
        # One target cannot belong to several root groups either.
        text=self.prediction.read_text().replace('\tTargetA__candidate','\t')
        self.prediction.write_text(text)

    def evaluate(self,name='benchmark'):
        return run(argparse.Namespace(results=self.results,panel=self.panel,reference_raw=self.raw,out=self.root/name,protocol='custom'))

    def test_partial_bundle_cannot_claim_original_benchmark_validation(self):
        summary=run(argparse.Namespace(results=self.results,panel=self.panel,reference_raw=self.raw,out=self.root/'benchmark'))
        self.assertEqual(summary['status'],'invalid_benchmark_inputs')
        self.assertIsNone(summary['validated_weight'])
        self.assertEqual(summary['original_dataset_integrity']['observed']['reference_proteomes'],2)

    def test_perfect_cohort_produces_weight_one_and_auditable_mapping(self):
        summary=self.evaluate()
        self.assertEqual(summary['status'],'validated')
        self.assertEqual(summary['validated_weight'],1.0)
        self.assertEqual(summary['reference_proteomes'],2)
        self.assertEqual(summary['expected_reference_proteins'],6)
        self.assertEqual(summary['mapped_reference_proteins'],6)
        self.assertTrue((self.root/'benchmark/identifier_mapping.tsv').exists())
        self.assertIn(str(self.prediction.resolve()),summary['inputs'])

    def test_merged_refogs_reduce_the_actual_supplied_scorer_metric(self):
        self.predict([('A__a1','B__b1','A__a2','B__b2'),('A__a0','B__b0')])
        summary=self.evaluate()
        self.assertEqual(summary['status'],'validated')
        self.assertAlmostEqual(summary['validated_weight'],1/3)
        self.assertEqual(summary['metrics']['recall'],1.0)

    def test_zero_f1_is_a_valid_zero_weight_not_missing_evidence(self):
        self.predict([(name,) for name in ('A__a1','A__a2','A__a0','B__b1','B__b2','B__b0')])
        summary=self.evaluate()
        self.assertEqual(summary['status'],'validated')
        self.assertEqual(summary['validated_weight'],0.0)

    def test_unrelated_target_members_are_not_false_positive_reference_genes(self):
        text=self.prediction.read_text()
        text=text.replace('A__a1\tB__b1\t','A__a1\tB__b1\tTargetA__candidate')
        self.prediction.write_text(text)
        self.assertEqual(self.evaluate()['validated_weight'],1.0)

    def test_missing_or_changed_input_reference_disables_weight(self):
        path=self.root/'references/A/prepared/A.faa'
        path.write_text(path.read_text().replace('>A__a0\nMAAA','>A__a0\nMCCC'))
        summary=self.evaluate()
        self.assertEqual(summary['status'],'incompatible_reference_cohort')
        self.assertIsNone(summary['validated_weight'])

    def test_candidate_only_predictions_cannot_supply_a_benchmark(self):
        self.prediction.unlink()
        summary=self.evaluate()
        self.assertEqual(summary['status'],'invalid_benchmark_inputs')
        self.assertIsNone(summary['validated_weight'])

    def test_duplicate_membership_disables_weight(self):
        self.predict([('A__a1','B__b1'),('A__a1','B__b2'),('A__a0','B__b0')])
        self.assertEqual(self.evaluate()['status'],'invalid_benchmark_inputs')

    def test_missing_ground_truth_is_explicit_and_nonfatal(self):
        summary=run(argparse.Namespace(results=self.results,panel=self.panel,reference_raw=self.root/'empty',out=self.root/'benchmark'))
        self.assertEqual(summary['status'],'unavailable')
        self.assertIsNone(summary['validated_weight'])

    def test_uncertain_members_and_singletons_are_safe(self):
        metrics,details,log=measure({'one':{'a','b','c'},'singleton':{'x'}},
            {'one':{'c'},'singleton':set()},[{'a','b','c'},{'x'}])
        self.assertEqual(metrics['f1'],1.0)
        self.assertEqual(metrics['evaluated_refogs'],1)


if __name__=='__main__':
    unittest.main()
