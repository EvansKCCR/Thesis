"""Synthetic acceptance and version/tree-level tests; no real inference."""
import json
from pathlib import Path
import tempfile
import unittest

from benchmark_configuration import acceptance,quality_policy,prediction_source,labelled_clades


class BenchmarkConfigurationTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.root=Path(self.temp.name)

    def tearDown(self):self.temp.cleanup()

    def table(self,name):
        path=self.root/name;path.parent.mkdir(parents=True,exist_ok=True)
        path.write_text('synthetic table');return path

    def test_reported_low_precision_fails_user_policy_without_changing_metrics(self):
        path,criteria=quality_policy()
        self.assertEqual(criteria,{'precision':.8,'f1':.7})
        metrics=dict(precision=.40038240663204905,recall=.8345526505954374,f1=.5411461870014839)
        original=dict(metrics)
        result=acceptance(metrics,criteria)
        self.assertEqual(result['status'],'rejected');self.assertEqual(metrics,original)
        self.assertEqual(acceptance(dict(precision=.8,f1=.7),criteria)['status'],'accepted')
        self.assertEqual(acceptance(dict(precision=.7999999,f1=.9),criteria)['status'],'rejected')

    def test_policy_rejects_nonfinite_thresholds_and_can_add_explicit_exact_criterion(self):
        path=self.root/'policy.json'
        for value in (True,-.1,1.1,float('nan'),'80%'):
            path.write_text(json.dumps(dict(metric='open_orthobench_weighted_pairwise',minimum_precision=value)))
            with self.assertRaises(ValueError):quality_policy(path)
        path.write_text(json.dumps(dict(metric='open_orthobench_weighted_pairwise',minimum_precision=.8,minimum_exact_fraction=.5)))
        self.assertEqual(quality_policy(path)[1],{'precision':.8,'exact_reconstruction_fraction':.5})

    def test_v315_uses_renamed_root_table_even_when_an_old_n0_file_is_present(self):
        self.table('Phylogenetic_Hierarchical_Orthogroups/N0.tsv')
        expected=self.table('Orthogroups/Orthogroups.tsv')
        (self.root/'Log.txt').write_text('Started OrthoFinder version 3.1.5\n')
        table,meta,inputs=prediction_source(self.root)
        self.assertEqual(table,expected);self.assertEqual(meta['selected_node'],'N0')
        self.assertIn(self.root/'Log.txt',inputs)

    def test_lca_is_selected_by_species_not_by_benchmark_scores(self):
        tree=self.table('Species_Tree/SpeciesTree_rooted_node_labels.txt')
        tree.write_text('(Out:1,(A:0.1,B:0.2)N2:1)N0;')
        expected=self.table('Phylogenetic_Hierarchical_Orthogroups/N2.tsv')
        self.table('Orthogroups/Orthogroups.tsv')
        table,meta,inputs=prediction_source(self.root,'cohort-lca',{'A','B'})
        self.assertEqual(table,expected);self.assertEqual(meta['selected_node'],'N2');self.assertIn(tree,inputs)
        expected.unlink()
        with self.assertRaisesRegex(ValueError,'no candidate-filtered or alternate-clade fallback'):
            prediction_source(self.root,'cohort-lca',{'A','B'})

    def test_missing_or_malformed_species_tree_cannot_choose_a_convenient_clade(self):
        with self.assertRaises(ValueError):prediction_source(self.root,'cohort-lca',{'A','B'})
        tree=self.table('Species_Tree/SpeciesTree_rooted_node_labels.txt')
        for text in ('(A:1,A:2)N0;','(A:1,B:2);','(A:1,B:2)N0;garbage'):
            tree.write_text(text)
            with self.assertRaises(ValueError):labelled_clades(tree)
        tree.write_text('(A:1,B:2)N0;')
        with self.assertRaises(ValueError):prediction_source(self.root,'cohort-lca',{'A','Missing'})


if __name__=='__main__':unittest.main()
