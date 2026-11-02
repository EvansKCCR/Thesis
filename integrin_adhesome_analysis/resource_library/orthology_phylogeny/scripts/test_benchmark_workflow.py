
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
import tempfile
import unittest
from unittest.mock import patch

from benchmark_workflow import run, measure
from pipeline import write_tsv
from prepare_reference import prepare


class BenchmarkWorkflowTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.root=_layout_resource(Path(self.temp.name))
        self.raw=_layout_resource(self.root/'reference_raw')
        self.inputs=_layout_resource(self.raw/'BENCHMARKS/Input')
        self.refogs=_layout_resource(self.raw/'BENCHMARKS/RefOGs')
        self.inputs.mkdir(parents=True);self.refogs.mkdir()
        (_layout_resource(self.refogs/'RefOG001.txt')).write_text('a1\nb1\n')
        (_layout_resource(self.refogs/'RefOG002.txt')).write_text('a2\nb2\n')
        self.panel=_layout_resource(self.root/'panel');self.panel.mkdir()
        (_layout_resource(self.panel/'targets.json')).write_text(json.dumps(['TargetA']))
        entries=[]
        for species,genes in [('A',['a1','a2','a0']),('B',['b1','b2','b0'])]:
            source=self.inputs/(species+'.fa')
            source.write_text(''.join('>'+gene+'\nMAAA\n' for gene in genes))
            folder=_layout_resource(self.root/'references'/species/'prepared');folder.mkdir(parents=True)
            prepared=folder/(species+'.faa')
            prepared.write_text(''.join('>'+species+'__'+gene+'\nMAAA\n' for gene in genes))
            write_tsv(_layout_resource(folder/'identifier_map.tsv'),['protein_id','prepared_id','selected'],[
                dict(protein_id=gene,prepared_id=species+'__'+gene,selected=1) for gene in genes])
            entries.append(dict(species=species,fasta=str(prepared),gene_level=False,unit_level='benchmark_protein',source='fixture'))
        (_layout_resource(self.panel/'panel.snapshot.json')).write_text(json.dumps(dict(proteomes=entries)))
        (_layout_resource(self.panel/'workflow_manifest.json')).write_text('{}')
        self.results=_layout_resource(self.root/'results')
        (_layout_resource(self.results/'Phylogenetic_Hierarchical_Orthogroups')).mkdir(parents=True)
        self.prediction=_layout_resource(self.results/'Phylogenetic_Hierarchical_Orthogroups/N0.tsv')
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
        text=_layout_resource(self.prediction).read_text().replace('\tTargetA__candidate','\t')
        self.prediction.write_text(text)

    def evaluate(self,name='benchmark'):
        return run(argparse.Namespace(results=self.results,panel=self.panel,reference_raw=self.raw,out=self.root/name,protocol='custom'))

    def test_partial_bundle_cannot_claim_original_benchmark_validation(self):
        summary=run(argparse.Namespace(results=self.results,panel=self.panel,reference_raw=self.raw,out=_layout_resource(self.root/'benchmark')))
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
        self.assertTrue((_layout_resource(self.root/'benchmark/identifier_mapping.tsv')).exists())
        self.assertIn(str(_layout_resource(self.prediction).resolve()),summary['inputs'])

    def test_merged_refogs_reduce_the_actual_supplied_scorer_metric(self):
        self.predict([('A__a1','B__b1','A__a2','B__b2'),('A__a0','B__b0')])
        summary=self.evaluate()
        self.assertEqual(summary['status'],'quality_rejected')
        self.assertIsNone(summary['validated_weight'])
        self.assertAlmostEqual(summary['metrics']['f1'],1/3)
        self.assertEqual(summary['metrics']['recall'],1.0)

    def test_zero_f1_is_reported_as_zero_but_does_not_pass_acceptance(self):
        self.predict([(name,) for name in ('A__a1','A__a2','A__a0','B__b1','B__b2','B__b0')])
        summary=self.evaluate()
        self.assertEqual(summary['status'],'quality_rejected')
        self.assertIsNone(summary['validated_weight'])
        self.assertEqual(summary['metrics']['f1'],0.0)

    def test_unrelated_target_members_are_not_false_positive_reference_genes(self):
        text=_layout_resource(self.prediction).read_text()
        text=text.replace('A__a1\tB__b1\t','A__a1\tB__b1\tTargetA__candidate')
        self.prediction.write_text(text)
        self.assertEqual(self.evaluate()['validated_weight'],1.0)

    def test_missing_or_changed_input_reference_disables_weight(self):
        path=_layout_resource(self.root/'references/A/prepared/A.faa')
        path.write_text(_layout_resource(path).read_text().replace('>A__a0\nMAAA','>A__a0\nMCCC'))
        summary=self.evaluate()
        self.assertEqual(summary['status'],'incompatible_reference_cohort')
        self.assertIsNone(summary['validated_weight'])

    def test_explicit_stop_masking_matches_full_benchmark_cohort_and_reports_changed_truth(self):
        source=_layout_resource(self.inputs/'A.fa')
        source.write_text(_layout_resource(source).read_text().replace('>a1\nMAAA','>a1\nMA*A*'))
        mapping=_layout_resource(self.root/'map.tsv')
        mapping.write_text('protein_id\tgene_id\na1\ta1\na2\ta2\na0\ta0\n')
        prepared=prepare('A',source,mapping,_layout_resource(self.root/'masked'),internal_stop_policy='mask')
        snapshot=json.loads((_layout_resource(self.panel/'panel.snapshot.json')).read_text())
        snapshot['proteomes'][0]['fasta']=str(prepared)
        (_layout_resource(self.panel/'panel.snapshot.json')).write_text(json.dumps(snapshot))
        summary=self.evaluate()
        self.assertEqual(summary['status'],'validated',summary['diagnostics'])
        self.assertEqual(summary['validated_weight'],1.0)
        self.assertEqual(summary['sequence_qc']['affected_refog_proteins'],['a1'])
        self.assertEqual(summary['sequence_qc']['internal_stop_sites'],1)
        self.assertIn(str((_layout_resource(prepared.parent/'sequence_qc.tsv')).resolve()),summary['inputs'])
        self.assertIn('processed cohort',' '.join(summary['diagnostics']))
        prepared.write_text(_layout_resource(prepared).read_text().replace('MAXA','MAYA'))
        changed=self.evaluate('changed_benchmark')
        self.assertEqual(changed['status'],'incompatible_reference_cohort')
        self.assertIsNone(changed['validated_weight'])

    def test_candidate_only_predictions_cannot_supply_a_benchmark(self):
        self.prediction.unlink()
        summary=self.evaluate()
        self.assertEqual(summary['status'],'invalid_benchmark_inputs')
        self.assertIsNone(summary['validated_weight'])

    def test_duplicate_membership_disables_weight(self):
        self.predict([('A__a1','B__b1'),('A__a1','B__b2'),('A__a0','B__b0')])
        self.assertEqual(self.evaluate()['status'],'invalid_benchmark_inputs')

    def test_missing_ground_truth_is_explicit_and_nonfatal(self):
        summary=run(argparse.Namespace(results=self.results,panel=self.panel,reference_raw=_layout_resource(self.root/'empty'),out=_layout_resource(self.root/'benchmark')))
        self.assertEqual(summary['status'],'unavailable')
        self.assertIsNone(summary['validated_weight'])

    def test_uncertain_members_and_singletons_are_safe(self):
        metrics,details,log=measure({'one':{'a','b','c'},'singleton':{'x'}},
            {'one':{'c'},'singleton':set()},[{'a','b','c'},{'x'}])
        self.assertEqual(metrics['f1'],1.0)
        self.assertEqual(metrics['evaluated_refogs'],1)
        self.assertEqual(metrics['exact_reconstructed_refogs'],1)
        self.assertEqual(metrics['exact_reconstruction_fraction'],1)

    def subset(self, rows=None, name='subset', protocol='custom'):
        selection=_layout_resource(self.raw/'adhesome_refog_selection.tsv')
        write_tsv(selection,['refog','family','curation_source','curation_status'],rows if rows is not None else [
            dict(refog='RefOG001',family='profilin',curation_source='independent synthetic fixture',curation_status='curated')])
        return run(argparse.Namespace(results=self.results,panel=self.panel,reference_raw=self.raw,
            out=self.root/name,protocol=protocol,scope='adhesome'))

    def test_subset_ignores_unselected_truth_performance_without_deleting_prediction_members(self):
        self.predict([('A__a1','B__b1'),('A__a2',),('B__b2',),('A__a0','B__b0')])
        summary=self.subset()
        self.assertEqual(summary['status'],'validated',summary['diagnostics'])
        self.assertEqual(summary['reference_orthogroups'],2)
        self.assertEqual(summary['metrics']['evaluated_refogs'],1)
        self.assertEqual(summary['metrics']['ppv'],1)
        self.assertEqual(summary['family_metrics']['profilin']['f1'],1)
        self.assertEqual(summary['selected_reference_orthogroups'],1)
        self.assertIn(str((_layout_resource(self.raw/'adhesome_refog_selection.tsv')).resolve()),summary['inputs'])
        self.assertTrue((_layout_resource(self.root/'subset/family_coverage.tsv')).is_file())

    def test_unselected_and_background_members_still_penalize_subset_precision(self):
        for suffix, extra in [('other_refog',('A__a2','B__b2')),('background',('A__a0','B__b0'))]:
            self.predict([('A__a1','B__b1',*extra),tuple(x for x in ('A__a0','B__b0','A__a2','B__b2') if x not in extra)])
            summary=self.subset(name=suffix)
            self.assertAlmostEqual(summary['metrics']['precision'],.2)
            self.assertAlmostEqual(summary['metrics']['f1'],1/3)
            self.assertEqual(summary['metrics']['recall'],1)
            self.assertEqual(summary['family_metrics']['profilin']['f1'],summary['metrics']['f1'])
            import csv
            with (_layout_resource(self.root/suffix/'refog_scores.tsv')).open() as handle:
                row=next(csv.DictReader(handle,delimiter='\t'))
            self.assertAlmostEqual(float(row['ppv']),.2)
            self.assertAlmostEqual(float(row['f1']),1/3)

    def test_multiple_refogs_in_one_family_do_not_merge_ground_truth(self):
        self.predict([('A__a1','B__b1','A__a2','B__b2'),('A__a0','B__b0')])
        summary=self.subset([dict(refog=name,family='profilin',curation_source='synthetic fixture',curation_status='curated')
            for name in ('RefOG001','RefOG002')])
        self.assertEqual(summary['family_metrics']['profilin']['evaluated_refogs'],2)
        self.assertAlmostEqual(summary['family_metrics']['profilin']['f1'],1/3)

    def test_fibronectin_uses_curated_selection_and_other_fibronectin_scopes_remain_untested(self):
        summary=self.subset([dict(refog='RefOG001',family='fibronectin_like',curation_source='synthetic curation only',curation_status='curated')])
        self.assertEqual(set(summary['family_metrics']),{'fibronectin_like'})
        with (_layout_resource(self.root/'subset/family_coverage.tsv')).open() as handle:
            import csv
            coverage={row['family']:row for row in csv.DictReader(handle,delimiter='\t')}
        self.assertEqual(coverage['fibronectin_like']['status'],'evaluated')
        self.assertEqual(coverage['fibronectin_III_ligand']['status'],'no_curated_ground_truth')
        self.assertEqual(coverage['fibronectin_like_candidate']['f1'],'')

    def test_invalid_or_empty_curated_selections_never_fallback_to_all_refogs(self):
        base=dict(refog='RefOG001',family='profilin',curation_source='synthetic fixture',curation_status='curated')
        cases=[[],[dict(base,refog='RefOG999')],[dict(base,family='not_a_family')],
               [dict(base,curation_source='')],[dict(base,curation_status='provisional')],[base,base]]
        for index,rows in enumerate(cases):
            summary=self.subset(rows,name='invalid_'+str(index))
            self.assertEqual(summary['status'],'invalid_benchmark_inputs')
            self.assertIsNone(summary['validated_weight'])
            self.assertNotIn('metrics',summary)

    def test_missing_selection_and_original_integrity_cannot_be_bypassed_by_subset_scope(self):
        missing=run(argparse.Namespace(results=self.results,panel=self.panel,reference_raw=self.raw,
            out=_layout_resource(self.root/'missing_selection'),protocol='custom',scope='adhesome'))
        self.assertEqual(missing['status'],'invalid_benchmark_inputs')
        summary=self.subset(protocol='open_orthobench')
        self.assertEqual(summary['status'],'invalid_benchmark_inputs')
        self.assertIsNone(summary['validated_weight'])
        self.assertEqual(summary['original_dataset_integrity']['observed']['reference_orthogroups'],2)

    def test_selected_uncertain_members_are_excluded_only_for_their_refog(self):
        folder=_layout_resource(self.refogs/'low_certainty_assignments');folder.mkdir()
        (_layout_resource(folder/'RefOG001.txt')).write_text('a0\n')
        self.predict([('A__a1','B__b1','A__a0'),('A__a2','B__b2'),('B__b0',)])
        self.assertEqual(self.subset()['metrics']['f1'],1)

    def test_explicit_selection_is_honoured_and_cannot_be_used_with_all_scope(self):
        selection=_layout_resource(self.root/'explicit.tsv')
        write_tsv(selection,['refog','family','curation_source','curation_status'],[
            dict(refog='RefOG002',family='profilin',curation_source='synthetic fixture',curation_status='curated')])
        args=argparse.Namespace(results=self.results,panel=self.panel,reference_raw=self.raw,
            out=_layout_resource(self.root/'explicit'),protocol='custom',scope='adhesome',selection=selection)
        summary=run(args)
        self.assertEqual(summary['selection_file'],str(_layout_resource(selection).resolve()))
        self.assertEqual(summary['family_metrics']['profilin']['f1'],1)
        args.scope='all';args.out=_layout_resource(self.root/'invalid_all')
        self.assertEqual(run(args)['status'],'invalid_benchmark_inputs')

    def test_separate_subset_entrypoint_defaults_to_adhesome_scope(self):
        from adhesome_benchmark_subset import main
        with patch('benchmark_workflow.run',return_value={'status':'synthetic'}) as execute:
            self.assertEqual(main(['--results','results','--panel','panel','--reference-raw','raw','--out','out']),0)
        self.assertEqual(execute.call_args.args[0].scope,'adhesome')

    def test_separate_subset_entrypoint_cannot_silently_run_the_original_scope(self):
        import contextlib
        import io
        from adhesome_benchmark_subset import main
        with patch('benchmark_workflow.run') as execute,contextlib.redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit) as error:
                main(['--results','results','--panel','panel','--reference-raw','raw','--out','out','--scope','all'])
        self.assertEqual(error.exception.code,2)
        execute.assert_not_called()

    def test_cohort_lca_preserves_all_benchmark_members_and_declares_the_tree_level(self):
        folder=_layout_resource(self.results/'Species_Tree');folder.mkdir()
        (_layout_resource(folder/'SpeciesTree_rooted_node_labels.txt')).write_text('(TargetA:1,(A:1,B:1)N2:1)N0;')
        (_layout_resource(self.results/'Phylogenetic_Hierarchical_Orthogroups/N2.tsv')).write_text(_layout_resource(self.prediction).read_text())
        args=argparse.Namespace(results=self.results,panel=self.panel,reference_raw=self.raw,
            out=_layout_resource(self.root/'lca'),protocol='custom',level='cohort-lca')
        summary=run(args)
        self.assertEqual(summary['status'],'validated',summary['diagnostics'])
        self.assertEqual(summary['prediction_configuration']['selected_node'],'N2')
        self.assertEqual(summary['quality_acceptance']['status'],'accepted')
        self.assertIn(str((_layout_resource(folder/'SpeciesTree_rooted_node_labels.txt')).resolve()),summary['inputs'])
        self.assertTrue((_layout_resource(args.out/'orthogroup_error_matrix.tsv')).exists())
        self.assertTrue((_layout_resource(args.out/'acceptance_checks.tsv')).exists())


if __name__=='__main__':
    unittest.main()
