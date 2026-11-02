"""Synthetic wiring/provenance checks; never start external inference or real runs."""

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
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0,str(_layout_resource(Path(__file__)).resolve().parents[2]))
import adhesome_workflow as w


class SubsetWiringTests(unittest.TestCase):
    def test_inference_tuning_is_explicit_and_cannot_modify_imported_results(self):
        args=argparse.Namespace(orthofinder_inflation=1.5,orthofinder_search='diamond_ultra_sens',orthofinder_method='msa')
        settings,options=w.orthofinder_configuration(args)
        self.assertEqual(options,['-I','1.5','-S','diamond_ultra_sens','-M','msa'])
        self.assertEqual(w.orthofinder_configuration(argparse.Namespace())[1],[])
        args.results=_layout_resource(Path('existing_results'))
        with self.assertRaisesRegex(ValueError,'cannot change imported results'):w.orthology(args)
        for invalid in (float('nan'),.5,11):
            args.orthofinder_inflation=invalid
            with self.assertRaises(ValueError):w.orthofinder_configuration(args)

    def test_recovery_restores_and_overrides_subset_settings_without_running_stages(self):
        with tempfile.TemporaryDirectory() as td:
            run=_layout_resource(Path(td)/'downstream');run.mkdir()
            selection=_layout_resource(Path(td)/'curated.tsv')
            state=dict(settings=dict(annotation=str(_layout_resource(Path(td)/'annotation')),out=str(run),reference_raw=str(_layout_resource(Path(td)/'raw')),
                benchmark_scope='adhesome',benchmark_selection=str(selection),benchmark_protocol='custom',
                threads=1,analysis_threads=1,iqtree_threads=1),stages={})
            (_layout_resource(run/'recovery.json')).write_text(json.dumps(state))
            args=argparse.Namespace(run=run,threads=None,analysis_threads=None,iqtree_threads=None,retry_benchmark=True)
            with patch.object(w,'downstream') as downstream:
                w.recover(args)
                restored=downstream.call_args.args[0]
                self.assertEqual(restored.benchmark_scope,'adhesome')
                self.assertEqual(restored.benchmark_selection,selection)
                self.assertTrue(restored.retry_benchmark)
                args.benchmark_scope='all'
                w.recover(args)
                self.assertIsNone(downstream.call_args.args[0].benchmark_selection)
                args.benchmark_selection=_layout_resource(Path(td)/'other.tsv')
                w.recover(args)
                self.assertEqual(downstream.call_args.args[0].benchmark_scope,'all')
                self.assertEqual(downstream.call_args.args[0].benchmark_selection,_layout_resource(Path(td)/'other.tsv'))

    def test_public_commands_accept_selection_without_executing_any_stage(self):
        for command, required, handler, scope, selection in (
            ('resume',['--handoff','handoff','--context-dir','context'],'resume','benchmark_scope','benchmark_selection'),
            ('continue',['--annotation','annotation'],'downstream','benchmark_scope','benchmark_selection'),
            ('orthology',['--panel','panel'],'orthology','benchmark_scope','benchmark_selection'),
            ('benchmark',['--orthology','orthology'],'benchmark','scope','selection'),
        ):
            prefix='' if command=='benchmark' else 'benchmark-'
            with patch.object(w,handler) as action:
                result=w.main([command,*required,'--out','output','--'+prefix+'scope','adhesome',
                    '--'+prefix+'selection','curated.tsv','--'+prefix+'level','cohort-lca',
                    '--'+prefix+'quality-policy','policy.json'])
                self.assertEqual(result,0)
                args=action.call_args.args[0]
                self.assertEqual(getattr(args,scope),'adhesome')
                self.assertEqual(getattr(args,selection),_layout_resource(Path('curated.tsv')))
                self.assertEqual(getattr(args,prefix.replace('-','_')+'level'),'cohort-lca')
                self.assertEqual(getattr(args,prefix.replace('-','_')+'quality_policy'),_layout_resource(Path('policy.json')))
        with patch.object(w,'recover') as action:
            self.assertEqual(w.main(['recover','--run','run','--retry-benchmark','--benchmark-scope','adhesome',
                '--benchmark-selection','curated.tsv']),0)
            self.assertTrue(action.call_args.args[0].retry_benchmark)

    def test_adapter_command_keeps_prediction_source_and_passes_scope(self):
        with patch.object(w,'execute') as execute:
            w.benchmark_command(_layout_resource(Path('full_results')),_layout_resource(Path('panel')),_layout_resource(Path('raw')),_layout_resource(Path('out')),_layout_resource(Path('logs')),
                scope='adhesome',selection=_layout_resource(Path('selection.tsv')))
        command=execute.call_args.args[0]
        self.assertEqual(command[1].name,'adhesome_benchmark_subset.py')
        self.assertNotIn('--workflow',command)
        self.assertEqual(command[command.index('--results')+1],_layout_resource(Path('full_results')))
        self.assertEqual(command[command.index('--scope')+1],'adhesome')
        self.assertEqual(command[command.index('--selection')+1],_layout_resource(Path('selection.tsv')))

    def test_original_benchmark_is_default_and_execution_failure_never_switches_to_subset(self):
        with patch.object(w,'execute',side_effect=ValueError('original benchmark failed')) as execute:
            with self.assertRaisesRegex(ValueError,'original benchmark failed'):
                w.benchmark_command(_layout_resource(Path('results')),_layout_resource(Path('panel')),_layout_resource(Path('raw')),_layout_resource(Path('out')),_layout_resource(Path('logs')))
        self.assertEqual(execute.call_count,1)
        command=execute.call_args.args[0]
        self.assertEqual(command[1].name,'benchmark.py')
        self.assertIn('--workflow',command)
        self.assertEqual(command[command.index('--scope')+1],'all')

    def test_scoped_summary_preserves_provenance_and_cannot_supply_global_audit_weight(self):
        with tempfile.TemporaryDirectory() as td:
            root=_layout_resource(Path(td));panel=_layout_resource(root/'panel');panel.mkdir()
            orth=_layout_resource(root/'orthology');orth.mkdir()
            results=_layout_resource(orth/'Results_fixture');results.mkdir()
            prediction=_layout_resource(results/'N0.tsv');prediction.write_text('synthetic full predictions')
            selection=_layout_resource(root/'curated.tsv');selection.write_text('synthetic curated selection')
            (_layout_resource(panel/'workflow_manifest.json')).write_text('{}')
            (_layout_resource(orth/'workflow_manifest.json')).write_text('{}')
            folder=_layout_resource(root/'benchmark');folder.mkdir();(_layout_resource(folder/'workflow_manifest.json')).write_text('{}')
            summary=dict(status='validated',validated_weight=.8,scope='adhesome',
                quality_acceptance=dict(status='accepted'),
                evaluated_scope='adhesome_subset_reference_orthogroup_inference',
                metrics=dict(f1=.8,precision=.8,recall=.8,evaluated_refogs=1),
                family_metrics={'fibronectin_like':dict(f1=.8)},selected_families=['fibronectin_like'],
                diagnostics=[],panel_manifest_sha256=w.digest(_layout_resource(panel/'workflow_manifest.json')),
                prediction_file=str(prediction),inputs={str(selection):w.digest(selection)})
            (_layout_resource(folder/'summary.json')).write_text(json.dumps(summary))
            with patch.object(w,'verify',return_value=dict(orthology_manifest=w.digest(_layout_resource(orth/'workflow_manifest.json')))):
                fields,_=w.benchmark_evidence(folder,orth,dict(results=str(results)),panel)
                self.assertEqual(fields['Orthology_benchmark_status'],'validated')
                self.assertEqual(fields['Orthology_benchmark_weight'],'')
                self.assertEqual(json.loads(fields['Orthology_benchmark_family_weights']),{'fibronectin_like':.8})
                summary['quality_acceptance']['status']='rejected';summary['status']='quality_rejected';summary['validated_weight']=None
                (_layout_resource(folder/'summary.json')).write_text(json.dumps(summary))
                fields,_=w.benchmark_evidence(folder,orth,dict(results=str(results)),panel)
                self.assertEqual(fields['Orthology_benchmark_status'],'quality_rejected')
                self.assertEqual(fields['Orthology_benchmark_weight'],'')
                self.assertNotIn('Orthology_benchmark_family_weights',fields)
                self.assertEqual(fields['Orthology_benchmark_f1'],.8)
                selection.write_text('changed selection')
                fields,_=w.benchmark_evidence(folder,orth,dict(results=str(results)),panel)
                self.assertEqual(fields['Orthology_benchmark_status'],'invalid_provenance')


if __name__=='__main__':
    unittest.main()
