"""Deployment/input routing tests; only temporary, explicitly synthetic fixtures."""
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from pipeline_paths import activate
activate()

import contextlib
import io
import json
import tempfile
import unittest
from argparse import Namespace
from unittest.mock import patch
from input_dependencies import InputRegistry, InputDependencyError
import deployment_validation as deployment
from schisto_adhesome import workflows


class ResolverTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory(prefix='deployment_resolver_test_')
        self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
        self.config=json.loads((ROOT/'config/input_registry.json').read_text())
        self.save()

    def save(self):
        path=self.root/'config/input_registry.json'
        path.parent.mkdir(exist_ok=True)
        path.write_text(json.dumps(self.config))

    def test_raw_aliases_resolve_project_relative(self):
        r=InputRegistry(self.root)
        self.assertEqual(r.path_for('reference_proteomes'),self.root/'input/reference_raw')
        self.assertEqual(r.path_for('target_proteomes'),self.root/'input/target_raw')

    def test_release_code_line_endings_normalize_without_weakening_input_hashes(self):
        from input_dependencies import deployment_source_digest
        first=self.root/'one.py';second=self.root/'two.py'
        first.write_bytes(b'print(1)\r\n');second.write_bytes(b'print(1)\n')
        self.assertEqual(deployment_source_digest(first),deployment_source_digest(second))
        first=self.root/'one.hmm';second=self.root/'two.hmm'
        first.write_bytes(b'SCIENTIFIC\r\n');second.write_bytes(b'SCIENTIFIC\n')
        self.assertNotEqual(deployment_source_digest(first),deployment_source_digest(second))

    def test_missing_mandatory_stops_requesting_stage(self):
        with self.assertRaises(InputDependencyError) as caught:
            InputRegistry(self.root).require('target_raw','target.detection')
        self.assertEqual(caught.exception.event['stage'],'target.detection')
        self.assertEqual(caught.exception.event['decision'],'stop_stage')

    def test_unknown_mandatory_dataset_is_not_substituted(self):
        with self.assertRaisesRegex(InputDependencyError,'unregistered_dataset'):
            InputRegistry(self.root).require('unknown_species','homology.mapping')

    def test_optional_is_skipped_with_reason(self):
        event=InputRegistry(self.root).assess('calibration_fasta','calibration',mandatory=False,execution=True)
        self.assertEqual(event['decision'],'skip_operation')
        self.assertEqual(event['reason'],'path_missing')

    def test_missing_selftest_never_disables_mandatory_exclusion_library(self):
        r=InputRegistry(self.root)
        self.assertEqual(r.assess('calibration_fasta','calibration',False)['decision'],'skip_operation')
        with self.assertRaises(InputDependencyError):r.require('negative_exclusion_fasta','target.audit')

    def test_optional_strict_policy_and_mandatory_safeguard(self):
        self.config['unidentified_input']['optional_missing']='stop';self.save()
        with self.assertRaises(InputDependencyError):
            InputRegistry(self.root).assess('calibration_fasta','calibration',False,execution=True)
        self.config['unidentified_input']['mandatory_missing']='skip';self.save()
        with self.assertRaisesRegex(ValueError,'Mandatory input failures must stop'):
            InputRegistry(self.root)

    def test_no_fallback_to_stale_proteome_outputs(self):
        old=self.root/'reference_proteomes';old.mkdir();(old/'Human.faa').write_text('>synthetic\nMACDE\n')
        with self.assertRaises(InputDependencyError):InputRegistry(self.root).require('reference_proteomes','homology')

    def test_selective_panel_does_not_require_other_proteomes(self):
        self.config['reference_panels']={'human_only':['Homo_sapiens.faa']};self.save()
        raw=self.root/'input/reference_raw';raw.mkdir(parents=True)
        (raw/'Homo_sapiens.faa').write_text('>synthetic\nMACDE\n')
        self.assertEqual(InputRegistry(self.root).select_references(panel='human_only'),[raw/'Homo_sapiens.faa'])
        with self.assertRaises(InputDependencyError):InputRegistry(self.root).select_references(files=['Missing.faa'])

    def test_unknown_panel_and_path_escape_rejected(self):
        r=InputRegistry(self.root)
        with self.assertRaises(ValueError):r.select_references(panel='unmapped')
        with self.assertRaises(ValueError):r.select_references(files=['../unrelated.faa'])

    def test_evaluation_paths_cannot_enter_target_root(self):
        self.config['datasets']['validation_fasta']['path']='input/target_raw/validation.faa';self.save()
        with self.assertRaisesRegex(ValueError,'Evaluation controls'):InputRegistry(self.root)

    def test_ambiguous_alias_rejected(self):
        self.config['datasets']['target_raw']['aliases'].append('reference_proteomes');self.save()
        with self.assertRaisesRegex(ValueError,'Ambiguous'):InputRegistry(self.root)

    def test_file_override_is_explicit_and_format_checked(self):
        p=self.root/'protein.faa';p.write_text('>synthetic\nMACDE\n')
        r=InputRegistry(self.root)
        self.assertEqual(r.require('target_raw','homology',p,allow_file=True),p)
        with self.assertRaises(InputDependencyError):r.require('target_raw','target',p)
        txt=self.root/'protein.txt';txt.write_text('not FASTA')
        with self.assertRaises(InputDependencyError):r.require('target_raw','homology',txt,allow_file=True)


class DeploymentTests(unittest.TestCase):
    def test_clean_state_never_opens_biological_inputs_or_launches_tools(self):
        original=Path.open
        def guarded(path,*args,**kwargs):
            if path.is_relative_to(ROOT/'validation') or path.is_relative_to(ROOT/'input') or path.is_relative_to(ROOT/'inputs'):
                raise AssertionError('Biological input opened during clean-state dry run: '+str(path))
            return original(path,*args,**kwargs)
        with patch.object(Path,'open',guarded),patch('subprocess.run',side_effect=AssertionError('Tool launched')):
            report=deployment.audit(ROOT,clean_state=True,phylogeny=True)
        self.assertEqual(report['structural_readiness'],'pass',report['errors'])
        self.assertFalse(report['execution_ready'])
        self.assertFalse(report['analysis_launched'])
        self.assertTrue(all(not e['available'] for e in report['input_decisions']))

    def test_graph_gate_isolation_audit_before_optional_trees(self):
        nodes=deployment.graph('all',True);index={n['id']:i for i,n in enumerate(nodes)}
        for gate in ('target','validation','calibration'):
            self.assertLess(index[gate+'.orthology'],index[gate+'.audit'])
            self.assertLess(index[gate+'.audit'],index[gate+'.phylogeny'])
            node=next(n for n in nodes if n['id']==gate+'.audit')
            self.assertNotIn(gate+'.phylogeny',node['after'])
            self.assertNotIn(gate+'.benchmark',node['after'])
        target=next(n for n in nodes if n['id']=='target.detection')
        self.assertEqual(target['inputs'],['target_raw','target_configuration'])

    def test_cycles_unknown_edges_and_duplicate_stages_fail(self):
        for nodes in [[{'id':'x','after':['x']}],[{'id':'x','after':['missing']}],[{'id':'x','after':[]},{'id':'x','after':[]}]]:
            with self.assertRaises(ValueError):deployment.ordered(nodes)

    def test_optional_calibration_skips_entire_simulated_gate(self):
        report=deployment.audit(ROOT,'target',True)
        states={n['id']:n['simulated_state'] for n in report['stages']}
        self.assertEqual(states['calibration.detection'],'skipped_optional_input')
        self.assertEqual(states['calibration.report'],'skipped_optional_gate')
        self.assertEqual(states['target.detection'],'blocked_on_data')

    def test_phylogeny_requires_generated_library_not_historical_run(self):
        report=deployment.audit(ROOT,'phylogeny',True)
        self.assertFalse(report['data_execution_readiness'])
        self.assertEqual(report['input_decisions'][0]['requested'],'phylogeny_library')

    def test_reference_only_mode_does_not_request_complete_proteomes(self):
        report=deployment.audit(ROOT,'homology',True,reference_only=True)
        self.assertNotIn('reference_raw',{e['requested'] for e in report['input_decisions']})

    def test_selected_inventory_excludes_unrelated_reference_data(self):
        from schisto_adhesome.datasets import reference_inventory
        with tempfile.TemporaryDirectory(prefix='selected_inventory_test_') as tmp:
            raw=Path(tmp);selected=raw/'Selected.faa';selected.write_text('>synthetic\nMACDE\n')
            (raw/'Other.faa').write_text('UNRELATED DATA')
            (raw/'Selected.protein_gene.tsv').write_text('protein_id\tgene_id\nsynthetic\tg1\n')
            result=reference_inventory(raw,[selected])
            self.assertEqual(set(result['files']),{'Selected.faa','Selected.protein_gene.tsv'})
            self.assertEqual(result['selection'],['Selected.faa'])

    def test_reports_only_written_when_requested(self):
        with tempfile.TemporaryDirectory(prefix='deployment_report_test_') as tmp:
            folder=Path(tmp)/'reports'
            with contextlib.redirect_stdout(io.StringIO()):
                code=deployment.main(['--dry-run','--clean-state','--report-dir',str(folder)])
            self.assertEqual(code,0)
            self.assertEqual({p.name for p in folder.iterdir()},{'deployment_audit.json','deployment_audit.md'})
            self.assertFalse(json.loads((folder/'deployment_audit.json').read_text())['analysis_launched'])

    def test_missing_calibration_does_not_materialize_or_launch_it(self):
        with tempfile.TemporaryDirectory(prefix='optional_selftest_test_') as tmp:
            root=Path(tmp);(root/'config').mkdir();(root/'config/input_registry.json').write_text((ROOT/'config/input_registry.json').read_text())
            raw=root/'input/reference_raw';raw.mkdir(parents=True);(raw/'Reference.faa').write_text('>synthetic_reference gene=synthetic_gene\nMACDE\n')
            val=root/'validation';val.mkdir();(val/'validation.faa').write_text('>synthetic_evaluation\nMACDE\n');(val/'validation_set.csv').write_text('fasta_id\nsynthetic_evaluation\n')
            args=Namespace(config=ROOT/'config/isolated_workflow.json',mode='validation',
                calibration_fasta=val/'calibration.faa',calibration_labels=val/'calibration_set.xlsx',
                validation_fasta=val/'validation.faa',validation_labels=val/'validation_set.csv',
                reference_raw=raw,select_labelled_validation=False,out=root/'run',resume=False)
            with patch.object(workflows,'ROOT',root),patch.object(workflows,'algorithm_inventory',return_value={}),\
                 patch.object(workflows,'load_dataset',return_value=([],{'excluded_unlabelled_ids':[]})),\
                 patch.object(workflows,'materialize') as materialize,patch.object(workflows,'reference_matches',return_value=[]):
                contract=workflows.prepare(args)
            self.assertEqual(set(contract['datasets']),{'validation'})
            self.assertEqual(contract['skipped_gates']['calibration']['status'],'skipped')
            self.assertEqual(materialize.call_args.args[0],'validation')
            self.assertFalse((root/'run/calibration').exists())

    def test_labels_only_control_ids_still_cannot_enter_target_stream(self):
        from schisto_adhesome import datasets
        with tempfile.TemporaryDirectory(prefix='labels_only_gate_test_') as tmp:
            root=Path(tmp);(root/'config').mkdir();(root/'config/input_registry.json').write_text((ROOT/'config/input_registry.json').read_text())
            target=root/'input/target_raw';target.mkdir(parents=True)
            (target/'TargetA.faa').write_text('>CONTROL gene=g1\nMACDEFGHIKLMNPQRSTVWY\n')
            labels=root/'validation';labels.mkdir();(labels/'calibration_set.xlsx').touch()
            config=json.loads((ROOT/'config/targets.example.json').read_text())
            config_path=root/'config/targets.json';config_path.write_text(json.dumps(config))
            args=Namespace(config=ROOT/'config/isolated_workflow.json',mode='target',targets_config=config_path,input_dir=target,
                calibration_fasta=labels/'calibration.faa',calibration_labels=labels/'calibration_set.xlsx',
                validation_fasta=labels/'validation.faa',validation_labels=labels/'validation_set.csv',
                reference_raw=root/'input/reference_raw',out=root/'run',resume=False)
            with patch.object(workflows,'ROOT',root),patch.object(datasets,'label_rows',return_value=[{'accession':'CONTROL'}]):
                with self.assertRaisesRegex(ValueError,'control input'):workflows.prepare(args)
            self.assertFalse((root/'run').exists())


if __name__=='__main__':unittest.main()
