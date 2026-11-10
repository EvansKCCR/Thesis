"""Synthetic annotation recovery tests; no biological search/inference tools."""

# Resolve source resources independently of the caller working directory.
import sys as _layout_sys
from pathlib import Path as _LayoutPath
_LAYOUT_BUNDLE_ROOT = next(p for p in _LayoutPath(__file__).resolve().parents
                          if (p/'config/layout_paths.json').is_file() and (p/'pipeline_paths.py').is_file())
if str(_LAYOUT_BUNDLE_ROOT) not in _layout_sys.path:
    _layout_sys.path.insert(0, str(_LAYOUT_BUNDLE_ROOT))
from pipeline_paths import activate as _layout_activate, resource_path as _layout_resource
_layout_activate()


import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0,str(_layout_resource(Path(__file__)).resolve().parents[2]))
import adhesome_workflow as w


class AnnotationRecoveryTests(unittest.TestCase):
    def fixture(self,root,context=False):
        handoff=_layout_resource(root/'handoff');(_layout_resource(handoff/'families')).mkdir(parents=True)
        text='>Target__one\nACDENPAY\n'
        (_layout_resource(handoff/'families/unassigned_domain_candidate.faa')).write_text(text)
        (_layout_resource(handoff/'external_candidates.faa')).write_text(text)
        (_layout_resource(handoff/'verified_regions.tsv')).write_text('family\tsequence_id\tregion\tstart\tend\n')
        w.seal(handoff,'handoff')
        args=['resume','--handoff',str(handoff),'--out',str(_layout_resource(root/'annotation')),'--annotation-only']
        if context:
            returned=_layout_resource(root/'returned');returned.mkdir()
            for name in ('DeepLoc_adhesome_candidate_list.csv','SignalP_adhesome_candidate_list.txt',
                         'NCBI_CDD_adhesome_hitdata.txt','deeptmhmm_predicted_topologies.3line'):
                (returned/name).write_text('synthetic context fixture\n')
            args+=['--context-dir',str(returned)]
        else:args+=['--allow-missing-context']
        return handoff,_layout_resource(root/'annotation'),args

    def fail(self,command,logs,label=None,**kwargs):
        logs.mkdir(parents=True,exist_ok=True)
        w.dump(_layout_resource(logs/'command_001.json'),dict(command=[str(x) for x in command],status='failed',exit_code=2))
        (_layout_resource(logs/'command_001.log')).write_text('Unknown or duplicate family: unassigned_domain_candidate\n')
        raise ValueError('simulated pre-output annotation failure')

    def test_checkpoint_retry_executes_only_synthetic_annotation_and_preserves_logs(self):
        with tempfile.TemporaryDirectory() as td:
            root=_layout_resource(Path(td));handoff,out,args=self.fixture(root)
            before={p.relative_to(handoff):w.digest(p) for p in _layout_resource(handoff).rglob('*') if _layout_resource(p).is_file()}
            with patch.object(w,'execute',side_effect=self.fail):self.assertEqual(w.main(args),1)
            failed_log=(_layout_resource(out/'logs/command_001.log')).read_bytes()
            with patch.object(w,'execute') as execute:
                self.assertEqual(w.main(args),1)
                execute.assert_not_called()
            self.assertEqual(w.main(args+['--retry-failed']),0)
            self.assertEqual(w.verify(out,'annotation')['handoff'],str(_layout_resource(handoff).resolve()))
            self.assertTrue(json.loads((_layout_resource(out/'motifs/manifest.json')).read_text())['completed'])
            self.assertEqual((_layout_resource(out/'logs/command_001.log')).read_bytes(),failed_log)
            self.assertTrue((_layout_resource(out/'logs/command_002.json')).exists())
            self.assertEqual(before,{p.relative_to(handoff):w.digest(p) for p in _layout_resource(handoff).rglob('*') if _layout_resource(p).is_file()})
            with patch.object(w,'execute') as execute:
                self.assertEqual(w.main(args+['--retry-failed']),1)
                execute.assert_not_called()
            self.assertFalse(any('backup' in p.name or '_attempt_' in p.name for p in _layout_resource(root).rglob('*')))

    def test_legacy_failed_command_adoption_checks_staged_context_before_writing(self):
        with tempfile.TemporaryDirectory() as td:
            root=_layout_resource(Path(td));handoff,out,args=self.fixture(root,context=True)
            with patch.object(w,'execute',side_effect=self.fail):self.assertEqual(w.main(args),1)
            (_layout_resource(out/'annotation_recovery.json')).unlink()  # Release before annotation checkpoints.
            before={p.relative_to(out):w.digest(p) for p in _layout_resource(out).rglob('*') if _layout_resource(p).is_file()}
            source=_layout_resource(root/'returned/DeepLoc_adhesome_candidate_list.csv');original=_layout_resource(source).read_bytes()
            source.write_text('changed prediction\n')
            with patch.object(w,'execute') as execute:
                self.assertEqual(w.main(args+['--retry-failed']),1)
                execute.assert_not_called()
            self.assertEqual(before,{p.relative_to(out):w.digest(p) for p in _layout_resource(out).rglob('*') if _layout_resource(p).is_file()})
            source.write_bytes(original)
            with patch.object(w,'execute') as execute:
                self.assertEqual(w.main(args+['--retry-failed']),0)
                execute.assert_called_once()
            self.assertEqual(w.digest(_layout_resource(out/'logs/command_001.json')),before[_layout_resource(Path('logs/command_001.json'))])
            self.assertEqual((_layout_resource(out/'context'/source.name)).read_bytes(),original)
            self.assertTrue((_layout_resource(out/'workflow_manifest.json')).exists())

    def test_retry_refuses_changed_settings_partial_output_or_unfailed_attempt(self):
        with tempfile.TemporaryDirectory() as td:
            root=_layout_resource(Path(td));handoff,out,args=self.fixture(root)
            with patch.object(w,'execute',side_effect=self.fail):self.assertEqual(w.main(args),1)
            checkpoint=_layout_resource(out/'annotation_recovery.json');saved=_layout_resource(checkpoint).read_bytes()
            with patch.object(w,'execute') as execute:
                self.assertEqual(w.main(args+['--retry-failed','--motif-scope','master']),1)
                state=json.loads(saved);state['status']='running';w.dump(checkpoint,state)
                self.assertEqual(w.main(args+['--retry-failed']),1)
                checkpoint.write_bytes(saved);(_layout_resource(out/'motifs')).mkdir()
                self.assertEqual(w.main(args+['--retry-failed']),1)
                execute.assert_not_called()


if __name__=='__main__':unittest.main()
