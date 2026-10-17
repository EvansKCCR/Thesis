"""Synthetic annotation recovery tests; no biological search/inference tools."""
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import adhesome_workflow as w


class AnnotationRecoveryTests(unittest.TestCase):
    def fixture(self,root,context=False):
        handoff=root/'handoff';(handoff/'families').mkdir(parents=True)
        text='>Target__one\nACDENPAY\n'
        (handoff/'families/unassigned_domain_candidate.faa').write_text(text)
        (handoff/'external_candidates.faa').write_text(text)
        (handoff/'verified_regions.tsv').write_text('family\tsequence_id\tregion\tstart\tend\n')
        w.seal(handoff,'handoff')
        args=['resume','--handoff',str(handoff),'--out',str(root/'annotation'),'--annotation-only']
        if context:
            returned=root/'returned';returned.mkdir()
            for name in ('DeepLoc_adhesome_candidate_list.csv','SignalP_adhesome_candidate_list.txt',
                         'NCBI_CDD_adhesome_hitdata.txt','deeptmhmm_predicted_topologies.3line'):
                (returned/name).write_text('synthetic context fixture\n')
            args+=['--context-dir',str(returned)]
        else:args+=['--allow-missing-context']
        return handoff,root/'annotation',args

    def fail(self,command,logs):
        logs.mkdir(parents=True,exist_ok=True)
        w.dump(logs/'command_001.json',dict(command=[str(x) for x in command],status='failed',exit_code=2))
        (logs/'command_001.log').write_text('Unknown or duplicate family: unassigned_domain_candidate\n')
        raise ValueError('simulated pre-output annotation failure')

    def test_checkpoint_retry_executes_only_synthetic_annotation_and_preserves_logs(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);handoff,out,args=self.fixture(root)
            before={p.relative_to(handoff):w.digest(p) for p in handoff.rglob('*') if p.is_file()}
            with patch.object(w,'execute',side_effect=self.fail):self.assertEqual(w.main(args),1)
            failed_log=(out/'logs/command_001.log').read_bytes()
            with patch.object(w,'execute') as execute:
                self.assertEqual(w.main(args),1)
                execute.assert_not_called()
            self.assertEqual(w.main(args+['--retry-failed']),0)
            self.assertEqual(w.verify(out,'annotation')['handoff'],str(handoff.resolve()))
            self.assertTrue(json.loads((out/'motifs/manifest.json').read_text())['completed'])
            self.assertEqual((out/'logs/command_001.log').read_bytes(),failed_log)
            self.assertTrue((out/'logs/command_002.json').exists())
            self.assertEqual(before,{p.relative_to(handoff):w.digest(p) for p in handoff.rglob('*') if p.is_file()})
            with patch.object(w,'execute') as execute:
                self.assertEqual(w.main(args+['--retry-failed']),1)
                execute.assert_not_called()
            self.assertFalse(any('backup' in p.name or '_attempt_' in p.name for p in root.rglob('*')))

    def test_legacy_failed_command_adoption_checks_staged_context_before_writing(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);handoff,out,args=self.fixture(root,context=True)
            with patch.object(w,'execute',side_effect=self.fail):self.assertEqual(w.main(args),1)
            (out/'annotation_recovery.json').unlink()  # Release before annotation checkpoints.
            before={p.relative_to(out):w.digest(p) for p in out.rglob('*') if p.is_file()}
            source=root/'returned/DeepLoc_adhesome_candidate_list.csv';original=source.read_bytes()
            source.write_text('changed prediction\n')
            with patch.object(w,'execute') as execute:
                self.assertEqual(w.main(args+['--retry-failed']),1)
                execute.assert_not_called()
            self.assertEqual(before,{p.relative_to(out):w.digest(p) for p in out.rglob('*') if p.is_file()})
            source.write_bytes(original)
            with patch.object(w,'execute') as execute:
                self.assertEqual(w.main(args+['--retry-failed']),0)
                execute.assert_called_once()
            self.assertEqual(w.digest(out/'logs/command_001.json'),before[Path('logs/command_001.json')])
            self.assertEqual((out/'context'/source.name).read_bytes(),original)
            self.assertTrue((out/'workflow_manifest.json').exists())

    def test_retry_refuses_changed_settings_partial_output_or_unfailed_attempt(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);handoff,out,args=self.fixture(root)
            with patch.object(w,'execute',side_effect=self.fail):self.assertEqual(w.main(args),1)
            checkpoint=out/'annotation_recovery.json';saved=checkpoint.read_bytes()
            with patch.object(w,'execute') as execute:
                self.assertEqual(w.main(args+['--retry-failed','--motif-scope','master']),1)
                state=json.loads(saved);state['status']='running';w.dump(checkpoint,state)
                self.assertEqual(w.main(args+['--retry-failed']),1)
                checkpoint.write_bytes(saved);(out/'motifs').mkdir()
                self.assertEqual(w.main(args+['--retry-failed']),1)
                execute.assert_not_called()


if __name__=='__main__':unittest.main()
