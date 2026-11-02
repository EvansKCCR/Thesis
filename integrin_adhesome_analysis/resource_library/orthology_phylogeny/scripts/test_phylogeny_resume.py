"""In-place resume tests with simulated external tools; no real phylogeny runs."""

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
import importlib.util
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import phylogeny as p


@unittest.skipUnless(importlib.util.find_spec('Bio'),'Biopython required')
class ResumeTests(unittest.TestCase):
    def fixture(self,root,status='running',checkpoint=True):
        group=_layout_resource(root/'OG1');work=_layout_resource(group/'inference');work.mkdir(parents=True)
        seqs={'a':'ACDE'*20,'b':'ACDF'*20,'c':'FGDE'*20,'d':'FGDF'*20}
        p.fasta(_layout_resource(group/'members.faa'),seqs)
        p.write_tsv(_layout_resource(group/'orthology_edges.tsv'),['candidate_id','orthologue_id'],[])
        p.write_tsv(_layout_resource(group/'duplication_evidence.tsv'),['Genes 1','Genes 2'],[])
        files={x.name:p.digest(x) for x in group.iterdir() if _layout_resource(x).is_file()}
        p.dump(_layout_resource(root/'manifest.json'),dict(groups=[dict(group='OG1',status='ready',files=files)]))
        args=argparse.Namespace(project=root,hog=None,outgroup=[],threads=2,iqtree_threads=2,seed=20260920,
            min_sites=50,max_gap=.9,mafft='mafft',trimal='trimal',iqtree='iqtree2',aligner='linsi',
            resume=True,retry_failed=False,allow_partial=True)
        p.fasta(_layout_resource(work/'aligned.faa'),seqs);p.fasta(_layout_resource(work/'trimmed.faa'),seqs)
        if checkpoint:(_layout_resource(work/'gene_tree.ckp.gz')).write_bytes(b'simulated IQ-TREE checkpoint')
        # Legacy interrupted run: no per-stage hashes and no process ID.
        p.dump(_layout_resource(work/'run.json'),dict(status=status,settings={k:str(v) if isinstance(v,Path) else v for k,v in vars(args).items()},
            input_hashes=files,commands=[['mafft'],['trimal'],['iqtree2']],versions={},failed_stage='trimmed_alignment_qc',
            failure_category='alignment_qc_failure',error='fixture failed QC requires review'))
        calls=[]
        def fake(command,**kwargs):
            if '--version' in command:return argparse.Namespace(stdout='fixture version',stderr='')
            calls.append(command)
            if command[0]=='mafft':kwargs['stdout'].write((_layout_resource(group/'members.faa')).read_text())
            elif command[0]=='trimal':
                _layout_resource(Path(command[command.index('-out')+1])).write_text(_layout_resource(Path(command[command.index('-in')+1])).read_text())
            else:
                self.assertEqual(command[0],'iqtree2')
                self.assertNotIn('-redo',command)
                self.assertEqual(command[command.index('--prefix')+1],str((_layout_resource(work/'gene_tree')).resolve()))
                (_layout_resource(work/'gene_tree.treefile')).write_text('(a:1,b:1,(c:1,d:1)90/99:1);')
            return argparse.Namespace(returncode=0)
        return group,work,args,calls,fake

    def test_legacy_iqtree_checkpoint_reuses_both_alignments_without_backups(self):
        with tempfile.TemporaryDirectory() as folder:
            group,work,args,calls,fake=self.fixture(_layout_resource(Path(folder)))
            before={n:(p.digest(work/n),(work/n).stat().st_mtime_ns) for n in ('aligned.faa','trimmed.faa','gene_tree.ckp.gz')}
            with patch.object(p.shutil,'which',side_effect=lambda x:x),patch.object(p.subprocess,'run',side_effect=fake), \
                 patch.object(p,'active_legacy_processes',return_value=[]):
                p.run(args)
                self.assertEqual([c[0] for c in calls],['iqtree2'])
                p.run(args)  # completed groups must be skipped, never recalculated
                self.assertEqual([c[0] for c in calls],['iqtree2'])
            record=json.loads((_layout_resource(work/'run.json')).read_text())
            self.assertEqual(record['status'],'complete')
            self.assertEqual(record['reused_stages'],['mafft','trimal'])
            self.assertEqual(record['resume_events'][0]['previous_status'],'running')
            self.assertFalse(list(_layout_resource(group).glob('inference_attempt_*')))
            for n,pair in before.items():self.assertEqual((p.digest(work/n),(work/n).stat().st_mtime_ns),pair)

    def test_failed_group_stays_flagged_without_explicit_retry(self):
        with tempfile.TemporaryDirectory() as folder:
            group,work,args,calls,fake=self.fixture(_layout_resource(Path(folder)),status='failed')
            before=(_layout_resource(work/'run.json')).read_bytes()
            with patch.object(p.shutil,'which',side_effect=lambda x:x),patch.object(p.subprocess,'run',side_effect=fake), \
                 patch.object(p,'active_legacy_processes',return_value=[]):p.run(args)
            self.assertEqual(calls,[])
            self.assertEqual((_layout_resource(work/'run.json')).read_bytes(),before)
            self.assertEqual(json.loads((_layout_resource(args.project/'run_summary.json')).read_text())['counts'],{'failed':1})

    def test_explicit_retry_resumes_failed_group_in_same_directory(self):
        with tempfile.TemporaryDirectory() as folder:
            group,work,args,calls,fake=self.fixture(_layout_resource(Path(folder)),status='failed')
            args.retry_failed=True
            with patch.object(p.shutil,'which',side_effect=lambda x:x),patch.object(p.subprocess,'run',side_effect=fake), \
                 patch.object(p,'active_legacy_processes',return_value=[]):p.run(args)
            self.assertEqual([c[0] for c in calls],['iqtree2'])
            self.assertFalse(list(_layout_resource(group).glob('inference_attempt_*')))

    def test_partial_mafft_output_restarts_only_the_incomplete_stage(self):
        with tempfile.TemporaryDirectory() as folder:
            group,work,args,calls,fake=self.fixture(_layout_resource(Path(folder)),checkpoint=False)
            (_layout_resource(work/'aligned.faa')).write_text('>a\nACDE\n')
            (_layout_resource(work/'trimmed.faa')).unlink()
            with patch.object(p.shutil,'which',side_effect=lambda x:x),patch.object(p.subprocess,'run',side_effect=fake), \
                 patch.object(p,'active_legacy_processes',return_value=[]):p.run(args)
            self.assertEqual([c[0] for c in calls],['mafft','trimal','iqtree2'])
            self.assertEqual(json.loads((_layout_resource(work/'run.json')).read_text())['status'],'complete')

    def test_changed_seed_or_checkpoint_hash_is_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            group,work,args,calls,fake=self.fixture(_layout_resource(Path(folder)))
            args.seed+=1
            with patch.object(p.shutil,'which',side_effect=lambda x:x),patch.object(p.subprocess,'run',side_effect=fake), \
                 patch.object(p,'active_legacy_processes',return_value=[]):
                with self.assertRaisesRegex(ValueError,'changed or unrecorded seed'):p.run(args)
                args.seed-=1
                record=json.loads((_layout_resource(work/'run.json')).read_text());record['checkpoints']={'trimmed.faa':'changed'}
                p.dump(_layout_resource(work/'run.json'),record)
                p.run(args)
            self.assertEqual(calls,[])
            self.assertIn('trimmed alignment changed',json.loads((_layout_resource(work/'run.json')).read_text())['error'])

    def test_active_old_process_refuses_resume(self):
        with tempfile.TemporaryDirectory() as folder:
            group,work,args,calls,fake=self.fixture(_layout_resource(Path(folder)))
            with patch.object(p,'active_legacy_processes',return_value=[123]):
                with self.assertRaisesRegex(ValueError,'still active'):p.run(args)
            self.assertEqual(calls,[])


if __name__=='__main__':unittest.main()
