
# Resolve the organized bundle independently of the caller working directory.
import sys as _bundle_sys
from pathlib import Path as _BundlePath
_BUNDLE_ROOT=next(p for p in _BundlePath(__file__).resolve().parents if (p/'pipeline_paths.py').is_file())
if str(_BUNDLE_ROOT) not in _bundle_sys.path:_bundle_sys.path.insert(0,str(_BUNDLE_ROOT))
from pipeline_paths import activate as _activate_bundle_paths
_activate_bundle_paths()
import argparse
import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from adhesome_core import digest,validate_protein_input,write_table
from adhesome_pipeline import batch,prepare,scan
from parse_adhesome_domtblout import FIELDS


class InPlaceResumeTests(unittest.TestCase):
    def setup_batch(self,root,codes):
        settings=dict(qc=dict(min_length=1,max_ambiguous_fraction=.1,max_ambiguous_run=5),
            proteomes=[dict(code=c,species=c,assembly='fixture',path=c+'.faa',gene_regex=r'gene=(\S+)') for c in codes])
        for code in codes:(root/(code+'.faa')).write_text('>one gene=one\nMAAA\n')
        config=root/'targets.json';config.write_text(json.dumps(settings))
        rules=root/'rules.json';rules.write_text('{}')
        manifest=root/'manifest.csv';manifest.write_text('Pfam accession\nPF00001\n')
        library=str(root/'library.hmm');Path(library+'.metadata.json').write_text('{}')
        return argparse.Namespace(config=config,input_dir=root,out=root/'run',prepare_only=False,library=library,
            rules=str(rules),manifest=str(manifest),crosswalk=None,reference_set=None,classification_mode='comparative',
            cpu=1,max_i=.001,max_c=.001,min_coverage=.35)

    def fake_scan(self,args):
        folder=Path(args.out);folder.mkdir(exist_ok=True)
        for name in ('domains.csv','domains.audit.csv'):write_table(folder/name,[],FIELDS)
        (folder/'combined.domtblout').write_text('# complete fixture\n')
        (folder/'scan.json').write_text(json.dumps(dict(proteome_sha256=digest(args.proteome),
            metadata_sha256=digest(args.library+'.metadata.json'),thresholds=dict(i=.001,c=.001,profile_coverage=.35),
            artifact_checksums={name:digest(folder/name) for name in ('domains.csv','domains.audit.csv','combined.domtblout')})))

    def fake_classify(self,args):
        folder=Path(args.out);folder.mkdir(exist_ok=True)
        (folder/'candidates.csv').write_text('protein_id,family,status\n'+Path(args.proteome).stem+'__one,test,provisional\n')
        (folder/'classification.json').write_text(json.dumps(dict(mode='comparative',
            checksums={str(Path(p).resolve()):digest(p) for p in (args.proteome,args.domains,args.rules,args.manifest)})))

    def tree(self,folder):return {str(p.relative_to(folder)):digest(p) for p in folder.rglob('*') if p.is_file()}

    def test_resume_keeps_completed_and_excluded_targets_without_backups(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);args=self.setup_batch(root,['Done','Bad','Next']);calls=[]
            def interrupted_scan(a):
                code=Path(a.proteome).stem;calls.append(code)
                if code=='Bad':
                    Path(a.out).mkdir(exist_ok=True);(Path(a.out)/'failed.log').write_text('original failed scan')
                    raise ValueError('simulated failed scan')
                self.fake_scan(a)
            with patch('adhesome_pipeline.shutil.which',return_value='fixture/tool'),patch('adhesome_pipeline.scan',side_effect=interrupted_scan),patch('adhesome_pipeline.classify',side_effect=self.fake_classify),contextlib.redirect_stdout(io.StringIO()):
                with self.assertRaisesRegex(ValueError,'failed scan'):batch(args)
                completed=self.tree(args.out/'Done');failed=self.tree(args.out/'Bad')
                args.resume=True;args.exclude_target=['Bad']
                batch(args)
                self.assertEqual(calls,['Done','Bad','Next'])
                self.assertEqual(self.tree(args.out/'Done'),completed)
                self.assertEqual(self.tree(args.out/'Bad'),failed)
                self.assertFalse((args.out/'archived_attempts').exists())
                summary=json.loads((args.out/'batch.json').read_text())
                self.assertEqual(set(summary['species']),{'Done','Next'})
                self.assertEqual(summary['excluded_targets']['Bad']['exclusion_reason'],'user_requested_exclusion')
                # Removing the failed target from targets.json is equivalent.
                settings=json.loads(args.config.read_text());settings['proteomes']=[r for r in settings['proteomes'] if r['code']!='Bad']
                args.config.write_text(json.dumps(settings));args.exclude_target=[]
                batch(args);self.assertEqual(calls,['Done','Bad','Next'])
                (root/'Done.faa').write_text('>one gene=one\nMCCC\n')
                with self.assertRaisesRegex(ValueError,'preserved input'):batch(args)

    def test_resume_after_classification_failure_reuses_preparation_and_scan(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);args=self.setup_batch(root,['TargetA']);attempts=[]
            def interrupted_classify(a):
                attempts.append(1)
                if len(attempts)==1:raise ValueError('classification failed')
                self.fake_classify(a)
            with patch('adhesome_pipeline.shutil.which',return_value='fixture/tool'),patch('adhesome_pipeline.scan',side_effect=self.fake_scan) as scans,patch('adhesome_pipeline.classify',side_effect=interrupted_classify),contextlib.redirect_stdout(io.StringIO()):
                with self.assertRaisesRegex(ValueError,'classification failed'):batch(args)
                prepared=self.tree(args.out/'TargetA/prepared');domains=self.tree(args.out/'TargetA/domains')
                args.resume=True;batch(args)
                self.assertEqual(scans.call_count,1)
                self.assertEqual(self.tree(args.out/'TargetA/prepared'),prepared)
                self.assertEqual(self.tree(args.out/'TargetA/domains'),domains)
                self.assertFalse((args.out/'archived_attempts').exists())

    def test_cannot_remove_a_completed_target_during_resume(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);args=self.setup_batch(root,['Done','Next'])
            with patch('adhesome_pipeline.shutil.which',return_value='fixture/tool'),patch('adhesome_pipeline.scan',side_effect=self.fake_scan),patch('adhesome_pipeline.classify',side_effect=self.fake_classify),contextlib.redirect_stdout(io.StringIO()):
                batch(args);before=self.tree(args.out/'Done');args.resume=True;args.exclude_target=['Done']
                with self.assertRaisesRegex(ValueError,'omits completed targets'):batch(args)
                self.assertEqual(self.tree(args.out/'Done'),before)

    def test_nucleotide_fasta_is_rejected_before_preparation_or_scanning(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);args=self.setup_batch(root,['Bad']);(root/'Bad.faa').write_text('>one gene=one\n'+'ACGT'*100+'\n')
            with self.assertRaisesRegex(ValueError,'nucleotide/transcript FASTA'):batch(args)
            self.assertFalse(args.out.exists())
            with self.assertRaisesRegex(ValueError,'nucleotide/transcript FASTA'):prepare(args.config,root/'prepared')
            self.assertFalse((root/'prepared').exists())
            with patch('adhesome_pipeline.run') as command:
                with self.assertRaisesRegex(ValueError,'nucleotide/transcript FASTA'):
                    scan(argparse.Namespace(proteome=root/'Bad.faa',out=root/'domains'))
                command.assert_not_called()

    def test_genuine_long_protein_gets_a_clear_pre_scan_error(self):
        with tempfile.TemporaryDirectory() as td:
            path=Path(td)/'protein.faa';path.write_text('>long\n'+'MELK'*25001+'\n')
            with self.assertRaisesRegex(ValueError,'100,000 residues'):validate_protein_input(path,hmmer_limit=True)


if __name__=='__main__':unittest.main()
