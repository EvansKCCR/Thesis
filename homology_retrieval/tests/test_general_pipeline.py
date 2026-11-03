"""General discovery/routing/checkpoint regressions, with no real BLAST calls."""
import contextlib
import csv
import gzip
import importlib.util
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

source=Path(__file__).resolve().parents[1]/'reciprocal_blast.py'
spec=importlib.util.spec_from_file_location('general_homology_tested',source)
p=importlib.util.module_from_spec(spec);sys.modules[spec.name]=p;spec.loader.exec_module(p)
A='MACDEFGHIKLMNPQRSTVWY'
B='MAEDEFGHIKLMNPQRSTVWY'
OTHER='MYYYYYYYYYYYYYYYYYYY'


class GeneralPipelineTests(unittest.TestCase):
    def test_selected_family_does_not_read_unrelated_reference_proteome(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);args,tool=self.fixture(root)
            (root/'reference_proteomes/Mus_musculus/primary.faa').write_text('UNRELATED INVALID INPUT; MUST NOT BE READ')
            status,text,_=self.invoke(root,args,tool,['--check-inputs','--reference-set','kinase'])
            self.assertEqual(status,0)
            report=json.loads(text)
            self.assertEqual(len(report['reference_sets']),1)
            self.assertEqual([r['reference_proteome'] for r in report['reference_proteomes']],['Homo_sapiens__primary'])

    def test_unlabelled_query_does_not_use_unrelated_sole_background(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);args,tool=self.fixture(root)
            (root/'reference_fasta/kinase.faa').write_text('>sp|unknown\nMFFFFFGHIKLMNPQRSTVWY\n')
            status,text,_=self.invoke(root,args,tool,['--check-inputs','--reference-set','kinase',
                                                     '--reference-proteome','Homo_sapiens__primary','--augment-reference-background'])
            self.assertEqual(status,1)
            self.assertEqual(json.loads(text)['input_decisions'][0]['decision'],'stop_stage')

    def fixture(self, root, targets=2):
        refs,bgs,tgts=root/'reference_fasta',root/'reference_proteomes',root/'target_proteomes'
        (refs/'nested').mkdir(parents=True)
        (bgs/'Homo_sapiens').mkdir(parents=True)
        (bgs/'Mus_musculus').mkdir(parents=True)
        (refs/'kinase.faa').write_text('>homo_sapiens|human_reference\n'+A+'\n')
        (refs/'nested/scaffold.faa').write_text('>mus_musculus|mouse_reference\n'+B+'\n')
        (bgs/'Homo_sapiens/primary.faa').write_text('>native\n'+A+'\n>other\n'+OTHER+'\n')
        (bgs/'Mus_musculus/primary.faa').write_text('>native\n'+B+'\n>other\n'+OTHER+'\n')
        for index in range(targets):
            folder=tgts/f'clade{index}';folder.mkdir(parents=True)
            (folder/'proteome.faa').write_text('>candidate\n'+A+'\n')
        tool=root/'mock_executable';tool.write_text('synthetic executable')
        args=['--reference-dir',str(refs),'--background-dir',str(bgs),
              '--target-dir',str(tgts),'--outdir',str(root/'out')]
        return args,tool

    def invoke(self, root, args, tool, extra=(), runner=None):
        out,err=io.StringIO(),io.StringIO()
        with patch.object(p.shutil,'which',return_value=str(tool)),\
             patch.object(p,'run',side_effect=runner or AssertionError('BLAST unexpectedly called')),\
             contextlib.redirect_stdout(out),contextlib.redirect_stderr(err):
            status=p.main([*args,*extra])
        return status,out.getvalue(),err.getvalue()

    def runner(self, calls, no_hits=False, fail_once=False):
        failed=False
        def run(command, log, stdout_path=None):
            nonlocal failed
            calls.append(command)
            log.write('MOCK '+str(command)+'\n')
            if stdout_path is None:
                prefix=Path(command[command.index('-out')+1])
                for suffix in ('pin','phr','psq'):Path(str(prefix)+'.'+suffix).write_text('mock database')
                return
            if no_hits:
                stdout_path.write_text('');return
            def hit(query,subject,score):
                return f'{query}\t{subject}\t90\t20\t20\t20\t1\t20\t1\t20\t1e-20\t{score}\n'
            if '.forward.' in stdout_path.name:
                stdout_path.write_text(hit('R0000001','P0000001',100)+hit('R0000002','P0000001',90))
            else:
                if fail_once and not failed:
                    failed=True;raise RuntimeError('synthetic interrupted reverse search')
                label=stdout_path.parent.name
                bg='B001_0000001' if label.startswith('Homo') else 'B002_0000001'
                stdout_path.write_text(hit('P0000001',bg,100))
        return run

    def test_recursive_inputs_arbitrary_target_count_and_same_basename(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);args,tool=self.fixture(root,targets=3)
            status,text,_=self.invoke(root,args,tool,['--check-inputs'])
            report=json.loads(text)
            self.assertEqual(status,0)
            self.assertEqual(len(report['target_proteomes']),3)
            self.assertEqual(len(report['reference_sets']),2)
            self.assertEqual(len({r['species'] for r in report['target_proteomes']}),3)
            self.assertFalse((root/'out').exists())

    def test_pooled_forward_searches_and_per_source_reverse_databases(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);args,tool=self.fixture(root);calls=[]
            self.assertEqual(self.invoke(root,args,tool,runner=self.runner(calls))[0],0)
            searches=[cmd for cmd in calls if cmd[0]=='blastp']
            self.assertEqual(sum(Path(cmd[cmd.index('-query')+1]).name=='curated_references.faa' for cmd in searches),2)
            self.assertEqual(len(searches),6)
            self.assertEqual(sum(cmd[0]=='makeblastdb' for cmd in calls),4)
            with (root/'out/homology_audit.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
            self.assertEqual(len(rows),4)
            self.assertTrue(all(r['status']=='RBH_candidate' for r in rows))
            for row in rows:
                self.assertEqual(row['reference_proteome'].split('__')[0].lower(),row['reference_id'].split('|')[0])
            manifest=json.loads((root/'out/run_manifest.json').read_text())
            self.assertEqual(manifest['summary']['unique_rbh_candidates'],2)
            ids=[line for line in (root/'out/rbh_candidates.faa').read_text().splitlines() if line.startswith('>')]
            self.assertEqual(len(ids),2)
            self.assertEqual(len(set(ids)),2)

    def test_per_reference_set_exports_keep_all_matching_queries(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);args,tool=self.fixture(root)
            self.invoke(root,args,tool,runner=self.runner([]))
            for label in ('kinase','nested__scaffold'):
                with (root/f'out/reference_sets/{label}/homology_audit.tsv').open() as f:
                    rows=list(csv.DictReader(f,delimiter='\t'))
                self.assertEqual(len(rows),2)
                self.assertTrue(all(r['reference_set']==label for r in rows))
                self.assertEqual((root/f'out/reference_sets/{label}/rbh_candidates.faa').read_text().count('>'),2)

    def test_completed_resume_executes_no_tools(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);args,tool=self.fixture(root)
            self.invoke(root,args,tool,runner=self.runner([]))
            self.assertEqual(self.invoke(root,args,tool,['--resume'])[0],0)

    def test_failed_resume_reuses_databases_and_forward_search(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);args,tool=self.fixture(root);calls=[]
            runner=self.runner(calls,fail_once=True)
            with self.assertRaisesRegex(RuntimeError,'synthetic interrupted'):
                self.invoke(root,args,tool,runner=runner)
            first_forward=sum(cmd[0]=='blastp' and 'curated_references.faa' in cmd[2] for cmd in calls)
            self.assertEqual(first_forward,1)
            self.assertEqual(self.invoke(root,args,tool,['--resume'],runner=runner)[0],0)
            self.assertEqual(sum(cmd[0]=='makeblastdb' for cmd in calls),4)
            self.assertEqual(sum(cmd[0]=='blastp' and 'curated_references.faa' in cmd[2] for cmd in calls),2)

    def test_changed_input_refuses_resume(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);args,tool=self.fixture(root)
            self.invoke(root,args,tool,runner=self.runner([]))
            path=root/'target_proteomes/clade0/proteome.faa';path.write_text(path.read_text()+'\n')
            with self.assertRaisesRegex(ValueError,'saved inputs/code/settings/tools'):
                self.invoke(root,args,tool,['--resume'])

    def test_no_hits_are_valid_and_cached(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);args,tool=self.fixture(root)
            self.assertEqual(self.invoke(root,args,tool,runner=self.runner([],no_hits=True))[0],0)
            with (root/'out/homology_audit.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
            self.assertEqual(len(rows),4)
            self.assertTrue(all(row['status']=='no_forward_hit' for row in rows))
            self.assertEqual((root/'out/rbh_candidates.faa').read_text(),'')
            self.assertEqual(self.invoke(root,args,tool,['--resume'])[0],0)

    def test_unique_exact_sequence_can_identify_unlabelled_source(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);args,tool=self.fixture(root)
            (root/'reference_fasta/kinase.faa').write_text('>sp|unlabelled\n'+A+'\n')
            status,text,_=self.invoke(root,args,tool,['--check-inputs'])
            report=json.loads(text);self.assertEqual(status,0)
            row=next(x for x in report['reference_mappings'] if x['reference_id']=='sp|unlabelled')
            self.assertEqual(row['source_assignment'],'unique_exact_sequence_source')

    def test_ambiguous_source_never_guessed(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);args,tool=self.fixture(root)
            (root/'reference_fasta/kinase.faa').write_text('>sp|unlabelled\n'+OTHER+'\n')
            status,text,_=self.invoke(root,args,tool,['--check-inputs'])
            self.assertEqual(status,1)
            self.assertTrue(json.loads(text)['errors'])
            self.assertFalse((root/'out').exists())

    def test_missing_explicit_species_cannot_route_to_another_species(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);args,tool=self.fixture(root)
            (root/'reference_fasta/kinase.faa').write_text('>gallus_gallus|missing_source\n'+A+'\n')
            status,text,_=self.invoke(root,args,tool,['--check-inputs'])
            self.assertEqual(status,1)
            self.assertIn('no matching proteome',str(json.loads(text)['errors']))

    def test_explicit_mapping_resolves_unlabelled_missing_query(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);args,tool=self.fixture(root)
            (root/'reference_fasta/kinase.faa').write_text('>candidate|unlabelled\nMWWWWWWWWWWWWWWWWWWW\n')
            mapping=root/'sources.tsv'
            mapping.write_text('reference_file\treference_id\treference_proteome\nkinase.faa\tcandidate|unlabelled\tHomo_sapiens/primary.faa\n')
            status,text,_=self.invoke(root,args,tool,['--source-map',str(mapping),'--augment-reference-background','--check-inputs'])
            self.assertEqual(status,0)
            report=json.loads(text)
            self.assertEqual(report['curated_sequences_added'],1)
            row=next(x for x in report['reference_mappings'] if x['reference_id']=='candidate|unlabelled')
            self.assertEqual(row['source_assignment'],'explicit_source_map')

    def test_gzipped_recursive_fasta_is_supported(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);folder=root/'nested';folder.mkdir()
            with gzip.open(folder/'species.faa.gz','wt') as f:f.write('>x\n'+A+'\n')
            base,files=p.collection(root)
            self.assertIn('nested__species',files)
            self.assertEqual(p.read_fasta(files['nested__species'],'P')[0].seq,A)

    def test_candidate_alias_is_human_for_all_queries(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);args,tool=self.fixture(root)
            (root/'reference_fasta/kinase.faa').write_text('>candidate|human_query\n'+A+'\n')
            (root/'reference_fasta/nested/scaffold.faa').write_text('>candidate|missing_human_query\nMWWWWWWWWWWWWWWWWWWW\n')
            status,text,_=self.invoke(root,args,tool,['--augment-reference-background','--check-inputs'])
            self.assertEqual(status,0)
            rows=json.loads(text)['reference_mappings']
            self.assertTrue(all(row['reference_proteome'].startswith('Homo_sapiens') for row in rows))
            self.assertTrue(all(row['source_assignment']=='user_confirmed_header_alias:candidate=homo_sapiens' for row in rows))
            self.assertEqual({row['reference_id'] for row in rows},{'candidate|human_query','candidate|missing_human_query'})

    def test_output_cannot_be_nested_under_recursive_input(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);args,tool=self.fixture(root)
            with self.assertRaisesRegex(ValueError,'outside recursively scanned input'):
                self.invoke(root,args,tool,['--outdir',str(root/'target_proteomes/results'),'--check-inputs'])

    def test_curated_reference_only_is_explicit_and_provisional(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);args,tool=self.fixture(root)
            index=args.index('--background-dir');del args[index:index+2]
            status,text,_=self.invoke(root,args,tool,['--reference-only','--check-inputs'])
            self.assertEqual(status,0)
            self.assertTrue(all(row['source_assignment']=='explicit_curated_reference_only'
                                for row in json.loads(text)['reference_mappings']))


if __name__=='__main__':unittest.main()
