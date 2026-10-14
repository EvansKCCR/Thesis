import tempfile
from collections import Counter
import unittest
import argparse
import importlib.util
import json
import subprocess
import sys
import types
from pathlib import Path
from unittest.mock import patch
from phylogeny import alignment_qc, fasta, load_sequences, annotate, run, trim_input, restore_trimmed, failure_category, family_plans
from pipeline import write_tsv, read_tsv, digest


class PhylogenyTests(unittest.TestCase):
    def test_current_resolved_broad_provisional_is_kept_as_broad_group(self):
        from openpyxl import Workbook
        with tempfile.TemporaryDirectory() as td:
            path=Path(td)/'audit.xlsx';book=Workbook();sheet=book.active;sheet.title='Resolved_candidates'
            sheet.append(['protein_id','assigned_family','broad_homology_group','classification_state','grade'])
            sheet.append(['T__a','unassigned','actin_superfamily','broad_only','C']);book.save(path)
            plans,_,coverage=family_plans([],None,'family',['T','R'],audit_workbook=path,
                target_species=['T'],sequences={'T__a':'MAAA'})
            self.assertEqual(plans[0][0],'family_actin_superfamily')
            self.assertEqual(plans[0][2],{'T__a'})
            self.assertEqual(coverage[0]['classification_state'],'broad_only')
            self.assertEqual(coverage[0]['assigned_family'],'unassigned')

    def test_audited_selection_excludes_other_targets_and_keeps_orphans(self):
        from openpyxl import Workbook
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            book=Workbook()
            sheet=book.active
            sheet.title='Collapsed_candidates'
            sheet.append(['protein_id','assigned_family','grade','family'])
            for row in [('old1','talin','A','wrong'),('T__b','talin','B','wrong'),
                        ('T__c','talin','C','wrong'),('T__d','talin','D','wrong'),
                        ('T__e','fibronectin_like','C','wrong')]:
                sheet.append(row)
            path=root/'audit.xlsx'
            book.save(path)
            write_tsv(root/'report.tsv',['candidate_id','prepared_id'],[dict(candidate_id='old1',prepared_id='T__a')])
            plans=[('OG1',{}, {'T__a','T__d','T__unselected','R__x'})]
            seqs={n:'MAAA' for n in ['T__a','T__b','T__c','T__d','T__e','T__unselected','R__x']}
            merged,links,coverage=family_plans(plans,None,'family',['T','R'],root/'report.tsv',
                audit_workbook=path,target_species=['T'],sequences=seqs)
            byname={g:(r,ids) for g,r,ids in merged}
            self.assertEqual(byname['family_talin'][1],{'T__a','T__b','T__c','R__x'})
            self.assertEqual(byname['family_fibronectin_like'][1],{'T__e'})
            self.assertEqual(byname['family_talin'][0]['_source_by_tip']['T__b'],[])
            self.assertEqual(Counter(r['grade'] for r in coverage),{'A':1,'B':1,'C':2})
            self.assertEqual(sum(r['status']=='included_without_source_group' for r in coverage),3)
            # Exercise CLI preparation when every audited target lacks a source group.
            write_tsv(root/'hogs.tsv',['HOG','T','R'],[])
            write_tsv(root/'orthologues.tsv',['candidate_id','orthologues','relationship','source_file'],[])
            (root/'species.json').write_text(json.dumps(['T','R']))
            (root/'targets.json').write_text(json.dumps(['T']))
            for sp in ('T','R'):
                fasta(root/f'{sp}.faa',{n:s for n,s in seqs.items() if n.startswith(sp+'__')})
            result=subprocess.run([sys.executable,str(Path(__file__).with_name('phylogeny.py')),'prepare',
                '--hogs',str(root/'hogs.tsv'),'--orthologues',str(root/'orthologues.tsv'),
                '--proteomes',str(root),'--species-list',str(root/'species.json'),
                '--target-species',str(root/'targets.json'),'--audit-workbook',str(path),
                '--candidate-report',str(root/'report.tsv'),'--allow-missing-duplications',
                '--output',str(root/'prepared')],capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stderr)
            _,summary=read_tsv(root/'prepared/groups.tsv')
            self.assertEqual(sum(int(r['candidates']) for r in summary),4)
            self.assertTrue(all(r['reference_sequences']=='0' for r in summary))
            self.assertTrue(all(r['status']=='too_few_taxa_for_supported_tree' for r in summary))
            del seqs['T__e']
            with self.assertRaisesRegex(ValueError,'missing from target proteomes'):
                family_plans(plans,None,'family',['T','R'],root/'report.tsv',
                    audit_workbook=path,target_species=['T'],sequences=seqs)

    def test_family_union_preserves_membership_and_ambiguity(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            write_tsv(root/'families.tsv',['protein_id','family'],[
                dict(protein_id='old1',family='talin'),dict(protein_id='T__b',family='talin'),
                dict(protein_id='old1',family='other'),dict(protein_id='T__missing',family='talin'),
                dict(protein_id='old1',family='talin')])
            write_tsv(root/'report.tsv',['candidate_id','prepared_id'],[dict(candidate_id='old1',prepared_id='T__a')])
            plans=[('OG1',{}, {'T__a','R__x'}),('OG2',{}, {'T__b','R__x','R__y'})]
            merged,links,coverage=family_plans(plans,root/'families.tsv','family',['T','R'],root/'report.tsv')
            byname={g:(r,ids) for g,r,ids in merged}
            row,ids=byname['family_talin']
            self.assertEqual(ids,{'T__a','T__b','R__x','R__y'})
            self.assertEqual(row['_source_by_tip']['R__x'],['OG1','OG2'])
            self.assertEqual(row['Matched_candidates'],'T__a, T__b')
            self.assertEqual(byname['family_other'][1],{'T__a','R__x'})
            self.assertEqual(len(links),3)
            self.assertTrue(any(r['prepared_id']=='T__missing' and r['status']=='not_in_selected_groups' for r in coverage))
            with self.assertRaisesRegex(ValueError,'Unknown families'):
                family_plans(plans,root/'families.tsv','family',['T','R'],selected=['absent'])
    @unittest.skipUnless(importlib.util.find_spec('Bio'), 'Biopython required')
    def test_rooting_and_duplication_evidence(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder=Path(tmp)
            work=folder/'inference'
            work.mkdir()
            tree=work/'tree'
            tree.write_text('(Cele__o:1,(Shae__a:1,Shae__b:1)99/100:1,Hsap__h:1);')
            write_tsv(folder/'orthology_edges.tsv',['candidate_id','orthologue_id'],
                      [dict(candidate_id='Shae__a',orthologue_id='Hsap__h')])
            write_tsv(folder/'duplication_evidence.tsv',['Genes 1','Genes 2'],
                      [{'Genes 1':'Shae__a','Genes 2':'Shae__b'}])
            annotate(tree,work,[])
            _,rows=read_tsv(work/'branch_evidence.tsv')
            self.assertTrue(all('unrooted split' in r['interpretation'] for r in rows))
            self.assertTrue(all(not r['exact_duplication_rows'] for r in rows))
            annotate(tree,work,['Cele__o'])
            _,rows=read_tsv(work/'branch_evidence.tsv')
            matched=[r for r in rows if r['exact_duplication_rows']=='1']
            self.assertEqual(len(matched),1)
            self.assertEqual(matched[0]['support_label'],'99/100')

    @unittest.skipUnless(importlib.util.find_spec('Bio'), 'Biopython required')
    def test_inference_driver_with_simulated_external_tools(self):
        # Tests orchestration only; this does not validate MAFFT/trimAl/IQ-TREE themselves.
        with tempfile.TemporaryDirectory() as tmp:
            project=Path(tmp)
            group=project/'OG1'
            group.mkdir()
            seqs={'Cele__o':'ACDE'*20,'Hsap__h':'ACDF'*20,'Shae__a':'FGDE'*20,'Shae__b':'FGDF'*20}
            fasta(group/'members.faa',seqs)
            write_tsv(group/'orthology_edges.tsv',['candidate_id','orthologue_id'],[])
            write_tsv(group/'duplication_evidence.tsv',['Genes 1','Genes 2'],[])
            (project/'manifest.json').write_text(json.dumps({'groups':[dict(group='OG1',status='ready',
                files={p.name:digest(p) for p in group.iterdir()})]}))
            args=argparse.Namespace(project=project,hog=None,outgroup=[],threads=1,seed=123,min_sites=50,
                max_gap=0.9,mafft='mafft',trimal='trimal',iqtree='iqtree2',aligner='linsi')
            def fake_run(command,**kwargs):
                if '--version' in command:
                    return argparse.Namespace(stdout='simulated version',stderr='',returncode=0)
                if command[0]=='mafft':
                    kwargs['stdout'].write((group/'members.faa').read_text())
                elif command[0]=='trimal':
                    Path(command[command.index('-out')+1]).write_text(Path(command[command.index('-in')+1]).read_text())
                else:
                    prefix=command[command.index('--prefix')+1]
                    Path(prefix+'.treefile').write_text('(Cele__o:1,Hsap__h:1,(Shae__a:1,Shae__b:1)90/99:1);')
                    self.assertIn('-bnni',command)
                    self.assertIn('-keep-ident',command)
                    self.assertEqual(command[command.index('-T')+1],'2')
                return argparse.Namespace(returncode=0)
            with patch('phylogeny.shutil.which',side_effect=lambda s:s),patch('phylogeny.subprocess.run',side_effect=fake_run):
                run(args)
            record=json.loads((group/'inference/run.json').read_text())
            self.assertEqual(record['status'],'complete')
            self.assertEqual(len(record['commands']),3)
            self.assertTrue((group/'inference/branch_evidence.tsv').is_file())

    def test_fasta_roundtrip_and_alignment_statistics(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'aligned.faa'
            seqs={'a':'ACDE','b':'ACDF','c':'FGDE','d':'FGDF'}
            fasta(path,seqs)
            self.assertEqual(load_sequences(path),seqs)
            _,_,qc=alignment_qc(path,seqs,4,0.9)
            self.assertEqual(qc['parsimony_informative_sites'],3)
            self.assertEqual(qc['distinct_aligned_sequences'],4)

    def test_trimming_aliases_and_qc_diagnostics(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            seqs={'Species__very_long_identifier_1':'ACDE','Species__very_long_identifier_2':'ACDF'}
            aliases=trim_input(seqs,root)
            restore_trimmed(root/'trimal_input.faa',aliases,root/'restored.faa')
            self.assertEqual(load_sequences(root/'restored.faa'),seqs)
            fasta(root/'bad.faa',{next(iter(aliases)):'ACDE'})
            with self.assertRaisesRegex(ValueError,'membership'):
                restore_trimmed(root/'bad.faa',aliases,root/'bad_restored.faa')
            self.assertTrue((root/'bad.membership.json').exists())
            fasta(root/'gaps.faa',{'x':'----','y':'ACDE'})
            with self.assertRaisesRegex(ValueError,'gaps'):
                alignment_qc(root/'gaps.faa',{'x','y'},1,0.9)
            self.assertTrue((root/'gaps.sequence_qc.tsv').exists())

    def test_tool_failure_retry_and_intentional_skip(self):
        with tempfile.TemporaryDirectory() as td:
            project=Path(td)
            group=project/'OG1'
            group.mkdir()
            fasta(group/'members.faa',{'a':'ACDE'*20,'b':'ACDF'*20,'c':'FGDE'*20,'d':'FGDF'*20})
            (project/'manifest.json').write_text(json.dumps({'groups':[
                dict(group='OG1',status='ready',files={'members.faa':digest(group/'members.faa')}),
                dict(group='OG2',status='too_few_taxa_for_supported_tree',files={})]}))
            args=argparse.Namespace(project=project,hog=None,outgroup=[],threads=8,iqtree_threads=2,seed=1,min_sites=50,
                max_gap=0.9,mafft='mafft',trimal='trimal',iqtree='iqtree2',aligner='linsi',allow_partial=True,retry_failed=False)
            bio=types.ModuleType('Bio')
            bio.Phylo=object()
            def fake(command,**kw):
                if '--version' in command:
                    return argparse.Namespace(stdout='test',stderr='')
                if command[0]=='mafft':
                    kw['stdout'].write((group/'members.faa').read_text())
                elif command[0]=='trimal':
                    Path(command[command.index('-out')+1]).write_text(Path(command[command.index('-in')+1]).read_text())
                else:
                    self.assertEqual(command[command.index('-T')+1],'2')
                    kw['stderr'].write('fixture IQ-TREE error; cause must be inspected')
                    raise subprocess.CalledProcessError(2,command)
            with patch.dict('sys.modules',{'Bio':bio}),patch('phylogeny.shutil.which',side_effect=lambda s:s),patch('phylogeny.subprocess.run',side_effect=fake):
                run(args)
                record=json.loads((group/'inference/run.json').read_text())
                self.assertEqual(record['failure_category'],'tool_failure')
                self.assertEqual(record['exit_code'],2)
                self.assertIn('fixture IQ-TREE error',record['stderr_tail'])
                summary=json.loads((project/'run_summary.json').read_text())
                self.assertEqual(summary['counts'],{'failed':1,'skipped':1})
                args.retry_failed=True
                args.allow_partial=False
                with self.assertRaisesRegex(ValueError,'groups failed'):
                    run(args)
                self.assertTrue((group/'inference_attempt_001/run.json').exists())
            self.assertEqual(failure_category(ValueError('excessive gaps'),'trimmed_alignment_qc'),'alignment_qc_failure')

    def test_rejects_changed_ids_lengths_and_poor_sequences(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'aligned.faa'
            for seqs,expected in [({'a':'AAAA'},['b']),({'a':'AAAA','b':'AAA'},['a','b']),
                                  ({'a':'----','b':'AAAA'},['a','b']),({'a':'XXXX'},['a'])]:
                fasta(path,seqs)
                with self.assertRaises(ValueError):
                    alignment_qc(path,expected,1,0.9)
            fasta(path,{'a':'AAA'})
            with self.assertRaises(ValueError):
                alignment_qc(path,['a'],50,0.9)


if __name__=='__main__':
    unittest.main()
