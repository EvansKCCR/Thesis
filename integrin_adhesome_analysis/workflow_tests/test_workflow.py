import argparse
import csv
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import adhesome_workflow as w


class WorkflowTests(unittest.TestCase):
    def test_failed_start_resumes_in_place_and_excludes_failed_target_everywhere(self):
        with tempfile.TemporaryDirectory() as td:
            handoff=Path(td)/'handoff';detection=handoff/'detection'
            for code in ('Good','Bad'):
                base=detection/code
                for name in ('prepared','domains','candidates'):(base/name).mkdir(parents=True,exist_ok=True)
                (base/'prepared'/f'{code}.faa').write_text('>Good__g\nMAAA\n' if code=='Good' else 'not a protein FASTA\n')
                w.write(base/'domains/domains.csv',[dict(query_name=code+'__g',pfam_accession='PF00041',ali_from=1,ali_to=4,decision='retained')])
                # Detection tables are CSV, unlike the handoff TSVs.
                for path in (base/'domains/domains.csv',):path.write_text(path.read_text().replace('\t',','))
                w.write(base/'candidates/candidates.csv',[dict(protein_id=code+'__g',family='fibronectin_like',status='provisional')])
                path=base/'candidates/candidates.csv';path.write_text(path.read_text().replace('\t',','))
            w.dump(detection/'batch.json',dict(species={'Good':dict(status='completed_domain_screen'),'Bad':dict(status='failed')}))
            (handoff/'families').mkdir();(handoff/'logs').mkdir()
            (handoff/'logs/command_001.json').write_text('{"status":"failed"}')
            original=(handoff/'logs/command_001.json').read_bytes()
            self.assertEqual(w.main(['start','--resume','--exclude-target','Bad','--detection',str(detection),'--out',str(handoff)]),0)
            metadata=w.verify(handoff,'handoff')
            self.assertEqual(metadata['target_species'],['Good']);self.assertEqual(metadata['excluded_target_species'],['Bad'])
            self.assertFalse(any('/Bad/' in path.replace('\\','/') for path in metadata['sources']))
            self.assertEqual((handoff/'logs/command_001.json').read_bytes(),original)
            self.assertEqual(set(w.fasta(handoff/'external_candidates.faa')),{'Good__g'})
            self.assertFalse(list(handoff.glob('*backup*')))

    def test_phylogeny_cli_is_opt_in(self):
        for flags,expected in [([],False),(['--skip-phylogeny'],False),(['--phylogeny'],True),(['--phylogeny-grouping','family'],False)]:
            with patch.object(w,'downstream') as action:
                self.assertEqual(w.main(['continue','--annotation','fixture','--out','new_output',*flags]),0)
                self.assertEqual(action.call_args.args[0].phylogeny,expected)

    def test_reference_benchmark_runs_automatically_before_audit_without_phylogeny(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            handoff,annotation,unused=self.fixture(root)
            target=root/'detection/TargetA/prepared/TargetA.faa'
            target.parent.mkdir(parents=True)
            w.write_fasta(target,w.fasta(handoff/'external_candidates.faa'))
            w.seal(handoff,'handoff',target_species=['TargetA'],detection=str(root/'detection'),sources={str(target.resolve()):w.digest(target)})
            w.seal(annotation,'annotation',handoff=str(handoff),handoff_manifest=w.digest(handoff/'workflow_manifest.json'))
            raw=root/'inputs/reference_raw'
            inputs=raw/'BENCHMARKS/Input';truth=raw/'BENCHMARKS/RefOGs'
            inputs.mkdir(parents=True);truth.mkdir()
            (truth/'RefOG001.txt').write_text('a1\nb1\n')
            (truth/'RefOG002.txt').write_text('a2\nb2\n')
            for species,prefix in [('A','a'),('B','b')]:
                (inputs/(species+'.fa')).write_text(''.join('>'+prefix+str(i)+'\n'+('M*AA*' if i==0 else 'MAAA')+'\n' for i in range(3)))
            actual=w.execute;commands=[];fail_reference=[True]
            def simulated(command,logs):
                commands.append([str(x) for x in command])
                if 'prepare_references.py' in str(command[1]) and fail_reference[0]:
                    fail_reference[0]=False
                    folder=Path(command[command.index('--out')+1]);folder.mkdir()
                    (folder/'partial.txt').write_text('failed reference attempt')
                    raise ValueError('simulated reference preparation failure')
                if str(command[0])=='orthofinder':
                    result=Path(command[command.index('-o')+1])/'Results_fixture'
                    (result/'Phylogenetic_Hierarchical_Orthogroups').mkdir(parents=True)
                    (result/'Orthologues').mkdir()
                    w.write(result/'Phylogenetic_Hierarchical_Orthogroups/N0.tsv',[
                        dict(HOG='HOG'+str(i),TargetA='TargetA__g' if i==1 else '',Bench_A='Bench_A__a'+str(i),Bench_B='Bench_B__b'+str(i)) for i in range(3)])
                    for species,gene in [('Bench_A','a1'),('Bench_B','b1')]:
                        w.write(result/f'Orthologues/{species}.tsv',[dict(Orthogroup='OG1',TargetA='TargetA__g',**{species:species+'__'+gene})])
                else:
                    actual(command,logs)
            out=root/'automatic'
            with patch.object(w.shutil,'which',return_value='/fixture/tool'),patch.object(w,'execute',side_effect=simulated):
                with self.assertRaisesRegex(ValueError,'reference preparation failure'):
                    w.downstream(argparse.Namespace(annotation=annotation,out=out,reference_raw=raw,threads=1,analysis_threads=1,compare=None,benchmark_protocol='custom'))
                annotation_hash=w.digest(annotation/'workflow_manifest.json')
                w.recover(argparse.Namespace(run=out,threads=1,analysis_threads=None,iqtree_threads=None))
            self.assertEqual(annotation_hash,w.digest(annotation/'workflow_manifest.json'))
            self.assertTrue((out/'references_attempt_001/partial.txt').exists())
            summary=json.loads((out/'benchmark/summary.json').read_text())
            self.assertEqual(summary['status'],'validated',summary.get('diagnostics'))
            self.assertEqual(summary['validated_weight'],1.0)
            self.assertEqual(summary['mapped_reference_proteins'],6)
            self.assertEqual(summary['sequence_qc']['internal_stop_proteins'],2)
            joined=w.read(out/'final/joined_evidence.tsv')[0]
            self.assertEqual(joined['Orthology_benchmark_status'],'validated')
            self.assertEqual(float(joined['Orthology_benchmark_weight']),1.0)
            self.assertFalse(any('phylogeny.py' in ' '.join(command) for command in commands))
            self.assertEqual(w.verify(out,'completed_workflow')['phylogeny_status'],'not_requested')
            panel=json.loads((out/'references/panel.auto.json').read_text())
            self.assertTrue(all(row.get('unit_level')=='benchmark_protein' and row['gene_level'] is False for row in panel['proteomes'] if row['species'].startswith('Bench_')))
            benchmark_index=next(i for i,c in enumerate(commands) if 'benchmark.py' in c[1])
            audit_index=next(i for i,c in enumerate(commands) if 'integrated_adhesome_pipeline.py' in c[1])
            self.assertLess(benchmark_index,audit_index)
            def failed_benchmark(command,logs):
                if 'benchmark.py' in str(command[1]):raise ValueError('simulated benchmark execution failure')
                return simulated(command,logs)
            failed=root/'benchmark_failure_run'
            with patch.object(w.shutil,'which',return_value='/fixture/tool'),patch.object(w,'execute',side_effect=failed_benchmark):
                w.downstream(argparse.Namespace(annotation=annotation,out=failed,reference_raw=raw,threads=1,analysis_threads=1,compare=None,benchmark_protocol='custom'))
            self.assertTrue((failed/'final/Stage_2_output.xlsx').exists())
            self.assertEqual(w.read(failed/'final/joined_evidence.tsv')[0]['Orthology_benchmark_status'],'failed')
            self.assertEqual(w.verify(failed,'completed_workflow')['completion_status'],'completed_with_optional_failures')
            boundary=len(commands)
            with patch.object(w.shutil,'which',return_value='/fixture/tool'),patch.object(w,'execute',side_effect=simulated):
                w.recover(argparse.Namespace(run=failed,threads=1,analysis_threads=None,iqtree_threads=None))
            self.assertFalse(any(command[0]=='orthofinder' for command in commands[boundary:]))
            self.assertEqual(w.read(failed/'final/joined_evidence.tsv')[0]['Orthology_benchmark_status'],'validated')
            self.assertEqual(w.verify(failed,'completed_workflow')['completion_status'],'complete')
            # A changed ground-truth file cannot keep supplying a validated weight.
            (truth/'RefOG001.txt').write_text('a1\nb2\n')
            fields,sha=w.benchmark_evidence(out/'benchmark',out/'orthology',json.loads((out/'orthology/workflow_manifest.json').read_text()),out/'panel')
            self.assertEqual(fields['Orthology_benchmark_status'],'invalid_provenance')
            self.assertEqual(fields['Orthology_benchmark_weight'],'')

    def test_start_preflight_checks_all_targets_before_creating_output(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);raw=root/'input/targets';raw.mkdir(parents=True)
            (raw/'TargetA.faa').write_text('>a gene=gene1\nMAAA\n')
            (raw/'validation.faa').write_text('>candidate|P001|GENE1\nMAAA\n>hard_negative|P002|GENE2\nMAAA\n')
            entries=[dict(code=code,species=code,assembly='fixture',path=code+'.faa',gene_regex=r'gene=(\S+)') for code in ('TargetA','validation')]
            cfg=root/'targets.json';w.dump(cfg,dict(qc=dict(min_length=1,max_ambiguous_fraction=.1,max_ambiguous_run=5),proteomes=entries))
            args=argparse.Namespace(detection=None,input_dir=raw,targets_config=cfg,out=root/'handoff',threads=1,check_inputs=True)
            with patch.object(w,'execute') as execute:
                with self.assertRaisesRegex(ValueError,'validation: gene_regex'):
                    w.start(args)
                self.assertFalse(args.out.exists());execute.assert_not_called()
                entries[1]['gene_regex']=r'^[^|]+\|([^|]+)\|'
                w.dump(cfg,dict(proteomes=entries))
                w.start(args)
                self.assertFalse(args.out.exists());execute.assert_not_called()
                args.input_dir=root/'inputs/targets'
                with self.assertRaisesRegex(ValueError,'use --input-dir'):
                    w.start(args)
                self.assertFalse(args.out.exists());execute.assert_not_called()

    def test_automatic_references_and_continuous_downstream(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            h,a,unused=self.fixture(root)
            target=root/'detection/TargetA/prepared/TargetA.faa'
            target.parent.mkdir(parents=True)
            w.write_fasta(target,w.fasta(h/'external_candidates.faa'))
            w.seal(h,'handoff',target_species=['TargetA'],detection=str(root/'detection'),sources={str(target.resolve()):w.digest(target)})
            w.seal(a,'annotation',handoff=str(h),handoff_manifest=w.digest(h/'workflow_manifest.json'))
            raw=root/'reference_raw'
            raw.mkdir()
            (raw/'RefA.faa').write_text('>p1\nMAAA\n>p2\nMA\n')
            (raw/'RefA.gff3').write_text('chr\tx\tCDS\t1\t9\t.\t+\t0\tprotein_id=p1;Dbxref=GeneID:1\nchr\tx\tCDS\t1\t9\t.\t+\t0\tprotein_id=p2;Dbxref=GeneID:1\n')
            (raw/'RefB.faa').write_text('>b gene=geneB\nMAAA\n')
            (raw/'RefC.faa').write_text('>c\nMAAA\n')
            w.dump(raw/'RefC.reference.json',dict(fasta_sha256=w.digest(raw/'RefC.faa'),complete_proteome=False,panel_scope='fixture_targeted_panel'))
            (raw/'RefC.protein_gene.tsv').write_text('protein_id\tgene_id\nc\tgeneC\n')
            actual_execute=w.execute
            commands=[]
            fail_tree=[True]
            fail_extract=[True]
            def simulated(command,logs):
                commands.append([str(x) for x in command])
                if str(command[0])=='orthofinder':
                    result=Path(command[command.index('-o')+1])/'Results_fixture'
                    (result/'Phylogenetic_Hierarchical_Orthogroups').mkdir(parents=True)
                    (result/'Orthologues').mkdir()
                    w.write(result/'Phylogenetic_Hierarchical_Orthogroups/N0.tsv',[
                        dict(HOG='HOG1',TargetA='TargetA__g',RefA='RefA__GeneID_1',RefB='RefB__geneB',RefC='RefC__geneC')])
                    for ref,gene in [('RefA','GeneID_1'),('RefB','geneB'),('RefC','geneC')]:
                        w.write(result/f'Orthologues/{ref}.tsv',[dict(Orthogroup='OG1',TargetA='TargetA__g',**{ref:ref+'__'+gene})])
                elif 'phylogeny.py' in str(command[1]) and command[2]=='run':
                    # External inference is simulated; preparation is exercised for real.
                    if fail_tree[0]:
                        fail_tree[0]=False
                        raise ValueError('simulated tree tool failure')
                    w.dump(Path(command[command.index('--project')+1])/'run_summary.json',dict(status='complete'))
                    return
                elif 'pipeline.py' in str(command[1]) and command[2]=='extract' and fail_extract[0]:
                    fail_extract[0]=False
                    raise ValueError('simulated extraction failure')
                else:
                    actual_execute(command,logs)
            out=root/'automatic'
            with patch.object(w.shutil,'which',return_value='/fixture/tool'),patch.object(w,'execute',side_effect=simulated):
                with self.assertRaisesRegex(ValueError,'simulated extraction failure'):
                    w.downstream(argparse.Namespace(annotation=a,out=out,reference_raw=raw,threads=1,analysis_threads=1,compare=None,phylogeny=True,phylogeny_grouping='family'))
                self.assertTrue((out/'orthology/inference_complete.json').exists())
                w.recover(argparse.Namespace(run=out,threads=2,analysis_threads=None,iqtree_threads=2))
                self.assertEqual(sum(c[0]=='orthofinder' for c in commands),1)
                state=json.loads((out/'recovery.json').read_text())
                self.assertEqual(state['stages']['orthology']['status'],'complete')
                self.assertEqual(state['stages']['phylogeny_run']['status'],'failed')
                self.assertEqual(state['stages']['audit']['status'],'complete')
                self.assertTrue((out/'final/Stage_2_output.xlsx').exists())
                self.assertEqual(w.verify(out,'completed_workflow')['completion_status'],'completed_with_optional_failures')
                self.assertEqual(json.loads((out/'benchmark/summary.json').read_text())['status'],'unavailable')
                audit_hash=w.digest(out/'final/Stage_2_output.xlsx')
                w.recover(argparse.Namespace(run=out,threads=2,analysis_threads=None,iqtree_threads=2,phylogeny=True))
                self.assertEqual(audit_hash,w.digest(out/'final/Stage_2_output.xlsx'))
                # Adopted recovery must refresh an audit built with superseded rules.
                w.archive_attempt(out/'workflow_manifest.json')
                audit_meta=out/'final/workflow_manifest.json'
                saved=json.loads(audit_meta.read_text());saved['rules_sha256']='superseded-rule-hash'
                w.dump(audit_meta,saved)
                boundary=len(commands)
                # Simulate an older run without a recovery file; adopt only verified stages.
                (out/'recovery.json').unlink()
                w.recover(argparse.Namespace(run=out,threads=2,analysis_threads=None,iqtree_threads=2,
                    adopt_existing=True,annotation=a,phylogeny=True))
                recovered=commands[boundary:]
                self.assertFalse(any(c[0]=='orthofinder' for c in recovered))
                self.assertFalse(any('prepare_references.py' in c[1] for c in recovered))
                self.assertTrue(any('--retry-failed' in c for c in recovered))
            self.assertTrue((out/'final/Stage_2_output.xlsx').exists())
            self.assertEqual(json.loads((out/'final/workflow_manifest.json').read_text())['rules_sha256'],w.digest(w.audit_configuration()[0]))
            self.assertEqual(w.fasta(out/'references/RefA/prepared/RefA.faa'),{'RefA__GeneID_1':'MAAA'})
            self.assertFalse(json.loads((out/'references/RefC/source.json').read_text())['reference_scope']['complete_proteome'])
            self.assertEqual(json.loads((out/'phylogeny/manifest.json').read_text())['duplication_evidence_status'],'unavailable_not_evidence_of_no_duplication')
            self.assertEqual(json.loads((out/'phylogeny/manifest.json').read_text())['grouping'],'family')
            self.assertTrue((out/'phylogeny/family_fibronectin_like/tips.tsv').exists())
            self.assertEqual(json.loads((out/'phylogeny/manifest.json').read_text())['selected_candidates'],1)
            self.assertTrue(any('phylogeny.py' in c[1] and c[2]=='run' for c in commands))
            self.assertEqual(w.verify(out,'completed_workflow')['kind'],'completed_workflow')

    def test_general_target_preparation(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            inputs=root/'inputs'
            inputs.mkdir()
            config = json.loads((w.ROOT.parent / "README" / "targets.example.json").read_text(encoding="utf-8"))
            config['proteomes']=[]
            for code in ('TargetX','TargetY'):
                filename=f'custom_{code}_shared.faa'
                (inputs/filename).write_text('>protein1 gene=gene1\n'+'M'+'A'*40+'\n')
                config['proteomes'].append(dict(code=code,species='Example shared',assembly='fixture',
                    path=filename,gene_regex=r'gene=(\S+)'))
            w.dump(root/'targets.json',config)
            w.execute([sys.executable,w.DETECT/'adhesome_pipeline.py','batch','--prepare-only',
                '--config',root/'targets.json','--input-dir',inputs,'--out',root/'prepared'],root/'logs')
            for code in ('TargetX','TargetY'):
                self.assertEqual(set(w.fasta(root/'prepared'/code/'prepared'/f'{code}.faa')),{code+'__gene1'})

    def test_fibronectin_screening_handoff_and_return(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            detection = root/'detection'
            for sp in ('TargetA',):
                base=detection/sp
                for folder in ('prepared','candidates','domains'):
                    (base/folder).mkdir(parents=True)
                w.write_fasta(base/'prepared'/f'{sp}.faa',{sp+'__g':'MAAA'})
                rows=[dict(protein_id=sp+'__g',family=f,status='unresolved',domain_match='True')
                      for f in ('fibronectin_like','fibronectin_III_ligand')]
                for row in rows:
                    row.update(similarity='no_qualifying_hit',similarity_reason='completed_search_no_qualifying_hit',
                               reciprocal='not_tested',reciprocal_reason='no_qualifying_forward_hit_for_reciprocal_assessment')
                rows.append(dict(rows[0],family='component_reference_candidate',status='provisional',similarity='supported',
                    similarity_reason='component_evidence_only_not_family_assignment',crosswalk_reference_components='FIXTURE'))
                for row in rows:
                    row.setdefault('crosswalk_reference_components','')
                with (base/'candidates/candidates.csv').open('w',newline='') as handle:
                    writer=csv.DictWriter(handle,fieldnames=list(rows[0]))
                    writer.writeheader(); writer.writerows(rows)
                with (base/'domains/domains.csv').open('w',newline='') as handle:
                    writer=csv.DictWriter(handle,fieldnames=['query_name','pfam_accession','ali_from','ali_to','decision'])
                    writer.writeheader()
                    writer.writerow(dict(query_name=sp+'__g',pfam_accession='PF00041',ali_from=1,ali_to=4,decision='retained'))
                with (base/'candidates/component_hypotheses.csv').open('w',newline='') as handle:
                    writer=csv.DictWriter(handle,fieldnames=['protein_id','component_symbol','interpretation'])
                    writer.writeheader()
                    writer.writerow(dict(protein_id=sp+'__g',component_symbol='FIXTURE',interpretation='fixture overlap'))
                with (base/'candidates/reference_coverage.csv').open('w',newline='') as handle:
                    writer=csv.DictWriter(handle,fieldnames=['family','reference_count','calibration_action'])
                    writer.writeheader()
                    writer.writerow(dict(family='fibronectin_like',reference_count=1,calibration_action=''))
            handoff=root/'handoff'
            w.start(argparse.Namespace(detection=detection,input_dir=None,out=handoff,threads=1))
            manifest=json.loads((handoff/'workflow_manifest.json').read_text())
            self.assertEqual(manifest['workflow_version'],'1.0')
            self.assertEqual(manifest['completed_stages'],['detection'])
            self.assertEqual(manifest['deferred_stages'],['topology','localization','phylogeny'])
            self.assertFalse((handoff/'downstream').exists())
            self.assertFalse((handoff/'motifs').exists())
            self.assertEqual(w.read(handoff/'component_hypotheses.tsv')[0]['component_symbol'],'FIXTURE')
            rows=w.read(handoff/'provisional_hypotheses.tsv')
            self.assertEqual(len(rows),3)
            self.assertEqual({r['similarity'] for r in rows},{'no_qualifying_hit','supported'})
            self.assertEqual(w.read(handoff/'reference_coverage.tsv')[0]['Specie'],'TargetA')
            self.assertEqual({r['status'] for r in rows},{'unresolved','provisional'})
            self.assertEqual(len(w.fasta(handoff/'external_candidates.faa')),1)
            self.assertEqual(len(w.fasta(handoff/'families/fibronectin_like_candidate.faa')),1)
            self.assertEqual(len(w.read(handoff/'fibronectin_like_hypotheses.tsv')),2)
            topology=root/'returned.3line'
            topology.write_text(''.join(f'>{sp}__g\nMAAA\nOOOO\n' for sp in ('TargetA',)))
            evidence=w.read(handoff/'external_evidence_template.tsv')
            evidence[0].update(DeepLoc_2_1='Extracellular',source='fixture prediction')
            w.write(root/'returned.tsv',evidence)
            annotation=root/'annotation'
            w.resume(argparse.Namespace(handoff=handoff,out=annotation,topology_file=topology,
                context_dir=None,external_evidence=root/'returned.tsv',evidence=None,allow_missing_context=False,annotation_only=True))
            audit=w.read(annotation/'motifs/topology_audit.tsv')
            self.assertEqual({r['status'] for r in audit},{'exact_sequence_verified'})
            self.assertEqual({r['family'] for r in audit},{'fibronectin_like_candidate','component_reference_candidate'})

    def test_fasta_rejects_duplicate_and_invalid(self):
        with tempfile.TemporaryDirectory() as td:
            p = Path(td)/'x.faa'
            for text in ('>x\nMAA\n>x\nMAA\n', '>x\nMA*A\n', 'MAA\n'):
                p.write_text(text)
                with self.assertRaises(ValueError):
                    w.fasta(p)

    def test_external_linkage_and_source(self):
        with tempfile.TemporaryDirectory() as td:
            p = Path(td)/'e.tsv'
            row = dict(protein_id='Shae__g', sequence_sha256=w.hashlib.sha256(b'MAAA').hexdigest(),
                       DeepLoc_2_1='Cytoplasm', DeepTMHMM='',SignalP_v6_0='',NCBI_CDD_architecture='',source='predictor run 1')
            w.write(p,[row])
            self.assertEqual(len(w.context_rows(p, {'Shae__g':'MAAA'})),1)
            for change in ({'source':''},{'sequence_sha256':'bad'},{'protein_id':'unknown'}):
                w.write(p,[dict(row,**change)])
                with self.assertRaises(ValueError):
                    w.context_rows(p, {'Shae__g':'MAAA'})

    def test_sealed_output_detects_changes(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            (root/'data').write_text('original')
            w.seal(root,'test')
            w.verify(root,'test')
            (root/'data').write_text('changed')
            with self.assertRaises(ValueError):
                w.verify(root,'test')

    def fixture(self, root):
        h = root/'handoff'
        h.mkdir()
        seqs = {sp+'__g':'M'+'A'*499 for sp in ('TargetA',)}
        w.write_fasta(h/'external_candidates.faa',seqs)
        w.write(h/'hypotheses.tsv', [dict(protein_id=k,Specie=k.split('__')[0],family='fibronectin_like',motif_family='fibronectin_like_candidate',
               assigned_family='fibronectin_like',status='supported_candidate',domain_match='True',domain_type_fraction='1',
               domain_copy_fraction='1',architecture='PF00039:1-80;PF00040:90-150;PF00041:160-250',length=500,
               reference_subfamily='fibronectin_like',subfamily_support_method='curated_reference',reference_accession='FIXTURE') for k in seqs])
        w.seal(h,'handoff')
        a = root/'annotation'
        (a/'motifs').mkdir(parents=True)
        w.write(a/'motifs/adhesion_motif_annotations.tsv',[],['family','sequence_id','candidate_tier'])
        w.write(a/'motifs/used_regions.tsv',[],['family','sequence_id','region','start','end'])
        w.seal(a,'annotation',handoff=str(h),handoff_manifest=w.digest(h/'workflow_manifest.json'))
        entries=[]
        for sp in ('TargetA','NewRef'):
            p=root/f'{sp}.faa'
            w.write_fasta(p,{sp+'__g':'M'+'A'*499})
            entries.append(dict(species=sp,fasta=p.name,gene_level=True,source='test fixture'))
        config=root/'panel.json'
        w.dump(config,dict(proteomes=entries))
        panel=root/'panel'
        w.calibrate(argparse.Namespace(handoff=h,panel=config,out=panel))
        return h,a,panel

    def test_custom_panel_extraction_join_and_audit(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            h,a,panel=self.fixture(root)
            results=root/'Results_test'
            (results/'Phylogenetic_Hierarchical_Orthogroups').mkdir(parents=True)
            (results/'Orthologues').mkdir()
            row=dict(HOG='N0.HOG0001',OG='OG0001',Gene_Tree_Parent_Clade='n0',
                     TargetA='TargetA__g',NewRef='NewRef__g')
            w.write(results/'Phylogenetic_Hierarchical_Orthogroups/N0.tsv',[row])
            for sp in ('TargetA',):
                w.write(results/f'Orthologues/{sp}__v__NewRef.tsv',[{'Orthogroup':'OG0001',sp:sp+'__g','NewRef':'NewRef__g'}])
            ort=root/'orthology'
            w.orthology(argparse.Namespace(panel=panel,out=ort,results=results,
                input_manifest=panel/'workflow_manifest.json',compare=None))
            report=w.read(ort/'candidates/candidate_report.tsv')
            self.assertEqual(report[0]['NewRef_orthologues'],'NewRef__g')
            final=root/'final'
            def evaluation_after_calls(audit_dir,handoff):
                self.assertTrue((audit_dir/'Stage_2_output.xlsx').exists())
                self.assertEqual(handoff,h)
                return dict(status='validation_target_not_present')
            with patch.object(w,'evaluate_controls',side_effect=evaluation_after_calls) as evaluate:
                w.audit(argparse.Namespace(annotation=a,orthology=ort,out=final))
            evaluate.assert_called_once_with(final,h)
            self.assertEqual(w.verify(final,'audit')['validation_evaluation']['status'],'validation_target_not_present')
            joined=w.read(final/'joined_evidence.tsv')
            self.assertEqual(joined[0]['Orthology_direct_reference_count'],'1')
            from openpyxl import load_workbook
            wb=load_workbook(final/'Stage_2_output.xlsx',read_only=True)
            rows=list(wb['Audited_family_assignment'].values)
            grade=rows[0].index('grade')
            # Complete family architecture plus orthology gives A. Missing
            # localization loses its bonus but does not override the family gate.
            self.assertEqual([r[grade] for r in rows[1:]],['A'])
            self.assertIn('domain_score_0_4',rows[0])
            self.assertIn('putative_audited_family',rows[0])
            self.assertIn('protein_name',rows[0])
            wb.close()
            failed_evaluation=root/'final_eval_error'
            with patch('evaluate_validation.evaluate',side_effect=ValueError('simulated invalid validation input')):
                w.audit(argparse.Namespace(annotation=a,orthology=ort,out=failed_evaluation))
            self.assertTrue((failed_evaluation/'Stage_2_output.xlsx').exists())
            self.assertEqual(w.verify(failed_evaluation,'audit')['validation_evaluation']['status'],'failed')
            # Add a reference, rerun extraction and compare by candidate ID.
            cfg=json.loads((root/'panel.json').read_text())
            w.write_fasta(root/'Extra.faa',{'Extra__g':'MAAA'})
            cfg['proteomes'].append(dict(species='Extra',fasta='Extra.faa',gene_level=True,source='extra fixture'))
            w.dump(root/'expanded.json',cfg)
            expanded=root/'expanded'
            w.calibrate(argparse.Namespace(handoff=h,panel=root/'expanded.json',out=expanded))
            row['Extra']='Extra__g'
            w.write(results/'Phylogenetic_Hierarchical_Orthogroups/N0.tsv',[row])
            for sp in ('TargetA',):
                w.write(results/f'Orthologues/{sp}__v__Extra.tsv',[{'Orthogroup':'OG0001',sp:sp+'__g','Extra':'Extra__g'}])
            rerun=root/'orthology_expanded'
            w.orthology(argparse.Namespace(panel=expanded,out=rerun,results=results,
                input_manifest=expanded/'workflow_manifest.json',compare=ort))
            comparison=w.read(rerun/'calibration_comparison.tsv')
            self.assertEqual(comparison[0]['previous_references'],'NewRef')
            self.assertEqual(comparison[0]['current_references'],'NewRef;Extra')
            # A changed panel may never authorize another audit.
            (panel/'proteomes/NewRef.faa').write_text('>NewRef__g\nMCCC\n')
            with self.assertRaises(ValueError):
                w.audit(argparse.Namespace(annotation=a,orthology=ort,out=root/'bad'))

    def test_panel_rejects_candidate_sequence_mismatch(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            h,a,panel=self.fixture(root)
            (root/'TargetA.faa').write_text('>TargetA__g\nMCCC\n')
            with self.assertRaisesRegex(ValueError,'sequence changed'):
                w.calibrate(argparse.Namespace(handoff=h,panel=root/'panel.json',out=root/'badpanel'))

    def test_target_without_candidates_still_required_in_panel(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            h,a,panel=self.fixture(root)
            w.seal(h,'handoff',target_species=['TargetA','EmptyTarget'])
            with self.assertRaisesRegex(ValueError,'every handoff target'):
                w.calibrate(argparse.Namespace(handoff=h,panel=root/'panel.json',out=root/'missing'))
            config=json.loads((root/'panel.json').read_text())
            w.write_fasta(root/'EmptyTarget.faa',{'EmptyTarget__g':'MAAA'})
            config['proteomes'].append(dict(species='EmptyTarget',fasta='EmptyTarget.faa',gene_level=True,source='fixture'))
            w.dump(root/'panel_complete.json',config)
            w.calibrate(argparse.Namespace(handoff=h,panel=root/'panel_complete.json',out=root/'complete'))
            self.assertEqual(json.loads((root/'complete/targets.json').read_text()),['EmptyTarget','TargetA'])


if __name__ == '__main__':
    unittest.main()
