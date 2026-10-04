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
    def test_automatic_references_and_continuous_downstream(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            h,a,unused=self.fixture(root)
            target=root/'detection/TargetA/prepared/TargetA.faa'
            target.parent.mkdir(parents=True)
            w.write_fasta(target,{'TargetA__g':'MAAA'})
            w.seal(h,'handoff',target_species=['TargetA'],detection=str(root/'detection'),sources={str(target.resolve()):w.digest(target)})
            w.seal(a,'annotation',handoff=str(h),handoff_manifest=w.digest(h/'workflow_manifest.json'))
            raw=root/'reference_raw'
            raw.mkdir()
            (raw/'RefA.faa').write_text('>p1\nMAAA\n>p2\nMA\n')
            (raw/'RefA.gff3').write_text('chr\tx\tCDS\t1\t9\t.\t+\t0\tprotein_id=p1;Dbxref=GeneID:1\nchr\tx\tCDS\t1\t9\t.\t+\t0\tprotein_id=p2;Dbxref=GeneID:1\n')
            (raw/'RefB.faa').write_text('>b gene=geneB\nMAAA\n')
            (raw/'RefC.faa').write_text('>c\nMAAA\n')
            (raw/'RefC.protein_gene.tsv').write_text('protein_id\tgene_id\nc\tgeneC\n')
            actual_execute=w.execute
            commands=[]
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
                    return
                else:
                    actual_execute(command,logs)
            out=root/'automatic'
            with patch.object(w.shutil,'which',return_value='/fixture/tool'),patch.object(w,'execute',side_effect=simulated):
                w.downstream(argparse.Namespace(annotation=a,out=out,reference_raw=raw,threads=1,analysis_threads=1,compare=None,skip_phylogeny=False))
            self.assertTrue((out/'final/audited_candidates.xlsx').exists())
            self.assertEqual(w.fasta(out/'references/RefA/prepared/RefA.faa'),{'RefA__GeneID_1':'MAAA'})
            self.assertEqual(json.loads((out/'phylogeny/manifest.json').read_text())['duplication_evidence_status'],'unavailable_not_evidence_of_no_duplication')
            self.assertTrue(any('phylogeny.py' in c[1] and c[2]=='run' for c in commands))
            self.assertEqual(w.verify(out,'completed_workflow')['kind'],'completed_workflow')

    def test_general_target_preparation(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            inputs=root/'inputs'
            inputs.mkdir()
            config=json.loads((w.ROOT/'targets.example.json').read_text())
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
                with (base/'candidates/candidates.csv').open('w',newline='') as handle:
                    writer=csv.DictWriter(handle,fieldnames=list(rows[0]))
                    writer.writeheader(); writer.writerows(rows)
                with (base/'domains/domains.csv').open('w',newline='') as handle:
                    writer=csv.DictWriter(handle,fieldnames=['query_name','pfam_accession','ali_from','ali_to','decision'])
                    writer.writeheader()
                    writer.writerow(dict(query_name=sp+'__g',pfam_accession='PF00041',ali_from=1,ali_to=4,decision='retained'))
            handoff=root/'handoff'
            w.start(argparse.Namespace(detection=detection,input_dir=None,out=handoff,threads=1))
            rows=w.read(handoff/'provisional_hypotheses.tsv')
            self.assertEqual(len(rows),2)
            self.assertEqual({r['status'] for r in rows},{'unresolved'})
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
            self.assertEqual({r['family'] for r in audit},{'fibronectin_like_candidate'})

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
        seqs = {sp+'__g':'MAAA' for sp in ('TargetA',)}
        w.write_fasta(h/'external_candidates.faa',seqs)
        w.write(h/'hypotheses.tsv', [dict(protein_id=k,Specie=k.split('__')[0],family='fibronectin_like',motif_family='fibronectin_like_candidate',
               assigned_family='fibronectin_like',status='supported_candidate',domain_match='True',domain_type_fraction='1',
               domain_copy_fraction='1',architecture='PF00001:1-4',length=4) for k in seqs])
        w.seal(h,'handoff')
        a = root/'annotation'
        (a/'motifs').mkdir(parents=True)
        w.write(a/'motifs/adhesion_motif_annotations.tsv',[],['family','sequence_id','candidate_tier'])
        w.write(a/'motifs/used_regions.tsv',[],['family','sequence_id','region','start','end'])
        w.seal(a,'annotation',handoff=str(h),handoff_manifest=w.digest(h/'workflow_manifest.json'))
        entries=[]
        for sp in ('TargetA','NewRef'):
            p=root/f'{sp}.faa'
            w.write_fasta(p,{sp+'__g':'MAAA'})
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
            w.audit(argparse.Namespace(annotation=a,orthology=ort,out=final))
            joined=w.read(final/'joined_evidence.tsv')
            self.assertEqual(joined[0]['Orthology_direct_reference_count'],'1')
            from openpyxl import load_workbook
            wb=load_workbook(final/'audited_candidates.xlsx',read_only=True)
            rows=list(wb['Collapsed_candidates'].values)
            grade=rows[0].index('grade')
            self.assertEqual([r[grade] for r in rows[1:]],['A'])
            wb.close()
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
