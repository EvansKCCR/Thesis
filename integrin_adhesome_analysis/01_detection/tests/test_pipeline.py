import argparse
import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from adhesome_core import fasta, table, write_table
from adhesome_pipeline import prepare, classify_strict as classify, architecture_matches, blast_hits, scan, homology, discover_proteomes, batch
from build_adhesome_pfam_library import build
from parse_adhesome_domtblout import parse, FIELDS

class PipelineTests(unittest.TestCase):
    def test_reference_catalogue_accepts_computational_evidence_without_experimental_metadata(self):
        from adhesome_pipeline import catalogue
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);refs=root/'refs.faa';metadata=root/'refs.csv'
            refs.write_text('>ref\nMAAA\n')
            metadata.write_text('reference_id,family,module,species,accession,evidence_source\nref,actin,actin_linking,human,P1,sequence and domain evidence\n')
            self.assertEqual(set(catalogue(metadata,refs)),{'ref'})

    def test_selected_target_retry_preserves_completed_outputs_and_batch_status(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);out=root/'results'
            qc=dict(min_length=1,max_ambiguous_fraction=.1,max_ambiguous_run=5)
            entries=[dict(code=c,species=c,assembly='fixture',path=c+'.faa',gene_regex=r'gene=(\S+)') for c in ('TargetA','validation')]
            (root/'TargetA.faa').write_text('>one gene=one\nMAAA\n')
            (root/'validation.faa').write_text('>candidate|P2|Two gene=two\nMAAA\n')
            config=root/'targets.json';config.write_text(json.dumps(dict(qc=qc,proteomes=entries)))
            rules=root/'rules.json';rules.write_text('{}')
            manifest=root/'manifest.csv';manifest.write_text('Pfam accession\nPF00001\n')
            library=str(root/'library.hmm');Path(library+'.metadata.json').write_text('{}')
            args=argparse.Namespace(config=config,input_dir=root,out=out,prepare_only=False,library=library,
                rules=str(rules),manifest=str(manifest),crosswalk=None,reference_set=None,classification_mode='comparative',
                cpu=1,max_i=.001,max_c=.001,min_coverage=.35)
            from adhesome_core import digest
            attempts=[]
            def fake_scan(a):
                folder=Path(a.out);folder.mkdir()
                write_table(folder/'domains.csv',[],FIELDS)
                (folder/'scan.json').write_text(json.dumps(dict(metadata_sha256=digest(library+'.metadata.json'),
                    thresholds=dict(i=.001,c=.001,profile_coverage=.35))))
            def fake_classify(a):
                code=Path(a.proteome).stem;attempts.append(code)
                if code=='validation' and attempts.count(code)==1:raise ValueError('simulated classifier failure')
                folder=Path(a.out);folder.mkdir()
                (folder/'candidates.csv').write_text('protein_id,family,status\n'+code+'__one,test,provisional\n')
                (folder/'classification.json').write_text(json.dumps(dict(mode='comparative',checksums={str(Path(p).resolve()):digest(p) for p in (a.proteome,a.domains,a.rules,a.manifest)})))
            with patch('adhesome_pipeline.shutil.which',return_value='fixture/tool'),patch('adhesome_pipeline.scan',side_effect=fake_scan),patch('adhesome_pipeline.classify',side_effect=fake_classify),contextlib.redirect_stdout(io.StringIO()):
                with self.assertRaisesRegex(ValueError,'simulated classifier failure'):batch(args)
                retained={str(p.relative_to(out/'TargetA')):digest(p) for p in (out/'TargetA').rglob('*') if p.is_file()}
                entries[1]['gene_regex']=r'^[^|]+\|([^|]+)\|'
                config.write_text(json.dumps(dict(qc=qc,proteomes=entries)));args.species=['validation']
                batch(args)
                self.assertEqual(attempts,['TargetA','validation','validation'])
                self.assertEqual(retained,{str(p.relative_to(out/'TargetA')):digest(p) for p in (out/'TargetA').rglob('*') if p.is_file()})
                summary=json.loads((out/'batch.json').read_text())
                self.assertEqual(set(summary['species']),{'TargetA','validation'})
                self.assertTrue(all(r['status']=='completed_domain_screen' for r in summary['species'].values()))
                self.assertTrue((out/'archived_attempts/validation_001/proteome.json').exists())
                self.assertEqual(list(fasta(out/'validation/prepared/validation.faa')),[('validation__P2','MAAA')])
                (root/'TargetA.faa').write_text('>one gene=one\nMCCC\n')
                with self.assertRaisesRegex(ValueError,'preserved input'):batch(args)

    def test_input_check_only_does_not_require_tools_or_write_results(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            (root/'TargetA.faa').write_text('>a gene=one\nMAAA\n  >b gene=one\nMAAA\n')
            cfg=root/'targets.json'
            cfg.write_text(json.dumps(dict(proteomes=[dict(code='TargetA',species='TargetA',path='TargetA.faa',gene_regex=r'gene=(\S+)')])) )
            args=argparse.Namespace(config=cfg,input_dir=root,out=root/'out',check_inputs=True)
            with patch('adhesome_pipeline.shutil.which',return_value=None),contextlib.redirect_stdout(io.StringIO()) as output:
                batch(args)
            report=json.loads(output.getvalue())
            self.assertEqual(report['species']['TargetA']['proteins'],2)
            self.assertEqual(report['species']['TargetA']['genes'],1)
            self.assertFalse((root/'out').exists())

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def put(self, name, text):
        path = self.root / name
        path.write_text(text, encoding='utf-8')
        return path

    def test_preparation_collapse_before_filter_and_stable_ties(self):
        self.put('p.fa', '>b gene=g1\nAAAA\n>a gene=g1\nCCCC\n>c gene=g2\nAXXXXX\n>d gene=g2\nAAAA\n>e gene=g3\nACDE*\n')
        config = self.put('config.json', json.dumps(dict(qc=dict(min_length=3,max_ambiguous_fraction=.2,max_ambiguous_run=2),
            proteomes=[dict(code='X',species='test',assembly='test',path='p.fa',gene_regex=r'gene=(\S+)')])))
        with contextlib.redirect_stdout(io.StringIO()):
            prepare(config, self.root / 'out')
        self.assertEqual(list(fasta(self.root / 'out/X.faa')), [('X__g1','CCCC'),('X__g3','ACDE')])
        rows = table(self.root / 'out/X.id_map.csv')
        self.assertIn('ambiguous_fraction', rows[2]['decision'])
        self.assertEqual(rows[3]['decision'], 'shorter_isoform_or_tie')

    def test_missing_gene_fails(self):
        self.put('p.fa', '>a\nAAA\n')
        config = self.put('config.json', json.dumps(dict(qc=dict(min_length=1,max_ambiguous_fraction=.2,max_ambiguous_run=2),proteomes=[dict(code='X',path='p.fa',gene_regex=r'gene=(\S+)')])))
        with self.assertRaisesRegex(ValueError, 'No gene'):
            prepare(config, self.root / 'out')

    def test_build_missing_model_and_ga_split(self):
        sheet = self.put('m.csv', 'Pfam accession,Domain or family\nPF00001,A\nPF00002,B\n')
        hmm = self.put('p.hmm', 'HMMER3/f\nNAME A\nACC PF00001.1\nLENG 100\nGA 20 10;\n//\n')
        with self.assertRaisesRegex(ValueError, 'Missing Pfam'):
            build(sheet, hmm, self.root / 'out.hmm')
        with open(hmm, 'a') as handle:
            handle.write('HMMER3/f\nNAME B\nACC PF00002.2\nLENG 50\n//\n')
        models = build(sheet, hmm, self.root / 'out.hmm')
        self.assertEqual(models['PF00001']['ga'], [20,10])
        self.assertIsNone(models['PF00002']['ga'])
        self.assertIn('PF00002.2', (self.root / 'out.hmm.fallback.hmm').read_text())

    def row(self, acc='PF00001.1', score=20, full=30, start=1, end=50, hto=50, ce='1e-5', ie='1e-5'):
        return f'M {acc} 100 protein - 300 1e-5 {full} 0 1 1 {ce} {ie} {score} 0 1 {hto} {start} {end} {start} {end} 0.9 description\n'

    def test_ga_both_scores_coverage_and_overlap(self):
        source = self.put('d', self.row() + self.row(score=25,start=30,end=80) + self.row(start=100,end=150)
                          + self.row(full=5,start=200,end=250) + self.row(hto=10,start=260,end=270))
        rows = parse(source, {'PF00001':'A'}, {'PF00001':dict(accession='PF00001.1',ga=[20,10])})
        self.assertEqual([r['decision'] for r in rows], ['same_family_overlap','retained','retained','threshold_rejected','threshold_rejected'])

    def test_fallback_requires_both_evalues(self):
        source = self.put('d', self.row(ce='0.1') + self.row(ie='0.1',start=100,end=150))
        rows = parse(source, {'PF00001':'A'}, {'PF00001':dict(accession='PF00001.1',ga=None)})
        self.assertTrue(all(r['decision'] == 'threshold_rejected' for r in rows))

    def test_version_mismatch_fails(self):
        source = self.put('d', self.row())
        with self.assertRaisesRegex(ValueError, 'version mismatch'):
            parse(source, {'PF00001':'A'}, {'PF00001':dict(accession='PF00001.2',ga=None)})

    def test_malformed_data_fails(self):
        with self.assertRaisesRegex(ValueError, 'malformed'):
            parse(self.put('d','broken line\n'), {}, {})

    def test_domain_order(self):
        domains = [dict(pfam_accession='B'),dict(pfam_accession='A')]
        rule = dict(domain_alternatives=[dict(A=1,B=1)],ordered_domains=['A','B'])
        self.assertFalse(architecture_matches(domains,rule))
        self.assertTrue(architecture_matches(list(reversed(domains)),rule))

    def test_generic_domains_not_accepted_without_evidence(self):
        proteome = self.put('p.fa', '>p\n'+'A'*300+'\n>q\n'+'A'*300+'\n')
        domains = self.root / 'domains.csv'
        write_table(domains,[dict(query_name=p,pfam_accession='PF00069',ali_from=1,ali_to=200) for p in ['p','q']],FIELDS)
        rule = self.put('rules.json',json.dumps(dict(status='test',families=dict(kinase=dict(module='downstream_signalling',min_length=200,
                         diagnostic=False,domain_alternatives=[{'PF00069':1}])))))
        evidence = self.put('e.csv','protein_id,family,similarity,reciprocal,source\np,kinase,true,true,test\n')
        args = argparse.Namespace(proteome=proteome,domains=domains,rules=rule,homology=evidence,evidence=None,out=self.root/'out')
        classify(args)
        rows = table(self.root/'out/candidates.csv')
        self.assertEqual([r['status'] for r in rows],['supported_candidate','excluded'])

    def test_blast_coverage_uses_coordinates(self):
        path = self.put('b.tsv','r p 50 500 1 10 1 10 1e-30 100 100 100\nr q 50 50 1 50 1 50 1e-30 100 100 100\n')
        self.assertEqual([h['subject'] for h in blast_hits(path,1e-5,.3)], ['q'])

    def test_detection_length_gate_and_external_evidence_rejection(self):
        proteome = self.put('p.fa', '>long\n'+'A'*300+'\n>short\n'+'A'*30+'\n')
        domains = self.root / 'domains.csv'
        write_table(domains,[],FIELDS)
        rules = self.put('r.json', json.dumps(dict(status='test',families=dict(receptor=dict(
            module='integrin_receptors',min_length=200,domain_alternatives=[{'PF00001':1}],requires_topology=True)))))
        evidence = self.put('e.csv','protein_id,family,similarity,reciprocal,source\nlong,receptor,true,true,test\nshort,receptor,true,true,test\n')
        args = argparse.Namespace(proteome=proteome,domains=domains,rules=rules,homology=evidence,evidence=None,out=self.root/'out')
        classify(args)
        self.assertEqual([r['status'] for r in table(self.root/'out/candidates.csv')], ['provisional','excluded'])
        extra = self.put('t.csv','protein_id,family,topology,source\nlong,receptor,true,prediction\n')
        args.evidence = [extra]
        with self.assertRaisesRegex(ValueError,'resume'):
            classify(args)

    def test_nonfinite_scores_fail(self):
        with self.assertRaisesRegex(ValueError, 'non-finite'):
            parse(self.put('d',self.row(score='nan')), {'PF00001':'A'}, {'PF00001':dict(accession='PF00001.1',ga=None)})

    def test_recursive_discovery_and_output_exclusion(self):
        settings = dict(proteomes=[dict(code='Sman',species='Schistosoma mansoni',path='Smansoni_proteins.fa'),
            dict(code='Shae',species='Schistosoma haematobium',path='Shaematobium_proteins.fa'),
            dict(code='Sjap',species='Schistosoma japonicum',path='Sjaponicum_proteins.fa')])
        for species in settings['proteomes']:
            folder = self.root / species['species'] / 'nested'
            folder.mkdir(parents=True)
            (folder/'proteins.FASTA').write_text('>a gene=a\nAAA\n')
        self.put('unrelated.faa','>a\nAAA\n')
        generated = self.root/'results'
        generated.mkdir()
        (generated/'Sman.faa').write_text('>a\nAAA\n')
        found, ignored = discover_proteomes(self.root,settings,generated)
        self.assertEqual(set(found),{'Sman','Shae','Sjap'})
        self.assertEqual(len(ignored),1)
        with self.assertRaisesRegex(ValueError,'found 2'):
            discover_proteomes(self.root,settings)

    def test_recursive_discovery_requires_each_species(self):
        settings = dict(proteomes=[dict(code='Sman',species='Schistosoma mansoni',path='Smansoni_proteins.fa')])
        with self.assertRaisesRegex(ValueError,'found 0'):
            discover_proteomes(self.root,settings)

    def test_scan_routes_ga_and_fallback_with_shared_search_space(self):
        sheet = self.put('m.csv','Pfam accession,Domain or family\nPF00001,A\nPF00002,B\n')
        hmm = self.put('p.hmm','HMMER3/f\nNAME A\nACC PF00001.1\nLENG 100\nGA 20 10;\n//\nHMMER3/f\nNAME B\nACC PF00002.1\nLENG 100\n//\n')
        library = str(self.root/'library.hmm')
        build(sheet,hmm,library)
        proteome = self.put('p.fa','>protein\n'+'A'*300+'\n')
        commands = []
        def fake_run(command, log):
            commands.append(command)
            if command[0] == 'hmmscan':
                accession = 'PF00001.1' if '--cut_ga' in command else 'PF00002.1'
                Path(command[command.index('--domtblout')+1]).write_text(self.row(acc=accession))
        args = argparse.Namespace(out=self.root/'out',library=library,manifest=sheet,proteome=str(proteome),cpu=1,max_i=.001,max_c=.001,min_coverage=.35)
        with patch('adhesome_pipeline.run', side_effect=fake_run):
            scan(args)
        scans = [c for c in commands if c[0]=='hmmscan']
        self.assertEqual(len(scans),2)
        self.assertIn('--cut_ga',scans[0])
        self.assertNotIn('--cut_ga',scans[1])
        self.assertTrue(all(c[c.index('-Z')+1]=='2' and c[c.index('--domZ')+1]=='2' for c in scans))
        self.assertEqual(len(table(self.root/'out/domains.csv')),2)

    def test_reciprocal_competitor_ties_prevent_family_support(self):
        refs = self.put('ref.fa','>r\n'+'A'*100+'\n')
        cat = self.put('ref.csv','reference_id,family,module,species,accession,evidence_source,experimental_evidence\nr,kinase,downstream_signalling,human,R,PMID:test,test\n')
        broad = self.put('broad.fa','>b\n'+'A'*100+'\n>n\n'+'A'*100+'\n')
        broad_cat = self.put('broad.csv','reference_id,family\nb,kinase\nn,non_adhesome\n')
        proteome = self.put('p.fa','>p\n'+'A'*100+'\n')
        def fake_run(command, log):
            path = Path(command[command.index('-out')+1])
            lines = 'r p 50 100 1 100 1 100 1e-30 100 100 100\n' if path.name=='forward.tsv' else (
                'p b 50 100 1 100 1 100 1e-30 100 100 100\np n 50 100 1 100 1 100 1e-30 99.5 100 100\n')
            path.write_text(lines)
        args = argparse.Namespace(references=str(refs),catalogue=str(cat),metazoan=str(broad),metazoan_catalogue=str(broad_cat),
            proteome=str(proteome),out=self.root/'out',cpu=1,evalue=1e-5,coverage=.3,tie_fraction=.01)
        with patch('adhesome_pipeline.run',side_effect=fake_run):
            homology(args)
        result = table(self.root/'out/homology.csv')[0]
        self.assertEqual(result['reciprocal'],'false')
        self.assertEqual(result['reverse_families'],'kinase;non_adhesome')

if __name__ == '__main__':
    unittest.main()
