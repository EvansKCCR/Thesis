"""Synthetic regressions. BLAST execution is mocked; project data is never used."""
import contextlib,csv,importlib.util,io,json,sys,tempfile,types,unittest
from pathlib import Path
from unittest.mock import patch
source=Path(__file__).resolve().parents[1]/'reciprocal_blast.py'
if not source.is_file():source=Path(__file__).resolve().parent/'reciprocal_blast.py'
spec=importlib.util.spec_from_file_location('parvin_tested',source)
p=importlib.util.module_from_spec(spec);sys.modules[spec.name]=p;spec.loader.exec_module(p)

class Tests(unittest.TestCase):
    def protein(self,name,seq,internal):return p.Protein(name,name,seq,internal)
    def fixture(self,root):
        (root/'mock_tool').write_text('synthetic executable')
        (root/'refs.faa').write_text('>curated\nMACDEFGHIKLMNPQRSTVWY\n')
        (root/'bg.faa').write_text('>competitor\nMYYYYYYYYYYYYYYYYYYY\n')
        (root/'sources.tsv').write_text('reference_file\treference_id\treference_proteome\nrefs.faa\tcurated\tbg\n')
        (root/'targets').mkdir()
        for i in range(7):(root/f'targets/Target{i}.faa').write_text('>target\nMACDEFGHIKLMNPQRSTVWY\n')
        return ['--references',str(root/'refs.faa'),'--proteomes',str(root/'targets'),
                '--reference-background',str(root/'bg.faa'),'--source-map',str(root/'sources.tsv'),'--outdir',str(root/'out')]
    def test_strict_exact_mapping(self):
        r=self.protein('r','MACDE','R1');bank=[self.protein('other','MACDE','B1')]
        mapping,methods=p.background_mapping([r],bank)
        self.assertEqual(mapping,{'R1':'B1'});self.assertEqual(methods['R1'],'identical_sequence');self.assertEqual(len(bank),1)
    def test_strict_reports_all_missing_without_mutation(self):
        refs=[self.protein('first','MACDE','R1'),self.protein('second','MWWWW','R2')]
        bank=[self.protein('bg','MYYYY','B1')];before=list(bank)
        with self.assertRaises(ValueError) as exc:p.background_mapping(refs,bank)
        self.assertIn('first',str(exc.exception));self.assertIn('second',str(exc.exception));self.assertEqual(bank,before)
    def test_augmentation_preserves_competitors(self):
        refs=[self.protein('exact','MACDE','R1'),self.protein('missing','MWWWW','R2')]
        bank=[self.protein('native','MACDE','B1'),self.protein('competitor','MYYYY','B2')];before=list(bank)
        mapping,methods=p.background_mapping(refs,bank,True)
        self.assertEqual(bank[:2],before);self.assertEqual(len(bank),3)
        self.assertEqual(bank[-1].seq,refs[-1].seq);self.assertEqual(methods['R2'],'curated_reference_added')
    def test_multiple_matches_still_fail_atomically(self):
        refs=[self.protein('missing','MWWWW','R1'),self.protein('ambiguous','MACDE','R2')]
        bank=[self.protein('one','MACDE','B1'),self.protein('two','MACDE','B2')];before=list(bank)
        with self.assertRaisesRegex(ValueError,'found 2 exact'):p.background_mapping(refs,bank,True)
        self.assertEqual(bank,before)
    def test_added_ids_do_not_collide(self):
        bank=[self.protein('curated_added|r','MYYYY','A0000001')]
        p.background_mapping([self.protein('r','MACDE','R1')],bank,True)
        self.assertEqual(len({x.original_id for x in bank}),2);self.assertEqual(len({x.internal_id for x in bank}),2)
    def test_empty_header(self):
        with tempfile.TemporaryDirectory() as td:
            f=Path(td)/'bad.faa';f.write_text('>\nMACDE\n')
            with self.assertRaisesRegex(ValueError,'Empty FASTA header'):p.read_fasta(f,'R')
    def test_stops_retained_and_summarized(self):
        with tempfile.TemporaryDirectory() as td:
            f=Path(td)/'stops.faa';f.write_text(''.join(f'>p{i}\nMA*CD*\n' for i in range(9)))
            qc=[];err=io.StringIO()
            with contextlib.redirect_stderr(err):rows=p.read_fasta(f,'P',qc,'target')
            self.assertEqual(len(err.getvalue().splitlines()),1);self.assertEqual(rows[0].seq,'MA*CD')
            self.assertEqual(len(qc),9);self.assertEqual(qc[0]['internal_stop_count'],1)
    def test_preflight_never_runs_tools_or_writes(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);args=self.fixture(root);stdout=io.StringIO()
            with patch.object(p.subprocess,'run',side_effect=AssertionError('BLAST executed')),\
                 patch.object(p.shutil,'which',side_effect=AssertionError('BLAST required')),contextlib.redirect_stdout(stdout):
                self.assertEqual(p.main(args+['--augment-reference-background','--check-inputs']),0)
            report=json.loads(stdout.getvalue());self.assertFalse(report['files_written'])
            self.assertEqual(report['curated_sequences_added'],1);self.assertEqual(len(report['target_proteomes']),7)
            self.assertFalse((root/'out').exists())
    def test_strict_failure_does_not_write(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);args=self.fixture(root)
            with patch.object(p.subprocess,'run',side_effect=AssertionError('BLAST executed')):
                with self.assertRaisesRegex(ValueError,'found 0 exact'):p.main(args)
            self.assertFalse((root/'out').exists())
    def test_augmentation_requires_background(self):
        with contextlib.redirect_stderr(io.StringIO()),self.assertRaises(SystemExit):
            p.main(['--references','x','--proteomes','y','--outdir','z','--reference-only','--augment-reference-background'])
    def test_nonfinite_arguments(self):
        for val in ('nan','inf'):
            with self.subTest(val=val),contextlib.redirect_stderr(io.StringIO()),self.assertRaises(SystemExit):
                p.main(['--references','x','--proteomes','y','--outdir','z','--evalue',val])
    def test_original_tie_and_filter_rules(self):
        args=types.SimpleNamespace(evalue=1e-5,min_identity=25,min_qcov=50,min_scov=50,min_bit_gap=0)
        hit={'bitscore':100,'evalue':1e-20,'pident':80,'qcov':100,'scov':100}
        self.assertTrue(p.passes(hit,args));self.assertTrue(p.ambiguous([hit,dict(hit)],args))
        self.assertFalse(p.passes(dict(hit,qcov=49),args))
    def fake_run(self,cmd,log,stdout_path=None,tie=False):
        log.write('MOCK synthetic command\n')
        if stdout_path is None:
            prefix=Path(cmd[cmd.index('-out')+1])
            for suffix in ('pin','phr','psq'):Path(str(prefix)+'.'+suffix).write_text('synthetic database')
            return
        def hit(q,s,score):return f'{q}\t{s}\t90\t20\t20\t20\t1\t20\t1\t20\t1e-20\t{score}\n'
        body=hit('R0000001','P0000001',100) if '.forward.' in stdout_path.name else hit('P0000001','A0000001',100)+hit('P0000001','B001_0000001',100 if tie else 80)
        stdout_path.write_text(body)
    def mocked_analysis(self,root,tie=False):
        args=self.fixture(root)
        with patch.object(p.shutil,'which',return_value=str(root/'mock_tool')),\
             patch.object(p.subprocess,'run',return_value=types.SimpleNamespace(stdout='synthetic version')),\
             patch.object(p,'run',side_effect=lambda *a,**kw:self.fake_run(*a,**kw,tie=tie)),contextlib.redirect_stderr(io.StringIO()):
            p.main(args+['--augment-reference-background'])
        return root/'out'
    def test_mocked_outputs_label_added_references(self):
        with tempfile.TemporaryDirectory() as td:
            out=self.mocked_analysis(Path(td))
            with (out/'homology_audit.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
            self.assertEqual(len(rows),7);self.assertTrue(all(x['status']=='RBH_candidate_augmented_reference' for x in rows))
            self.assertTrue(all(x['reverse_database_mode']=='full_reference_proteome_plus_curated' for x in rows))
            self.assertEqual((out/'rbh_candidates.faa').read_text().count('>'),7)
            manifest=json.loads((out/'run_manifest.json').read_text());self.assertEqual(manifest['input_report']['curated_sequences_added'],1)
            self.assertIn('B001_0000001',(out/'work/reference_databases/bg/proteins.faa').read_text());self.assertTrue((out/'sequence_qc.tsv').is_file())
    def test_mocked_reverse_ties_remain_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            out=self.mocked_analysis(Path(td),True)
            with (out/'homology_audit.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
            self.assertTrue(all(x['status']=='ambiguous_reverse_best' for x in rows));self.assertEqual((out/'rbh_candidates.faa').read_text(),'')
    def test_existing_output_not_overwritten(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);args=self.fixture(root);(root/'out').mkdir();saved=root/'out/run_manifest.json';saved.write_text('retained')
            with patch.object(p.shutil,'which',return_value=str(root/'mock_tool')),patch.object(p.subprocess,'run',side_effect=AssertionError('tool executed')):
                with self.assertRaisesRegex(ValueError,'Output is not empty'):p.main(args+['--augment-reference-background'])
            self.assertEqual(saved.read_text(),'retained')

if __name__=='__main__':unittest.main()
