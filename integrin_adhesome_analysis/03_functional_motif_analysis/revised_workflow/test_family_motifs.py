import json
import tempfile
import unittest
from pathlib import Path
from collections import defaultdict
from unittest.mock import patch
import family_motifs as fm


class MotifTests(unittest.TestCase):
    def test_overlaps_and_ambiguity(self):
        self.assertEqual(list(fm.occurrences('P..P', 'PPPPP')), [(1, 4, 'PPPP'), (2, 5, 'PPPP')])
        self.assertEqual(list(fm.occurrences('NP.Y', 'NPXY')), [])

    def test_fasta_validation(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)/'a.fa'
            p.write_text('>a description\nnpay*\n')
            self.assertEqual(fm.read_fasta(p)[0]['sequence'], 'NPAY')
            for content in ['>a\nAA\n>a\nAA\n', '>a\nAA*AA\n', '>a\n', 'AAA\n']:
                p.write_text(content)
                with self.assertRaises(ValueError):
                    fm.read_fasta(p)

    def test_context_coordinates(self):
        spec = dict(context='c_terminal_window', window=10, regions=['cytoplasmic_tail'])
        self.assertEqual(fm.context_state(spec, 90, 94, 100, []), 'outside_context')
        self.assertEqual(fm.context_state(spec, 91, 94, 100, []), 'tail_window_only')
        self.assertEqual(fm.context_state(spec, 91, 94, 100, [('cytoplasmic_tail', 95, 100)]), 'outside_context')
        self.assertEqual(fm.context_state(spec, 95, 99, 100, [('cytoplasmic_tail', 95, 100)]), 'context_supported')
        self.assertEqual(fm.context_state(dict(context='domain', regions=['PF1']), 1, 3, 100, []), 'context_unverified')

    def test_dedup_background_and_singleton(self):
        rec = lambda i,s: dict(id=i, sequence=s)
        fam = {'a': [rec('a', 'ACDEFG'), rec('b', 'ACDEFG')], 'b': [rec('a', 'ACDEFG'), rec('c', 'LMNPQR')]}
        self.assertEqual(len(fm.unique(fam['a'])), 1)
        lib = dict(families={'a': {'motifs':['m']}, 'b': {'motifs':[]}}, motifs={'m':dict(pattern='ACD', interpretation='',source='')})
        _, _, p = fm.inspect(fam, lib, defaultdict(list))
        self.assertEqual((p[0]['target_total'],p[0]['other_total'],p[0]['other_positive']), (1,1,0))
        words, _ = fm.discover(fam, 3, 4, .5, 5)
        self.assertFalse(any(r['family']=='a' for r in words))

    def test_discovery_coordinates(self):
        fam = {'a':[dict(id='1',sequence='ACDEFG'),dict(id='2',sequence='ACDEFA')],
               'b':[dict(id='3',sequence='LMNPQR')]}
        words, sites = fm.discover(fam, 4, 5, 1, 10)
        self.assertTrue(any(r['motif']=='ACDEF' for r in words))
        for row in sites:
            seq = next(r['sequence'] for r in fam[row['family']] if r['id']==row['sequence_id'])
            self.assertEqual(seq[row['start']-1:row['end']], row['motif'])

    def test_region_validation_and_pwm_alphabet(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'r.tsv'
            p.write_text('family\tsequence_id\tregion\tstart\tend\na\tx\tPF1\t0\t4\n')
            with self.assertRaises(ValueError):
                fm.load_regions(p, {'a':[dict(id='x',sequence='ACDE')]})
            p.write_text('MEME version 5\nALPHABET= ACGT\nMOTIF x\n')
            with self.assertRaises(ValueError):
                fm.check_meme(p)

    def test_end_to_end_and_missing_tools(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d); inputs=root/'inputs'; inputs.mkdir()
            (inputs/'integrin_beta.fasta').write_text('>one\nACDENPAY\n>two\nACDENPVY\n')
            out=root/'out'
            with patch('family_motifs.shutil.which', return_value=None):
                self.assertEqual(fm.main(['--input-dir',str(inputs),'--out',str(out),'--external','run','--min-streme','2']),0)
            manifest=json.loads((out/'manifest.json').read_text())
            self.assertTrue(manifest['completed'])
            self.assertIn('executable unavailable',(out/'external_stages.tsv').read_text())
            self.assertIn('SLIM017',(out/'motif_annotations.tsv').read_text())
            with self.assertRaises(SystemExit):
                fm.main(['--input-dir',str(inputs),'--out',str(out)])

    def test_domain_counts_deduplicate_intervals(self):
        fam={'a':[dict(id='x',sequence='ACDEFG')]}
        lib={'families':{'a':{'domain_alternatives':[{'PF1':2}]}}}
        regions=defaultdict(list, {('a','x'):[('PF1',1,2),('PF1',1,2)]})
        self.assertEqual(fm.inspect_domains(fam,lib,regions)[0]['state'],'domain_count_rule_not_met')
        regions['a','x'].append(('PF1',4,6))
        self.assertEqual(fm.inspect_domains(fam,lib,regions)[0]['state'],'domain_count_rule_supported')

    def test_external_failure_blocks_dependents(self):
        import argparse
        import subprocess
        with tempfile.TemporaryDirectory() as d:
            args=argparse.Namespace(external='run',timeout=5,min_width=6,max_width=15,
                top=2,seed=1729,min_streme=2,qvalue=.05,known_meme=None)
            fam={'a':[dict(id='1',sequence='ACDEFG'),dict(id='2',sequence='LMNPQR')]}
            with patch('family_motifs.shutil.which', side_effect=lambda x:x), patch(
                    'family_motifs.subprocess.run', return_value=subprocess.CompletedProcess([],1,'','bad')):
                stages,_,_=fm.external(fam,Path(d),args)
            self.assertEqual(stages[0]['status'],'failed')
            self.assertEqual(stages[1]['status'],'skipped')
            commands=json.loads((Path(d)/'commands.json').read_text())
            self.assertIn('--protein',commands[0]['argv'])
            self.assertNotIn('--dna',commands[0]['argv'])


if __name__ == '__main__':
    unittest.main()
