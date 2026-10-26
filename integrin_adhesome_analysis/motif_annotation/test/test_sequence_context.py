
# Resolve the organized bundle independently of the caller working directory.
import sys as _bundle_sys
from pathlib import Path as _BundlePath
_BUNDLE_ROOT=next(p for p in _BundlePath(__file__).resolve().parents if (p/'pipeline_paths.py').is_file())
if str(_BUNDLE_ROOT) not in _bundle_sys.path:_bundle_sys.path.insert(0,str(_BUNDLE_ROOT))
from pipeline_paths import activate as _activate_bundle_paths
_activate_bundle_paths()
import tempfile
import unittest
from collections import defaultdict
from pathlib import Path

import import_sequence_context as sc
import family_motifs as fm


class SequenceContextTests(unittest.TestCase):
    def test_exact_sequence_gate_and_cdd_deduplication(self):
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            (base/'DeepLoc_adhesome_candidate_list.csv').write_text(
                'Protein_ID,Localizations,Cytoplasm,Extracellular\na,Extracellular,0.1,0.9\nb,Cytoplasm,0.9,0.1\n')
            (base/'SignalP_adhesome_candidate_list.txt').write_text(
                '# ID\tPrediction\tOTHER\tSP(Sec/SPI)\tCS Position\n'
                'a\tSP\t0.1\t0.9\tCS pos: 3-4. Pr: 0.99\n'
                'b\tOTHER\t1\t0\t\n')
            hit='Q#1 - >a\tspecific\t1\t5\t8\t1e-9\t40\tcd123\tExample\t - \tcl1\n'
            (base/'NCBI_CDD_adhesome_hitdata.txt').write_text(
                'Query\tHit type\tPSSM-ID\tFrom\tTo\tE-Value\tBitscore\tAccession\tShort name\tIncomplete\tSuperfamily\n'
                + hit + hit + hit.replace('specific','non-specific')
                + 'Q#2 - >b\tspecific\t1\t2\t4\t1e-9\t40\tcd456\tWrong\t - \tcl1\n')
            families={'first':[dict(id='a',sequence='ABCDEFGHIJ')],
                      'second':[dict(id='a',sequence='ABCDEFGHIJ'),dict(id='b',sequence='ABCDEFGHIJ')]}
            audit=[dict(sequence_id='a',status='exact_sequence_verified'),
                   dict(sequence_id='b',status='sequence_mismatch')]
            regions=defaultdict(list)
            context,hits=sc.load(base,families,audit,regions)
            self.assertEqual(len(hits),2)
            self.assertEqual(context['first','a']['signalp_cleavage_after'],3)
            self.assertIn(('signal_peptide',1,3),regions['first','a'])
            self.assertIn(('CDD:cd123',5,8),regions['second','a'])
            self.assertEqual(context['second','b']['verification'],'no_exact_topology_sequence')
            self.assertEqual(context['second','b']['cdd_specific_hits'],0)

    def test_new_family_order_screen(self):
        rules={'families':{'Shc':{'domain_alternatives':[{'PF00640':1,'PF00017':1}],
                    'ordered_domain_groups':[['PF00640'],['PF00017']], 'support_requires_order':True}}}
        families={'Shc':[dict(id='x',sequence='A'*300)]}
        regions=defaultdict(list,{('Shc','x'):[('PF00017',10,100),('PF00640',150,250)]})
        self.assertEqual(fm.inspect_domains(families,rules,regions)[0]['state'],'domain_count_order_conflict')
        regions['Shc','x']=[('PF00640',10,100),('PF00017',150,250)]
        self.assertEqual(fm.inspect_domains(families,rules,regions)[0]['state'],'domain_count_order_rule_supported')


if __name__ == '__main__': unittest.main()
