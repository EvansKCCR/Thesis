import csv,hashlib,json,tempfile,unittest
from pathlib import Path
from collections import defaultdict
import master_library as ml
import family_motifs as fm


class MasterTests(unittest.TestCase):
    def test_current_scopes_and_workbook_parity(self):
        root=fm.HERE/'Motif_library'
        core=ml.load(root/'adhesion_core.json');full=ml.load(root/'motif_master.json','master')
        excel=ml.load(root/'adhesion_motif_master_evidence_based.xlsx')
        self.assertEqual(len(core['rules']),53);self.assertEqual(len(full['rules']),347)
        self.assertEqual({r['rule_id']:r['regex_python'] for r in core['rules']},
                         {r['rule_id']:r['regex_python'] for r in excel['rules']})
        self.assertTrue({r['rule_id'] for r in core['rules']}<={r['rule_id'] for r in full['rules']})
        self.assertEqual(len({r['regex_python'] for r in full['rules']}),347)
        self.assertFalse(core['context_filters_implemented'])

    def test_checksums_and_empty_patterns_are_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            path=Path(td)/'library.json'
            data=json.loads((fm.HERE/'Motif_library/adhesion_core.json').read_text())
            data['rules'][0]['regex_python']='CHANGED'
            path.write_text(json.dumps(data))
            with self.assertRaisesRegex(ValueError,'checksum'):ml.load(path)
            data['rules'][0]['regex_python']='A*'
            data['rules'][0]['regex_sha256']=hashlib.sha256(b'A*').hexdigest()
            path.write_text(json.dumps(data))
            with self.assertRaisesRegex(ValueError,'Empty-matching'):ml.load(path)

    def test_one_regex_scan_keeps_full_sequence_anchors_and_candidate_tier(self):
        with tempfile.TemporaryDirectory() as td:
            path=Path(td)/'library.json'
            data=json.loads((fm.HERE/'Motif_library/adhesion_core.json').read_text())
            r=data['rules'][0];r.update(regex_python='P..P$',regex_sha256=hashlib.sha256(b'P..P$').hexdigest())
            data['rules']=[r];path.write_text(json.dumps(data));snapshot=ml.load(path)
            family={'ARP7':[dict(id='x',sequence='PPPPAPPXPAPPPP')]}
            ml.scan(snapshot,family,defaultdict(list),[],Path(td),fm.write_tsv)
            with (Path(td)/'adhesion_motif_annotations.tsv').open() as h:rows=list(csv.DictReader(h,delimiter='\t'))
            self.assertEqual(len(rows),1);self.assertEqual(rows[0]['peptide'],'PPPP')
            self.assertEqual(rows[0]['end'],'14')
            self.assertEqual(rows[0]['candidate_tier'],'sequence_only_candidate')
            self.assertNotIn('confidence_score',rows[0])
