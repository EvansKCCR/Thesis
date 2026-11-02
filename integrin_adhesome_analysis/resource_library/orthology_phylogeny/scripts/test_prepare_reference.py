
# Resolve source resources independently of the caller working directory.
import sys as _layout_sys
from pathlib import Path as _LayoutPath
_LAYOUT_BUNDLE_ROOT = next(p for p in _LayoutPath(__file__).resolve().parents
                          if (p/'config/layout_paths.json').is_file() and (p/'pipeline_paths.py').is_file())
if str(_LAYOUT_BUNDLE_ROOT) not in _layout_sys.path:
    _layout_sys.path.insert(0, str(_LAYOUT_BUNDLE_ROOT))
from pipeline_paths import activate as _layout_activate, resource_path as _layout_resource
_layout_activate()

import tempfile
import json
import unittest
from pathlib import Path
from prepare_reference import prepare, normalize_sequence
from pipeline import read_tsv, digest


class ReferencePreparationTests(unittest.TestCase):
    def test_collapse_and_mapping_failure(self):
        with tempfile.TemporaryDirectory() as td:
            root=_layout_resource(Path(td))
            fa=_layout_resource(root/'raw.faa')
            fa.write_text('>p2\nMAAA*\n>p1\nMAAA\n>p3\nMA\n')
            mapping=_layout_resource(root/'map.tsv')
            mapping.write_text('protein_id\tgene_id\np2\tg1\np1\tg1\np3\tg2\n')
            out=prepare('RefA',fa,mapping,_layout_resource(root/'prepared'))
            self.assertEqual(_layout_resource(out).read_text(),'>RefA__g1\nMAAA\n>RefA__g2\nMA\n')
            self.assertIn('p1\tg1\tRefA__g1\tp1\t1', (_layout_resource(root/'prepared/identifier_map.tsv')).read_text())
            mapping.write_text('protein_id\tgene_id\np1\tg1\n')
            with self.assertRaisesRegex(ValueError,'missing gene mapping'):
                prepare('RefA',fa,mapping,_layout_resource(root/'bad'))
            self.assertFalse((_layout_resource(root/'bad')).exists())

    def test_benchmark_stop_masking_preserves_ids_positions_and_provenance(self):
        with tempfile.TemporaryDirectory() as td:
            root=_layout_resource(Path(td))
            fa=_layout_resource(root/'raw.fa');original='>p1\nMA**G*\n>p2\nMCGG\n'
            fa.write_text(original)
            mapping=_layout_resource(root/'map.tsv');mapping.write_text('protein_id\tgene_id\np1\tp1\np2\tp2\n')
            out=prepare('Bench_A',fa,mapping,_layout_resource(root/'prepared'),internal_stop_policy='mask')
            self.assertEqual(_layout_resource(out).read_text(),'>Bench_A__p1\nMAXXG\n>Bench_A__p2\nMCGG\n')
            self.assertEqual(_layout_resource(fa).read_text(),original)
            fields,rows=read_tsv(_layout_resource(out.parent/'sequence_qc.tsv'))
            self.assertEqual(len(rows),1)
            self.assertEqual(rows[0]['masked_positions_1based'],'3;4')
            self.assertEqual(rows[0]['terminal_stop_removed'],'1')
            meta=json.loads((_layout_resource(out.parent/'preparation.json')).read_text())
            self.assertEqual(meta['sequence_qc']['internal_stop_sites'],2)
            self.assertEqual(meta['sequence_qc']['report_sha256'],digest(_layout_resource(out.parent/'sequence_qc.tsv')))
            self.assertEqual(meta['input_sha256'],digest(fa))

    def test_ordinary_references_still_reject_internal_stops_and_other_invalid_characters(self):
        with self.assertRaisesRegex(ValueError,r"Invalid sequence: p1.*unsupported characters.*\*"):
            normalize_sequence('MA*GG','p1')
        for sequence in ('MA-GG','MA?GG','*',''):
            with self.subTest(sequence=sequence),self.assertRaisesRegex(ValueError,'Invalid sequence: p1'):
                normalize_sequence(sequence,'p1',internal_stop_policy='mask')
