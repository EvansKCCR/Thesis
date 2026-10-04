import tempfile
import unittest
from pathlib import Path
from prepare_reference import prepare


class ReferencePreparationTests(unittest.TestCase):
    def test_collapse_and_mapping_failure(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            fa=root/'raw.faa'
            fa.write_text('>p2\nMAAA*\n>p1\nMAAA\n>p3\nMA\n')
            mapping=root/'map.tsv'
            mapping.write_text('protein_id\tgene_id\np2\tg1\np1\tg1\np3\tg2\n')
            out=prepare('RefA',fa,mapping,root/'prepared')
            self.assertEqual(out.read_text(),'>RefA__g1\nMAAA\n>RefA__g2\nMA\n')
            self.assertIn('p1\tg1\tRefA__g1\tp1\t1', (root/'prepared/identifier_map.tsv').read_text())
            mapping.write_text('protein_id\tgene_id\np1\tg1\n')
            with self.assertRaisesRegex(ValueError,'missing gene mapping'):
                prepare('RefA',fa,mapping,root/'bad')
            self.assertFalse((root/'bad').exists())
