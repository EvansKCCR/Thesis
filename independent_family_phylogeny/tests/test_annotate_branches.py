"""Synthetic checks only; no saved family trees or analytical tools are used."""
import importlib.util
import tempfile
import unittest
from pathlib import Path

MODULE = Path(__file__).resolve().parents[1] / "annotate_branches.py"
spec = importlib.util.spec_from_file_location("annotate_branches", MODULE)
annotation = importlib.util.module_from_spec(spec)
spec.loader.exec_module(annotation)


class Node:
    def __init__(self, name=None, children=(), confidence=None):
        self.name, self.clades, self.confidence = name, list(children), confidence

    def get_terminals(self):
        return [self] if not self.clades else [tip for child in self.clades for tip in child.get_terminals()]


class Tree:
    def __init__(self):
        self.root = Node(children=[Node(children=[Node("a"), Node("b")], confidence=90),
                                   Node(children=[Node("c"), Node("d")], confidence=90)])

    def get_nonterminals(self, order="preorder"):
        return [self.root] + self.root.clades


class BranchEvidenceTests(unittest.TestCase):
    def branch_rows(self, rooted=False, events=()):
        return annotation.branch_rows("example", Tree(), set("abcd"),
                                      {"a": "Species1", "b": "Species2", "c": "Species3", "d": "validation"},
                                      dict.fromkeys("abcd", "OG1"),
                                      [{"sequence_id_1": "a", "sequence_id_2": "b"},
                                       {"sequence_id_1": "a", "sequence_id_2": "c"}], list(events), rooted)

    def event(self, first, second, name="D1"):
        return {"genes_1": set(first), "genes_2": set(second), "event_id": name}

    def test_unrooted_split_deduplicated_and_counts(self):
        result = self.branch_rows()
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0]["tip_ids"], "a,b")
        self.assertEqual(result[0]["support_label"], "90")
        self.assertEqual(result[0]["direct_orthology_edges_within"], 1)
        self.assertEqual(result[0]["orthology_edges_across_split"], 1)

    def test_unrooted_compatibility_is_not_exact_duplication_assignment(self):
        row = self.branch_rows(events=[self.event("ab", "cd")])[0]
        self.assertEqual(row["duplication_split_compatible_rows"], "D1")
        self.assertEqual(row["exact_duplication_rows"], "")

    def test_rooted_exact_child_sets_only(self):
        result = self.branch_rows(rooted=True, events=[self.event("a", "b"), self.event("a", "bc", "D2")])
        self.assertEqual(len(result), 2)
        self.assertEqual(result[0]["exact_duplication_rows"], "D1")

    def test_partial_event_never_counts_as_complete_match(self):
        row = self.branch_rows(events=[self.event("abx", "cd")])[0]
        self.assertEqual(row["duplication_split_compatible_rows"], "")
        self.assertEqual(row["incomplete_duplication_rows"], "D1")

    def test_control_panel_not_reported_as_species(self):
        row = self.branch_rows(rooted=True)[1]
        self.assertEqual(row["species"], "Species3")
        self.assertEqual(row["mixed_control_tip_ids"], "d")

    def test_explicit_identifier_aliases_and_no_suffix_guessing(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "WorkingDirectory").mkdir()
            (root / "WorkingDirectory/SequenceIDs.txt").write_text("0_0: Species1__a\n")
            mapper = annotation.Identifiers({"OG1": {"Species1__a"}}, {"Species1__a": "Species1"}, root)
            for alias in ("Species1__a", "Species1_Species1__a", "0_0"):
                self.assertEqual(mapper.resolve(alias, "OG1", "fixture"), {"Species1__a"})
            with self.assertRaises(ValueError):
                mapper.resolve("a", "OG1", "fixture")

    def test_native_full_evidence_duplicate_pairs_and_nonempty_events(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            table = root / "Orthologues/Orthologues_A/A__v__B.tsv"
            table.parent.mkdir(parents=True)
            table.write_text("Orthogroup\tA\tB\nOG1\tA__a\tB__b, B__c\nOG1\tA__a\tB__b\n")
            duplication = root / "Gene_Duplication_Events/Duplications.tsv"
            duplication.parent.mkdir()
            duplication.write_text("Orthogroup\tGenes 1\tGenes 2\tSupport\nOG1\tA_A__a\tB_B__b, B_B__c\t0.8\n")
            mapper = annotation.Identifiers({"OG1": {"A__a", "B__b", "B__c"}},
                                            {"A__a": "A", "B__b": "B", "B__c": "B"}, root)
            edges, events, sources = annotation.evidence(root, mapper, {"OG1"}, {"A", "B"})
            self.assertEqual(len(edges), 2)
            self.assertEqual(edges[0]["source_rows"], "2,3")
            self.assertEqual(events[0]["genes_1"], {"A__a"})
            self.assertEqual(events[0]["genes_2"], {"B__b", "B__c"})
            self.assertEqual(len(sources), 2)

    def test_changed_tree_membership_rejected(self):
        try:
            from Bio import Phylo  # noqa: F401
        except ImportError:
            self.skipTest("Biopython is not installed in this test environment")
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "tree.treefile").write_text("((a,b)90,(c,d)90);\n")
            with self.assertRaises(ValueError):
                annotation.load_tree(root, "tree.treefile", set("abce"))


if __name__ == "__main__":
    unittest.main()
