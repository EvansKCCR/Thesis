"""Synthetic regression checks for checkpoint recovery; never launch tools."""
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

MODULE = Path(__file__).resolve().parents[1] / "run_phylogeny.py"
spec = importlib.util.spec_from_file_location("family_run_resume", MODULE)
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)


class ResumeTests(unittest.TestCase):
    def test_legacy_redo_receipt_matches_without_repeating_redo(self):
        current = {"command": ["iqtree2", "-m", "MFP"], "input_sha256": {"trimmed.fasta": "abc"},
                   "executable_sha256": "tool"}
        previous = dict(current, command=current["command"] + ["--redo"])
        self.assertTrue(runner.matching_fingerprint("tree", previous, current))
        previous["input_sha256"] = {"trimmed.fasta": "changed"}
        self.assertFalse(runner.matching_fingerprint("tree", previous, current))

    def test_changed_inference_settings_remain_rejected(self):
        previous = {"command": ["iqtree2", "-m", "LG", "--redo"]}
        self.assertFalse(runner.matching_fingerprint("tree", previous,
                                                   {"command": ["iqtree2", "-m", "MFP"]}))

    def test_completed_legacy_tree_is_reused_without_launch(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            binary, source, output = root / "fake_iqtree", root / "trimmed.fasta", root / "tree.treefile"
            binary.write_text("synthetic binary placeholder")
            source.write_text(">a\nAAAA\n")
            output.write_text("(a,b,c);\n")
            argv = [str(binary), "-s", str(source)]
            runner.dump(root / "tree.complete.json", {
                "fingerprint": {"command": argv + ["--redo"],
                                "input_sha256": {source.name: runner.digest(source)},
                                "executable_sha256": runner.digest(binary)},
                "output_sha256": {output.name: runner.digest(output)}})
            with patch.object(runner.subprocess, "run", side_effect=AssertionError("No tool launch allowed")):
                runner.stage(root, "tree", argv, [source], [output], lambda: None)
            self.assertEqual(json.loads((root / "tree.complete.json").read_text())["fingerprint"]["command"][-1], "--redo")

    def test_unfinished_resume_keeps_model_checkpoint_and_appends_logs(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            binary, source, output = root / "fake_iqtree", root / "trimmed.fasta", root / "tree.treefile"
            binary.write_text("synthetic binary placeholder")
            source.write_text(">a\nAAAA\n")
            model = root / "tree.model.gz"
            model.write_bytes(b"synthetic saved model data")
            argv = [str(binary), "-s", str(source)]
            runner.dump(root / "tree.command.json", {"argv": argv + ["--redo"], "started_utc": "old"})
            (root / "tree.stdout.log").write_text("old log\n")

            def fake_tool(command, **kwargs):
                self.assertNotIn("--redo", command)
                self.assertNotIn("-redo", command)
                self.assertEqual(model.read_bytes(), b"synthetic saved model data")
                kwargs["stdout"].write("synthetic continuation\n")
                output.write_text("(a,b,c);\n")
                return type("Result", (), {"returncode": 0})()

            with patch.object(runner.subprocess, "run", side_effect=fake_tool):
                runner.stage(root, "tree", argv, [source], [output], lambda: None)
            self.assertEqual((root / "tree.stdout.log").read_text(), "old log\nsynthetic continuation\n")
            self.assertEqual(len(json.loads((root / "tree.command.json").read_text())["attempts"]), 2)
            self.assertEqual(model.read_bytes(), b"synthetic saved model data")


if __name__ == "__main__":
    unittest.main()
