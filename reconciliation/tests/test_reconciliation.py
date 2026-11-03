"""Synthetic data/provenance regressions; no analytical tools are invoked."""
import contextlib
import importlib.util
import io
import json
import tempfile
import unittest
from argparse import Namespace
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("reconcile_candidates", ROOT / "reconcile_candidates.py")
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)
spec = importlib.util.spec_from_file_location("independent_phylogeny", ROOT.parent / "independent_family_phylogeny/run_phylogeny.py")
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)


def dump(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value), encoding="utf-8")


def write_fasta(path, seqs):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(f">{name}\n{seq}\n" for name, seq in seqs.items()), encoding="ascii")


class HandoffTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="adhesome_reconciliation_test_")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.hr, self.inputs = self.root / "homology/results/run", self.root / "homology"
        self.audit, self.handoff = self.root / "audit", self.root / "handoff"
        self.hr.mkdir(parents=True)
        self.audit.mkdir()
        self.raw = {"g1.1": "ACDEFGHIK", "g1.2": "ACDEFGHIK", "g1.3": "ACDEYGHIK",
                    "g2.1": "ACDEFGHIL", "g3.1": "ACDEFGHIM", "g4.1": "ACDEFGHIN"}
        self.prepared = {f"test__g{n}": self.raw[f"g{n}.1"] for n in range(1, 5)}
        self.references = {f"homo_sapiens|R{n}": "ACDEFGHI" + end for n, end in enumerate("PQRS", 1)}
        self.hits = [f"g{n}.1" for n in range(1, 5)]
        self.main = [{"Specie": "test", "protein_id": f"test__g{n}", "putative_audited_family": "parvin_like",
                      "protein_name": "Putative parvin-like protein", "family_decision": "family_supported" if n < 4 else "architecture_gate_failed",
                      "assignment_confidence": "high", "stage1_family_hypothesis": "parvin_like", "Orthology_Orthogroup": "OG0000818"}
                     for n in range(1, 5)]
        self.main.append({"Specie": "validation", "protein_id": "validation__R1", "family_decision": "family_supported",
                          "putative_audited_family": "parvin_like"})
        self.negative = [{"Specie": "test", "protein_id": "test__g3", "exclusion_reason": "explicit background exclusion"}]
        self.args = Namespace(homology_run=self.hr, homology_inputs=self.inputs, audit=self.audit,
                              handoff=self.handoff, crosswalk=self.root / "crosswalk.json", out=self.root / "library")
        dump(self.args.crosswalk, {"main_accepted_decisions": ["family_supported", "provisional_family"], "families": [
            {"family": "parvin_like", "reference_sets": ["parvin_references"], "main_labels": ["parvin_like"]}]})
        self.refresh()

    def refresh(self):
        target = self.inputs / "target_proteomes/test.faa"
        ref = self.inputs / "reference_fasta/parvin_references.faa"
        write_fasta(target, self.raw)
        write_fasta(ref, self.references)
        folder = self.handoff / "detection/test/prepared"
        write_fasta(folder / "test.faa", self.prepared)
        dump(folder / "preparation.json", {"species": {"test": {"input_sha256": mod.digest(target),
                                                                  "output_sha256": mod.digest(folder / "test.faa")}}})
        maps = [{"protein_id": "test__" + rid.split(".")[0], "original_id": rid,
                 "decision": "retained" if rid.endswith(".1") else "shorter_isoform_or_tie"} for rid in self.raw]
        mod.write_table(folder / "map.tsv", maps, [])
        (folder / "test.id_map.csv").write_text((folder / "map.tsv").read_text().replace("\t", ","))
        rows = [{"reference_set": "parvin_references", "species": "test", "reference_id": list(self.references)[n % len(self.references)],
                 "best_target_id": rid, "status": "RBH_candidate_augmented_reference"} for n, rid in enumerate(self.hits)]
        mapping = [{"reference_set": "parvin_references", "reference_id": rid, "reference_proteome": "Homo_sapiens",
                    "source_assignment": "reference_header_species", "mapping_method": "identical_sequence"} for rid in self.references]
        mod.write_table(self.hr / "homology_audit.tsv", rows, ["reference_set", "species", "reference_id", "best_target_id", "status"])
        mod.write_table(self.hr / "reference_mapping.tsv", mapping, [])
        rbh = {f"test|{rid}": self.raw[rid] for rid in self.hits}
        write_fasta(self.hr / "rbh_candidates.faa", rbh)
        write_fasta(self.hr / "reference_sets/parvin_references/rbh_candidates.faa", rbh)
        dump(self.hr / "run_config.json", {"input_sha256": {str(target): mod.digest(target), str(ref): mod.digest(ref)}})
        dump(self.hr / "run_manifest.json", {"status": "complete", "input_report": {
            "target_proteomes": [{"species": "test", "file": str(target)}],
            "reference_sets": [{"reference_set": "parvin_references", "file": str(ref), "queries": len(self.references)}]},
            "output_sha256": {str(p.relative_to(self.hr)).replace("\\", "/"): mod.digest(p) for p in self.hr.rglob("*") if p.is_file()}})
        mod.write_table(self.audit / "audited_family_assignments.tsv", self.main, [])
        mod.write_table(self.audit / "negative_control_exclusions.tsv", self.negative, ["Specie", "protein_id", "exclusion_reason"])
        dump(self.audit / "workflow_manifest.json", {"kind": "audit", "annotation": str(self.root / "annotation"), "files": {
            name: mod.digest(self.audit / name) for name in ["audited_family_assignments.tsv", "negative_control_exclusions.tsv"]}})

    def test_union_keeps_main_and_rbh_confidence_separate(self):
        data = mod.reconcile(self.args)
        self.assertEqual(data["summary"]["candidate_sources"], {"both": 2, "rbh_only": 1})
        self.assertEqual(data["summary"]["audit_rows_outside_target_cohort"], 1)
        row = next(r for r, s in data["records"] if r["original_ids"] == "g4.1")
        self.assertFalse(row["main_accepted"])
        self.assertEqual(row["main_decision"], "architecture_gate_failed")
        self.assertEqual(row["evidence_source"], "rbh_only")

    def test_explicit_negative_is_excluded_from_both_sources(self):
        data = mod.reconcile(self.args)
        self.assertEqual(data["summary"]["rbh_negative_excluded"], 1)
        self.assertFalse(any("g3" in r["original_ids"] for r, s in data["records"]))
        self.assertEqual({r["source"] for r in data["excluded"]}, {"main", "rbh"})

    def test_identical_alternative_isoform_merges_with_explicit_gene_map(self):
        self.hits.append("g1.2")
        self.refresh()
        data = mod.reconcile(self.args)
        self.assertEqual(len(data["records"]), 3)
        row = next(r for r, s in data["records"] if r["main_protein_id"] == "test__g1")
        self.assertEqual(row["original_ids"], "g1.1;g1.2")

    def test_different_isoform_retained_as_separate_sequence(self):
        self.hits.append("g1.3")
        self.refresh()
        data = mod.reconcile(self.args)
        self.assertEqual(len(data["records"]), 4)
        self.assertEqual(data["conflicts"][0]["relationship"], "alternative_isoform_sequence")
        variant = next(r for r, s in data["records"] if r["original_ids"] == "g1.3")
        self.assertFalse(variant["main_accepted"])

    def test_distinct_genes_with_identical_residues_not_collapsed(self):
        self.raw["g2.1"] = self.raw["g1.1"]
        self.prepared["test__g2"] = self.raw["g2.1"]
        self.refresh()
        self.assertEqual(len(mod.reconcile(self.args)["records"]), 3)

    def test_family_conflict_visible_without_overwriting_main_name(self):
        self.main[0]["putative_audited_family"] = "filamin"
        self.refresh()
        data = mod.reconcile(self.args)
        row = next(r for r, s in data["records"] if r["original_ids"] == "g1.1")
        self.assertTrue(row["family_conflict"])
        self.assertEqual(row["main_family"], "filamin")
        self.assertEqual(row["phylogeny_families"], "parvin_like")

    def test_main_only_calls_retained_when_rbh_has_no_hits(self):
        self.hits = []
        self.refresh()
        data = mod.reconcile(self.args)
        self.assertEqual(data["summary"]["candidate_sources"], {"main_only": 2})
        self.assertEqual(len(data["prepared_families"]), 1)

    def test_provisional_main_call_preserves_its_tier(self):
        self.main[0]["family_decision"] = "provisional_family"
        self.refresh()
        row = next(r for r, s in mod.reconcile(self.args)["records"] if r["original_ids"] == "g1.1")
        self.assertEqual(row["main_decision"], "provisional_family")

    def test_modified_saved_inputs_are_rejected(self):
        with (self.inputs / "target_proteomes/test.faa").open("a") as handle:
            handle.write("\n")
        with self.assertRaisesRegex(ValueError, "source hash changed"):
            mod.reconcile(self.args)

    def test_saved_rbh_membership_and_residues_are_verified(self):
        path = self.hr / "rbh_candidates.faa"
        path.write_text(">test|g1.1\nACDEFGHIK\n")
        manifest = mod.load_json(self.hr / "run_manifest.json")
        manifest["output_sha256"]["rbh_candidates.faa"] = mod.digest(path)
        dump(self.hr / "run_manifest.json", manifest)
        with self.assertRaisesRegex(ValueError, "Global RBH FASTA does not match"):
            mod.reconcile(self.args)

    def test_unknown_family_requires_explicit_crosswalk(self):
        config = mod.load_json(self.args.crosswalk)
        config["families"][0]["reference_sets"] = ["other"]
        dump(self.args.crosswalk, config)
        with self.assertRaisesRegex(ValueError, "explicit family crosswalk"):
            mod.reconcile(self.args)

    def test_handoff_auto_resolution_from_annotation(self):
        dump(self.root / "annotation/workflow_manifest.json", {"handoff": str(self.handoff)})
        self.args.handoff = None
        self.assertEqual(mod.reconcile(self.args)["summary"]["reconciled_candidates"], 3)

    def test_prepare_is_compatible_and_refuses_overwrite(self):
        data = mod.reconcile(self.args)
        with contextlib.redirect_stdout(io.StringIO()):
            mod.prepare(self.args.out, data)
            manifest = mod.check(self.args.out)
        self.assertEqual(manifest["families"][0]["records"], 7)
        self.assertTrue(manifest["preserve_identical_tips"])
        with self.assertRaises(FileExistsError):
            mod.prepare(self.args.out, data)
        (self.args.out / "parvin_like/input.fasta").write_text(">changed\nACDE\n")
        with self.assertRaisesRegex(ValueError, "Prepared file changed"):
            mod.check(self.args.out)

    def test_plan_writes_nothing(self):
        before = {p.relative_to(self.root): mod.digest(p) for p in self.root.rglob("*") if p.is_file()}
        with contextlib.redirect_stdout(io.StringIO()):
            code = mod.main(["--homology-run", str(self.hr), "--homology-inputs", str(self.inputs),
                             "--audit", str(self.audit), "--handoff", str(self.handoff), "--crosswalk", str(self.args.crosswalk),
                             "--out", str(self.args.out)])
        self.assertEqual(code, 0)
        after = {p.relative_to(self.root): mod.digest(p) for p in self.root.rglob("*") if p.is_file()}
        self.assertEqual(before, after)

    def test_reference_only_families_are_not_prepared_as_target_trees(self):
        self.hits = []
        for row in self.main:
            row["family_decision"] = "architecture_gate_failed"
        self.refresh()
        data = mod.reconcile(self.args)
        self.assertEqual(data["prepared_families"], [])
        self.assertEqual(data["family_summary"][0]["status"], "not_prepared_no_target_candidates")
        with contextlib.redirect_stdout(io.StringIO()):
            mod.prepare(self.args.out, data)
            mod.check(self.args.out)

    def test_prepared_library_is_independent_of_source_installations(self):
        data = mod.reconcile(self.args)
        with contextlib.redirect_stdout(io.StringIO()):
            mod.prepare(self.args.out, data)
        self.audit.rename(self.root / "moved_audit")
        self.handoff.rename(self.root / "moved_handoff")
        self.inputs.rename(self.root / "moved_homology")
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(mod.check(self.args.out)["summary"]["reconciled_candidates"], 3)

    def test_source_change_before_export_is_rejected_without_output(self):
        data = mod.reconcile(self.args)
        with (self.audit / "audited_family_assignments.tsv").open("a") as handle:
            handle.write("\n")
        with self.assertRaisesRegex(ValueError, "Source changed during reconciliation"):
            mod.prepare(self.args.out, data)
        self.assertFalse(self.args.out.exists())

    def test_accepted_main_residue_mismatch_is_not_silently_replaced(self):
        self.prepared["test__g1"] = "ACDEFGHIV"
        self.refresh()
        with self.assertRaisesRegex(ValueError, "Accepted main sequence disagrees"):
            mod.reconcile(self.args)

    def test_output_inside_raw_inputs_is_refused(self):
        with contextlib.redirect_stderr(io.StringIO()):
            code = mod.main(["--prepare", "--homology-run", str(self.hr), "--homology-inputs", str(self.inputs),
                             "--audit", str(self.audit), "--handoff", str(self.handoff), "--crosswalk", str(self.args.crosswalk),
                             "--out", str(self.inputs / "target_proteomes/should_not_exist")])
        self.assertEqual(code, 1)
        self.assertFalse((self.inputs / "target_proteomes/should_not_exist").exists())


class RunnerTests(unittest.TestCase):
    def test_reconciled_plan_retains_identical_tips_and_old_plan_unchanged(self):
        with tempfile.TemporaryDirectory(prefix="phylogeny_plan_test_") as tmp:
            root = Path(tmp)
            write_fasta(root / "library/demo/input.fasta", {f"C{n:06d}": "ACDEFG" + aa for n, aa in enumerate("HIKL", 1)})
            args = Namespace(exclude_panel=[], mafft_threads=4, alignment="genafpair", trim="gappyout",
                             iqtree_threads="4", seed=123, stop_after="tree", execute=False, reconciled_library=True)
            stdout = io.StringIO()
            with contextlib.redirect_stdout(stdout):
                runner.run_family(args, "demo", root / "library", root / "results", dict(mafft="mafft", trimal="trimal", iqtree="iqtree"))
            self.assertIn("-keep-ident", stdout.getvalue())
            self.assertIn("--ep 0 --genafpair --maxiterate 1000", stdout.getvalue())
            self.assertIn("-gappyout", stdout.getvalue())
            self.assertFalse((root / "results").exists())
            args.reconciled_library = False
            stdout = io.StringIO()
            with contextlib.redirect_stdout(stdout):
                runner.run_family(args, "demo", root / "library", root / "results", dict(mafft="mafft", trimal="trimal", iqtree="iqtree"))
            self.assertNotIn("-keep-ident", stdout.getvalue())

    def test_tree_validator_checks_exact_tip_membership(self):
        with tempfile.TemporaryDirectory(prefix="newick_test_") as tmp:
            path = Path(tmp) / "tree.treefile"
            path.write_text("(C000001:0.1,(R000001:0.2,R000002:0.3)99:0.4);")
            runner.check_tree_membership(path, {"C000001", "R000001", "R000002"})
            with self.assertRaisesRegex(ValueError, "membership"):
                runner.check_tree_membership(path, {"C000001", "R000001", "R000003"})

    def test_trimming_must_keep_global_columns(self):
        with tempfile.TemporaryDirectory(prefix="trim_test_") as tmp:
            aligned, trimmed = Path(tmp) / "aligned.fasta", Path(tmp) / "trimmed.fasta"
            aligned.write_text(">C000001\nA-CDE\n>R000001\nAB-DE\n")
            trimmed.write_text(">C000001\nADE\n>R000001\nADE\n")
            runner.check_trimmed_columns(trimmed, aligned)
            trimmed.write_text(">C000001\nACDE\n>R000001\nABDE\n")
            with self.assertRaisesRegex(ValueError, "columns or residues"):
                runner.check_trimmed_columns(trimmed, aligned)


if __name__ == "__main__":
    unittest.main()
