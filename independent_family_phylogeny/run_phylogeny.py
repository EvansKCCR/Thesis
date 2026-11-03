#!/usr/bin/env python3
"""Standalone family phylogeny. Planning is default; --execute launches tools."""
import argparse
import hashlib
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def fasta(path):
    records, name, chunks = {}, None, []
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                if name in records:
                    raise ValueError(f"Duplicate identifier: {name}")
                records[name] = "".join(chunks).upper()
            name, chunks = line[1:].split()[0], []
        elif line.strip():
            if name is None:
                raise ValueError(f"Sequence before header: {path}")
            chunks.append(line.strip())
    if name is not None:
        if name in records:
            raise ValueError(f"Duplicate identifier: {name}")
        records[name] = "".join(chunks).upper()
    if not records or any(not seq.replace("-", "").replace(".", "") for seq in records.values()):
        raise ValueError(f"Empty FASTA or sequence with no residues: {path}")
    return records


def check_alignment(path, expected, preserve_residues=False):
    actual = fasta(path)
    if set(actual) != set(expected):
        raise ValueError(f"Alignment changed sequence identifiers/membership: {path}")
    if len({len(seq) for seq in actual.values()}) != 1:
        raise ValueError(f"Unequal alignment row lengths: {path}")
    if preserve_residues and any(seq.replace("-", "").replace(".", "") != expected[name]
                                 for name, seq in actual.items()):
        raise ValueError(f"Alignment changed residues: {path}")


def check_trimmed_columns(path, aligned):
    """A reconciled library must retain the same global alignment columns."""
    source, trimmed = fasta(aligned), fasta(path)
    names = list(source)
    column = 0
    for n in range(len(next(iter(trimmed.values())))):
        wanted = tuple(trimmed[name][n].replace(".", "-") for name in names)
        while column < len(source[names[0]]):
            actual = tuple(source[name][column].replace(".", "-") for name in names)
            column += 1
            if actual == wanted:
                break
        else:
            raise ValueError(f"Trimming changed alignment columns or residues: {path}")


def check_tree_membership(path, expected):
    text = path.read_text().strip()
    if not text.endswith(";"):
        raise ValueError(f"Incomplete Newick tree: {path}")
    leaves = [quoted.replace("''", "'") if quoted else plain for quoted, plain in re.findall(
        r"(?:\(|,)\s*(?:'((?:[^']|'')+)'|([^():,;\s]+))\s*(?=:|,|\))", text)]
    if len(leaves) != len(expected) or set(leaves) != set(expected):
        raise ValueError(f"Tree changed tip identifiers/membership: {path}")


def tool(name, execute):
    candidates = [name] if name != "iqtree" else ["iqtree", "iqtree2", "iqtree3"]
    for candidate in candidates:
        found = shutil.which(candidate)
        if found:
            return str(Path(found).resolve())
    if execute:
        raise ValueError(f"Required executable is not on PATH: {name}")
    return name


def dump(path, value):
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


def canonical_command(name, argv):
    """Accept old completed tree receipts without repeating their unsafe --redo."""
    return [part for part in argv if part not in ("-redo", "--redo")] if name == "tree" else list(argv)


def matching_fingerprint(name, saved, current):
    previous = dict(saved)
    previous["command"] = canonical_command(name, previous["command"])
    return previous == current


def stage(folder, name, argv, inputs, outputs, validator, stdout_path=None):
    argv = canonical_command(name, argv)
    checkpoint = folder / f"{name}.complete.json"
    fingerprint = {"command": argv,
                   "input_sha256": {p.name: digest(p) for p in inputs},
                   "executable_sha256": digest(argv[0])}
    if checkpoint.exists():
        saved = json.loads(checkpoint.read_text())
        if not matching_fingerprint(name, saved["fingerprint"], fingerprint):
            raise ValueError(f"{name}: input/options/tool changed; choose a new --run-name")
        if any(not p.is_file() or saved["output_sha256"].get(p.name) != digest(p) for p in outputs):
            raise ValueError(f"{name}: completed output changed/missing; choose a new --run-name")
        validator()
        print(f"REUSE {folder.name}/{name}", flush=True)
        return
    command_path = folder / f"{name}.command.json"
    prior = json.loads(command_path.read_text()) if command_path.exists() else None
    if prior and canonical_command(name, prior["argv"]) != argv:
        raise ValueError(f"{name}: unfinished command settings/tool path changed; restore the original environment/options")
    if name == "tree":
        native = [folder / "tree.ckp.gz", folder / "tree.model.gz"]
        if any(path.exists() for path in native) and prior is None:
            raise ValueError("tree: native checkpoint exists without original command provenance; resume refused")
        for path in native:
            if path.exists():
                print(f"RESUME {folder.name}/tree: reuse {path.name}", flush=True)
    started = datetime.now(timezone.utc).isoformat()
    attempts = list(prior.get("attempts", [{"argv": prior["argv"], "started_utc": prior.get("started_utc")}])) if prior else []
    attempts.append({"argv": argv, "started_utc": started})
    dump(command_path, {"argv": argv, "cwd": str(folder), "started_utc": started, "attempts": attempts})
    print(shlex.join(argv) + (" > " + shlex.quote(str(stdout_path)) if stdout_path else ""), flush=True)
    stdout_log = stdout_path or folder / f"{name}.stdout.log"
    log_mode = "a" if name == "tree" and prior else "w"
    with stdout_log.open(log_mode) as out, (folder / f"{name}.stderr.log").open(log_mode) as err:
        process = subprocess.run(argv, cwd=folder, stdout=out, stderr=err, check=False)
    if process.returncode:
        raise ValueError(f"{name} exited {process.returncode}; inspect {folder / (name + '.stderr.log')}")
    validator()
    if any(not p.is_file() or p.stat().st_size == 0 for p in outputs):
        raise ValueError(f"{name} did not produce expected nonempty outputs")
    dump(checkpoint, {"completed_utc": datetime.now(timezone.utc).isoformat(),
                      "fingerprint": fingerprint,
                      "output_sha256": {p.name: digest(p) for p in outputs}})


def run_family(args, family, library, out, binaries):
    source = library / family / "input.fasta"
    expected = fasta(source)
    excluded = set(args.exclude_panel)
    # Panel labels are explicit in the manifest, rather than inferred from gene names.
    if excluded:
        import csv
        with (library / family / "tips.tsv").open(newline="") as handle:
            omitted = {row["sequence_id"] for row in csv.DictReader(handle, delimiter="\t")
                       if row["panel_label"] in excluded}
        expected = {name: seq for name, seq in expected.items() if name not in omitted}
    if len(expected) < 4 or len(set(expected.values())) < 4:
        raise ValueError(f"{family}: fewer than four taxa or four distinct sequences after selection")
    folder = out / family
    selected = folder / "selected_input.fasta"
    aligned, trimmed = folder / "aligned.fasta", folder / "trimmed.fasta"
    align_cmd = [binaries["mafft"], "--thread", str(args.mafft_threads)]
    align_cmd += ["--ep", "0", "--genafpair"] if args.alignment == "genafpair" else ["--localpair"]
    align_cmd += ["--maxiterate", "1000", str(selected)]
    trim_cmd = [binaries["trimal"], "-in", str(aligned), "-out", str(trimmed)]
    trim_cmd += ["-gappyout"] if args.trim == "gappyout" else ["-gt", "0.2"]
    tree_cmd = [binaries["iqtree"], "-s", str(trimmed), "-m", "MFP", "-B", "1000",
                "-T", args.iqtree_threads, "-seed", str(args.seed), "--prefix", str(folder / "tree")]
    # The new handoff uses short mapped IDs and explicitly retains identical tips.
    # Existing orthogroup libraries retain their original command/checkpoints.
    reconciled = getattr(args, "reconciled_library", False)
    if reconciled:
        tree_cmd += ["-st", "AA", "-keep-ident"]
    commands = [("align", align_cmd), ("trim", trim_cmd), ("tree", tree_cmd)]
    stop = ["align", "trim", "tree"].index(args.stop_after) + 1
    if not args.execute:
        print(f"\n{family}: {len(expected)} records, {len(set(expected.values()))} unique sequences")
        for name, argv in commands[:stop]:
            print(shlex.join(argv) + (" > " + shlex.quote(str(aligned)) if name == "align" else ""))
        return
    folder.mkdir(parents=True, exist_ok=True)
    # Freeze run-level settings even if the first external command fails.
    selection = {"family": family, "library_input_sha256": digest(source),
                 "excluded_panels": sorted(excluded), "selected_ids": list(expected),
                 "alignment": args.alignment, "trim": args.trim,
                 "mafft_threads": args.mafft_threads, "iqtree_threads": args.iqtree_threads,
                 "seed": args.seed}
    config_path = folder / "run_config.json"
    if config_path.exists() and json.loads(config_path.read_text()) != selection:
        raise ValueError(f"{family}: settings changed; choose a new --run-name")
    if not config_path.exists():
        dump(config_path, selection)
    text = "".join(f">{name}\n{seq}\n" for name, seq in expected.items())
    if selected.exists() and selected.read_text() != text:
        raise ValueError(f"{family}: selected input changed; choose a new --run-name")
    if not selected.exists():
        selected.write_text(text, encoding="ascii")
    stage(folder, "align", align_cmd, [selected], [aligned],
          lambda: check_alignment(aligned, expected, True), stdout_path=aligned)
    if stop < 2:
        return
    def validate_trim():
        check_alignment(trimmed, expected)
        if reconciled:
            check_trimmed_columns(trimmed, aligned)
    stage(folder, "trim", trim_cmd, [aligned], [trimmed], validate_trim)
    if stop < 3:
        return
    if len({seq.replace("-", "").replace(".", "") for seq in fasta(trimmed).values()}) < 4:
        raise ValueError(f"{family}: trimming leaves fewer than four distinct sequences")
    def validate_tree():
        if reconciled:
            for name in ("tree.treefile", "tree.contree"):
                check_tree_membership(folder / name, expected)
    stage(folder, "tree", tree_cmd, [trimmed],
          [folder / "tree.treefile", folder / "tree.iqtree", folder / "tree.contree"], validate_tree)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--execute", action="store_true", help="Actually launch MAFFT, trimAl and IQ-TREE")
    parser.add_argument("--dry-run",action="store_true",help="Structural planning without a prepared library or biological data")
    parser.add_argument("--clean-state",action="store_true")
    parser.add_argument("--report-dir",type=Path)
    parser.add_argument("--resume", action="store_true", help="Require an existing run; reuse completed stages and native IQ-TREE checkpoints")
    parser.add_argument("--library", type=Path, default=ROOT / "library")
    parser.add_argument("--run-name", default="genafpair_gappyout")
    parser.add_argument("--family", action="append", help="Repeat to select families; default: all prepared families")
    parser.add_argument("--alignment", choices=["genafpair", "localpair"], default="genafpair")
    parser.add_argument("--trim", choices=["gappyout", "gt02"], default="gappyout")
    parser.add_argument("--mafft-threads", type=int, default=4)
    parser.add_argument("--iqtree-threads", default="AUTO", help="AUTO or positive integer")
    parser.add_argument("--seed", type=int, default=20261028)
    parser.add_argument("--stop-after", choices=["align", "trim", "tree"], default="tree")
    parser.add_argument("--exclude-panel", action="append", default=[],
                        help="Optional explicit panel exclusion, e.g. validation; default: retain all")
    args = parser.parse_args(argv)
    try:
        if args.dry_run:
            if args.execute:raise ValueError('Choose --dry-run or --execute, not both')
            sys.path.insert(0,str(ROOT.parent/'integrin_adhesome_analysis/resource_library/scripts'))
            from deployment_validation import main as validate_deployment
            flags=['--dry-run','--workflow','phylogeny','--library',str(args.library)]
            if args.clean_state:flags.append('--clean-state')
            if args.report_dir:flags+=['--report-dir',str(args.report_dir)]
            return validate_deployment(flags)
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]*", args.run_name):
            raise ValueError("--run-name must contain only letters, numbers, underscores or hyphens")
        if args.mafft_threads < 1 or (args.iqtree_threads != "AUTO" and
                                     (not args.iqtree_threads.isdigit() or int(args.iqtree_threads) < 1)):
            raise ValueError("Thread counts must be positive integers (IQ-TREE also accepts AUTO)")
        if args.execute and os.name == "nt":
            raise ValueError("Execute this bundle inside WSL/Linux with the bioinformatics tools on PATH")
        library = args.library.resolve()
        manifest = json.loads((library / "library_manifest.json").read_text())
        args.reconciled_library = manifest.get("library_kind") == "reconciled_candidates"
        families = [row["family"] for row in manifest["families"]]
        requested = list(dict.fromkeys(args.family or families))
        if set(requested) - set(families):
            raise ValueError(f"Unknown families: {set(requested) - set(families)}")
        import csv
        panels = set()
        for family in requested:
            for name in ("input.fasta", "tips.tsv"):
                path = library / family / name
                key = f"{family}/{name}"
                if digest(path) != manifest["files"][key]:
                    raise ValueError(f"Prepared library changed: {path}")
            with (library / family / "tips.tsv").open(newline="") as handle:
                panels.update(row["panel_label"] for row in csv.DictReader(handle, delimiter="\t"))
        if set(args.exclude_panel) - panels:
            raise ValueError(f"Unknown excluded panels in selected families: {set(args.exclude_panel) - panels}")
        binaries = {name: tool(name, args.execute) for name in ("mafft", "trimal", "iqtree")}
        out = ROOT / "results" / args.run_name
        if args.resume and not out.is_dir():
            raise ValueError(f"--resume requires an existing run directory: {out}")
        failures = []
        for family in requested:
            try:
                run_family(args, family, library, out, binaries)
            except (OSError, ValueError, KeyError) as error:
                failures.append({"family": family, "error": str(error)})
                print(f"ERROR {family}: {error}", file=sys.stderr, flush=True)
        if args.execute:
            out.mkdir(parents=True, exist_ok=True)
            dump(out / "run_status.json", {"updated_utc": datetime.now(timezone.utc).isoformat(),
                                           "requested_families": requested, "stop_after": args.stop_after,
                                           "failures": failures,
                                           "status": "failed" if failures else "requested_stages_complete"})
        else:
            print("\nPLAN ONLY: no commands executed; no results written. Add --execute to run locally.")
        return 1 if failures else 0
    except (OSError, ValueError, KeyError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
