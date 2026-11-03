#!/usr/bin/env python3
"""Annotate existing family trees using the corresponding full OrthoFinder evidence.

No alignment, inference, rooting, reconciliation, subprocesses or workflow imports.
Requires Biopython. Output is written separately; existing output is never replaced.
"""
import argparse
import csv
import hashlib
import json
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def digest(path):
    value = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def rows(path, required=()):
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or set(required) - set(reader.fieldnames):
            raise ValueError(f"Missing required columns {required}: {path}")
        for line, row in enumerate(reader, 2):
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f"Malformed TSV row: {path}:{line}")
            yield line, row


def members(value):
    return [item.strip() for item in value.split(",") if item.strip()]


def write_tsv(path, records, fields):
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)


def write_json(path, data):
    with path.open("x", encoding="utf-8") as handle:
        json.dump(data, handle, indent=2)
        handle.write("\n")


class Identifiers:
    """Resolve only explicit, unique aliases, never an arbitrary suffix match."""
    def __init__(self, by_group, panel_by_id, results):
        self.by_group = by_group
        self.panel_by_id = panel_by_id
        self.aliases = defaultdict(set)
        for group_ids in by_group.values():
            for protein_id in group_ids:
                panel = panel_by_id[protein_id]
                self.aliases[protein_id].add(protein_id)
                # Full OrthoFinder duplication files can add species + '_' again.
                self.aliases[panel + "_" + protein_id].add(protein_id)
        sequence_ids = results / "WorkingDirectory/SequenceIDs.txt"
        self.mapping_file = sequence_ids if sequence_ids.exists() else None
        if self.mapping_file:
            with sequence_ids.open(encoding="utf-8-sig") as handle:
                for line in handle:
                    numeric_id, sep, original_id = line.strip().partition(": ")
                    if sep and original_id in panel_by_id:
                        self.aliases[numeric_id].add(original_id)

    def resolve(self, value, group, source):
        resolved = set()
        for raw in members(value):
            candidates = self.aliases.get(raw, set()) & self.by_group[group]
            if len(candidates) != 1:
                raise ValueError(f"Unresolved/ambiguous evidence identifier {raw!r}, {group}, {source}")
            resolved.add(next(iter(candidates)))
        if not resolved:
            raise ValueError(f"Empty gene set in evidence: {source}")
        return resolved


def source_membership(results, selected_groups):
    by_group, panels = {}, {}
    table = results / "Orthogroups/Orthogroups.tsv"
    for _, row in rows(table, ["Orthogroup"]):
        group = row["Orthogroup"]
        if group not in selected_groups:
            continue
        if group in by_group:
            raise ValueError(f"Duplicate source orthogroup row: {group}")
        ids = set()
        for panel, value in row.items():
            if panel == "Orthogroup":
                continue
            for protein_id in members(value):
                if protein_id in panels:
                    raise ValueError(f"Nonunique source sequence identifier: {protein_id}")
                ids.add(protein_id)
                panels[protein_id] = panel
        by_group[group] = ids
    if set(by_group) != selected_groups:
        raise ValueError(f"Missing source orthogroups: {selected_groups - set(by_group)}")
    return by_group, panels


def pair_tables(results):
    """Choose one direction per species pair from native OrthoFinder pair tables."""
    chosen = {}
    for path in sorted((results / "Orthologues").glob("Orthologues_*/*.tsv")):
        with path.open(encoding="utf-8-sig", newline="") as handle:
            header = next(csv.reader(handle, delimiter="\t"), [])
        if len(header) != 3 or header[0] != "Orthogroup" or header[1] == header[2]:
            raise ValueError(f"Unexpected native pair-table header: {path}")
        chosen.setdefault(tuple(sorted(header[1:])), (path, header[1:]))
    if not chosen:
        raise ValueError("No native Orthologues/Orthologues_*/*.tsv tables found; "
                         "candidate-only exports are deliberately not substituted")
    return list(chosen.values())


def evidence(results, mapping, selected_groups, needed_panels):
    edges, events, sources = {}, [], []
    for path, panel_columns in pair_tables(results):
        if not set(panel_columns) <= needed_panels:
            continue
        relative = path.relative_to(results).as_posix()
        matching = 0
        for line, row in rows(path, ["Orthogroup"] + panel_columns):
            group = row["Orthogroup"]
            if group not in selected_groups:
                continue
            matching += 1
            source = f"{relative}:{line}"
            left = mapping.resolve(row[panel_columns[0]], group, source)
            right = mapping.resolve(row[panel_columns[1]], group, source)
            if any(mapping.panel_by_id[p] != panel_columns[0] for p in left) or any(
                    mapping.panel_by_id[p] != panel_columns[1] for p in right):
                raise ValueError(f"Pair-table species/identifier disagreement: {source}")
            relation = "one-to-one" if len(left) == len(right) == 1 else (
                "one-to-many" if min(len(left), len(right)) == 1 else "many-to-many")
            for first in left:
                for second in right:
                    a, b = sorted((first, second))
                    key = (group, a, b)
                    if key not in edges:
                        edges[key] = {"orthogroup": group, "sequence_id_1": a, "sequence_id_2": b,
                                      "relationship": set(), "source_file": relative,
                                      "source_rows": set()}
                    edges[key]["relationship"].add(relation)
                    edges[key]["source_rows"].add(str(line))
        sources.append({"path": relative, "sha256": digest(path), "selected_group_rows": matching})
    duplication_file = results / "Gene_Duplication_Events/Duplications.tsv"
    for line, row in rows(duplication_file, ["Orthogroup", "Genes 1", "Genes 2", "Support"]):
        group = row["Orthogroup"]
        if group not in selected_groups:
            continue
        source = f"Gene_Duplication_Events/Duplications.tsv:{line}"
        left = mapping.resolve(row["Genes 1"], group, source)
        right = mapping.resolve(row["Genes 2"], group, source)
        if left & right:
            raise ValueError(f"Duplication child gene sets overlap: {source}")
        events.append({"event_id": f"D{line:08d}", "orthogroup": group,
                       "species_tree_node": row.get("Species Tree Node", ""),
                       "gene_tree_node": row.get("Gene Tree Node", ""),
                       "orthofinder_support": row["Support"], "type": row.get("Type", ""),
                       "source_file": "Gene_Duplication_Events/Duplications.tsv", "source_row": line,
                       "genes_1": left, "genes_2": right})
    sources.append({"path": "Gene_Duplication_Events/Duplications.tsv", "sha256": digest(duplication_file),
                    "selected_group_rows": len(events)})
    for edge in edges.values():
        edge["source_rows"] = ",".join(sorted(edge["source_rows"], key=int))
        edge["relationship"] = ",".join(sorted(edge["relationship"]))
    return [edges[key] for key in sorted(edges)], events, sources


def load_tree(folder, filename, expected_ids):
    try:
        from Bio import Phylo
    except ImportError as error:
        raise ValueError("Biopython is required: conda install -c conda-forge biopython") from error
    tree_path = folder / filename
    tree = Phylo.read(tree_path, "newick")
    tips = [tip.name for tip in tree.get_terminals()]
    if any(not tip for tip in tips) or len(tips) != len(set(tips)):
        raise ValueError(f"Missing or duplicate tree-tip identifiers: {tree_path}")
    ids = set(tips)
    if ids != expected_ids:
        raise ValueError(f"Tree/config membership mismatch: {tree_path}; "
                         f"missing={sorted(expected_ids - ids)}, unexpected={sorted(ids - expected_ids)}")
    if len(ids) < 4:
        raise ValueError(f"Too few tree tips: {tree_path}")
    return tree, ids


def label(node):
    return node.name or (f"{node.confidence:.15g}" if node.confidence is not None else "")


BRANCH_FIELDS = ["branch_id", "support_label", "tip_count", "species", "tip_ids",
                 "direct_orthology_edges_within", "exact_duplication_rows", "interpretation",
                 "family", "split_id", "rooting_mode", "complement_tip_count", "complement_tip_ids",
                 "source_orthogroups", "mixed_control_tip_ids", "orthology_edges_across_split",
                 "orthology_edges_in_complement", "duplication_rows_within",
                 "duplication_rows_within_count", "duplication_split_compatible_rows",
                 "duplication_split_compatible_count", "incomplete_duplication_rows"]


def branch_rows(family, tree, ids, panel_by_id, group_by_id, edges, events, rooted):
    """Unrooted results describe bipartitions, never ancestral duplication timing."""
    clades, seen = [], {}
    for node in tree.get_nonterminals(order="preorder"):
        if node is tree.root:
            continue
        displayed = {tip.name for tip in node.get_terminals()}
        other = ids - displayed
        if not displayed or not other:
            continue
        if not rooted and min(len(displayed), len(other)) < 2:
            continue
        smaller = min((tuple(sorted(displayed)), tuple(sorted(other))), key=lambda x: (len(x), x))
        split = hashlib.sha256(json.dumps(smaller).encode()).hexdigest()[:16]
        side = displayed if rooted else set(smaller)
        if not rooted and split in seen:
            previous = seen[split]
            if label(node) and label(previous) and label(node) != label(previous):
                raise ValueError(f"Conflicting support labels on a duplicate unrooted split: {family}")
            if label(node) and not label(previous):
                seen[split] = node
                for index, item in enumerate(clades):
                    if item[0] == split:
                        clades[index] = (split, node, side)
            continue
        seen[split] = node
        clades.append((split, node, side))
    output = []
    for index, (split, node, side) in enumerate(clades, 1):
        complement = ids - side
        within = sum({e["sequence_id_1"], e["sequence_id_2"]} <= side for e in edges)
        across = sum(bool({e["sequence_id_1"], e["sequence_id_2"]} & side) and
                     bool({e["sequence_id_1"], e["sequence_id_2"]} & complement) for e in edges)
        opposite = sum({e["sequence_id_1"], e["sequence_id_2"]} <= complement for e in edges)
        contained, compatible, incomplete, exact = [], [], [], []
        children = [{tip.name for tip in child.get_terminals()} for child in node.clades]
        for event in events:
            a, b = event["genes_1"], event["genes_2"]
            present = (a | b) & ids
            if not (a | b) <= ids:
                if present & side:
                    incomplete.append(event["event_id"])
                continue
            if (a | b) <= side:
                contained.append(event["event_id"])
            if (a <= side and b <= complement) or (b <= side and a <= complement):
                compatible.append(event["event_id"])
            if rooted and len(children) == 2 and {frozenset(a), frozenset(b)} == set(map(frozenset, children)):
                exact.append(event["event_id"])
        interpretation = ("User-declared rooted clade; duplication matches compare complete original child sets"
                          if rooted else "Unrooted split; compatible duplication partitions are topological context, "
                          "not a rooted duplication assignment or timing inference")
        output.append(dict(zip(BRANCH_FIELDS, [
            f"B{index:04d}", label(node), len(side),
            ",".join(sorted({panel_by_id[p] for p in side if panel_by_id[p] != "validation"})),
            ",".join(sorted(side)), within, ",".join(exact), interpretation,
            family, split, "user_declared_rooted" if rooted else "unrooted",
            len(complement), ",".join(sorted(complement)),
            ",".join(sorted({group_by_id[p] for p in side})),
            ",".join(sorted(p for p in side if panel_by_id[p] == "validation")), across, opposite,
            ",".join(contained), len(contained), ",".join(compatible), len(compatible), ",".join(incomplete)])))
    return output


def prepare_family(args, family, declared, manifest, by_group, panel_by_id):
    folder = args.run_dir / family
    if not (folder / args.tree_name).is_file():
        raise FileNotFoundError(f"Missing family tree: {folder / args.tree_name}")
    config_path = folder / "run_config.json"
    config = json.loads(config_path.read_text())
    if config["family"] != family or len(config["selected_ids"]) != len(set(config["selected_ids"])):
        raise ValueError(f"Invalid family/selected identifiers in {config_path}")
    fasta_key = f"{family}/input.fasta"
    if (config["library_input_sha256"] != manifest["files"][fasta_key] or
            digest(args.library / fasta_key) != manifest["files"][fasta_key]):
        raise ValueError(f"Run input does not match the prepared family library: {family}")
    tips_path = args.library / family / "tips.tsv"
    if digest(tips_path) != manifest["files"][f"{family}/tips.tsv"]:
        raise ValueError(f"Prepared tip metadata changed: {family}")
    groups = set(declared[family]["orthogroups"].split(","))
    group_by_id = {protein_id: group for group in groups for protein_id in by_group[group]}
    tips = {row["sequence_id"]: row for _, row in rows(tips_path, ["sequence_id", "source_orthogroup", "panel_label"])}
    if set(tips) != set(group_by_id) or any(tips[p]["source_orthogroup"] != group_by_id[p] or
                                           tips[p]["panel_label"] != panel_by_id[p] for p in tips):
        raise ValueError(f"Prepared tips disagree with original orthogroup membership: {family}")
    expected = {p for p in group_by_id if panel_by_id[p] not in set(config["excluded_panels"])}
    if expected != set(config["selected_ids"]):
        raise ValueError(f"Run selection does not match declared panel exclusions: {family}")
    tree, ids = load_tree(folder, args.tree_name, expected)
    # If the input tree is an original runner output, validate its recorded checksum.
    checkpoint = folder / "tree.complete.json"
    if checkpoint.exists():
        recorded = json.loads(checkpoint.read_text()).get("output_sha256", {}).get(args.tree_name)
        if recorded and recorded != digest(folder / args.tree_name):
            raise ValueError(f"Completed tree changed: {folder / args.tree_name}")
    if args.rooted and len(tree.root.clades) != 2:
        raise ValueError(f"--rooted requires an already rooted bifurcating input: {family}")
    return {"family": family, "folder": folder, "groups": groups, "tree": tree, "ids": ids,
            "group_by_id": group_by_id, "config": config,
            "input_receipts": {"tree_sha256": digest(folder / args.tree_name),
                               "run_config_sha256": digest(config_path), "tips_sha256": digest(tips_path)}}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True, type=Path, help="Existing results/<run-name> directory")
    parser.add_argument("--orthofinder-results", required=True, type=Path, help="Matching original Results_* directory")
    parser.add_argument("--out", required=True, type=Path, help="New separate annotation directory; must not exist")
    parser.add_argument("--library", type=Path, default=ROOT / "library")
    parser.add_argument("--family", action="append", help="Repeat to limit families; default: all prepared library families")
    parser.add_argument("--tree-name", default="tree.treefile", help="Existing Newick basename, e.g. gene_tree.treefile")
    parser.add_argument("--skip-incomplete", action="store_true", help="Explicitly skip missing trees, recording them in summary")
    parser.add_argument("--rooted", action="store_true", help="Declare that ALL supplied trees are already biologically rooted; does not reroot")
    args = parser.parse_args(argv)
    try:
        args.run_dir, args.library = args.run_dir.resolve(), args.library.resolve()
        results, out = args.orthofinder_results.resolve(), args.out.resolve()
        if not args.run_dir.is_dir():
            raise ValueError(f"Run directory does not exist: {args.run_dir}")
        if out.exists():
            raise ValueError(f"Output already exists; select a new --out directory: {out}")
        if Path(args.tree_name).name != args.tree_name:
            raise ValueError("--tree-name must be a filename, not a path")
        manifest_path = args.library / "library_manifest.json"
        manifest = json.loads(manifest_path.read_text())
        if digest(results / "Orthogroups/Orthogroups.tsv") != manifest["orthogroups_table_sha256"]:
            raise ValueError("Source orthogroup table differs from the prepared library; wrong Results_* directory")
        declared = {row["family"]: row for row in manifest["families"]}
        requested = list(dict.fromkeys(args.family or declared))
        if set(requested) - set(declared):
            raise ValueError(f"Unknown families: {set(requested) - set(declared)}")
        selected_groups = {g for name in requested for g in declared[name]["orthogroups"].split(",")}
        by_group, panel_by_id = source_membership(results, selected_groups)
        ready, skipped = [], []
        for name in requested:
            try:
                ready.append(prepare_family(args, name, declared, manifest, by_group, panel_by_id))
            except FileNotFoundError as error:
                if not args.skip_incomplete or (args.run_dir / name / args.tree_name).is_file():
                    raise
                skipped.append({"family": name, "status": "skipped_missing_tree", "reason": str(error)})
                print(f"SKIP {name}: tree not available", file=sys.stderr)
        if not ready:
            raise ValueError("No family trees available to annotate")
        selected_groups = set.union(*(item["groups"] for item in ready))
        subset = {group: by_group[group] for group in selected_groups}
        mapping = Identifiers(subset, panel_by_id, results)
        needed_panels = {panel_by_id[p] for item in ready for p in item["ids"]}
        # Include excluded-panel comparisons in the evidence sidecar too.
        needed_panels.update(panel_by_id[p] for group in subset.values() for p in group)
        print("Reading full native orthology/duplication evidence (no inference is launched).", flush=True)
        edges, events, sources = evidence(results, mapping, selected_groups, needed_panels)
        source_manifest = results.parent.parent / "workflow_manifest.json"
        if not source_manifest.is_file() or digest(source_manifest) != manifest["source_workflow_manifest_sha256"]:
            raise ValueError("Source orthology manifest does not match the prepared library")
        recorded = json.loads(source_manifest.read_text())["files"]
        if mapping.mapping_file:
            sources.append({"path": mapping.mapping_file.relative_to(results).as_posix(),
                            "sha256": digest(mapping.mapping_file), "purpose": "explicit_identifier_aliases"})
        for receipt in sources:
            key = (results.relative_to(source_manifest.parent) / receipt["path"]).as_posix()
            if recorded.get(key) != receipt["sha256"]:
                raise ValueError(f"Source evidence differs from saved orthology manifest: {receipt['path']}")
        receipts = {"library_manifest_sha256": digest(manifest_path),
                    "orthofinder_results": str(results), "orthology_manifest_sha256": digest(source_manifest),
                    "source_evidence": sources,
                    "sequence_ids_sha256": digest(mapping.mapping_file) if mapping.mapping_file else None,
                    "identifier_policy": "Exact original ID, exact species-prefixed ID, or explicit SequenceIDs alias"}
        prepared = []
        for item in ready:
            family_edges = [e for e in edges if e["orthogroup"] in item["groups"]]
            family_events = [e for e in events if e["orthogroup"] in item["groups"]]
            active_edges = [e for e in family_edges if {e["sequence_id_1"], e["sequence_id_2"]} <= item["ids"]]
            branches = branch_rows(item["family"], item["tree"], item["ids"], panel_by_id,
                                   item["group_by_id"], active_edges, family_events, args.rooted)
            prepared.append((item, family_edges, family_events, branches))
        # All source/identifier checks and annotation computation finish before publishing output.
        out.mkdir(parents=True)
        summary = []
        for item, family_edges, family_events, branches in prepared:
            family = item["family"]
            folder = out / family
            folder.mkdir()
            write_tsv(folder / "branch_evidence.tsv", branches, BRANCH_FIELDS)
            edge_fields = ["orthogroup", "sequence_id_1", "sequence_id_2", "relationship", "source_file", "source_rows", "both_tips_in_tree"]
            edge_rows = [dict(e, both_tips_in_tree=int({e["sequence_id_1"], e["sequence_id_2"]} <= item["ids"])) for e in family_edges]
            write_tsv(folder / "orthology_evidence.tsv", edge_rows, edge_fields)
            event_fields = ["event_id", "orthogroup", "species_tree_node", "gene_tree_node", "orthofinder_support", "type",
                            "source_file", "source_row", "genes_1", "genes_2", "complete_in_tree", "missing_tip_ids"]
            event_rows = []
            for event in family_events:
                missing = (event["genes_1"] | event["genes_2"]) - item["ids"]
                event_rows.append(dict(event, genes_1=",".join(sorted(event["genes_1"])),
                                       genes_2=",".join(sorted(event["genes_2"])), complete_in_tree=int(not missing),
                                       missing_tip_ids=",".join(sorted(missing))))
            write_tsv(folder / "duplication_evidence.tsv", event_rows, event_fields)
            provenance = {"created_utc": datetime.now(timezone.utc).isoformat(), "family": family,
                          "input_tree": str(item["folder"] / args.tree_name), "inputs": item["input_receipts"],
                          "selected_tip_ids": sorted(item["ids"]), "excluded_panels": item["config"]["excluded_panels"],
                          "source_orthogroups": sorted(item["groups"]), "rooting": "user_declared_rooted" if args.rooted else "unrooted",
                          "orthology_scope": "full native pairwise tables; unique unordered protein pairs; no edges inferred across orthogroups",
                          "duplication_policy": "Only complete original events contribute within/split/exact matches; incomplete events remain explicit",
                          "evidence_provenance_file": "../evidence_sources.json",
                          "output_sha256": {p.name: digest(p) for p in sorted(folder.glob("*.tsv"))}}
            write_json(folder / "branch_evidence.provenance.json", provenance)
            summary.append({"family": family, "status": "complete", "tree_tips": len(item["ids"]),
                            "branches": len(branches), "source_orthology_pairs": len(family_edges),
                            "pairs_in_tree": sum(e["both_tips_in_tree"] for e in edge_rows),
                            "source_duplication_events": len(family_events),
                            "complete_duplication_events": sum(e["complete_in_tree"] for e in event_rows), "reason": ""})
            print(f"WROTE {folder / 'branch_evidence.tsv'}", flush=True)
        summary.extend(skipped)
        write_tsv(out / "annotation_summary.tsv", summary,
                  ["family", "status", "tree_tips", "branches", "source_orthology_pairs", "pairs_in_tree",
                   "source_duplication_events", "complete_duplication_events", "reason"])
        write_json(out / "evidence_sources.json", receipts)
        return 0
    except (OSError, ValueError, KeyError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
