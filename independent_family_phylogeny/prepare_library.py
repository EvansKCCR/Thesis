#!/usr/bin/env python3
"""Export exact, full orthogroup membership; never launch analytical tools."""
import argparse
import collections
import csv
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_fasta(path):
    records = {}
    header = None
    chunks = []
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if header is not None:
                add_record(records, header, "".join(chunks), path)
            header, chunks = line[1:].strip(), []
        elif line.strip():
            if header is None:
                raise ValueError(f"Sequence before first header: {path}")
            chunks.append(line.strip())
    if header is not None:
        add_record(records, header, "".join(chunks), path)
    if not records:
        raise ValueError(f"Empty FASTA: {path}")
    return records


def add_record(records, header, sequence, path):
    protein_id = header.split()[0]
    if protein_id in records or not sequence:
        raise ValueError(f"Duplicate ID or empty sequence: {path}: {protein_id}")
    if set(sequence) - set("ACDEFGHIKLMNPQRSTVWYXBZJUO"):
        raise ValueError(f"Unexpected amino-acid characters: {path}: {protein_id}")
    records[protein_id] = (header, sequence)


def write_tsv(path, rows, fields):
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def prepare(results, out, append=False):
    results, out = results.resolve(), out.resolve()
    if out.exists() and not append:
        raise ValueError(f"Refusing to overwrite prepared library: {out}")
    config = json.loads((ROOT / "families.json").read_text())
    previous = None
    protected = {}
    previous_manifest_sha256 = None
    if append:
        manifest_path = out / "library_manifest.json"
        previous = json.loads(manifest_path.read_text())
        previous_manifest_sha256 = sha256(manifest_path)
        for relative, expected in previous["files"].items():
            if sha256(out / relative) != expected:
                raise ValueError(f"Existing library changed; append refused: {relative}")
            if "/" in relative:
                protected[relative] = expected
        for entry in previous["families"]:
            family = entry["family"]
            if family not in config or config[family]["orthogroups"] != entry["orthogroups"].split(",") or any(
                    config[family][key] != entry[key] for key in ("label", "category")):
                raise ValueError(f"Append cannot change or remove an existing family: {family}")
        old_names = {entry["family"] for entry in previous["families"]}
        additions = [family for family in config if family not in old_names]
        if not additions:
            if sha256(ROOT / "families.json") != previous["families_config_sha256"]:
                raise ValueError("Family configuration changed without an additive extension")
            print("Library already contains all configured families; nothing changed.")
            return
        for family in additions:
            if (out / family).exists():
                raise ValueError(f"New-family destination already exists; append refused: {out / family}")
    selected = [og for entry in config.values() for og in entry["orthogroups"]]
    if len(selected) != len(set(selected)):
        raise ValueError("An orthogroup is assigned more than once")
    table = results / "Orthogroups/Orthogroups.tsv"
    upstream_manifest_path = results.parent.parent / "workflow_manifest.json"
    upstream_manifest = json.loads(upstream_manifest_path.read_text()) if upstream_manifest_path.exists() else None
    if previous and (sha256(table) != previous["orthogroups_table_sha256"] or
                     not upstream_manifest or sha256(upstream_manifest_path) != previous["source_workflow_manifest_sha256"]):
        raise ValueError("Append requires the exact original source orthogroup table and orthology manifest")
    with table.open(encoding="utf-8-sig", newline="") as handle:
        membership = {row["Orthogroup"]: row for row in csv.DictReader(handle, delimiter="\t")
                      if row["Orthogroup"] in selected}
    pending, source_rows = {}, []
    all_ids = set()
    for family, entry in config.items():
        records, tips = {}, []
        for og in entry["orthogroups"]:
            path = results / "Orthogroup_Sequences" / f"{og}.fa"
            if upstream_manifest is not None:
                key = str(path.relative_to(upstream_manifest_path.parent)).replace("\\", "/")
                if upstream_manifest["files"].get(key) != sha256(path):
                    raise ValueError(f"Source FASTA does not match saved orthology manifest: {path}")
            members = read_fasta(path)
            expected = {}
            for panel, cell in membership[og].items():
                if panel == "Orthogroup":
                    continue
                for protein_id in filter(None, (x.strip() for x in cell.split(","))):
                    if protein_id in expected:
                        raise ValueError(f"Ambiguous species mapping: {og}: {protein_id}")
                    expected[protein_id] = panel
            if set(members) != set(expected):
                raise ValueError(f"FASTA/table membership mismatch: {og}")
            for protein_id, (header, sequence) in members.items():
                if protein_id in all_ids:
                    raise ValueError(f"Repeated orthogroup member: {protein_id}")
                all_ids.add(protein_id)
                records[protein_id] = (header, sequence)
                panel = expected[protein_id]
                tips.append({"sequence_id": protein_id, "source_orthogroup": og,
                             "panel_label": panel,
                             "record_role": "mixed_control_panel" if panel == "validation" else "proteome",
                             "length_aa": len(sequence),
                             "sequence_sha256": hashlib.sha256(sequence.encode("ascii")).hexdigest(),
                             "original_header": header})
            source_rows.append({"family": family, "orthogroup": og, "records": len(members),
                                "source_fasta": str(path), "sha256": sha256(path)})
        pending[family] = (entry, records, tips)
    # Finish all source checks before writing any library output.
    out.mkdir(parents=True, exist_ok=append)
    summary = []
    for family, (entry, records, tips) in pending.items():
        folder = out / family
        existing = previous is not None and family in old_names
        if not existing:
            folder.mkdir()
            with (folder / "input.fasta").open("w", encoding="ascii", newline="\n") as handle:
                for header, sequence in records.values():
                    handle.write(f">{header}\n")
                    for start in range(0, len(sequence), 80):
                        handle.write(sequence[start:start + 80] + "\n")
            write_tsv(folder / "tips.tsv", tips, list(tips[0]))
        counts = collections.Counter(tip["panel_label"] for tip in tips)
        if not existing:
            write_tsv(folder / "panel_counts.tsv",
                      [{"panel_label": panel, "records": count} for panel, count in sorted(counts.items())],
                      ["panel_label", "records"])
        lengths = [len(seq) for _, seq in records.values()]
        summary.append({"family": family, "label": entry["label"], "category": entry["category"],
                        "orthogroups": ",".join(entry["orthogroups"]), "records": len(records),
                        "unique_sequences": len({seq for _, seq in records.values()}),
                        "proteome_panels": sum(panel != "validation" for panel in counts),
                        "mixed_control_records": counts["validation"],
                        "min_length_aa": min(lengths), "max_length_aa": max(lengths)})
    write_tsv(out / "family_summary.tsv", summary, list(summary[0]))
    write_tsv(out / "source_orthogroups.tsv", source_rows, list(source_rows[0]))
    files = {str(path.relative_to(out)).replace("\\", "/"): sha256(path)
             for path in sorted(out.rglob("*")) if path.is_file() and path.name != "library_manifest.json"}
    manifest = {"bundle_version": "1.0", "created_utc": datetime.now(timezone.utc).isoformat(),
                "status": "prepared_only_no_analysis", "source_results": str(results),
                "orthogroups_table_sha256": sha256(table),
                "source_workflow_manifest": str(upstream_manifest_path) if upstream_manifest else None,
                "source_workflow_manifest_sha256": sha256(upstream_manifest_path) if upstream_manifest else None,
                "families_config_sha256": sha256(ROOT / "families.json"),
                "selection": f"Complete membership of the {len(selected)} explicitly requested orthogroups",
                "sequence_policy": "No subsampling, sequence deduplication, domain extraction or renaming",
                "mixed_control_panel": "validation is a mixed reference collection, not a species",
                "families": summary, "source_orthogroups": source_rows, "files": files}
    if previous:
        for relative, expected in protected.items():
            if sha256(out / relative) != expected:
                raise ValueError(f"Protected family file changed during append: {relative}")
        manifest["extensions"] = list(previous.get("extensions", [])) + [{
            "created_utc": datetime.now(timezone.utc).isoformat(), "added_families": additions,
            "previous_manifest_sha256": previous_manifest_sha256,
            "preserved_family_file_sha256": protected,
            "interpretation": "Additive library preparation only; no run files modified or tools launched"}]
    publish_path = out / "library_manifest.publish.tmp"
    with publish_path.open("x", encoding="utf-8") as handle:
        handle.write(json.dumps(manifest, indent=2) + "\n")
    os.replace(publish_path, out / "library_manifest.json")
    print(f"Prepared {len(summary)} families, {len(selected)} orthogroups, {len(all_ids)} records: {out}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True, help="Existing OrthoFinder Results_* directory")
    parser.add_argument("--out", type=Path, default=ROOT / "library")
    parser.add_argument("--append", action="store_true", help="Add only new configured families; preserve existing family files byte for byte")
    args = parser.parse_args()
    prepare(args.results, args.out, append=args.append)


if __name__ == "__main__":
    main()
