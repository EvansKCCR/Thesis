#!/usr/bin/env python3
"""Read-only integrity and provenance checks; does not invoke analytical tools."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

from prepare_library import read_fasta, sha256

ROOT = Path(__file__).resolve().parent


def check(library, source_results=None):
    manifest = json.loads((library / "library_manifest.json").read_text())
    reconciled = manifest.get("library_kind") == "reconciled_candidates"
    if reconciled and source_results:
        raise ValueError("Reconciled inputs use candidate-level provenance; omit --source-results")
    if not reconciled and sha256(ROOT / "families.json") != manifest["families_config_sha256"]:
        raise ValueError("Family selection configuration changed")
    for relative, expected in manifest["files"].items():
        if sha256(library / relative) != expected:
            raise ValueError(f"Prepared file changed: {relative}")
    total = 0
    for family in manifest["families"]:
        folder = library / family["family"]
        records = read_fasta(folder / "input.fasta")
        with (folder / "tips.tsv").open(newline="") as handle:
            tips = list(csv.DictReader(handle, delimiter="\t"))
        if len(tips) != len(records) or {row["sequence_id"] for row in tips} != set(records):
            raise ValueError(f"Tip mapping changed: {family['family']}")
        if len(records) != family["records"]:
            raise ValueError(f"Record count changed: {family['family']}")
        for tip in tips:
            header, seq = records[tip["sequence_id"]]
            if ((not reconciled and header != tip["original_header"]) or len(seq) != int(tip["length_aa"]) or
                    hashlib.sha256(seq.encode("ascii")).hexdigest() != tip["sequence_sha256"]):
                raise ValueError(f"Tip sequence/header changed: {tip['sequence_id']}")
        if source_results:
            originals = {}
            for og in family["orthogroups"].split(","):
                originals.update(read_fasta(source_results / "Orthogroup_Sequences" / f"{og}.fa"))
            if originals != records:
                raise ValueError(f"Prepared family differs from exact source sequences: {family['family']}")
        total += len(records)
    if source_results:
        if sha256(source_results / "Orthogroups/Orthogroups.tsv") != manifest["orthogroups_table_sha256"]:
            raise ValueError("Source membership table changed")
        for row in manifest["source_orthogroups"]:
            if sha256(source_results / "Orthogroup_Sequences" / (row["orthogroup"] + ".fa")) != row["sha256"]:
                raise ValueError(f"Source orthogroup changed: {row['orthogroup']}")
    print(f"PASS: {len(manifest['families'])} families; {len(manifest.get('source_orthogroups', []))} source orthogroups; "
          f"{total} records; file hashes, identifiers, sequences and tip mappings verified.")
    if source_results:
        print("PASS: complete source FASTA membership and original headers/residues preserved exactly.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--library", type=Path, default=ROOT / "library")
    parser.add_argument("--source-results", type=Path)
    args = parser.parse_args()
    check(args.library, args.source_results)
