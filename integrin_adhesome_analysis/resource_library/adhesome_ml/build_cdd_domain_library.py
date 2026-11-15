#!/usr/bin/env python3
"""Extract traceable, deduplicated CDD domain sequences from labelled controls."""
import argparse
import csv
import hashlib
import json
import re
from collections import defaultdict
from pathlib import Path

HIT_CLASS_MAP = {"specific": "specific", "non-specific": "non-specific",
                 "superfamily": "broad", "broad": "broad"}
HIT_CLASSES = set(HIT_CLASS_MAP.values())
FIELDS = ["library_id", "sequence_sha256", "source_protein_id", "origin", "cdd_accession",
          "pssm_id", "domain_annotation", "hit_classification", "superfamily_accession",
          "start_1based", "end_1based", "evalue", "bitscore", "domain_length_aa"]


def read_fasta(path):
    result = {}
    name, chunks = None, []
    def save():
        if name is not None:
            seq = "".join(chunks).upper().replace("*", "X").replace("-", "").replace(".", "")
            if not seq:
                raise ValueError(f"Empty protein sequence: {name} ({path})")
            if name in result and result[name] != seq:
                raise ValueError(f"Conflicting duplicate FASTA identifier {name}: {path}")
            result[name] = seq
    with path.open(encoding="utf-8-sig") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                save()
                name, chunks = line[1:].split()[0], []
            elif name is None:
                raise ValueError(f"Sequence precedes FASTA header: {path}")
            else:
                chunks.append(line)
    save()
    return result


def tsv(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def cdd_rows(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        lines = (line for line in handle if line.strip() and not line.startswith("#"))
        reader = csv.DictReader(lines, delimiter="\t")
        if not reader.fieldnames or "Query" not in reader.fieldnames:
            raise ValueError(f"No NCBI CDD hit table header in {path}")
        yield from reader


def query_id(value):
    value = re.sub(r"^Q#\d+\s*-\s*>", "", str(value or "").strip())
    return value.split()[0] if value else ""


def prepared_sequences(branch):
    found = {}
    for path in (branch / "handoff/detection").glob("*/prepared/*.faa"):
        for ident, seq in read_fasta(path).items():
            if ident in found and found[ident] != seq:
                raise ValueError(f"Protein ID appears with conflicting sequences: {ident}")
            found[ident] = seq
    if not found:
        raise ValueError(f"No prepared FASTA files in {branch}/handoff/detection/*/prepared")
    return found


def build(args):
    if args.out.exists():
        raise ValueError(f"Output exists; select a new output directory: {args.out}")
    records, dedup = [], set()
    stats = defaultdict(int)
    panels = (("calibration", args.calibration_run, args.calibration_cdd),
              ("validation", args.validation_run, args.validation_cdd))
    for origin, branch, cdd in panels:
        labels = {row["protein_id"]: row for row in tsv(branch / "reports/predictions.tsv")
                  if str(row.get("control_class", "")).strip().casefold() in {"positive", "positive_control", "adhesome_associated"}}
        proteins = prepared_sequences(branch)
        missing = set(labels) - proteins.keys()
        if missing:
            raise ValueError(f"{origin}: {len(missing)} positive controls lack prepared sequences; examples: {sorted(missing)[:5]}")
        stats[origin + "_positive_controls"] = len(labels)
        for hit in cdd_rows(cdd):
            ident = query_id(hit.get("Query"))
            if ident not in labels:
                continue
            source_classification = str(hit.get("Hit type", "")).strip().casefold()
            classification = HIT_CLASS_MAP.get(source_classification)
            if not classification:
                continue
            try:
                start, end = int(hit["From"]), int(hit["To"])
            except (KeyError, TypeError, ValueError):
                stats[origin + "_invalid_coordinates"] += 1
                continue
            protein = proteins[ident]
            if start < 1 or end < start or end > len(protein):
                raise ValueError(f"Out-of-range CDD coordinates for {ident}: {start}-{end}, sequence length {len(protein)}")
            domain = protein[start - 1:end]
            digest = hashlib.sha256(domain.encode()).hexdigest()
            accession = str(hit.get("Accession", "")).strip()
            annotation = str(hit.get("Short name", "")).strip()
            key = (origin, ident, accession, classification, start, end, digest)
            if key in dedup:
                stats[origin + "_duplicate_hits_skipped"] += 1
                continue
            dedup.add(key)
            records.append({"library_id": "CDDSEQ_" + digest[:20], "sequence_sha256": digest,
                "source_protein_id": ident, "origin": origin, "cdd_accession": accession,
                "pssm_id": str(hit.get("PSSM-ID", "")).strip(), "domain_annotation": annotation,
                "hit_classification": classification,
                "superfamily_accession": str(hit.get("Superfamily", "")).strip(),
                "start_1based": start, "end_1based": end,
                "evalue": str(hit.get("E-Value", "")).strip(),
                "bitscore": str(hit.get("Bitscore", "")).strip(), "domain_length_aa": len(domain),
                "_sequence": domain})
            stats[origin + "_" + classification + "_hit_records"] += 1

    if not records:
        raise ValueError("No specific, non-specific or broad CDD hits mapped to positive controls")
    args.out.mkdir(parents=True)
    sequences = {}
    for row in records:
        sequences.setdefault(row["library_id"], row["_sequence"])
    with (args.out / "cdd_domain_sequences.faa").open("w", encoding="utf-8", newline="\n") as handle:
        for ident, sequence in sorted(sequences.items()):
            handle.write(f">{ident}\n")
            for i in range(0, len(sequence), 60):
                handle.write(sequence[i:i+60] + "\n")
    with (args.out / "cdd_domain_metadata.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({key: value for key, value in row.items() if key != "_sequence"} for row in records)
    manifest = {"format_version": "1.0", "source": "NCBI CD-Search returned hit tables joined to labelled positive controls",
        "coordinate_convention": "1-based inclusive input coordinates; FASTA extraction uses [From-1:To]",
        "unique_domain_sequences": len(sequences), "traceable_hit_records": len(records),
        "classification_levels": {key: sum(r["hit_classification"] == key for r in records) for key in sorted(HIT_CLASSES)},
        "sequence_deduplication": "Identical extracted amino-acid strings share one FASTA record; all source hits remain as metadata rows",
        "counts": dict(stats)}
    (args.out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--calibration-run", type=Path, required=True)
    parser.add_argument("--validation-run", type=Path, required=True)
    parser.add_argument("--calibration-cdd", type=Path, required=True)
    parser.add_argument("--validation-cdd", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    try:
        build(args)
    except (OSError, ValueError, KeyError, csv.Error) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()
