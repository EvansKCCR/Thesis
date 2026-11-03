#!/usr/bin/env python3
"""Independent, standard-library handoff: completed RBH + audit -> phylogeny inputs.

No BLAST, audit, alignment or tree inference is invoked. Default is a read-only
plan. Use --prepare to write a NEW library, or --check to validate that library.
"""
import argparse
import csv
import gzip
import hashlib
import json
import os
import re
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent
BUNDLE = ROOT.parent
ACCEPTED_RBH = {"RBH_candidate", "RBH_candidate_augmented_reference", "provisional_reference_only_RBH"}
AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def sequence_hash(seq):
    return hashlib.sha256(seq.encode("ascii")).hexdigest()


def load_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


def read_table(path, delimiter="\t"):
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def fasta(path, allow_empty=False):
    result, name, chunks = {}, None, []
    def save():
        if name is not None:
            if name in result:
                raise ValueError(f"Duplicate FASTA ID {name}: {path}")
            seq = "".join(chunks).upper().rstrip("*")
            if not seq or not re.fullmatch(r"[A-Z*]+", seq):
                raise ValueError(f"Empty/non-protein FASTA sequence {name}: {path}")
            result[name] = seq
    opener = gzip.open if str(path).lower().endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8-sig") as handle:
        for line in handle:
            if line.startswith(">"):
                save()
                if not line[1:].strip():
                    raise ValueError(f"Empty FASTA header: {path}")
                name, chunks = line[1:].split()[0], []
            elif line.strip():
                if name is None:
                    raise ValueError(f"Sequence before FASTA header: {path}")
                chunks.append("".join(line.split()))
    save()
    if not result and not allow_empty:
        raise ValueError(f"Empty FASTA: {path}")
    return result


def native_path(value):
    """Translate recorded Windows/WSL paths; never guess biological IDs."""
    value = str(value).replace("\\", "/")
    if os.name == "nt":
        match = re.match(r"^/mnt/([a-zA-Z])/(.*)$", value)
        if match:
            value = f"{match[1].upper()}:/{match[2]}"
    elif re.match(r"^[a-zA-Z]:/", value):
        value = f"/mnt/{value[0].lower()}/{value[3:]}"
    return Path(value)


def safe_relative(value):
    path = Path(value)
    if path.is_absolute() or ".." in path.parts or "\\" in value or ":" in value:
        raise ValueError(f"Unsafe library-relative path: {value}")
    return path


class Sources:
    def __init__(self):
        self.records = []

    def add(self, path, role, expected=None, recorded_path=None):
        path = Path(path).resolve()
        actual = digest(path)
        if expected is not None and actual != expected:
            raise ValueError(f"Saved source hash changed ({role}): {path}")
        self.records.append({"role": role, "path": str(path), "sha256": actual,
                             "recorded_path": recorded_path, "saved_hash_verified": expected is not None})
        return path


def verified_output(folder, manifest, name, sources, role):
    expected = manifest.get("files", manifest.get("output_sha256", {})).get(name)
    if not expected:
        raise ValueError(f"Saved manifest has no hash for {name}: {folder}")
    return sources.add(folder / safe_relative(name), role, expected)


def input_path(recorded, inputs, anchors, hashes, sources, role):
    """Prefer the relocated bundle copy, checking the original run's exact hash."""
    normalized = recorded.replace("\\", "/")
    expected = hashes.get(recorded)
    if not expected:
        raise ValueError(f"No saved input hash for {recorded}")
    candidates = []
    if role.startswith('target:') and '/input/target_raw/' in normalized:
        sys.path.insert(0,str(BUNDLE/'integrin_adhesome_analysis/resource_library/scripts'))
        from input_dependencies import raw_input
        candidates.append(raw_input('target_raw')/safe_relative(normalized.rsplit('/input/target_raw/',1)[1]))
    for anchor in anchors:
        marker = f"/{anchor}/"
        if marker in normalized:
            suffix = normalized.rsplit(marker, 1)[1]
            for local in anchors:
                candidates.append(inputs / local / safe_relative(suffix))
    candidates.append(native_path(recorded))
    for candidate in candidates:
        if candidate.is_file():
            return sources.add(candidate, role, expected, recorded)
    raise ValueError(f"Missing {role}: {recorded}; supply the matching --homology-inputs bundle")


def crosswalk(path):
    config = load_json(path)
    by_set, by_label, families = {}, {}, {}
    for row in config["families"]:
        family = row["family"]
        slug = family.lower()
        if not re.fullmatch(r"[a-z0-9][a-z0-9_]*", slug) or slug in families:
            raise ValueError(f"Invalid/duplicate family slug: {family}")
        families[slug] = row
        for reference_set in row["reference_sets"]:
            if reference_set in by_set:
                raise ValueError(f"Reference set mapped twice: {reference_set}")
            by_set[reference_set] = slug
        for label in row["main_labels"]:
            if label.casefold() in by_label:
                raise ValueError(f"Audit label mapped twice: {label}")
            by_label[label.casefold()] = slug
    return config, families, by_set, by_label


def handoff_path(audit_manifest, explicit, sources):
    if explicit:
        return explicit.resolve()
    annotation = native_path(audit_manifest["annotation"])
    path = sources.add(annotation / "workflow_manifest.json", "annotation_provenance")
    handoff = load_json(path).get("handoff")
    if not handoff:
        raise ValueError("Annotation manifest has no handoff; supply --handoff explicitly")
    return native_path(handoff).resolve()


def make_record(species, raw_ids, seq, main=None, main_accepted=False):
    return {"species": species, "raw_ids": set(raw_ids), "seq": seq, "main": main or {},
            "main_accepted": main_accepted, "rbh": [], "rbh_families": set(), "relationships": set()}


def reconcile(args):
    sources = Sources()
    sources.add(Path(__file__), "reconciliation_code")
    sources.add(args.crosswalk, "family_crosswalk")
    config, families, by_set, by_label = crosswalk(args.crosswalk)
    hr, audit = args.homology_run.resolve(), args.audit.resolve()
    hm = load_json(sources.add(hr / "run_manifest.json", "homology_manifest"))
    hc = load_json(sources.add(hr / "run_config.json", "homology_configuration"))
    if hm.get("status") != "complete":
        raise ValueError("A completed homology run is required")
    am = load_json(sources.add(audit / "workflow_manifest.json", "audit_manifest"))
    if am.get("kind") != "audit":
        raise ValueError("--audit must point to a completed audit directory")
    main_rows = read_table(verified_output(audit, am, "audited_family_assignments.tsv", sources, "audit_calls"))
    exclusions = read_table(verified_output(audit, am, "negative_control_exclusions.tsv", sources, "negative_gates"))
    rbhs = read_table(verified_output(hr, hm, "homology_audit.tsv", sources, "rbh_evidence"))
    mappings = read_table(verified_output(hr, hm, "reference_mapping.tsv", sources, "reference_mapping"))
    global_rbh = fasta(verified_output(hr, hm, "rbh_candidates.faa", sources, "rbh_fasta"), True)
    report = hm["input_report"]
    sets = {row["reference_set"]: row for row in report["reference_sets"]}
    if set(sets) - set(by_set):
        raise ValueError(f"Add explicit family crosswalk entries for: {sorted(set(sets) - set(by_set))}")
    target_sources = {}
    for row in report["target_proteomes"]:
        species = row["species"]
        if species in target_sources:
            raise ValueError(f"Duplicate target species: {species}")
        path = input_path(row["file"], args.homology_inputs, ["target_proteomes", "schistosome_proteomes"],
                          hc["input_sha256"], sources, f"target:{species}")
        target_sources[species] = (path, fasta(path))
    species_scope = set(target_sources)
    main, negative = {}, {}
    for row in main_rows:
        key = (row["Specie"], row["protein_id"])
        if key[0] not in species_scope:
            continue
        if key in main:
            raise ValueError(f"Duplicate audit key: {key}")
        main[key] = row
    for row in exclusions:
        if row["Specie"] in species_scope:
            negative[(row["Specie"], row["protein_id"])] = row
    handoff = handoff_path(am, args.handoff, sources)
    id_maps, retained, prepared, preparation_checks = {}, {}, {}, []
    for species, (target_path, raw) in target_sources.items():
        folder = handoff / "detection" / species / "prepared"
        prep_path = sources.add(folder / "preparation.json", f"preparation:{species}")
        prep = load_json(prep_path)["species"][species]
        seq_path = sources.add(folder / f"{species}.faa", f"prepared_sequences:{species}", prep["output_sha256"])
        map_path = sources.add(folder / f"{species}.id_map.csv", f"id_map:{species}")
        prepared[species] = fasta(seq_path)
        id_maps[species], retained[species] = {}, {}
        for row in read_table(map_path, ","):
            if row["original_id"] in id_maps[species]:
                raise ValueError(f"Duplicate original ID in map: {species}/{row['original_id']}")
            id_maps[species][row["original_id"]] = row
            if row["decision"] == "retained":
                pid = row["protein_id"]
                if pid in retained[species] or pid not in prepared[species]:
                    raise ValueError(f"Invalid retained ID in map: {species}/{pid}")
                retained[species][pid] = row["original_id"]
        preparation_checks.append({"species": species, "same_raw_file_sha256": digest(target_path) == prep["input_sha256"]})
    accepted = set(config["main_accepted_decisions"])
    catalogue, excluded, conflicts = {}, [], []
    for (species, pid), row in main.items():
        if row["family_decision"] not in accepted:
            continue
        if (species, pid) in negative:
            excluded.append({"species": species, "protein_id": pid, "source": "main",
                             "reason": negative[(species, pid)]["exclusion_reason"]})
            continue
        if pid not in retained[species]:
            raise ValueError(f"Accepted audit protein lacks retained mapping: {species}/{pid}")
        raw_id = retained[species][pid]
        seq = prepared[species][pid]
        if target_sources[species][1].get(raw_id) != seq:
            raise ValueError(f"Accepted main sequence disagrees with current target: {species}/{raw_id}")
        catalogue[(species, "main", pid)] = make_record(species, [raw_id], seq, row, True)
    accepted_rbh = [r for r in rbhs if r["status"] in ACCEPTED_RBH]
    expected_global, expected_by_set = {}, defaultdict(dict)
    details = defaultdict(list)
    for row in accepted_rbh:
        species, raw_id, reference_set = row["species"], row["best_target_id"], row["reference_set"]
        if species not in species_scope or reference_set not in sets:
            raise ValueError("RBH audit refers to an unknown target/reference set")
        full_id = f"{species}|{raw_id}"
        seq = target_sources[species][1].get(raw_id)
        if seq is None:
            raise ValueError(f"RBH candidate missing from target proteome: {full_id}")
        expected_global[full_id] = seq
        expected_by_set[reference_set][full_id] = seq
        details[(species, raw_id)].append(row)
    if global_rbh != expected_global:
        raise ValueError("Global RBH FASTA does not match accepted audit hits and original target residues")
    for reference_set in sets:
        name = f"reference_sets/{reference_set}/rbh_candidates.faa"
        actual = fasta(verified_output(hr, hm, name, sources, f"rbh_set:{reference_set}"), True)
        if actual != expected_by_set[reference_set]:
            raise ValueError(f"Per-set RBH FASTA disagrees with audit evidence: {reference_set}")
    reconciliation_rows = []
    for (species, raw_id), rows in sorted(details.items()):
        seq = target_sources[species][1][raw_id]
        mapped = id_maps[species].get(raw_id, {})
        pid = mapped.get("protein_id", "")
        main_row = main.get((species, pid), {})
        match = "unmapped_original_id"
        key = (species, "rbh", raw_id)
        if pid in prepared[species]:
            if prepared[species][pid] == seq:
                key = (species, "main", pid)
                match = "exact_sequence_mapped_gene"
            else:
                match = "alternative_isoform_sequence" if mapped.get("decision") != "retained" else "sequence_release_conflict"
                conflicts.append({"species": species, "original_id": raw_id, "protein_id": pid,
                                  "relationship": match, "rbh_sequence_sha256": sequence_hash(seq),
                                  "prepared_sequence_sha256": sequence_hash(prepared[species][pid])})
        neg = negative.get((species, pid))
        rec = {"species": species, "original_id": raw_id, "main_protein_id": pid,
               "sequence_relationship": match, "id_map_decision": mapped.get("decision", ""),
               "main_audit_present": bool(main_row), "main_family": main_row.get("putative_audited_family", ""),
               "main_decision": main_row.get("family_decision", ""),
               "main_stage1_hypothesis": main_row.get("stage1_family_hypothesis", ""),
               "rbh_reference_sets": ";".join(sorted({r["reference_set"] for r in rows})),
               "rbh_statuses": ";".join(sorted({r["status"] for r in rows})),
               "negative_excluded": bool(neg), "negative_reason": neg["exclusion_reason"] if neg else ""}
        reconciliation_rows.append(rec)
        if neg:
            excluded.append({"species": species, "protein_id": raw_id, "source": "rbh", "reason": neg["exclusion_reason"]})
            continue
        if key not in catalogue:
            catalogue[key] = make_record(species, [raw_id], seq, main_row)
        candidate = catalogue[key]
        candidate["raw_ids"].add(raw_id)
        candidate["rbh"].extend(rows)
        candidate["rbh_families"].update(by_set[r["reference_set"]] for r in rows)
        candidate["relationships"].add(match)
    records, family_members = [], defaultdict(list)
    for n, (key, candidate) in enumerate(sorted(catalogue.items()), 1):
        row, seq = candidate["main"], candidate["seq"]
        main_family = by_label.get(row.get("putative_audited_family", "").casefold()) if candidate["main_accepted"] else None
        claims = set(candidate["rbh_families"])
        if main_family:
            claims.add(main_family)
        evidence = "both" if candidate["main_accepted"] and candidate["rbh"] else "main_only" if candidate["main_accepted"] else "rbh_only"
        conflict = bool(main_family and candidate["rbh_families"] - {main_family}) or bool(
            candidate["main_accepted"] and candidate["rbh"] and main_family is None)
        record = {"sequence_id": f"C{n:06d}", "panel_label": candidate["species"], "species": candidate["species"],
                  "original_ids": ";".join(sorted(candidate["raw_ids"])), "evidence_source": evidence,
                  "main_protein_id": row.get("protein_id", key[2] if key[1] == "main" else ""),
                  "main_audit_present": bool(row), "main_accepted": candidate["main_accepted"],
                  "main_family": row.get("putative_audited_family", ""), "main_protein_name": row.get("protein_name", ""),
                  "main_decision": row.get("family_decision", ""), "main_confidence": row.get("assignment_confidence", ""),
                  "main_stage1_hypothesis": row.get("stage1_family_hypothesis", ""),
                  "source_orthogroup": row.get("Orthology_Orthogroup", ""),
                  "rbh_reference_sets": ";".join(sorted({r["reference_set"] for r in candidate["rbh"]})),
                  "rbh_statuses": ";".join(sorted({r["status"] for r in candidate["rbh"]})),
                  "rbh_supporting_queries": len({(r["reference_set"], r["reference_id"]) for r in candidate["rbh"]}),
                  "sequence_relationship": ";".join(sorted(candidate["relationships"])),
                  "family_conflict": conflict, "phylogeny_families": ";".join(sorted(claims)),
                  "record_role": "target_candidate", "length_aa": len(seq), "sequence_sha256": sequence_hash(seq),
                  "original_header": ";".join(f"{candidate['species']}|{rid}" for rid in sorted(candidate["raw_ids"]))}
        records.append((record, seq))
        for family in claims:
            family_members[family].append((record, seq))
    refs, ref_rows = {}, {}
    for reference_set, row in sorted(sets.items()):
        path = input_path(row["file"], args.homology_inputs, ["reference_fasta", "refrence_fasta"],
                          hc["input_sha256"], sources, f"references:{reference_set}")
        refs[reference_set] = fasta(path)
    for row in mappings:
        key = (row["reference_set"], row["reference_id"])
        if key in ref_rows:
            raise ValueError(f"Duplicate curated reference mapping: {key}")
        if key[0] not in refs or key[1] not in refs[key[0]]:
            raise ValueError(f"Reference mapping missing from original FASTA: {key}")
        ref_rows[key] = row
    for reference_set, seqs in refs.items():
        if {rid for s, rid in ref_rows if s == reference_set} != set(seqs):
            raise ValueError(f"Incomplete curated reference mapping: {reference_set}")
    for row in accepted_rbh:
        if (row["reference_set"], row["reference_id"]) not in ref_rows:
            raise ValueError("Accepted RBH has no curated query mapping")
    summaries, prepared_families = [], []
    for family, definition in sorted(families.items()):
        members = family_members[family]
        family_refs = {}
        for reference_set in definition["reference_sets"]:
            for rid, seq in refs.get(reference_set, {}).items():
                mapping = ref_rows[(reference_set, rid)]
                key = (mapping["reference_proteome"], rid, sequence_hash(seq))
                if key not in family_refs:
                    family_refs[key] = (mapping, seq, set())
                family_refs[key][2].add(reference_set)
        summary = {"family": family, "canonical_family": definition["family"],
                   "reference_sets": ";".join(definition["reference_sets"]), "targets": len(members),
                   "rbh_targets": sum(bool(r["rbh_reference_sets"]) for r, seq in members),
                   "main_accepted_targets": sum(r["main_accepted"] for r, seq in members),
                   "curated_references": len(family_refs),
                   "status": "prepared" if members else "not_prepared_no_target_candidates"}
        summaries.append(summary)
        if not members:
            continue
        members = list(members)
        for n, ((species, rid, sha), (mapping, seq, reference_sets)) in enumerate(sorted(family_refs.items()), 1):
            tip = {"sequence_id": f"R{n:06d}", "panel_label": f"reference_{species}", "species": species,
                   "original_ids": rid, "evidence_source": "curated_reference", "record_role": "curated_reference",
                   "reference_sets": ";".join(sorted(reference_sets)), "reference_source_assignment": mapping["source_assignment"],
                   "reference_mapping_method": mapping["mapping_method"], "source_orthogroup": "",
                   "length_aa": len(seq), "sequence_sha256": sha, "original_header": rid}
            members.append((tip, seq))
        for tip, seq in members:
            if "*" in seq or not AMINO_ACIDS.intersection(seq):
                raise ValueError(f"Phylogeny input has internal stops/no informative residues: {family}/{tip['sequence_id']}")
        if len(members) < 4 or len({seq for tip, seq in members}) < 4:
            summary["status"] = "not_prepared_too_few_distinct_sequences"
            continue
        prepared_families.append((family, members))
    counts = Counter(row["evidence_source"] for row, seq in records)
    summary = {"target_species": sorted(species_scope), "audit_rows_outside_target_cohort": len(main_rows) - len(main),
               "rbh_unique_input_candidates": len(global_rbh), "reconciled_candidates": len(records),
               "candidate_sources": dict(counts), "rbh_main_decisions": dict(Counter(r["main_decision"] for r in reconciliation_rows)),
               "rbh_negative_excluded": sum(r["negative_excluded"] for r in reconciliation_rows),
               "sequence_conflicts_or_alternative_isoforms": len(conflicts),
               "prepared_families": len(prepared_families), "prepared_phylogeny_records": sum(len(m) for f, m in prepared_families)}
    return {"records": records, "prepared_families": prepared_families, "family_summary": summaries,
            "reconciliation": reconciliation_rows, "excluded": excluded, "conflicts": conflicts,
            "negative_gates": [r for r in exclusions if r["Specie"] in species_scope],
            "accepted_rbh_evidence": accepted_rbh, "sources": sources.records, "summary": summary,
            "preparation_checks": preparation_checks}


def write_table(path, rows, empty_fields):
    fields = list(dict.fromkeys(k for row in rows for k in row)) or empty_fields
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_fasta(path, records):
    with path.open("w", encoding="ascii", newline="\n") as handle:
        for row, seq in records:
            handle.write(f">{row['sequence_id']}\n")
            for n in range(0, len(seq), 80):
                handle.write(seq[n:n + 80] + "\n")


def prepare(out, data):
    out = out.resolve()
    for source in data["sources"]:
        if digest(source["path"]) != source["sha256"]:
            raise ValueError(f"Source changed during reconciliation: {source['path']}")
    # Never overwrite an existing library, current tree run, input or source file.
    out.mkdir(parents=True, exist_ok=False)
    write_fasta(out / "reconciled_candidates.faa", data["records"])
    write_table(out / "candidates.tsv", [r for r, s in data["records"]], ["sequence_id"])
    for name, key, fields in [
        ("family_summary.tsv", "family_summary", ["family", "status"]),
        ("rbh_main_reconciliation.tsv", "reconciliation", ["species", "original_id"]),
        ("excluded_candidates.tsv", "excluded", ["species", "protein_id", "source", "reason"]),
        ("sequence_conflicts.tsv", "conflicts", ["species", "original_id", "relationship"]),
        ("negative_exclusion_gates.tsv", "negative_gates", ["Specie", "protein_id", "exclusion_reason"]),
        ("rbh_supporting_evidence.tsv", "accepted_rbh_evidence", ["species", "reference_set", "reference_id"])]:
        write_table(out / name, data[key], fields)
    families = []
    for family, members in data["prepared_families"]:
        folder = out / family
        folder.mkdir()
        write_fasta(folder / "input.fasta", members)
        write_table(folder / "tips.tsv", [r for r, s in members], ["sequence_id", "panel_label"])
        families.append({"family": family, "records": len(members),
                         "distinct_sequences": len({seq for row, seq in members})})
    manifest = {"schema_version": "1.0", "library_kind": "reconciled_candidates",
                "created_utc": datetime.now(timezone.utc).isoformat(), "preserve_identical_tips": True,
                "families": families, "summary": data["summary"], "sources": data["sources"],
                "main_preparation_checks": data["preparation_checks"],
                "interpretation": "Source annotations are unchanged. RBH-only inputs are exploratory phylogenetic candidates. Family conflicts are retained explicitly; membership is not an identity assignment.",
                "files": {str(p.relative_to(out)).replace("\\", "/"): digest(p)
                          for p in sorted(out.rglob("*")) if p.is_file()}}
    (out / "library_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    check(out)


def check(out):
    manifest = load_json(out / "library_manifest.json")
    if manifest.get("library_kind") != "reconciled_candidates":
        raise ValueError("Expected a reconciled-candidate library")
    for name, sha in manifest["files"].items():
        if digest(out / safe_relative(name)) != sha:
            raise ValueError(f"Prepared file changed: {name}")
    for folder, fasta_name, tips_name in [(out, "reconciled_candidates.faa", "candidates.tsv")] + [
            (out / row["family"], "input.fasta", "tips.tsv") for row in manifest["families"]]:
        seqs = fasta(folder / fasta_name, allow_empty=True)
        tips = read_table(folder / tips_name)
        if len(tips) != len(seqs) or {r["sequence_id"] for r in tips} != set(seqs):
            raise ValueError(f"Tip membership changed: {folder}")
        for row in tips:
            seq = seqs[row["sequence_id"]]
            if sequence_hash(seq) != row["sequence_sha256"] or len(seq) != int(row["length_aa"]):
                raise ValueError(f"Tip residues changed: {folder}/{row['sequence_id']}")
    print(f"PASS: {len(manifest['families'])} reconciled phylogeny families; hashes, sequences and tip mappings verified")
    return manifest


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--prepare", action="store_true", help="Write a new, independent library; existing outputs are refused")
    mode.add_argument("--check", action="store_true", help="Validate the prepared --out library without opening source workflows")
    mode.add_argument('--dry-run',action='store_true',help='Dataset-free structural deployment validation')
    parser.add_argument('--clean-state',action='store_true')
    parser.add_argument('--report-dir',type=Path)
    parser.add_argument("--homology-run", type=Path, default=BUNDLE / "homology_retrieval/results/homology_run")
    parser.add_argument("--homology-inputs", type=Path, default=BUNDLE / "homology_retrieval")
    parser.add_argument("--audit", type=Path, help='Explicit completed main audit directory; no historical audit is assumed')
    parser.add_argument("--handoff", type=Path, help="Main workflow handoff containing detection/<species>/prepared; auto-read from annotation provenance")
    parser.add_argument("--crosswalk", type=Path, default=ROOT / "config/family_crosswalk.json")
    parser.add_argument("--out", type=Path, default=BUNDLE / "independent_family_phylogeny/library_reconciled")
    args = parser.parse_args(argv)
    try:
        if args.dry_run:
            sys.path.insert(0,str(BUNDLE/'integrin_adhesome_analysis/resource_library/scripts'))
            from deployment_validation import main as validate_deployment
            flags=['--dry-run','--workflow','reconciliation']
            if args.clean_state:flags.append('--clean-state')
            if args.report_dir:flags+=['--report-dir',str(args.report_dir)]
            return validate_deployment(flags)
        if args.check:
            check(args.out.resolve())
            return 0
        if args.audit is None:
            raise ValueError('unidentified_input: stage=reconciliation.audit; supply --audit with a completed current audit. Use --dry-run --clean-state for dataset-free validation.')
        if args.prepare and args.out.exists():
            raise ValueError(f"Output exists; use --check or select a new --out: {args.out}")
        if args.prepare:
            out = args.out.resolve()
            protected = [args.homology_run.resolve(), args.audit.resolve()]
            if args.handoff:
                protected.append(args.handoff.resolve())
            protected += [args.homology_inputs.resolve() / name for name in (
                "target_proteomes", "schistosome_proteomes", "reference_fasta", "refrence_fasta",
                "reference_proteomes", "reference_proteome")]
            if any(out == path or path in out.parents for path in protected):
                raise ValueError("--out must be outside source inputs, source runs and the main handoff")
        data = reconcile(args)
        print(json.dumps(data["summary"], indent=2))
        if args.prepare:
            prepare(args.out, data)
            print(f"Prepared independent library: {args.out.resolve()}")
        else:
            print("PLAN ONLY: sources checked; no files written or analyses launched. Add --prepare to export.")
        return 0
    except (OSError, ValueError, KeyError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
