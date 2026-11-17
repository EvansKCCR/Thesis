"""Separate, evidence-feature machine-learning classifier for adhesome candidates.

Training uses calibration control and reference-family labels. FASTA-derived
motifs and verified raw scan features are identical at training and inference.
Validation is evaluated separately and is never used for fitting.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import sys
from pathlib import Path

import joblib
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_extraction import DictVectorizer
from sklearn.metrics import (accuracy_score, confusion_matrix, f1_score,
                             precision_score, recall_score)
from sklearn.pipeline import Pipeline

VERSION = "5.8"
FEATURE_SCHEMA = "sequence_family_gate_v5_8_membership_calibrated_raw_vote_family"
_BUNDLE = Path(__file__).resolve().parents[2]
_SETTINGS_PATH = _BUNDLE / "config/ml_classifier/classifier.json"
_SETTINGS = json.loads(_SETTINGS_PATH.read_text(encoding="utf-8"))
MAX_PREDICTION_RECORDS = int(_SETTINGS["maximum_prediction_records_per_invocation"])
CATEGORIES = ("specific", "broad", "domain_only", "unresolved", "non_adhesome")
TOKEN = re.compile(r"[A-Za-z][A-Za-z0-9_.:+/-]*|\d+(?:\.\d+)?")
EVIDENCE_KEY = re.compile(
    r"pfam|cdd|motif|topolog|localiz|signalp|deeptmhmm|deeploc|architecture|"
    r"domain|orthogroup|orthology|evidence|sequence_length|protein_length|length|"
    r"composition|tm_helix|transmembrane|signal_peptide|posthoc_(?:pfam|cdd|evidence)", re.I)
LEAK_KEY = re.compile(
    r"(^|_)(control_class|expected|truth|label|supported|provisional|predicted|"
    r"assigned|putative|family|protein_name|grade|confidence|status|proposal|"
    r"exclusion|recruited|annotated|identity|accession|protein_id|original_id|"
    r"posthoc_evidence_tier)(_|$)", re.I)
ID_FIELDS = ("protein_id", "original_id", "query_id", "sequence_id", "query", "seqid", "ID", "id")
TIER_TARGET_FIELDS = {"family_support_level", "stage2_category", "annotation_category", "grade", "supported"}
SEQUENCE_ALIASES = ("original_id", "protein_id", "query_id", "sequence_id", "query", "seqid", "ID", "id")
AA = "ACDEFGHIKLMNPQRSTVWY"
DOMAIN_COORD = re.compile(r"(?P<name>PF\d{5}(?:\.\d+)?|[A-Za-z][A-Za-z0-9_.+/-]*)\s*:\s*(?P<start>\d+)\s*-\s*(?P<end>\d+)")
PFAM_ID = re.compile(r"PF\d{5}(?:\.\d+)?", re.I)
NEGATIVE_LABELS = {"negative", "hard_negative", "background_negative", "non_adhesome", "negative_control"}
POSITIVE_LABELS = {"positive", "adhesome_associated", "positive_control"}


def read_fasta(path: Path) -> list[dict[str, str]]:
    """Read protein FASTA without external tools; preserve the first header token as ID."""
    records = []
    header = None
    chunks = []
    def emit():
        if header is None:
            return
        sequence = "".join(chunks).upper().replace("*", "X").replace("-", "").replace(".", "")
        if not sequence:
            raise ValueError(f"Empty protein sequence in FASTA record {header!r}")
        invalid = sorted(set(sequence) - set(AA + "BXZJUO"))
        if invalid:
            raise ValueError(f"Invalid amino-acid character(s) {invalid} in FASTA record {header!r}")
        records.append({"protein_id": header.split()[0], "sequence": sequence})
    with path.open(encoding="utf-8-sig") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                emit()
                header = line[1:].strip()
                if not header:
                    raise ValueError("FASTA contains an empty header")
                chunks = []
            else:
                if header is None:
                    raise ValueError("FASTA sequence found before the first header")
                chunks.append(line)
    emit()
    ids = [r["protein_id"] for r in records]
    if not records:
        raise ValueError(f"No FASTA records found: {path}")
    if len(ids) != len(set(ids)):
        raise ValueError("FASTA input contains duplicate first-token identifiers")
    return records


def row_aliases(row: dict) -> set[str]:
    return {str(row.get(key, "")).strip() for key in SEQUENCE_ALIASES if str(row.get(key, "")).strip()}


def read_tsv(path: Path) -> list[dict[str, str]]:
    # Evidence JSON cells can exceed Python's 128 KiB CSV default.
    csv.field_size_limit(2**31 - 1)
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        fields = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def row_id(row: dict) -> str:
    for key in ID_FIELDS:
        value = str(row.get(key, "")).strip()
        if value:
            return value
    return ""


def tokens(text: str) -> list[str]:
    return [m.group(0).strip(".,;()[]{}") for m in TOKEN.finditer(text or "") if m.group(0).strip(".,;()[]{}")]


def _json_records(value) -> list[dict]:
    if isinstance(value, list):
        return [x for x in value if isinstance(x, dict)]
    try:
        parsed = json.loads(str(value or ""))
    except (TypeError, ValueError):
        return []
    if isinstance(parsed, dict):
        return [parsed]
    return [x for x in parsed if isinstance(x, dict)] if isinstance(parsed, list) else []


def _safe_float(value, default=None):
    try:
        number = float(value)
        return number if math.isfinite(number) else default
    except (TypeError, ValueError):
        return default


def _norm_feature_token(value: str) -> str:
    return re.sub(r"[^a-z0-9_.+-]+", "_", str(value or "").strip().lower()).strip("_")


def _domain_records(row: dict) -> list[dict]:
    hits = []
    for record in _json_records(row.get("ml_domain_hits_json")):
        name = str(record.get("id") or record.get("name") or "").strip()
        if not name:
            continue
        hits.append({"source": str(record.get("source") or "pfam").lower(), "name": name,
                     "start": _safe_float(record.get("start")), "end": _safe_float(record.get("end")),
                     "sequence_length": _safe_float(record.get("sequence_length")),
                     "score": _safe_float(record.get("score")),
                     "evalue": _safe_float(record.get("evalue")),
                     "profile_coverage": _safe_float(record.get("profile_coverage"))})
    have_profile_hits = any(PFAM_ID.fullmatch(str(hit["name"])) for hit in hits)
    text_fields = [(key, str(row.get(key, "") or "")) for key in (
        "Pfam_architecture", "NCBI_CDD_architecture", "architecture", "matched_domains",
        "NCBI_CDD_selected_hits", "CDD_architecture", "domain_architecture")]
    text_fields.extend(("ml_hypothesis_architecture", str(record.get("architecture", "")))
                       for record in _json_records(row.get("ml_hypothesis_architectures_json")))
    for key, value in text_fields:
        value = value.strip()
        if not value:
            continue
        source = "cdd" if re.search(r"cdd", key, re.I) else "pfam" if "pfam" in key.lower() else "domain"
        for match in DOMAIN_COORD.finditer(value):
            name = match.group("name")
            if source == "domain" and PFAM_ID.fullmatch(name):
                source = "pfam"
            if source == "pfam" and have_profile_hits and PFAM_ID.fullmatch(name):
                continue
            hits.append({"source": source, "name": name,
                         "start": float(match.group("start")), "end": float(match.group("end"))})
        if key == "matched_domains" and not have_profile_hits:
            for name in PFAM_ID.findall(value):
                hits.append({"source": "pfam", "name": name, "start": None, "end": None})
    unique = {}
    for hit in hits:
        name = hit["name"].strip()
        source = hit["source"]
        if source == "domain":
            source = "pfam" if PFAM_ID.fullmatch(name) else "cdd"
        start, end = hit.get("start"), hit.get("end")
        key = (source, name.lower(), start, end)
        unique.setdefault(key, {**hit, "source": source, "name": name, "start": start, "end": end})
    return sorted(unique.values(), key=lambda hit: (hit["source"], hit["start"] is None,
        hit["start"] if hit["start"] is not None else 0, hit["end"] if hit["end"] is not None else 0,
        hit["name"].lower()))


def _motif_records(row: dict) -> list[dict]:
    hits = []
    for record in _json_records(row.get("ml_motif_hits_json")):
        motif = str(record.get("motif") or record.get("pattern") or "").strip()
        if motif:
            hits.append({"motif": str(record.get("motif") or "").strip(),
                         "pattern": str(record.get("pattern") or "").strip(),
                         "state": str(record.get("state") or "").strip(),
                         "context": str(record.get("context") or "").strip(),
                         "start": _safe_float(record.get("start")), "end": _safe_float(record.get("end"))})
    return hits


def _homology_records(row: dict) -> list[dict]:
    records = []
    for record in _json_records(row.get("ml_homology_hits_json")):
        direction = str(record.get("direction", "")).lower()
        if direction not in {"forward", "reverse"}:
            continue
        records.append({"direction": direction,
                        "identity": _safe_float(record.get("identity")),
                        "alignment_length": _safe_float(record.get("alignment_length")),
                        "evalue": _safe_float(record.get("evalue")),
                        "bitscore": _safe_float(record.get("bitscore")),
                        "query_coverage": _safe_float(record.get("query_coverage")),
                        "subject_coverage": _safe_float(record.get("subject_coverage"))})
    return records


def _domain_features(row: dict) -> tuple[dict[str, float], list[dict]]:
    hits = _domain_records(row)
    features: dict[str, float] = {"domain:total_hit_count": float(len(hits))}
    by_source_name: dict[tuple[str, str], list[dict]] = {}
    by_source: dict[str, list[dict]] = {}
    for hit in hits:
        source, name = hit["source"], hit["name"]
        norm = _norm_feature_token(name)
        by_source_name.setdefault((source, norm), []).append(hit)
        by_source.setdefault(source, []).append(hit)
        features[f"domain:present:{source}:{norm}"] = 1.0
    for (source, name), copies in by_source_name.items():
        features[f"domain:copy_number:{source}:{name}"] = float(len(copies))
        scores = [hit["score"] for hit in copies if hit.get("score") is not None]
        evalues = [hit["evalue"] for hit in copies if hit.get("evalue") is not None]
        coverages = [hit["profile_coverage"] for hit in copies if hit.get("profile_coverage") is not None]
        if scores:
            features[f"domain:best_score:{source}:{name}"] = max(scores)
        if evalues:
            features[f"domain:best_neglog10_evalue:{source}:{name}"] = max(
                -math.log10(max(min(value, 1.0), 1e-300)) for value in evalues)
        if coverages:
            features[f"domain:best_profile_coverage:{source}:{name}"] = max(coverages)
    for source, source_hits in by_source.items():
        positioned = [x for x in source_hits if x.get("start") is not None]
        ordered = positioned if positioned else source_hits
        ordered = sorted(ordered, key=lambda x: (x.get("start") is None,
                         x.get("start") if x.get("start") is not None else source_hits.index(x)))
        names = [_norm_feature_token(x["name"]) for x in ordered]
        features[f"domain:source_hit_count:{source}"] = float(len(source_hits))
        for left, right in zip(names, names[1:]):
            pair = f"domain:order_pair:{source}:{left}>{right}"
            features[pair] = features.get(pair, 0.0) + 1.0
        for hit in positioned:
            length = next((_safe_float(row.get(key)) for key in ("length", "sequence_length", "protein_length", "query_length", "ml_seq_len_aa")
                           if _safe_float(row.get(key)) and _safe_float(row.get(key)) > 0), None)
            length = length or hit.get("sequence_length")
            if length:
                fraction = float(hit["start"]) / length
                region = "n_term" if fraction < 0.33 else "middle" if fraction < 0.67 else "c_term"
                features[f"domain:position:{source}:{_norm_feature_token(hit['name'])}:{region}"] = 1.0
    return features, hits


def _motif_features(row: dict) -> tuple[dict[str, float], list[dict]]:
    hits = _motif_records(row)
    features: dict[str, float] = {"motif:site_count": float(len(hits))}
    counts: dict[str, int] = {}
    for hit in hits:
        for key in ("motif", "pattern"):
            token = _norm_feature_token(hit.get(key, ""))
            if token:
                feature = f"motif:{key}:{token}"
                features[feature] = 1.0
                counts[feature] = counts.get(feature, 0) + 1
        for key in ("state", "context"):
            token = _norm_feature_token(hit.get(key, ""))
            if token:
                features[f"motif:{key}:{token}"] = 1.0
    for key, count in counts.items():
        features[f"{key}:copy_number"] = float(count)
    for key, value in row.items():
        low_key = key.lower()
        if "motif" not in low_key or "family" in low_key:
            continue
        if not str(value or "").strip():
            continue
        numeric = _safe_float(value)
        normalized = _norm_feature_token(key)
        if numeric is not None:
            features[f"motif:metric:{normalized}"] = numeric
        elif low_key in {"motif_check", "motif_evidence", "motif_details"}:
            for token in tokens(str(value)):
                features[f"motif:raw:{normalized}:{_norm_feature_token(token)}"] = 1.0
    return features, hits


def _homology_features(row: dict) -> tuple[dict[str, float], list[dict]]:
    hits = _homology_records(row)
    features: dict[str, float] = {}
    for direction in ("forward", "reverse"):
        directional = [x for x in hits if x["direction"] == direction]
        features[f"homology:{direction}:hit_count"] = float(len(directional))
        features[f"homology:{direction}:tested"] = float(str(row.get(f"ml_{direction}_tested", "")).lower() in {"1", "true", "yes"})
        if directional:
            features[f"homology:{direction}:best_identity"] = max(x["identity"] for x in directional if x.get("identity") is not None) if any(x.get("identity") is not None for x in directional) else 0.0
            features[f"homology:{direction}:best_bitscore"] = max(x["bitscore"] for x in directional if x.get("bitscore") is not None) if any(x.get("bitscore") is not None for x in directional) else 0.0
            features[f"homology:{direction}:best_alignment_length"] = max(x["alignment_length"] for x in directional if x.get("alignment_length") is not None) if any(x.get("alignment_length") is not None for x in directional) else 0.0
            evalues = [x["evalue"] for x in directional if x.get("evalue") is not None]
            if evalues:
                features[f"homology:{direction}:best_neglog10_evalue"] = max(-math.log10(max(min(e, 1.0), 1e-300)) for e in evalues)
            for metric in ("query_coverage", "subject_coverage"):
                values = [x[metric] for x in directional if x.get(metric) is not None]
                if values:
                    features[f"homology:{direction}:best_{metric}"] = max(values)
    fwd = [x for x in hits if x["direction"] == "forward"]
    rev = [x for x in hits if x["direction"] == "reverse"]
    features["homology:directional_concordance"] = float(bool(fwd) and bool(rev))
    reciprocal = _safe_float(row.get("ml_reciprocal_pair_count"), 0.0) or 0.0
    features["homology:reciprocal_pair_count"] = reciprocal
    status_records = _json_records(row.get("ml_homology_status_json"))
    for status in status_records:
        for field in ("similarity", "reciprocal", "homology_check", "related_similarity", "related_reciprocal"):
            value = _norm_feature_token(status.get(field, ""))
            if value and value not in {"not_tested", "", "none", "na"}:
                features[f"homology:status:{field}:{value}"] = 1.0
        for field in ("reference_count", "related_reference_count"):
            value = _safe_float(status.get(field))
            if value is not None:
                key = f"homology:{field}"
                features[key] = max(features.get(key, 0.0), value)
    return features, hits


def evidence_groups(row: dict) -> dict[str, bool]:
    _, domains = _domain_features(row)
    _, motifs = _motif_features(row)
    hfeatures, homology = _homology_features(row)
    homology_present = bool(homology) or any(key.startswith("homology:status:") for key in hfeatures)
    return {"domain": bool(domains), "motif": bool(motifs), "homology": homology_present}


def feature_dict(row: dict, config_vocab: set[str]) -> dict[str, float]:
    """Encode only raw motif, domain architecture and directional homology evidence."""
    features, domains = _domain_features(row)
    motif_features, motifs = _motif_features(row)
    homology_features, homology = _homology_features(row)
    features.update(motif_features)
    features.update(homology_features)
    groups = evidence_groups(row)
    for group, present in groups.items():
        features[f"evidence_group:{group}_present"] = float(present)
    features["evidence_group:informative_group_count"] = float(sum(groups.values()))
    observed = {_norm_feature_token(hit["name"]) for hit in domains}
    observed.update(_norm_feature_token(hit.get(key, "")) for hit in motifs for key in ("motif", "pattern"))
    for accession in config_vocab:
        if _norm_feature_token(accession) in observed:
            features[f"configured_model_observed:{_norm_feature_token(accession)}"] = 1.0
    features["configured_model_match_count"] = float(sum(_norm_feature_token(x) in observed for x in config_vocab))
    return features


def configured_vocab(bundle: Path) -> tuple[set[str], dict[str, str]]:
    files = {
        "family_architecture": bundle / "config/annotation/family_rules.json",
        "family_identity": bundle / "config/annotation/family_identity_rules.json",
        "stage2_evidence": bundle / "config/annotation/stage2_evidence_rules.json",
        "detection": bundle / "config/detection/family_rules.json",
        "motif_master": bundle / "config/motif_annotation/motif_master.json",
        "motif_context": bundle / "config/motif_annotation/context_rules.json",
        "adhesion_core_motifs": bundle / "config/motif_annotation/adhesion_core.json",
    }
    vocab: set[str] = set()
    hashes: dict[str, str] = {}
    for name, path in files.items():
        if not path.is_file():
            continue
        hashes[name] = sha(path)
        raw = path.read_text(encoding="utf-8-sig")
        # Capture profile IDs and diagnostic CDD names from architecture/identity rules.
        vocab.update(x.lower() for x in re.findall(r"PF\d{5}(?:\.\d+)?", raw, re.I))
        vocab.update(x.lower() for x in re.findall(r"\b(?:CH_[A-Za-z0-9_.-]+|PTKc_[A-Za-z0-9_.-]+|"
                                                  r"SH[23]_[A-Za-z0-9_.-]+|FERM_[A-Za-z0-9_.-]+|"
                                                  r"[A-Z][A-Za-z0-9_.-]{3,})\b", raw))
    if _SETTINGS_PATH.is_file():
        hashes["ml_classifier_settings"] = sha(_SETTINGS_PATH)
    # The post-hoc rules are currently implemented as code, so fingerprint that
    # source and consume its saved evidence columns when present.
    posthoc = bundle.parent / "schisto_pfam_cdd_reclassify/schisto_pfam_cdd_reclassify.py"
    if posthoc.is_file():
        hashes["posthoc_rules_code"] = sha(posthoc)
    return vocab, hashes


def merge_rows(base: list[dict], extra: list[dict]) -> list[dict]:
    by_id: dict[str, dict] = {}
    for row in extra:
        for ident in row_aliases(row):
            target = by_id.setdefault(ident, {})
            for key, value in row.items():
                if value is None or not str(value).strip() or (LEAK_KEY.search(key) and key not in TIER_TARGET_FIELDS):
                    continue
                # Preserve repeated evidence rather than silently replacing a hit.
                if target.get(key) and target[key] != value:
                    target[key] = f"{target[key]};{value}"
                else:
                    target[key] = value
    output = []
    for source in base:
        merged = dict(source)
        ident = row_id(source)
        extra = {}
        for alias in row_aliases(source):
            for key, value in by_id.get(alias, {}).items():
                if not extra.get(key):
                    extra[key] = value
                elif extra[key] != value:
                    extra[key] = f"{extra[key]};{value}"
        for key, value in extra.items():
            if not merged.get(key):
                merged[key] = value
            elif merged[key] != value:
                merged[key] = f"{merged[key]};{value}"
        output.append(merged)
    return output


def read_returned(directory: Path | None) -> list[dict]:
    if directory is None:
        return []
    if not directory.is_dir():
        raise ValueError(f"Returned-evidence directory not found: {directory}")
    records = []
    for path in sorted([*directory.rglob("*.tsv"), *directory.rglob("*.csv")]):
        try:
            if path.suffix.lower() == ".tsv":
                records.extend(read_tsv(path))
            else:
                with path.open(encoding="utf-8-sig", newline="") as handle:
                    records.extend(csv.DictReader(handle))
        except (OSError, csv.Error, UnicodeError):
            continue
    return records


def _append_json_record(container: dict, field: str, record: dict) -> None:
    values = container.setdefault(field, [])
    marker = json.dumps(record, sort_keys=True, separators=(",", ":"))
    if all(json.dumps(value, sort_keys=True, separators=(",", ":")) != marker for value in values):
        values.append(record)


def _reference_alias(value: str) -> str:
    """Normalize the harmless naming difference used by the two BLAST directions."""
    return re.sub(r"^([A-Za-z0-9]+)_component_", r"\1_", str(value or ""))


def read_stage1_evidence(directory: Path | None) -> tuple[list[dict], dict[str, str], list[dict]]:
    """Read Stage 1 Pfam hits, motif observations and both directional BLAST tables."""
    if directory is None:
        return [], {}, []
    if not directory.is_dir():
        raise ValueError(f"Stage 1 evidence directory not found: {directory}")
    by_id: dict[str, dict] = {}
    hashes: dict[str, str] = {}
    scan_status: list[dict] = []

    def record_for(ident: str) -> dict:
        ident = str(ident or "").strip()
        if not ident:
            return {}
        return by_id.setdefault(ident, {"protein_id": ident})

    for species_dir in sorted(path for path in directory.iterdir() if path.is_dir()):
        species = species_dir.name
        domain_file = species_dir / "domains/domains.csv"
        motif_file = species_dir / "candidates/motif_observations.csv"
        forward_file = species_dir / "homology/forward.tsv"
        reverse_file = species_dir / "homology/reverse.tsv"
        scan_status.append({"species": species,
            "domain_tested": (species_dir / "domains/scan.json").is_file() or domain_file.is_file(),
            "motif_tested": motif_file.is_file(),
            "forward_tested": forward_file.is_file(), "reverse_tested": reverse_file.is_file()})
        for path in (domain_file, motif_file, forward_file, reverse_file):
            if path.is_file():
                hashes[path.relative_to(directory).as_posix()] = sha(path)
        if domain_file.is_file():
            with domain_file.open(encoding="utf-8-sig", newline="") as handle:
                for hit in csv.DictReader(handle):
                    ident = str(hit.get("query_name") or hit.get("protein_id") or "").strip()
                    if not ident:
                        continue
                    accession = str(hit.get("pfam_accession") or hit.get("target_accession") or "").strip()
                    name = str(hit.get("curated_domain_family") or hit.get("description") or hit.get("target_name") or accession).strip()
                    _append_json_record(record_for(ident), "ml_domain_hits_json", {
                        "source": "pfam" if PFAM_ID.fullmatch(accession) else "domain",
                        "id": accession or name, "name": name,
                        "start": _safe_float(hit.get("ali_from") or hit.get("env_from")),
                        "end": _safe_float(hit.get("ali_to") or hit.get("env_to")),
                        "sequence_length": _safe_float(hit.get("query_length")),
                        "evalue": _safe_float(hit.get("i_evalue") or hit.get("full_evalue")),
                        "score": _safe_float(hit.get("domain_score") or hit.get("full_score")),
                        "profile_coverage": _safe_float(hit.get("profile_coverage"))})
                    sequence_length = _safe_float(hit.get("query_length"))
                    if sequence_length and sequence_length > 0:
                        record_for(ident)["ml_seq_len_aa"] = str(max(
                            sequence_length, _safe_float(record_for(ident).get("ml_seq_len_aa"), 0.0) or 0.0))
        if motif_file.is_file():
            with motif_file.open(encoding="utf-8-sig", newline="") as handle:
                for hit in csv.DictReader(handle):
                    ident = str(hit.get("protein_id") or "").strip()
                    if not ident:
                        continue
                    _append_json_record(record_for(ident), "ml_motif_hits_json", {
                        "motif": hit.get("motif", ""), "pattern": hit.get("pattern", ""),
                        "state": hit.get("state", ""), "context": hit.get("context", ""),
                        "start": _safe_float(hit.get("start")), "end": _safe_float(hit.get("end"))})
        directional_refs: dict[str, dict[str, set[str]]] = {}
        for direction, path in (("forward", forward_file), ("reverse", reverse_file)):
            if not path.is_file():
                continue
            with path.open(encoding="utf-8-sig", newline="") as handle:
                for values in csv.reader(handle, delimiter="\t"):
                    # Stage 1 BLAST_FIELDS is exactly:
                    # qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen
                    if len(values) != 12:
                        raise ValueError(f"Malformed Stage 1 {direction} BLAST row in {path}: expected 12 fields, found {len(values)}")
                    identity = _safe_float(values[2])
                    alignment_length = _safe_float(values[3])
                    evalue = _safe_float(values[8])
                    bitscore = _safe_float(values[9])
                    if identity is None or alignment_length is None or evalue is None or bitscore is None:
                        continue
                    query, subject = values[0], values[1]
                    qcov = scov = None
                    qstart, qend = _safe_float(values[4]), _safe_float(values[5])
                    sstart, send = _safe_float(values[6]), _safe_float(values[7])
                    qlen, slen = _safe_float(values[10]), _safe_float(values[11])
                    if qstart is not None and qend is not None and qlen and qlen > 0:
                        qcov = min(1.0, (abs(qend - qstart) + 1.0) / qlen)
                    if sstart is not None and send is not None and slen and slen > 0:
                        scov = min(1.0, (abs(send - sstart) + 1.0) / slen)
                    ident = subject if direction == "forward" else query
                    reference = query if direction == "forward" else subject
                    _append_json_record(record_for(ident), "ml_homology_hits_json", {
                        "direction": direction, "identity": identity,
                        "alignment_length": alignment_length, "evalue": evalue, "bitscore": bitscore,
                        "query_coverage": qcov, "subject_coverage": scov})
                    directional_refs.setdefault(ident, {"forward": set(), "reverse": set()})[direction].add(_reference_alias(reference))
        for ident, directions in directional_refs.items():
            count = len(directions["forward"] & directions["reverse"])
            record_for(ident)["ml_reciprocal_pair_count"] = str(count)

    # The Stage 1 family-hypothesis export has ordered profile coordinates and
    # explicit reciprocal/similarity statuses; keep only those raw evidence fields.
    hypothesis_file = directory.parent / "hypotheses.tsv"
    if hypothesis_file.is_file():
        hashes[hypothesis_file.relative_to(directory.parent).as_posix()] = sha(hypothesis_file)
        for hypothesis in read_tsv(hypothesis_file):
            ident = str(hypothesis.get("protein_id", "")).strip()
            if not ident:
                continue
            architecture = str(hypothesis.get("architecture", "")).strip()
            if architecture:
                _append_json_record(record_for(ident), "ml_hypothesis_architectures_json", {"architecture": architecture})
            status = {field: hypothesis.get(field, "") for field in (
                "similarity", "reciprocal", "homology_check", "related_similarity",
                "related_reciprocal", "reference_count", "related_reference_count") if str(hypothesis.get(field, "")).strip()}
            if status:
                _append_json_record(record_for(ident), "ml_homology_status_json", status)
    rows = []
    for ident, record in sorted(by_id.items()):
        for field in ("ml_domain_hits_json", "ml_motif_hits_json", "ml_homology_hits_json", "ml_hypothesis_architectures_json", "ml_homology_status_json"):
            if isinstance(record.get(field), list):
                record[field] = json.dumps(record[field], sort_keys=True, separators=(",", ":"))
        rows.append(record)
    return rows, hashes, scan_status


def apply_scan_status(rows: list[dict], scan_status: list[dict]) -> list[dict]:
    """Mark whether each evidence block was actually searched for a protein."""
    output = []
    for source in rows:
        row = dict(source)
        ident = row_id(row)
        for status in scan_status:
            if ident.startswith(status["species"] + "__"):
                for field in ("domain_tested", "motif_tested", "forward_tested", "reverse_tested"):
                    row[f"ml_{field}"] = str(bool(status[field])).lower()
                break
        output.append(row)
    return output


def load_branch(branch: Path, returned: Path | None, bundle: Path) -> tuple[list[dict], dict[str, str]]:
    report = branch / "reports/predictions.tsv"
    if not report.is_file():
        raise ValueError(f"Missing completed branch report: {report}")
    vocab, hashes = configured_vocab(bundle)
    rows = read_tsv(report)
    # Merge evidence exports using exact protein IDs. Do not merge label tables.
    candidates = {branch / "downstream/final/audited_family_assignments.tsv",
                  branch / "posthoc_refinement/Stage_2_posthoc_refinement.tsv"}
    candidates.update(branch.rglob("Stage_2_posthoc_refinement.tsv"))
    hashes[f"{branch.name}_predictions"] = sha(report)
    for path in sorted(candidates):
        if path.is_file():
            rows = merge_rows(rows, read_tsv(path))
            hashes[f"{branch.name}_{path.name}_{len(hashes)}"] = sha(path)
    stage1_rows, stage1_hashes, scan_status = read_stage1_evidence(branch / "handoff/detection")
    rows = apply_scan_status(merge_rows(rows, stage1_rows), scan_status)
    for name, digest in stage1_hashes.items():
        hashes[f"{branch.name}_stage1_{name}"] = digest
    returned_rows = read_returned(returned)
    rows = merge_rows(rows, returned_rows)
    if returned:
        for path in sorted([*returned.rglob("*.tsv"), *returned.rglob("*.csv")]):
            hashes[f"{branch.name}_returned_{path.relative_to(returned).as_posix()}"] = sha(path)
    return rows, hashes


def positive_tier(row: dict) -> str:
    """Pipeline-derived tier target; saved as pseudo-label, never gold truth."""
    supported = str(row.get("supported", "")).lower() in {"true", "1", "yes"}
    broad_values = " ".join(str(row.get(k, "")) for k in ("family_support_level", "annotation_category", "stage2_category")).lower()
    if supported:
        return "broad" if "broad" in broad_values or "superfamily" in broad_values else "specific"
    raw = " ".join(str(v) for k, v in row.items()
                   if EVIDENCE_KEY.search(k) and not LEAK_KEY.search(k)
                   and not re.search(r"length|composition", k, re.I)).strip()
    stage1 = str(row.get("stage1_recruited", "")).lower() in {"true", "1", "yes"}
    if stage1 or raw:
        return "domain_only"
    return "unresolved"


def make_feature_matrix(rows: list[dict], vocab: set[str]):
    return [feature_dict(row, vocab) for row in rows]


def random_forest(settings: dict, jobs: int) -> RandomForestClassifier:
    """Remove descriptive config metadata before calling sklearn's estimator."""
    parameters = dict(settings)
    algorithm = parameters.pop("algorithm", "RandomForestClassifier")
    parameters.pop("classes", None)  # Tier label inventory is documentation only.
    if algorithm != "RandomForestClassifier":
        raise ValueError(f"Unsupported configured classifier algorithm: {algorithm}")
    return RandomForestClassifier(**parameters, n_jobs=jobs)


def class_metrics(truth: list[str], predicted: list[str], labels: list[str]) -> dict:
    return {
        "records": len(truth),
        "accuracy": accuracy_score(truth, predicted) if truth else None,
        "macro_precision": precision_score(truth, predicted, labels=labels, average="macro", zero_division=0) if truth else None,
        "macro_recall": recall_score(truth, predicted, labels=labels, average="macro", zero_division=0) if truth else None,
        "macro_f1": f1_score(truth, predicted, labels=labels, average="macro", zero_division=0) if truth else None,
        "confusion_matrix_labels": labels,
        "confusion_matrix": confusion_matrix(truth, predicted, labels=labels).tolist() if truth else [],
    }


def normalized_control_class(row: dict) -> str | None:
    label = re.sub(r"[^a-z0-9]+", "_", str(row.get("control_class", "")).strip().lower()).strip("_")
    if label in POSITIVE_LABELS:
        return "adhesome_associated"
    if label in NEGATIVE_LABELS:
        return "non_adhesome"
    return None


def coverage_by_class(rows: list[dict]) -> dict:
    coverage = {}
    for label in ("positive", "negative"):
        control = "adhesome_associated" if label == "positive" else "non_adhesome"
        group = [row for row in rows if normalized_control_class(row) == control]
        counts = {name: sum(evidence_groups(row)[name] for row in group)
                  for name in ("domain", "motif", "homology")}
        coverage[label] = {"records": len(group), "records_with_feature": counts,
                           "feature_fraction": {name: count / len(group) if group else None
                                                for name, count in counts.items()}}
    return coverage


def membership_calls(binary, vectors: list[dict], rows: list[dict], threshold: float,
                     calibration: dict | None = None) -> tuple[list[str], list[float], list[float]]:
    """Apply saved calibration; return calls, calibrated probabilities, raw diagnostics."""
    if calibration is None:
        raise ValueError("Membership calls require the saved v5.8 calibrator; raw votes are not calibrated probabilities")
    from probability_calibration import apply_binary_sigmoid
    classes = list(binary.named_steps["model"].classes_)
    raw_probabilities = binary.predict_proba(vectors)
    calibrated_probabilities = apply_binary_sigmoid(calibration, raw_probabilities, classes)
    calls, calibrated, raw = [], [], []
    for index, row in enumerate(rows):
        probability = float(calibrated_probabilities[index])
        raw_vote = float(raw_probabilities[index][classes.index("adhesome_associated")])
        calibrated.append(probability)
        raw.append(raw_vote)
        if not any(evidence_groups(row).values()):
            calls.append("unresolved")
        elif probability < threshold:
            calls.append("non_adhesome")
        else:
            calls.append("adhesome_associated")
    return calls, calibrated, raw


def train(args) -> int:
    from family_gate import train_family_gate
    return train_family_gate(args, sys.modules[__name__])


def rows_from_fasta(fasta_records: list[dict[str, str]], source_rows: list[dict] | None = None):
    source_rows = source_rows or []
    rows = []
    for fasta_row in fasta_records:
        matches = [row for row in source_rows if row_aliases(fasta_row) & row_aliases(row)]
        row = merge_rows([dict(matches[0])], matches[1:])[0] if matches else {"protein_id": fasta_row["protein_id"]}
        row.pop("sequence", None)
        row.pop("protein_sequence", None)
        row["input_fasta_id"] = fasta_row["protein_id"]
        if not row_id(row):
            raise ValueError(f"FASTA ID {fasta_row['protein_id']!r} cannot be used as a row identifier")
        rows.append(row)
    return rows


def prediction_evidence_summary(row: dict) -> dict[str, str]:
    domains = _domain_records(row)
    motifs = _motif_records(row)
    homology = _homology_records(row)
    copies: dict[tuple[str, str], int] = {}
    ordered: dict[str, list[dict]] = {}
    for hit in domains:
        key = (hit["source"], _norm_feature_token(hit["name"]))
        copies[key] = copies.get(key, 0) + 1
        ordered.setdefault(hit["source"], []).append(hit)
    order_text = []
    for source, hits in sorted(ordered.items()):
        hits.sort(key=lambda hit: (hit.get("start") is None,
                                   hit.get("start") if hit.get("start") is not None else 0))
        order_text.append(f"{source}:" + ">".join(hit["name"] for hit in hits))
    directions = {}
    for direction in ("forward", "reverse"):
        hits = [hit for hit in homology if hit["direction"] == direction]
        identities = [hit["identity"] for hit in hits if hit.get("identity") is not None]
        evalues = [hit["evalue"] for hit in hits if hit.get("evalue") is not None]
        bitscores = [hit["bitscore"] for hit in hits if hit.get("bitscore") is not None]
        qcoverages = [hit["query_coverage"] for hit in hits if hit.get("query_coverage") is not None]
        scoverages = [hit["subject_coverage"] for hit in hits if hit.get("subject_coverage") is not None]
        directions[direction] = (f"tested={row.get('ml_' + direction + '_tested', 'unknown')};"
            f"hits={len(hits)};best_identity={max(identities):.3f}" if identities else
            f"tested={row.get('ml_' + direction + '_tested', 'unknown')};hits=0")
        if evalues:
            directions[direction] += f";best_evalue={min(evalues):.4g}"
        if bitscores:
            directions[direction] += f";best_bitscore={max(bitscores):.3f}"
        if qcoverages:
            directions[direction] += f";best_query_coverage={max(qcoverages):.3f}"
        if scoverages:
            directions[direction] += f";best_subject_coverage={max(scoverages):.3f}"
    group_presence = evidence_groups(row)
    motif_text = ";".join(f"{hit.get('motif') or hit.get('pattern')}:{hit.get('state') or 'detected'}:{hit.get('context') or 'context_unknown'}"
                           for hit in motifs)
    return {"domain_architecture": ";".join(f"{hit['source']}:{hit['name']}" +
            (f"@{int(hit['start'])}-{int(hit['end'])}" if hit.get("start") is not None and hit.get("end") is not None else "")
            for hit in domains),
            "domain_copy_number": ";".join(f"{source}:{name}={count}" for (source, name), count in sorted(copies.items())),
            "domain_order": ";".join(order_text), "motif_evidence": motif_text,
            "forward_homology": directions["forward"], "reverse_homology": directions["reverse"],
            "reciprocal_homology_pair_count": str(_safe_float(row.get("ml_reciprocal_pair_count"), 0.0) or 0),
            "feature_groups_present": ";".join(name for name, present in group_presence.items() if present),
            "core_feature_group_count": str(sum(group_presence.values()))}


def predict(args) -> int:
    from family_gate import predict_family_gate
    return predict_family_gate(args, sys.modules[__name__])


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    fit = commands.add_parser("train", help="Fit on a completed labelled calibration branch; optionally evaluate on validation")
    fit.add_argument("--calibration", type=Path, required=True, help="Calibration branch directory with reports/predictions.tsv")
    fit.add_argument("--validation", type=Path, help="Optional independent validation branch; metrics only, never fitted")
    fit.add_argument("--calibration-returned", type=Path, help="Legacy option: ignored; CDD fragments come from the frozen library")
    fit.add_argument("--validation-returned", type=Path, help="Legacy option: ignored; validation fragments never enter the reference bank")
    fit.add_argument("--baseline-validation-predictions", type=Path,
        help="Optional predictions.tsv from the prior classifier; comparison only and must cover the identical validation IDs")
    fit.add_argument("--out", type=Path, required=True, help="New model output directory")
    fit.add_argument("--jobs", type=int, default=-1)
    fit.set_defaults(func=train)
    pred = commands.add_parser("predict", help="Predict families from FASTA and verified Stage 1 evidence; 250 proteins per batch")
    pred.add_argument("--model", type=Path, required=True, help="Saved adhesome_ml.joblib or its model directory")
    pred.add_argument("--fasta", type=Path, help="Protein FASTA for sequence motif scans; must match the verified Stage 1 scan")
    pred.add_argument("--input-tsv", type=Path, help="Legacy option; unsupported by the sequence family gate")
    pred.add_argument("--evidence-dir", type=Path, help="Stage 1 detection directory with per-species domains, motif observations and forward/reverse BLAST tables")
    pred.add_argument("--returned", type=Path, help="Legacy option: returned CDD hits are not used as query features")
    pred.add_argument("--ids", type=Path, help="Optional newline-delimited selection from the input FASTA")
    pred.add_argument("--batch", action="store_true", help="Process any input size in sequential chunks of at most 250 proteins")
    pred.add_argument("--out", type=Path, required=True)
    pred.add_argument("--threshold", type=float, default=float(_SETTINGS["default_membership_probability_cutoff"]), help="Cutoff on the calibrated adhesome-membership probability")
    pred.add_argument("--family-vote-fraction-cutoff",type=float,default=None,
        help="Override the top-family raw-vote cutoff (default: configured 0.30; see validation sensitivity grid)")
    pred.add_argument("--family-vote-margin",type=float,default=None,
        help="Override the minimum top-versus-runner-up raw-vote margin (default: configured 0.30; see validation sensitivity grid)")
    pred.add_argument("--family-similarity-profile", choices=sorted(_SETTINGS["family_gate"]["family_similarity"]["profiles"]),
        default=None,
        help="Family-resolution E-value/coverage tier; defaults to the saved model policy and does not alter adhesome-membership calls")
    pred.set_defaults(func=predict)
    evaluation = commands.add_parser('evaluate', help='Evaluate a frozen CDD-bank family classifier on validation; no fitting')
    evaluation.add_argument('--model', type=Path, required=True)
    evaluation.add_argument('--validation', type=Path, required=True)
    evaluation.add_argument('--out', type=Path, required=True, help='New validation report directory')
    evaluation.add_argument('--family-vote-fraction-cutoff',type=float,default=None,
        help='Evaluate with an explicit top-family raw-vote cutoff; default comes from configuration')
    evaluation.add_argument('--family-vote-margin',type=float,default=None,
        help='Evaluate with an explicit minimum raw-vote margin; default comes from configuration')
    evaluation.add_argument('--family-similarity-profile', choices=sorted(_SETTINGS['family_gate']['family_similarity']['profiles']),
        default=None,help='Evaluate this family-resolution profile; defaults to the saved model policy and does not change membership calls')
    from family_gate import evaluate_family_gate
    evaluation.set_defaults(func=lambda args: evaluate_family_gate(args, sys.modules[__name__]))
    args = parser.parse_args(argv)
    try:
        if args.command == "train":
            if args.jobs == 0:
                raise ValueError("--jobs cannot be zero")
            if args.out.exists():
                raise ValueError("Choose a new output directory; existing model files are never overwritten")
            if args.validation and not (args.validation / "reports/predictions.tsv").is_file():
                raise ValueError("Validation input has no reports/predictions.tsv")
            if args.baseline_validation_predictions and not args.validation:
                raise ValueError("--baseline-validation-predictions requires --validation")
            if args.baseline_validation_predictions and not args.baseline_validation_predictions.is_file():
                raise ValueError("Baseline predictions TSV does not exist")
        elif args.command == "predict":
            if args.input_tsv is not None:
                raise ValueError("The current family gate requires --fasta and --evidence-dir; --input-tsv is unsupported")
            if args.input_tsv is None and args.returned is None and args.fasta is None:
                raise ValueError("Supply candidate proteins with --fasta and --evidence-dir")
            if args.evidence_dir is not None and not args.evidence_dir.is_dir():
                raise ValueError("--evidence-dir must be an existing Stage 1 detection directory")
            if not 0.0 < args.threshold < 1.0:
                raise ValueError("--threshold must be between 0 and 1")
            if args.family_vote_fraction_cutoff is not None and not 0.0 <= args.family_vote_fraction_cutoff <= 1.0:
                raise ValueError("--family-vote-fraction-cutoff must be between 0 and 1")
            if args.family_vote_margin is not None and not 0.0 <= args.family_vote_margin <= 1.0:
                raise ValueError("--family-vote-margin must be between 0 and 1")
            if args.out.exists():
                raise ValueError("Choose a new prediction output path; existing files are never overwritten")
            if not args.model.is_file() and not (args.model.is_dir() and (args.model / "adhesome_ml.joblib").is_file()):
                raise ValueError("--model must be a saved adhesome_ml.joblib file or a directory containing it")
        elif args.command == "evaluate":
            if args.family_vote_fraction_cutoff is not None and not 0.0 <= args.family_vote_fraction_cutoff <= 1.0:
                raise ValueError("--family-vote-fraction-cutoff must be between 0 and 1")
            if args.family_vote_margin is not None and not 0.0 <= args.family_vote_margin <= 1.0:
                raise ValueError("--family-vote-margin must be between 0 and 1")
        return args.func(args)
    except (OSError, ValueError, KeyError, csv.Error) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
