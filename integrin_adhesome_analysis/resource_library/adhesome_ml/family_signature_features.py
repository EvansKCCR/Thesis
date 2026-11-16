"""Crosswalk-derived CDD plus family-architecture combination features.

Only sequence observations and the configured family identity-rule predicates are
used. Crosswalk labels and validation outcomes are not read as query evidence.
"""
from collections import defaultdict
from functools import lru_cache
import hashlib
import json
import math
import re
from pathlib import Path


def _token(value):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value or "")).strip("_").casefold()


def _split(value):
    if value is None:
        return []
    return [part.strip() for part in re.split(r"[;|]", str(value))
            if part and part.strip() and str(part).strip().casefold() not in {"nan", "none", "null"}]


def _canonical(value, aliases):
    current = str(value or "").strip()
    seen = set()
    while current.casefold() in aliases and current.casefold() not in seen:
        seen.add(current.casefold())
        current = str(aliases[current.casefold()]).strip()
    return current


def _sheet_rows(worksheet):
    iterator = worksheet.iter_rows(values_only=True)
    headers = next(iterator, None)
    if not headers:
        return [], set()
    names = [str(value).strip() if value is not None else "" for value in headers]
    records = []
    for values in iterator:
        if not values or not any(value is not None and str(value).strip() for value in values):
            continue
        records.append({name: values[index] if index < len(values) else None
                        for index, name in enumerate(names) if name})
    return records, set(names)


@lru_cache(maxsize=4)
def _load_library_cached(path_text, rules_json, aliases_json):
    path = Path(path_text)
    if not path.is_file():
        raise ValueError(f"Adhesome family-rule crosswalk workbook is missing: {path}")
    try:
        from openpyxl import load_workbook
    except ImportError as exc:
        raise ValueError("Crosswalk features require the workflow's openpyxl dependency") from exc

    rules = json.loads(rules_json)
    aliases = json.loads(aliases_json)
    rule_by_id = {str(rule.get("id", "")).strip(): rule for rule in rules
                  if str(rule.get("id", "")).strip()}
    rules_by_family = defaultdict(list)
    for rule in rules:
        family = _token(_canonical(rule.get("family", ""), aliases))
        if family and rule.get("required_groups"):
            rules_by_family[family].append(rule)

    workbook = load_workbook(path, read_only=True, data_only=True)
    try:
        required_sheets = {"CDD_Review", "Full_Catalogue_5570", "Alias_Proposals"}
        missing_sheets = required_sheets - set(workbook.sheetnames)
        if missing_sheets:
            raise ValueError("Crosswalk workbook is missing sheets: " + ", ".join(sorted(missing_sheets)))
        review, review_headers = _sheet_rows(workbook["CDD_Review"])
        catalogue, catalogue_headers = _sheet_rows(workbook["Full_Catalogue_5570"])
        aliases_rows, alias_headers = _sheet_rows(workbook["Alias_Proposals"])
    finally:
        workbook.close()

    required_review = {"CDD accession", "CDD annotation", "Evidence level", "Review priority",
                       "Family-rule IDs", "Registry family IDs", "Family-identity Pfam anchors",
                       "Existing ML feature"}
    required_catalogue = {"CDD accession", "Annotation", "Original category",
                          "Original family pattern", "Pfam core anchor", "CDD classifier feature"}
    if missing := required_review - review_headers:
        raise ValueError("CDD_Review lacks required columns: " + ", ".join(sorted(missing)))
    if missing := required_catalogue - catalogue_headers:
        raise ValueError("Full_Catalogue_5570 lacks required columns: " + ", ".join(sorted(missing)))
    if "Status" not in alias_headers:
        raise ValueError("Alias_Proposals lacks its Status column")

    signatures = {}

    def add_mapping(accession, annotation, family_name, tier, evidence_level, source, explicit_rule_ids=()):
        family = _canonical(family_name, aliases)
        family_key = _token(family)
        accession_key = _token(accession)
        if not family_key or not accession_key:
            return
        family_rules = rules_by_family.get(family_key, [])
        selected_rules = [rule_by_id[rid] for rid in explicit_rule_ids
                          if rid in rule_by_id and
                          _token(_canonical(rule_by_id[rid].get("family", ""), aliases)) == family_key]
        if not selected_rules:
            selected_rules = family_rules
        for rule in selected_rules:
            rule_id = str(rule.get("id", "")).strip()
            if not rule_id or not rule.get("required_groups"):
                continue
            key = (family_key, rule_id, accession_key)
            item = signatures.setdefault(key, {
                "family": family, "family_token": family_key, "rule_id": rule_id,
                "rule_token": _token(rule_id), "accession": str(accession).strip(),
                "accession_token": accession_key, "annotation": str(annotation or "").strip(),
                "tiers": set(), "evidence_levels": set(), "sources": set()})
            item["tiers"].add(tier)
            item["evidence_levels"].add(str(evidence_level or "").strip())
            item["sources"].add(source)

    # The full catalogue supplies configured family-pattern/accession links. They
    # are only used in conjunction with a complete family-identity architecture.
    for record in catalogue:
        category = str(record.get("Original category") or "").strip().casefold()
        if category not in {"family_pattern_candidate", "family_pattern_ambiguous"}:
            continue
        accession = record.get("CDD accession")
        annotation = record.get("Annotation")
        tier = "configured_pattern" if category == "family_pattern_candidate" else "ambiguous_pattern"
        families = _split(record.get("Original family pattern"))
        anchors = {_token(_canonical(value, aliases)) for value in _split(record.get("Pfam core anchor"))}
        if anchors:
            families = [family for family in families if _token(_canonical(family, aliases)) in anchors]
        for family in families:
            add_mapping(accession, annotation, family, tier, category,
                        "Full_Catalogue_5570")

    allowed_priorities = {"P1_named_model_review", "P2_competing_identity_review",
                          "P3_architecture_review"}
    excluded_levels = {"broad_domain_class_only", "broad_pfam_domain_class_only",
                       "no_defensible_configured_family_link",
                       "cooccurrence_multiple_families_context_only",
                       "multiple_named_candidates_unresolved"}
    for record in review:
        priority = str(record.get("Review priority") or "").strip()
        evidence_level = str(record.get("Evidence level") or "").strip()
        if priority not in allowed_priorities or evidence_level in excluded_levels:
            continue
        accession = record.get("CDD accession")
        annotation = record.get("CDD annotation")
        listed_rule_ids = [rule_id for rule_id in _split(record.get("Family-rule IDs"))
                           if rule_id in rule_by_id]
        registry_families = _split(record.get("Registry family IDs"))
        anchor_families = _split(record.get("Family-identity Pfam anchors"))
        if anchor_families:
            anchors = {_token(_canonical(value, aliases)) for value in anchor_families}
            registry_families = [family for family in registry_families
                                 if _token(_canonical(family, aliases)) in anchors]
            if not registry_families:
                registry_families = anchor_families
        if registry_families:
            family_rules = defaultdict(list)
            for family in registry_families:
                key = _token(_canonical(family, aliases))
                for rule_id in listed_rule_ids:
                    if _token(_canonical(rule_by_id[rule_id].get("family", ""), aliases)) == key:
                        family_rules[key].append(rule_id)
                if not family_rules[key]:
                    family_rules[key] = [rule["id"] for rule in rules_by_family.get(key, [])]
            for family in registry_families:
                key = _token(_canonical(family, aliases))
                add_mapping(accession, annotation, family, priority, evidence_level,
                            "CDD_Review", family_rules.get(key, []))
        elif listed_rule_ids:
            for rule_id in listed_rule_ids:
                rule = rule_by_id[rule_id]
                add_mapping(accession, annotation, rule.get("family"), priority,
                            evidence_level, "CDD_Review", [rule_id])

    compiled = []
    for item in signatures.values():
        compiled.append({**item,
            "tiers": tuple(sorted(item["tiers"])),
            "evidence_levels": tuple(sorted(x for x in item["evidence_levels"] if x)),
            "sources": tuple(sorted(item["sources"]))})
    compiled.sort(key=lambda x: (x["family_token"], x["rule_id"], x["accession_token"]))
    if not compiled:
        raise ValueError("Crosswalk contains no family-rule/CDD combinations that match configured identity rules")

    alias_statuses = [str(row.get("Status") or "").strip() for row in aliases_rows]
    metadata = {
        "path": path.name,
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "review_rows": len(review),
        "catalogue_rows": len(catalogue),
        "signature_combinations": len(compiled),
        "mapped_families": len({item["family_token"] for item in compiled}),
        "mapped_identity_rules": len({item["rule_id"] for item in compiled}),
        "signature_combinations_by_review_priority": {
            priority: sum(priority in item["tiers"] for item in compiled)
            for priority in sorted(allowed_priorities)},
        "alias_proposal_rows_ignored": len(alias_statuses),
        "alias_proposal_statuses_ignored": dict(sorted((status, alias_statuses.count(status))
                                                        for status in set(alias_statuses))),
        "feature_policy": "CDD accession similarity is combined with configured family identity-rule architecture predicates; proposed family/class, same-protein family context, calibration counts, and alias proposals are not query features."}
    return tuple(compiled), metadata


def load_library(crosswalk_path, identity_rules, aliases):
    rules_json = json.dumps(identity_rules, sort_keys=True, separators=(",", ":"))
    aliases_json = json.dumps(aliases, sort_keys=True, separators=(",", ":"))
    return _load_library_cached(str(Path(crosswalk_path).resolve()), rules_json, aliases_json)


def extract(row, vector, signatures, search_max_evalue):
    """Return learned-model features and traceable observed CDD/architecture combos."""
    cutoff = float(search_max_evalue)
    best_by_accession = {}
    for hit in row.get("_family_similarity_hsps", []) or []:
        try:
            evalue = float(hit.get("evalue", 1.0))
        except (TypeError, ValueError):
            continue
        if evalue > cutoff:
            continue
        for accession in hit.get("accessions", []) or []:
            key = _token(accession)
            if not key:
                continue
            previous = best_by_accession.get(key)
            rank = (evalue, -float(hit.get("bitscore", 0.0)))
            if previous is None or rank < previous[0]:
                best_by_accession[key] = (rank, hit)

    output = {}
    evidence = []
    family_stats = defaultdict(lambda: {"accessions_by_tier": defaultdict(set),
        "full": 0, "partial_sum": 0.0, "combos": 0, "rules_with_hits": set(), "best_strength": 0.0})
    family_totals = defaultdict(lambda: {"accessions_by_tier": defaultdict(set), "combos": 0})

    for signature in signatures:
        family = signature["family_token"]
        rule_id = signature["rule_id"]
        accession = signature["accession_token"]
        prefix = f"identity_rule:{rule_id}:"
        observed = best_by_accession.get(accession)
        hit = observed[1] if observed else None
        has_cdd = bool(hit)
        required_fraction = float(vector.get(prefix + "required_satisfied_fraction", 0.0))
        required_all = float(vector.get(prefix + "required_all_satisfied", 0.0)) >= 0.5
        pattern_defined = float(vector.get(prefix + "specific_cdd_pattern_defined", 0.0)) >= 0.5
        pattern_ok = (not pattern_defined or
                      float(vector.get(prefix + "specific_cdd_pattern_present", 0.0)) >= 0.5)
        order_count = float(vector.get(prefix + "ordered_group_count", 0.0))
        order_ok = order_count == 0 or float(vector.get(prefix + "ordered_groups_satisfied", 0.0)) >= 0.5
        max_ok = float(vector.get(prefix + "maximum_constraints_satisfied", 1.0)) >= 0.5
        length_defined = prefix + "length_range_min" in vector
        length_ok = (not length_defined or float(vector.get(prefix + "length_in_range", 0.0)) >= 0.5)
        repeat_count = float(vector.get(prefix + "required_cdd_repeat_count", 0.0))
        repeat_ok = (repeat_count == 0 or
                     float(vector.get(prefix + "required_cdd_repeats_satisfied", 0.0)) >= 0.5)
        full = has_cdd and required_all and pattern_ok and order_ok and max_ok and length_ok and repeat_ok
        base = f"family_signature:combo:{family}:{signature['rule_token']}:{accession}:"
        output[base + "full_architecture"] = float(full)
        output[base + "partial_architecture_score"] = (
            required_fraction * float(has_cdd) if has_cdd else 0.0)

        stat = family_stats[family]
        total = family_totals[family]
        stat["combos"] += 1
        total["combos"] += 1
        if has_cdd:
            stat["rules_with_hits"].add(rule_id)
            strength = min(100.0, -math.log10(max(float(hit.get("evalue", 1.0)), 1e-100)))
            stat["best_strength"] = max(stat["best_strength"], strength)
            for tier in signature["tiers"]:
                stat["accessions_by_tier"][tier].add(accession)
                total["accessions_by_tier"][tier].add(accession)
            stat["partial_sum"] += required_fraction
            if full:
                stat["full"] += 1
            evidence.append({
                "family": signature["family"], "rule_id": rule_id,
                "cdd_accession": signature["accession"],
                "cdd_annotation": signature["annotation"],
                "crosswalk_tiers": list(signature["tiers"]),
                "crosswalk_evidence_levels": list(signature["evidence_levels"]),
                "cdd_evalue": float(hit.get("evalue", 1.0)),
                "cdd_bitscore": float(hit.get("bitscore", 0.0)),
                "cdd_query_coverage": float(hit.get("query_coverage", 0.0)),
                "cdd_subject_coverage": float(hit.get("subject_coverage", 0.0)),
                "identity_rule_required_fraction": required_fraction,
                "identity_rule_required_all_satisfied": required_all,
                "identity_rule_specific_pattern_satisfied": pattern_ok,
                "identity_rule_order_satisfied": order_ok,
                "identity_rule_constraints_satisfied": bool(max_ok and length_ok and repeat_ok),
                "full_architecture_combination": bool(full),
                "crosswalk_sources": list(signature["sources"])})

    for family, total in family_totals.items():
        stat = family_stats[family]
        prefix = f"family_signature:family:{family}:"
        for tier in ("P1_named_model_review", "P2_competing_identity_review",
                     "P3_architecture_review", "configured_pattern", "ambiguous_pattern"):
            output[prefix + tier + "_accession_count"] = float(
                len(total["accessions_by_tier"].get(tier, set())))
        output[prefix + "full_combo_count"] = float(stat["full"])
        output[prefix + "full_combo_fraction"] = stat["full"] / total["combos"] if total["combos"] else 0.0
        output[prefix + "partial_architecture_mean"] = stat["partial_sum"] / total["combos"] if total["combos"] else 0.0
        output[prefix + "rules_with_cdd_hits"] = float(len(stat["rules_with_hits"]))
        output[prefix + "best_cdd_neglog10_evalue"] = stat["best_strength"]

    evidence.sort(key=lambda item: (not item["full_architecture_combination"],
                                    item["family"], item["rule_id"], item["cdd_accession"]))
    return output, evidence
