"""Turn sequence predicates in the independent Pfam/CDD post-hoc rules into ML features.

Only rule conditions and the domain tokens used by those conditions are read.
Return values (family names, tiers, proposals), orthogroup arguments, and final
annotation decisions are intentionally excluded.
"""
import ast
import hashlib
import json
import re
from functools import lru_cache
from pathlib import Path


def _direct_tuple_return(statements):
    """Whether this block directly returns a tuple, without descending into ifs."""
    for statement in statements:
        if isinstance(statement, ast.Return) and isinstance(statement.value, ast.Tuple):
            return True
        if isinstance(statement, ast.If):
            continue
        for child in ast.iter_child_nodes(statement):
            if isinstance(child, ast.Return) and isinstance(child.value, ast.Tuple):
                return True
            if isinstance(child, ast.If):
                continue
            if any(isinstance(x, ast.Return) and isinstance(x.value, ast.Tuple) for x in ast.walk(child)):
                return True
    return False


def _condition_nodes(path):
    source = Path(path).read_text(encoding="utf-8-sig")
    tree = ast.parse(source, filename=str(path))
    functions = {n.name: n for n in tree.body if isinstance(n, ast.FunctionDef)}
    rules = []
    pfam_tokens, cdd_tokens = set(), set()
    for function_name in ("specific", "broad"):
        function = functions.get(function_name)
        if function is None:
            continue
        for node in ast.walk(function):
            if not isinstance(node, ast.If):
                continue
            if not (_direct_tuple_return(node.body) or _direct_tuple_return(node.orelse)):
                continue
            expression = ast.unparse(node.test)
            rules.append({"feature_id": f"{function_name}_line_{node.lineno}",
                          "expression": expression, "source_function": function_name})
            for call in (x for x in ast.walk(node.test) if isinstance(x, ast.Call)):
                name = call.func.id if isinstance(call.func, ast.Name) else ""
                destination = cdd_tokens if name == "a" else pfam_tokens if name in {"b", "f", "pf_count"} else None
                if destination is not None:
                    destination.update(arg.value for arg in call.args
                                       if isinstance(arg, ast.Constant) and isinstance(arg.value, str))
    if not rules:
        raise ValueError(f"No sequence-based post-hoc rule conditions found in {path}")
    allowed = {"p", "cc", "og", "a", "b", "f", "pf_count"}
    for rule in rules:
        expression_tree = ast.parse(rule["expression"], mode="eval")
        names = {n.id for n in ast.walk(expression_tree) if isinstance(n, ast.Name)}
        unsupported = names - allowed
        if unsupported:
            raise ValueError(f"Post-hoc condition {rule['feature_id']} uses unsupported fields: {sorted(unsupported)}")
        if "og" in names:
            raise ValueError(f"Orthogroup-dependent post-hoc rule excluded: {rule['feature_id']}")
        rule["_code"] = compile(expression_tree, f"{path}:{rule['feature_id']}", "eval")
    return source, rules, sorted(pfam_tokens), sorted(cdd_tokens)


@lru_cache(maxsize=4)
def load_rules(path):
    source, rules, pfam_tokens, cdd_tokens = _condition_nodes(path)
    return hashlib.sha256(source.encode()).hexdigest(), rules, pfam_tokens, cdd_tokens


def _annotation_counts(cdd_features):
    counts = {}
    prefix, suffix = "cdd:annotation:", ":hit_count"
    for key, value in cdd_features.items():
        if key.startswith(prefix) and key.endswith(suffix):
            try:
                counts[key[len(prefix):-len(suffix)]] = max(0.0, float(value))
            except (TypeError, ValueError):
                continue
    return counts


def extract(row, cdd_features, source_path):
    """Compute rule, partial-signature, and copy-count features from raw sequence evidence."""
    source_sha, rules, pfam_tokens, cdd_tokens = load_rules(str(source_path))
    domain_hits = json.loads(row.get("ml_domain_hits_json") or "[]")
    pfams = [str(hit.get("id", "")).split(".")[0] for hit in domain_hits
             if str(hit.get("id", "")).upper().startswith("PF")]
    cdd_counts = _annotation_counts(cdd_features)
    cc = [token for token, count in cdd_counts.items() if count > 0]
    a = lambda *terms: any(term in cc for term in terms)
    b = lambda *terms: any(term in pfams for term in terms)
    f = lambda *terms: all(term in pfams for term in terms)
    def pf_count(values, token):
        return values.count(token)

    env = {"p": pfams, "cc": cc, "a": a, "b": b, "f": f, "pf_count": pf_count}
    features = {}
    hits = []
    for rule in rules:
        hit = bool(eval(rule["_code"], {"__builtins__": {}}, env))
        features["posthoc_rule:" + rule["feature_id"]] = float(hit)
        hits.append((rule["source_function"], hit))

    pfam_present = 0
    for token in pfam_tokens:
        count = pf_count(pfams, token)
        features[f"posthoc_pfam:{token}:copy_count"] = float(count)
        features[f"posthoc_pfam:{token}:present"] = float(count > 0)
        pfam_present += count > 0
    cdd_present = 0
    for token in cdd_tokens:
        # cdd_features uses the same normalized model-name token convention.
        normalized = re.sub(r"[^A-Za-z0-9_.-]+", "_", token).strip("_") or "unknown"
        count = cdd_counts.get(normalized, 0.0)
        features[f"posthoc_cdd:{normalized}:hit_count"] = float(count)
        features[f"posthoc_cdd:{normalized}:present"] = float(count > 0)
        cdd_present += count > 0

    specific_signals = sum(hit for kind, hit in hits if kind == "specific")
    broad_signals = sum(hit for kind, hit in hits if kind == "broad")
    features["posthoc_signature:pfam_present_count"] = float(pfam_present)
    features["posthoc_signature:cdd_present_count"] = float(cdd_present)
    features["posthoc_signature:specific_rule_signal_count"] = float(specific_signals)
    features["posthoc_signature:broad_rule_signal_count"] = float(broad_signals)
    features["posthoc_signature:multiple_specific_rule_signals"] = float(specific_signals > 1)
    features["posthoc_signature:specific_and_broad_signal"] = float(specific_signals > 0 and broad_signals > 0)
    features["posthoc_signature:rule_signal_count"] = float(specific_signals + broad_signals)
    features["posthoc_signature:rule_count"] = float(len(rules))
    return features, source_sha
