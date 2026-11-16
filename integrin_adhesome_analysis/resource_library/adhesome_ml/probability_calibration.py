"""Leakage-aware sigmoid calibration for Random Forest score outputs.

Only out-of-fold predictions from labelled calibration controls are used to fit
the sigmoid mappings. Exact-sequence groups are kept intact across folds; these
group identifiers are provenance only and never enter the feature matrix.
"""
from collections import Counter
import hashlib
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedGroupKFold


class CalibrationUnavailable(ValueError):
    """Raised when calibration data cannot support an honest out-of-fold fit."""


def sequence_cluster_groups(sequences, executable="cd-hit", identity=0.98,
                            shorter_coverage=0.90, longer_coverage=0.90,
                            max_memory_mb=1024):
    """Group exact duplicates and, when available, near-identical full-length proteins.

    CD-HIT cluster labels are used only to keep related records in the same OOF fold.
    The returned group IDs are hashes of sequence digests, not protein identifiers.
    If CD-HIT is missing or fails, exact-sequence SHA-256 groups remain the safe fallback.
    """
    sequences = [re.sub(r"\s+", "", str(sequence)).upper() for sequence in sequences]
    if any(not sequence for sequence in sequences):
        raise CalibrationUnavailable("Cannot form leakage groups for an empty protein sequence")
    sequence_hashes = [hashlib.sha256(sequence.encode("utf-8")).hexdigest()
                       for sequence in sequences]

    def exact_fallback(reason):
        return ["exact:" + digest for digest in sequence_hashes], {
            "method": "exact_sequence_sha256_fallback",
            "executable": executable,
            "identity": None,
            "minimum_shorter_sequence_coverage": None,
            "minimum_longer_sequence_coverage": None,
            "fallback_reason": reason}

    tool = shutil.which(executable) if executable else None
    if not tool:
        return exact_fallback("CD-HIT executable not found")
    try:
        with tempfile.TemporaryDirectory(prefix="adhesome_oof_groups_") as temp:
            temp = Path(temp)
            input_fasta = temp / "sequences.faa"
            output_fasta = temp / "clusters.faa"
            with input_fasta.open("w", encoding="ascii", newline="\n") as handle:
                for index, sequence in enumerate(sequences):
                    handle.write(f">seq_{index:08d}\n")
                    for offset in range(0, len(sequence), 80):
                        handle.write(sequence[offset:offset + 80] + "\n")
            command = [tool, "-i", str(input_fasta), "-o", str(output_fasta),
                       "-c", str(float(identity)), "-aS", str(float(shorter_coverage)),
                       "-aL", str(float(longer_coverage)), "-g", "1", "-d", "0",
                       "-T", "1", "-M", str(int(max_memory_mb))]
            completed = subprocess.run(command, capture_output=True, text=True, timeout=600)
            cluster_file = Path(str(output_fasta) + ".clstr")
            if completed.returncode != 0 or not cluster_file.is_file():
                reason = (completed.stderr or completed.stdout or
                          f"CD-HIT exited with status {completed.returncode}").strip()[-1000:]
                return exact_fallback(reason)
            cluster_members = []
            current = []
            for line in cluster_file.read_text(encoding="utf-8").splitlines():
                if line.startswith(">Cluster "):
                    if current:
                        cluster_members.append(current)
                    current = []
                    continue
                match = re.search(r">seq_(\d{8})\.\.\.", line)
                if match:
                    current.append(int(match.group(1)))
            if current:
                cluster_members.append(current)
            if sorted(index for cluster in cluster_members for index in cluster) != list(range(len(sequences))):
                return exact_fallback("CD-HIT cluster report did not cover each sequence exactly once")
            groups = [""] * len(sequences)
            for members in cluster_members:
                digest_payload = "\n".join(sorted({sequence_hashes[index] for index in members}))
                group_id = "cluster:" + hashlib.sha256(digest_payload.encode("utf-8")).hexdigest()
                for index in members:
                    groups[index] = group_id
            return groups, {"method": "cdhit_grouped_out_of_fold",
                "executable": tool, "identity": float(identity),
                "minimum_shorter_sequence_coverage": float(shorter_coverage),
                "minimum_longer_sequence_coverage": float(longer_coverage),
                "max_memory_mb": int(max_memory_mb), "threads": 1,
                "fallback_reason": ""}
    except (OSError, subprocess.SubprocessError, UnicodeError, ValueError) as exc:
        return exact_fallback(f"CD-HIT grouping failed: {exc}")


def grouped_splits(labels, groups, max_splits=5, random_state=3817):
    labels = np.asarray(labels, dtype=str)
    groups = np.asarray(groups, dtype=str)
    if len(labels) != len(groups) or not len(labels):
        raise CalibrationUnavailable("Calibration labels/groups are empty or have different lengths")
    classes = sorted(set(labels))
    group_counts = {label: len(set(groups[labels == label])) for label in classes}
    possible = min(group_counts.values(), default=0)
    for n_splits in range(min(int(max_splits), possible), 1, -1):
        splitter = StratifiedGroupKFold(n_splits=n_splits, shuffle=True,
                                        random_state=int(random_state))
        try:
            folds = list(splitter.split(np.zeros(len(labels)), labels, groups))
        except ValueError:
            continue
        valid = True
        for train_idx, test_idx in folds:
            if set(labels[train_idx]) != set(classes) or set(labels[test_idx]) != set(classes):
                valid = False
                break
            if set(groups[train_idx]) & set(groups[test_idx]):
                valid = False
                break
        if valid:
            return folds
    raise CalibrationUnavailable(
        "Cannot create at least two group-disjoint folds with every class represented "
        f"in train and held-out partitions; distinct sequence groups by class: {group_counts}")


def cross_fitted_probabilities(estimator_factory, matrix, labels, groups, max_splits=5,
                               random_state=3817):
    """Return complete OOF predictions, exact fold assignments and split summary."""
    labels = np.asarray(labels, dtype=str)
    groups = np.asarray(groups, dtype=str)
    folds = grouped_splits(labels, groups, max_splits, random_state)
    classes = sorted(set(labels))
    oof = np.full((len(labels), len(classes)), np.nan, dtype=float)
    assignments = np.full(len(labels), -1, dtype=int)
    summary = []
    for fold_number, (train_idx, test_idx) in enumerate(folds, start=1):
        estimator = estimator_factory()
        estimator.fit([matrix[i] for i in train_idx], labels[train_idx])
        fitted_classes = list(estimator.classes_)
        if fitted_classes != classes:
            raise CalibrationUnavailable(
                f"Fold {fold_number} trained classes differ from the full calibration classes")
        oof[test_idx, :] = estimator.predict_proba([matrix[i] for i in test_idx])
        assignments[test_idx] = fold_number
        summary.append({
            "fold": fold_number,
            "train_records": int(len(train_idx)),
            "heldout_records": int(len(test_idx)),
            "train_groups": int(len(set(groups[train_idx]))),
            "heldout_groups": int(len(set(groups[test_idx]))),
            "train_class_counts": dict(sorted(Counter(labels[train_idx]).items())),
            "heldout_class_counts": dict(sorted(Counter(labels[test_idx]).items()))})
    if np.isnan(oof).any() or (assignments < 1).any():
        raise CalibrationUnavailable("Out-of-fold predictions do not cover every calibration record")
    return classes, oof, assignments, summary


def _logit(values, epsilon=1e-6):
    values = np.clip(np.asarray(values, dtype=float), epsilon, 1.0 - epsilon)
    return np.log(values / (1.0 - values)).reshape(-1, 1)


def _fit_sigmoid(scores, labels):
    model = LogisticRegression(C=1e6, solver="lbfgs", max_iter=2000)
    model.fit(_logit(scores), np.asarray(labels, dtype=int))
    return model


def fit_binary_sigmoid(oof_probabilities, labels, classes, groups,
                       positive_label="adhesome_associated", minimum_groups_per_class=5):
    labels = np.asarray(labels, dtype=str)
    groups = np.asarray(groups, dtype=str)
    classes = list(classes)
    if positive_label not in classes or len(classes) != 2:
        raise CalibrationUnavailable("Binary probability calibration requires both membership classes")
    class_group_counts = {label: len(set(groups[labels == label])) for label in classes}
    if min(class_group_counts.values()) < int(minimum_groups_per_class):
        raise CalibrationUnavailable(
            "Insufficient independent sequence groups for membership calibration: "
            f"{class_group_counts}; require at least {minimum_groups_per_class} per class")
    positive_index = classes.index(positive_label)
    binary_truth = (labels == positive_label).astype(int)
    sigmoid = _fit_sigmoid(np.asarray(oof_probabilities)[:, positive_index], binary_truth)
    return {"method": "sigmoid_on_logit_oof_vote", "positive_label": positive_label,
            "classes": classes, "class_group_counts": class_group_counts,
            "minimum_groups_per_class": int(minimum_groups_per_class),
            "parameters": {"estimator": "sklearn.linear_model.LogisticRegression",
                           "C": 1e6, "solver": "lbfgs", "max_iter": 2000,
                           "input": "logit-transformed random-forest positive-class probability",
                           "target": "binary calibration label"},
            "estimator": sigmoid, "oof_records": int(len(labels))}


def fit_family_sigmoids(oof_probabilities, labels, classes, groups,
                        minimum_positive_groups=10, minimum_negative_groups=5,
                        minimum_eligible_classes=2):
    """Fit one-vs-rest sigmoid mappings for adequately represented family labels.

    Calibrated per-class estimates are renormalized over the eligible classes.
    Classes below the independent-group floor remain in raw diagnostic rankings,
    but are excluded from probability rankings and evidence-gated calls.
    """
    labels = np.asarray(labels, dtype=str)
    groups = np.asarray(groups, dtype=str)
    classes = list(classes)
    raw = np.asarray(oof_probabilities, dtype=float)
    models, support, unavailable = {}, {}, {}
    for class_index, label in enumerate(classes):
        positive = labels == label
        negative = ~positive
        positive_groups = len(set(groups[positive]))
        negative_groups = len(set(groups[negative]))
        support[label] = {"records": int(np.count_nonzero(positive)),
                          "independent_groups": positive_groups,
                          "negative_records": int(np.count_nonzero(negative)),
                          "negative_independent_groups": negative_groups}
        if positive_groups < int(minimum_positive_groups):
            unavailable[label] = (f"only {positive_groups} independent positive sequence groups; "
                                  f"requires {minimum_positive_groups}")
            continue
        if negative_groups < int(minimum_negative_groups):
            unavailable[label] = (f"only {negative_groups} independent negative sequence groups; "
                                  f"requires {minimum_negative_groups}")
            continue
        models[label] = _fit_sigmoid(raw[:, class_index], positive.astype(int))
    if len(models) < int(minimum_eligible_classes):
        reason=(f"only {len(models)} families met calibration support; "
                f"at least {minimum_eligible_classes} eligible families are required to form a non-degenerate distribution")
        for label in list(models):
            unavailable[label]=reason
        models={}
    return {"method": "one_vs_rest_sigmoid_on_logit_oof_vote_renormalized",
            "classes": classes, "eligible_classes": list(models), "estimators": models,
            "class_support": support, "unavailable_classes": unavailable,
            "minimum_positive_groups": int(minimum_positive_groups),
            "minimum_negative_groups": int(minimum_negative_groups),
            "minimum_eligible_classes": int(minimum_eligible_classes),
            "parameters": {"estimator": "sklearn.linear_model.LogisticRegression",
                           "C": 1e6, "solver": "lbfgs", "max_iter": 2000,
                           "input": "per-family logit-transformed random-forest probability",
                           "mapping": "one-vs-rest sigmoid, then normalize across eligible families"},
            "oof_records": int(len(labels))}


def apply_binary_sigmoid(calibration, raw_probabilities, classes):
    classes = list(classes)
    positive_index = classes.index(calibration["positive_label"])
    raw = np.asarray(raw_probabilities, dtype=float)
    return calibration["estimator"].predict_proba(_logit(raw[:, positive_index]))[:, 1]


def apply_family_sigmoids(calibration, raw_probabilities, classes):
    classes = list(classes)
    raw = np.asarray(raw_probabilities, dtype=float)
    eligible = calibration.get("eligible_classes", [])
    if not eligible:
        return eligible, np.zeros((len(raw), 0), dtype=float)
    calibrated = np.column_stack([
        calibration["estimators"][label].predict_proba(_logit(raw[:, classes.index(label)]))[:, 1]
        for label in eligible])
    row_sums = calibrated.sum(axis=1, keepdims=True)
    calibrated = np.divide(calibrated, row_sums, out=np.zeros_like(calibrated), where=row_sums > 0)
    return eligible, calibrated


def reliability_bins(truth, probabilities, n_bins=10):
    truth = np.asarray(truth, dtype=float)
    probabilities = np.clip(np.asarray(probabilities, dtype=float), 0.0, 1.0)
    edges = np.linspace(0.0, 1.0, int(n_bins) + 1)
    bins = []
    weighted_error = 0.0
    total = len(truth)
    for index in range(int(n_bins)):
        include = ((probabilities >= edges[index]) &
                   ((probabilities <= edges[index + 1]) if index == n_bins - 1
                    else (probabilities < edges[index + 1])))
        count = int(np.count_nonzero(include))
        mean_probability = float(probabilities[include].mean()) if count else None
        observed_fraction = float(truth[include].mean()) if count else None
        if count:
            weighted_error += (count / total) * abs(mean_probability - observed_fraction)
        bins.append({"lower": float(edges[index]), "upper": float(edges[index + 1]),
                     "count": count, "mean_probability": mean_probability,
                     "observed_fraction": observed_fraction})
    return {"bins": bins, "expected_calibration_error": weighted_error if total else None,
            "records": total}


def probability_metrics_binary(truth, probabilities, raw_probabilities=None, n_bins=10):
    truth = np.asarray(truth, dtype=int)
    probabilities = np.clip(np.asarray(probabilities, dtype=float), 1e-12, 1.0 - 1e-12)
    if not len(truth):
        return {"records": 0, "brier_score": None, "log_loss": None,
                "reliability": reliability_bins([], [], n_bins)}
    from sklearn.metrics import log_loss
    result = {"records": int(len(truth)),
              "brier_score": float(np.mean((probabilities - truth) ** 2)),
              "log_loss": float(log_loss(truth, np.column_stack((1.0 - probabilities, probabilities)),
                                          labels=[0, 1])),
              "reliability": reliability_bins(truth, probabilities, n_bins)}
    if raw_probabilities is not None:
        raw = np.clip(np.asarray(raw_probabilities, dtype=float), 1e-12, 1.0 - 1e-12)
        result["raw"] = {"brier_score": float(np.mean((raw - truth) ** 2)),
                         "log_loss": float(log_loss(truth, np.column_stack((1.0 - raw, raw)),
                                                     labels=[0, 1])),
                         "reliability": reliability_bins(truth, raw, n_bins)}
    return result


def probability_metrics_multiclass(truth, probabilities, classes, raw_probabilities=None, n_bins=10):
    classes = list(classes)
    truth = np.asarray(truth, dtype=str)
    probabilities = np.asarray(probabilities, dtype=float)
    if not len(truth):
        return {"records": 0, "classes": classes, "brier_score": None, "log_loss": None,
                "top_label_reliability": reliability_bins([], [], n_bins),
                "per_class_reliability": {}}
    if probabilities.shape != (len(truth), len(classes)):
        raise ValueError("Multiclass probability matrix does not match its labels/classes")
    if not np.allclose(probabilities.sum(axis=1), 1.0, atol=1e-6):
        raise ValueError("Calibrated family probabilities must sum to one")
    if set(truth) - set(classes):
        raise ValueError("Multiclass probability cohort contains a family outside calibrated classes")
    from sklearn.metrics import log_loss
    encoded = np.column_stack([(truth == label).astype(float) for label in classes])
    top_index = probabilities.argmax(axis=1)
    result = {"records": int(len(truth)), "classes": classes,
              "brier_score": float(np.mean(np.sum((probabilities - encoded) ** 2, axis=1))),
              "log_loss": float(log_loss(truth, probabilities, labels=classes)),
              "top_label_reliability": reliability_bins(
                  (truth == np.asarray(classes, dtype=object)[top_index]).astype(int),
                  probabilities[np.arange(len(truth)), top_index], n_bins),
              "per_class_reliability": {
                  label: reliability_bins(encoded[:, index], probabilities[:, index], n_bins)
                  for index, label in enumerate(classes)}}
    if raw_probabilities is not None:
        raw = np.asarray(raw_probabilities, dtype=float)
        if raw.shape != probabilities.shape:
            raise ValueError("Raw and calibrated family scores do not share the same scored cohort")
        raw_top = raw.argmax(axis=1)
        raw_top_classes = np.asarray(classes, dtype=object)[raw_top]
        result["raw"] = {"brier_score": float(np.mean(np.sum((raw - encoded) ** 2, axis=1))),
                         "log_loss": float(log_loss(truth, raw, labels=classes)),
                         "top_label_reliability": reliability_bins(
                             (truth == raw_top_classes).astype(int),
                             raw[np.arange(len(truth)), raw_top], n_bins)}
    return result
