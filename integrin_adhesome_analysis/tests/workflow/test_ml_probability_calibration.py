"""Focused synthetic tests for v5.8 probability calibration and leakage control."""
import sys
import json
import tempfile
from pathlib import Path
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np

ML = Path(__file__).resolve().parents[2] / 'resource_library/adhesome_ml'
sys.path.insert(0, str(ML))
try:
    import joblib
    import sklearn
except ImportError:
    raise unittest.SkipTest('ML optional dependencies are not installed')

import probability_calibration as calibration
import family_gate


class ProbabilityCalibrationTests(unittest.TestCase):
    def test_group_folds_keep_duplicate_sequences_together_and_cover_records(self):
        labels = np.array(['positive'] * 20 + ['negative'] * 20)
        groups = np.array([f'p{i // 2}' for i in range(20)] +
                          [f'n{i // 2}' for i in range(20)])
        folds = calibration.grouped_splits(labels, groups, max_splits=5, random_state=7)
        held_out = []
        for train, test in folds:
            self.assertFalse(set(groups[train]) & set(groups[test]))
            self.assertEqual(set(labels[train]), {'positive', 'negative'})
            self.assertEqual(set(labels[test]), {'positive', 'negative'})
            held_out.extend(test.tolist())
        self.assertEqual(sorted(held_out), list(range(len(labels))))

    def test_missing_cdhit_falls_back_to_exact_sequence_hash_groups(self):
        groups, metadata = calibration.sequence_cluster_groups(
            ['MKWVTF', 'MKWVTF', 'GHHHHH'], executable='definitely-not-cd-hit')
        self.assertEqual(groups[0], groups[1])
        self.assertNotEqual(groups[0], groups[2])
        self.assertEqual(metadata['method'], 'exact_sequence_sha256_fallback')
        self.assertIn('not found', metadata['fallback_reason'])

    def test_cdhit_cluster_report_defines_fold_groups(self):
        def fake_cdhit(command, **_kwargs):
            output = Path(command[command.index('-o') + 1])
            Path(str(output) + '.clstr').write_text(
                '>Cluster 0\n0 6aa, >seq_00000000... *\n1 6aa, >seq_00000001... at 99%\n'
                '>Cluster 1\n0 6aa, >seq_00000002... *\n', encoding='utf-8')
            return SimpleNamespace(returncode=0, stderr='', stdout='')

        with patch.object(calibration.shutil, 'which', return_value='fake-cd-hit'), \
             patch.object(calibration.subprocess, 'run', side_effect=fake_cdhit):
            groups, metadata = calibration.sequence_cluster_groups(
                ['MKWVTF', 'MKWVTA', 'GHHHHH'])
        self.assertEqual(groups[0], groups[1])
        self.assertNotEqual(groups[0], groups[2])
        self.assertEqual(metadata['method'], 'cdhit_grouped_out_of_fold')

    def test_oof_estimator_is_fit_without_held_out_sequence_groups(self):
        from sklearn.ensemble import RandomForestClassifier
        from sklearn.feature_extraction import DictVectorizer
        from sklearn.pipeline import make_pipeline

        labels = np.array(['positive'] * 20 + ['negative'] * 20)
        groups = np.array([f'p{i // 2}' for i in range(20)] +
                          [f'n{i // 2}' for i in range(20)])
        matrix = [{'signal': float(label == 'positive'), 'noise': float(i % 3)}
                  for i, label in enumerate(labels)]
        classes, probabilities, assignment, summary = calibration.cross_fitted_probabilities(
            lambda: make_pipeline(DictVectorizer(), RandomForestClassifier(
                n_estimators=8, random_state=11, n_jobs=1)), matrix, labels, groups,
            max_splits=4, random_state=7)
        self.assertEqual(set(classes), {'positive', 'negative'})
        self.assertEqual(probabilities.shape, (40, 2))
        self.assertTrue(np.all((probabilities >= 0) & (probabilities <= 1)))
        self.assertEqual(set(assignment), {1, 2, 3, 4})
        self.assertEqual(sum(item['heldout_records'] for item in summary), 40)

    def test_membership_and_family_calibrators_return_valid_probabilities(self):
        rng = np.random.default_rng(42)
        membership_labels = np.array(['adhesome_associated'] * 30 + ['non_adhesome'] * 30)
        membership_groups = np.array([f'm{i}' for i in range(60)])
        membership_raw = np.column_stack((np.full(60, .5), np.full(60, .5)))
        membership_raw[:30, 1] = rng.uniform(.55, .95, 30)
        membership_raw[:30, 0] = 1 - membership_raw[:30, 1]
        membership_raw[30:, 1] = rng.uniform(.05, .45, 30)
        membership_raw[30:, 0] = 1 - membership_raw[30:, 1]
        membership = calibration.fit_binary_sigmoid(
            membership_raw, membership_labels,
            ['non_adhesome', 'adhesome_associated'], membership_groups,
            minimum_groups_per_class=5)
        calibrated_membership = calibration.apply_binary_sigmoid(
            membership, membership_raw, ['non_adhesome', 'adhesome_associated'])
        self.assertTrue(np.all((calibrated_membership >= 0) & (calibrated_membership <= 1)))

        family_labels = np.array(['A'] * 20 + ['B'] * 20 + ['C'] * 20)
        family_groups = np.array([f'f{i}' for i in range(60)])
        family_raw = rng.uniform(.05, .2, (60, 3))
        for index, label in enumerate(family_labels):
            family_raw[index, ['A', 'B', 'C'].index(label)] = rng.uniform(.6, .9)
        family_raw /= family_raw.sum(axis=1, keepdims=True)
        family = calibration.fit_family_sigmoids(
            family_raw, family_labels, ['A', 'B', 'C'], family_groups,
            minimum_positive_groups=10, minimum_negative_groups=5,
            minimum_eligible_classes=2)
        names, calibrated_family = calibration.apply_family_sigmoids(
            family, family_raw, ['A', 'B', 'C'])
        self.assertEqual(names, ['A', 'B', 'C'])
        self.assertTrue(np.all((calibrated_family >= 0) & (calibrated_family <= 1)))
        np.testing.assert_allclose(calibrated_family.sum(axis=1), np.ones(60), atol=1e-7)

    def test_insufficient_family_support_is_unavailable_not_raw_votes(self):
        labels = np.array(['common'] * 20 + ['rare'] * 4)
        groups = np.array([f'c{i}' for i in range(20)] + [f'r{i}' for i in range(4)])
        raw = np.tile([.7, .3], (24, 1))
        result = calibration.fit_family_sigmoids(
            raw, labels, ['common', 'rare'], groups, minimum_positive_groups=10,
            minimum_negative_groups=5, minimum_eligible_classes=2)
        self.assertEqual(result['eligible_classes'], [])
        self.assertIn('rare', result['unavailable_classes'])
        names, values = calibration.apply_family_sigmoids(result, raw, ['common', 'rare'])
        self.assertEqual(names, [])
        self.assertEqual(values.shape, (24, 0))

    def test_membership_calibration_refuses_too_few_independent_groups(self):
        raw = np.tile([.4, .6], (8, 1))
        labels = np.array(['adhesome_associated'] * 4 + ['non_adhesome'] * 4)
        groups = np.array(['p1', 'p1', 'p2', 'p2', 'n1', 'n1', 'n2', 'n2'])
        with self.assertRaisesRegex(calibration.CalibrationUnavailable, 'Insufficient independent'):
            calibration.fit_binary_sigmoid(raw, labels,
                ['non_adhesome', 'adhesome_associated'], groups,
                minimum_groups_per_class=5)

    def test_probability_metrics_report_scores_and_reliability_counts(self):
        binary = calibration.probability_metrics_binary(
            [1, 0, 1, 0], [.8, .2, .6, .4], [.7, .3, .7, .3], n_bins=4)
        self.assertIn('brier_score', binary)
        self.assertIn('log_loss', binary)
        self.assertEqual(sum(x['count'] for x in binary['reliability']['bins']), 4)
        family = calibration.probability_metrics_multiclass(
            ['A', 'B', 'A'], [[.8, .2], [.1, .9], [.6, .4]], ['A', 'B'],
            raw_probabilities=[[.7, .3], [.2, .8], [.5, .5]], n_bins=3)
        self.assertIn('brier_score', family)
        self.assertIn('per_class_reliability', family)
        self.assertEqual(sum(x['count'] for x in
                             family['per_class_reliability']['A']['bins']), 3)
        with self.assertRaisesRegex(ValueError, 'sum to one'):
            calibration.probability_metrics_multiclass(
                ['A'], [[.8, .8]], ['A', 'B'])

    def test_calibrator_save_load_and_old_model_rejection(self):
        raw = np.tile([.3, .7], (20, 1))
        labels = np.array(['adhesome_associated'] * 10 + ['non_adhesome'] * 10)
        groups = np.array([f'g{i}' for i in range(20)])
        membership = calibration.fit_binary_sigmoid(raw, labels,
            ['non_adhesome', 'adhesome_associated'], groups, minimum_groups_per_class=5)
        saved = {'version': '5.8',
                 'feature_schema': family_gate.SCHEMA,
                 'probability_calibration': {
                     'membership': membership,
                     'family': {'eligible_classes': [], 'estimators': {}}}}
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'model.joblib'
            joblib.dump(saved, path)
            loaded = joblib.load(path)
            family_gate._require_calibrated_model(loaded, 'Synthetic test')
            output = calibration.apply_binary_sigmoid(
                loaded['probability_calibration']['membership'], raw,
                ['non_adhesome', 'adhesome_associated'])
            self.assertEqual(output.shape, (20,))
            self.assertEqual(loaded['version'], '5.8')
        with self.assertRaisesRegex(ValueError, 'must be retrained'):
            family_gate._require_calibrated_model(
                {'version': '5.7', 'feature_schema': 'sequence_family_gate_v5_7'}, 'Evaluation')

    def test_parent_family_probability_records_unrepresented_descendants(self):
        result = family_gate._family_probability_assignment_details(
            'parent', {'parent': .2, 'A': .7, 'B': .3},
            {'A': 'parent', 'B': 'parent', 'C': 'parent'}, ['parent', 'A', 'B', 'C'])
        self.assertEqual(result['basis'], 'sum_of_eligible_calibrated_descendant_classes')
        self.assertAlmostEqual(result['probability'], 1.0)
        self.assertEqual(result['eligible_descendant_classes'], ['A', 'B'])
        self.assertEqual(result['unrepresented_descendant_classes'], ['C'])

    def test_validation_report_writes_raw_and_calibrated_comparisons(self):
        membership_metrics = calibration.probability_metrics_binary(
            [1, 0], [.8, .2], [.7, .3], n_bins=2)
        family_metrics = calibration.probability_metrics_multiclass(
            ['A', 'B'], [[.8, .2], [.2, .8]], ['A', 'B'],
            raw_probabilities=[[.7, .3], [.3, .7]], n_bins=2)
        report = {
            'note': 'Validation is evaluation only.',
            'class_counts': {'specific': 2},
            'membership': {
                'records': 2, 'positive_precision': 1.0, 'positive_recall': 1.0,
                'positive_f1': 1.0,
                'precalibration_policy_comparison': {
                    'positive_precision': .5, 'positive_recall': .5, 'positive_f1': .5,
                    'confusion_matrix': [[1, 0], [1, 0]]},
                'probability_calibration': membership_metrics},
            'family': {
                'labelled_positives': 2, 'unmapped_positives': 0,
                'correct_accepted': 2, 'accepted': 2, 'strict_recall': 1.0,
                'top_candidate_correct': 2, 'top_three_contains_truth': 2,
                'precalibration_policy_comparison': {
                    'macro_precision': .5, 'macro_recall': .5, 'macro_f1': .5,
                    'confusion_matrix': [[1, 0], [1, 0]]},
                'probability_calibration': {
                    'family_labelled_records': 2, 'scored_records': 2,
                    'unavailable_true_family_records': 0, 'coverage': 1.0,
                    'calibration_availability_by_family': {}, 'metrics': family_metrics},
                'per_family': {
                    'A': {'records': 1, 'correct': 1, 'precalibration_correct': 0},
                    'B': {'records': 1, 'correct': 1, 'precalibration_correct': 0}}}}
        with tempfile.TemporaryDirectory() as temp:
            family_gate.write_validation_report(Path(temp), report)
            html = (Path(temp) / 'internal_validation_report.html').read_text(encoding='utf-8')
            report_json = json.loads((Path(temp) / 'internal_validation_report.json').read_text(encoding='utf-8'))
        self.assertIn('Previous versus calibrated-policy calls', html)
        self.assertIn('Membership probability calibration', html)
        self.assertIn('Family probability calibration', html)
        self.assertEqual(report_json['membership']['records'], 2)


if __name__ == '__main__':
    unittest.main()
