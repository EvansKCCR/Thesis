"""Paired expansion regression safeguards; no biological analysis is launched."""
import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import evaluate_family_expansion as expansion
import evaluate_validation as evaluation


class ExpansionTests(unittest.TestCase):
    def pair(self):
        families = {'pipeline_version': '1.0', 'families': {'parent': {}}}
        identity = dict(pipeline_version='1.0', pfam_cdd_aliases={}, additional_cdd_domains={}, rules=[
            dict(id='test', parents=['parent'], family='child', protein_name='Child-family protein',
                 required_groups=[dict(models=['PF00017'], count=1)], specific_cdd_pattern='^SH2_child$',
                 ordered_groups=[0], allow_diagnostic_dispatch=True)])
        return families, identity

    def test_supplied_pair_uses_only_supported_fields_and_routes(self):
        expansion.validate_pair(*(json.loads((expansion.AUDIT / n).read_text()) for n in
            ('family_rules_expanded.json', 'family_identity_rules_expanded.json')))

    def test_unknown_gate_cannot_be_silently_ignored(self):
        f, i = self.pair(); i['rules'][0]['forbidden_cdd_pattern'] = 'competing_family'
        with self.assertRaisesRegex(ValueError, 'Unsupported/missing'):
            expansion.validate_pair(f, i)

    def test_bad_route_order_duplicate_and_empty_core_fail(self):
        for mutation in (lambda r: r.update(parents=['missing']),
                         lambda r: r.update(ordered_groups=[1]),
                         lambda r: r.update(required_groups=[])):
            f, i = self.pair(); mutation(i['rules'][0])
            with self.assertRaises(ValueError):
                expansion.validate_pair(f, i)
        f, i = self.pair(); i['rules'] *= 2
        with self.assertRaisesRegex(ValueError, 'Duplicate'):
            expansion.validate_pair(f, i)

    def test_unsupported_protocol_version_rejected(self):
        f, i = self.pair(); i['pipeline_version'] = '2.0'
        with self.assertRaisesRegex(ValueError, 'version 1.0'):
            expansion.validate_pair(f, i)

    def test_source_hash_mismatch_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / 'source.tsv'; path.write_text('source\n')
            with self.assertRaisesRegex(ValueError, 'SHA256'):
                expansion.read_verified(path, 'changed', {})

    def context(self, folder):
        raw = Path(folder) / 'validation.faa'
        raw.write_text('>candidate|POS|POS\nMALL\n>candidate|NEG|NEG\nMFLL\n>candidate|MISS|MISS\nMPLL\n')
        mapping = {f'candidate|{a}|{a}': dict(protein_id='validation__'+a, decision='retained') for a in ('POS', 'NEG', 'MISS')}
        sequences = {'validation__POS': 'MALL', 'validation__NEG': 'MFLL', 'validation__MISS': 'MPLL'}
        initial = {'validation__POS': {'parent'}, 'validation__NEG': {'parent'}}
        gold = {a: dict(control_class='positive' if a != 'NEG' else 'hard_negative', expected_families={'child'},
                        expected_broad_families=set(), expected_specific_families={'child'}) for a in ('POS', 'NEG', 'MISS')}
        trusted = {a: evaluation.sequence_hash(s) for a, s in [('POS', 'MALL'), ('NEG', 'MFLL'), ('MISS', 'MPLL')]}
        return raw, mapping, sequences, initial, gold, trusted, set(), set(), set()

    def call(self):
        return dict(assigned_family='child', domain_score_0_4=4, grade='A', family_confidence='strong',
                    family_support_level='specific_supported', annotation_resolution='family_lineage')

    def test_full_cohort_and_exclusion_ablation_are_independent(self):
        with tempfile.TemporaryDirectory() as folder:
            context = self.context(folder)
            before = {'validation__POS': self.call(), 'validation__NEG': self.call()}
            final = {'validation__POS': self.call()}
            excluded = {'validation__NEG': dict(exclusion_gate='exact_negative_control_sequence',
                                              exclusion_gates='exact_negative_control_sequence')}
            result = expansion.predictions(*context, before, final, excluded, {})
            by_id = {r['accession']: r for r in result}
            self.assertEqual(len(result), 3)
            self.assertFalse(by_id['MISS']['stage1_recruitment'])
            self.assertFalse(by_id['MISS']['dual_gate_supported'])
            self.assertTrue(by_id['NEG']['audit_before_negative_exclusion'])
            self.assertTrue(by_id['NEG']['audit_without_exact_negative_exclusion'])
            self.assertFalse(by_id['NEG']['dual_gate_supported'])
            excluded['validation__NEG']['exclusion_gates'] += ';diagnostic_negative_subfamily_and_complete_architecture'
            result = expansion.predictions(*context, before, final, excluded, {})
            self.assertFalse(next(r for r in result if r['accession'] == 'NEG')['audit_without_exact_negative_exclusion'])

    def test_expected_label_does_not_make_a_wrong_supported_family_correct(self):
        with tempfile.TemporaryDirectory() as folder:
            context = self.context(folder)
            calls = {'validation__POS': self.call()}
            baseline = expansion.predictions(*context, calls, calls, {}, {'POS': 'identity_proposal'})
            wrong = copy.deepcopy(calls); wrong['validation__POS']['assigned_family'] = 'wrong'
            expanded = expansion.predictions(*context, wrong, wrong, {}, {'POS': 'identity_proposal'})
            row = next(r for r in expansion.compare(baseline, expanded) if r['accession'] == 'POS')
            self.assertTrue(row['baseline_correct_supported'])
            self.assertFalse(row['expanded_correct_supported'])
            self.assertTrue(row['expanded_supported'])


if __name__ == '__main__':
    unittest.main()
