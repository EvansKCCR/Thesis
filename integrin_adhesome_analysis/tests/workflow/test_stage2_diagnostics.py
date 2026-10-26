"""Synthetic saved-output tests; no biological inference or calibration fitting."""

# Resolve the organized bundle independently of the caller working directory.
import sys as _bundle_sys
from pathlib import Path as _BundlePath
_BUNDLE_ROOT=next(p for p in _BundlePath(__file__).resolve().parents if (p/'pipeline_paths.py').is_file())
if str(_BUNDLE_ROOT) not in _bundle_sys.path:_bundle_sys.path.insert(0,str(_BUNDLE_ROOT))
from pipeline_paths import activate as _activate_bundle_paths
_activate_bundle_paths()
import importlib
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
import evaluate_validation as ev
import validation_reporting as vr


class Stage2DiagnosticTests(unittest.TestCase):
    def candidate(self,**changes):
        row=dict(putative_audited_family='coronin_like',grade='B',domain_score_0_4=4,
            assignment_confidence='moderate',annotation_resolution='family',family_support_level='specific_supported')
        row.update(changes);return row

    def prediction(self,key='P1',**changes):
        row=dict(accession=key,protein_id=key,sequence_sha256=key,evaluation_status='scored',
            label_sequence_verification='trusted_fasta_match',control_class='positive',expected_families='coronin_like',
            expected_specific_families='coronin_like',expected_broad_families='',family_rule_coverage='configured',
            identity_calibration_training_overlap=False,exact_forward_reference_overlap=False,exact_reverse_reference_overlap=False,
            stage1_recruitment=True,stage2_candidate_retention=True,audit_before_negative_exclusion=True,
            audit_without_exact_negative_exclusion=True,dual_gate_supported=True,dual_gate_strong=False,
            broad_family_supported=True,specific_family_supported=True,provisional_annotation=False,
            ambiguous_annotation=False,audit_including_provisional=True,audited_family='coronin_like',
            pre_exclusion_family='coronin_like',negative_exclusion=False,exclusion_gate='',preparation_decision='retained')
        row.update(changes);return row

    def test_supported_resolution_and_legacy_unknown_specificity(self):
        self.assertEqual(ev.resolution_decisions(self.candidate())['specific_family_supported'],True)
        broad=ev.resolution_decisions(self.candidate(family_support_level='broad_supported',annotation_resolution='broad_family'))
        self.assertTrue(broad['broad_family_supported']);self.assertFalse(broad['specific_family_supported'])
        legacy=ev.resolution_decisions(self.candidate(family_support_level='',annotation_resolution=''))
        self.assertTrue(legacy['broad_family_supported']);self.assertIsNone(legacy['specific_family_supported'])
        self.assertFalse(ev.resolution_decisions({})['specific_family_supported'])

    def test_resolution_labels_cannot_promote_grade_c_or_ambiguity(self):
        for changes in [dict(grade='C',domain_score_0_4=2),dict(assignment_confidence='ambiguous'),dict(domain_score_0_4=2)]:
            with self.subTest(changes=changes):
                result=ev.resolution_decisions(self.candidate(**changes))
                self.assertFalse(result['broad_family_supported']);self.assertFalse(result['specific_family_supported'])
        result=ev.resolution_decisions(self.candidate(grade='C',family_support_level='provisional',assignment_confidence='provisional'))
        self.assertTrue(result['provisional_annotation']);self.assertFalse(result['broad_family_supported'])

    def test_macro_recall_weights_rare_families_and_keeps_unknown_wrong_names(self):
        rows=[self.prediction('A'+str(i)) for i in range(10)]
        rows.append(self.prediction('B',expected_families='ERM_like',expected_specific_families='ERM_like',
            dual_gate_supported=False,audited_family='unassigned'))
        result,groups=vr.macro_family_metrics(rows,'dual_gate_supported')
        self.assertEqual(result['macro_recall'],.5);self.assertEqual(result['macro_F1'],.5)
        self.assertEqual(result['macro_precision_defined_groups'],1)
        rows.append(self.prediction('N',control_class='hard_negative',audited_family='unknown_wrong_family'))
        result,groups=vr.macro_family_metrics(rows,'dual_gate_supported')
        self.assertEqual(result['macro_precision'],.5)
        self.assertTrue(any(r['group_kind']=='predicted_only' and r['FP']==1 for r in groups))

    def test_broad_ancestor_credit_is_separate_from_specific_and_legacy_gold(self):
        row=self.prediction(expected_families='arp2_like',expected_specific_families='arp2_like',expected_broad_families='actin',
            audited_family='actin',family_support_level='broad_supported',specific_family_supported=False)
        _,_,families,_=ev.aggregate([row]);by_policy={r['policy']:r for r in families}
        self.assertEqual(by_policy['broad_family_supported']['TP'],1)
        self.assertEqual(by_policy['specific_family_supported']['FN'],1)
        self.assertEqual(by_policy['dual_gate_supported']['TP'],0)
        self.assertNotIn('stage2_candidate_retention',by_policy)
        self.assertNotIn('ambiguous_annotation',by_policy)

    def test_gate_codes_require_saved_trace_or_observed_disposition(self):
        row=self.prediction(identity_gate_trace=json.dumps({'rules':[{'rule':'x','reason_codes':['MISSING_REQUIRED_DOMAIN','DOMAIN_ORDER_CONFLICT']}]}))
        reason=vr.gate_diagnostic(row)
        self.assertEqual(reason['reason_codes'],'MISSING_REQUIRED_DOMAIN;DOMAIN_ORDER_CONFLICT')
        self.assertEqual(reason['gate_trace_status'],'available')
        row=self.prediction(review_architecture_evidence='family-specific CDD absent; no alias recognition',identity_gate_trace='')
        self.assertNotIn('CDD_ALIAS_UNRECOGNIZED',vr.gate_diagnostic(row)['reason_codes'])
        row['identity_gate_trace']=json.dumps({'limitation_codes':['CDD_ALIAS_UNRECOGNIZED']})
        self.assertEqual(vr.gate_diagnostic(row)['reason_codes'],'CDD_ALIAS_UNRECOGNIZED')

    def test_missing_stage1_exclusion_and_ambiguity_are_observed_codes(self):
        row=self.prediction(stage1_recruitment=False,negative_exclusion=True,assignment_confidence='ambiguous')
        result=vr.gate_diagnostic(row)
        self.assertEqual(result['reason_codes'],'NO_STAGE1_HYPOTHESIS;COMPETING_FAMILY;NEGATIVE_EXCLUSION')
        self.assertEqual(result['gate_trace_status'],'unavailable')
        row['identity_gate_trace']='bad JSON'
        self.assertEqual(vr.gate_diagnostic(row)['gate_trace_status'],'invalid_json')

    def test_unmapped_family_stays_in_membership_denominator_and_coverage_is_explicit(self):
        rows=[self.prediction(),self.prediction('P2',expected_families='',expected_specific_families='',family_rule_coverage='needs_mapping_review')]
        metrics,_,family,_=ev.aggregate(rows);reports=vr.reports(rows,metrics)
        coverage=reports['family_label_coverage_summary.tsv'][0]
        self.assertEqual((coverage['positive_controls'],coverage['mapped_expected_family_controls'],coverage['unmapped_expected_family_controls']),(2,1,1))
        self.assertEqual(coverage['mapped_fraction'],.5)
        main=next(r for r in metrics if r['policy']=='dual_gate_supported' and r['population']=='all_controls')
        self.assertEqual(main['n'],2)
        self.assertEqual(next(r for r in family if r['policy']=='dual_gate_supported')['positive_controls_without_mapped_family'],1)

    def test_saved_output_integration_exports_trace_and_unknown_legacy_resolution(self):
        fixture=importlib.import_module('tests.workflow.test_validation_evaluation').ValidationEvaluationTests()
        with tempfile.TemporaryDirectory() as td:
            args=fixture.fixture(Path(td));path=args['audit_dir']/'audited_family_assignments.tsv'
            rows=ev.read_table(path)
            for row in rows:
                row.update(identity_gate_trace='',family_support_level='',annotation_resolution='')
                if row['protein_id']=='validation__P1':
                    row.update(family_support_level='specific_supported',annotation_resolution='family',
                        identity_gate_trace=json.dumps({'rules':[{'reason_codes':['ORTHOLOGY_UNRESOLVED']}]}))
            fixture.write(path,rows);original=ev.sha(path)
            ev.evaluate(**args)
            predictions={r['accession']:r for r in ev.read_table(args['out']/'control_predictions.tsv')}
            self.assertEqual(predictions['P1']['specific_family_supported'],'True')
            self.assertEqual(predictions['P2']['specific_family_supported'],'')
            self.assertEqual(predictions['P4']['provisional_annotation'],'True')
            diagnostics={r['accession']:r for r in ev.read_table(args['out']/'stage2_gate_diagnostics.tsv')}
            self.assertEqual(diagnostics['P1']['reason_codes'],'ORTHOLOGY_UNRESOLVED')
            self.assertEqual(diagnostics['N1']['reason_codes'],'NEGATIVE_EXCLUSION')
            self.assertEqual(diagnostics['P3']['reason_codes'],'NO_STAGE1_HYPOTHESIS')
            self.assertTrue((args['out']/'family_macro_metrics.tsv').exists())
            self.assertTrue((args['out']/'family_label_coverage_summary.tsv').exists())
            self.assertEqual(ev.sha(path),original)

    def test_sequence_deduplication_does_not_mix_unknown_and_known_specificity(self):
        rows=[self.prediction('A',sequence_sha256='same',specific_family_supported=None),
              self.prediction('B',sequence_sha256='same',specific_family_supported=False)]
        populations=dict(vr.populations(rows))
        self.assertEqual(populations['unique_sequence_consistent_controls'],[])

    def test_broad_only_gold_does_not_establish_specific_family_correctness(self):
        from openpyxl import Workbook
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);labels=root/'labels.xlsx';fasta=root/'control.faa'
            book=Workbook();sheet=book.active;sheet.title='control_identity'
            sheet.append(['class','accession','genes','protein_name','trusted_family'])
            sheet.append(['Positive control','BROAD','','Actin-family protein','actin']);book.save(labels)
            fasta.write_text('>BROAD\nMAAA\n')
            gold,*_=ev.load_gold(labels,fasta)
            self.assertEqual(gold['BROAD']['expected_families'],{'actin'})
            self.assertEqual(gold['BROAD']['expected_specific_families'],set())


if __name__=='__main__':unittest.main()
