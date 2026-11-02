"""Synthetic checks of heuristic context scaling and paired-ablation accounting."""

# Resolve source resources independently of the caller working directory.
import sys as _layout_sys
from pathlib import Path as _LayoutPath
_LAYOUT_BUNDLE_ROOT = next(p for p in _LayoutPath(__file__).resolve().parents
                          if (p/'config/layout_paths.json').is_file() and (p/'pipeline_paths.py').is_file())
if str(_LAYOUT_BUNDLE_ROOT) not in _layout_sys.path:
    _layout_sys.path.insert(0, str(_LAYOUT_BUNDLE_ROOT))
from pipeline_paths import activate as _layout_activate, resource_path as _layout_resource
_layout_activate()


import argparse
import importlib.util
import json
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

ROOT=_LAYOUT_BUNDLE_ROOT;sys.path.insert(0,str(ROOT))
import adhesome_workflow as workflow
import identity_probability
from orthogroup_context import validate_weight,configuration
from maintenance.evaluate_orthogroup_ablation import invariants,changes_summary

spec=importlib.util.spec_from_file_location('context_test_audit',_layout_resource(ROOT/'annotation/integrated_adhesome_pipeline.py'))
audit=importlib.util.module_from_spec(spec);spec.loader.exec_module(audit)


class OrthogroupContextTests(unittest.TestCase):
    def setUp(self):
        self.rules=json.loads((_layout_resource(ROOT/'detection/config/family_rules.json')).read_text())
        self.row=dict(protein_id='Target__p',Specie='Target',family='profilin',length=150,
            architecture='PF00235:10-130',Orthology_HOG_reference_species_count=3,Orthology_HOG_status='shared HOG',
            MotifPrediction_context_supported_count=1,**{'DeepLoc_2.1':'Cytoplasm'})

    def call(self,weight=None,row=None):
        rules=dict(self.rules,_orthogroup_context_weight=weight)
        return audit.audit_record(row or self.row,rules)

    def test_0541_changes_only_context_score_and_is_not_applied_twice(self):
        row=dict(self.row,Orthology_benchmark_status='validated',Orthology_benchmark_weight=.5)
        base=self.call(row=row);weighted=self.call(.541,row)
        self.assertEqual(base['orthology_score_0_3'],1.5)
        self.assertEqual(weighted['orthology_score_0_3'],1.623)
        self.assertEqual(weighted['total_score_0_10'],8.623)
        for field in ('domain_score_0_4','motif_score_0_2','localization_score_0_1','audit_family','identity_rule','architecture_evidence'):
            self.assertEqual(base.get(field),weighted.get(field))
        self.assertEqual(weighted['orthology_benchmark_weight_used'],'')
        self.assertEqual(weighted['orthogroup_context_policy'],'heuristic')
        self.assertFalse(configuration(.541)['probability_interpretation'])

    def test_default_baseline_and_rejected_benchmark_remain_reproducible(self):
        baseline=audit.audit_record(self.row,self.rules)
        explicit_none=self.call()
        self.assertEqual(baseline,explicit_none)
        rejected=dict(self.row,Orthology_benchmark_status='quality_rejected',Orthology_benchmark_weight='')
        self.assertEqual(self.call(row=rejected)['orthology_score_0_3'],3)
        self.assertEqual(self.call(.541,rejected)['orthology_score_0_3'],1.623)

    def test_seven_species_and_all_sensitivity_coefficients_use_the_same_scaling(self):
        for sp in ('Sbovis','Scura','Shae','Sinter','Sjap','Sman','Srhoda'):
            for weight in (.25,.4,.541,.7,1):
                row=self.call(weight,dict(self.row,Specie=sp,protein_id=sp+'__p'))
                self.assertEqual(row['orthology_score_0_3'],round(3*weight,6))
                self.assertEqual(row['orthogroup_context_weight'],weight)
                self.assertEqual(row['grade'],'A')

    def test_zero_is_an_ablation_but_complete_domain_identity_is_retained(self):
        result=self.call(0)
        self.assertEqual(result['orthology_score_0_3'],0)
        self.assertEqual(result['grade'],'B')
        self.assertEqual(result['assigned_family'],'profilin')
        self.assertTrue(result['domain_gate_passed'])
        self.assertEqual(result['orthogroup_context_policy'],'disabled')

    def test_missing_domains_cannot_be_rescued_by_context_or_direct_reference_claim(self):
        row=dict(self.row,architecture='',reference_subfamily='profilin',subfamily_support_method='curated_reference')
        self.assertEqual(self.call(1,row)['assigned_family'],'unassigned')
        partial=dict(self.row,family='alpha_actinin',length=900,architecture='PF00307:10-110')
        self.assertEqual(self.call(.541,partial)['grade'],'D')
        legacy=dict(self.rules,_stage2_profile='original',_orthogroup_context_weight=.541)
        self.assertEqual(audit.audit_record(partial,legacy)['grade'],'C')
        self.assertEqual(self.call(0,partial)['grade'],'D')

    def test_invalid_weights_and_probability_misinterpretation_are_rejected(self):
        for value in (True,-1,1.01,float('nan'),float('inf'),'0.541'):
            with self.assertRaises(ValueError):validate_weight(value)
        row={'assignment_probability':.541,'assigned_family':'profilin'}
        status=identity_probability.annotate([row],orthogroup_context_weight=.541)
        self.assertNotIn('assignment_probability',row)
        self.assertFalse(status['probability_interpretation'])

    def test_workflow_flags_are_optional_and_zero_is_not_lost(self):
        for name,required,handler in [('audit',['--annotation','ann','--orthology','ort'],'audit'),
            ('resume',['--handoff','hand','--context-dir','ctx'],'resume'),('continue',['--annotation','ann'],'downstream')]:
            with patch.object(workflow,handler) as action:
                self.assertEqual(workflow.main([name,*required,'--out','out']),0)
                self.assertIsNone(action.call_args.args[0].orthogroup_context_weight)
                self.assertEqual(workflow.main([name,*required,'--out','out','--orthogroup-context-weight','0']),0)
                self.assertEqual(action.call_args.args[0].orthogroup_context_weight,0)

    def test_independent_evidence_invariance_excludes_orthology_dependent_explanations(self):
        base=self.call();zero=self.call(0)
        for row in (base,zero):row.update(row_origin='original_stage1',original_hypothesis_row=1)
        self.assertNotEqual(base['naming_evidence'],zero['naming_evidence'])
        self.assertEqual(invariants([base]),invariants([zero]))


if __name__=='__main__':unittest.main()
