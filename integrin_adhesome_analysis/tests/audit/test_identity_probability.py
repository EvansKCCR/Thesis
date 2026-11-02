
# Resolve source resources independently of the caller working directory.
import sys as _layout_sys
from pathlib import Path as _LayoutPath
_LAYOUT_BUNDLE_ROOT = next(p for p in _LayoutPath(__file__).resolve().parents
                          if (p/'config/layout_paths.json').is_file() and (p/'pipeline_paths.py').is_file())
if str(_LAYOUT_BUNDLE_ROOT) not in _layout_sys.path:
    _layout_sys.path.insert(0, str(_LAYOUT_BUNDLE_ROOT))
from pipeline_paths import activate as _layout_activate, resource_path as _layout_resource
_layout_activate()



import sys as _sys
from pathlib import Path as _Path
_PIPELINE_ROOT=_layout_resource(_Path(__file__)).resolve().parents[2]
_sys.path.insert(0,str(_PIPELINE_ROOT))
_sys.path.insert(0,str(_layout_resource(_PIPELINE_ROOT/"annotation")))
import json
from pathlib import Path
import tempfile
import unittest
from identity_probability import annotate,blocks,estimate,fingerprint
from calibrate_identity import expected_families


class IdentityProbabilityTests(unittest.TestCase):
    def test_isotonic_pooling_is_monotone_and_data_driven(self):
        samples=[dict(score=x,correct=y) for x,y in [(1,0),(2,1),(3,0),(4,1)]]
        fitted=blocks(samples,min_bin=1)
        self.assertEqual([(b['n'],b['correct']) for b in fitted],[(1,0),(2,1),(1,1)])
        rates=[b['correct']/b['n'] for b in fitted];self.assertEqual(rates,sorted(rates))

    def test_probability_does_not_extrapolate_beyond_control_scores(self):
        samples=[dict(score=7,correct=1,stratum='same') for i in range(20)]
        self.assertEqual(estimate(samples,dict(score=9,stratum='same'))['status'],'outside_observed_calibration_score_range')

    def test_query_sequence_labels_are_excluded_in_full_and_do_not_set_names(self):
        with tempfile.TemporaryDirectory() as td:
            path=_layout_resource(Path(td)/'model.json')
            samples=[dict(score=7,correct=0,stratum='strong|family',sequence_sha256=str(i)) for i in range(20)]
            samples.extend([dict(score=7,correct=1,stratum='strong|family',sequence_sha256='SELF')]*2)
            path.write_text(json.dumps(dict(code_fingerprint=fingerprint(),input_hashes={},minimum_controls=20,
                context={'reference_overlap':True},samples=samples,cross_validation={})))
            rows=[dict(assigned_family='talin',total_score_0_10=7,family_confidence='strong',protein_name_scope='family',
                sequence_sha256='SELF',protein_name='Putative Talin'),
                dict(assigned_family='talin',total_score_0_10=7,family_confidence='strong',protein_name_scope='family',
                sequence_sha256='NEW',protein_name='Putative Talin')]
            result=annotate(rows,path)
            self.assertEqual(result['applied'],2)
            self.assertEqual(rows[0]['assignment_probability'],0)
            self.assertEqual(rows[0]['probability_calibration_n'],20)
            self.assertGreater(rows[0]['assignment_probability_upper_95'],0)
            self.assertIn('leave_one_sequence_out',rows[0]['probability_method'])
            self.assertAlmostEqual(rows[1]['assignment_probability'],2/22,places=6)
            self.assertTrue(all(r['protein_name']=='Putative Talin' for r in rows))

    def test_negative_adhesome_control_can_have_correct_hip1_identity(self):
        truth=expected_families(dict(genes='HIP1',protein_name='Huntingtin-interacting protein 1',
            **{'class':'hard_negative control'}),{})
        self.assertIn('HIP1_HIP1R_like',truth)
        self.assertNotIn('talin',truth)

    def test_stale_probability_model_does_not_supply_a_number(self):
        with tempfile.TemporaryDirectory() as td:
            path=_layout_resource(Path(td)/'model.json');path.write_text(json.dumps(dict(code_fingerprint={})))
            row=dict(assigned_family='talin')
            self.assertEqual(annotate([row],path)['status'],'probability_model_stale_after_scoring_change')
            self.assertNotIn('assignment_probability',row)

    def test_experimental_metadata_never_changes_calibrated_assignment_probability(self):
        with tempfile.TemporaryDirectory() as td:
            path=_layout_resource(Path(td)/'model.json')
            samples=[dict(score=7,correct=int(i<15),stratum='strong|family',sequence_sha256=str(i)) for i in range(20)]
            path.write_text(json.dumps(dict(code_fingerprint=fingerprint(),input_hashes={},minimum_controls=20,
                context={'positive_family_scope':['talin'],'positive_control_only':True},samples=samples,cross_validation={})))
            base=dict(assigned_family='talin',total_score_0_10=7,family_confidence='strong',protein_name_scope='family',
                sequence_sha256='NEW',protein_name='Putative Talin')
            rows=[dict(base),dict(base,experimental_evidence='',adhesome_experiment_verified=False),
                dict(base,experimental_evidence='experiment',adhesome_experiment_verified=True)]
            self.assertEqual(annotate(rows,path)['applied'],3)
            self.assertEqual({row['assignment_probability'] for row in rows},{0.75})
            self.assertEqual(len({row['assignment_probability_lower_95'] for row in rows}),1)


if __name__=='__main__':unittest.main()
