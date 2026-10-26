
# Resolve the organized bundle independently of the caller working directory.
import sys as _bundle_sys
from pathlib import Path as _BundlePath
_BUNDLE_ROOT=next(p for p in _BundlePath(__file__).resolve().parents if (p/'pipeline_paths.py').is_file())
if str(_BUNDLE_ROOT) not in _bundle_sys.path:_bundle_sys.path.insert(0,str(_BUNDLE_ROOT))
from pipeline_paths import activate as _activate_bundle_paths
_activate_bundle_paths()
import unittest
from collections import defaultdict
from family_identity_review import review, interval_count


class IdentityReviewTests(unittest.TestCase):
    def test_overlapping_domains_do_not_create_tandem_architecture(self):
        self.assertEqual(interval_count([('PF00307', 1, 80), ('PF00307', 20, 90)], 'PF00307'), 1)

    def test_ch_architecture_preserves_original_family_and_competing_evidence(self):
        families = {'actin': [dict(id='x', sequence='A'*350)]}
        regions = defaultdict(list, {('actin', 'x'): [('PF00307', 1, 80), ('PF00307', 100, 180)]})
        row = review(families, regions)[0]
        self.assertEqual(row['review_priority'], 'architecture_followup')
        self.assertEqual(row['family'], 'actin')
        self.assertEqual(row['identity_status'], 'unresolved_original_family_retained')
        regions['actin', 'x'].append(('PF00022', 1, 300))
        self.assertEqual(review(families, regions)[0]['review_priority'], 'unresolved_review')


if __name__ == '__main__':
    unittest.main()
