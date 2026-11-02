"""Observed Pfam accessions remain usable beside their architecture aliases."""

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
import sys
from pathlib import Path
import unittest

sys.path.insert(0,str(_layout_resource(Path(__file__)).resolve().parents[2]))
from family_identity import combined_domains,distinct_footprints,ordered_core


class ObservedModelPreservationTests(unittest.TestCase):
    def test_ig_and_pet_aliases_preserve_observed_accessions(self):
        for accession,alias in [('PF07679','IG'),('PF13927','IG'),('PF06297','PET')]:
            with self.subTest(accession=accession):
                hits=combined_domains({'length':300},[(accession,20,90)])
                self.assertEqual({h['model'] for h in hits},{accession,alias})
                self.assertTrue(all(h['sources']=={'Pfam'} for h in hits))
                self.assertEqual(len(distinct_footprints(hits,{accession,alias})),1)

    def test_alias_cannot_supply_a_second_copy_or_ordered_domain(self):
        hits=combined_domains({'length':300},[('PF07679',20,90)])
        rule=dict(required_groups=[dict(models=['PF07679','IG'],count=2)],ordered_groups=[0])
        self.assertFalse(ordered_core(hits,rule))
        two=combined_domains({'length':300},[('PF07679',20,90),('PF07679',110,180)])
        self.assertEqual(len(distinct_footprints(two,{'PF07679','IG'})),2)
        self.assertTrue(ordered_core(two,rule))

    def test_invalid_or_coordinate_free_observation_does_not_gain_coordinates(self):
        hits=combined_domains({'length':300},[('PF13927',None,None)])
        self.assertEqual(len(distinct_footprints(hits,{'PF13927','IG'})),0)


if __name__=='__main__':unittest.main()
