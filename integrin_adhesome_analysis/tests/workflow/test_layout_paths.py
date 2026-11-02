"""Layout regressions use synthetic inputs; no stored run is executed."""

# Resolve source resources independently of the caller working directory.
import sys as _layout_sys
from pathlib import Path as _LayoutPath
_LAYOUT_BUNDLE_ROOT = next(p for p in _LayoutPath(__file__).resolve().parents
                          if (p/'config/layout_paths.json').is_file() and (p/'pipeline_paths.py').is_file())
if str(_LAYOUT_BUNDLE_ROOT) not in _layout_sys.path:
    _layout_sys.path.insert(0, str(_LAYOUT_BUNDLE_ROOT))
from pipeline_paths import activate as _layout_activate, resource_path as _layout_resource
_layout_activate()

from pathlib import Path
import hashlib
import json
import tempfile
import unittest
from unittest.mock import patch

import pipeline_paths
from calibrate_identity import calibration_inputs
import control_exclusion
import identity_probability
from schisto_adhesome.workflows import algorithm_inventory


class LayoutPathsTests(unittest.TestCase):
    def test_every_recorded_destination_exists_and_is_idempotent(self):
        layout=json.loads((pipeline_paths.CONFIG/'layout_paths.json').read_text())
        for previous,current in layout['file_map'].items():
            target=pipeline_paths.ROOT/current
            with self.subTest(previous=previous):
                self.assertTrue(target.is_file())
                self.assertEqual(pipeline_paths.recorded_resource(previous),target)
                self.assertEqual(pipeline_paths.resource_path(target),target)

    def test_layout_receipt_accepts_only_the_exact_recorded_transition(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);(root/'main').mkdir();(root/'maintenance').mkdir()
            target=root/'main/runtime.py';target.write_text('relocated source')
            old=hashlib.sha256(b'old source layout').hexdigest()
            new=hashlib.sha256(target.read_bytes()).hexdigest()
            receipt={'files':[{'old_path':'old.py','new_path':'main/runtime.py',
                              'old_sha256':old,'new_sha256':new}]}
            (root/'maintenance/layout_migration.json').write_text(json.dumps(receipt))
            with patch.object(pipeline_paths,'ROOT',root),patch.object(pipeline_paths,'MAINTENANCE',root/'maintenance'),\
                 patch.object(pipeline_paths,'_LAYOUT',{'file_map':{'old.py':'main/runtime.py'},'directory_map':{}}),\
                 patch.object(pipeline_paths,'_RECEIPT_CACHE',None):
                self.assertTrue(pipeline_paths.recorded_digest_matches('old.py',old))
                self.assertFalse(pipeline_paths.recorded_digest_matches('old.py','0'*64))
                target.write_text('later unapproved code change')
                self.assertFalse(pipeline_paths.recorded_digest_matches('old.py',old))

    def test_external_and_run_paths_are_not_relocated(self):
        with tempfile.TemporaryDirectory() as td:
            external=Path(td)/'annotation/family_rules.json'
            self.assertEqual(pipeline_paths.resource_path(external),external)
        run=pipeline_paths.ROOT/'runs/annotation/family_rules.json'
        self.assertEqual(pipeline_paths.resource_path(run),run)
        self.assertEqual(pipeline_paths.resource_path(Path('runs/new/annotation')),Path('runs/new/annotation'))

    def test_historical_model_keys_resolve_without_rewriting_model(self):
        expected={
            '01_detection/config/family_rules.json':'config/detection/family_rules.json',
            '02_calibration_set/calibration_set.xlsx':'validation/calibration_set.xlsx',
            'calibrate_identity.py':'resource_library/scripts/calibrate_identity.py',
        }
        for recorded,current in expected.items():
            self.assertEqual(pipeline_paths.recorded_resource(recorded),pipeline_paths.ROOT/current)
        absolute=_layout_resource(pipeline_paths.ROOT/'saved_run'/'unchanged.json')
        self.assertEqual(pipeline_paths.recorded_resource(absolute),absolute)

    def test_frozen_inventory_covers_relocated_classifier_modules(self):
        inventory=algorithm_inventory()
        for name in ('resource_library/scripts/audit_join.py','resource_library/scripts/stage2_evidence.py',
                     'resource_library/scripts/control_exclusion.py','main/pipeline_paths.py','config/annotation/family_rules.json'):
            self.assertEqual(inventory[name],hashlib.sha256((_layout_resource(pipeline_paths.ROOT/name)).read_bytes()).hexdigest())

    def test_shared_folder_selects_calibration_not_validation_or_auxiliary_workbook(self):
        from openpyxl import Workbook
        with tempfile.TemporaryDirectory() as td:
            root=_layout_resource(Path(td));folder=_layout_resource(root/'validation');folder.mkdir()
            book=Workbook();sheet=book.active;sheet.title='control_identity'
            sheet.append(['class','accession','trusted_family','protein_name'])
            sheet.append(['hard_negative control','N1','test_family','Negative'])
            book.save(_layout_resource(folder/'calibration_set.xlsx'));book.close()
            (_layout_resource(folder/'calibration.faa')).write_text('>N1\nMACE\n')
            (_layout_resource(folder/'validation.faa')).write_text('>N1\nMWWW\n')
            (_layout_resource(folder/'components_source.xlsx')).write_text('auxiliary workbook, never opened')
            self.assertEqual(calibration_inputs(folder),(_layout_resource(folder/'calibration_set.xlsx'),_layout_resource(folder/'calibration.faa')))
            control_exclusion.controls.cache_clear()
            try:
                with patch.object(control_exclusion,'ROOT',root):bank=control_exclusion.controls()
                self.assertIn(hashlib.sha256(b'MACE').hexdigest(),bank['negative_by_hash'])
                self.assertNotIn(hashlib.sha256(b'MWWW').hexdigest(),bank['negative_by_hash'])
            finally:control_exclusion.controls.cache_clear()

    def test_incomplete_named_pair_cannot_select_other_panel(self):
        with tempfile.TemporaryDirectory() as td:
            root=_layout_resource(Path(td));folder=_layout_resource(root/'validation');folder.mkdir()
            (_layout_resource(folder/'calibration_set.xlsx')).touch();(_layout_resource(folder/'validation.faa')).touch()
            with self.assertRaisesRegex(ValueError,'Both calibration_set.xlsx and calibration.faa'):
                calibration_inputs(folder)
            control_exclusion.controls.cache_clear()
            try:
                with patch.object(control_exclusion,'ROOT',root):
                    with self.assertRaisesRegex(ValueError,'Missing calibration exclusion inputs'):
                        control_exclusion.controls()
            finally:control_exclusion.controls.cache_clear()

    def test_validation_inputs_do_not_change_calibration_exclusion_hashes(self):
        with tempfile.TemporaryDirectory() as td:
            root=_layout_resource(Path(td));(_layout_resource(root/'scripts')).mkdir();(_layout_resource(root/'validation')).mkdir()
            for name in ('identity_probability.py','stage2_resolution.py','control_exclusion.py','cdd_evidence.py','family_catalogue.py','phylogenetic_identity.py'):
                (_layout_resource(root/'scripts')/name).write_text('synthetic source')
            for name in ('calibration_set.xlsx','calibration.faa','validation.faa','components_source.xlsx'):
                (_layout_resource(root/'validation')/name).write_text('original fixture')
            with patch.object(identity_probability,'ROOT',root),patch.object(identity_probability,'DEFAULT_MODEL',_layout_resource(root/'absent_model.json')):
                before=identity_probability.resource_hashes()
                (_layout_resource(root/'validation/validation.faa')).write_text('changed validation input')
                (_layout_resource(root/'validation/components_source.xlsx')).write_text('changed auxiliary workbook')
                self.assertEqual(before,identity_probability.resource_hashes())
                (_layout_resource(root/'validation/calibration.faa')).write_text('changed exclusion input')
                self.assertNotEqual(before,identity_probability.resource_hashes())


if __name__=='__main__':unittest.main()
