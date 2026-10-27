"""Layout regressions use synthetic inputs; no stored run is executed."""
from pathlib import Path
import hashlib
import tempfile
import unittest
from unittest.mock import patch

import pipeline_paths
from calibrate_identity import calibration_inputs
import control_exclusion
import identity_probability
from schisto_adhesome.workflows import algorithm_inventory


class LayoutPathsTests(unittest.TestCase):
    def test_historical_model_keys_resolve_without_rewriting_model(self):
        expected={
            '01_detection/config/family_rules.json':'detection/config/family_rules.json',
            '02_calibration_set/calibration_set.xlsx':'validation/calibration_set.xlsx',
            'calibrate_identity.py':'scripts/calibrate_identity.py',
        }
        for recorded,current in expected.items():
            self.assertEqual(pipeline_paths.recorded_resource(recorded),pipeline_paths.ROOT/current)
        absolute=pipeline_paths.ROOT/'saved_run'/'unchanged.json'
        self.assertEqual(pipeline_paths.recorded_resource(absolute),absolute)

    def test_frozen_inventory_covers_relocated_classifier_modules(self):
        inventory=algorithm_inventory()
        for name in ('scripts/audit_join.py','scripts/stage2_evidence.py',
                     'scripts/control_exclusion.py','pipeline_paths.py'):
            self.assertEqual(inventory[name],hashlib.sha256((pipeline_paths.ROOT/name).read_bytes()).hexdigest())

    def test_shared_folder_selects_calibration_not_validation_or_auxiliary_workbook(self):
        from openpyxl import Workbook
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);folder=root/'validation';folder.mkdir()
            book=Workbook();sheet=book.active;sheet.title='control_identity'
            sheet.append(['class','accession','trusted_family','protein_name'])
            sheet.append(['hard_negative control','N1','test_family','Negative'])
            book.save(folder/'calibration_set.xlsx');book.close()
            (folder/'calibration.faa').write_text('>N1\nMACE\n')
            (folder/'validation.faa').write_text('>N1\nMWWW\n')
            (folder/'components_source.xlsx').write_text('auxiliary workbook, never opened')
            self.assertEqual(calibration_inputs(folder),(folder/'calibration_set.xlsx',folder/'calibration.faa'))
            control_exclusion.controls.cache_clear()
            try:
                with patch.object(control_exclusion,'ROOT',root):bank=control_exclusion.controls()
                self.assertIn(hashlib.sha256(b'MACE').hexdigest(),bank['negative_by_hash'])
                self.assertNotIn(hashlib.sha256(b'MWWW').hexdigest(),bank['negative_by_hash'])
            finally:control_exclusion.controls.cache_clear()

    def test_incomplete_named_pair_cannot_select_other_panel(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);folder=root/'validation';folder.mkdir()
            (folder/'calibration_set.xlsx').touch();(folder/'validation.faa').touch()
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
            root=Path(td);(root/'scripts').mkdir();(root/'validation').mkdir()
            for name in ('identity_probability.py','stage2_resolution.py','control_exclusion.py','cdd_evidence.py','family_catalogue.py','phylogenetic_identity.py'):
                (root/'scripts'/name).write_text('synthetic source')
            for name in ('calibration_set.xlsx','calibration.faa','validation.faa','components_source.xlsx'):
                (root/'validation'/name).write_text('original fixture')
            with patch.object(identity_probability,'ROOT',root),patch.object(identity_probability,'DEFAULT_MODEL',root/'absent_model.json'):
                before=identity_probability.resource_hashes()
                (root/'validation/validation.faa').write_text('changed validation input')
                (root/'validation/components_source.xlsx').write_text('changed auxiliary workbook')
                self.assertEqual(before,identity_probability.resource_hashes())
                (root/'validation/calibration.faa').write_text('changed exclusion input')
                self.assertNotEqual(before,identity_probability.resource_hashes())


if __name__=='__main__':unittest.main()
