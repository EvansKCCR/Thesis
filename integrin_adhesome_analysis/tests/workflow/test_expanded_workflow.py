"""Expanded rules reach the public workflow; synthetic inputs only."""

# Resolve the organized bundle independently of the caller working directory.
import sys as _bundle_sys
from pathlib import Path as _BundlePath
_BUNDLE_ROOT=next(p for p in _BundlePath(__file__).resolve().parents if (p/'pipeline_paths.py').is_file())
if str(_BUNDLE_ROOT) not in _bundle_sys.path:_bundle_sys.path.insert(0,str(_BUNDLE_ROOT))
from pipeline_paths import activate as _activate_bundle_paths
_activate_bundle_paths()
import argparse
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
import adhesome_workflow as w
import family_identity
import identity_probability
import evaluate_validation as ev
from tests.workflow import test_workflow as fixtures


class ExpandedWorkflowTests(unittest.TestCase):
    def test_cli_selector_is_available_for_audit_resume_and_continue(self):
        cases=[('audit',w.audit,['--annotation','a','--orthology','o']),
               ('resume',w.resume,['--handoff','h']),('continue',w.downstream,['--annotation','a'])]
        for command,fn,extra in cases:
            with self.subTest(command=command),patch.object(w,fn.__name__) as run:
                self.assertEqual(w.main([command,*extra,'--out','new','--family-coverage','expanded']),0)
                self.assertEqual(run.call_args.args[0].family_coverage,'expanded')

    def test_expanded_is_the_default_for_all_main_commands(self):
        cases=[('audit',w.audit,['--annotation','a','--orthology','o']),
               ('resume',w.resume,['--handoff','h']),('continue',w.downstream,['--annotation','a'])]
        for command,fn,extra in cases:
            with self.subTest(command=command),patch.object(w,fn.__name__) as run:
                self.assertEqual(w.main([command,*extra,'--out','new']),0)
                self.assertEqual(run.call_args.args[0].family_coverage,'expanded')
        self.assertEqual(w.audit_configuration(),w.audit_configuration('expanded'))

    def test_new_rule_is_reached_and_default_library_is_restored(self):
        spec=importlib.util.spec_from_file_location('expanded_test_audit',w.AUDIT)
        module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
        from classification_contract import stage_one_workbook
        original=family_identity.LIBRARY
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);source=root/'input.xlsx'
            stage_one_workbook([dict(Specie='AnyTarget',protein_id='AnyTarget__test',family='component_reference_candidate',
                length=500,architecture='PF07679:100-180',NCBI_CDD_architecture='IgI_1_Palladin_C:100-180',
                Orthology_HOG_reference_species_count=1)],source)
            for mode,expected in [('original','unassigned'),('expanded','palladin_like')]:
                rules_path,identity_path=w.audit_configuration(mode)
                rules=json.loads(rules_path.read_text());rules['_sha256']=w.digest(rules_path)
                rules['_source_path']=str(rules_path.resolve())
                output=root/mode;output.mkdir()
                meta=module.run_xlsx(source,output/'Stage_2_output.xlsx','candidates',rules,identity_library=identity_path)
                rows=w.read(output/'audited_family_assignments.tsv')
                self.assertEqual(rows[0]['putative_audited_family'],expected)
                self.assertEqual(meta['identity_library_sha256'],w.digest(identity_path))
                self.assertEqual(ev.audit_rule_paths(output),(rules_path.resolve(),identity_path.resolve()))
                if mode=='expanded':self.assertEqual(rows[0]['grade'],'A')
                self.assertEqual(family_identity.LIBRARY,original)

    def test_old_probability_fit_is_not_transferred_to_expanded_rules(self):
        with tempfile.TemporaryDirectory() as folder:
            model=Path(folder)/'model.json'
            model.write_text(json.dumps(dict(code_fingerprint=identity_probability.fingerprint())))
            _,expanded=w.audit_configuration('expanded')
            with family_identity.using_library(expanded):
                row=dict(assigned_family='palladin_like')
                self.assertEqual(identity_probability.annotate([row],model)['status'],'probability_model_stale_after_scoring_change')
                self.assertNotIn('assignment_probability',row)

    def test_recorded_configuration_cannot_change_silently(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);f=root/'families.json';i=root/'identities.json'
            f.write_text('{}');i.write_text('{}')
            w.dump(root/'Stage_2_output.manifest.json',dict(family_rules_path=str(f),identity_library_path=str(i),
                rule_sha256=w.digest(f),identity_library_sha256=w.digest(i)))
            self.assertEqual(ev.audit_rule_paths(root),(f,i))
            i.write_text('{"changed":true}')
            with self.assertRaisesRegex(ValueError,'configuration differs'):
                ev.audit_rule_paths(root)

    def test_public_audit_selects_and_seals_both_expanded_files(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);h,annotation,panel=fixtures.WorkflowTests().fixture(root)
            ort=root/'orthology';(ort/'candidates').mkdir(parents=True)
            w.write(ort/'candidates/candidate_report.tsv',[],['protein_id'])
            w.dump(ort/'candidates/provenance.json',{})
            w.seal(ort,'orthology',panel=str(panel),panel_manifest=w.digest(panel/'workflow_manifest.json'))
            output=root/'expanded'
            with patch('audit_join.join',return_value=(w.read(h/'hypotheses.tsv'),{})), \
                 patch.object(w,'evaluate_controls',return_value=dict(status='validation_target_not_present')):
                self.assertEqual(w.main(['audit','--annotation',str(annotation),'--orthology',str(ort),
                    '--out',str(output)]),0)
            meta=w.verify(output,'audit')
            families,identities=w.audit_configuration('expanded')
            self.assertEqual(meta['family_coverage'],'expanded')
            self.assertEqual(meta['rules_sha256'],w.digest(families))
            self.assertEqual(meta['identity_library_sha256'],w.digest(identities))
            self.assertEqual((output/'configuration/family_rules.json').read_bytes(),families.read_bytes())
            self.assertEqual((output/'configuration/family_identity_rules.json').read_bytes(),identities.read_bytes())
            self.assertTrue((output/'Stage_2_output.xlsx').is_file())
            self.assertEqual(json.loads((output/'logs/command_001.json').read_text())['status'],'complete')


if __name__=='__main__':unittest.main()
