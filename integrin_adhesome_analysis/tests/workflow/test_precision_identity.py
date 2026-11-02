"""Discriminating CDD evidence and optional evolutionary identity review."""

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
import csv
import json
import tempfile
import unittest

from cdd_evidence import interpretation,diagnostic
from family_identity import using_configuration,scoped_refinement,combined_domains
from phylogenetic_identity import apply,sha
from family_catalogue import report,canonical_family


class PrecisionIdentityTests(unittest.TestCase):
    def test_recovery_can_explicitly_disable_optional_identity_review(self):
        from unittest.mock import patch
        import adhesome_workflow as core
        with patch.object(core,'recover') as recover:
            self.assertEqual(core.main(['recover','--run','fixture','--skip-phylogenetic-identity-resolution']),0)
            self.assertFalse(recover.call_args.args[0].phylogenetic_identity_resolution)
    def test_src_aliases_do_not_register_the_retracted_typo(self):
        for name in ('src','Src','SRC'):self.assertEqual(canonical_family(name),'Src')
        self.assertEqual(canonical_family('SCR'),'SCR')

    def test_specificity_is_reported_category_not_model_name(self):
        for category in ('non-specific','superfamily','multidomain',''):
            hit=dict(short_name='PTKc_FAK',accession='cd00123',hit_type=category,start=1,end=90)
            self.assertFalse(diagnostic(hit,'FAK'))
        self.assertTrue(diagnostic(dict(short_name='PTKc_FAK',accession='cd00123',hit_type='specific'),'FAK',True))

    def test_model_threshold_failure_overrides_category(self):
        hit=dict(short_name='FAK',accession='cd00123',hit_type='specific',bitscore=50,specific_bitscore_cutoff=60)
        self.assertFalse(interpretation(hit)['specific'])
        hit.update(bitscore=60);self.assertTrue(interpretation(hit)['specific'])
        hit.update(accession='pfam00069');self.assertFalse(diagnostic(hit,'FAK',True))
        hit.update(accession='cd123',evalue='nan');self.assertFalse(diagnostic(hit,'FAK'))
        hit.update(accession='cl00123',evalue=1e-30);self.assertFalse(diagnostic(hit,'FAK'))

    def test_exact_legacy_specific_export_marker_is_traceable(self):
        hit=dict(short_name='LAR',accession='cd123',interpretation='CDD_specific_domain_prediction_not_family_or_motif_validation')
        self.assertEqual(interpretation(hit)['category_source'],'legacy_specific_only_export')
        self.assertTrue(interpretation(hit)['specific'])
        del hit['interpretation'];self.assertFalse(interpretation(hit)['specific'])

    def test_all_complete_diagnostic_alternatives_are_compared(self):
        entries=[dict(id=f,parents=['parent'] if f=='A' else ['another_dispatch'],family=f,protein_name=f,
            required_groups=[dict(models=['PF00069'],count=1)],specific_cdd_pattern=f,ordered_groups=[])
            for f in ('A','B')]
        config=dict(precision_safeguards=True,pfam_cdd_aliases={},additional_cdd_domains={},rules=entries)
        row=dict(length=100,NCBI_CDD_hits=[dict(short_name=f,start=1,end=90,accession='cd123',hit_type='specific') for f in ('A','B')])
        with using_configuration(config):
            hits=combined_domains(row,[('PF00069',1,90)])
            result=scoped_refinement(row,'parent',hits)
            self.assertIsNone(result['rule']);self.assertEqual(result['competing_families'],['A','B'])
            entries[1]['no_auto_promotion']=True
            result=scoped_refinement(row,'parent',hits)
            self.assertEqual(result['rule']['id'],'A')

    def test_all_stage1_families_have_explicit_handling(self):
        result=report()
        self.assertEqual(len(result['stage1_policies']),result['stage1_families'])
        self.assertFalse(result['duplicate_rule_ids']);self.assertFalse(result['missing_registry_routes'])

    def test_lar_merlin_and_trio_boundaries_with_reported_specific_hits(self):
        import sys
        from pipeline_paths import ROOT
        sys.path.insert(0,str(_layout_resource(ROOT/'annotation')))
        from integrated_adhesome_pipeline import classify_rows
        rules=json.loads((_layout_resource(ROOT/'annotation/family_rules.json')).read_text())
        rules['_stage2_profile']='precision_v1'
        cases=[('lar','LAR_RPTP','PF07679:10-90;PF00041:110-180;PF00102:250-450','PTPc_LAR',250,450),
            ('merlin','NF2_merlin_like','PF00373:10-200;PF00769:500-600','FERM_F1_Merlin',10,100),
            ('trio','RhoGEF_like','PF00621:100-250;PF00435:300-360;PF00435:400-460;PF00435:500-560;PF00169:600-700','PH1_Kalirin_Trio_like',600,700)]
        for pid,expected,architecture,name,start,end in cases:
            row=dict(protein_id=pid,Specie='novel',family='component_reference_candidate',length=1000,
                architecture=architecture,Orthology_HOG_reference_species_count=2,
                NCBI_CDD_hits=[dict(accession='cd12345',short_name=name,start=start,end=end,hit_type='specific')])
            _,calls=classify_rows([row],rules)
            self.assertEqual(calls[0]['assigned_family'],expected)
            self.assertNotIn(calls[0]['assigned_family'],('ERM_like','spectrin_like'))

    def test_missing_phylogeny_does_not_turn_supported_family_review_into_failure(self):
        import sys
        from pipeline_paths import ROOT
        sys.path.insert(0,str(_layout_resource(ROOT/'annotation')))
        from family_set_review_gate import gate
        registry=dict(family_sets={'test':dict(review_flags=[],user_inspection_required=False)})
        row=dict(family='test',stage2_architecture_compatible=True,stage2_cdd_pattern_match=True,
            stage2_named_curated_cd_match=True,ncbi_cdd_specific_hit_threshold_met=True)
        self.assertEqual(gate(row,registry)['user_inspection_required'],'false')
        self.assertEqual(gate(row,registry,True)['user_inspection_required'],'true')

    def test_phylogenetic_identity_disabled_does_not_read_missing_evidence(self):
        row=dict(protein_id='target',assigned_family='FAK',grade='A',identity_outcome='identity_supported')
        status=apply([row],_layout_resource(Path('does_not_exist.tsv')),False)
        self.assertEqual(status['status'],'not_requested');self.assertEqual(row['grade'],'A')
        self.assertEqual(row['phylogenetic_identity_status'],'not_requested')

    def test_optional_failure_cannot_block_family_call(self):
        row=dict(protein_id='target',assigned_family='FAK',grade='A')
        self.assertEqual(apply([row],_layout_resource(Path('absent.tsv')),True)['status'],'unavailable')
        self.assertEqual(row['assigned_family'],'FAK');self.assertEqual(row['grade'],'A')

    def test_conflicting_tree_triggers_review_without_relabelling(self):
        with tempfile.TemporaryDirectory() as td:
            root=_layout_resource(Path(td));tree=_layout_resource(root/'tree.nwk');tree.write_text('(target,ref1,ref2);')
            hit=dict(protein_id='target',sequence_sha256='sequence',tree_file=str(tree),tree_sha256=sha(tree),
                family='Src',ufboot=99,sh_alrt=95,reference_count=2,reference_taxa_count=2,
                reference_label_source='curated reference labels',clade_id='clade1',placement_status='pure_supported_clade')
            evidence=_layout_resource(root/'placements.tsv')
            with _layout_resource(evidence).open('w',newline='') as h:
                w=csv.DictWriter(h,fieldnames=list(hit),delimiter='\t');w.writeheader();w.writerow(hit)
            row=dict(protein_id='target',sequence_sha256='sequence',assigned_family='FAK',grade='A')
            apply([row],evidence,True)
            self.assertEqual(row['identity_outcome'],'conflicting_identity')
            self.assertEqual(row['assigned_family'],'FAK');self.assertEqual(row['grade'],'A')


if __name__=='__main__':unittest.main()
