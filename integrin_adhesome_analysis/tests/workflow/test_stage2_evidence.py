"""Synthetic tests of observed evidence, hierarchical gates and experiment freeze."""
import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
import family_identity as identity
import stage2_evidence as evidence
from tools import optimize_stage2 as experiment


class EvidenceTests(unittest.TestCase):
    def row(self,**changes):
        row=dict(protein_id='target__unknown',length=200,family='unknown',
                 NCBI_CDD_hits=json.dumps([dict(accession='pfam01146',short_name='Caveolin',start=30,end=150)]))
        row.update(changes);return row

    def test_exact_accession_import_retains_provenance_and_deduplicates(self):
        config=copy.deepcopy(identity.library());config['import_observed_pfam_accessions']=True
        for accession in ('pfam01146','PF01146.5','PFAM01146'):
            row=self.row(NCBI_CDD_hits=json.dumps([dict(accession=accession,short_name='observed_model',start=30,end=150)]))
            with identity.using_configuration(config):
                hits=identity.combined_domains(row,[])
                model=next(h for h in hits if h['model']=='PF01146')
                self.assertEqual(model['sources'],{'CDD'})
                self.assertTrue(model['cdd_pfam_accession'])
                merged=identity.combined_domains(row,[('PF01146',32,148)])
                self.assertEqual(len(identity.distinct_footprints(merged,{'PF01146'})),1)
                self.assertEqual(next(h for h in merged if h['model']=='PF01146')['sources'],{'Pfam','CDD'})

    def test_no_import_from_unverified_text_or_invalid_coordinates(self):
        config=copy.deepcopy(identity.library());config['import_observed_pfam_accessions']=True
        with identity.using_configuration(config):
            for accession,start,end in [('PF01146 extra',30,150),('pfam01146',0,150),('pfam01146',30,250)]:
                row=self.row(NCBI_CDD_hits=json.dumps([dict(accession=accession,short_name='unmapped',start=start,end=end)]))
                self.assertNotIn('PF01146',{h['model'] for h in identity.combined_domains(row,[])})
        with identity.using_configuration(dict(config,import_observed_pfam_accessions=False)):
            self.assertNotIn('PF01146',{h['model'] for h in identity.combined_domains(self.row(NCBI_CDD_hits=json.dumps([
                dict(accession='pfam01146',short_name='unmapped',start=30,end=150)])),[])})

    def test_recruitment_keeps_ids_and_requires_specific_model(self):
        rules={'_stage2_profile':'calibrated_v1'}
        with evidence.annotation_context(rules):
            proposals=evidence.candidate_hypotheses(self.row(),{'unknown'},rules,lambda row:[])
            self.assertEqual([r['family'] for r in proposals],['caveolin_like'])
            self.assertEqual(proposals[0]['protein_id'],'target__unknown')
            for changes in (dict(NCBI_CDD_hits='',NCBI_CDD_architecture='PH:30-150'),
                            dict(NCBI_CDD_hits='',Orthology_HOG_reference_species_count=12,protein_name='Caveolin')):
                self.assertEqual(evidence.candidate_hypotheses(self.row(**changes),{'unknown'},rules,lambda row:[]),[])

    def test_profile_excludes_whole_confounded_family_and_restores_library(self):
        base=identity.library()
        with evidence.annotation_context({'_stage2_profile':'calibrated_v1'}):
            self.assertNotIn('evidence_limk',{r['id'] for r in identity.library()['rules']})
            self.assertTrue(identity.library()['import_observed_pfam_accessions'])
        self.assertEqual(identity.library(),base)
        with evidence.annotation_context({'_stage2_profile':'evidence_v1'}):
            self.assertIn('evidence_limk',{r['id'] for r in identity.library()['rules']})

    def test_original_does_not_recruit_and_unknown_profile_fails(self):
        self.assertEqual(evidence.candidate_hypotheses(self.row(),set(),{'_stage2_profile':'original'},lambda row:[]),[])
        with self.assertRaises(ValueError):evidence.profile({'_stage2_profile':'invalid'})

    def test_generic_partial_is_unresolved_and_full_core_is_unchanged(self):
        result=dict(score=2,cdd=False,note='Shared PH only',rule_trace=[])
        rules={'_stage2_profile':'calibrated_v1'}
        self.assertEqual(evidence.partial_confirmation({},'unknown',result,rules)['score'],0)
        full=dict(result,score=4)
        self.assertIs(evidence.partial_confirmation({},'unknown',full,rules),full)
        self.assertIs(evidence.partial_confirmation({},'unknown',result,{'_stage2_profile':'original'}),result)

    def test_partial_core_remains_partial_and_length_failure_cannot_rescue(self):
        trace=dict(family='PDE4_like',diagnostic_cdd_matched=True,pfam_anchor_passed=True,
                   coordinates_passed=True,length_passed=True,maximum_counts_passed=True,
                   domain_groups=[dict(observed=1,required=1),dict(observed=0,required=1)],rule_id='evidence_pde4')
        rules={'_stage2_profile':'calibrated_v1'}
        with evidence.annotation_context(rules):
            result=dict(score=0,cdd=False,note='Incomplete PDE4',rule_trace=[trace])
            partial=evidence.partial_confirmation({},'PDE4_like',result,rules)
            self.assertEqual(partial['score'],2);self.assertFalse(partial['complete'])
            trace['length_passed']=False
            self.assertEqual(evidence.partial_confirmation({},'PDE4_like',result,rules)['score'],0)

    def test_validation_start_blocks_calibration_extension_even_if_incomplete(self):
        with tempfile.TemporaryDirectory() as directory:
            folder=Path(directory)
            self.assertFalse(experiment.validation_started(folder))
            (folder/'validation_started.json').write_text('{}')
            self.assertTrue(experiment.validation_started(folder))
            with self.assertRaisesRegex(ValueError,'after validation'):
                experiment.finalize_freeze(None,None,folder)

    def test_selection_rejects_added_false_positive(self):
        def pred(pid,positive,accepted):return dict(protein_id=pid,control_class='positive' if positive else 'hard_negative',
            supported=accepted,predicted_family='family',expected_families='family' if positive else '')
        rows=[pred('P',True,True),pred('N',False,False)]
        baseline=dict(arm='original',profile='original',context_weight=1.,TP=1,FP=0,F1=1.,
            hard_negative_FP=0,background_negative_FP=0,family=dict(FP=0,F1=1.))
        unsafe=dict(baseline,arm='unsafe',profile='evidence_v1',FP=1,hard_negative_FP=1)
        safe=dict(baseline,arm='safe',profile='calibrated_v1')
        with tempfile.TemporaryDirectory() as directory:
            selected=experiment.freeze_selection(Path(directory),[baseline,unsafe,safe],
                dict(original=rows,unsafe=[rows[0],pred('N',False,True)],safe=rows))
        self.assertEqual(selected['selected_arm'],'safe')
        self.assertFalse(unsafe['eligible_by_calibration'])


if __name__=='__main__':unittest.main()
