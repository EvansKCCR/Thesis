"""Synthetic checks for automatic post hoc reporting and benchmark RBHs."""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'resource_library'))
from schisto_adhesome import negative_homology
from schisto_adhesome.stage_performance_report import generate


def write_tsv(path, rows, fields=None):
    fields=fields or list(dict.fromkeys(key for row in rows for key in row))
    with Path(path).open('w',encoding='utf-8',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=fields,delimiter='\t')
        writer.writeheader();writer.writerows(rows)


class IntegratedRefinementReportTests(unittest.TestCase):
    def test_stage_report_keeps_cohort_and_exports_evidence_and_deltas(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);reports=root/'reports';reports.mkdir();posthoc=root/'posthoc';posthoc.mkdir()
            refinement=root/'refinement';refinement.mkdir();out=root/'stage'
            source=[dict(Specie='Val_Test',protein_id='p1',control_class='positive',expected_families='actin',
                stage1_recruited='True',stage1_family_hypotheses='actin',stage2_annotated='False',
                stage2_family='unassigned',supported='False',predicted_family='unassigned'),
                dict(Specie='Val_Test',protein_id='p2',control_class='negative',expected_families='',
                stage1_recruited='True',stage1_family_hypotheses='background',stage2_annotated='True',
                stage2_family='background',supported='False',predicted_family='background'),
                dict(Specie='Val_Test',protein_id='p3',control_class='positive',expected_families='integrin',
                stage1_recruited='False',stage1_family_hypotheses='',stage2_annotated='False',
                stage2_family='unassigned',supported='False',predicted_family='unassigned')]
            write_tsv(reports/'predictions.tsv',source)
            proposed=[dict(Specie='Val_Test',protein_id='p1',control_class='positive',expected_families='actin',
                posthoc_AB='True',posthoc_evaluation_family='actin',posthoc_AB_family_correct='True',
                original_supported='False',posthoc_status='reviewed',posthoc_tier='A'),
                dict(Specie='Val_Test',protein_id='p2',control_class='negative',expected_families='',
                posthoc_AB='False',posthoc_evaluation_family='',posthoc_AB_family_correct='False',
                original_supported='False',posthoc_status='reviewed',posthoc_tier='C'),
                dict(Specie='Val_Test',protein_id='p3',control_class='positive',expected_families='integrin',
                posthoc_AB='False',posthoc_evaluation_family='',posthoc_AB_family_correct='False',
                original_supported='False',posthoc_status='absent_from_audit',posthoc_tier='')]
            write_tsv(posthoc/'predictions.tsv',proposed)
            write_tsv(refinement/'Stage_2_posthoc_refinement.tsv',[
                {'Species':'Val_Test','Protein ID':'p1','posthoc_proposed_family_or_class':'actin',
                 'posthoc_evidence_rationale':'PF00022 plus CDD actin-family architecture','posthoc_evidence_tier':'A'},
                {'Species':'Val_Test','Protein ID':'p2','posthoc_proposed_family_or_class':'background_candidate',
                 'posthoc_evidence_rationale':'broad domain only','posthoc_evidence_tier':'C'}])
            (refinement/'run_manifest.json').write_text(json.dumps({'column_mapping':{'Specie':'Species','protein_id':'Protein ID'}}))
            write_tsv(refinement/'potential_false_positive_candidates.tsv',[
                {'Specie':'Val_Test','protein_id':'p3','audited_family':'ERM','domain_proposed_family':'Merlin',
                 'evidence_tier':'A','domain_evidence':'FERM architecture conflicts with canonical ERM signature'}])
            summary=generate(reports/'predictions.tsv',posthoc,refinement,out,'validation')
            metrics={row['stage']:row for row in csv.DictReader((out/'stage_performance.tsv').open(encoding='utf-8'),delimiter='\t')}
            changes=list(csv.DictReader((out/'candidate_changes.tsv').open(encoding='utf-8'),delimiter='\t'))
            self.assertEqual(metrics['Stage 1 detection']['n'],'3')
            self.assertEqual(metrics['Stage 3 plus post hoc A/B proposals']['TP'],'1')
            self.assertEqual(summary['posthoc_candidates_rescued'],1)
            self.assertEqual(summary['supported_assignments_flagged_as_potential_false_positive'],1)
            self.assertIn('PF00022',next(row for row in changes if row['protein_id']=='p1')['evidence'])
            self.assertTrue((out/'stage_changes.tsv').is_file())

    def test_negative_homology_adds_only_exact_unique_rbh_with_source_evidence(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);panel=root/'Benchmark'/'Input';panel.mkdir(parents=True)
            target_sequence='ACDEFGHIKLMNPQRSTVWY'*4
            query_sequence='YWVTSRQPNMLKIHGFEDCA'*4
            target=panel/'Schisto.faa';target.write_text(f'>target1\n{target_sequence}\n',encoding='ascii')
            out=root/'homology'

            def fake_homology(command,**kwargs):
                output=Path(command[command.index('--outdir')+1]);output.mkdir(parents=True)
                write_tsv(output/'targets_summary.tsv',[{'species':'Schisto','proteins_in_proteome':'1'}])
                write_tsv(output/'homology_audit.tsv',[{'reference_id':'candidate|NEG0000001','species':'Schisto',
                    'status':'RBH_candidate','mapping_method':'identical_sequence','source_assignment':'explicit_source_map',
                    'reference_proteome':'homo_sapiens_negative_controls','best_target_id':'target1',
                    'reverse_best_id':'candidate|NEG0000001','forward_pident':'45.0','forward_qcov':'95.0',
                    'forward_scov':'92.0','reverse_pident':'45.0','reverse_qcov':'95.0','reverse_scov':'92.0'}])
                return type('Result',(),{'returncode':0})()

            query=dict(original_id='human_negative|Q1',species='Homo sapiens',control_class='negative',
                accession='Q1',sequence=query_sequence,sequence_sha256=hashlib.sha256(query_sequence.encode()).hexdigest())
            with patch.object(negative_homology.subprocess,'run',side_effect=fake_homology):
                added,summary=negative_homology.run([query],root,out,threads=2)
            self.assertEqual(len(added),1)
            self.assertEqual(summary['confirmed_one_to_one'],1)
            self.assertEqual(added[0]['negative_homology_source_id'],'human_negative|Q1')
            matches=list(csv.DictReader((out/'homology_matches.tsv').open(encoding='utf-8'),delimiter='\t'))
            self.assertEqual(matches[0]['one_to_one_decision'],'confirmed_unique_reciprocal_best_hit_to_exact_human_source')


if __name__=='__main__':unittest.main()
