"""Evidence joins preserve the stage 1 pool and measured motif context."""

# Resolve source resources independently of the caller working directory.
import sys as _layout_sys
from pathlib import Path as _LayoutPath
_LAYOUT_BUNDLE_ROOT = next(p for p in _LayoutPath(__file__).resolve().parents
                          if (p/'config/layout_paths.json').is_file() and (p/'pipeline_paths.py').is_file())
if str(_LAYOUT_BUNDLE_ROOT) not in _layout_sys.path:
    _layout_sys.path.insert(0, str(_LAYOUT_BUNDLE_ROOT))
from pipeline_paths import activate as _layout_activate, resource_path as _layout_resource
_layout_activate()


import csv
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0,str(_layout_resource(Path(__file__)).resolve().parents[2]))
from audit_join import join


class EvidenceJoinTests(unittest.TestCase):
    def write(self,path,rows):
        path.parent.mkdir(parents=True,exist_ok=True)
        with _layout_resource(path).open('w',newline='') as handle:
            writer=csv.DictWriter(handle,fieldnames=list(rows[0]),delimiter='\t')
            writer.writeheader();writer.writerows(rows)

    def test_hog_references_and_motif_context_without_replacing_hypotheses(self):
        with tempfile.TemporaryDirectory() as td:
            root=_layout_resource(Path(td));ann=_layout_resource(root/'annotation');panel=_layout_resource(root/'panel');ort=_layout_resource(root/'orthology')
            panel.mkdir();(_layout_resource(panel/'targets.json')).write_text(json.dumps(['Target']))
            (_layout_resource(panel/'species.json')).write_text(json.dumps(['Target','RefA','RefB']))
            self.write(_layout_resource(ort/'candidates/candidate_report.tsv'),[dict(candidate_id='Target__p',HOG='H1',HOG_status='shared HOG',mapping_status='mapped')])
            self.write(_layout_resource(ort/'candidates/candidate_HOGs.tsv'),[dict(HOG='H1',Target='Target__p',RefA='RefA__a, RefA__b',RefB='RefB__a')])
            self.write(_layout_resource(ann/'motifs/current_motif_rules.tsv'),[dict(rule_id='CURRENT')])
            self.write(_layout_resource(ann/'motifs/motif_annotations.tsv'),[
                dict(sequence_id='Target__p',motif_id='CURRENT',start=1,end=3,peptide='AAA',candidate_tier='context_supported_candidate'),
                dict(sequence_id='Target__p',motif_id='CURRENT',start=4,end=6,peptide='AAA',candidate_tier='sequence_only_candidate'),
                dict(sequence_id='Target__p',motif_id='CURRENT',start=7,end=9,peptide='AAA',candidate_tier='context_conflict'),
                dict(sequence_id='Target__p',motif_id='OLD',start=10,end=12,peptide='AAA',candidate_tier='context_supported_candidate')])
            # The join accepts both sealed annotations and older handoffs.
            # Include the minimal handoff manifest so either path is covered.
            (_layout_resource(root/'workflow_manifest.json')).write_text(json.dumps({'target_species':['Target']}))
            rows=[dict(protein_id='Target__p',Specie='Target',family='profilin',assigned_family='profilin'),
                  dict(protein_id='Target__p',Specie='Target',family='Nck',assigned_family='')]
            before=json.dumps(rows,sort_keys=True)
            result,support=join(ann,root,panel,ort,rows,{'Orthology_benchmark_weight':.5})
            self.assertEqual(before,json.dumps(rows,sort_keys=True))
            self.assertEqual([r['family'] for r in result],['profilin','Nck'])
            self.assertEqual({r['protein_id'] for r in result},{'Target__p'})
            row=result[0]
            self.assertEqual(row['Orthology_HOG_reference_species_count'],2)
            self.assertEqual(row['Orthology_direct_reference_count'],0)
            self.assertEqual(row['Orthology_RefA_HOG_copy_count'],2)
            self.assertEqual(row['MotifPrediction_context_supported_count'],1)
            self.assertEqual(row['MotifPrediction_sequence_only_count'],1)
            self.assertEqual(row['MotifPrediction_context_conflict_count'],1)
            self.assertEqual(row['MotifPrediction_priority_site_count'],2)
            self.assertEqual(row['MotifPrediction_ignored_superseded_count'],1)
            self.assertEqual(len(support['tables']['Motif_annotations']),3)
            self.assertEqual(support['reference_species'],['RefA','RefB'])


if __name__=='__main__':unittest.main()
