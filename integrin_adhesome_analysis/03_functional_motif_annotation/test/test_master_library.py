import json
import unittest
from unittest.mock import patch
import master_library as ml
import family_motifs as fm


class MasterTests(unittest.TestCase):
    def test_supplied_master_counts_and_exclusions_gap(self):
        rules=json.loads((fm.HERE/'Motif_library/context_rules.json').read_text())
        curated,sheets,_,generic,audit,meta=ml.load(fm.HERE/'Motif_library/exhaustive_comprehensive_adhesion_motif_library_PINCH_Parvin_updated.xlsx',rules)
        self.assertEqual(meta['master_records'],691)
        self.assertEqual(len(curated),37)
        self.assertEqual(sum(r['enabled'] for r in generic),352)
        self.assertFalse(meta['qa_sheet_present'])
        self.assertEqual(len(meta['summary_discrepancies']),1)
        entries={r['motif_id']:r for r in curated}
        self.assertTrue(entries['SLIM037']['reader_only'])
        for mid in ('SLIM038','SLIM039'):
            self.assertEqual(rules['motifs'][mid]['target_families'],['paxillin_like'])
            self.assertEqual(entries[mid]['screening_recommendation'],'Review only')
        self.assertTrue(all(not r['enabled'] for r in generic if r['source_collection']=='Legacy ELM candidate library'))

    def test_flags_qa_invalid_and_reader(self):
        header=['library_id','source_collection','family','scan_enabled','regex_valid','regex_python','screening_recommendation','motif_type']
        rows=[['a','ELM class catalogue','ELM','Yes','Yes','RGD','Include with ELM context filters','LIG'],
              ['b','Legacy ELM candidate library','Legacy candidate','No','Yes','RGD','Review only','LIG'],
              ['c','ELM class catalogue','ELM','Yes','Yes','[','Include with ELM context filters','LIG'],
              ['d','ELM class catalogue','ELM','Yes','Yes','RGD','Annotation only','motif-reader annotation']]
        sheets={'Master_motif_library':[{'row':1,'values':header}]+[{'row':i+2,'values':r} for i,r in enumerate(rows)],
                'QA_exclusions':[{'row':1,'values':['library_id']},{'row':2,'values':['a']}]}
        with patch('master_library.ca.read_workbook',return_value=sheets):
            _,_,_,generic,audit,meta=ml.load('unused',{'motifs':{}})
        self.assertTrue(all(not r['enabled'] for r in generic))
        self.assertEqual(audit[0]['effective_state'],'excluded_by_QA')
        self.assertEqual(audit[2]['issue'],'invalid_or_ambiguous_regex')
        self.assertEqual(audit[3]['effective_state'],'reader_annotation_only')


if __name__=='__main__':unittest.main()
