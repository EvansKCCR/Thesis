
# Resolve the organized bundle independently of the caller working directory.
import sys as _bundle_sys
from pathlib import Path as _BundlePath
_BUNDLE_ROOT=next(p for p in _BundlePath(__file__).resolve().parents if (p/'pipeline_paths.py').is_file())
if str(_BUNDLE_ROOT) not in _bundle_sys.path:_bundle_sys.path.insert(0,str(_BUNDLE_ROOT))
from pipeline_paths import activate as _activate_bundle_paths
_activate_bundle_paths()

import sys as _sys
from pathlib import Path as _Path
_PIPELINE_ROOT=_Path(__file__).resolve().parents[2]
_sys.path.insert(0,str(_PIPELINE_ROOT))
_sys.path.insert(0,str(_PIPELINE_ROOT/"annotation"))
import csv
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from control_exclusion import apply
from calibrate_identity import fit


class ControlExclusionTests(unittest.TestCase):
    def bank(self):
        negative={'class':'hard_negative control','accession':'NEG',
            'negative_specific_models':{'spectrin_specific'}}
        return dict(negative_by_hash={'NEG_HASH':[negative]},positive_hashes={'POS_HASH'},
            negative_families={'spectrin_like':[negative]})

    def row(self,**values):
        record=dict(sequence_sha256='NEW_HASH',assigned_family='spectrin_like',
            protein_name='Putative spectrin',identity_rule='configured_family_architecture',
            domain_score_0_4=4,family_specific_CDD_support=True,length=1000,
            NCBI_CDD_architecture='spectrin_specific:20-200')
        record.update(values);return record

    def test_exact_negative_is_excluded_even_when_classifier_abstains(self):
        with patch('control_exclusion.controls',return_value=self.bank()):
            kept,removed=apply([self.row(sequence_sha256='NEG_HASH',assigned_family='unassigned',domain_score_0_4=0)])
        self.assertEqual(kept,[]);self.assertEqual(len(removed),1)
        self.assertEqual(removed[0]['exclusion_gate'],'exact_negative_control_sequence')

    def test_complete_core_and_negative_specific_observed_cdd_excludes(self):
        with patch('control_exclusion.controls',return_value=self.bank()):
            kept,removed=apply([self.row()])
        self.assertEqual(kept,[])
        self.assertEqual(removed[0]['exclusion_gate'],'diagnostic_negative_subfamily_and_complete_architecture')

    def test_positive_controls_shared_domains_partial_cores_and_inferences_are_not_excluded(self):
        rows=[self.row(sequence_sha256='POS_HASH'),self.row(NCBI_CDD_architecture='SPEC:20-200'),
            self.row(domain_score_0_4=2),self.row(NCBI_CDD_architecture='[Inferred: spectrin_specific:20-200]')]
        with patch('control_exclusion.controls',return_value=self.bank()):
            kept,removed=apply(rows)
        self.assertEqual(len(kept),4);self.assertEqual(removed,[])

    def test_final_export_drops_negatives_and_preserves_accounting(self):
        from openpyxl import Workbook,load_workbook
        from integrated_adhesome_pipeline import run_xlsx,ROOT
        rules=json.loads((ROOT/'detection/config/family_rules.json').read_text())
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);source=root/'input.xlsx';output=root/'output.xlsx'
            rows=[dict(Specie='test',protein_id=pid,family='actin',length=400,
                architecture='PF00022:10-390',sequence_sha256=digest)
                for pid,digest in [('positive','POS_HASH'),('negative','NEG_HASH')]]
            book=Workbook();sheet=book.active;sheet.title='candidates';sheet.append(list(rows[0]))
            for row in rows:sheet.append(list(row.values()))
            book.save(source)
            with patch('control_exclusion.controls',return_value=self.bank()):
                result=run_xlsx(source,output,'candidates',rules)
            self.assertEqual(result['input_unique_proteins'],2)
            self.assertEqual(result['unique_proteins'],1);self.assertEqual(result['excluded_negative_controls'],1)
            self.assertTrue(result['protein_accounting_complete'])
            book=load_workbook(output,read_only=True);values=iter(book['Audited_family_assignment'].values)
            headings=next(values);final=[dict(zip(headings,r)) for r in values];book.close()
            self.assertEqual([r['protein_id'] for r in final],['positive'])
            with (root/'negative_control_exclusions.tsv').open(newline='') as handle:
                excluded=list(csv.DictReader(handle,delimiter='\t'))
            self.assertEqual([r['protein_id'] for r in excluded],['negative'])

    def test_identity_calibration_uses_only_positive_controls(self):
        from openpyxl import Workbook
        from identity_probability import ROOT
        # fit records provenance relative to the installed package.
        with tempfile.TemporaryDirectory(dir=ROOT) as td:
            root=Path(td);data=root/'controls';data.mkdir();feature_path=root/'features.tsv'
            labels=[('Positive control','P1','MAAA'),('Positive control','P2','MCCC'),
                ('hard_negative control','N1','MDDD'),('background_negative control','N2','MEEE')]
            book=Workbook();sheet=book.active;sheet.title='control_identity'
            sheet.append(['class','accession','genes','protein_name'])
            for category,accession,seq in labels:sheet.append([category,accession,'TLN1','Talin-1'])
            book.save(data/'controls.xlsx')
            (data/'controls.faa').write_text(''.join(f'>{acc}\n{seq}\n' for cat,acc,seq in labels))
            with feature_path.open('w',newline='') as handle:
                writer=csv.DictWriter(handle,delimiter='\t',fieldnames=['sequence_sha256','assigned_family',
                    'total_score_0_10','family_confidence','protein_name_scope','orthology_panel_sha256'])
                writer.writeheader()
                for cat,acc,seq in labels:writer.writerow(dict(sequence_sha256=hashlib.sha256(seq.encode()).hexdigest(),
                    assigned_family='talin',total_score_0_10=7,family_confidence='strong',protein_name_scope='family',
                    orthology_panel_sha256='panel'))
            result=fit(data,feature_path,root/'model.json',root/'review',minimum=2)
            model=json.loads((root/'model.json').read_text())
            self.assertEqual({s['accession'] for s in model['samples']},{'P1','P2'})
            self.assertTrue(result['context']['positive_control_only'])
            self.assertEqual(result['context']['excluded']['negative_exclusion_control_not_identity_training'],2)

    def test_an_all_negative_final_table_has_no_placeholder_protein(self):
        from audit_report import prepare_sheet
        records,fields,removed=prepare_sheet('Audited_family_assignment',[])
        self.assertEqual(records,[])
        self.assertIn('protein_id',fields);self.assertIn('putative_audited_family',fields)


if __name__=='__main__':unittest.main()
