import csv
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import evaluate_validation as ev


class ValidationEvaluationTests(unittest.TestCase):
    def write(self,path,rows,delimiter='\t'):
        path.parent.mkdir(parents=True,exist_ok=True)
        with path.open('w',newline='',encoding='utf-8') as handle:
            writer=csv.DictWriter(handle,fieldnames=list(rows[0]),delimiter=delimiter);writer.writeheader();writer.writerows(rows)

    def fixture(self,root):
        from openpyxl import Workbook
        detection=root/'detection';unit=detection/'validation';prepared=unit/'prepared';prepared.mkdir(parents=True)
        cohort=[('P1','MAAA'),('P2','MCCC'),('P3','MDDD'),('P4','MEEE'),('N1','MFFF'),('N2','MGGG'),
            ('U1','MHHH'),('P5','MIII'),('N3','MKKK')]
        raw=root/'validation.faa';raw.write_text(''.join(f'>candidate|{key}|GENE\n{seq}\n' for key,seq in cohort))
        fasta=prepared/'validation.faa';fasta.write_text(''.join(f'>validation__{key}\n{seq}\n' for key,seq in cohort))
        (unit/'proteome.json').write_text(json.dumps(dict(proteomes=[dict(code='validation',species='Pipeline validation',path=str(raw))])))
        (prepared/'preparation.json').write_text(json.dumps(dict(species={'validation':dict(input_sha256=ev.sha(raw),output_sha256=ev.sha(fasta))})))
        mapping=prepared/'validation.id_map.csv'
        self.write(mapping,[dict(original_id=f'candidate|{key}|GENE',protein_id='validation__'+key,decision='retained') for key,seq in cohort],',')
        labels=root/'labels.xlsx';book=Workbook();sheet=book.active;sheet.title='control_identity'
        sheet.append(['class','accession','genes','protein_name','trusted_family'])
        for key,seq in cohort:
            if key=='U1':continue
            category='Positive control' if key.startswith('P') else 'hard_negative control' if key=='N1' else 'background_negative control'
            sheet.append([category,key,'','Talin' if key in ('P1','P2') else 'Actin','talin' if key in ('P1','P2') else 'actin'])
        book.save(labels)
        trusted=root/'trusted.faa'
        trusted.write_text(''.join(f'>{key}\n'+('MWWW' if key=='N3' else seq)+'\n' for key,seq in cohort if key not in ('U1','P5')))
        selected={'P1':'talin','P2':'Syk','P4':'actin','N1':'talin','U1':'talin'}
        handoff=root/'handoff';handoff.mkdir()
        self.write(handoff/'hypotheses.tsv',[dict(Specie='validation',protein_id='validation__'+key,family=family) for key,family in selected.items()])
        source=unit/'candidates/candidates.csv';self.write(source,[dict(protein_id='validation__'+key,status='provisional') for key in selected],',')
        (handoff/'workflow_manifest.json').write_text(json.dumps(dict(kind='handoff',detection=str(detection),target_species=['validation'],
            sources={str(p.resolve()):ev.sha(p) for p in (fasta,source)})))
        audited=root/'audit';audited.mkdir();screen=[];final=[]
        for key,family in selected.items():
            grade='C' if key=='P4' else 'A';score=2 if grade=='C' else 4
            row=dict(Specie='validation',protein_id='validation__'+key,family=family,preaudit_family=family,audit_family=family,
                assigned_family='',assignment_status='provisional' if grade=='C' else 'assigned',grade=grade,
                family_decision='provisional_family' if grade=='C' else 'family_supported',domain_score_0_4=score,
                family_specific_CDD_support=True,orthology_score_0_3=2,motif_score_0_2=0,localization_score_0_1=0,
                total_score_0_10=score+2,row_origin='original_stage1',status='provisional',
                family_confidence='provisional' if grade=='C' else 'strong',protein_name='Putative '+family,
                protein_name_scope='family',naming_evidence='fixture computational evidence')
            screen.append(row)
            if key!='N1':final.append(dict(Specie='validation',protein_id='validation__'+key,putative_audited_family=family,
                assignment_confidence=row['family_confidence'],domain_score_0_4=score,grade=grade,
                protein_name=row['protein_name'],protein_name_scope='family'))
        book=Workbook();sheet=book.active;sheet.title='Screening_hypothesis';sheet.append(list(screen[0]))
        for row in screen:sheet.append(list(row.values()))
        book.save(audited/'Stage_2_output.xlsx')
        self.write(audited/'audited_family_assignments.tsv',final)
        self.write(audited/'negative_control_exclusions.tsv',[dict(Specie='validation',protein_id='validation__N1',exclusion_gate='exact_negative_control_sequence')])
        return dict(audit_dir=audited,handoff=handoff,out=root/'report',labels=labels,control_fasta=trusted)

    def test_binary_metrics_and_undefined_denominators(self):
        result=ev.binary_metrics(3,1,5,2)
        self.assertEqual(result['PPV'],.75);self.assertAlmostEqual(result['F1'],6/9)
        self.assertEqual(result['n'],11)
        empty=ev.binary_metrics(0,0,4,0)
        self.assertIsNone(empty['PPV']);self.assertIsNone(empty['recall']);self.assertIsNone(empty['F1'])
        self.assertEqual(empty['specificity'],1)

    def test_end_to_end_denominators_missed_positives_and_negative_gate_diagnostics(self):
        with tempfile.TemporaryDirectory() as td:
            args=self.fixture(Path(td));before=ev.sha(args['audit_dir']/'audited_family_assignments.tsv')
            result=ev.evaluate(**args)
            self.assertEqual(result['status'],'complete');self.assertEqual(result['scored_controls'],7)
            metrics=ev.read_table(args['out']/'metrics.tsv')
            get=lambda policy:next(r for r in metrics if r['policy']==policy and r['population']=='all_controls')
            main=get('dual_gate_supported');self.assertEqual(tuple(int(main[k]) for k in ('TP','FP','TN','FN')),(2,0,2,3))
            self.assertAlmostEqual(float(main['F1']),4/7);self.assertEqual(float(main['PPV']),1)
            prior=get('audit_before_negative_exclusion');self.assertEqual(int(prior['FP']),1)
            self.assertAlmostEqual(float(prior['PPV']),2/3)
            stage1=get('stage1_recruitment');self.assertEqual(tuple(int(stage1[k]) for k in ('TP','FP','TN','FN')),(3,1,1,2))
            self.assertAlmostEqual(float(stage1['F1']),2/3)
            family=next(r for r in ev.read_table(args['out']/'family_metrics.tsv') if r['policy']=='dual_gate_supported')
            self.assertEqual(tuple(int(family[k]) for k in ('TP','FP','FN')),(1,1,4))
            self.assertAlmostEqual(float(family['F1']),2/7)
            predictions={r['accession']:r for r in ev.read_table(args['out']/'control_predictions.tsv')}
            self.assertEqual(predictions['P3']['stage1_recruitment'],'False')
            self.assertEqual(predictions['U1']['evaluation_status'],'unlabelled_control')
            self.assertEqual(predictions['P5']['evaluation_status'],'scored')
            self.assertEqual(predictions['P5']['label_sequence_verification'],'accession_label_only')
            self.assertEqual(predictions['N3']['evaluation_status'],'labelled_sequence_mismatch')
            summary=json.loads((args['out']/'summary.json').read_text())
            self.assertEqual(summary['input_sequences'],9);self.assertEqual(summary['negative_control_exclusions'],1)
            self.assertEqual(summary['sequence_verified_controls'],6);self.assertEqual(summary['accession_label_only_controls'],1)
            verified=next(r for r in metrics if r['policy']=='dual_gate_supported' and r['population']=='sequence_verified_controls')
            self.assertEqual(int(verified['n']),6);self.assertAlmostEqual(float(verified['F1']),2/3)
            self.assertIn('not independent',summary['interpretation'])
            self.assertEqual(ev.sha(args['audit_dir']/'audited_family_assignments.tsv'),before)
            self.assertTrue((args['out']/'performance_report.html').is_file())

    def test_changed_raw_input_cannot_be_scored_as_the_original_control_cohort(self):
        with tempfile.TemporaryDirectory() as td:
            args=self.fixture(Path(td));raw=Path(td)/'validation.faa';raw.write_text(raw.read_text().replace('MAAA','MAAC'))
            with self.assertRaisesRegex(ValueError,'changed after detection'):ev.evaluate(**args)
            self.assertFalse(args['out'].exists())

    def test_positive_lost_during_preparation_stays_in_the_metric_denominator(self):
        from openpyxl import load_workbook
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);args=self.fixture(root);unit=root/'detection/validation';prepared=unit/'prepared/validation.faa'
            mapping=unit/'prepared/validation.id_map.csv';rows=ev.read_table(mapping,',')
            for row in rows:
                if row['protein_id']=='validation__P1':row['decision']='below_global_min_length'
            self.write(mapping,rows,',')
            prepared.write_text(prepared.read_text().replace('>validation__P1\nMAAA\n',''))
            receipt=unit/'prepared/preparation.json';meta=json.loads(receipt.read_text());meta['species']['validation']['output_sha256']=ev.sha(prepared);receipt.write_text(json.dumps(meta))
            paths=[args['handoff']/'hypotheses.tsv',args['audit_dir']/'audited_family_assignments.tsv']
            for path in paths:self.write(path,[r for r in ev.read_table(path) if r['protein_id']!='validation__P1'])
            source=unit/'candidates/candidates.csv';self.write(source,[r for r in ev.read_table(source,',') if r['protein_id']!='validation__P1'],',')
            manifest=args['handoff']/'workflow_manifest.json';meta=json.loads(manifest.read_text())
            for path in (prepared,source):meta['sources'][str(path.resolve())]=ev.sha(path)
            manifest.write_text(json.dumps(meta))
            book=load_workbook(args['audit_dir']/'Stage_2_output.xlsx');sheet=book['Screening_hypothesis']
            for index in range(sheet.max_row,1,-1):
                if sheet.cell(index,2).value=='validation__P1':sheet.delete_rows(index)
            book.save(args['audit_dir']/'Stage_2_output.xlsx');book.close()
            ev.evaluate(**args)
            main=next(r for r in ev.read_table(args['out']/'metrics.tsv') if r['policy']=='dual_gate_supported' and r['population']=='all_controls')
            self.assertEqual(int(main['n']),7);self.assertEqual(int(main['FN']),4)

    def test_an_excluded_control_cannot_reappear_in_final_assignments(self):
        with tempfile.TemporaryDirectory() as td:
            args=self.fixture(Path(td));path=args['audit_dir']/'audited_family_assignments.tsv';rows=ev.read_table(path)
            rows.append(dict(rows[0],protein_id='validation__N1'));self.write(path,rows)
            with self.assertRaisesRegex(ValueError,'remains in final assignments'):ev.evaluate(**args)

    def test_ambiguous_and_provisional_calls_cannot_pass_the_full_dual_gate(self):
        self.assertFalse(ev.supported(dict(putative_audited_family='talin',domain_score_0_4=4,grade='A',assignment_confidence='ambiguous')))
        self.assertFalse(ev.supported(dict(putative_audited_family='talin',domain_score_0_4=2,grade='C',assignment_confidence='provisional')))
        row=dict(putative_audited_family='talin',domain_score_0_4=4,grade='B',assignment_confidence='moderate')
        self.assertTrue(ev.supported(row));self.assertFalse(ev.supported(row,strong=True))

    def test_no_validation_target_is_an_explicit_non_applicable_status(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);handoff=root/'handoff';handoff.mkdir();(handoff/'workflow_manifest.json').write_text('{"target_species":["TargetA"]}')
            result=ev.evaluate(root/'audit',handoff,root/'report')
            self.assertEqual(result['status'],'validation_target_not_present');self.assertFalse((root/'report').exists())

    def test_evaluation_failure_preserves_audit_and_reports_the_error(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);audit=root/'audit';audit.mkdir();sentinel=audit/'calls.tsv';sentinel.write_text('original calls')
            with patch('evaluate_validation.evaluate',side_effect=ValueError('invalid control provenance')):
                result=ev.after_audit(audit,root/'handoff')
            self.assertEqual(result['status'],'failed');self.assertEqual(sentinel.read_text(),'original calls')
            self.assertIn('invalid control provenance',(audit/'validation_evaluation/error.log').read_text())

    def test_existing_report_is_not_overwritten(self):
        with tempfile.TemporaryDirectory() as td:
            args=self.fixture(Path(td));args['out'].mkdir();(args['out']/'summary.json').write_text('original')
            with self.assertRaises(FileExistsError):ev.evaluate(**args)
            self.assertEqual((args['out']/'summary.json').read_text(),'original')


if __name__=='__main__':unittest.main()
