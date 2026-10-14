import json,tempfile,unittest
from pathlib import Path
from integrated_adhesome_pipeline import audit_record,run_xlsx,checksum

ROOT=Path(__file__).resolve().parents[1]


class AuditTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rules=json.loads((ROOT/'01_detection/config/family_rules.json').read_text())

    def row(self,**extra):
        return dict(Specie='AnyTarget',protein_id='AnyTarget__g',family='actin',length=380,
                    architecture='PF00022:10-370',**extra)

    def test_domain_fold_and_hog_do_not_resolve_actin_paralog(self):
        result=audit_record(self.row(Orthology_direct_reference_count=5),self.rules)
        self.assertEqual(result['classification_state'],'broad_only')
        self.assertEqual(result['assigned_family'],'unassigned')
        self.assertEqual(result['broad_homology_group'],'actin_superfamily')
        self.assertEqual(result['orthology_score_0_3'],1)

    def test_independent_reference_resolves_conventional_actin_and_arp(self):
        for label,family in [('ACTB','actin_conventional'),('ARP6','ARP6')]:
            with self.subTest(label=label):
                result=audit_record(self.row(reference_subfamily=label,subfamily_support_method='curated_reference',reference_accession='FIXTURE'),self.rules)
                self.assertEqual(result['assigned_family'],family)
                self.assertEqual(result['grade'],'A')

    def test_reference_string_and_unconcordant_phylogeny_do_not_resolve(self):
        for extra in [dict(reference_subfamily='ACTB'),dict(reference_subfamily='ACTB',subfamily_support_method='phylogeny',phylogeny_concordant='no')]:
            self.assertEqual(audit_record(self.row(**extra),self.rules)['classification_state'],'broad_only')

    def test_motif_or_reference_cannot_rescue_missing_architecture(self):
        row=self.row(reference_subfamily='ACTB',subfamily_support_method='curated_reference',reference_accession='FIXTURE',MotifPrediction_context_supported_count=9)
        row['architecture']=''
        result=audit_record(row,self.rules)
        self.assertEqual(result['grade'],'D');self.assertEqual(result['assigned_family'],'unassigned')

    def test_new_workbook_schema_and_source_preservation(self):
        from openpyxl import Workbook,load_workbook
        with tempfile.TemporaryDirectory() as td:
            source,dest=Path(td)/'source.xlsx',Path(td)/'result.xlsx'
            book=Workbook();sheet=book.active;sheet.title='candidates'
            row=self.row(reference_subfamily='ACTB',subfamily_support_method='curated_reference',reference_accession='FIXTURE')
            sheet.append(list(row));sheet.append(list(row.values()));sheet.append(list(row.values()))
            book.create_sheet('Notes').append(['Preserve this source worksheet']);book.save(source)
            source_hash=checksum(source)
            meta=run_xlsx(source,dest,'candidates',self.rules)
            result=load_workbook(dest,read_only=True)
            self.assertEqual(result['Resolved_candidates'].max_row,2)
            self.assertEqual(result['All_hypotheses_audited'].max_row,3)
            self.assertEqual(result['Notes']['A1'].value,'Preserve this source worksheet')
            result.close();self.assertEqual(checksum(source),source_hash)
            self.assertEqual(meta['pipeline_version'],'1.0')
