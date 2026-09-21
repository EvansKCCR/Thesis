"""Small regression fixtures; no external software or real results required."""
import argparse
import gzip
import tempfile
import unittest
from pathlib import Path
import pipeline as p
from refseq_gene_map import parse_gff


class PipelineTests(unittest.TestCase):
    def test_refseq_mapping_repeated_cds_and_conflicting_genes(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'annotation.gff'
            def cds(gene):
                return 'chr1\tRefSeq\tCDS\t1\t9\t.\t+\t0\tprotein_id=NP_123.2;Dbxref=GeneID:'+gene+'\n'
            path.write_text('#!genome-build test\n'+cds('123')+cds('123'))
            mapping, metadata = parse_gff(path)
            self.assertEqual(mapping['NP_123.2'], {'123'})
            self.assertEqual(metadata, ['#!genome-build test'])
            path.write_text(cds('123')+cds('456'))
            mapping, _ = parse_gff(path)
            self.assertEqual(mapping['NP_123.2'], {'123','456'})

    def test_compressed_fasta_and_candidate_header(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            f = root/'input.fa'
            with gzip.open(f,'wt') as out:
                out.write('>x gene=g\nMAAA\n')
            self.assertEqual(list(p.records(f)), [('x gene=g','MAAA')])
            f.write_text('protein_id\nShae__g\nShae__g\n')
            self.assertEqual(p.candidates(f), (['Shae__g'],2))

    def test_pairwise_direction_and_hog_membership(self):
        previous = p.ROOT
        try:
            with tempfile.TemporaryDirectory() as tmp:
                p.ROOT = Path(tmp)
                result = p.ROOT/'results'
                p.write_tsv(p.ROOT/'03_qc/candidate_audit.tsv',
                    ['candidate_set','candidate_id','prepared_id','status'],
                    [dict(candidate_set='test',candidate_id=c,prepared_id=c,status='matched')
                     for c in ['Shae__a','Shae__b','Shae__absent']])
                fields = ['HOG','OG',*p.SPECIES]
                p.write_tsv(result/'Phylogenetic_Hierarchical_Orthogroups/N0.tsv',fields,
                    [dict(HOG='N0.HOG1',OG='OG1',Hsap='Hsap__h',Cele='Cele__c',
                          Shae='Shae__a, Shae__b',Sjap='',Sman='')])
                # Reverse column order: multiplicity must remain candidate -> human.
                p.write_tsv(result/'Orthologues/Hsap/Hsap__v__Shae.tsv', ['HOG','Hsap','Shae'],
                    [dict(HOG='N0.HOG1',Hsap='Hsap__h',Shae='Shae__a, Shae__b')])
                out = p.ROOT/'filtered'
                p.extract(argparse.Namespace(results=result,hog='N0.tsv',candidate_set='test',output=out))
                _, rows = p.read_tsv(out/'candidate_report.tsv')
                self.assertEqual(rows[0]['Hsap_relationship'],'many-to-one')
                self.assertEqual(rows[0]['Cele_relationship'],'orthogroup member without direct pairwise orthologue')
                self.assertEqual(rows[2]['HOG_status'],'not found in selected HOG table')
                # Modern root table and combined tables, including duplicate pairwise records.
                p.write_tsv(result/'Orthogroups/Orthogroups.tsv', ['Orthogroup',*p.SPECIES],
                    [dict(Orthogroup='OG1',Hsap='Hsap__h',Cele='Cele__c',Shae='Shae__a, Shae__b',Sjap='',Sman='')])
                (result/'Phylogenetic_Hierarchical_Orthogroups/N0.tsv').unlink()
                (result/'Log.txt').write_text('Started OrthoFinder version 3.1.5\n')
                p.write_tsv(result/'Orthologues/Shae.tsv', ['Orthogroup','Species','Shae','Orthologs'],
                    [dict(Orthogroup='OG1',Species='Hsap',Shae='Shae__a, Shae__b',Orthologs='Hsap__h')])
                modern = p.ROOT/'modern'
                p.extract(argparse.Namespace(results=result,hog='auto',candidate_set='test',output=modern))
                _, rows = p.read_tsv(modern/'candidate_report.tsv')
                self.assertEqual(rows[0]['HOG'],'OG1')
                self.assertEqual(rows[0]['Hsap_relationship'],'many-to-one')
                _, pairs = p.read_tsv(modern/'direct_orthologues.tsv')
                self.assertEqual(len(pairs),2)
                (result/'Log.txt').write_text('Started OrthoFinder version 2.5.5\n')
                with self.assertRaises(ValueError):
                    p.resolve_group_table(result,'auto')
        finally:
            p.ROOT = previous


if __name__ == '__main__':
    unittest.main()
