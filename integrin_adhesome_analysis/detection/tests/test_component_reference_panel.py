
# Resolve the organized bundle independently of the caller working directory.
import sys as _bundle_sys
from pathlib import Path as _BundlePath
_BUNDLE_ROOT=next(p for p in _BundlePath(__file__).resolve().parents if (p/'pipeline_paths.py').is_file())
if str(_BUNDLE_ROOT) not in _bundle_sys.path:_bundle_sys.path.insert(0,str(_BUNDLE_ROOT))
from pipeline_paths import activate as _activate_bundle_paths
_activate_bundle_paths()
import argparse
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch, Mock

from adhesome_core import digest, table, write_table, write_fasta, sequence_dict
from adhesome_pipeline import catalogue, reciprocal_agreement
from expand_component_references import component_identity, resolve_component, expand
from build_reference_sets import UniProtClient, build as build_base_panel


def entry(accession, gene, taxon=9606):
    return dict(primaryAccession=accession, genes=[dict(geneName=dict(value=gene))],
        organism=dict(taxonId=taxon,scientificName='Homo sapiens'),
        entryType='UniProtKB reviewed (Swiss-Prot)',sequence=dict(value='A'*100),
        entryAudit=dict(sequenceVersion=1,entryVersion=1))


class ComponentPanelTests(unittest.TestCase):
    def test_mixed_pins_require_matching_component_and_taxon(self):
        component={'Official Symbol':'ONE','Gene ID':'1','reference_accession':'P00001, P00002'}
        entries={'P00001':entry('P00001','ONE'),'P00002':entry('P00002','OTHER')}
        resolved,seq,rejected=resolve_component(component,entries)
        self.assertEqual(resolved['primaryAccession'],'P00001')
        self.assertEqual(rejected[0]['accession'],'P00002')
        entries['P00001']['organism']['taxonId']=10090
        with self.assertRaisesRegex(ValueError,'expected one matching'):
            resolve_component(component,entries)

    def test_geneid_allows_verified_renamed_symbol_but_not_arbitrary_accession(self):
        component={'Official Symbol':'OLD','Gene ID':'1','reference_accession':'P00001'}
        value=entry('P00001','NEW')
        self.assertFalse(component_identity(value,component))
        value['uniProtKBCrossReferences']=[dict(database='GeneID',id='1')]
        self.assertTrue(component_identity(value,component))

    def test_inactive_pin_can_only_be_repaired_by_unique_gene_taxon_match(self):
        component={'Official Symbol':'ONE','Gene ID':'1','reference_accession':'P00001'}
        entries={'P00001':dict(primaryAccession='P00001',entryType='Inactive')}
        client=Mock()
        client.get.side_effect=[dict(results=[entry('P00002','ONE')]),entry('P00002','ONE')]
        resolved,seq,rejected=resolve_component(component,entries,client,True)
        self.assertEqual(resolved['primaryAccession'],'P00002')
        self.assertEqual(rejected[0]['accession'],'P00001')
        client.get.side_effect=[dict(results=[entry('P00002','ONE'),entry('P00003','ONE')])]
        with self.assertRaisesRegex(ValueError,'Ambiguous'):
            resolve_component(component,entries,client,True)

    def test_hierarchical_duplicate_references_do_not_break_reciprocity(self):
        rows=[dict(reference_id='family',family='broad',accession='P1'),
              dict(reference_id='subfamily',family='specific',accession='P1',component_symbol='ONE'),
              dict(reference_id='component',family='component_reference',accession='P1',component_symbol='ONE')]
        self.assertTrue(reciprocal_agreement(dict(family='broad'),rows))
        self.assertTrue(reciprocal_agreement(dict(family='specific'),rows))
        self.assertTrue(reciprocal_agreement(dict(family='component_reference',component_symbol='ONE'),rows))
        rows.append(dict(reference_id='other',family='component_reference',accession='P2',component_symbol='TWO'))
        self.assertFalse(reciprocal_agreement(dict(family='broad'),rows))
        self.assertFalse(reciprocal_agreement(dict(family='component_reference',component_symbol='ONE'),rows))

    def test_complete_paired_expansion_and_rebuild_do_not_accumulate_duplicates(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);base=root/'base';base.mkdir()
            fields=['reference_id','family','module','species','taxon_id','gene','accession','evidence_source','experimental_evidence']
            original=dict(zip(fields,['old','broad','actin_linking','Homo sapiens','9606','ONE','P00001','fixture','seed claim']))
            competitor=dict(original,reference_id='negative',family='non_adhesome',accession='P00003',gene='NEG')
            for label,rows in [('search',[original]),('metazoan',[original,competitor])]:
                write_table(base/(label+'_catalogue.csv'),rows,fields)
                write_fasta(base/(label+'_references.faa'),[(r['reference_id'],'A'*100) for r in rows])
            (base/'reference_lock.json').write_text(json.dumps(dict(files={p.name:digest(p) for p in base.iterdir()})))
            crosswalk=root/'crosswalk.csv'
            cf=['Official Symbol','Gene ID','Protein name','Swiss-Prot ID','Functional Category','FA','reference_accession','reference_pfams']
            write_table(crosswalk,[dict(zip(cf,[gene,str(i),gene,accession,'Adaptor','fixture',accession,'PF00001']))
                for i,(gene,accession) in enumerate([('ONE','P00001'),('TWO','P00002')],1)],cf)
            rules=root/'rules.json'
            rules.write_text(json.dumps(dict(families=dict(broad=dict(module='actin_linking'),
                specific=dict(module='actin_linking',reference_labels=['ONE'])))))
            args=argparse.Namespace(base=base,out=root/'expanded',crosswalk=crosswalk,rules=rules,cache=root/'cache',workers=2,offline=True)
            def get(client,url):
                accession=url.rsplit('/',1)[-1].split('.')[0]
                return entry(accession,dict(P00001='ONE',P00002='TWO')[accession])
            with patch.object(UniProtClient,'get',get):
                expand(args)
            self.assertEqual(len(table(args.out/'component_coverage.csv')),2)
            self.assertEqual(len(sequence_dict(args.out/'search_references.faa')),5)
            self.assertEqual(len(sequence_dict(args.out/'metazoan_references.faa')),6)
            rows=catalogue(args.out/'search_catalogue.csv',args.out/'search_references.faa',{'broad':{},'specific':{}})
            self.assertIn('Hsa_component_TWO',rows)
            self.assertEqual(rows['Hsa_component_TWO']['family'],'component_reference')
            self.assertEqual(rows['Hsa_component_ONE__specific']['family'],'specific')
            self.assertEqual(rows['Hsa_component_ONE']['adhesome_experiment_verified'],'false')
            lock=json.loads((args.out/'reference_lock.json').read_text())
            self.assertTrue(all(digest(args.out/name)==sha for name,sha in lock['files'].items()))
            with self.assertRaisesRegex(ValueError,'expanded crosswalk panel'):
                build_base_panel(argparse.Namespace(out=args.out))
            old={r['reference_id'] for r in table(args.out/'search_catalogue.csv')}
            args.base=args.out;args.out=root/'rebuilt'
            with patch.object(UniProtClient,'get',get):
                expand(args)
            self.assertEqual(old,{r['reference_id'] for r in table(args.out/'search_catalogue.csv')})


if __name__=='__main__':
    unittest.main()
