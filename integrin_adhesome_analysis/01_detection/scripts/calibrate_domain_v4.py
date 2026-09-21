"""Check requested architectures against positive references and shared-domain competitors."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import argparse,json
from collections import defaultdict,Counter
from adhesome_pipeline import scan
from adhesome_classification import architecture_evidence
from adhesome_core import table,write_table,digest,save_json

FAMILIES=('fibronectin_like','vitronectin_like','ILK','integrin_alpha','FAK')

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reuse-scan',action='store_true')
    args=parser.parse_args()
    root=Path(__file__).resolve().parents[1]
    out=root/'results/reference_calibration_v4'
    refs=root/'references/schistosome_v4'
    library=root/'databases/pfam/adhesome_v4.hmm'
    if args.reuse_scan:
        meta=json.loads((out/'domains/scan.json').read_text())
        if meta['proteome_sha256']!=digest(refs/'metazoan_references.faa') or meta['metadata_sha256']!=digest(str(library)+'.metadata.json'):
            raise ValueError('Calibration input changed; rescan required')
    else:
        scan(argparse.Namespace(out=str(out/'domains'),library=str(library),
            manifest=str(root/'config/adhesome_pfam_manifest.csv'),proteome=str(refs/'metazoan_references.faa'),
            cpu=2,max_i=0.001,max_c=0.001,min_coverage=0.35))
    rules=json.loads((root/'config/family_rules.json').read_text())['families']
    hits=defaultdict(list)
    for row in table(out/'domains/domains.csv'):
        hits[row['query_name']].append(row)
    rows=[]
    for reference in table(refs/'metazoan_catalogue.csv'):
        rid=reference['reference_id']
        domains=sorted(hits[rid],key=lambda row:int(row['ali_from']))
        for family in FAMILIES:
            arch=architecture_evidence(domains,rules[family])
            if reference['family']!=family and not arch['partial']:
                continue
            rows.append(dict(reference_id=rid,accession=reference['accession'],annotated_family=reference['family'],
                tested_family=family,expected_positive=reference['family']==family,
                complete=arch['complete'],specific=arch['specific'],full_order=arch['order'],support_order=arch['support_order'],
                observed=arch['observed'],missing=arch['missing']))
    write_table(out/'architecture_calibration.csv',rows,list(rows[0]))
    counts={family:dict(positive_references=sum(r['expected_positive'] for r in rows if r['tested_family']==family),
        positive_complete=sum(r['expected_positive'] and r['complete'] for r in rows if r['tested_family']==family),
        competitor_complete=sum(not r['expected_positive'] and r['complete'] for r in rows if r['tested_family']==family)) for family in FAMILIES}
    save_json(out/'calibration.json',dict(counts=counts,rules_sha256=digest(root/'config/family_rules.json'),
        library_sha256=digest(library),references_sha256=digest(refs/'metazoan_references.faa'),
        interpretation='Same-panel architecture calibration, not an independent accuracy or orthology benchmark; complete architecture alone does not assign a family.'))
    print(json.dumps(counts,indent=2),flush=True)

if __name__=='__main__':
    main()
