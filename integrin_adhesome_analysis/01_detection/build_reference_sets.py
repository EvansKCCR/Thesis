#!/usr/bin/env python3
"""Build traceable paired FASTA/catalogues from the supplied adhesome seed bundle."""
import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import time
import urllib.parse
import urllib.request
from adhesome_core import table, write_table, write_fasta, save_json, digest, sequence_dict

ROOT = Path(__file__).resolve().parent
DEFAULT_SEEDS = ROOT / 'references' / 'seeds' / 'schistosome_v5'

class UniProtClient:
    def __init__(self, cache, offline=False):
        self.cache = Path(cache)
        self.cache.mkdir(parents=True, exist_ok=True)
        self.offline = offline
        self.requests = []

    def get(self, url):
        path = self.cache / (hashlib.sha256(url.encode()).hexdigest() + '.json')
        if path.exists():
            record = json.loads(path.read_text(encoding='utf-8'))
            if record['url'] != url:
                raise ValueError('Cache URL mismatch')
        elif self.offline:
            raise ValueError(f'Offline cache missing {url}')
        else:
            for attempt in range(3):
                try:
                    request = urllib.request.Request(url, headers={'User-Agent':'adhesome-reference-builder/2.0','Accept':'application/json'})
                    with urllib.request.urlopen(request, timeout=40) as response:
                        record = dict(url=url, retrieved_at=datetime.now(timezone.utc).isoformat(),
                            uniprot_release=response.headers.get('X-UniProt-Release'),
                            release_date=response.headers.get('X-UniProt-Release-Date'),data=json.load(response))
                    save_json(path,record)
                    break
                except Exception:
                    if attempt == 2:
                        raise
                    time.sleep(attempt + 1)
        self.requests.append(dict(url=url,cache_file=str(path),sha256=digest(path)))
        return record['data']

def gene_names(entry):
    names = []
    for gene in entry.get('genes',[]):
        if gene.get('geneName'):
            names.append(gene['geneName']['value'])
        for key in ('synonyms','orderedLocusNames','orfNames'):
            names += [n['value'] for n in gene.get(key,[])]
    return names

def valid_identity(entry, row):
    if row.get('identity_mode') == 'accession_taxon' and row.get('role') == 'phylogenetic_anchor':
        return (str(entry.get('organism',{}).get('taxonId')) == row['taxon_id'] and
                row.get('accession') in [entry.get('primaryAccession'),*entry.get('secondaryAccessions',[])])
    return (str(entry.get('organism',{}).get('taxonId')) == row['taxon_id'] and
            row['gene'].casefold() in {g.casefold() for g in gene_names(entry)})

def resolve(row, client, repair=False, ambiguity_policy='error'):
    requested = row.get('accession','').strip()
    if requested:
        entry = client.get('https://rest.uniprot.org/uniprotkb/' + urllib.parse.quote(requested) + '.json')
        if valid_identity(entry,row):
            return entry, 'accession', []
        if not repair:
            raise ValueError(f'Identity mismatch: {requested} -> {gene_names(entry)} taxon {entry.get("organism",{}).get("taxonId")}')
    query = f'(gene_exact:{row["gene"]}) AND (organism_id:{row["taxon_id"]})'
    url = 'https://rest.uniprot.org/uniprotkb/search?' + urllib.parse.urlencode(dict(query=query,format='json',size=500))
    results = client.get(url).get('results',[])
    if len(results) >= 500:
        raise ValueError('Gene query potentially truncated; pin an accession')
    exact = [e for e in results if valid_identity(e,row)]
    reviewed = [e for e in exact if 'Swiss-Prot' in e.get('entryType','')]
    choices = reviewed or exact
    if len(choices) > 1 and ambiguity_policy == 'prefer-complete':
        def rank(entry):
            existence = entry.get('proteinExistence','5: Uncertain')
            match = re.match(r'(\d+):',existence)
            pe = int(match.group(1)) if match else 5
            return (entry.get('proteinDescription',{}).get('flag') == 'Fragment',pe,
                    -entry.get('sequence',{}).get('length',0),entry['primaryAccession'])
        choices = sorted(choices,key=rank)[:1]
        method = 'ranked_gene_taxon_match'
    else:
        method = 'unique_gene_taxon_match'
    if len(choices) != 1:
        raise ValueError(f'Ambiguous/missing gene query: {row["gene"]}; pin one of {[e["primaryAccession"] for e in choices]}')
    accession = choices[0]['primaryAccession']
    entry = client.get('https://rest.uniprot.org/uniprotkb/' + accession + '.json')
    return entry, 'repaired_identity' if requested else method, [e['primaryAccession'] for e in exact]

def experimental_annotations(entry):
    """Preserve exact annotated claims: ECO:0000269 is experimental assertion evidence."""
    result = []
    def walk(value, path):
        if isinstance(value,dict):
            for evidence in value.get('evidences',[]):
                if evidence.get('evidenceCode') == 'ECO:0000269':
                    result.append(dict(path=path,claim=value.get('value',value.get('description','')),
                                       code=evidence['evidenceCode'],source=evidence.get('source',''),id=evidence.get('id','')))
            for key, item in value.items():
                if key != 'evidences':
                    walk(item,path+'/'+key)
        elif isinstance(value,list):
            for i,item in enumerate(value):
                walk(item,path+'/'+str(i))
    walk(entry,'')
    return result

def build(args):
    out = Path(args.out)
    out.mkdir(parents=True,exist_ok=True)
    client = UniProtClient(args.cache,args.offline)
    rules = json.loads(Path(args.rules).read_text(encoding='utf-8'))['families']
    sets, audit, errors, evidence_rows = {}, [], [], []
    datasets = [('adhesome',args.curated_seed)]
    if getattr(args,'anchor_seed',None):
        datasets.append(('anchors',args.anchor_seed))
    if getattr(args,'supplemental_seed',None):
        datasets.append(('supplemental',args.supplemental_seed))
    datasets.append(('metazoan',args.metazoan_seed))
    for label, seed in datasets:
        records, rows, seen = [], [], set()
        for row in table(seed):
            rid = row['reference_id']
            try:
                if not re.fullmatch(r'[A-Za-z0-9_.-]+',rid) or rid in seen:
                    raise ValueError('Invalid or duplicate reference ID')
                seen.add(rid)
                if row['family'] not in rules and row['family'] != 'non_adhesome':
                    raise ValueError('Family not present in classification rules')
                if rules.get(row['family'],{}).get('screening_only') and label != 'metazoan':
                    raise ValueError('Generic screening buckets cannot label forward reference families')
                if label == 'adhesome' and any(not row.get(k) for k in ('module','evidence_source','experimental_evidence')):
                    raise ValueError('Missing curated evidence fields')
                entry,method,choices = resolve(row,client,args.repair_identity,args.ambiguity_policy)
                if not valid_identity(entry,row):
                    raise ValueError('Resolved record identity mismatch')
                seq = entry.get('sequence',{}).get('value','').upper()
                if not seq or set(seq)-set('ACDEFGHIKLMNPQRSTVWYXBZJUO'):
                    raise ValueError('Invalid sequence')
                reviewed = 'Swiss-Prot' in entry.get('entryType','')
                if args.strict_reviewed_curated and label=='adhesome' and not reviewed:
                    raise ValueError('Curated reference is not reviewed')
                exp = experimental_annotations(entry)
                accession = entry['primaryAccession']
                name = entry.get('proteinDescription',{}).get('recommendedName',{}).get('fullName',{}).get('value',entry.get('uniProtkbId',''))
                meta = dict(reference_id=rid,family=row['family'],module=row.get('module',''),species=entry['organism']['scientificName'],
                    taxon_id=entry['organism']['taxonId'],gene=';'.join(gene_names(entry)),accession=accession,requested_accession=row.get('accession',''),
                    protein_name=name,length=len(seq),sequence_sha256=hashlib.sha256(seq.encode()).hexdigest(),
                    sequence_version=entry.get('entryAudit',{}).get('sequenceVersion',''),entry_version=entry.get('entryAudit',{}).get('entryVersion',''),
                    reviewed=reviewed,resolve_method=method,source_url='https://www.uniprot.org/uniprotkb/'+accession+'/entry',
                    evidence_source=row.get('evidence_source',''),experimental_evidence=row.get('experimental_evidence',''),
                    evidence_status='seed_claim_with_experimental_annotations' if exp else 'seed_claim_requires_review',
                    adhesome_experiment_verified='false',experimental_annotation_count=len(exp),
                    experimental_pubmed_ids=';'.join(sorted({e['id'] for e in exp if e['source']=='PubMed'})),role=row.get('role','adhesome_seed'),
                    reference_tier='computational_anchor' if row.get('role')=='phylogenetic_anchor' else 'curated_seed' if label=='adhesome' else 'panel_reference',
                    phylogenetic_group=row.get('phylogenetic_group',''),family_assignment_basis=row.get('family_assignment_basis','seed assertion'),
                    anchor_limitations=row.get('anchor_limitations',''))
                rows.append(meta)
                records.append((rid,seq))
                audit.append(dict(dataset=label,**meta,query_candidates=choices))
                evidence_rows.extend(dict(dataset=label,reference_id=rid,accession=accession,**item) for item in exp)
                print(f'{label}: {rid} -> {accession} ({len(seq)} aa; {method})',flush=True)
            except Exception as error:
                errors.append(dict(dataset=label,reference_id=rid,error=str(error)))
                print(f'REVIEW {label}: {rid}: {error}',flush=True)
        sets[label] = (rows,records)
    save_json(out/'resolution_audit.json',dict(records=audit,errors=errors,requests=client.requests))
    write_table(out/'resolution_errors.csv',errors,['dataset','reference_id','error'])
    if errors:
        raise ValueError(f'{len(errors)} reference identity/resolution errors; see {out}/resolution_errors.csv. No new reference set published.')
    if not any(r['family']=='non_adhesome' for r in sets['metazoan'][0]):
        raise ValueError('Reverse-search panel lacks competitors')
    # Check shared records before publishing either pair.
    forward_rows = sets['adhesome'][0] + sets.get('anchors',([],[]))[0] + sets.get('supplemental',([],[]))[0]
    if len({r['reference_id'] for r in forward_rows}) != len(forward_rows):
        raise ValueError('Duplicate curated/anchor forward reference IDs')
    curated = {r['reference_id']:r for r in forward_rows}
    broad = {r['reference_id']:r for r in sets['metazoan'][0]}
    for rid,row in curated.items():
        if rid not in broad or any(row[k]!=broad[rid][k] for k in ('family','accession','sequence_sha256')):
            raise ValueError(f'Curated reference missing/inconsistent in reciprocal panel: {rid}')
    fields = list(sets['adhesome'][0][0])
    for label,(rows,records) in sets.items():
        fasta_name = {'adhesome':'adhesome_references.faa','anchors':'phylogenetic_anchors.faa','metazoan':'metazoan_references.faa','supplemental':'supplemental_references.faa'}[label]
        csv_name = {'adhesome':'reference_catalogue.csv','anchors':'anchor_catalogue.csv','metazoan':'metazoan_catalogue.csv','supplemental':'supplemental_catalogue.csv'}[label]
        write_fasta(out/fasta_name,records)
        write_table(out/csv_name,rows,fields)
        individual = out / 'individual' / label
        individual.mkdir(parents=True,exist_ok=True)
        for rid,seq in records:
            write_fasta(individual/(rid+'.faa'),[(rid,seq)])
        assert set(sequence_dict(out/fasta_name)) == {r['reference_id'] for r in table(out/csv_name)}
    write_table(out/'experimental_annotations.csv',evidence_rows,['dataset','reference_id','accession','path','claim','code','source','id'])
    data_files=['adhesome_references.faa','reference_catalogue.csv','metazoan_references.faa','metazoan_catalogue.csv','experimental_annotations.csv']
    if 'supplemental' in sets:
        data_files += ['supplemental_references.faa','supplemental_catalogue.csv']
    if 'anchors' in sets:
        write_fasta(out/'search_references.faa',sets['adhesome'][1]+sets['anchors'][1]+sets.get('supplemental',([],[]))[1])
        write_table(out/'search_catalogue.csv',forward_rows,fields)
        data_files += ['phylogenetic_anchors.faa','anchor_catalogue.csv','search_references.faa','search_catalogue.csv']
        coverage=[]
        for family,rule in rules.items():
            if rule.get('screening_only'):
                continue
            counts=Counter(r['phylogenetic_group'] for r in sets['anchors'][0] if r['family']==family)
            coverage.append(dict(family=family,curated_seeds=sum(r['family']==family for r in sets['adhesome'][0]),
                supplemental_references=sum(r['family']==family for r in sets.get('supplemental',([],[]))[0]),
                platyhelminth_anchors=counts['Platyhelminthes'],molluscan_anchors=counts['Mollusca'],annelid_anchors=counts['Annelida'],
                reciprocal_references=sum(r['family']==family for r in sets['metazoan'][0]),
                coverage_gap='' if sum(counts.values()) else 'no_selected_spiralian_anchor'))
        write_table(out/'anchor_coverage.csv',coverage,list(coverage[0]))
        data_files.append('anchor_coverage.csv')
    save_json(out/'reference_lock.json',dict(schema_version=1,created_at=datetime.now(timezone.utc).isoformat(),
        interpretation='Targeted metazoan competitor panel; seed assertions are not automatically verified adhesome experiments.',
        ambiguity_policy=args.ambiguity_policy,repair_identity=args.repair_identity,
        builder_sha256=digest(__file__),rules_sha256=digest(args.rules),
        seeds={str(p):digest(p) for label,p in datasets},
        files={name:digest(out/name) for name in data_files},
        species_counts={label:dict(Counter(r['species'] for r in rows)) for label,(rows,records) in sets.items()},
        references=audit,requests=client.requests))
    print(f'Published paired reference sets in {out}',flush=True)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--curated-seed',default=str(DEFAULT_SEEDS/'curated_reference_seed.csv'))
    p.add_argument('--metazoan-seed',default=str(DEFAULT_SEEDS/'metazoan_reference_seed.csv'))
    p.add_argument('--anchor-seed',default=str(DEFAULT_SEEDS/'phylogenetic_anchor_seed.csv'))
    p.add_argument('--supplemental-seed',default=str(DEFAULT_SEEDS/'supplemental_reference_seed.csv'))
    p.add_argument('--rules',default=str(ROOT/'config/family_rules.json'))
    p.add_argument('--out',default=str(ROOT/'references/schistosome_v5'))
    p.add_argument('--cache',default=str(ROOT/'references/cache/uniprot'))
    p.add_argument('--offline',action='store_true')
    p.add_argument('--repair-identity',action='store_true',help='Resolve incorrect supplied accessions by unique gene/taxon match, retaining provenance')
    p.add_argument('--strict-reviewed-curated',action='store_true')
    p.add_argument('--ambiguity-policy',choices=['error','prefer-complete'],default='error',
                   help='For matching genes only: prefer reviewed/nonfragment/protein evidence/longest; record all alternatives')
    build(p.parse_args())

if __name__=='__main__':
    main()
