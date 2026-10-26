#!/usr/bin/env python3
"""Retrieve real UniProt sequences; never substitute missing sequences.
Run: python build_calibration_set.py --components components.csv --out output
Requires: requests. Internet access required.
"""
import argparse,csv,hashlib,json,pathlib,random,re,time
import requests

BASE='https://rest.uniprot.org/uniprotkb'
HARD={
  "talin": [
    "HIP1",
    "HIP1R",
    "EZR",
    "RDX",
    "MSN",
    "NF2",
    "EPB41",
    "EPB41L1",
    "EPB41L2",
    "EPB41L3",
    "EPB41L4A",
    "EPB41L4B",
    "EPB41L5",
    "FRMD3",
    "FRMD4A",
    "FRMD4B",
    "FRMD5",
    "FRMD6",
    "FRMD7",
    "FRMD8"
  ],
  "alpha-actinin": [
    "SPTA1",
    "SPTAN1",
    "SPTB",
    "SPTBN1",
    "SPTBN2",
    "SPTBN4",
    "SPTBN5",
    "DMD",
    "UTRN",
    "DRP2",
    "PLEC",
    "MACF1",
    "DST",
    "FLNA",
    "FLNB",
    "FLNC",
    "CNN1",
    "CNN2",
    "CNN3",
    "TAGLN",
    "TAGLN2",
    "TAGLN3",
    "IQGAP1",
    "IQGAP2",
    "IQGAP3",
    "VAV1",
    "VAV2",
    "VAV3",
    "SPECC1",
    "SPECC1L",
    "GAS2",
    "GAS2L1",
    "GAS2L2",
    "GAS2L3"
  ],
  "vinculin": [
    "CTNNA1",
    "CTNNA2",
    "CTNNA3",
    "CTNNB1",
    "CTNND1",
    "CTNND2",
    "JUP",
    "APC",
    "APC2",
    "ARMC1",
    "ARMC2",
    "ARMC3",
    "ARMC4",
    "ARMC5",
    "ARMC6",
    "ARMC7",
    "ARMC8",
    "ARMC9",
    "ARMC10",
    "ARMC12",
    "ARMCX1",
    "ARMCX2",
    "ARMCX3",
    "ARMCX4",
    "ARMCX5",
    "ARMCX6",
    "CTNNBL1",
    "KPNA1",
    "KPNA2",
    "KPNA3",
    "KPNA4",
    "KPNA5",
    "KPNA6",
    "KPNA7",
    "KPNB1",
    "IPO4",
    "IPO5",
    "IPO7",
    "IPO8",
    "IPO9",
    "IPO11",
    "IPO13",
    "XPO1",
    "XPO4",
    "XPO5",
    "XPO6",
    "XPO7",
    "CSE1L",
    "TNPO1",
    "TNPO2",
    "TNPO3"
  ],
  "PINCH": [
    "LHX1",
    "LHX2",
    "LHX3",
    "LHX4",
    "LHX5",
    "LHX6",
    "LHX8",
    "LHX9",
    "ISL1",
    "ISL2",
    "LMO1",
    "LMO2",
    "LMO3",
    "LMO4",
    "PDLIM1",
    "PDLIM2",
    "PDLIM3",
    "PDLIM4",
    "PDLIM5",
    "PDLIM7",
    "LIMK1",
    "LIMK2",
    "LMO7",
    "CSRP1",
    "CSRP2",
    "CSRP3",
    "CRIP1",
    "CRIP2",
    "LASP1",
    "FHL1",
    "FHL2",
    "FHL3",
    "FHL5",
    "ZYX",
    "LPP",
    "TRIP6",
    "AJUBA",
    "LIMD1",
    "WTIP"
  ]
}

def rows(path):
 return list(csv.DictReader(open(path,encoding='utf-8-sig')))
def accessions(row):
 return [a.strip() for a in row['Swiss-Prot ID'].split(',') if a.strip() and a.strip().lower()!='null']
def table(path,data,fields):
 with open(path,'w',newline='') as f:
  w=csv.DictWriter(f,fields,extrasaction='ignore');w.writeheader();w.writerows(data)
def fetch(session,url,params=None):
 for attempt in range(5):
  r=session.get(url,params=params,timeout=120)
  if r.status_code==429 or r.status_code>=500:
   time.sleep(2**attempt);continue
  r.raise_for_status();return r
 raise RuntimeError('Repeated HTTP failure: '+url)
def search(session,query):
 url=BASE+'/search';params={'query':query,'format':'json','size':500}
 while url:
  r=fetch(session,url,params);yield from r.json()['results']
  url=r.links.get('next',{}).get('url');params=None

def genes(e):
 return {v['value'] for g in e.get('genes',[]) for k in ('geneName','synonyms','orderedLocusNames','orfNames') for v in ([g[k]] if k=='geneName' and k in g else g.get(k,[]))}
def pfams(e):
 return {x['id'] for x in e.get('uniProtKBCrossReferences',[]) if x['database']=='Pfam'}
def name(e):
 p=e.get('proteinDescription',{});return p.get('recommendedName',{}).get('fullName',{}).get('value',str(p))
def seq(e):return e['sequence']['value']
def fasta(path,entries,label):
 with open(path,'w') as f:
  for e in entries:
   s=seq(e);assert s and re.fullmatch('[ACDEFGHIKLMNPQRSTVWYBXZJUO]+',s)
   g=','.join(sorted(genes(e))) or 'NA'
   f.write('>'+label+'|'+e['primaryAccession']+'|'+g+' taxid='+str(e['organism']['taxonId'])+' '+name(e)+'\n')
   f.write('\n'.join(s[i:i+60] for i in range(0,len(s),60))+'\n')

def main():
 p=argparse.ArgumentParser();p.add_argument('--components',default='components.csv');p.add_argument('--out',default='output');p.add_argument('--background-count',type=int,default=232);p.add_argument('--seed',type=int,default=20261014);a=p.parse_args()
 out=pathlib.Path(a.out);out.mkdir(parents=True,exist_ok=True)
 source=rows(a.components);s=requests.Session();s.headers['User-Agent']='components-fasta-builder/1.0'
 mapping={};issues=[]
 for i,row in enumerate(source,1):
  aa=accessions(row)
  if not aa:issues.append({'row':i,'symbol':row['Official Symbol'],'issue':'No accession; not inferred'})
  for ac in aa:mapping.setdefault(ac,[]).append((i,row))
 candidates={};raw={};audit=[]
 for ac,rr in mapping.items():
  try:
   r=fetch(s,BASE+'/'+ac+'.json');e=r.json();raw[ac]=e
   if e['organism']['taxonId']!=9606:raise ValueError('Nonhuman accession requires manual review')
   candidates[e['primaryAccession']]=e
   for i,row in rr:
    audit.append({'row':i,'source_symbol':row['Official Symbol'],'requested_accession':ac,'resolved_accession':e['primaryAccession'],'resolved_genes':';'.join(sorted(genes(e))),'symbol_match':row['Official Symbol'] in genes(e),'length':len(seq(e))})
  except Exception as ex:issues.append({'row':';'.join(str(i) for i,_ in rr),'symbol':';'.join(x['Official Symbol'] for _,x in rr),'issue':ac+': '+str(ex)})
 if not candidates:raise RuntimeError('No sequences retrieved; no FASTA files generated')
 fasta(out/'candidates.fasta',list(candidates.values()),'candidate')
 # Keep all supplied accessions, including cross-gene entries, with explicit audit.
 # These are library candidates, not experimentally confirmed positives.
 hard={};hard_audit=[];source_genes={r['Official Symbol'] for r in source}
 for target,gg in HARD.items():
  for gene in gg:
   if gene in source_genes:continue
   found=list(search(s,'(gene_exact:'+gene+') AND (organism_id:9606) AND (reviewed:true)'))
   found=[e for e in found if gene in genes(e) and e['primaryAccession'] not in candidates]
   if len(found)!=1:
    issues.append({'row':'','symbol':gene,'issue':'Hard-negative query returned '+str(len(found))+' exact matches; skipped'});continue
   e=found[0];hard[e['primaryAccession']]=e
   hard_audit.append({'gene':gene,'accession':e['primaryAccession'],'negative_for':target,'scope':'family-specific challenge, not global non-adhesome claim'})
 fasta(out/'hard_negatives.fasta',list(hard.values()),'hard_negative')
 # Complete human reference-proteome scope; reviewed subset provides independent annotations.
 pool=list(search(s,'(proteome:UP000005640) AND (reviewed:true)'))
 target_pfams=set().union(*(pfams(e) for e in candidates.values()))
 target_gene_names=source_genes|set().union(*(genes(e) for e in candidates.values()))
 excluded=set(candidates)|set(hard);excluded_sequences={seq(e) for e in list(candidates.values())+list(hard.values())}
 eligible=[]
 for e in pool:
  if e['primaryAccession'] in excluded or seq(e) in excluded_sequences:continue
  if genes(e)&target_gene_names or pfams(e)&target_pfams:continue
  # Require independent curated enzyme annotation. This conservative background is
  # enzyme-enriched, not a representative random sample of the entire proteome.
  pd=e.get('proteinDescription',{}).get('recommendedName',{})
  if not pd.get('ecNumbers'):continue
  if e.get('proteinDescription',{}).get('flag')=='Fragment':continue
  eligible.append(e)
 eligible.sort(key=lambda e:e['primaryAccession']);rng=random.Random(a.seed);rng.shuffle(eligible)
 background=[];seen=set(excluded_sequences)
 for e in eligible:
  if seq(e) in seen:continue
  background.append(e);seen.add(seq(e))
  if len(background)>=a.background_count:break
 if len(background)<a.background_count:issues.append({'row':'','symbol':'background','issue':'Only '+str(len(background))+' eligible unique sequences'})
 fasta(out/'background_negatives.fasta',background,'background_negative_provisional')
 with open(out/'components_with_negatives.fasta','w') as f:
  for n in ['candidates.fasta','hard_negatives.fasta','background_negatives.fasta']:f.write((out/n).read_text())
 table(out/'candidate_mapping_audit.csv',audit,['row','source_symbol','requested_accession','resolved_accession','resolved_genes','symbol_match','length'])
 table(out/'hard_negative_manifest.csv',hard_audit,['gene','accession','negative_for','scope'])
 table(out/'issues.csv',issues,['row','symbol','issue'])
 manifest=[]
 for label,entries in [('candidate',list(candidates.values())),('hard_negative',list(hard.values())),('background_negative_provisional',background)]:
  for e in entries:manifest.append({'class':label,'accession':e['primaryAccession'],'genes':';'.join(sorted(genes(e))),'protein_name':name(e),'length':len(seq(e)),'sha256':hashlib.sha256(seq(e).encode()).hexdigest(),'source':'https://www.uniprot.org/uniprotkb/'+e['primaryAccession'],'pfam':';'.join(sorted(pfams(e)))})
 table(out/'sequence_manifest.csv',manifest,['class','accession','genes','protein_name','length','sha256','source','pfam'])
 (out/'retrieved_records.json').write_text(json.dumps({'candidates':raw,'hard_negatives':list(hard.values()),'background_negatives':background}))
 summary={'candidate_sequences':len(candidates),'hard_negative_sequences':len(hard),'background_sequences':len(background),'seed':a.seed,'background_pool':'UP000005640 reviewed enzyme subset; no candidate Pfam overlap','warning':'Background labels provisional; full family catalogue not supplied. Hard negatives are target-family-specific. Do not treat all as globally non-adhesome.'}
 (out/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
