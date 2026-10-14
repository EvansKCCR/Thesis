#!/usr/bin/env python3
"""Integrated conservative adhesome audit and duplicate-collapse pipeline.

Usage:
  python integrated_adhesome_pipeline.py INPUT.xlsx OUTPUT.xlsx [--sheet candidates] [--keep-all-audited]

Scoring: architecture 0-4, orthology 0-3, motif 0-2, localisation/topology 0-1.
Architecture is a mandatory family gate. A high total cannot rescue a failed gate.
"""
from __future__ import annotations
import argparse, hashlib, json, re
from collections import Counter, defaultdict
from pathlib import Path
from openpyxl import load_workbook
from openpyxl.styles import Font, PatternFill, Alignment
from openpyxl.utils import get_column_letter

VERSION='1.0'
TARGET={'PINCH_like','parvin_like','Csk','Syk','Nck','Grb2','Crk','Shc'}
GRADE_RANK={'A':4,'B':3,'C':2,'D':1}
DECISION_RANK={'Retain family assignment':4,'Retain provisionally':3,'Family-level provisional only':2,'Ambiguous; manual review':1,'Do not assign at family level':0}
STATUS_RANK={'supported_candidate':3,'supported candidate':3,'provisional':2,'unresolved':1}
REF=('H. sapiens','C. elegans','D. melanogaster','M. musculus','X. laevis')


def sha256(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1048576),b''): h.update(b)
 return h.hexdigest()
def t(r,k): return str(r.get(k) or '')
def pos(a,x): return [(int(i),int(j)) for i,j in re.findall(rf'{re.escape(x)}:(\d+)-(\d+)',str(a or ''))]
def cnt(a,x): return len(pos(a,x))
def ordered(a, classes):
 last=-1
 for cl in classes:
  opts=cl if isinstance(cl,tuple) else (cl,)
  hits=sorted((s,e) for x in opts for s,e in pos(a,x))
  hit=next((z for z in hits if z[0]>last),None)
  if not hit:return False
  last=hit[0]
 return True
def ortho(r):
 s=t(r,'Orthology_reference_support'); n=int(r['Orthology_direct_reference_count']) if r.get('Orthology_direct_reference_count') not in (None,'') else sum(x in s for x in REF); shared=t(r,'Orthology_HOG_status')=='shared HOG'
 return 3 if shared and n>=4 else 2 if shared and n else 1 if shared else 0
def motif(r):
 tier=t(r,'MotifPrediction_candidate_tiers')
 if 'context_supported_candidate' in tier or float(r.get('MotifPrediction_context_supported_count') or 0)>0:return 2
 if 'sequence_only_candidate' in tier or float(r.get('MotifPrediction_priority_site_count') or 0)>0:return 1
 return 0
def loc(r,cytosolic=False):
 dl,tm,sp=t(r,'DeepLoc_2.1'),t(r,'DeepTMHMM'),t(r,'SignalP.v6.0')
 intra=any(x in dl for x in ('Cytoplasm','Nucleus','Cell membrane')); secret='Extracellular' in dl or sp.startswith('SP') or bool(re.search(r'(^|;)\s*SP\b',tm)); trans=bool(re.search(r'(^|;)\s*TM\b',tm)) and 'no TM' not in tm
 return int(intra and not secret and (not cytosolic or not trans))
def specific(r,f):
 s=';'.join(t(r,x) for x in ('NCBI_CDD_architecture','Orthology_reference_support','Orthology_interpretation'))
 p={'PINCH_like':(r'LIM[1-5]_PINCH',r'\bPINCH\b',r'\bLIMS[12]?\b'),'parvin_like':(r'PARVIN',r'PARVA',r'PARVB',r'PARVG',r'ACTOPAXIN',r'AFFIXIN'),'Csk':(r'PTKc_Csk',r'C-terminal Src kinase',r'(?:^|[_ ;])CSK(?:[_ ;]|$)',r'(?:^|[_ ;])CHK(?:[_ ;]|$)'),'Syk':(r'PTKc_Syk',r'(?:^|[_ ;])SYK(?:[_ ;]|$)',r'ZAP[-_ ]?70'),'Nck':(r'NCK[12]?(?:_|\b)',),'Grb2':(r'GRB2(?:_|\b)',r'SEM[-_ ]?5'),'Crk':(r'CRK(?:L)?(?:_|\b)',),'Shc':(r'PTB_Shc',r'SH2_SHC',r'SHC[1-4]?(?:_|\b)')}
 return any(re.search(x,s,re.I) for x in p.get(f,()))

def family_domain(r,f):
 if f=='fibronectin_III_ligand': return 0,'generic FNIII screening bucket; not a fibronectin-family assignment'
 a=t(r,'architecture'); c=t(r,'NCBI_CDD_architecture'); L=float(r.get('length') or 0); kin=('PF00069','PF07714')
 if f=='PINCH_like':
  ok=cnt(a,'PF00412')==5 and 250<=L<=500 and len(set(re.findall(r'LIM[1-5]_PINCH',c,re.I)))>=2
  return (4,'five LIM domains plus PINCH-specific CDD') if ok else (0,'PINCH gate failed')
 if f=='parvin_like':
  ok=cnt(a,'PF00307')>=2 and 250<=L<=700 and specific(r,f)
  return (4,'tandem CH domains plus parvin-specific evidence') if ok else (0,'parvin gate failed')
 if f=='Syk':
  comp=cnt(a,'PF00017')>=2 and max(cnt(a,x) for x in kin)>=1
  if comp and ordered(a,['PF00017','PF00017',kin]) and specific(r,f): return 4,'SH2-SH2-kinase plus Syk-specific evidence'
  if comp and specific(r,f): return 3,'Syk-compatible composition; order unresolved'
  return 0,'Syk gate failed'
 if f=='Csk':
  arch=cnt(a,'PF00018')>=1 and cnt(a,'PF00017')>=1 and max(cnt(a,x) for x in kin)>=1 and ordered(a,['PF00018','PF00017',kin])
  if arch and specific(r,f):return 4,'SH3-SH2-kinase plus Csk/Chk-specific evidence'
  if arch:return 2,'SH3-SH2-kinase is compatible but shared with Src'
  return 0,'Csk gate failed'
 if f=='Nck':
  ok=cnt(a,'PF00018')==3 and cnt(a,'PF00017')==1 and ordered(a,['PF00018','PF00018','PF00018','PF00017']) and specific(r,f)
  return (4,'SH3-SH3-SH3-SH2 plus Nck-specific evidence') if ok else (0,'Nck gate failed')
 if f=='Grb2':
  ok=cnt(a,'PF00018')==2 and cnt(a,'PF00017')==1 and ordered(a,['PF00018','PF00017','PF00018']) and specific(r,f)
  return (4,'SH3-SH2-SH3 plus Grb2-specific evidence') if ok else (0,'Grb2 gate failed')
 if f=='Crk':
  ok=cnt(a,'PF00017')==1 and 1<=cnt(a,'PF00018')<=2 and ordered(a,['PF00017','PF00018']) and specific(r,f)
  return (4,'SH2 followed by one/two SH3 domains plus Crk-specific evidence') if ok else (0,'Crk gate failed')
 if f=='Shc':
  ok=cnt(a,'PF00640')>=1 and cnt(a,'PF00017')>=1 and ordered(a,['PF00640','PF00017']) and specific(r,f)
  return (4,'PTB-SH2 plus Shc-specific evidence') if ok else (0,'Shc gate failed')
 # Existing-family generic architecture rule.
 dm=r.get('domain_match'); tf=float(r.get('domain_type_fraction') or 0); cf=float(r.get('domain_copy_fraction') or 0)
 if dm in (True,1,'1'): return 4,'configured family architecture matched'
 frac=min(tf,cf)
 return (2,'partial configured architecture') if frac>=.5 else (1,'weak partial architecture') if frac>0 else (0,'no configured architecture support')

def audit(r,f=None,discovery='original_hypothesis'):
 f=f or t(r,'family'); d,ev=family_domain(r,f); o=ortho(r); m=motif(r); l=loc(r,f in TARGET); total=d+o+m+l
 gate=d==4
 if gate and o>=2:grade,dec='A','Retain family assignment'
 elif gate:grade,dec='B','Retain provisionally'
 elif d>=2 and o>=2:grade,dec='C','Family-level provisional only'
 else:grade,dec='D','Do not assign at family level'
 assigned=f if dec!='Do not assign at family level' else 'unassigned'
 z=dict(r);z.update({'preaudit_family':t(r,'family'),'assigned_family':assigned,'audit_family':assigned,'architecture_evidence':ev,'domain_score_0_4':d,'orthology_score_0_3':o,'motif_score_0_2':m,'localization_score_0_1':l,'total_score_0_10':total,'grade':grade,'family_decision':dec,'rationale':ev if gate else ev+'; supporting evidence cannot rescue the failed family gate','discovery_rule':discovery,'pipeline_version':VERSION})
 return z

def global_recovery(base):
 out=[]
 for f in ('PINCH_like','parvin_like'):
  d,_=family_domain(base,f)
  if d==4:out.append(audit(base,f,'global_architecture_recovery'))
 return out

def rank(r):
 assigned=int(t(r,'assigned_family') not in ('','unassigned','unresolved'))
 decision=DECISION_RANK.get(t(r,'family_decision'),0); grade=GRADE_RANK.get(t(r,'grade'),0); d=float(r.get('domain_score_0_4') or 0); o=float(r.get('orthology_score_0_3') or 0);m=float(r.get('motif_score_0_2') or 0);l=float(r.get('localization_score_0_1') or 0);tot=float(r.get('total_score_0_10') or 0);specificity=int(t(r,'discovery_rule')=='global_architecture_recovery')
 return (assigned,decision,grade,d,specificity,o,m,l,tot,STATUS_RANK.get(t(r,'status'),0))

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--version',action='version',version='adhesome-audit '+VERSION);ap.add_argument('input',type=Path);ap.add_argument('output',type=Path);ap.add_argument('--sheet',default='candidates');ap.add_argument('--keep-all-audited',action='store_true');a=ap.parse_args()
 wb=load_workbook(a.input);src=wb[a.sheet];headers=[c.value for c in src[1] if c.value not in (None,'')]
 raw=[]
 for vals in src.iter_rows(min_row=2,max_col=len(headers),values_only=True):
  if any(v not in (None,'') for v in vals):raw.append(dict(zip(headers,vals)))
 audited=[audit(r) for r in raw];groups=defaultdict(list)
 for r in audited:groups[(r.get('Specie'),r.get('protein_id'))].append(r)
 recovered=[]
 for key,rs in list(groups.items()):
  new=global_recovery(max(rs,key=rank));groups[key].extend(new);recovered.extend(new)
 collapsed=[]
 for key,rs in groups.items():
  win=max(rs,key=rank);alts=sorted({t(x,'preaudit_family') for x in rs if t(x,'preaudit_family') and x is not win});win=dict(win);win['alternative_family_hypotheses']='; '.join(alts);win['collapsed_row_count']=sum(t(x,'discovery_rule')=='original_hypothesis' for x in rs);collapsed.append(win)
 extras=['preaudit_family','audit_family','architecture_evidence','domain_score_0_4','orthology_score_0_3','motif_score_0_2','localization_score_0_1','total_score_0_10','grade','family_decision','rationale','discovery_rule','pipeline_version','alternative_family_hypotheses','collapsed_row_count']
 outheaders=headers+[x for x in ['assigned_family']+extras if x not in headers]
 def write_sheet(name,rows):
  if name in wb.sheetnames:del wb[name]
  sh=wb.create_sheet(name,0);sh.append(outheaders)
  for r in rows:sh.append([r.get(h) for h in outheaders])
  fill=PatternFill('solid',fgColor='1F4E78')
  for c in sh[1]:c.fill=fill;c.font=Font(color='FFFFFF',bold=True);c.alignment=Alignment(wrap_text=True)
  sh.freeze_panes='A2';sh.auto_filter.ref=sh.dimensions;sh.sheet_view.showGridLines=False
  for i,h in enumerate(outheaders,1):sh.column_dimensions[get_column_letter(i)].width=34 if any(x in h.lower() for x in ('rationale','architecture','alternative','support')) else 16
 if a.keep_all_audited:write_sheet('All_hypotheses_audited',audited)
 write_sheet('Collapsed_candidates',sorted(collapsed,key=lambda r:(t(r,'Specie'),t(r,'protein_id'))))
 if 'Pipeline_summary' in wb.sheetnames:del wb['Pipeline_summary']
 sm=wb.create_sheet('Pipeline_summary',1);fam=Counter(t(r,'assigned_family') for r in collapsed);grades=Counter(t(r,'grade') for r in collapsed);pinch=[t(r,'protein_id') for r in collapsed if t(r,'assigned_family')=='PINCH_like'];parvin=[t(r,'protein_id') for r in collapsed if t(r,'assigned_family')=='parvin_like']
 summary=[['Metric','Value'],['Pipeline version',VERSION],['Family-hypothesis rows audited',len(raw)],['Unique species-protein pairs retained',len(collapsed)],['Rows collapsed',len(raw)-len(collapsed)],['PINCH-like assignments',len(pinch)],['Parvin-like assignments',len(parvin)],['PINCH-like protein IDs',', '.join(pinch) or 'None'],['Parvin-like protein IDs',', '.join(parvin) or 'None'],['Input SHA-256',sha256(a.input)]]+[[f'Grade {g}',grades[g]] for g in 'ABCD']+[[f'Family: {k}',v] for k,v in sorted(fam.items())]
 for r in summary:sm.append(r)
 for c in sm[1]:c.fill=PatternFill('solid',fgColor='1F4E78');c.font=Font(color='FFFFFF',bold=True)
 sm.column_dimensions['A'].width=38;sm.column_dimensions['B'].width=90;sm.freeze_panes='A2';sm.sheet_view.showGridLines=False
 wb.save(a.output);manifest=a.output.with_suffix('.manifest.json');manifest.write_text(json.dumps({'version':VERSION,'input':a.input.name,'input_sha256':sha256(a.input),'output':a.output.name,'output_sha256':sha256(a.output),'rows_audited':len(raw),'unique_pairs':len(collapsed),'rows_collapsed':len(raw)-len(collapsed),'global_recovery_rows':len(recovered),'pinch_like':pinch,'parvin_like':parvin,'keep_all_audited':a.keep_all_audited},indent=2)+'\n')
 print(json.dumps({'output':str(a.output),'manifest':str(manifest),'rows_audited':len(raw),'unique_pairs':len(collapsed),'rows_collapsed':len(raw)-len(collapsed),'pinch_like':pinch,'parvin_like':parvin},indent=2))

if __name__=='__main__':main()
