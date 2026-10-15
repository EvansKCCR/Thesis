#!/usr/bin/env python3
"""Reference-gated protein-family adjudication for comparative adhesome candidates.

USAGE
 python integrated_adhesome_pipeline.py INPUT.xlsx OUTPUT.xlsx [--sheet candidates] [--rules ../01_detection/config/family_rules.json]

Expected additional columns for protein/subfamily calls (optional; without them output stays broad):
 reference_subfamily, subfamily_support_method, reciprocal_reference_subfamily,
 phylogeny_concordant, subfamily_profile_calibrated, top_profile_subfamily,
 profile_margin_pass, reference_accession, reference_dataset_version.

Accepted methods: curated_reference, phylogeny, reciprocal_best_hit, calibrated_hmm.
WARNING: Reference-supported is not experimentally demonstrated adhesome localization.
The gene symbol is not used as a classifier unless separately supplied as curated reference evidence.
"""
from __future__ import annotations
import argparse, hashlib, json, re
from collections import Counter, defaultdict
from pathlib import Path

VERSION='1.0'
DOM_RE=re.compile(r'(?<!\w)(PF\d{5})(?:\s*:\s*(\d+)\s*-\s*(\d+))?',re.I)
TRUTH={'1','TRUE','YES','Y','PASS','SUPPORTED'}
BAD={'CONFLICT','CONTRADICT','UNSUPPORTED','FAILED','UNVERIFIED'}
ORDER_KIN={'PF00069','PF07714'}


def value(r,*names):
    for name in names:
        if name in r and r[name] is not None and str(r[name]).strip():return str(r[name]).strip()
    return ''

def num(v,default=0.0):
    try:return float(str(v).replace(',',''))
    except (ValueError,TypeError):return default

def boolean(v):return str(v or '').strip().upper() in TRUTH

def norm(v):return re.sub(r'[^A-Z0-9]','',str(v or '').upper())

def checksum(path):
    h=hashlib.sha256()
    with open(path,'rb') as handle:
        for part in iter(lambda:handle.read(1<<20),b''):h.update(part)
    return h.hexdigest()


def domains(row):
    """Parse semicolon-delimited Pfam calls, de-duplicating same-model overlapping hits.

    Unknown positions are retained as separate unpositioned calls, and order is unknown.
    """
    raw=value(row,'architecture','Pfam_architecture','pfam_architecture')
    hits=[]
    for match in DOM_RE.finditer(raw):
        pf=match.group(1).upper()
        if match.group(2):
            st,en=int(match.group(2)),int(match.group(3))
            if en<st:continue
            overlap=any(pf==old[0] and old[1] is not None and min(en,old[2])-max(st,old[1])+1>=0.5*min(en-st+1,old[2]-old[1]+1) for old in hits)
            if not overlap:hits.append((pf,st,en))
        else:hits.append((pf,None,None))
    return hits


def counts(hits):return Counter(x[0] for x in hits)


def position_order(hits, ordered_groups):
    """Each requirement group must occupy a different strictly succeeding interval.
    Distinct domain positions are mandatory; absence of positions is unknown, not false.
    """
    if not ordered_groups:return True
    last=-1;used=set()
    for group in ordered_groups:
        opts=[group] if isinstance(group,str) else group
        choices=sorted((start,end,i) for i,(pf,start,end) in enumerate(hits) if pf in opts and start is not None and i not in used and start>last)
        if not choices:return False
        st,en,i=choices[0];last=st;used.add(i)
    return True


def compatible_order(hits,rule):
    if not rule.get('support_requires_order') and not rule.get('ordered_domains'):return True
    choices=rule.get('support_order_alternatives') or [rule.get('ordered_domain_groups') or rule.get('ordered_domains',[])]
    return any(position_order(hits,groups) for groups in choices if groups)


def architecture(row,rule):
    hits=domains(row);available=counts(hits)
    alternatives=rule.get('domain_alternatives',[])
    if not alternatives:return 0,'no architecture configured',0.0
    ratios=[]
    for required in alternatives:
        need=sum(int(n) for n in required.values())
        found=sum(min(int(n),available.get(pf,0)) for pf,n in required.items())
        ratios.append(found/max(need,1))
    frac=max(ratios,default=0.0)
    if frac==0:return 0,'no required Pfam anchors',frac
    full=frac>=1-1e-9
    ordered_ok=compatible_order(hits,rule) if full else False
    min_len=num(rule.get('min_length'));length=num(value(row,'length','Length'))
    length_partial=bool(min_len and length and length<min_len)
    if full and ordered_ok and not length_partial:score=4
    elif full:score=3
    elif frac>=0.5:score=2
    else:score=1
    notes=[f'domain-match fraction={frac:.2f}']
    if full and not ordered_ok:notes.append('order absent/contradicted')
    if length_partial:notes.append(f'below reference-like length ({int(min_len)} aa); may be fragment')
    if not any(h[1] is not None for h in hits):notes.append('positionless architecture; order cannot be verified')
    return score,'; '.join(notes),frac


def family_labels(rule):
    return set(norm(x) for x in rule.get('reference_labels',[]) if x) | {norm(rule.get('label',''))}


def identity(row,rule):
    """Explicit, non-circular species-aware evidence only; HOG membership never qualifies.

    Row family/gene-symbol/NCBI CDD descriptions are not used to establish paralog identity.
    """
    target=norm(value(row,'reference_subfamily','Reference_subfamily','orthology_resolved_family'))
    labels=family_labels(rule)
    if not target or target not in labels:return 0,'no matched independently supported reference label'
    method=norm(value(row,'subfamily_support_method','reference_evidence_method'))
    if method=='CURATEDREFERENCE':
        if value(row,'reference_accession','curated_reference_accession'):
            return 3,'curated reference accession + explicit subfamily identity'
        return 1,'curated_reference method missing reference_accession'
    if method=='PHYLOGENY':
        return (3,'explicit concordant phylogenetic clade') if boolean(value(row,'phylogeny_concordant')) else (1,'phylogeny label without concordance')
    if method=='RECIPROCALBESTHIT':
        reciprocal=norm(value(row,'reciprocal_reference_subfamily'))
        return (3,'reciprocal reference identity concordant') if reciprocal==target else (1,'reciprocal hit not independently concordant')
    if method=='CALIBRATEDHMM':
        top=norm(value(row,'top_profile_subfamily'))
        calibrated=boolean(value(row,'subfamily_profile_calibrated'))
        margin=boolean(value(row,'profile_margin_pass'))
        return (3,'benchmark-calibrated competing-profile winner') if calibrated and margin and top==target else (1,'profile lacking benchmark calibration/margin')
    return 1,'reference_subfamily string alone is not orthology proof'


def context(row,rule):
    """Topology has no power to discriminate paralogs; annotate contradictions only."""
    tm=value(row,'DeepTMHMM').upper();sp=value(row,'SignalP.v6.0').upper()
    dl=value(row,'DeepLoc_2.1').upper()
    sp_pred=bool(re.search(r'\bSP\b',sp)) and not 'NO SP' in sp or bool(re.search(r'\bSP\b',tm)) and not 'NO SP' in tm
    trans=bool(re.search(r'\bTM\s+\d+[-–]\d+',tm)) or bool(re.search(r'^\s*TM\b',tm)) or 'SP+TM' in tm
    if 'NO TM' in tm:trans=False
    group=rule.get('module','')
    if group in ('actin_linking','focal_adhesion_scaffolds','membrane_proximal_adaptors','downstream_signalling') and (sp_pred or trans):
        return -1,'secretory/transmembrane prediction inconsistent with typical cytosolic candidate (review exceptions)'
    if rule.get('requires_topology') and not any([tm,sp,dl]):return 0,'topology missing; not demonstrated'
    if group=='integrin_receptors':
        if trans and sp_pred:return 1,'SP and TM compatible with type-I integrin'
        if tm or sp:return -1,'integrin SP/TM architecture not confirmed'
    if group=='extracellular_matrix' and rule.get('requires_topology'):
        if sp_pred:return 1,'signal peptide supports secretory compartment'
        if tm or sp:return -1,'expected extracellular routing not supported (manual review)'
    return 0,'no discriminatory topology evidence'


def orthology_score(row,identity_strength):
    hog=value(row,'Orthology_HOG_status','HOG_status')
    shared='SHARED HOG' in hog.upper()
    # Count orthology links as context, not subfamily identity.
    linked=sum(bool(v) for k,v in row.items() if k.endswith('_orthologues'))
    linked=max(linked,int(num(value(row,'Orthology_direct_reference_count'))))
    if identity_strength>=3:return 3
    if shared or linked:return 1
    return 0


def benchmark_weight(row):
    """A measured reference-cohort F1 is a context weight, never identity proof."""
    if value(row,'Orthology_benchmark_status')!='validated':
        return 1.0,'unavailable_neutral_weight'
    import math
    try:
        weight=float(value(row,'Orthology_benchmark_weight'))
    except (TypeError,ValueError):
        return 1.0,'invalid_weight_ignored'
    if not math.isfinite(weight) or not 0<=weight<=1:
        return 1.0,'invalid_weight_ignored'
    return weight,'validated_reference_cohort_weight'


def motif_score(row):
    s=value(row,'motif_candidate_tiers','MotifPrediction_candidate_tiers','motif_context_summary').lower()
    if any(w in s for w in ('context_conflict','context_conflicted','contradict','unverified')):return 0
    if 'context_supported' in s or num(value(row,'MotifPrediction_context_supported_count'))>0:return 2
    if 'sequence_only' in s or num(value(row,'MotifPrediction_priority_site_count'))>0:return 1
    return 0


def ancestors(name,families):
    seen=set()
    while name in families and name not in seen:
        seen.add(name);name=families[name].get('parent_family')
        if name:yield name
        else:break


def row_species(row):
    explicit=value(row,'Specie','species','Species','taxon','taxon_id','Taxon')
    if explicit:return explicit
    pid=value(row,'protein_id','candidate_id','Protein_ID')
    if '__' in pid:return pid.split('__',1)[0]
    return 'species_unknown'


def assess(row,families):
    result=[]
    for name,rule in families.items():
        d,architecture_note,frac=architecture(row,rule)
        ident,identity_note=identity(row,rule)
        if d==0 and ident==3 and rule.get('allow_profile_rescue') and boolean(value(row,'profile_domain_gate_verified')):
            d=2
            architecture_note='explicit independently validated profile/fold rescue; Pfam anchor undetected; manual review'
        if d==0:continue
        loc,top_note=context(row,rule)
        screening=bool(rule.get('screening_only') or rule.get('classification_level')=='screen')
        identity_needed=bool(rule.get('requires_identity_evidence'))
        valid_arch=d==4 or (d>=2 and ident==3)
        eligible=not screening and valid_arch and (ident==3 or not identity_needed) and loc>=0
        entry={'name':name,'family_group':rule.get('family_group',name),'classification_level':rule.get('classification_level','family'),
               'score_architecture':d,'architecture_fraction':round(frac,3),'architecture_note':architecture_note,'score_identity':ident,'identity_note':identity_note,
               'topology_score':loc,'topology_note':top_note,'eligible':eligible,'screen':screening,
               'partial':d<4,'parent':rule.get('parent_family') or ''}
        result.append(entry)
    return result


def audit_record(row,rules):
    """Audit independent of XLSX I/O; deterministic and unit-testable."""
    fam=rules['families']
    evaluations=assess(row,fam)
    viable=[v for v in evaluations if v['eligible']]
    leaf=[x for x in viable if not any(x['name'] in set(ancestors(y['name'],fam)) for y in viable if y['name']!=x['name'])]
    leaf.sort(key=lambda z:(-z['score_identity'],-z['score_architecture'],z['name']))
    state='unresolved';assigned='unassigned';broad='unassigned';label='';note='';winner=None
    if len(leaf)==1:
        winner=leaf[0]
        assigned=winner['name'];broad=winner['family_group']
        label=assigned if winner['classification_level']=='subfamily' else ''
        state='reference_supported' if winner['score_architecture']==4 else 'provisional_fragment'
        note=f'{winner["architecture_note"]}; {winner["identity_note"]}'
    elif len(leaf)>1:
        state='ambiguous';note='Competing eligible non-ancestor hypotheses: '+', '.join(z['name'] for z in leaf)
    else:
        strong_arch=[x for x in evaluations if x['score_architecture']>=3]
        scored_groups=defaultdict(int)
        for item in strong_arch:scored_groups[item['family_group']]=max(scored_groups[item['family_group']],item['score_architecture'])
        # In no-label cases report broad fold only, not a fabricated protein identity.
        supported=[x for x in scored_groups if x in {'actin_superfamily','ADF_H_superfamily','LIM_scaffolds','nonreceptor_SH3_SH2_kinases'}]
        if len(supported)==1:
            broad=supported[0];state='broad_only';note=f'Broad {broad} conserved domain architecture; no validated resolving evidence'
        elif len(supported)>1:
            state='ambiguous';note='Multiple unrelated broad architecture hypotheses: '+', '.join(sorted(supported))
        else:
            note='No unique family hypothesis passes reference/topology/architecture gates'
    # Flag contradictory strong evidence if a matched reference belongs to a different family group.
    strong_other=[e for e in evaluations if e['score_identity']==3 and e['score_architecture']>=2 and (not winner or e['family_group']!=winner['family_group'])]
    if winner and strong_other:
        state='ambiguous';assigned='unassigned';label='';note='Conflicting independent reference identities across superfamilies'
    ao=winner['score_architecture'] if winner else max((x['score_architecture'] for x in evaluations),default=0)
    oi=winner['score_identity'] if winner else 0
    raw_ortho=orthology_score(row,oi)
    weight,weight_status=benchmark_weight(row)
    # Curated identity, RBH identity, profiles and candidate-specific phylogeny
    # are not calibrated by a benchmark of orthogroup partitions.
    independent_identity=oi>=3
    ortho=raw_ortho if independent_identity else round(raw_ortho*weight,6)
    mo=motif_score(row)
    top=winner['topology_score'] if winner else 0
    total=ao+ortho+mo+max(top,0)
    if state=='reference_supported':grade='A' if winner and winner['classification_level']=='subfamily' else 'B'
    elif state=='provisional_fragment':grade='C'
    elif state=='broad_only':grade='C'
    else:grade='D'
    rank= {'reference_supported':3,'provisional_fragment':2,'broad_only':1,'ambiguous':0,'unresolved':0}[state]
    # Original family assertion never treated as independent identity evidence.
    out=dict(row)
    out.update({'preaudit_family':value(row,'family'),'assigned_family':assigned,'broad_homology_group':broad,
      'assigned_subfamily':label,'classification_state':state,'grade':grade,'family_decision':state.replace('_',' '),
      'architecture_evidence':note,'domain_score_0_4':ao,'orthology_score_0_3':ortho,'motif_score_0_2':mo,
      'orthology_raw_score_0_3':raw_ortho,'orthology_benchmark_weight_used':1.0 if independent_identity else weight,
      'orthology_benchmark_weight_status':'independent_identity_preserved' if independent_identity else weight_status,
      'localization_score_0_1':max(top,0),'total_score_0_10':total,'protein_identity_score_0_3':oi,
      'adhesome_membership':'not_inferred','confidence_statement':'Heuristic evidence tier; not a calibrated probability',
      'conflict_hypotheses':'; '.join(e['name'] for e in leaf[1:]) if len(leaf)>1 else '',
      'candidate_hypotheses':'; '.join(f'{e["name"]}:{e["score_architecture"]}/{e["score_identity"]}' for e in sorted(evaluations,key=lambda a:(-a['score_architecture'],-a['score_identity'],a['name']))[:20]),
      'family_rule_version':rules.get('schema_version','unknown'),'pipeline_version':VERSION,
      'resolution_rank':rank,'rationale':note})
    return out


def unique_key(row,n):
    pid=value(row,'protein_id','candidate_id','Protein_ID')
    if row_species(row)=='species_unknown':return ('unknown_taxon_do_not_collapse',str(n))
    return (row_species(row),pid) if pid else ('missing_protein_id',str(n))


def collapse_rows(rows):
    grouped=defaultdict(list)
    for i,r in enumerate(rows):grouped[unique_key(r,i)].append(r)
    result=[]
    order={'reference_supported':4,'provisional_fragment':3,'broad_only':2,'ambiguous':1,'unresolved':0}
    for key,group in grouped.items():
        confident={x['assigned_family'] for x in group if x['classification_state']=='reference_supported'}
        broad_confident={x['broad_homology_group'] for x in group if x['classification_state']=='reference_supported'}
        group.sort(key=lambda x:(-order.get(x['classification_state'],0),-x['protein_identity_score_0_3'],-x['domain_score_0_4'],x['assigned_family']))
        out=dict(group[0])
        if len(confident)>1 or len(broad_confident)>1:
            out['classification_state']='ambiguous';out['assigned_family']='unassigned';out['assigned_subfamily']='';out['grade']='D'
            out['family_decision']='ambiguous';out['resolution_rank']=0
            out['rationale']='Conflicting row-level verified labels; manual review of isoforms/hypotheses required'
            out['conflict_hypotheses']='; '.join(sorted(confident))
        out['original_row_count']=len(group)
        out['alternative_family_hypotheses']='; '.join(sorted(set(value(g,'preaudit_family') for g in group)-{out['assigned_family']}))
        result.append(out)
    return sorted(result,key=lambda r:(row_species(r),value(r,'protein_id','candidate_id')))


def run_xlsx(source:Path,dest:Path,sheet_name:str,rules):
    """Keeps Excel workbook structure; import lazily to let tests run without spreadsheet dependency."""
    from openpyxl import load_workbook
    from openpyxl.styles import Font, PatternFill, Alignment
    from openpyxl.utils import get_column_letter
    if source.resolve()==dest.resolve():raise ValueError('Input and output must be different files')
    wb=load_workbook(source)
    names={s.lower():s for s in wb.sheetnames}
    key=names.get(sheet_name.lower())
    if key is None:raise ValueError(f'Sheet {sheet_name!r} not found; available: {wb.sheetnames}')
    ws=wb[key]
    headings=[str(c.value).strip() if c.value is not None else '' for c in ws[1]]
    if not headings or any(not x for x in headings[:2]):raise ValueError('Empty or invalid header row')
    if len(set(headings))!=len(headings):raise ValueError('Duplicate input column names')
    raw=[dict(zip(headings,v)) for v in ws.iter_rows(min_row=2,max_col=len(headings),values_only=True) if any(x is not None and str(x).strip() for x in v)]
    audited=[audit_record(row,rules) for row in raw]
    resolved=collapse_rows(audited)
    extra=['preaudit_family','assigned_family','broad_homology_group','assigned_subfamily','classification_state','grade','family_decision','architecture_evidence','protein_identity_score_0_3','domain_score_0_4','orthology_score_0_3','motif_score_0_2','localization_score_0_1','total_score_0_10','adhesome_membership','confidence_statement','candidate_hypotheses','conflict_hypotheses','rationale','original_row_count','alternative_family_hypotheses','family_rule_version','pipeline_version']
    fields=headings+[h for h in extra if h not in headings]
    def add(name,data):
        if name in wb.sheetnames:del wb[name]
        sh=wb.create_sheet(name,0);sh.append(fields)
        for item in data:sh.append([item.get(h) for h in fields])
        for cell in sh[1]:
            cell.fill=PatternFill('solid',fgColor='1F4E78');cell.font=Font(color='FFFFFF',bold=True)
            cell.alignment=Alignment(wrap_text=True)
        sh.freeze_panes='A2';sh.auto_filter.ref=sh.dimensions;sh.sheet_view.showGridLines=False
        for col,h in enumerate(fields,1):sh.column_dimensions[get_column_letter(col)].width=35 if any(s in h.lower() for s in ['evidence','hypotheses','rationale','confidence']) else 21
    add('All_hypotheses_audited',audited)
    add('Resolved_candidates',resolved)
    issues=[r for r in resolved if r['classification_state'] in {'ambiguous','unresolved','broad_only','provisional_fragment'}]
    add('Manual_review_queue',issues)
    if 'Pipeline_summary' in wb.sheetnames:del wb['Pipeline_summary']
    sh=wb.create_sheet('Pipeline_summary',0)
    sh.append(['Metric','Value'])
    metrics=[('Pipeline version',VERSION),('Rule schema',rules.get('schema_version')),('Raw hypotheses',len(raw)),('Distinct species/protein',len(resolved)),('Collapsed hypotheses',len(raw)-len(resolved)),('Rule file SHA-256',rules.get('_sha256','')),('Input SHA-256',checksum(source)),('Reference-supported',sum(x['classification_state']=='reference_supported' for x in resolved)),('Broad-only',sum(x['classification_state']=='broad_only' for x in resolved)),('Ambiguous',sum(x['classification_state']=='ambiguous' for x in resolved)),('Unresolved',sum(x['classification_state']=='unresolved' for x in resolved)),('Provisional fragments',sum(x['classification_state']=='provisional_fragment' for x in resolved)),('Subfamily annotations',sum(bool(x['assigned_subfamily']) for x in resolved))]
    for row in metrics:sh.append(row)
    sh.column_dimensions['A'].width=32;sh.column_dimensions['B'].width=85
    for cell in sh[1]:cell.fill=PatternFill('solid',fgColor='1F4E78');cell.font=Font(color='FFFFFF',bold=True)
    sh.freeze_panes='A2'
    wb.save(dest)
    manifest={'pipeline_version':VERSION,'rule_schema':rules.get('schema_version'),'input':source.name,'input_sha256':checksum(source),'rule_sha256':rules.get('_sha256'),'output':dest.name,'output_sha256':checksum(dest),'counts':dict(Counter(x['classification_state'] for x in resolved)),'unique_proteins':len(resolved),'input_hypotheses':len(raw)}
    dest.with_suffix('.manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    return manifest


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input',type=Path);parser.add_argument('output',type=Path)
    parser.add_argument('--sheet',default='candidates')
    parser.add_argument('--keep-all-audited',action='store_true',help='Compatibility flag: all audited hypotheses are always preserved')
    parser.add_argument('--rules',type=Path,default=Path(__file__).resolve().parents[1]/'01_detection/config/family_rules.json')
    parser.add_argument('--version',action='version',version='adhesome-audit '+VERSION)
    args=parser.parse_args()
    rules=json.loads(args.rules.read_text());rules['_sha256']=checksum(args.rules)
    meta=run_xlsx(args.input,args.output,args.sheet,rules)
    print(json.dumps(meta,indent=2))

if __name__=='__main__':main()
