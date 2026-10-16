"""Observed Pfam/CDD architecture and scoped putative protein naming (version 1.0).

CDD imports of Pfam models corroborate coordinates, not independent votes.
Inferred regions and target identifiers never supply family-specific evidence.
"""
from collections import Counter
from functools import lru_cache
import json
from pathlib import Path
import re

ROOT=Path(__file__).resolve().parent
LIBRARY=ROOT/'05_conservative_adhesome_family_assignment_audit/family_identity_rules.json'


def observed_text(value):
    text=re.sub(r'\[\s*Inferred\s*:.*?\]','',str(value or ''),flags=re.I)
    return ';'.join(part for part in text.split(';') if not re.search(r'\binferred\b|not.detected|\babsent\b',part,re.I))


@lru_cache(maxsize=1)
def library():return json.loads(LIBRARY.read_text(encoding='utf-8'))


def cdd_hits(row):
    raw=row.get('NCBI_CDD_hits')
    if raw:
        records=json.loads(raw) if isinstance(raw,str) else raw
    else:
        text=observed_text(row.get('NCBI_CDD_architecture') or row.get('CDD_architecture'))
        records=[dict(short_name=m[1].strip(),start=m[2],end=m[3],accession='') for m in
            re.finditer(r'([^;:]+):\s*(\d+)\s*-\s*(\d+)',text)]
    seen=set();result=[]
    for raw in records:
        item=dict(raw);name=str(item.get('short_name') or '')
        if not name or re.search('inferred|not.detected|absent',name,re.I):continue
        try:start,end=int(item['start']),int(item['end'])
        except (KeyError,TypeError,ValueError):continue
        if start<1 or end<start:continue
        length=int(float(row.get('length') or row.get('Length') or 0))
        if length and end>length:continue
        key=(name,start,end,item.get('accession',''))
        if key in seen:continue
        seen.add(key);item.update(short_name=name,start=start,end=end);result.append(item)
    return result


def overlaps(a,b):
    return a['start'] is not None and b['start'] is not None and min(a['end'],b['end'])-max(a['start'],b['start'])+1>=.5*min(a['end']-a['start']+1,b['end']-b['start']+1)


def combined_domains(row,pfam_hits):
    """Complement missing domain models using observed CDD footprints; merge copies."""
    config=library();combined=[]
    pfam_extra={'PF07679':'IG','PF00047':'IG','PF13927':'IG','PF13895':'IG'}
    for model,start,end in pfam_hits:
        combined.append(dict(model=pfam_extra.get(model,model),start=start,end=end,sources={'Pfam'},cdd_models=set()))
    for hit in cdd_hits(row):
        name=hit['short_name']
        models=[model for model,pattern in {**config['pfam_cdd_aliases'],**config['additional_cdd_domains']}.items() if re.search(pattern,name,re.I)]
        for model in models:
            probe=dict(model=model,start=hit['start'],end=hit['end'],sources={'CDD'},cdd_models={name})
            matches=[old for old in combined if old['model']==model and overlaps(old,probe)]
            if matches:
                for old in matches:old['sources'].add('CDD');old['cdd_models'].add(name)
            else:combined.append(probe)
    return sorted(combined,key=lambda hit:(hit['start'] is None,hit['start'] or 0,hit['model']))


def cdd_text(row):
    names=[hit['short_name'] for hit in cdd_hits(row)]
    # Coordinate-free returned annotations may support a label, never missing copies.
    text=observed_text(row.get('NCBI_CDD_architecture') or row.get('CDD_architecture'))
    return ';'.join(names+[text])


def ordered_core(hits,rule):
    last=0
    for index in rule.get('ordered_groups',[]):
        group=rule['required_groups'][index]
        choices=sorted((hit['start'],hit['end']) for hit in hits if hit['model'] in group['models'] and hit['start'] is not None and hit['start']>last)
        if len(choices)<group['count']:return False
        # Select distinct nonoverlapping footprints, not overlapping alternative models.
        accepted=0
        for start,end in choices:
            if start<=last:continue
            last=end;accepted+=1
            if accepted==group['count']:break
        if accepted<group['count']:return False
    return True


def scoped_refinement(row,parent,hits):
    counts=Counter(hit['model'] for hit in hits);text=cdd_text(row);matches=[]
    for rule in library()['rules']:
        allowed=rule['parents']
        if parent not in allowed and not (parent=='fibronectin_like_candidate' and 'fibronectin_III_ligand' in allowed):continue
        if not re.search(rule['specific_cdd_pattern'],text,re.I):continue
        if rule['id']=='actin_conventional' and re.search(r'ARP\d',text,re.I):continue
        if not all(sum(counts[model] for model in group['models'])>=group['count'] for group in rule['required_groups']):continue
        if not ordered_core(hits,rule):continue
        core={model for group in rule['required_groups'] for model in group['models']}
        if not any(hit['model'] in core and 'Pfam' in hit['sources'] for hit in hits):continue
        matches.append(rule)
    # A named collagen subtype refines the broader fibrillar architecture.
    if any(rule['id']=='collagen_v' for rule in matches):matches=[r for r in matches if r['id']!='collagen_fibrillar']
    families={rule['family'] for rule in matches}
    return dict(rule=matches[0] if len(families)==1 else None,
        competing_families=sorted(families) if len(families)>1 else [],
        matches=matches)


NAMES={'talin':'Talin','kindlin':'Kindlin-like protein','actin':'Actin-family protein',
    'cofilin':'Cofilin/ADF-family protein','vinculin':'Vinculin-family protein',
    'alpha_actinin':'Alpha-actinin-like protein','filamin':'Filamin-like protein',
    'paxillin_like':'Paxillin-like protein','zyxin_like':'Zyxin-like protein',
    'PDLIM_family':'PDZ-LIM family protein','PINCH_like':'PINCH-like protein','parvin_like':'Parvin-like protein',
    'FAK':'Focal adhesion kinase-like protein','ILK':'Integrin-linked kinase-like protein',
    'Src_Csk_like':'Src/Csk-related tyrosine kinase','Src':'Src-family tyrosine kinase',
    'Csk':'Csk/Chk-family tyrosine kinase','Syk':'Syk/ZAP70-family tyrosine kinase',
    'Nck':'Nck-like adaptor protein','Grb2':'Grb2-like adaptor protein','Crk':'Crk-like adaptor protein','Shc':'Shc-like adaptor protein',
    'integrin_alpha':'Integrin alpha subunit-like protein','integrin_beta':'Integrin beta subunit-like protein',
    'fibronectin_like':'Fibronectin-like protein','collagen_like':'Collagen-like protein',
    'laminin_like':'Laminin-like protein','vitronectin_like':'Vitronectin-like protein',
    'PTP_PEST':'PEST-type protein tyrosine phosphatase-like protein',
    'LAR_RPTP':'LAR-type receptor protein tyrosine phosphatase-like protein',
    'RPTP_alpha_epsilon':'Receptor protein tyrosine phosphatase alpha/epsilon-like protein'}


def annotate_name(row,result):
    family=row['assigned_family'];chosen=result.get('identity_rule')
    if family=='unassigned':
        models={hit['model'] for hit in result.get('combined_hits',[])}
        description=next((label for model,label in [('PF00041','fibronectin type III domain-containing protein'),
            ('PF00102','protein tyrosine phosphatase domain-containing protein'),('PF00069','protein kinase domain-containing protein'),
            ('PF07714','protein kinase domain-containing protein'),('PF00412','LIM domain-containing protein'),
            ('PF00307','calponin-homology domain-containing protein')] if model in models),'uncharacterized protein')
        name='Putative '+description;confidence='unresolved';scope='domain_only' if models else 'unresolved'
    else:
        label=chosen['protein_name'] if chosen else NAMES.get(family,family.replace('_',' ')+' family protein')
        name='Putative '+label
        concordant=result.get('cdd_concordant',False)
        confidence='strong' if row['grade']=='A' and concordant else 'moderate' if row['grade'] in ('A','B') else 'provisional'
        scope='family_lineage' if chosen else 'family'
    conflicts=result.get('identity_competitors',[])
    if conflicts:confidence='ambiguous';scope='ambiguous_family';name='Putative '+NAMES.get(family,'unresolved family protein')+' (competing CDD identities)'
    hits=result.get('combined_hits',[])
    row.update(protein_name=name,family_confidence=confidence,protein_identity_confidence=confidence,
        protein_name_scope=scope,identity_rule=chosen['id'] if chosen else 'configured_family_architecture',
        identity_competitors='; '.join(conflicts),
        observed_domain_architecture='; '.join(f"{h['model']}:{h['start']}-{h['end']}[{'+'.join(sorted(h['sources']))}]" for h in hits),
        naming_evidence=('Observed diagnostic CDD and compatible ordered domain core' if chosen else 'Configured family architecture')+
            ('; positive orthology support' if row['orthology_score_0_3']>0 else '; orthology support unavailable/insufficient'),
        probability_status='not_calibrated_for_family_or_protein_identity')
    return row
