"""Observed Pfam/CDD architecture and scoped putative protein naming (version 1.0).

CDD imports of Pfam models corroborate coordinates, not independent votes.
Inferred regions and target identifiers never supply family-specific evidence.
"""

# Resolve the organized bundle independently of the caller working directory.
import sys as _bundle_sys
from pathlib import Path as _BundlePath
_BUNDLE_ROOT=next(p for p in _BundlePath(__file__).resolve().parents if (p/'pipeline_paths.py').is_file())
if str(_BUNDLE_ROOT) not in _bundle_sys.path:_bundle_sys.path.insert(0,str(_BUNDLE_ROOT))
from pipeline_paths import activate as _activate_bundle_paths
_activate_bundle_paths()
from collections import Counter
from contextlib import contextmanager
from functools import lru_cache
import json
from pathlib import Path
import re

ROOT=Path(__file__).resolve().parents[1]
LIBRARY=ROOT/'annotation/family_identity_rules.json'
_CONFIGURATION=None


def observed_text(value):
    text=re.sub(r'\[\s*Inferred\s*:.*?\]','',str(value or ''),flags=re.I)
    return ';'.join(part for part in text.split(';') if not re.search(r'\binferred\b|not.detected|\babsent\b',part,re.I))


@lru_cache(maxsize=1)
def library():return _CONFIGURATION if _CONFIGURATION is not None else json.loads(LIBRARY.read_text(encoding='utf-8'))


@contextmanager
def using_configuration(configuration):
    """Temporary evidence-rule overlay; never rewrites the reference library."""
    global _CONFIGURATION
    previous=_CONFIGURATION
    try:
        _CONFIGURATION=configuration;library.cache_clear();yield
    finally:
        _CONFIGURATION=previous;library.cache_clear()


@contextmanager
def using_library(path=None):
    """Select one audit's identity library and restore the default afterwards."""
    global LIBRARY
    previous=LIBRARY
    selected=Path(path).resolve() if path is not None else previous
    # Validate readability before changing process-wide state.
    json.loads(selected.read_text(encoding='utf-8'))
    try:
        LIBRARY=selected;library.cache_clear()
        yield selected
    finally:
        LIBRARY=previous;library.cache_clear()


def cdd_hits(row):
    raw=row.get('NCBI_CDD_hits')
    if raw:
        records=json.loads(raw) if isinstance(raw,str) else raw
    else:
        text=observed_text(row.get('NCBI_CDD_architecture') or row.get('CDD_architecture'))
        records=[]
        for part in text.split(';'):
            match=re.fullmatch(r'\s*(.*?)\s*(?::\s*(\d+)\s*-\s*(\d+)|\[\s*(\d+)\s*-\s*(\d+)\s*\])\s*',part)
            if not match:continue
            name=re.sub(r'^(?:cd\d+|smart\d+|pfam\d+|COG\d+)\s+','',match[1].strip(),flags=re.I)
            records.append(dict(short_name=name,start=match[2] or match[4],end=match[3] or match[5],accession=''))
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
    pfam_extra={'PF07679':'IG','PF00047':'IG','PF13927':'IG','PF13895':'IG','PF06297':'PET'}
    for model,start,end in pfam_hits:
        # Broad architecture aliases must not erase the original observed model:
        # expanded rules may require its more specific Pfam accession. These
        # entries share a footprint, so distinct_footprints never counts an
        # accession and its alias as two domain copies within a model group.
        for label in dict.fromkeys((model,pfam_extra.get(model,model))):
            combined.append(dict(model=label,start=start,end=end,sources={'Pfam'},cdd_models=set()))
    for hit in cdd_hits(row):
        name=hit['short_name']
        models=[model for model,pattern in {**config['pfam_cdd_aliases'],**config['additional_cdd_domains']}.items() if re.search(pattern,name,re.I)]
        accession=re.fullmatch(r'(?:pfam|PF)(\d{5})(?:\.\d+)?',str(hit.get('accession') or ''),re.I)
        verified_model='PF'+accession[1] if accession and config.get('import_observed_pfam_accessions') else None
        if verified_model:models=list(dict.fromkeys([*models,verified_model]))
        for model in models:
            probe=dict(model=model,start=hit['start'],end=hit['end'],sources={'CDD'},cdd_models={name})
            if model==verified_model:probe['cdd_pfam_accession']=True
            matches=[old for old in combined if old['model']==model and overlaps(old,probe)]
            if matches:
                for old in matches:
                    old['sources'].add('CDD');old['cdd_models'].add(name)
                    if model==verified_model:old['cdd_pfam_accession']=True
            else:combined.append(probe)
    return sorted(combined,key=lambda hit:(hit['start'] is None,hit['start'] or 0,hit['model']))


def cdd_text(row):
    names=[hit['short_name'] for hit in cdd_hits(row)]
    # Coordinate-free returned annotations may support a label, never missing copies.
    text=observed_text(row.get('NCBI_CDD_architecture') or row.get('CDD_architecture'))
    return ';'.join(names+[text])


def observed_cdd_names(row):
    """Read model labels without treating workbook summary counts as coordinates."""
    names=[hit['short_name'] for hit in cdd_hits(row)]
    # Coordinate-free returned labels can corroborate identity but never add
    # missing domain copies. Exclude inferred/absent entries as elsewhere.
    for part in observed_text(row.get('NCBI_CDD_architecture') or row.get('CDD_architecture')).split(';'):
        name=re.sub(r'^(?:cd\d+|smart\d+|pfam\d+|COG\d+)\s+','',part.strip(),flags=re.I)
        # Located observations were already validated in cdd_hits. Do not
        # resurrect rejected out-of-sequence hits through their text labels.
        if re.search(r'(?::\s*\d+\s*-\s*\d+|\[\s*\d+\s*-\s*\d+\s*\])\s*$',name):continue
        name=re.sub(r'\s*(?::\s*\d+\s*-\s*\d+|\[.*?\])\s*$','',name).strip()
        if name:names.append(name)
    return list(dict.fromkeys(names))


def diagnostic_cdd(row,pattern):
    """Match each observed model separately, including anchored model patterns."""
    return any(re.search(pattern,name,re.I) for name in observed_cdd_names(row))


def distinct_footprints(hits,models):
    """Alternative domain models on one footprint count as a single copy."""
    accepted=[]
    # Selecting the shortest-ending footprint also avoids broad CDD superfamily
    # spans consuming several independently observed repeat copies.
    for hit in sorted((h for h in hits if h['model'] in models and h['start'] is not None),key=lambda h:(h['end'],h['start'])):
        if not any(overlaps(hit,old) for old in accepted):accepted.append(hit)
    return accepted


def model_counts(hits):
    return Counter({model:len(distinct_footprints(hits,{model})) for model in {h['model'] for h in hits}})


def ordered_core(hits,rule):
    last=0
    for index in rule.get('ordered_groups',[]):
        group=rule['required_groups'][index]
        choices=sorted(((hit['start'],hit['end']) for hit in hits if hit['model'] in group['models'] and hit['start'] is not None and hit['start']>last),key=lambda pair:(pair[1],pair[0]))
        if len(choices)<group['count']:return False
        # Select distinct nonoverlapping footprints, not overlapping alternative models.
        accepted=0
        for start,end in choices:
            if start<=last:continue
            last=end;accepted+=1
            if accepted==group['count']:break
        if accepted<group['count']:return False
    return True


def refinement_trace(row,parent,hits):
    """Explain each reviewed rule without adding proteins or relaxing its core.

    Cross-parent dispatch is opt-in per rule and requires both a diagnostic CDD
    model and a coordinate-bearing Pfam anchor belonging to that rule's core.
    Generic domains, gene names, and calibration labels cannot dispatch a rule.
    """
    counts=model_counts(hits);text=cdd_text(row);names=observed_cdd_names(row);traces=[]
    length=int(float(row.get('length') or row.get('Length') or 0))
    for rule in library()['rules']:
        allowed=rule['parents']
        parent_scope=parent in allowed or (parent=='fibronectin_like_candidate' and 'fibronectin_III_ligand' in allowed)
        core={model for group in rule['required_groups'] for model in group['models']}
        diagnostic=[name for name in names if re.search(rule['specific_cdd_pattern'],name,re.I)]
        anchor=any(hit['model'] in core and hit['start'] is not None and (
            'Pfam' in hit['sources'] or library().get('accept_cdd_pfam_anchor') and hit.get('cdd_pfam_accession')) for hit in hits)
        dispatched=bool(rule.get('allow_diagnostic_dispatch') and diagnostic and anchor)
        reached=parent_scope or dispatched
        groups=[dict(models=group['models'],required=group['count'],observed=len(distinct_footprints(hits,set(group['models'])))) for group in rule['required_groups']]
        missing=[group for group in groups if group['observed']<group['required']]
        order=ordered_core(hits,rule)
        coordinates=all(hit['start'] is not None and hit['start']>=1 and hit['end']>=hit['start'] and (not length or hit['end']<=length) for hit in hits if hit['model'] in core)
        length_ok=not rule.get('length_range') or rule['length_range'][0]<=length<=rule['length_range'][1]
        maximum=all(counts[model]<=limit for model,limit in rule.get('maximum_model_counts',{}).items())
        competitors=rule['id']=='actin_conventional' and bool(re.search(r'ARP\d',text,re.I))
        codes=[]
        if not reached:codes.append('NO_IDENTITY_RULE_REACHED')
        if not diagnostic:codes.append('DIAGNOSTIC_CDD_UNCONFIRMED')
        if not anchor:codes.append('PFAM_ANCHOR_MISSING')
        if missing:codes.append('MISSING_REQUIRED_DOMAIN')
        if not order:codes.append('DOMAIN_ORDER_CONFLICT')
        if not coordinates:codes.append('DOMAIN_COORDINATES_UNCONFIRMED')
        if not length_ok:codes.append('LENGTH_OUTSIDE_RULE_RANGE')
        if not maximum:codes.append('EXCESS_DOMAIN_COPIES')
        if competitors:codes.append('COMPETING_FAMILY')
        traces.append(dict(rule_id=rule['id'],family=rule['family'],parent=parent,parent_scope_passed=parent_scope,
            routing_scope='parent' if parent_scope else 'diagnostic_dispatch' if dispatched else 'out_of_scope',
            reached=reached,evaluated=reached,diagnostic_cdd_matched=bool(diagnostic),matched_cdd_models=diagnostic,
            pfam_anchor_passed=anchor,domain_groups=groups,missing_domains=missing,order_passed=order,
            coordinates_passed=coordinates,length_passed=length_ok,maximum_counts_passed=maximum,
            matched=not codes,limitation_codes=codes,reason_codes=codes))
    return traces


def scoped_refinement(row,parent,hits):
    trace=refinement_trace(row,parent,hits)
    matched={entry['rule_id'] for entry in trace if entry['matched']}
    matches=[rule for rule in library()['rules'] if rule['id'] in matched]
    # A named collagen subtype refines the broader fibrillar architecture.
    if any(rule['id']=='collagen_v' for rule in matches):matches=[r for r in matches if r['id']!='collagen_fibrillar']
    families={rule['family'] for rule in matches}
    return dict(rule=matches[0] if len(families)==1 else None,
        competing_families=sorted(families) if len(families)>1 else [],
        matches=matches,trace=trace)


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
        # Grade A already means the required family gate and orthology passed.
        # Optional CDD corroboration and experimental annotations are not extra gates.
        confidence='strong' if row['grade']=='A' else 'moderate' if row['grade']=='B' else 'provisional'
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
