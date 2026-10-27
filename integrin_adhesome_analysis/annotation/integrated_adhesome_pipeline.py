#!/usr/bin/env python3
"""Version 1.0: architecture-gated refinement of stage 1 family hypotheses."""

# Resolve the organized bundle independently of the caller working directory.
import sys as _bundle_sys
from pathlib import Path as _BundlePath
_BUNDLE_ROOT=next(p for p in _BundlePath(__file__).resolve().parents if (p/'pipeline_paths.py').is_file())
if str(_BUNDLE_ROOT) not in _bundle_sys.path:_bundle_sys.path.insert(0,str(_BUNDLE_ROOT))
from pipeline_paths import activate as _activate_bundle_paths
_activate_bundle_paths()
import argparse
import csv
from collections import Counter,defaultdict
import hashlib
import json
import math
from pathlib import Path
import re
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from classification_contract import general_family,general_rules,write_sheet,schema
from family_identity import combined_domains,scoped_refinement,cdd_text,annotate_name,model_counts,LIBRARY as IDENTITY_LIBRARY
import family_identity
from orthogroup_context import validate_weight,configuration as context_configuration
from stage2_resolution import broad_architecture,annotate_resolution,common_broad_family,BROAD_LABELS

VERSION='1.0'
DOM_RE=re.compile(r'(?<!\w)(PF\d{5})(?:\s*:\s*(\d+)\s*-\s*(\d+))?',re.I)
TRUTH={'1','TRUE','YES','Y','PASS','SUPPORTED','CONSISTENT'}


def value(row,*keys):
    for key in keys:
        if row.get(key) is not None and str(row[key]).strip():return str(row[key]).strip()
    return ''


def num(value,default=0.0):
    try:
        result=float(str(value).replace(',',''))
        return result if math.isfinite(result) else default
    except (TypeError,ValueError):return default


def boolean(value):return str(value).strip().upper() in TRUTH
def norm(value):return re.sub('[^A-Z0-9]','',str(value).upper())
def checksum(path):
    with Path(path).open('rb') as handle:return hashlib.file_digest(handle,'sha256').hexdigest()


def identity(row,rule):
    """Compatibility API for optional detection evidence; never bypasses stage 2 gates."""
    target=norm(value(row,'reference_subfamily','orthology_resolved_family'))
    labels={norm(x) for x in rule.get('reference_labels',[])}|{norm(rule.get('label',''))}
    if not target or target not in labels:return 0,'no independent identity label'
    method=norm(value(row,'subfamily_support_method'))
    passed=(method=='CURATEDREFERENCE' and bool(value(row,'reference_accession')) or
        method=='RECIPROCALBESTHIT' and norm(value(row,'reciprocal_reference_subfamily'))==target or
        method=='PHYLOGENY' and boolean(value(row,'phylogeny_concordant')) or
        method=='CALIBRATEDHMM' and boolean(value(row,'subfamily_profile_calibrated')) and boolean(value(row,'profile_margin_pass')) and norm(value(row,'top_profile_subfamily'))==target)
    return (3,'independent identity annotation') if passed else (1,'unresolved identity annotation')


def policy_for(rules):
    return rules.get('_audit_policy') or json.loads(Path(__file__).with_name('audit_policy.json').read_text())


def domains(row):
    hits=[]
    for match in DOM_RE.finditer(value(row,'architecture','Pfam_architecture','pfam_architecture')):
        pf=match[1].upper();start=int(match[2]) if match[2] else None;end=int(match[3]) if match[3] else None
        if start is not None:
            if start<1 or end<start:continue
            if any(pf==old[0] and old[1] is not None and min(end,old[2])-max(start,old[1])+1>=.5*min(end-start+1,old[2]-old[1]+1) for old in hits):continue
        hits.append((pf,start,end))
    return hits


def position_order(hits,groups):
    last=0;used=set()
    for group in groups:
        choices=[group] if isinstance(group,str) else group
        found=sorted((start,end,i) for i,(pf,start,end) in enumerate(hits) if pf in choices and start is not None and start>last and i not in used)
        if not found:return False
        start,end,i=found[0];last=end;used.add(i)
    return True


def cdd_support(row,family,policy):
    if family_identity.library().get('precision_safeguards'):
        pattern=policy['specific_cdd_patterns'].get(family)
        return bool(pattern and family_identity.diagnostic_cdd(row,pattern))
    text=re.sub(r'[_:/;(),-]+',' ',cdd_text(row).upper())
    pattern=policy['specific_cdd_patterns'].get(family)
    if not pattern:return False
    return bool(re.search(r'(?<![A-Z0-9])(?:'+pattern+r')(?![A-Z0-9])',text,re.I))


def _architecture(row,family,rules):
    """Only the supplied hypothesis is tested; no catalogue-wide replacement search."""
    policy=policy_for(rules);catalogue=rules['families'];general=general_family(family,catalogue)
    aliases={'Csk_family':'Csk','Src_family':'Src','paxillin_family':'paxillin_like','zyxin_family':'zyxin_like'}
    general=aliases.get(general,general)
    rule=dict(general_rules(rules).get(general,catalogue.get(general,catalogue.get(family,{}))))
    combined=combined_domains(row,domains(row))
    hits=[(hit['model'],hit['start'],hit['end']) for hit in combined]
    refinement=scoped_refinement(row,general,combined)
    if refinement['rule']:
        named=refinement['rule'];core={model for group in named['required_groups'] for model in group['models']}
        length=num(value(row,'length','Length'))
        valid=all(h['start'] is not None and (not length or h['end']<=length) for h in combined if h['model'] in core)
        if valid:
            return dict(score=4,fraction=1,family=named['family'],note='Observed Pfam/CDD core and diagnostic CDD support scoped '+named['id']+' refinement',
                cdd=True,complete=True,identity_rule=named,combined_hits=combined,cdd_concordant=True,identity_competitors=[],
                rule_trace=refinement.get('trace',[]),baseline_gate=dict(configured=True,complete=True,reason_codes=[]))
    if rule.get('screening_only') or rule.get('classification_level') in ('screen','screening','component_hypothesis'):
        return dict(score=0,fraction=0,family=general,note='Generic screening bucket cannot establish a family assignment',cdd=False,complete=False,
                    combined_hits=combined,identity_competitors=refinement['competing_families'],rule_trace=refinement.get('trace',[]),baseline_gate=dict(configured=False,complete=False,reason_codes=['NO_IDENTITY_RULE_REACHED']))
    counts=model_counts(combined)
    length=num(value(row,'length','Length'))
    specific=general in policy['specific_cdd_patterns']
    if general=='PINCH_like':
        rule.update(domain_alternatives=[{'PF00412':5}],ordered_domain_groups=[['PF00412']]*5,min_length=250)
    elif general=='parvin_like':
        rule.update(domain_alternatives=[{'PF00307':2}],ordered_domain_groups=[['PF00307']]*2,min_length=0)
    elif general=='Syk':
        rule.update(domain_alternatives=[{'PF00017':2,'PF07714':1},{'PF00017':2,'PF00069':1}],ordered_domain_groups=[['PF00017'],['PF00017'],['PF07714','PF00069']])
    elif general in ('Csk','Src'):
        rule.update(domain_alternatives=[{'PF00018':1,'PF00017':1,'PF07714':1},{'PF00018':1,'PF00017':1,'PF00069':1}],ordered_domain_groups=[['PF00018'],['PF00017'],['PF07714','PF00069']])
    elif general=='alpha_actinin':
        rule.update(domain_alternatives=[{'PF00307':2,'PF00435':2}])
    requirements=rule.get('domain_alternatives',[])
    if not requirements:return dict(score=0,fraction=0,family=general,note='No family-specific architecture configured for this stage 1 hypothesis',cdd=False,complete=False,
        combined_hits=combined,identity_competitors=refinement['competing_families'],rule_trace=refinement.get('trace',[]),baseline_gate=dict(configured=False,complete=False,reason_codes=['NO_IDENTITY_RULE_REACHED']))
    fractions=[sum(min(counts[p],n) for p,n in group.items())/sum(group.values()) for group in requirements]
    fraction=max(fractions)
    if fraction==0:return dict(score=0,fraction=0,family=general,note='Required family domain anchors absent',cdd=False,complete=False,
        combined_hits=combined,identity_competitors=refinement['competing_families'],rule_trace=refinement.get('trace',[]),baseline_gate=dict(configured=True,complete=False,reason_codes=['MISSING_REQUIRED_DOMAIN']))
    complete=fraction>=1
    groups=rule.get('ordered_domain_groups') or rule.get('ordered_domains')
    alternatives=rule.get('support_order_alternatives') or ([groups] if groups else [])
    ordered=any(position_order(hits,order) for order in alternatives) if alternatives else True
    valid_coordinates=all(start is not None and (not length or end<=length) for pf,start,end in hits if any(pf in group for group in requirements))
    length_ok=not rule.get('min_length') or length>=rule['min_length']
    if general=='PINCH_like':length_ok=policy['pinch_length_range_inclusive'][0]<=length<=policy['pinch_length_range_inclusive'][1]
    if general=='PINCH_like' and counts['PF00412']!=5:complete=False
    if general=='alpha_actinin' and (counts['PF00435']>4 or length>1500):complete=False
    cdd=cdd_support(row,general,policy)
    score=4 if complete and ordered and valid_coordinates and length_ok and (not specific or cdd) else 3 if complete else 2 if fraction>=.5 else 1
    strict=family_identity.library().get('precision_safeguards',False)
    proposal=strict and (rule.get('no_auto_promotion') or rule.get('user_inspection_required'))
    shared_only=strict and not cdd and all(p in {'PF00069','PF07714','PF00373','PF00412','PF00041','PF00169','PF00018','PF00017'} for group in requirements for p in group)
    if proposal or shared_only:score=min(score,3)
    notes=[f'family domain recovery={fraction:.2f}']
    if not ordered:notes.append('required domain order unconfirmed')
    if not valid_coordinates:notes.append('required domain coordinates unavailable or invalid')
    if not length_ok:notes.append('outside configured family length range')
    if specific:notes.append('family-specific CDD present' if cdd else 'required family-specific CDD absent')
    if proposal:notes.append('supplied family proposal requires review; no automatic supported promotion')
    if shared_only:notes.append('shared domain core alone cannot establish this specific family')
    if refinement['competing_families']:
        score=min(score,3);notes.append('Competing diagnostic CDD identities: '+', '.join(refinement['competing_families']))
    call=general
    if general=='Src':
        if cdd_support(row,'Csk',policy):call='Csk';cdd=True;notes.append('Csk/Chk-specific CDD refines the SH3-SH2-kinase hypothesis')
        else:call='Src_Csk_like';notes.append('compatible Src/Csk core; no Csk-specific call')
    elif general=='Csk' and not cdd:call='Src_Csk_like'
    # Baseline support is retained as evidence, not as an exemption from the gate.
    if value(row,'status') in ('supported','supported_candidate'):
        notes.append('stage 1 supported baseline rechecked against domain types/copies')
    core_models={model for group in requirements for model in group}
    concordant=any(h['model'] in core_models and 'CDD' in h['sources'] for h in combined)
    reasons=[]
    if not complete:reasons.append('MISSING_REQUIRED_DOMAIN')
    if not ordered:reasons.append('DOMAIN_ORDER_CONFLICT')
    if not valid_coordinates:reasons.append('DOMAIN_COORDINATES_UNCONFIRMED')
    if not length_ok:reasons.append('LENGTH_OUTSIDE_RULE_RANGE')
    if specific and not cdd:reasons.append('DIAGNOSTIC_CDD_UNCONFIRMED')
    if proposal:reasons.append('RULE_PROPOSAL_REVIEW_REQUIRED')
    if shared_only:reasons.append('SHARED_CORE_NOT_FAMILY_DIAGNOSTIC')
    if not complete and (str(row.get('sequence_completeness','')).lower() in ('partial','fragment','fragmented') or
            any(str(h.get('incomplete','')).upper() in ('N','C','NC') for h in family_identity.cdd_hits(row))):
        reasons.append('MISSING_REGION_UNASSESSABLE_PARTIAL_SEQUENCE')
    if refinement['competing_families']:reasons.append('COMPETING_FAMILY')
    return dict(score=score,fraction=fraction,family=call,note='; '.join(notes),cdd=cdd,complete=score==4,
        combined_hits=combined,cdd_concordant=cdd or concordant,identity_competitors=refinement['competing_families'],
        rule_trace=refinement.get('trace',[]),baseline_gate=dict(configured=True,complete=score==4,required_domain_groups=requirements,
            observed_counts=dict(counts),ordered=ordered,coordinates_passed=valid_coordinates,length_passed=length_ok,
            diagnostic_cdd_required=specific,diagnostic_cdd_passed=cdd,reason_codes=reasons))


def architecture(row,family,rules):
    from stage2_evidence import partial_confirmation
    return partial_confirmation(row,family,broad_architecture(row,_architecture(row,family,rules)),rules)


def benchmark_weight(row,family=None):
    if value(row,'Orthology_benchmark_status')!='validated':return 1.0,'unavailable_neutral_weight'
    if value(row,'Orthology_benchmark_scope')=='adhesome_subset_reference_orthogroup_inference':
        try:
            weights=json.loads(value(row,'Orthology_benchmark_family_weights') or '{}')
            if not isinstance(weights,dict):raise ValueError('Family weights must be a mapping')
            family=family or value(row,'family')
            if family not in weights:return 1.0,'family_not_benchmarked_neutral_weight'
            weight=float(weights[family])
            return (weight,'validated_adhesome_family_weight') if math.isfinite(weight) and 0<=weight<=1 else (1.0,'invalid_weight_ignored')
        except (ValueError,TypeError):
            return 1.0,'invalid_weight_ignored'
    weight=num(value(row,'Orthology_benchmark_weight'),-1)
    return (weight,'validated_reference_cohort_weight') if 0<=weight<=1 else (1.0,'invalid_weight_ignored')


def orthology_score(row,policy,family=None,context_weight=None):
    references=set()
    for key,item in row.items():
        if key.startswith('Orthology_') and key.endswith('_HOG_copy_count') and num(item)>0:references.add(key[len('Orthology_'):-len('_HOG_copy_count')])
    targets=set(value(row,'Orthology_target_species').split(';'))-{''}
    references-=targets
    count=max(len(references),int(num(value(row,'Orthology_HOG_reference_species_count','Orthology_direct_reference_count'))))
    if count==0 and value(row,'Orthology_reference_support'):
        count=len({x.strip() for x in re.split('[;,]',value(row,'Orthology_reference_support')) if x.strip()})
    score=next((score for minimum,score in policy['orthology_reference_species_bands'] if count>=minimum),0)
    context_weight=validate_weight(context_weight)
    weight,status=benchmark_weight(row,family) if context_weight is None else (context_weight,'heuristic_orthogroup_context' if context_weight else 'orthogroup_context_disabled')
    return score,round(score*weight,6),count,weight,status


def motif_score(row):
    if num(value(row,'MotifPrediction_context_supported_count'))>0:return 2
    tiers=value(row,'MotifPrediction_candidate_tiers','motif_candidate_tiers').lower()
    if 'context_supported' in tiers and 'context_conflict' not in tiers:return 2
    if num(value(row,'MotifPrediction_sequence_only_count'))>0:return 1
    if 'sequence_only' in tiers:return 1
    if 'context_conflict' not in tiers and num(value(row,'MotifPrediction_priority_site_count'))>0:return 1
    return 0


def localization(row,family,rules):
    rule=rules['families'].get(family,{})
    tm=value(row,'DeepTMHMM').upper();sp=value(row,'SignalP.v6.0').upper();loc=value(row,'DeepLoc_2.1').upper()
    secreted=bool(re.search(r'\bSP\b',sp+' '+tm)) and 'NO SP' not in sp+' '+tm
    trans=bool(re.search(r'\bTM\b|TRANSMEMBRANE',tm)) and 'NO TM' not in tm
    if rule.get('module')=='integrin_receptors':return (1,'compatible SP and transmembrane topology') if secreted and trans else (0,'integrin topology unconfirmed')
    if rule.get('module')=='extracellular_matrix':return (1,'compatible extracellular routing') if secreted or 'EXTRACELLULAR' in loc else (0,'extracellular routing unconfirmed')
    if trans or secreted:return 0,'predicted secretory/TM context conflicts with typical intracellular family context'
    if any(x in loc for x in ('CYTOPLASM','CYTOSOL','CELL MEMBRANE','NUCLEUS')) or tm in ('GLOB','GLOBULAR','SOLUBLE','NO TM'):
        return 1,'compatible intracellular context'
    return 0,'localization not informative'


def row_species(row):
    explicit=value(row,'Specie','species','Species')
    if explicit:return explicit
    identifier=value(row,'protein_id','candidate_id')
    return identifier.split('__',1)[0] if '__' in identifier else 'species_unknown'


def unique_key(row,index=0):
    pid=value(row,'protein_id','candidate_id','Protein_ID')
    if not pid:raise ValueError('Stage 1 hypothesis is missing protein_id')
    return row_species(row),pid


def audit_record(row,rules):
    policy=policy_for(rules);baseline=value(row,'family')
    result=architecture(row,baseline,rules)
    context_weight=validate_weight(rules.get('_orthogroup_context_weight'))
    raw,orth,ref_count,weight,weight_status=orthology_score(row,policy,result['family'],context_weight)
    motif=motif_score(row);loc,loc_note=localization(row,result['family'],rules)
    domain=result['score'];supported_orthology=orth>0
    grade='A' if domain==4 and supported_orthology else 'B' if domain==4 else 'C' if domain>0 and supported_orthology else 'D'
    tested=boolean(value(row,'Orthology_analysis_completed')) or value(row,'Orthology_HOG_status').lower() in ('shared hog','target-only in sampled hog')
    if grade=='A':decision='family_supported';assignment='assigned'
    elif grade=='B':decision='target_panel_restricted' if tested and ref_count==0 else 'orthology_unavailable' if not tested else 'orthology_support_insufficient';assignment='assigned'
    elif grade=='C':decision='provisional_family';assignment='provisional'
    else:decision='architecture_gate_failed' if domain==0 else 'partial_architecture_without_orthology';assignment='unassigned'
    assigned=result['family'] if grade in ('A','B','C') else 'unassigned'
    out=dict(row)
    out['_precision_safeguards']=family_identity.library().get('precision_safeguards',False)
    out['_competition_review']=family_identity.library().get('identity_ablation',{}).get('competition') is not False
    out.update(preaudit_family=baseline,audit_family=result['family'],assigned_family=assigned,
        score_interpretation='Evidence summary; not an identity probability or proof of physical adhesome membership',
        functional_association_status='computational_family_and_context_candidate',
        stage1_assigned_family=value(row,'assigned_family'),Specie=row_species(row),architecture_evidence=result['note'],
        domain_score_0_4=domain,orthology_raw_score_0_3=raw,orthology_score_0_3=orth,motif_score_0_2=motif,
        localization_score_0_1=loc,total_score_0_10=round(domain+orth+motif+loc,6),grade=grade,family_decision=decision,
        assignment_status=assignment,classification_state='family_supported' if grade in ('A','B') else 'provisional_family' if grade=='C' else 'unresolved',
        domain_gate_passed=domain==4,family_specific_CDD_support=result['cdd'],localization_evidence=loc_note,
        orthology_reference_species_count=ref_count,
        orthology_benchmark_weight_used=weight if context_weight is None else '',
        orthology_benchmark_weight_status=weight_status if context_weight is None else 'not_applied_heuristic_context_override',
        orthogroup_context_weight=weight,orthogroup_context_policy='baseline' if context_weight is None else 'disabled' if weight==0 else 'heuristic',
        restriction_scope='sampled target/reference panel' if decision=='target_panel_restricted' else '',
        discovery_rule=value(row,'discovery_rule') or 'stage1_hypothesis',pipeline_version=VERSION,
        rationale=result['note']+f'; reference species={ref_count}; orthology={orth:g}; motif={motif}; localization={loc}; '+decision)
    if context_weight is None and value(row,'Orthology_benchmark_scope')=='adhesome_subset_reference_orthogroup_inference':
        out['Orthology_benchmark_weight']=weight if weight_status=='validated_adhesome_family_weight' else ''
        if weight_status=='family_not_benchmarked_neutral_weight':out['Orthology_benchmark_status']='family_not_benchmarked'
    if 'combined_hits' not in result:result['combined_hits']=combined_domains(row,domains(row))
    out['_rule_gate_trace']=result.get('rule_trace',[])
    return annotate_resolution(annotate_name(out,result),result)


def rank(row):
    return ({'assigned':3,'provisional':2,'unassigned':0}.get(row['assignment_status'],0),
        ({'family_supported':3,'target_panel_restricted':2,'orthology_support_insufficient':2,'orthology_unavailable':1,'provisional_family':1}.get(row['family_decision'],0),int(row.get('family_specific_CDD_support',False))),
        {'A':4,'B':3,'C':2,'D':1}[row['grade']],row['domain_score_0_4'],
        2 if row.get('family_support_level')=='specific_supported' else 1 if row.get('family_support_level')=='broad_supported' else 0,
        int(row.get('row_origin')=='stage2_recovery'),row['orthology_score_0_3'],row['motif_score_0_2'],
        row['localization_score_0_1'],row['total_score_0_10'],{'supported_candidate':3,'supported':3,'provisional':2,'conflicting':1}.get(value(row,'status'),0))


def collapse_rows(rows):
    grouped=defaultdict(list)
    for row in rows:grouped[unique_key(row)].append(row)
    resolved=[]
    for key,group in grouped.items():
        ordered=sorted(group,key=rank,reverse=True);chosen=dict(ordered[0])
        initial=[r for r in group if r.get('row_origin','original_stage1')!='stage2_recovery']
        chosen.update(alternative_family_hypotheses='; '.join(sorted({value(r,'family') for r in initial}-{chosen['preaudit_family']})),
            collapsed_row_count=len(group),original_row_count=len(initial),recovered_hypothesis_count=len(group)-len(initial),
            equally_ranked_hypotheses='; '.join(sorted({r['audit_family'] for r in ordered if rank(r)==rank(ordered[0])})))
        tied={r['assigned_family'] for r in ordered if rank(r)==rank(ordered[0])}-{'unassigned'}
        if any(r.get('_precision_safeguards') and r.get('_competition_review',True) for r in ordered):
            tied|={r['assigned_family'] for r in ordered if r.get('domain_gate_passed') and r.get('family_specific_CDD_support') and r['assigned_family']!='unassigned'}
        if len(tied)>1:
            broad=common_broad_family(tied) if all(r['domain_score_0_4']==4 for r in ordered if rank(r)==rank(ordered[0])) else None
            if broad:
                chosen.update(assigned_family=broad,audit_family=broad,family_support_level='broad_supported',
                    annotation_resolution='broad_family',protein_name_scope='broad_family',
                    protein_name='Putative '+BROAD_LABELS[broad],protein_identity_confidence='supported_at_broad_family_scope',
                    specific_identity_status='specific_lineage_or_paralogue_unresolved',unresolved_specific_families=';'.join(sorted(tied)))
                trace=json.loads(chosen.get('identity_gate_trace') or '{}')
                trace.update(tested_family=broad,unresolved_specific_families=sorted(tied),
                             reason_codes=sorted(set(trace.get('reason_codes',[]))|{'SPECIFIC_IDENTITY_UNRESOLVED'}))
                chosen['identity_gate_trace']=json.dumps(trace,sort_keys=True)
            else:
                chosen.update(family_confidence='ambiguous',protein_identity_confidence='ambiguous',protein_name_scope='ambiguous_family',
                    family_support_level='ambiguous',annotation_resolution='ambiguous_family',
                    protein_name='Putative unresolved '+' / '.join(sorted(tied))+' family protein',
                    naming_evidence=chosen['naming_evidence']+'; equally ranked competing families remain unresolved')
                trace=json.loads(chosen.get('identity_gate_trace') or '{}')
                trace['reason_codes']=sorted(set(trace.get('reason_codes',[]))|{'COMPETING_FAMILY'})
                chosen['identity_gate_trace']=json.dumps(trace,sort_keys=True)
        if chosen.get('family_confidence')=='ambiguous':chosen['identity_outcome']='conflicting_identity'
        resolved.append(chosen)
    return sorted(resolved,key=unique_key)


def classify_rows(rows,rules):
    from stage2_evidence import annotation_context
    with annotation_context(rules):return _classify_rows(rows,rules)


def _classify_rows(rows,rules):
    grouped=defaultdict(list);original=[]
    for index,raw in enumerate(rows,1):
        row=dict(raw,original_hypothesis_row=index,row_origin='original_stage1')
        grouped[unique_key(row)].append(row);original.append(row)
    recovered=[]
    for key,group in grouped.items():
        hashes={value(row,'sequence_sha256') for row in group}-{''}
        if len(hashes)>1:raise ValueError('Conflicting sequences within stage 1 species/protein group: '+str(key))
        cdd=';'.join(sorted({value(row,'NCBI_CDD_architecture') for row in group}-{''}))
        structured={}
        for row in group:
            raw=row.get('NCBI_CDD_hits') or '[]'
            try:records=json.loads(raw) if isinstance(raw,str) else raw
            except json.JSONDecodeError as exc:
                raise ValueError('CDD hit JSON is incomplete (Excel cells are limited to 32767 characters); supply the full --supporting-data sidecar') from exc
            for hit in records:
                structured[(hit.get('accession'),hit.get('short_name'),str(hit.get('start')),str(hit.get('end')))]=hit
        for row in group:
            row['NCBI_CDD_architecture']=cdd
            if structured:row['NCBI_CDD_hits']=json.dumps(list(structured.values()),sort_keys=True)
        existing={general_family(value(row,'family'),rules['families']) for row in group}
        from stage2_evidence import candidate_hypotheses
        recovered.extend(candidate_hypotheses(group[0],existing,rules,domains))
        for family in policy_for(rules)['recoveries']:
            if family in existing:continue
            probe=dict(group[0],family=family,assigned_family='',row_origin='stage2_recovery',
                discovery_rule='recovered_'+family+'_complete_architecture',recovered_from_families='; '.join(sorted({value(r,'family') for r in group})))
            if architecture(probe,family,rules)['score']==4:recovered.append(probe)
    audited=[audit_record(row,rules) for row in original+recovered]
    resolved=collapse_rows(audited)
    if {unique_key(r) for r in resolved}!=set(grouped):raise ValueError('Stage 2 changed the stage 1 protein pool')
    from control_exclusion import apply as exclude_controls
    resolved,excluded=exclude_controls(resolved)
    exclusions={unique_key(row):row for row in excluded}
    for row in audited:
        match=exclusions.get(unique_key(row))
        if match:
            for key in ('stage2_excluded','exclusion_gate','exclusion_gates','control_class','control_accession','exclusion_reason','excluded_putative_family','excluded_protein_name'):
                row[key]=match[key]
    if {unique_key(row) for row in resolved}|set(exclusions)!=set(grouped):raise ValueError('Retained and excluded protein accounting is incomplete')
    from identity_probability import annotate as calibrate_probabilities
    calibrate_probabilities(resolved,family_rules_sha256=rules.get('_sha256'),orthogroup_context_weight=rules.get('_orthogroup_context_weight'))
    return audited,resolved


def run_xlsx(source,dest,sheet_name,rules,supporting=None,identity_library=None):
    with family_identity.using_library(identity_library):
        return _run_xlsx(source,dest,sheet_name,rules,supporting)


def _run_xlsx(source,dest,sheet_name,rules,supporting=None):
    from openpyxl import load_workbook,Workbook
    if source.resolve()==dest.resolve():raise ValueError('Choose a new output workbook')
    input_book=load_workbook(source,read_only=True,data_only=True)
    names={s.lower():s for s in input_book.sheetnames}
    selected=names.get(sheet_name.lower())
    if not selected:raise ValueError('Input hypothesis sheet not found: '+sheet_name)
    values=input_book[selected].values;headings=next(values)
    if len(set(headings))!=len(headings):raise ValueError('Duplicate input column names')
    rows=[dict(zip(headings,row)) for row in values if any(item is not None for item in row)]
    input_book.close()
    supporting=supporting or {}
    cdd_by_id=defaultdict(dict)
    for hit in supporting.get('tables',{}).get('NCBI_CDD_selected_hits',[]):
        cdd_by_id[hit['sequence_id']][(hit.get('accession'),hit.get('short_name'),str(hit['start']),str(hit['end']))]=hit
    for row in rows:
        if row['protein_id'] in cdd_by_id:row['NCBI_CDD_hits']=list(cdd_by_id[row['protein_id']].values())
    audited,resolved=classify_rows(rows,rules)
    from phylogenetic_identity import apply as phylogenetic_resolution
    phylogenetic_status=phylogenetic_resolution(resolved,rules.get('_phylogenetic_identity_evidence'),
        rules.get('_phylogenetic_identity_resolution',False))
    from family_catalogue import registry_path,report as catalogue_report
    from family_set_review_gate import gate as family_review
    from cdd_evidence import interpretation as cdd_interpretation
    registry=json.loads(registry_path().read_text());review=[];cdd_audit=[]
    for row in resolved:
        hits=family_identity.cdd_hits(row)
        specifics=[h for h in hits if cdd_interpretation(h)['specific']]
        for hit in hits:cdd_audit.append({**hit,**cdd_interpretation(hit),'protein_id':row['protein_id'],'Specie':row['Specie']})
        matches=[r for r in family_identity.library()['rules'] if r['id']==row.get('identity_rule')]
        family_set=matches[0].get('family_set',row['assigned_family']) if matches else row['assigned_family']
        assessment=family_review(dict(row,family_set=family_set,
            stage2_architecture_compatible=row['domain_gate_passed'],stage2_cdd_pattern_match=row.get('family_specific_CDD_support',False),
            stage2_named_curated_cd_match=any(cdd_interpretation(h)['curated_model'] for h in specifics),
            ncbi_cdd_specific_hit_threshold_met=bool(specifics)),registry,
            rules.get('_phylogenetic_identity_resolution',False))
        assessment['registry_policy']='Coverage/QC metadata only; historical calibration inspection flags do not restrict the main calibrated_v1 classifier'
        from stage2_evidence import profile as selected_profile
        if selected_profile(rules)[1].get('precision_safeguards'):
            for key in ('family_set_acceptance_status','family_set_review_flags','user_inspection_required'):
                row[key]=assessment[key]
        review.append(assessment)
    def sidecar(name,records,empty=('protein_id','status')):
        fields=list(dict.fromkeys(k for r in records for k in r if not k.startswith('_'))) or list(empty)
        with (dest.parent/name).open('w',encoding='utf-8',newline='') as h:
            writer=csv.DictWriter(h,fieldnames=fields,delimiter='\t',extrasaction='ignore');writer.writeheader()
            writer.writerows({k:json.dumps(v,sort_keys=True) if isinstance(v,(dict,list)) else v for k,v in r.items() if not k.startswith('_')} for r in records)
    sidecar('family_set_review.tsv',review)
    sidecar('cdd_specificity_audit.tsv',cdd_audit)
    decision_fields=('Specie','protein_id','preaudit_family','assigned_family','protein_name','identity_outcome',
        'supported_identity_level','family_support_level','annotation_resolution','grade','family_adjudication',
        'identity_gate_trace','phylogenetic_identity_status','phylogenetic_identity_reason')
    sidecar('annotation_decisions.tsv',[{k:r.get(k,'') for k in decision_fields} for r in resolved])
    (dest.parent/'family_configuration_reconciliation.json').write_text(json.dumps(catalogue_report(),indent=2)+'\n')
    # Detailed rule explanations remain in TSV, never in an Excel cell whose
    # 32767-character limit could truncate the structured evidence.
    trace_fields=['Specie','protein_id','stage1_family','rule_id','family','reached','routing_scope',
                  'diagnostic_cdd_matched','matched_cdd_models','pfam_anchor_passed','domain_groups','missing_domains',
                  'order_passed','coordinates_passed','length_passed','maximum_counts_passed','matched','reason_codes']
    with (dest.parent/'stage2_rule_gate_trace.tsv').open('w',encoding='utf-8',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=trace_fields,delimiter='\t',extrasaction='ignore');writer.writeheader()
        for row in audited:
            for trace in row.get('_rule_gate_trace',[]):
                if not trace.get('reached') and not trace.get('diagnostic_cdd_matched'):continue
                entry=dict(trace,Specie=row['Specie'],protein_id=row['protein_id'],stage1_family=row['preaudit_family'])
                writer.writerow({k:json.dumps(v,sort_keys=True) if isinstance(v,(dict,list)) else v for k,v in entry.items()})
    candidates=collapse_rows(audited)
    with (dest.parent/'stage2_candidate_gate_trace.tsv').open('w',encoding='utf-8',newline='') as handle:
        fields=['Specie','protein_id','assigned_family','family_support_level','annotation_resolution','grade',
                'domain_gate_passed','stage2_excluded','exclusion_gate','identity_gate_trace']
        writer=csv.DictWriter(handle,fieldnames=fields,delimiter='\t',extrasaction='ignore');writer.writeheader()
        for row in candidates:
            if row.get('stage2_excluded'):
                trace=json.loads(row['identity_gate_trace']);trace['reason_codes']=sorted(set(trace.get('reason_codes',[]))|{'NEGATIVE_EXCLUSION'})
                row=dict(row,identity_gate_trace=json.dumps(trace,sort_keys=True))
            writer.writerow(row)
    screening=[]
    for row in audited:
        item={key:value for key,value in row.items() if not key.startswith('_')};item['assigned_family']=row['stage1_assigned_family'];screening.append(item)
    tables=dict(supporting.get('tables',{}))
    tables.update(Screening_hypothesis=screening,Audited_family_assignment=resolved)
    counts=Counter(row['grade'] for row in resolved)
    exclusions=list({unique_key(row):row for row in audited if row.get('stage2_excluded')}.values())
    exclusion_fields=['Specie','protein_id','control_class','control_accession','exclusion_gate','exclusion_gates','exclusion_reason','excluded_putative_family','excluded_protein_name']
    with (dest.parent/'negative_control_exclusions.tsv').open('w',encoding='utf-8',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=exclusion_fields,delimiter='\t',extrasaction='ignore')
        writer.writeheader();writer.writerows(exclusions)
    from identity_probability import annotate as probability_context
    calibration=probability_context(resolved,family_rules_sha256=rules.get('_sha256'),orthogroup_context_weight=rules.get('_orthogroup_context_weight'))
    metrics=[('Pipeline version',VERSION),('Original stage 1 hypotheses',len(rows)),('Stage 2 recovery hypotheses',sum(r['row_origin']=='stage2_recovery' for r in audited)),
        ('Stage 1 species/protein pairs',len({unique_key(r) for r in rows})),('Final species/protein pairs',len(resolved)),
        ('Negative controls excluded',len(exclusions)),('All input proteins accounted for',len(resolved)+len(exclusions)==len({unique_key(r) for r in rows})),
        ('Collapsed hypothesis records',len(audited)-len(resolved)-len(exclusions)),
        *[(f'Grade {grade}',counts[grade]) for grade in 'ABCD'],('Domain score maximum',4),('Orthology score maximum',3),('Motif score maximum',2),('Localization score maximum',1),
        ('Family/protein probability status',calibration['status'])]
    metrics.append(('Orthogroup context policy',json.dumps(context_configuration(rules.get('_orthogroup_context_weight')),sort_keys=True)))
    metrics.extend(('Benchmark '+key,value) for key,value in supporting.get('benchmark',{}).items() if not isinstance(value,(dict,list)))
    metrics.extend(('Identity calibration '+key,json.dumps(item,sort_keys=True) if isinstance(item,(dict,list)) else item) for key,item in calibration.items())
    tables['Pipeline_summary']=[dict(Metric=key,Value=item) for key,item in metrics]
    from audit_report import prepare_sheet
    cleanup={};book=Workbook();book.remove(book.active)
    for spec in schema(2):
        records,fields,removed=prepare_sheet(spec['name'],tables.get(spec['name'],[]))
        cleanup[spec['name']]=removed
        write_sheet(book,2,spec['name'],records,fields=fields)
        if spec['name']=='Audited_family_assignment':
            with (dest.parent/'audited_family_assignments.tsv').open('w',encoding='utf-8',newline='') as handle:
                writer=csv.DictWriter(handle,fieldnames=fields,delimiter='\t',extrasaction='ignore')
                writer.writeheader();writer.writerows(records)
    book.save(dest)
    calibration_fields=['protein_id','Specie','sequence_sha256','assigned_family','protein_name','grade','family_confidence',
        'protein_name_scope','total_score_0_10','orthology_panel_sha256']
    with (dest.parent/'identity_calibration_features.tsv').open('w',encoding='utf-8',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=calibration_fields,delimiter='\t',extrasaction='ignore')
        writer.writeheader();writer.writerows(resolved)
    manifest=dict(pipeline_version=VERSION,input_sha256=checksum(source),rule_sha256=rules.get('_sha256'),
        policy_sha256=rules.get('_policy_sha256'),output_sha256=checksum(dest),counts=dict(counts),
        identity_library_sha256=checksum(family_identity.LIBRARY),identity_library_path=str(family_identity.LIBRARY.resolve()),
        family_rules_path=rules.get('_source_path'),identity_probability_calibration=calibration,column_cleanup=cleanup,
        orthogroup_context=context_configuration(rules.get('_orthogroup_context_weight')),
        input_hypotheses=len(rows),audited_hypotheses=len(audited),unique_proteins=len(resolved),
        input_unique_proteins=len({unique_key(r) for r in rows}),excluded_negative_controls=len(exclusions),
        protein_pool_preserved=not exclusions,protein_accounting_complete=True)
    from stage2_evidence import profile,CONFIG
    manifest['stage2_profile']=profile(rules)[0]
    manifest['stage2_evidence_rules_sha256']=checksum(CONFIG)
    manifest['phylogenetic_identity_resolution']=phylogenetic_status
    dest.with_suffix('.manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    return manifest


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input',type=Path);parser.add_argument('output',type=Path);parser.add_argument('--sheet',default='candidates')
    parser.add_argument('--rules',type=Path,default=Path(__file__).with_name('family_rules.json'))
    parser.add_argument('--policy',type=Path,default=Path(__file__).with_name('audit_policy.json'))
    parser.add_argument('--identity-rules',type=Path,default=Path(__file__).with_name('family_identity_rules.json'),
                        help='Identity rules paired with --rules; the workflow selects both together')
    parser.add_argument('--supporting-data',type=Path);parser.add_argument('--keep-all-audited',action='store_true')
    parser.add_argument('--phylogenetic-identity-resolution',action='store_true',help='Optional tree concordance review; missing/failed evidence never blocks the family audit')
    parser.add_argument('--phylogenetic-identity-evidence',type=Path,help='Supported-clade TSV; read only when optional resolution is explicitly enabled')
    parser.add_argument('--orthogroup-context-weight',type=float,default=None,help='Optional heuristic context coefficient (0-1); omitted preserves baseline; 0 disables context; never a correctness probability')
    from stage2_evidence import configuration,default_profile
    parser.add_argument('--stage2-profile',choices=list(configuration()['profiles']),default=default_profile(),
                        help='Explicit evidence recruitment/confirmation profile; no threshold fitting at runtime')
    parser.add_argument('--version',action='version',version='adhesome-audit '+VERSION)
    args=parser.parse_args();rules=json.loads(args.rules.read_text());rules['_sha256']=checksum(args.rules)
    rules['_orthogroup_context_weight']=validate_weight(args.orthogroup_context_weight)
    rules['_stage2_profile']=args.stage2_profile
    rules['_phylogenetic_identity_resolution']=args.phylogenetic_identity_resolution
    rules['_phylogenetic_identity_evidence']=args.phylogenetic_identity_evidence
    rules['_source_path']=str(args.rules.resolve())
    from evaluate_family_expansion import validate_pair
    validate_pair(rules,json.loads(args.identity_rules.read_text()))
    rules['_audit_policy']=json.loads(args.policy.read_text());rules['_policy_sha256']=checksum(args.policy)
    support=json.loads(args.supporting_data.read_text()) if args.supporting_data else None
    print(json.dumps(run_xlsx(args.input,args.output,args.sheet,rules,support,args.identity_rules),indent=2))


if __name__=='__main__':main()
