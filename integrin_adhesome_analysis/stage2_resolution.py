"""Explicit annotation resolution and compact gate explanations, version 1.0.

Broad support requires its own complete observed core. Partial specific-family
gates are never promoted merely because orthologues or calibration controls exist.
"""
import json

BROAD_LABELS={'actin':'Actin-family protein','Src_Csk_like':'Src/Csk-related tyrosine kinase',
              'PET_LIM_like':'PET-LIM family protein','collagen_like':'Collagen-like protein',
              'laminin_like':'Laminin-related protein','fibronectin_like':'Fibronectin-like protein',
              'fibronectin_like_candidate':'Fibronectin-like protein'}
ANCESTORS={'actin':{'actin_conventional',*(f'arp{i}_like' for i in (1,2,3,4,5,6,8,10))},
           'Src_Csk_like':{'Src','Csk'},'PET_LIM_like':{'Testin_like','Prickle_like'}}


def common_broad_family(families):
    names=set(families)-{'','unassigned'}
    return next((parent for parent,members in ANCESTORS.items() if len(names)>1 and names<=members),None)


def broad_architecture(row,result):
    """Resolve an explicitly shared full core, retaining narrower limitations."""
    from family_identity import distinct_footprints,ordered_core,diagnostic_cdd
    hits=result.get('combined_hits',[])
    competitors=result.get('identity_competitors',[])
    parent=common_broad_family(competitors)
    if parent=='actin':
        groups=[dict(models=['PF00022'],count=1)];order=[]
    elif parent=='Src_Csk_like':
        groups=[dict(models=['PF00018'],count=1),dict(models=['PF00017'],count=1),dict(models=['PF07714','PF00069'],count=1)];order=[0,1,2]
    else:
        # PET+three LIM domains establishes the broad PET-LIM architecture;
        # five distinct LIM copies may also fit it, without a Testin name.
        permitted={'PET_LIM_like','Testin_like','Prickle_like','paxillin_like','zyxin_like','PDLIM_family',
                   'PINCH_like','LIM_scaffold_screen','FHL_screen','component_reference_candidate','unassigned_domain_candidate'}
        if result.get('family') not in permitted or result.get('identity_rule'):return result
        if competitors and not set(competitors)<=ANCESTORS['PET_LIM_like']:return result
        if not diagnostic_cdd(row,r'^(?:PET)(?:_|$)|TESTIN|PRICKLE'):return result
        parent='PET_LIM_like';groups=[dict(models=['PET'],count=1),dict(models=['PF00412'],count=3)];order=[0,1]
        lim=len(distinct_footprints(hits,{'PF00412'}))
        if lim>5:return result
    core={m for group in groups for m in group['models']}
    if not all(len(distinct_footprints(hits,set(group['models'])))>=group['count'] for group in groups):return result
    if not any(h['model'] in core and 'Pfam' in h['sources'] for h in hits):return result
    length=int(float(row.get('length') or row.get('Length') or 0))
    if any(h['start'] is None or (length and h['end']>length) for h in hits if h['model'] in core):return result
    if not ordered_core(hits,dict(required_groups=groups,ordered_groups=order)):return result
    return dict(result,score=4,fraction=1,family=parent,complete=True,cdd=True,
                note='Complete observed '+parent+' broad-family core; specific lineage/paralogue unresolved',
                identity_rule=None,identity_competitors=[],unresolved_specific_families=competitors,
                broad_gate=dict(required_groups=groups,order=order,passed=True),
                baseline_gate=dict(configured=True,complete=True,reason_codes=[]))


def trace_for(row,result,orthology_score):
    rules=result.get('rule_trace',[])
    reached=[r for r in rules if r.get('reached')]
    diagnostic=[r for r in rules if r.get('diagnostic_cdd_matched')]
    reasons=[];baseline=result.get('baseline_gate',{})
    if result['score']<4:
        reasons.extend(baseline.get('reason_codes',[]))
        if reached:reasons=[code for code in reasons if code!='NO_IDENTITY_RULE_REACHED']
        if not reached and not baseline.get('configured'):reasons.append('NO_IDENTITY_RULE_REACHED')
        # Report evidence limitations for applicable diagnostic rules, rather
        # than the missing cores of every unrelated catalogue entry.
        relevant=diagnostic or [r for r in reached if not r.get('missing_domains')] or reached
        for rule in relevant:
            reasons.extend(rule.get('reason_codes',[]))
        if not reasons:reasons.append('MISSING_REQUIRED_DOMAIN')
    if result.get('identity_competitors'):reasons.append('COMPETING_FAMILY')
    if orthology_score<=0:reasons.append('ORTHOLOGY_UNRESOLVED')
    if result.get('unresolved_specific_families') or result.get('broad_gate'):reasons.append('SPECIFIC_IDENTITY_UNRESOLVED')
    return dict(stage1_hypothesis=row.get('family',''),tested_family=result['family'],
                reason_codes=sorted(set(reasons)),complete_family_gate=result['score']==4,
                reached_rule_ids=[r['rule_id'] for r in reached],matched_rule_ids=[r['rule_id'] for r in rules if r.get('matched')],
                diagnostic_rule_ids=[r['rule_id'] for r in diagnostic],baseline_gate=baseline,
                broad_gate=result.get('broad_gate'),unresolved_specific_families=result.get('unresolved_specific_families',[]))


def annotate_resolution(row,result):
    complete=row['domain_gate_passed'];ambiguous=row.get('family_confidence')=='ambiguous'
    broad=row['assigned_family'] in BROAD_LABELS
    level=('ambiguous' if ambiguous else 'broad_supported' if complete and broad else
           'specific_supported' if complete and row['grade'] in ('A','B') else
           'provisional' if row['grade']=='C' else 'unresolved')
    row.update(family_support_level=level,
               annotation_resolution='broad_family' if level=='broad_supported' else row.get('protein_name_scope','unresolved'),
               specific_identity_status='specific_lineage_or_paralogue_unresolved' if level=='broad_supported' else
               'named_family_supported_exact_gene_not_evaluated' if level=='specific_supported' else 'unresolved',
               identity_gate_trace=json.dumps(trace_for(row,result,row['orthology_score_0_3']),sort_keys=True),
               unresolved_specific_families=';'.join(result.get('unresolved_specific_families',[])))
    if level=='broad_supported':
        row.update(protein_name='Putative '+BROAD_LABELS[row['assigned_family']],protein_name_scope='broad_family',
                   protein_identity_confidence='supported_at_broad_family_scope')
    return row
