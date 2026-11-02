"""Saved-output validation diagnostics; no inference, fitting or call changes."""

# Resolve source resources independently of the caller working directory.
import sys as _layout_sys
from pathlib import Path as _LayoutPath
_LAYOUT_BUNDLE_ROOT = next(p for p in _LayoutPath(__file__).resolve().parents
                          if (p/'config/layout_paths.json').is_file() and (p/'pipeline_paths.py').is_file())
if str(_LAYOUT_BUNDLE_ROOT) not in _layout_sys.path:
    _layout_sys.path.insert(0, str(_LAYOUT_BUNDLE_ROOT))
from pipeline_paths import activate as _layout_activate, resource_path as _layout_resource
_layout_activate()


from collections import Counter,defaultdict
import json


def yes(value):return str(value).strip().lower() in ('true','1','yes')


def family_labels(value):return {v.strip() for v in str(value or '').split(';') if v.strip()}


def policy_expected_families(row,policy):
    labels=family_labels(row.get('expected_families'))
    if policy=='broad_family_supported':return labels|family_labels(row.get('expected_broad_families'))
    if policy=='specific_family_supported' and 'expected_specific_families' in row:
        return family_labels(row['expected_specific_families'])
    return labels


def policy_prediction(row,policy):
    if policy=='audit_before_negative_exclusion':return row.get('pre_exclusion_family','unassigned')
    if policy=='audit_without_exact_negative_exclusion':return row.get('without_exact_exclusion_family',row.get('audited_family','unassigned'))
    return row.get('audited_family','unassigned')


def macro_family_metrics(scored,policy):
    """Equal weight per admissible gold-family set; alternatives stay together.

    Wrong predictions outside the mapped gold vocabulary form predicted-only
    groups, so macro precision does not silently discard novel wrong labels.
    Undefined denominators are excluded from the corresponding mean and counted.
    """
    evaluated=[r for r in scored if r.get(policy) is not None]
    positives=[r for r in evaluated if r['control_class']=='positive' and policy_expected_families(r,policy)]
    groups={tuple(sorted(policy_expected_families(r,policy))) for r in positives}
    vocabulary={family for group in groups for family in group}
    called=[r for r in evaluated if yes(r.get(policy)) and (r['control_class']!='positive' or policy_expected_families(r,policy))]
    extra={policy_prediction(r,policy) for r in called if policy_prediction(r,policy) not in vocabulary}
    rows=[]
    for labels,predicted_only in [(g,False) for g in sorted(groups)]+[((f,),True) for f in sorted(extra)]:
        members=set(labels)
        actual=[r for r in positives if tuple(sorted(policy_expected_families(r,policy)))==labels] if not predicted_only else []
        tp=sum(yes(r.get(policy)) and policy_prediction(r,policy) in members for r in actual)
        fp=sum(policy_prediction(r,policy) in members and
               (r['control_class']!='positive' or policy_prediction(r,policy) not in policy_expected_families(r,policy)) for r in called)
        fn=len(actual)-tp
        rows.append(dict(policy=policy,expected_admissible_families=';'.join(labels),group_kind='predicted_only' if predicted_only else 'mapped_gold',
            positive_controls=len(actual),TP=tp,FP=fp,FN=fn,precision=tp/(tp+fp) if tp+fp else None,
            recall=tp/(tp+fn) if tp+fn else None,F1=2*tp/(2*tp+fp+fn) if 2*tp+fp+fn else None))
    mean=lambda field:sum(r[field] for r in rows if r[field] is not None)/sum(r[field] is not None for r in rows) if any(r[field] is not None for r in rows) else None
    return dict(macro_precision=mean('precision'),macro_recall=mean('recall'),macro_F1=mean('F1'),macro_family_groups=len(rows),
        macro_precision_defined_groups=sum(r['precision'] is not None for r in rows),
        macro_recall_defined_groups=sum(r['recall'] is not None for r in rows),
        macro_evaluable_positive_controls=len(positives),macro_family_basis='equal_admissible_gold_family_sets_plus_predicted_only_groups'),rows


REASON_CODES=('NO_STAGE1_HYPOTHESIS','NO_IDENTITY_RULE_REACHED','MISSING_REQUIRED_DOMAIN','DOMAIN_ORDER_CONFLICT',
              'CDD_ALIAS_UNRECOGNIZED','ORTHOLOGY_UNRESOLVED','COMPETING_FAMILY','NEGATIVE_EXCLUSION',
              'DIAGNOSTIC_CDD_UNCONFIRMED','DOMAIN_COORDINATES_UNCONFIRMED','LENGTH_OUTSIDE_RULE_RANGE',
              'EXCESS_DOMAIN_COPIES','PFAM_ANCHOR_MISSING','SPECIFIC_IDENTITY_UNRESOLVED')


def gate_diagnostic(row):
    """Report established failures only; unavailable traces are not diagnoses."""
    raw=row.get('identity_gate_trace') or row.get('review_identity_gate_trace') or ''
    trace=None;status='unavailable'
    if raw:
        if isinstance(raw,(dict,list)):trace=raw;status='available'
        else:
            try:trace=json.loads(raw);status='available'
            except (TypeError,ValueError):status='invalid_json'
    reasons=set()
    def visit(node,reason_field=False):
        if isinstance(node,dict):
            for key,value in node.items():visit(value,key in ('reason_codes','limitation_codes','reasons','reason_code','reason'))
        elif isinstance(node,(list,tuple)):
            for value in node:visit(value,reason_field)
        elif reason_field and isinstance(node,str):
            for value in node.replace(',',';').split(';'):
                if value.strip() in REASON_CODES:reasons.add(value.strip())
    if trace is not None:visit(trace)
    if not yes(row.get('stage1_recruitment')):reasons.add('NO_STAGE1_HYPOTHESIS')
    if yes(row.get('negative_exclusion')):reasons.add('NEGATIVE_EXCLUSION')
    if row.get('assignment_confidence')=='ambiguous' or row.get('review_assignment_confidence')=='ambiguous' or row.get('family_support_level')=='ambiguous':reasons.add('COMPETING_FAMILY')
    ordered=[code for code in REASON_CODES if code in reasons]
    return dict(row,reason_codes=';'.join(ordered),primary_reason_code=ordered[0] if ordered else '',
        gate_trace_status=status,diagnostic_scope='saved_gate_trace_and_observed_disposition; no inferred missing CDD aliases')


def sequence_groups(scored):
    groups=defaultdict(list)
    for row in scored:groups[row.get('sequence_sha256') or row['protein_id']].append(row)
    return groups


def populations(scored):
    yield 'all_controls',scored
    for name,allowed in [('positive_vs_hard_negative',{'positive','hard_negative'}),
                         ('positive_vs_background_negative',{'positive','background_negative'})]:
        yield name,[r for r in scored if r['control_class'] in allowed]
    yield 'sequence_verified_controls',[r for r in scored if r['label_sequence_verification']=='trusted_fasta_match']
    for label,key in [('identity_calibration','identity_calibration_training_overlap'),
                      ('forward_reference','exact_forward_reference_overlap'),('reverse_reference','exact_reverse_reference_overlap')]:
        for overlap in (True,False):
            yield ('' if overlap else 'no_')+label+'_overlap',[r for r in scored if yes(r.get(key,False))==overlap]
    keys=('exact_forward_reference_overlap','exact_reverse_reference_overlap','identity_calibration_training_overlap')
    yield 'no_reference_or_calibration_overlap',[r for r in scored if not any(yes(r.get(key,False)) for key in keys)]
    for train in (True,False):
        for reference in (True,False):
            name=f'calibration_{"overlap" if train else "absent"}__forward_{"overlap" if reference else "absent"}'
            yield name,[r for r in scored if yes(r.get(keys[2],False))==train and yes(r.get(keys[0],False))==reference]
    # Exact duplicates count once only when their labels AND saved decisions
    # agree. Conflicts are exported separately, never silently selected.
    policies=('stage1_recruitment','dual_gate_supported','dual_gate_strong','audit_including_provisional',
              'audit_before_negative_exclusion','audit_without_exact_negative_exclusion','stage2_candidate_retention',
              'broad_family_supported','specific_family_supported','provisional_annotation','ambiguous_annotation')
    unique=[]
    for rows in sequence_groups(scored).values():
        signatures={(r['control_class'],r.get('expected_families',''),r.get('expected_broad_families',''),
                     r.get('expected_specific_families',''),r.get('audited_family',''),
                     *((None if r.get(p) is None else yes(r[p])) for p in policies)) for r in rows}
        if len(signatures)==1:unique.append(rows[0])
    yield 'unique_sequence_consistent_controls',unique


def review_reason(row):
    if not yes(row['stage1_recruitment']):
        return 'not_recruited' if row.get('preparation_decision')=='retained' else 'lost_during_preparation'
    if yes(row['negative_exclusion']):return 'excluded_positive_requires_control_label_review'
    if row.get('review_assignment_confidence')=='ambiguous':return 'competing_family_hypotheses'
    diagnostic=gate_diagnostic(row)
    if diagnostic['primary_reason_code']:return diagnostic['primary_reason_code'].lower()
    note=str(row.get('review_architecture_evidence','')).lower()
    for text,reason in [('no family-specific architecture','family_gate_not_configured'),
                        ('generic screening bucket','screening_bucket_without_diagnostic_refinement'),
                        ('family-specific cdd absent','missing_required_diagnostic_cdd'),
                        ('order unconfirmed','domain_order_not_established'),
                        ('coordinates unavailable','domain_coordinates_missing'),
                        ('length range','outside_family_length_range')]:
        if text in note:return reason
    return 'partial_or_incompatible_architecture'


def reports(predictions,metrics):
    scored=[r for r in predictions if r['evaluation_status']=='scored']
    positives=[r for r in scored if r['control_class']=='positive']
    review=[];misclassified=[];mapping=[];tiers=Counter();coverage=Counter();family_groups=defaultdict(list)
    for row in scored:
        level=row.get('family_support_level','')
        tier=('not_recruited' if not yes(row['stage1_recruitment']) else
              'negative_excluded' if yes(row['negative_exclusion']) else
              level if yes(row['dual_gate_supported']) and level in ('broad_supported','specific_supported') else
              'supported' if yes(row['dual_gate_supported']) else
              'ambiguous_family' if row.get('assignment_confidence')=='ambiguous' or level=='ambiguous' else
              'provisional_only' if yes(row['audit_including_provisional']) else 'recruited_unresolved')
        tiers[row['control_class'],tier]+=1
    for row in positives:
        expected=set(row['expected_families'].split(';'))-{''}
        wrong=yes(row['dual_gate_supported']) and bool(expected) and row['audited_family'] not in expected
        if not yes(row['dual_gate_supported']) or wrong:
            review.append(dict(row,review_reason='incorrect_accepted_family' if wrong else review_reason(row)))
        if wrong:misclassified.append(dict(row,review_reason='incorrect_accepted_family'))
        mapping.append({k:row.get(k,'') for k in ('accession','protein_id','trusted_protein_name','trusted_genes','expected_families',
                                                 'expected_broad_families','expected_specific_families','family_label_source','family_rule_coverage',
                                                 'family_label_review_required','family_label_review_reason','legacy_catalogue_families_ignored')})
        coverage[row.get('family_rule_coverage') or 'needs_mapping_review']+=1
        family_groups[row['expected_families'] or '__unmapped__'].append(row)
    per_family=[]
    for family,rows in sorted(family_groups.items()):
        accepted=[r for r in rows if yes(r['dual_gate_supported'])]
        correct=sum(r['audited_family'] in set(r['expected_families'].split(';')) for r in accepted) if family!='__unmapped__' else None
        per_family.append(dict(expected_admissible_families=family,positive_controls=len(rows),
            recruited=sum(yes(r['stage1_recruitment']) for r in rows),supported=len(accepted),correct_supported=correct,
            wrong_supported=len(accepted)-correct if correct is not None else None,
            not_supported=len(rows)-len(accepted),provisional_only=sum(yes(r['audit_including_provisional']) and not yes(r['dual_gate_supported']) for r in rows)))
    ablation=[dict(row) for row in metrics if row['policy'] in
        ('audit_before_negative_exclusion','audit_without_exact_negative_exclusion','dual_gate_supported') and
        row['population'] in ('all_controls','positive_vs_hard_negative','positive_vs_background_negative')]
    impact=[]
    for gate in sorted({r['exclusion_gate'] for r in scored if yes(r['negative_exclusion'])}):
        for category in ('positive','hard_negative','background_negative'):
            rows=[r for r in scored if r['exclusion_gate']==gate and r['control_class']==category and yes(r['negative_exclusion'])]
            impact.append(dict(exclusion_gate=gate,control_class=category,excluded=len(rows),
                               passed_family_gate_before_exclusion=sum(yes(r['audit_before_negative_exclusion']) for r in rows)))
    conflicts=[]
    for digest,rows in sequence_groups(scored).items():
        if len({(r['control_class'],r.get('expected_families',''),r.get('audited_family',''),
                 yes(r['dual_gate_supported'])) for r in rows})>1:
            conflicts.append(dict(sequence_sha256=digest,accessions=';'.join(r['accession'] for r in rows),
                                  reason='exact_sequence_labels_or_saved_decisions_disagree'))
    mapped=sum(bool(family_labels(r.get('expected_families'))) for r in positives)
    coverage_summary=dict(positive_controls=len(positives),mapped_expected_family_controls=mapped,
        unmapped_expected_family_controls=len(positives)-mapped,mapped_fraction=mapped/len(positives) if positives else None,
        mapped_specific_family_controls=sum(bool(policy_expected_families(r,'specific_family_supported')) for r in positives),
        broad_ancestor_rollup_controls=sum(bool(family_labels(r.get('expected_broad_families'))) for r in positives),
        source_scope='trusted labels/name ontology; broad ancestors are declared ontology roll-ups, not independent gold or sequence-linked holdouts')
    diagnostics=[gate_diagnostic(r) for r in predictions]
    reason_counts=Counter((r.get('control_class',''),code) for r in diagnostics for code in r['reason_codes'].split(';') if code)
    macro_rows=[];per_group=[]
    for policy in sorted({r['policy'] for r in metrics}-{'stage1_recruitment','stage2_candidate_retention','ambiguous_annotation'}):
        macro,groups=macro_family_metrics(scored,policy);macro_rows.append(dict(policy=policy,**macro));per_group.extend(groups)
    return {'exclusion_ablation.tsv':ablation,'exclusion_impact.tsv':impact,
            'confidence_tiers.tsv':[dict(control_class=c,tier=t,count=n) for (c,t),n in sorted(tiers.items())],
            'family_label_coverage.tsv':[dict(family_rule_coverage=k,positive_controls=n) for k,n in sorted(coverage.items())],
            'family_mapping_review.tsv':mapping,'missed_positive_review.tsv':review,
            'family_misclassifications.tsv':misclassified,'per_expected_family.tsv':per_family,
            'sequence_group_conflicts.tsv':conflicts,'family_label_coverage_summary.tsv':[coverage_summary],
            'stage2_gate_diagnostics.tsv':diagnostics,
            'stage2_reason_code_counts.tsv':[dict(control_class=c,reason_code=code,count=n) for (c,code),n in sorted(reason_counts.items())],
            'family_macro_metrics.tsv':macro_rows,'family_per_group_metrics.tsv':per_group}
