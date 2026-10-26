"""Compact version 1.0 audit views; complete raw evidence remains on disk."""

# Resolve the organized bundle independently of the caller working directory.
import sys as _bundle_sys
from pathlib import Path as _BundlePath
_BUNDLE_ROOT=next(p for p in _BundlePath(__file__).resolve().parents if (p/'pipeline_paths.py').is_file())
if str(_BUNDLE_ROOT) not in _bundle_sys.path:_bundle_sys.path.insert(0,str(_BUNDLE_ROOT))
from pipeline_paths import activate as _activate_bundle_paths
_activate_bundle_paths()
import re

FINAL_COLUMNS=['Specie','protein_id','putative_audited_family','protein_name','assignment_confidence',
    'family_support_level','annotation_resolution','specific_identity_status','unresolved_specific_families','identity_gate_trace',
    'assignment_probability','assignment_probability_lower_95','assignment_probability_upper_95',
    'probability_calibration_n','probability_method','probability_calibration_scope',
    'protein_name_scope','grade','family_decision','length','Pfam_architecture','NCBI_CDD_architecture',
    'architecture_evidence','naming_evidence','identity_rule','identity_competitors',
    'domain_score_0_4','orthology_raw_score_0_3','orthology_score_0_3','motif_score_0_2',
    'orthogroup_context_weight','orthogroup_context_policy',
    'localization_score_0_1','total_score_0_10','family_specific_CDD_support',
    'Orthology_Orthogroup','Orthology_HOG','Orthology_HOG_status',
    'orthology_reference_species_count','Orthology_direct_reference_count','Orthology_reference_support',
    'Orthology_benchmark_status','Orthology_benchmark_scope','Orthology_benchmark_acceptance','Orthology_benchmark_weight','orthology_benchmark_weight_status',
    'MotifPrediction_distinct_site_count','MotifPrediction_context_supported_count',
    'MotifPrediction_sequence_only_count','MotifPrediction_context_conflict_count',
    'DeepLoc_2.1','DeepTMHMM','SignalP.v6.0','NCBI_CDD_selected_hit_count','NCBI_CDD_best_evalue',
    'stage1_family_hypothesis','alternative_family_hypotheses','equally_ranked_hypotheses',
    'collapsed_row_count','recovered_hypothesis_count','restriction_scope']

UNUSED={'NCBI_CDD_hits','topology_check','phylogeny','probability_status',
    'reference_subfamily','subfamily_support_method','reciprocal_reference_subfamily',
    'subfamily_profile_calibrated','top_profile_subfamily','profile_margin_pass','profile_domain_gate_verified'}
# Only known aliases are deduplicated. Equal zeros in different scores are meaningful.
ALIASES={'library_id':'motif_id','reported_candidate_tier':'candidate_tier',
    'order_check':'domain_order_check','orthology_benchmark_weight_used':'Orthology_benchmark_weight',
    'orthology_reference_species_count':'Orthology_HOG_reference_species_count'}


def blank(value):return value is None or isinstance(value,str) and not value.strip()


def prepare_sheet(name,rows):
    records=[dict(row) for row in rows]
    if name=='Audited_family_assignment':
        if not records:
            # An all-negative batch has zero assignments, not a placeholder protein.
            return [],['Specie','protein_id','putative_audited_family','protein_name','assignment_confidence','grade'],{}
        for row in records:
            row.update(putative_audited_family=row.get('assigned_family','unassigned'),
                assignment_confidence=row.get('family_confidence','unresolved'),
                Pfam_architecture=row.get('architecture') or row.get('Pfam_architecture',''),
                stage1_family_hypothesis=row.get('preaudit_family') or row.get('family',''))
        possible=FINAL_COLUMNS
    else:possible=list(dict.fromkeys(key for row in records for key in row))
    fields=[];removed={}
    for key in possible:
        if re.sub('[^a-z]','',key.lower()).startswith('mulocdeep') or key in UNUSED:
            removed[key]='unused_or_removed_by_request';continue
        if all(blank(row.get(key)) for row in records):
            removed[key]='entirely_blank';continue
        canonical=ALIASES.get(key)
        if canonical in possible and all(row.get(key)==row.get(canonical) for row in records):
            removed[key]='duplicate_alias_of_'+canonical;continue
        fields.append(key)
    if name=='Audited_family_assignment':
        for key in {key for row in rows for key in row}-set(possible):removed[key]='detail_in_screening_or_supporting_evidence'
    if not records or not fields:return [dict(evidence_status='No records supplied for this evidence type')],['evidence_status'],removed
    return records,fields,removed
