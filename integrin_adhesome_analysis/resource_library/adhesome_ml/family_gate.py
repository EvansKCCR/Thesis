"""FASTA-derived family learning. Validation labels never enter fitting or thresholds."""
import hashlib
import html
import json
import re
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

SCHEMA = 'sequence_family_gate_v5_6_family_similarity_gates'
LEGACY_SCHEMAS = {'sequence_family_gate_v5_4_verified_cdd_bank_coverage'}
from cdd_features import build_calibration_bank, similarity_features
from identity_rule_features import features as identity_rule_features, load_rules as load_identity_rules
from posthoc_rule_features import extract as posthoc_sequence_features, load_rules as load_posthoc_rules


def sequence_digest(sequence):
    return hashlib.sha256(sequence.encode()).hexdigest()


def load_resources(c):
    rules = json.loads((c._BUNDLE / 'config/annotation/family_rules.json').read_text(encoding='utf-8'))
    aliases = {}
    for field in ('family_aliases', 'canonical_family_aliases'):
        for key, value in rules.get(field, {}).items():
            if isinstance(value, str):
                aliases[key.casefold()] = value
    parents = {key: value['parent_family'] for key, value in rules['families'].items()
               if value.get('parent_family') and value['parent_family'] != key}
    library = c._BUNDLE / 'config/motif_annotation/motif_master.json'
    motifs = [(x['rule_id'], x['regex_python']) for x in json.loads(library.read_text(encoding='utf-8'))['rules']]
    for _, pattern in motifs:
        re.compile(pattern)
    return aliases, parents, motifs


def canonical(label, aliases):
    label = label.strip()
    visited = set()
    while label.casefold() in aliases and label.casefold() not in visited:
        visited.add(label.casefold())
        label = aliases[label.casefold()]
    return label


def ancestors(label, parents):
    found = []
    while label in parents and parents[label] not in found and parents[label] != label:
        label = parents[label]
        found.append(label)
    return found


def trusted_family(row, aliases, parents):
    # Expected families are reporting reference labels, never Stage 1 hypotheses or predictions.
    labels = {canonical(x, aliases) for x in re.split(r'[;|]', str(row.get('expected_families') or row.get('trusted_family') or '')) if x.strip()}
    if len(labels) == 1:
        return next(iter(labels))
    if labels:
        common = set.intersection(*(set([x] + ancestors(x, parents)) for x in labels))
        for candidate in [min(labels)] + ancestors(min(labels), parents):
            if candidate in common:
                return candidate
    return ''


def verified_sequences(directory, c):
    """Verify saved Pfam evidence against the exact FASTA scanned, not a labelled match lookup."""
    if not directory or not directory.is_dir():
        raise ValueError('Supply the completed Stage 1 detection directory with prepared FASTAs and domains/scan.json')
    sequences, hashes = {}, {}
    for species in sorted(p for p in directory.iterdir() if p.is_dir()):
        scan = species / 'domains/scan.json'
        prepared = species / 'prepared' / (species.name + '.faa')
        if not scan.is_file():
            continue
        meta = json.loads(scan.read_text(encoding='utf-8'))
        if not prepared.is_file() or c.sha(prepared) != meta.get('proteome_sha256'):
            raise ValueError(f'Prepared FASTA does not match the saved Pfam scan: {prepared}')
        domains = species / 'domains/domains.csv'
        if c.sha(domains) != meta.get('artifact_checksums', {}).get('domains.csv'):
            raise ValueError(f'Pfam hit table differs from scan receipt: {domains}')
        hashes[str(scan)] = c.sha(scan)
        for record in c.read_fasta(prepared):
            ident = record['protein_id']
            if ident in sequences:
                raise ValueError(f'Duplicate prepared FASTA ID: {ident}')
            sequences[ident] = record['sequence']
    if not sequences:
        raise ValueError(f'No verified prepared FASTAs found in {directory}')
    return sequences, hashes


def extract(rows, sequences, evidence, motifs, vocab, c, cached=None, cdd=None):
    raw, hashes, status = cached if cached is not None else c.read_stage1_evidence(evidence)
    raw = {r['protein_id']: r for r in raw}
    posthoc_script=c._BUNDLE.parent/'schisto_pfam_cdd_reclassify/schisto_pfam_cdd_reclassify.py'
    posthoc_rules_sha,_,_,_=load_posthoc_rules(str(posthoc_script))
    compiled = [(name, re.compile(pattern)) for name, pattern in motifs]
    output = []
    for i, original in enumerate(rows):
        ident = c.row_id(original)
        if ident not in sequences:
            raise ValueError(f'No verified sequence for {ident}; candidate cohort cannot be silently reduced')
        seq = sequences[ident]
        # Only raw scan evidence enters features. Audit scores, Stage 1 calls and label fields are excluded.
        row = {'protein_id': ident, 'ml_seq_len_aa': str(len(seq))}
        for field in ('ml_domain_hits_json', 'ml_homology_hits_json', 'ml_reciprocal_pair_count'):
            if field in raw.get(ident, {}):
                row[field] = raw[ident][field]
        sites = []
        for name, pattern in compiled:
            for match in re.finditer('(?=(' + pattern.pattern + '))', seq):
                if match.end(1) > match.start(1):
                    sites.append({'motif': name, 'start': match.start(1) + 1, 'end': match.end(1),
                                  'state': 'sequence_match', 'context': 'not_assessed'})
        row['ml_motif_hits_json'] = json.dumps(sites)
        row = c.apply_scan_status([row], status)[0]
        vector = c.feature_dict(row, vocab)
        cdd_row=cdd.get(ident,{}) if cdd else {}
        if cdd_row:
            row['_family_similarity_hsps']=cdd_row.get('_family_similarity_hsps',[])
            vector.update({key:value for key,value in cdd_row.items() if not key.startswith('_')})
        vector.update(identity_rule_features(row,cdd_row,
            c._BUNDLE/'config/annotation/family_identity_rules.json',c._domain_records))
        posthoc_features,source_sha=posthoc_sequence_features(row,cdd_row,posthoc_script)
        if source_sha!=posthoc_rules_sha: raise ValueError('Post-hoc sequence-rule source changed during feature extraction')
        vector.update(posthoc_features)
        for site in sites:
            region = 'n' if site['start'] / len(seq) < .33 else 'c' if site['start'] / len(seq) > .67 else 'middle'
            key = 'motif:position:' + site['motif'] + ':' + region
            vector[key] = vector.get(key, 0) + 1
        output.append((row, vector))
        if (i + 1) % 250 == 0:
            print(f'FASTA features: {i + 1}/{len(rows)}', flush=True)
    hashes['posthoc_sequence_rules_sha256']=posthoc_rules_sha
    return output, hashes


def read_cohort(branch, c):
    rows = c.read_tsv(branch / 'reports/predictions.tsv')
    if any(not c.normalized_control_class(r) for r in rows):
        raise ValueError('Unrecognized control label; evaluation denominators must include every labelled record')
    ids = [c.row_id(r) for r in rows]
    if len(ids) != len(set(ids)) or any(not x for x in ids):
        raise ValueError('Control cohort IDs must be unique and populated')
    sequences, hashes = verified_sequences(branch / 'handoff/detection', c)
    manifest = {r['protein_id']: r for r in c.read_tsv(branch / 'inputs/sequence_manifest.tsv')}
    for row in rows:
        ident = c.row_id(row)
        if ident not in manifest or ident not in sequences:
            raise ValueError(f'Control sequence absent from verified input cohort: {ident}')
        if sequence_digest(sequences[ident]) != manifest[ident]['sequence_sha256']:
            raise ValueError(f'Control FASTA differs from labelled input: {ident}')
    return rows, sequences, hashes


def domain_set(row, c):
    return {x['name'].split('.')[0].upper() for x in c._domain_records(row)}


def membership_model_vector(vector):
    """Keep the membership model's established evidence contract unchanged."""
    return {key:value for key,value in vector.items() if not key.startswith('family_similarity:')}


def _family_token(value):
    return re.sub(r'[^A-Za-z0-9_.-]+','_',str(value or '')).strip('_').casefold()


def _family_similarity_config(c, saved=None):
    return (saved or {}).get('family_similarity_settings') or c._SETTINGS['family_gate']['family_similarity']


def _family_similarity_profile(settings, requested=None):
    name=requested or settings.get('default_profile','balanced')
    profiles=settings.get('profiles',{})
    if name not in profiles:
        raise ValueError(f"Unknown family similarity profile {name!r}; choose from {', '.join(sorted(profiles))}")
    selected=profiles[name]
    search_max=float(settings.get('search_max_evalue',0))
    base_max=float(settings.get('legacy_model_feature_max_evalue',0))
    if search_max<=0 or base_max<=0 or search_max<base_max:
        raise ValueError('Family similarity search_max_evalue must be positive and at least the legacy model-feature E-value')
    for level in ('specific','broad'):
        values=selected.get(level,{})
        if float(values.get('max_evalue',0))<=0 or float(values.get('max_evalue',0))>search_max:
            raise ValueError(f"Family similarity {name}/{level} max_evalue must be positive and no greater than search_max_evalue")
        if not 0<=float(values.get('min_query_coverage',0))<=1 or not 0<=float(values.get('min_subject_coverage',0))<=1:
            raise ValueError(f"Family similarity {name}/{level} coverage limits must be between 0 and 1")
    return name,selected


def _rule_similarity_thresholds(family, level, selected, identity_rules, aliases):
    defaults=dict(selected[level]); thresholds=[(defaults,'selected_profile')]
    target=_family_token(canonical(str(family),aliases))
    for rule in identity_rules:
        if _family_token(canonical(str(rule.get('family','')),aliases))!=target:
            continue
        override=(rule.get('sequence_similarity_policy',{}) or {}).get(level)
        if override:
            if not isinstance(override,dict):
                raise ValueError(f"Family identity rule {rule.get('id','unknown')} {level} similarity policy must be an object")
            unknown=set(override)-{'max_evalue','min_query_coverage','min_subject_coverage'}
            if unknown:
                raise ValueError(f"Unknown family similarity fields in rule {rule.get('id','unknown')}: {', '.join(sorted(unknown))}")
            thresholds.append(({**defaults,**override},'family_identity_rule:'+str(rule.get('id','unknown'))))
    return thresholds


def _family_similarity_match(row, family, level, settings, selected, identity_rules, aliases):
    family_key=_family_token(canonical(str(family),aliases))
    hits=row.get('_family_similarity_hsps',[]) or []
    matched=[]
    for threshold,threshold_source in _rule_similarity_thresholds(family,level,selected,identity_rules,aliases):
        maximum=float(threshold['max_evalue'])
        min_query=float(threshold.get('min_query_coverage',0))
        min_subject=float(threshold.get('min_subject_coverage',0))
        if maximum<=0 or maximum>float(settings['search_max_evalue']):
            raise ValueError(f'Family-specific E-value threshold for {family} must be positive and not exceed the CDD search limit')
        if not 0<=min_query<=1 or not 0<=min_subject<=1:
            raise ValueError(f'Family-specific coverage thresholds for {family} must be between 0 and 1')
        for hit in hits:
            if not any(_family_token(canonical(str(label),aliases))==family_key for label in hit.get('families',[])):
                continue
            if (float(hit.get('evalue',1.0))<=maximum and
                float(hit.get('query_coverage',0.0))>=min_query and
                float(hit.get('subject_coverage',0.0))>=min_subject):
                matched.append({**hit,'_matched_threshold':threshold,'_matched_threshold_source':threshold_source})
    return min(matched,key=lambda hit:(float(hit.get('evalue',1.0)),-float(hit.get('bitscore',0)))) if matched else None


def decide(membership, ranking, row, profiles, counts, parents, policy, c,
           vector=None, identity_rules=None, aliases=None, rule_families=None,
           similarity_settings=None, similarity_selection=None):
    vector=vector or {}; identity_rules=identity_rules or []; aliases=aliases or {}; rule_families=set(rule_families or [])
    observed = domain_set(row, c)
    if not any(c.evidence_groups(row).values()):
        return 'unresolved', '', 'no_informative_features'
    if membership < policy['membership_cutoff']:
        return 'non_adhesome', '', 'membership_below_cutoff'
    top, score = ranking[0]
    margin = score - (ranking[1][1] if len(ranking) > 1 else 0)
    shared = observed & set(profiles[top])
    specific_similarity=(_family_similarity_match(row,top,'specific',similarity_settings,similarity_selection,
        identity_rules,aliases) if similarity_settings and similarity_selection else None)
    rule_specific=top in rule_families
    if (shared or specific_similarity or rule_specific) and counts[top] >= policy['minimum_family_examples'] and score >= policy['specific_vote_cutoff'] and margin >= policy['specific_margin']:
        support_reason=('family_vote_margin_and_identity_rule_support' if rule_specific else
                        'family_vote_margin_and_evalue_similarity_support' if specific_similarity else
                        'family_vote_margin_and_shared_domain_support')
        if any(top in ancestors(x, parents) for x in counts if x != top):
            return 'broad', top, support_reason
        return 'specific', top, support_reason
    aggregated = defaultdict(float)
    for label, vote in ranking:
        for parent in ancestors(label, parents):
            aggregated[parent] += vote
    for parent, vote in sorted(aggregated.items(), key=lambda x: (-x[1], x[0])):
        related = [name for name, _ in ranking if parent in ancestors(name, parents)]
        domain_support=any(observed & set(profiles[name]) for name in related)
        similarity_support=bool(similarity_settings and similarity_selection and any(
            _family_similarity_match(row,name,'broad',similarity_settings,similarity_selection,
                identity_rules,aliases) for name in related))
        rule_support=any(name in rule_families for name in related)
        if vote >= policy['broad_vote_cutoff'] and (domain_support or similarity_support or rule_support):
            reason=('shared_parent_vote_and_identity_rule_support' if rule_support else
                    'shared_parent_vote_and_evalue_similarity_support' if similarity_support else
                    'shared_parent_vote_and_domain_support')
            return 'broad', parent, reason
    return ('domain_only' if observed else 'unresolved'), '', 'insufficient_family_resolution'


def infer(saved, data, c, membership_cutoff=None, family_similarity_profile=None):
    rule_path=c._BUNDLE/'config/annotation/family_identity_rules.json'
    expected_rule_sha=saved.get('provenance',{}).get('family_identity_rules')
    if expected_rule_sha and c.sha(rule_path)!=expected_rule_sha:
        raise ValueError('Family identity-rule source differs from the one frozen into this classifier')
    identity_rules,_,_=load_identity_rules(str(rule_path))
    vectors = [x[1] for x in data]
    binary = saved['binary']
    votes = binary.predict_proba([membership_model_vector(vector) for vector in vectors])[:, list(binary.classes_).index('adhesome_associated')]
    family_votes = saved['family'].predict_proba(vectors)
    families = list(saved['family'].classes_)
    policy = dict(saved['policy'])
    similarity_settings=_family_similarity_config(c,saved)
    similarity_profile_name,similarity_selection=_family_similarity_profile(similarity_settings,family_similarity_profile)
    aliases=saved.get('aliases',{})
    if membership_cutoff is not None:
        policy['membership_cutoff'] = membership_cutoff
    output = []
    for i, (row, vector) in enumerate(data):
        ranking = sorted(zip(families, map(float, family_votes[i])), key=lambda x: (-x[1], x[0]))
        proposals=family_identity_rule_proposals(saved,vector,identity_rules)
        proposal_families={x['family'] for x in proposals}
        category, assignment, reason = decide(float(votes[i]), ranking, row, saved['profiles'], saved['counts'], saved['parents'], policy, c,
            vector=vector,identity_rules=identity_rules,aliases=aliases,rule_families=proposal_families,
            similarity_settings=similarity_settings,similarity_selection=similarity_selection)
        candidate_family=assignment if category=='specific' else ranking[0][0]
        candidate_level='specific'
        matched_similarity=_family_similarity_match(row,candidate_family,candidate_level,similarity_settings,
            similarity_selection,identity_rules,aliases)
        if not matched_similarity:
            candidate_level='broad'
            matched_similarity=_family_similarity_match(row,candidate_family,candidate_level,similarity_settings,
                similarity_selection,identity_rules,aliases)
        if category=='broad' and assignment:
            related_candidates=[name for name,_ in ranking if assignment in ancestors(name,saved['parents'])]
            family_hits=[]
            for name in related_candidates:
                hit=_family_similarity_match(row,name,'broad',similarity_settings,similarity_selection,
                    identity_rules,aliases)
                if hit:family_hits.append(hit)
            if family_hits:
                matched_similarity=min(family_hits,key=lambda hit:(float(hit['evalue']),-float(hit['bitscore'])))
                candidate_level='broad'
        groups=c.evidence_groups(row)
        membership_call=('unresolved' if not any(groups.values()) else
            'adhesome_associated' if float(votes[i])>=policy['membership_cutoff'] else 'non_adhesome')
        prediction={'protein_id': row['protein_id'], 'classifier_call': category,
            'predicted_family': assignment, 'best_family_candidate': ranking[0][0],
            'family_vote_fraction': ranking[0][1], 'family_vote_margin': ranking[0][1] - (ranking[1][1] if len(ranking)>1 else 0),
            'adhesome_vote_fraction': float(votes[i]), 'family_ranking_json': json.dumps(ranking),
            'assignment_reason': reason, 'exact_sequence_match_used': 'no',
            'score_interpretation': 'uncalibrated model vote; not a biological probability',
            'membership_call':membership_call,
            'random_forest_call':category,'random_forest_family':assignment,
            'family_similarity_profile':similarity_profile_name,
            'family_similarity_search_max_evalue':float(similarity_settings['search_max_evalue']),
            'family_similarity_assignment_policy_json':json.dumps(similarity_selection,sort_keys=True),
            'family_similarity_support_level':candidate_level if matched_similarity else '',
            'family_similarity_support_evalue':matched_similarity.get('evalue','') if matched_similarity else '',
            'family_similarity_support_query_coverage':matched_similarity.get('query_coverage','') if matched_similarity else '',
            'family_similarity_support_subject_coverage':matched_similarity.get('subject_coverage','') if matched_similarity else '',
            'family_similarity_support_reference_id':matched_similarity.get('reference_id','') if matched_similarity else '',
            'family_similarity_support_annotations':';'.join(matched_similarity.get('annotations',[])) if matched_similarity else '',
            'family_similarity_support_accessions':';'.join(matched_similarity.get('accessions',[])) if matched_similarity else '',
            'family_similarity_support_identity':matched_similarity.get('identity','') if matched_similarity else '',
            'family_similarity_support_bitscore':matched_similarity.get('bitscore','') if matched_similarity else '',
            'family_similarity_support_threshold_source':matched_similarity.get('_matched_threshold_source','') if matched_similarity else '',
            'family_similarity_support_threshold_json':json.dumps(matched_similarity.get('_matched_threshold',{}),sort_keys=True) if matched_similarity else '',
            **c.prediction_evidence_summary(row)}
        families_proposed=sorted({x['family'] for x in proposals})
        prediction.update(rule_family_candidate=families_proposed[0] if len(families_proposed)==1 else '',
            rule_candidate_families=';'.join(families_proposed),
            rule_candidate_status=('unique_rule_candidate_provisional' if len(families_proposed)==1 else
                                   'ambiguous_multiple_family_rules' if families_proposed else 'no_strong_rule_match'),
            rule_candidate_rule_ids=';'.join(sorted({x['rule_id'] for x in proposals})),
            rule_candidate_evidence_json=json.dumps(proposals,sort_keys=True),
            rule_candidate_family_untrained=(len(families_proposed)==1 and families_proposed[0] not in families),
            rule_conflicts_with_random_forest=(len(families_proposed)==1 and bool(assignment) and families_proposed[0]!=assignment),
            prediction_method='random_forest_with_separate_identity_rule_review',
            prediction_policy_version='5.5_identity_rule_fallback_v1')
        output.append(prediction)
    return output


def family_identity_rule_proposals(saved, vector, rules):
    """Return strict sequence-rule candidates separately from RF assignments."""
    proposals=[]
    for rule in rules:
        groups=rule.get('required_groups',[]) or []
        pattern=str(rule.get('specific_cdd_pattern','')).strip()
        if not groups or not pattern or pattern=='(?!)':
            continue
        rid=str(rule.get('id','')).strip(); family=canonical(str(rule.get('family','')).strip(),saved.get('aliases',{}))
        if not rid or not family:
            continue
        prefix=f'identity_rule:{rid}:'
        if vector.get(prefix+'required_all_satisfied',0.0)<0.5:
            continue
        if vector.get(prefix+'specific_cdd_pattern_present',0.0)<0.5:
            continue
        if rule.get('ordered_groups') and vector.get(prefix+'ordered_groups_satisfied',0.0)<0.5:
            continue
        if rule.get('maximum_model_counts') and vector.get(prefix+'maximum_constraints_satisfied',0.0)<0.5:
            continue
        if rule.get('length_range') and vector.get(prefix+'length_in_range',0.0)<0.5:
            continue
        if rule.get('required_cdd_repeat_indices') and vector.get(prefix+'required_cdd_repeats_satisfied',0.0)<0.5:
            continue
        proposals.append({'rule_id':rid,'family':family,'protein_name':rule.get('protein_name',''),
            'required_groups_satisfied':int(vector.get(prefix+'required_satisfied_count',0.0)),
            'required_group_count':len(groups),'specific_cdd_pattern':pattern,
            'specific_cdd_pattern_hit_count':int(vector.get(prefix+'specific_cdd_pattern_hit_count',0.0)),
            'ordered_groups_required':bool(rule.get('ordered_groups')),
            'ordered_groups_satisfied':bool(rule.get('ordered_groups') and vector.get(prefix+'ordered_groups_satisfied',0.0)>=0.5),
            'evidence_status':'sequence-rule-supported; provisional pending independent family-level resolution'})
    return proposals


def evaluate(rows, predictions, sequences, saved, c, family_similarity_profile=None):
    truth = [c.normalized_control_class(r) for r in rows]
    # Score membership from its own output; family-resolution profile changes must not
    # turn an unchanged membership call into a different membership metric.
    calls = ['adhesome_associated' if p.get('membership_call')=='adhesome_associated' else 'non_adhesome'
             for p in predictions]
    binary = c.class_metrics(truth, calls, ['adhesome_associated','non_adhesome'])
    binary['positive_precision'] = c.precision_score(truth, calls, pos_label='adhesome_associated', zero_division=0)
    binary['positive_recall'] = c.recall_score(truth, calls, pos_label='adhesome_associated', zero_division=0)
    binary['positive_f1'] = c.f1_score(truth, calls, pos_label='adhesome_associated', zero_division=0)
    labelled = [(i, trusted_family(r, saved['aliases'], saved['parents'])) for i,r in enumerate(rows) if truth[i]=='adhesome_associated']
    labelled = [(i,label) for i,label in labelled if label]
    expected = [label for _,label in labelled]
    assigned = [predictions[i]['predicted_family'] or '__abstained__' for i,_ in labelled]
    family = c.class_metrics(expected, assigned, sorted(set(expected) | set(assigned)))
    family['labelled_positives'] = len(labelled)
    family['unmapped_positives'] = truth.count('adhesome_associated') - len(labelled)
    family['unseen_label_records'] = sum(label not in saved['counts'] for label in expected)
    family['correct_accepted'] = sum(a==b for a,b in zip(expected,assigned))
    family['accepted'] = sum(x!='__abstained__' for x in assigned)
    family['accepted_precision'] = family['correct_accepted'] / family['accepted'] if family['accepted'] else 0.0
    family['strict_recall'] = family['correct_accepted'] / len(labelled) if labelled else 0.0
    family['top_candidate_correct'] = sum(predictions[i]['best_family_candidate']==label for i,label in labelled)
    family['top_three_contains_truth'] = sum(label in [x[0] for x in json.loads(predictions[i]['family_ranking_json'])[:3]] for i,label in labelled)
    family['correct_broad_or_specific'] = sum(predictions[i]['predicted_family'] in [label]+ancestors(label,saved['parents']) for i,label in labelled)
    family['per_family'] = {label: {'records': expected.count(label), 'correct': sum(a==label and b==label for a,b in zip(expected,assigned)),
                                  'missed_or_wrong': sum(a==label and b!=label for a,b in zip(expected,assigned))} for label in sorted(set(expected))}
    rule_families=[p.get('rule_family_candidate','') for p in predictions]
    rule_positive=[i for i,p in enumerate(predictions) if p.get('rule_candidate_status')=='unique_rule_candidate_provisional']
    rule_ambiguous=sum(p.get('rule_candidate_status')=='ambiguous_multiple_family_rules' for p in predictions)
    rule_labelled=[(i,trusted_family(row,saved['aliases'],saved['parents'])) for i,row in enumerate(rows)
        if truth[i]=='adhesome_associated' and trusted_family(row,saved['aliases'],saved['parents'])]
    rule_covered=sum(bool(rule_families[i]) for i,_ in rule_labelled)
    rule_exact=sum(bool(rule_families[i]) and rule_families[i]==label for i,label in rule_labelled)
    untrained=[(i,label) for i,label in rule_labelled if label not in saved['family'].classes_]
    rule_fallback={'unique_rule_candidate_records':len(rule_positive),
        'ambiguous_rule_candidate_records':rule_ambiguous,
        'positive_controls_with_unique_rule_candidate':sum(truth[i]=='adhesome_associated' for i in rule_positive),
        'non_adhesome_controls_with_unique_rule_candidate':sum(truth[i]=='non_adhesome' for i in rule_positive),
        'membership_precision':sum(truth[i]=='adhesome_associated' for i in rule_positive)/len(rule_positive) if rule_positive else 0.0,
        'membership_recall':sum(truth[i]=='adhesome_associated' for i in rule_positive)/truth.count('adhesome_associated') if truth.count('adhesome_associated') else 0.0,
        'family_labelled_positive_controls':len(rule_labelled),
        'family_labelled_positive_controls_with_unique_rule_candidate':rule_covered,
        'exact_family_correct':rule_exact,
        'exact_family_precision_among_rule_candidates':rule_exact/rule_covered if rule_covered else 0.0,
        'strict_exact_family_recall':rule_exact/len(rule_labelled) if rule_labelled else 0.0,
        'untrained_family_positive_controls':len(untrained),
        'untrained_family_exact_rule_recoveries':sum(bool(rule_families[i]) and rule_families[i]==label for i,label in untrained),
        'rule_rf_conflicts':sum(bool(p.get('rule_conflicts_with_random_forest')) for p in predictions),
        'interpretation':'Separate sequence-rule-supported provisional proposals, not Random Forest calls; no acceptance threshold was tuned on validation.'}
    overlap = [sequence_digest(sequences[c.row_id(row)]) in saved['calibration_sequence_hashes'] for row in rows]
    strata = {}
    for state in (True, False):
        idx = [i for i,x in enumerate(overlap) if x==state]
        strata['sequence_overlap' if state else 'no_exact_sequence_overlap'] = c.class_metrics([truth[i] for i in idx], [calls[i] for i in idx], ['adhesome_associated','non_adhesome'])
    similarity_settings=_family_similarity_config(c,saved)
    profile_name,profile_values=_family_similarity_profile(similarity_settings,family_similarity_profile)
    return {'membership': binary, 'family': family, 'rule_fallback':rule_fallback,
            'rule_fallback_policy_version':'5.5_identity_rule_fallback_v1',
            'family_similarity_policy':{'profile':profile_name,'search_max_evalue':similarity_settings['search_max_evalue'],
                'thresholds':profile_values,'membership_model_uses_similarity_extension':False},
            'class_counts': dict(Counter(p['classifier_call'] for p in predictions)),
            'calibration_overlap_strata': strata, 'note': 'Validation never fitted or used to choose thresholds. Exact sequence overlap is reporting only; evolutionary homology overlap may remain. Specific and broad votes are uncalibrated.'}


def compare_validation(rows, baseline, current, saved, c):
    """Same-cohort pre/post comparison; baseline calls are outcomes, never features."""
    truth=[c.normalized_control_class(r) for r in rows]
    expected=[trusted_family(r,saved['aliases'],saved['parents']) for r in rows]
    def summarize(predictions):
        calls=['adhesome_associated' if p.get('classifier_call') in ('specific','broad','domain_only') else 'non_adhesome' for p in predictions]
        tp=sum(y=='adhesome_associated' and z=='adhesome_associated' for y,z in zip(truth,calls))
        fn=sum(y=='adhesome_associated' and z=='non_adhesome' for y,z in zip(truth,calls))
        fp=sum(y=='non_adhesome' and z=='adhesome_associated' for y,z in zip(truth,calls))
        tn=sum(y=='non_adhesome' and z=='non_adhesome' for y,z in zip(truth,calls))
        precision=tp/(tp+fp) if tp+fp else 0.0; recall=tp/(tp+fn) if tp+fn else 0.0
        f1=2*precision*recall/(precision+recall) if precision+recall else 0.0
        fam=[i for i,(y,e) in enumerate(zip(truth,expected)) if y=='adhesome_associated' and e]
        exact=sum(predictions[i].get('predicted_family','')==expected[i] for i in fam)
        specific=[i for i in fam if predictions[i].get('classifier_call')=='specific']
        broad=[i for i in fam if predictions[i].get('classifier_call')=='broad']
        broad_ok=sum(predictions[i].get('predicted_family','') in [expected[i]]+ancestors(expected[i],saved['parents']) for i in broad)
        return {'tp':tp,'fn':fn,'fp':fp,'tn':tn,'sensitivity_recall':recall,'precision':precision,'f1':f1,
            'positive_controls_correctly_recovered':tp,'labelled_family_positive_count':len(fam),
            'exact_family_correct':exact,'exact_family_accuracy':exact/len(fam) if fam else 0.0,
            'specific_assignments_with_family_label':len(specific),
            'specific_family_correct':sum(predictions[i].get('predicted_family','')==expected[i] for i in specific),
            'specific_family_accuracy':sum(predictions[i].get('predicted_family','')==expected[i] for i in specific)/len(specific) if specific else 0.0,
            'broad_assignments_with_family_label':len(broad),'broad_family_or_ancestor_correct':broad_ok,
            'broad_family_accuracy':broad_ok/len(broad) if broad else 0.0}
    old=summarize(baseline); new=summarize(current)
    return {'cohort_records':len(rows),'same_denominator':True,'baseline':old,'updated_classifier':new,
        'delta_updated_minus_baseline':{k:new[k]-v for k,v in old.items() if isinstance(v,(int,float)) and not isinstance(v,bool)},
        'baseline_role':'Comparison only; baseline calls were not used as training features or labels.'}


def read_baseline_predictions(path, rows, c):
    baseline_rows=c.read_tsv(path)
    ids=[c.row_id(row) for row in baseline_rows]
    expected=[c.row_id(row) for row in rows]
    if len(ids)!=len(set(ids)) or set(ids)!=set(expected):
        raise ValueError('Baseline validation predictions must cover the identical unique validation cohort')
    by_id={c.row_id(row):row for row in baseline_rows}
    predictions=[]
    for ident in expected:
        row=by_id[ident]
        call=str(row.get('classifier_call') or row.get('category') or '').strip().lower()
        if call in {'adhesome_associated','positive','positive_control'}:
            call='specific'
        elif call in {'non-adhesome','nonadhesome','negative','negative_control'}:
            call='non_adhesome'
        if call not in {'specific','broad','domain_only','unresolved','non_adhesome'}:
            raise ValueError(f'Unrecognized baseline classifier call for {ident}: {call!r}')
        family=str(row.get('predicted_family') or row.get('family_assignment') or '').strip()
        predictions.append({'protein_id':ident,'classifier_call':call,'predicted_family':family})
    return predictions


def write_validation_report(directory, report):
    (directory/'internal_validation_report.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    membership, family = report['membership'], report['family']
    rows = [('Validation records',membership['records']),
        ('Membership precision',f"{membership['positive_precision']:.1%}"),
        ('Membership recall',f"{membership['positive_recall']:.1%}"),
        ('Membership F1',f"{membership['positive_f1']:.1%}"),
        ('Family-labelled positives',family['labelled_positives']),
        ('Correct / accepted family calls',f"{family['correct_accepted']} / {family['accepted']}"),
        ('Strict family recall (abstentions count as missed)',f"{family['strict_recall']:.1%}"),
        ('Top-ranked candidate correct (before acceptance)',family['top_candidate_correct']),
        ('Truth within top three candidates',family['top_three_contains_truth']),
        ('Positive controls without family labels',family['unmapped_positives'])]
    def table(items):
        return '<table>'+''.join('<tr><th>'+html.escape(str(k))+'</th><td>'+html.escape(str(v))+'</td></tr>' for k,v in items)+'</table>'
    document = '<!doctype html><meta charset="utf-8"><title>Family classifier validation</title><style>body{font:16px system-ui;max-width:1000px;margin:40px auto}td,th{padding:8px;border:1px solid #ccc;text-align:left}table{border-collapse:collapse}pre{white-space:pre-wrap}</style><h1>Independent panel evaluation</h1>'
    document += '<p>'+html.escape(report['note'])+'</p>'+table(rows)
    similarity=report.get('family_similarity_policy')
    if similarity:
        policy_rows=[('Family similarity profile',similarity['profile']),
            ('CDD search ceiling',similarity['search_max_evalue']),
            ('Specific-family threshold',json.dumps(similarity['thresholds']['specific'],sort_keys=True)),
            ('Broad-family threshold',json.dumps(similarity['thresholds']['broad'],sort_keys=True)),
            ('Used to alter membership model',similarity['membership_model_uses_similarity_extension'])]
        document += '<h2>Family-resolution similarity policy</h2><p>CDD-fragment sequence matches support family resolution by E-value and coverage. The membership model continues to use its existing feature set and decision cutoff.</p>'+table(policy_rows)
    comparison=report.get('classifier_version_comparison')
    if comparison:
        metric_rows=[]
        old,new,delta=(comparison['baseline'],comparison['updated_classifier'],comparison['delta_updated_minus_baseline'])
        for key,label in (("sensitivity_recall","Sensitivity / recall"),("fn","False negatives"),
                          ("precision","Precision"),("fp","False positives"),("f1","F1"),
                          ("positive_controls_correctly_recovered","Correctly recovered positives")):
            metric_rows.append((label,f"{old[key]} | {new[key]} | Δ {delta[key]}"))
        for label,correct,count,accuracy in (
            ("Exact-family accuracy","exact_family_correct","labelled_family_positive_count","exact_family_accuracy"),
            ("Specific-family accuracy","specific_family_correct","specific_assignments_with_family_label","specific_family_accuracy"),
            ("Broad-family accuracy","broad_family_or_ancestor_correct","broad_assignments_with_family_label","broad_family_accuracy")):
            old_text=f"{old[correct]}/{old[count]} ({old[accuracy]:.1%})"
            new_text=f"{new[correct]}/{new[count]} ({new[accuracy]:.1%})"
            metric_rows.append((label,f"{old_text} | {new_text} | Δ {delta[accuracy]:+.1%}"))
        document += '<h2>Baseline vs current sequence features</h2><p>Same unique validation cohort and denominators; baseline predictions are comparison-only.</p>'+table([('Metric','Baseline | v5.4 | Change'),*metric_rows])
    document += '<p>Membership positives include specific, broad and domain-only calls. Unresolved calls count as negative for membership metrics. Family recall includes abstentions; membership performance does not establish family identity.</p><h2>Resolution tiers</h2>'+table(report['class_counts'].items())
    fallback=report.get('rule_fallback')
    if fallback:
        fallback_rows=[('Unique provisional rule candidates',fallback['unique_rule_candidate_records']),
            ('Ambiguous rule candidate records',fallback['ambiguous_rule_candidate_records']),
            ('Non-adhesome controls with a unique rule candidate',fallback['non_adhesome_controls_with_unique_rule_candidate']),
            ('Rule membership precision',f"{fallback['membership_precision']:.1%}"),
            ('Rule membership recall',f"{fallback['membership_recall']:.1%}"),
            ('Family-labelled positives with a unique rule candidate',fallback['family_labelled_positive_controls_with_unique_rule_candidate']),
            ('Exact-family rule matches',fallback['exact_family_correct']),
            ('Strict exact-family rule recall',f"{fallback['strict_exact_family_recall']:.1%}"),
            ('Untrained-family exact rule recoveries',fallback['untrained_family_exact_rule_recoveries']),
            ('Conflicts with RF family calls',fallback['rule_rf_conflicts'])]
        document += '<h2>Separate provisional family-identity-rule pathway</h2><p>'+html.escape(fallback['interpretation'])+'</p>'+table(fallback_rows)
    document += '<h2>Per-family counts</h2>'+table((label,f"{v['correct']} correct / {v['records']} labelled") for label,v in family['per_family'].items())
    document += '<details><summary>Full metrics and overlap strata</summary><pre>'+html.escape(json.dumps(report,indent=2))+'</pre></details>'
    (directory/'internal_validation_report.html').write_text(document,encoding='utf-8')


def train_family_gate(args, c):
    print('[classifier] Loading family rules, calibration labels, and verified Stage 1 evidence',flush=True)
    aliases, parents, motifs = load_resources(c)
    similarity_settings=c._SETTINGS['family_gate']['family_similarity']
    similarity_profile_name,similarity_selection=_family_similarity_profile(similarity_settings)
    vocab, provenance = c.configured_vocab(c._BUNDLE)
    rows, sequences, receipts = read_cohort(args.calibration, c)
    library=c._BUNDLE/'resource_library/adhesome_ml/cdd_domain_library_v1'
    with tempfile.TemporaryDirectory(prefix='adhesome_cdd_train_') as temp:
        print('[classifier] Building SHA-256-verified calibration-positive CDD reference bank',flush=True)
        family_labels={c.row_id(row):trusted_family(row,aliases,parents) for row in rows}
        bank_fasta,bank_meta=build_calibration_bank(library,args.calibration,Path(temp),family_labels=family_labels)
        bank_manifest=json.loads((Path(temp)/'cdd_calibration_references.manifest.json').read_text(encoding='utf-8'))
        print(f"[classifier] CDD bank ready: {bank_manifest['unique_reference_fragments']} fragments; "
              f"{bank_manifest['traceable_calibration_hit_records']} calibration hit records; "
              f"{bank_manifest['unique_cdd_accessions']} accessions; validation references excluded",flush=True)
        calibration_queries=[{'protein_id':c.row_id(r),'sequence':sequences[c.row_id(r)]} for r in rows]
        print(f'[classifier] Extracting calibration CDD features for {len(calibration_queries)} controls',flush=True)
        train_cdd=similarity_features(calibration_queries,bank_fasta,bank_meta,exclude_self=True,
            search_evalue=similarity_settings['search_max_evalue'],
            base_feature_evalue=similarity_settings['legacy_model_feature_max_evalue'])
        print('[classifier] Extracting calibration motif, Pfam, architecture, identity-rule, post-hoc-rule, and homology features',flush=True)
        data, scan_hashes = extract(rows, sequences, args.calibration/'handoff/detection', motifs, vocab, c,cdd=train_cdd)
        frozen_cdd_fasta=bank_fasta.read_bytes()
        frozen_cdd_meta=(Path(temp)/'cdd_calibration_references.json').read_bytes()
        frozen_cdd_manifest=(Path(temp)/'cdd_calibration_references.manifest.json').read_bytes()
        val_rows=val_sequences=val_receipts=val_hashes=val_data=None
        if args.validation:
            val_rows,val_sequences,val_receipts=read_cohort(args.validation,c)
            val_queries=[{'protein_id':c.row_id(r),'sequence':val_sequences[c.row_id(r)]} for r in val_rows]
            print(f'[classifier] Extracting independent validation CDD features for {len(val_queries)} controls; validation does not enter the reference bank',flush=True)
            val_cdd=similarity_features(val_queries,bank_fasta,bank_meta,
                search_evalue=similarity_settings['search_max_evalue'],
                base_feature_evalue=similarity_settings['legacy_model_feature_max_evalue'])
            print('[classifier] Extracting independent validation sequence features',flush=True)
            val_data,val_hashes=extract(val_rows,val_sequences,args.validation/'handoff/detection',motifs,vocab,c,cdd=val_cdd)
    labels = [c.normalized_control_class(r) for r in rows]
    if set(labels) != {'adhesome_associated','non_adhesome'}:
        raise ValueError('Both positive and negative calibration controls are required')
    def estimator(seed):
        settings = dict(c._SETTINGS['binary_model'], random_state=seed)
        return c.Pipeline([('vectorizer',c.DictVectorizer()),('model',c.random_forest(settings,args.jobs))])
    print('[classifier] Fitting adhesome membership model',flush=True)
    binary = estimator(1729).fit([membership_model_vector(v) for _,v in data], labels)
    indices = [(i,trusted_family(row,aliases,parents)) for i,row in enumerate(rows) if labels[i]=='adhesome_associated']
    indices = [(i,label) for i,label in indices if label and any(c.evidence_groups(data[i][0]).values())]
    if len({label for _,label in indices}) < 2:
        raise ValueError('At least two trusted, feature-supported calibration families are required')
    counts = Counter(label for _,label in indices)
    model_indices=[(i,label) for i,label in indices if counts[label]>=2]
    print(f'[classifier] Fitting family model from {len(model_indices)} labelled, feature-supported controls',flush=True)
    family = estimator(2718).fit([data[i][1] for i,_ in model_indices], [label for _,label in model_indices])
    profiles = {label: sorted(set().union(*(domain_set(data[i][0],c) for i,f in indices if f==label))) for label in counts}
    saved = dict(version='5.6',feature_schema=SCHEMA,binary=binary,family=family,config_vocab=sorted(vocab),
        aliases=aliases,parents=parents,motifs=motifs,counts=dict(counts),profiles=profiles,
        family_model_counts=dict(Counter(label for _,label in model_indices)),
        single_example_family_profiles=sorted(label for label,n in counts.items() if n==1),
        policy=c._SETTINGS['family_gate'],family_similarity_settings=similarity_settings,
        calibration_sequence_hashes={sequence_digest(sequences[c.row_id(r)]) for r in rows},
        provenance={**provenance,**receipts,**scan_hashes,'labels':c.sha(args.calibration/'reports/predictions.tsv'),
            'cdd_domain_fasta':c.sha(library/'cdd_domain_sequences.faa'),
            'cdd_domain_metadata':c.sha(library/'cdd_domain_metadata.tsv'),
            'family_identity_rules':c.sha(c._BUNDLE/'config/annotation/family_identity_rules.json'),
            'posthoc_sequence_rules_sha256':scan_hashes['posthoc_sequence_rules_sha256']})
    summary = {'classifier_version':'5.6','feature_schema':SCHEMA,'calibration_records':len(rows),
        'family_training_records':len(model_indices),'family_profile_records':len(indices),
        'family_model_counts':dict(Counter(label for _,label in model_indices)),
        'single_example_family_profiles':sorted(label for label,n in counts.items() if n==1),
        'family_profile_counts':dict(counts),
        'unmapped_or_ambiguous_positive_labels':sum(x=='adhesome_associated' for x in labels)-len(indices),
        'family_labels_excluded_from_multiclass_fit':sum(1 for n in counts.values() if n<2),
        'policy':saved['policy'],'validation_used_for_fitting':False,
        'feature_contract':'Identical FASTA motif scans, verified Pfam hits, calibration CDD local similarities, domain order/copy/quality, directional homology, and sequence-only boolean/token features parsed from the post-hoc rule conditions; no post-hoc proposals, audit scores, or hypothesis labels.',
        'identity_rule_policy':'Only measurable sequence fields from family_identity_rules.json are added: required domain-group counts/satisfaction, specific CDD annotation pattern matches, ordered groups, maximum domain counts, repeat indices and protein length ranges. Topology, orthology, localization and signal peptide are excluded from this feature block.',
        'posthoc_rule_feature_policy':'Parse only the boolean conditions in specific() and broad() from schisto_pfam_cdd_reclassify.py. Include rule-condition hits, individual Pfam/CDD signature presence and counts, and concurrent-signal counts. Never read rule return values, proposed family labels, tiers, orthogroup arguments, post-hoc assignment tables, or validation outcomes as features.',
        'cdd_policy':'All SHA-256-verified CDD fragments from labelled calibration-positive controls are retained, including fragments without a mapped family label. Validation-origin fragments are excluded. Calibration queries exclude same-protein references. Query and subject coverage use inclusive BLAST HSP coordinate spans divided by full sequence lengths, capped at 1.0.',
        'optional_returned_policy':'Returned CDD hits are not required; candidate sequences are searched against the frozen calibration-only CDD fragment bank.'}
    summary['family_similarity_policy']={'selected_profile':similarity_profile_name,
        'search_max_evalue':similarity_settings['search_max_evalue'],
        'legacy_model_feature_max_evalue':similarity_settings['legacy_model_feature_max_evalue'],
        'profile_thresholds':similarity_selection,
        'membership_model_uses_family_similarity_extension':False,
        'rule_override_key':similarity_settings.get('rule_override_key')}
    val_predictions = []
    if args.validation:
        print('[classifier] Evaluating the frozen calibration-trained model on independent validation',flush=True)
        val_predictions = infer(saved,val_data,c,family_similarity_profile=similarity_profile_name)
        summary['validation'] = evaluate(val_rows,val_predictions,val_sequences,saved,c,
            family_similarity_profile=similarity_profile_name)
        summary['validation_provenance'] = {**val_receipts,**val_hashes,'labels':c.sha(args.validation/'reports/predictions.tsv')}
        if args.baseline_validation_predictions:
            baseline_predictions=read_baseline_predictions(args.baseline_validation_predictions,val_rows,c)
            summary['validation']['classifier_version_comparison']=compare_validation(
                val_rows,baseline_predictions,val_predictions,saved,c)
            summary['baseline_validation_predictions_sha256']=c.sha(args.baseline_validation_predictions)
        for row,prediction in zip(val_rows,val_predictions):
            prediction.update(expected_family=trusted_family(row,aliases,parents),control_class=row['control_class'],
                calibration_sequence_overlap=sequence_digest(val_sequences[c.row_id(row)]) in saved['calibration_sequence_hashes'])
    args.out.mkdir(parents=True,exist_ok=False)
    (args.out/'cdd_calibration_references.faa').write_bytes(frozen_cdd_fasta)
    (args.out/'cdd_calibration_references.json').write_bytes(frozen_cdd_meta)
    (args.out/'cdd_calibration_references.manifest.json').write_bytes(frozen_cdd_manifest)
    saved['provenance']['frozen_cdd_reference_fasta']=c.sha(args.out/'cdd_calibration_references.faa')
    saved['provenance']['frozen_cdd_reference_metadata']=c.sha(args.out/'cdd_calibration_references.json')
    saved['provenance']['frozen_cdd_reference_manifest']=c.sha(args.out/'cdd_calibration_references.manifest.json')
    summary['cdd_reference_bank']=json.loads(frozen_cdd_manifest.decode('utf-8'))
    (args.out/'model_summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
    c.joblib.dump(saved,args.out/'adhesome_ml.joblib')
    (args.out/'model_summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
    c.write_tsv(args.out/'family_training_labels.tsv',[{'protein_id':c.row_id(rows[i]),'family':label} for i,label in indices])
    if val_predictions:
        c.write_tsv(args.out/'validation_predictions.tsv',val_predictions)
        report = summary['validation']
        write_validation_report(args.out, report)
    print(f'[classifier] COMPLETE: model and independent validation report written to {args.out}',flush=True)
    print(json.dumps({'model_directory':str(args.out),'calibration_records':len(rows),
        'learned_families':len(counts),'validation_records':len(val_predictions),
        'reports':'internal_validation_report.json and internal_validation_report.html'},indent=2))
    return 0


def evaluate_family_gate(args, c):
    """Evaluate a frozen model; never fit or select thresholds on validation."""
    if args.out.exists():
        raise ValueError('Choose a new evaluation output directory')
    path = args.model/'adhesome_ml.joblib' if args.model.is_dir() else args.model
    saved = c.joblib.load(path)
    if saved.get('feature_schema') not in ({SCHEMA}|LEGACY_SCHEMAS):
        raise ValueError('Evaluation requires a supported frozen CDD-bank family classifier (v5.4 or later)')
    similarity_settings=_family_similarity_config(c,saved)
    similarity_profile_name,_=_family_similarity_profile(
        similarity_settings,getattr(args,'family_similarity_profile',None))
    rows, sequences, receipts = read_cohort(args.validation, c)
    bank_path=args.model if args.model.is_dir() else path.parent
    bank_fasta=bank_path/'cdd_calibration_references.faa'; meta_path=bank_path/'cdd_calibration_references.json'
    manifest_path=bank_path/'cdd_calibration_references.manifest.json'
    if not bank_fasta.is_file() or not meta_path.is_file() or not manifest_path.is_file(): raise ValueError('Model lacks its frozen, verified calibration-only CDD reference bank; retrain v5.4')
    if c.sha(bank_fasta)!=saved.get('provenance',{}).get('frozen_cdd_reference_fasta') or c.sha(meta_path)!=saved.get('provenance',{}).get('frozen_cdd_reference_metadata') or c.sha(manifest_path)!=saved.get('provenance',{}).get('frozen_cdd_reference_manifest'):
        raise ValueError('Frozen CDD reference bank differs from the bank recorded in the model')
    bank_meta=json.loads(meta_path.read_text(encoding='utf-8'))
    bank_manifest=json.loads(manifest_path.read_text(encoding='utf-8'))
    if bank_manifest.get('validation_reference_records_included') is not False:
        raise ValueError('CDD bank contains validation-origin references; refusing to use it')
    queries=[{'protein_id':c.row_id(r),'sequence':sequences[c.row_id(r)]} for r in rows]
    cdd=similarity_features(queries,bank_fasta,bank_meta,
        search_evalue=similarity_settings['search_max_evalue'],
        base_feature_evalue=similarity_settings['legacy_model_feature_max_evalue'])
    data, hashes = extract(rows, sequences, args.validation/'handoff/detection', saved['motifs'], set(saved['config_vocab']), c,cdd=cdd)
    if hashes.get('posthoc_sequence_rules_sha256') != saved.get('provenance',{}).get('posthoc_sequence_rules_sha256'):
        raise ValueError('Post-hoc sequence-rule source differs from the one used to train this model')
    predictions = infer(saved, data, c,family_similarity_profile=similarity_profile_name)
    report = evaluate(rows, predictions, sequences, saved, c,family_similarity_profile=similarity_profile_name)
    for row, prediction in zip(rows, predictions):
        prediction.update(expected_family=trusted_family(row,saved['aliases'],saved['parents']),
            control_class=row['control_class'],
            calibration_sequence_overlap=sequence_digest(sequences[c.row_id(row)]) in saved['calibration_sequence_hashes'])
    args.out.mkdir(parents=True,exist_ok=False)
    c.write_tsv(args.out/'validation_predictions.tsv', predictions)
    write_validation_report(args.out, report)
    receipt = {'model_sha256':c.sha(path),'labels_sha256':c.sha(args.validation/'reports/predictions.tsv'),
        'family_similarity_profile':similarity_profile_name,
        'family_similarity_policy':report['family_similarity_policy'],
        'records':len(rows),'validation_used_for_fitting':False,'sources':{**receipts,**hashes}}
    (args.out/'evaluation_manifest.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
    print(f"Validation report: {args.out/'internal_validation_report.html'}; records={len(rows)}; model unchanged")
    return 0


def predict_family_gate(args, c):
    print('[classifier] Starting family prediction; loading model and checking inputs',flush=True)
    if not args.fasta or not args.evidence_dir:
        raise ValueError('Family prediction requires --fasta and --evidence-dir for verified sequence-derived features')
    saved = c.joblib.load(args.model/'adhesome_ml.joblib' if args.model.is_dir() else args.model)
    if saved.get('feature_schema') not in ({SCHEMA}|LEGACY_SCHEMAS):
        raise ValueError('Prediction requires a supported frozen CDD-bank family classifier (v5.4 or later)')
    similarity_settings=_family_similarity_config(c,saved)
    similarity_profile_name,similarity_selection=_family_similarity_profile(similarity_settings,args.family_similarity_profile)
    records = c.read_fasta(args.fasta)
    if args.ids:
        requested = set(args.ids.read_text().splitlines())
        if requested - {r['protein_id'] for r in records}:
            raise ValueError('Requested IDs are absent from candidate FASTA')
        records = [r for r in records if r['protein_id'] in requested]
    if len(records)>250 and not args.batch:
        raise ValueError('Maximum 250 sequences; specify --batch for consecutive 250-record batches')
    print(f'[classifier] Model ready; {len(records)} sequences queued',flush=True)
    sequences, _ = verified_sequences(args.evidence_dir,c)
    for index,row in enumerate(records,1):
        if sequences.get(row['protein_id']) != row['sequence']:
            raise ValueError(f"Candidate sequence differs from its scanned FASTA: {row['protein_id']}")
        if index%100==0 or index==len(records):
            print(f'[classifier] Verified Stage 1 sequence evidence: {index}/{len(records)}',flush=True)
    output = []
    cached = c.read_stage1_evidence(args.evidence_dir)
    posthoc_script=c._BUNDLE.parent/'schisto_pfam_cdd_reclassify/schisto_pfam_cdd_reclassify.py'
    current_posthoc_sha,_,_,_=load_posthoc_rules(str(posthoc_script))
    if current_posthoc_sha != saved.get('provenance',{}).get('posthoc_sequence_rules_sha256'):
        raise ValueError('Post-hoc sequence-rule source differs from the one used to train this model')
    bank_path=args.model if args.model.is_dir() else args.model.parent
    bank_fasta=bank_path/'cdd_calibration_references.faa'; meta_path=bank_path/'cdd_calibration_references.json'
    manifest_path=bank_path/'cdd_calibration_references.manifest.json'
    if not bank_fasta.is_file() or not meta_path.is_file() or not manifest_path.is_file(): raise ValueError('Model lacks its frozen, verified calibration-only CDD reference bank; retrain v5.4')
    if c.sha(bank_fasta)!=saved.get('provenance',{}).get('frozen_cdd_reference_fasta') or c.sha(meta_path)!=saved.get('provenance',{}).get('frozen_cdd_reference_metadata') or c.sha(manifest_path)!=saved.get('provenance',{}).get('frozen_cdd_reference_manifest'):
        raise ValueError('Frozen CDD reference bank differs from the bank recorded in the model')
    bank_meta=json.loads(meta_path.read_text(encoding='utf-8'))
    bank_manifest=json.loads(manifest_path.read_text(encoding='utf-8'))
    if bank_manifest.get('validation_reference_records_included') is not False:
        raise ValueError('CDD bank contains validation-origin references; refusing to use it')
    hypotheses = defaultdict(set)
    hypothesis_file = args.evidence_dir.parent/'hypotheses.tsv'
    if hypothesis_file.is_file():
        for hypothesis in c.read_tsv(hypothesis_file):
            hypotheses[hypothesis.get('protein_id','')].add(hypothesis.get('family',''))
    total_batches=(len(records)+249)//250
    for batch_index,start in enumerate(range(0,len(records),250),1):
        chunk=records[start:start+250]
        batch_end=start+len(chunk)
        print(f'[classifier] Batch {batch_index}/{total_batches} ({start+1}-{batch_end}/{len(records)}): comparing CDD domain profiles',flush=True)
        cdd=similarity_features(chunk,bank_fasta,bank_meta,
            search_evalue=similarity_settings['search_max_evalue'],
            base_feature_evalue=similarity_settings['legacy_model_feature_max_evalue'])
        print(f'[classifier] Batch {batch_index}/{total_batches}: extracting motif, domain, identity-rule, post-hoc-rule and homology features',flush=True)
        data,_ = extract(chunk,sequences,args.evidence_dir,saved['motifs'],set(saved['config_vocab']),c,cached,cdd)
        print(f'[classifier] Batch {batch_index}/{total_batches}: assigning adhesome membership and family calls',flush=True)
        output.extend(infer(saved,data,c,args.threshold,similarity_profile_name))
        print(f'[classifier] Progress: {len(output)}/{len(records)} sequences classified ({100*len(output)/len(records):.1f}%)',flush=True)
    if len(output)!=len(records):
        raise ValueError('Prediction cohort accounting failed')
    for prediction in output:
        prediction['stage1_hypotheses_context_only'] = ';'.join(sorted(hypotheses[prediction['protein_id']] - {''}))
    c.write_tsv(args.out,output)
    receipt = {'classifier_version':'5.6','model_version':saved.get('version','unknown'),
               'rule_fallback_policy_version':'5.5_identity_rule_fallback_v1',
               'family_similarity_profile':similarity_profile_name,
               'family_similarity_search_max_evalue':similarity_settings['search_max_evalue'],
               'family_similarity_assignment_thresholds':similarity_selection,
               'membership_model_uses_family_similarity_extension':False,
               'feature_schema':SCHEMA,'records':len(output),'class_counts':dict(Counter(r['classifier_call'] for r in output)),
               'rule_candidate_status_counts':dict(Counter(r['rule_candidate_status'] for r in output)),
               'untrained_family_rule_candidates':sum(bool(r['rule_candidate_family_untrained']) for r in output),
               'rule_rf_conflicts':sum(bool(r['rule_conflicts_with_random_forest']) for r in output),
               'model_sha256':c.sha(args.model/'adhesome_ml.joblib' if args.model.is_dir() else args.model),
               'input_sha256':c.sha(args.fasta),'threshold':args.threshold,'exact_sequence_match_used':False}
    args.out.with_suffix('.summary.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
    print(f'[classifier] COMPLETE: wrote {len(output)} predictions to {args.out}',flush=True)
    print(json.dumps(receipt,indent=2))
    return 0
