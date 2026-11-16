"""FASTA-derived family learning. Validation labels never enter fitting or thresholds."""
import hashlib
import html
import json
import re
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

SCHEMA = 'sequence_family_gate_v5_8_calibrated_probabilities'
from cdd_features import build_calibration_bank, similarity_features
from family_signature_features import extract as crosswalk_signature_features, load_library as load_crosswalk_library
from identity_rule_features import features as identity_rule_features, load_rules as load_identity_rules
from posthoc_rule_features import extract as posthoc_sequence_features, load_rules as load_posthoc_rules
from probability_calibration import (CalibrationUnavailable, apply_binary_sigmoid,
    apply_family_sigmoids, cross_fitted_probabilities, fit_binary_sigmoid,
    fit_family_sigmoids, probability_metrics_binary, sequence_cluster_groups,
    probability_metrics_multiclass)


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
    family_rule_config=json.loads((c._BUNDLE/'config/annotation/family_rules.json').read_text(encoding='utf-8'))
    family_aliases={}
    for field in ('family_aliases','canonical_family_aliases'):
        for key,value in family_rule_config.get(field,{}).items():
            if isinstance(value,str):family_aliases[key.casefold()]=value
    crosswalk_signatures,crosswalk_metadata=_family_signature_library(c,family_aliases)
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
        signature_features,signature_evidence=crosswalk_signature_features(
            row,vector,crosswalk_signatures,c._SETTINGS['family_gate']['family_similarity']['search_max_evalue'])
        vector.update(signature_features)
        row['_family_signature_evidence_json']=json.dumps(signature_evidence,sort_keys=True)
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
    hashes['adhesome_family_rule_crosswalk_sha256']=crosswalk_metadata['sha256']
    hashes['family_rules']=c.sha(c._BUNDLE/'config/annotation/family_rules.json')
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
    return {key:value for key,value in vector.items()
            if not key.startswith(('family_similarity:','family_signature:'))}


def _family_token(value):
    return re.sub(r'[^A-Za-z0-9_.-]+','_',str(value or '')).strip('_').casefold()


def _family_similarity_config(c, saved=None):
    return (saved or {}).get('family_similarity_settings') or c._SETTINGS['family_gate']['family_similarity']


def _family_signature_library(c, aliases):
    path=c._BUNDLE/c._SETTINGS['family_signature_crosswalk']['path']
    identity_rule_path=c._BUNDLE/'config/annotation/family_identity_rules.json'
    identity_rules,_,_=load_identity_rules(str(identity_rule_path))
    return load_crosswalk_library(path,identity_rules,aliases)


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


def _identity_core_observed(rule, row, c):
    """Check the shared configured domain core, including repeat counts and order."""
    records=c._domain_records(row)
    group_hits=[]
    required_groups=rule.get('core_required_groups') or rule.get('required_groups',[])
    for group in required_groups or []:
        models={_family_policy_token(x) for x in group.get('models',[]) if str(x).strip()}
        matches=[record for record in records
                 if _family_policy_token(str(record.get('name','')).split('.')[0]) in models]
        required=int(group.get('count',1))
        coordinates={(record.get('start'),record.get('end')) for record in matches
                     if record.get('start') not in (None,'') and record.get('end') not in (None,'')}
        observed_count=len(coordinates) if coordinates else len(matches)
        if observed_count<required:
            return False
        group_hits.append(matches)
    order=rule.get('ordered_groups',[]) or []
    if order:
        if not all(isinstance(index,int) and 0<=index<len(group_hits) and group_hits[index]
                   for index in order):
            return False
        positions=[]
        for index in order:
            starts=[float(hit['start']) for hit in group_hits[index]
                    if hit.get('start') not in (None,'')]
            positions.append(min(starts) if starts else None)
        if not all(value is not None for value in positions) or not all(
                left<=right for left,right in zip(positions,positions[1:])):
            return False
    return True


def _family_label_tokens(label, aliases, assignment_policy):
    """Keep raw and configured alias labels for rule/class reconciliation."""
    values={_family_policy_token(label)}
    values.add(_family_policy_token(canonical(str(label),aliases or {})))
    _,_,policy_aliases,_,_=_family_policy_maps(assignment_policy or {})
    for token in list(values):
        if token in policy_aliases:
            values.add(_family_policy_token(policy_aliases[token]))
    return {value for value in values if value}


def _shared_core_family_resolution(top, ranking, row, c, assignment_policy,
                                   identity_rules, aliases, rule_source_families,
                                   similarity_settings, similarity_selection):
    """Prefer a configured broad parent when shared-core families lack a unique CDD call."""
    settings=(assignment_policy or {}).get('shared_core_resolution',{})
    if not settings.get('enabled',False):
        return None
    top_tokens=_family_label_tokens(top,aliases,assignment_policy)
    raw_thresholds=(assignment_policy or {}).get('raw_vote_thresholds',{})
    broad_cutoff=float(settings.get('minimum_broad_vote_fraction',
        raw_thresholds.get('broad_parent_vote_fraction_cutoff',0.6)))
    allowed,_,_,_,_=_family_policy_maps(assignment_policy or {})
    proposal_tokens={_family_policy_token(x) for x in (rule_source_families or set())}
    for group in settings.get('groups',[]) or []:
        broad_label=str(group.get('broad_label','')).strip()
        broad_token=_family_policy_token(broad_label)
        member_tokens={_family_policy_token(x) for x in group.get('members',[]) if str(x).strip()}
        if not broad_label or broad_token not in allowed or not member_tokens:
            continue
        if not (top_tokens & member_tokens):
            continue
        core_rule={'core_required_groups':group.get('core_required_groups',[]),
                   'ordered_groups':group.get('ordered_groups',[])}
        if not core_rule['core_required_groups'] or not _identity_core_observed(core_rule,row,c):
            continue
        collision_members=member_tokens

        # Identity-rule proposals encode their family-specific CDD pattern. A
        # single supported child can resolve the shared domain core; zero or
        # multiple child proposals remain broad.
        cdd_supported={token for token in collision_members if token in proposal_tokens}
        if similarity_settings and similarity_selection:
            for token in collision_members:
                family=next((str(rule.get('family') or rule.get('family_set') or '').strip()
                             for rule in (identity_rules or [])
                             if _family_policy_token(rule.get('family') or rule.get('family_set'))==token),
                            token)
                if family and _family_similarity_match(row,family,'specific',similarity_settings,
                        similarity_selection,identity_rules,aliases):
                    cdd_supported.add(token)
        top_members=top_tokens & collision_members
        if len(cdd_supported)==1 and top_members & cdd_supported:
            continue

        group_vote=sum(float(vote) for label,vote in ranking
                       if _family_label_tokens(label,aliases,assignment_policy)&collision_members)
        reason=('shared_core_cdd_signature_conflicts_with_top_family' if cdd_supported else
                'shared_core_cdd_signature_not_differentiating')
        detail=(f"shared core {','.join(sorted(collision_members))}; "
                f"CDD-supported={','.join(sorted(cdd_supported)) or 'none'}; "
                f"group raw-vote mass={group_vote:.3f}")
        return {'suppress_specific':True,
                'broad_call':('broad',broad_label,f'{reason}; {detail}')
                    if group_vote>=broad_cutoff else None}
    return None


def decide(membership, ranking, row, profiles, counts, parents, policy, c,
           vector=None, identity_rules=None, aliases=None, rule_families=None,
           similarity_settings=None, similarity_selection=None,
           assignment_policy=None, rule_source_families=None):
    vector=vector or {}; identity_rules=identity_rules or []; aliases=aliases or {}; rule_families=set(rule_families or [])
    observed = domain_set(row, c)
    if not any(c.evidence_groups(row).values()):
        return 'unresolved', '', 'no_informative_features'
    if membership < policy['membership_probability_cutoff']:
        return 'non_adhesome', '', 'membership_below_cutoff'
    if not ranking:
        return ('domain_only' if observed else 'unresolved'), '', 'family_probability_unavailable'
    top, score = ranking[0]
    margin = score - (ranking[1][1] if len(ranking) > 1 else 0)
    shared = observed & set(profiles[top])
    shared_core_resolution=_shared_core_family_resolution(top,ranking,row,c,assignment_policy,
        identity_rules,aliases,rule_source_families,similarity_settings,similarity_selection)
    if shared_core_resolution and shared_core_resolution.get('broad_call'):
        return shared_core_resolution['broad_call']
    specific_similarity=(_family_similarity_match(row,top,'specific',similarity_settings,similarity_selection,
        identity_rules,aliases) if similarity_settings and similarity_selection else None)
    rule_specific=top in rule_families
    if (not shared_core_resolution or not shared_core_resolution.get('suppress_specific')) and (shared or specific_similarity or rule_specific) and counts[top] >= policy['minimum_family_examples'] and score >= policy['specific_family_probability_cutoff'] and margin >= policy['specific_family_probability_margin']:
        support_reason=('family_vote_margin_and_identity_rule_support' if rule_specific else
                        'family_vote_margin_and_evalue_similarity_support' if specific_similarity else
                        'family_vote_margin_and_shared_domain_support')
        if any(top in ancestors(x, parents) for x in counts if x != top):
            return 'broad', top, support_reason
        return 'specific', top, support_reason
    aggregated = defaultdict(float)
    policy_allowed,_,policy_aliases,policy_children,policy_broad_labels=_family_policy_maps(assignment_policy or {})
    def parent_child_in_scope(parent,label):
        parent_token=_family_policy_token(policy_aliases.get(_family_policy_token(parent),parent))
        if not assignment_policy or parent_token not in policy_broad_labels:
            return True
        children=policy_children.get(parent_token,set())
        return bool(children and (_family_label_tokens(label,aliases,assignment_policy)&children))
    for label, vote in ranking:
        for parent in ancestors(label, parents):
            if parent_child_in_scope(parent,label):
                aggregated[parent] += vote
    for parent, vote in sorted(aggregated.items(), key=lambda x: (-x[1], x[0])):
        parent_token=_family_policy_token(policy_aliases.get(_family_policy_token(parent),parent))
        if assignment_policy and parent_token in policy_broad_labels and parent_token not in policy_allowed:
            continue
        related = [name for name, _ in ranking if parent in ancestors(name, parents)
                   and parent_child_in_scope(parent,name)]
        domain_support=any(observed & set(profiles[name]) for name in related)
        similarity_support=bool(similarity_settings and similarity_selection and any(
            _family_similarity_match(row,name,'broad',similarity_settings,similarity_selection,
                identity_rules,aliases) for name in related))
        rule_support=any(name in rule_families for name in related)
        if vote >= policy['broad_parent_probability_cutoff'] and (domain_support or similarity_support or rule_support):
            reason=('shared_parent_vote_and_identity_rule_support' if rule_support else
                    'shared_parent_vote_and_evalue_similarity_support' if similarity_support else
                    'shared_parent_vote_and_domain_support')
            return 'broad', parent, reason
    return ('domain_only' if observed else 'unresolved'), '', 'insufficient_family_resolution'


def _family_policy_token(value):
    return re.sub(r'[^a-z0-9]+', '', str(value or '').casefold())


def _policy_digest(value):
    payload=json.dumps(value,sort_keys=True,separators=(',',':')).encode('utf-8')
    return hashlib.sha256(payload).hexdigest()


def _family_policy_maps(policy):
    allowed, tiers = {}, {}
    for tier, labels in (policy.get('allowed_families') or {}).items():
        for label in labels:
            token = _family_policy_token(label)
            allowed[token] = str(label)
            tiers[token] = str(tier)
    aliases = {}
    for alias, target in (policy.get('label_aliases') or {}).items():
        target_token = _family_policy_token(target)
        if target_token in allowed:
            aliases[_family_policy_token(alias)] = allowed[target_token]
    children = {}
    for parent, labels in (policy.get('broad_family_children') or {}).items():
        parent_token = _family_policy_token(parent)
        children[parent_token] = {_family_policy_token(label) for label in labels}
    broad_labels = {_family_policy_token(label) for label in policy.get('broad_family_labels', [])}
    return allowed, tiers, aliases, children, broad_labels


def _family_call_compatible(predicted, expected, parents, policy):
    if not predicted or not expected:
        return False
    _, _, aliases, children, _ = _family_policy_maps(policy)
    pred_token = _family_policy_token(aliases.get(_family_policy_token(predicted), predicted))
    exp_token = _family_policy_token(aliases.get(_family_policy_token(expected), expected))
    return (pred_token == exp_token or exp_token in children.get(pred_token, set()) or
            predicted in ancestors(expected, parents))


def _resolve_family_call(category, assignment, reason, parents, policy):
    """Constrain the evidence-gated raw-vote call to the configured family scope."""
    if category not in {'specific', 'broad'} or not assignment:
        return None
    allowed, tiers, aliases, children, broad_labels = _family_policy_maps(policy)
    original_token = _family_policy_token(assignment)
    token = _family_policy_token(aliases.get(original_token, assignment))
    if token in allowed:
        label = allowed[token]
        selected_category = 'broad' if token in broad_labels else category
        return selected_category, label, reason, tiers[token], assignment
    for parent_token, child_tokens in children.items():
        if token in child_tokens and parent_token in allowed:
            label = allowed[parent_token]
            return 'broad', label, f'{reason}; collapsed_to_allowed_broad_family', tiers[parent_token], assignment
    for parent in ancestors(assignment, parents):
        parent_token = _family_policy_token(aliases.get(_family_policy_token(parent), parent))
        if parent_token in allowed:
            label = allowed[parent_token]
            return 'broad', label, f'{reason}; collapsed_to_allowed_parent_family', tiers[parent_token], assignment
    return None


def _family_call_vote(assignment, ranking, parents, policy):
    """Return the raw RF vote mass represented by a selected family or broad parent."""
    if not assignment:
        return 0.0
    _, _, aliases, children, _ = _family_policy_maps(policy)
    assignment_token = _family_policy_token(aliases.get(_family_policy_token(assignment), assignment))
    child_tokens = children.get(assignment_token, set())
    matching = []
    for label, vote in ranking:
        label_token = _family_policy_token(aliases.get(_family_policy_token(label), label))
        if (label_token == assignment_token or label_token in child_tokens or
                assignment in ancestors(label, parents)):
            matching.append(float(vote))
    return sum(matching) if matching else 0.0


def _family_policy_metrics(rows, predictions, truth, expected_labels, saved, c, policy):
    """Report scoped family calls without treating raw votes as probabilities."""
    allowed, tiers, aliases, children, _ = _family_policy_maps(policy)
    child_to_parent = {child: parent for parent, child_set in children.items() for child in child_set}

    def canonical_scope_token(label):
        token = _family_policy_token(label)
        token = _family_policy_token(aliases.get(token, label))
        return child_to_parent.get(token, token)

    def tier_for_label(label):
        token = canonical_scope_token(label)
        return tiers.get(token)

    labelled = [(i, label) for i, label in enumerate(expected_labels)
                if label and truth[i]=='adhesome_associated']
    scope_labelled = [(i, label) for i, label in labelled if tier_for_label(label)]
    family_calls = [(i, p.get('predicted_family', '')) for i, p in enumerate(predictions)
                    if p.get('predicted_family')]

    compatible_scope = [(i, label) for i, label in scope_labelled
                        if _family_call_compatible(predictions[i].get('predicted_family', ''),
                                                   label,saved['parents'],policy)]
    exact_scope = [(i, label) for i, label in scope_labelled
                   if _family_policy_token(aliases.get(_family_policy_token(predictions[i].get('predicted_family', '')),
                                                       predictions[i].get('predicted_family', ''))) ==
                      _family_policy_token(aliases.get(_family_policy_token(label), label))]
    incorrect_calls = [(i, label) for i, label in family_calls
                       if truth[i] != 'adhesome_associated' or
                       (expected_labels[i] and not _family_call_compatible(label,expected_labels[i],saved['parents'],policy)) or
                       (truth[i] == 'adhesome_associated' and not expected_labels[i])]
    results = {
        'version': policy.get('version', ''),
        'mode': policy.get('mode', ''),
        'score_interpretation': policy.get('score_interpretation', ''),
        'labelled_family_controls': len(labelled),
        'in_scope_labelled_family_controls': len(scope_labelled),
        'selected_family_calls': len(family_calls),
        'exact_correct_in_scope': len(exact_scope),
        'family_compatible_correct_in_scope': len(compatible_scope),
        'incorrect_or_unlabelled_family_calls': len(incorrect_calls),
        'negative_control_family_calls': sum(truth[i] == 'non_adhesome' for i, _ in family_calls),
        'out_of_scope_raw_assignments_suppressed': sum(
            bool(p.get('family_call_policy_suppressed_outside_scope')) for p in predictions),
        'scope_precision_family_compatible': len(compatible_scope) / len(family_calls) if family_calls else 0.0,
        'scope_recall_family_compatible': len(compatible_scope) / len(scope_labelled) if scope_labelled else 0.0,
        'overall_labelled_exact_recall': len(exact_scope) / len(labelled) if labelled else 0.0,
        'tiers': {}
    }
    for tier in (policy.get('allowed_families') or {}):
        tier_labels = [(i, label) for i, label in scope_labelled if tier_for_label(label) == tier]
        tier_calls = [(i, label) for i, label in family_calls
                      if tiers.get(_family_policy_token(aliases.get(_family_policy_token(label), label))) == tier]
        tier_correct = [(i, label) for i, label in tier_labels
                        if _family_call_compatible(predictions[i].get('predicted_family', ''),
                                                   label,saved['parents'],policy)]
        results['tiers'][tier] = {
            'labelled_family_controls': len(tier_labels),
            'family_calls': len(tier_calls),
            'family_compatible_correct': len(tier_correct),
            'negative_control_family_calls': sum(truth[i] == 'non_adhesome' for i, _ in tier_calls),
            'precision_family_compatible': len(tier_correct) / len(tier_calls) if tier_calls else 0.0,
            'recall_family_compatible': len(tier_correct) / len(tier_labels) if tier_labels else 0.0
        }
    return results


def _require_calibrated_model(saved, purpose):
    calibration=saved.get('probability_calibration') or {}
    membership=calibration.get('membership') or {}
    family=calibration.get('family') or {}
    if (saved.get('version')=='5.8' and saved.get('feature_schema')==SCHEMA and
            membership.get('estimator') is not None and
            isinstance(family.get('eligible_classes'),list) and
            isinstance(family.get('estimators'),dict)):
        return
    version=saved.get('version','unknown')
    raise ValueError(f'{purpose} requires a v5.8 calibrated-probability model; '
        f'the supplied model is v{version} ({saved.get("feature_schema","unknown")}) and must be retrained')


def _sequence_groups(rows, sequences, c, settings):
    return sequence_cluster_groups([sequences[c.row_id(row)] for row in rows],**settings)


def _family_probability_assignment_details(assignment, probability_map, parents, all_classes, policy=None):
    """Explain direct family probabilities and hierarchy-based parent sums."""
    if not assignment:
        return {"probability": None, "basis": "unavailable",
                "eligible_descendant_classes": [], "unrepresented_descendant_classes": []}
    descendants = sorted(name for name in all_classes
                         if assignment in ancestors(name, parents))
    if not descendants and policy:
        _, _, aliases, children, _ = _family_policy_maps(policy)
        parent_token = _family_policy_token(aliases.get(_family_policy_token(assignment),assignment))
        child_tokens = children.get(parent_token,set())
        if child_tokens:
            descendants = sorted(name for name in all_classes
                if _family_policy_token(aliases.get(_family_policy_token(name),name)) in child_tokens)
            if assignment in all_classes:
                descendants.append(assignment)
            descendants=sorted(set(descendants))
            eligible=[name for name in descendants if name in probability_map]
            if eligible:
                return {"probability": float(sum(probability_map[name] for name in eligible)),
                        "basis": "sum_of_eligible_configured_broad_family_classes",
                        "eligible_descendant_classes": eligible,
                        "unrepresented_descendant_classes": [name for name in descendants if name not in probability_map]}
            return {"probability": None, "basis": "unavailable",
                    "eligible_descendant_classes": [],
                    "unrepresented_descendant_classes": descendants}
    if descendants:
        eligible = [name for name in descendants if name in probability_map]
        if not eligible:
            return {"probability": None, "basis": "unavailable",
                    "eligible_descendant_classes": [],
                    "unrepresented_descendant_classes": descendants}
        return {"probability": float(sum(probability_map[name] for name in eligible)),
                "basis": "sum_of_eligible_calibrated_descendant_classes",
                "eligible_descendant_classes": eligible,
                "unrepresented_descendant_classes": [name for name in descendants if name not in probability_map]}
    if assignment in probability_map:
        return {"probability": float(probability_map[assignment]),
                "basis": "conditional_family_class_probability",
                "eligible_descendant_classes": [assignment],
                "unrepresented_descendant_classes": []}
    return {"probability": None, "basis": "unavailable",
            "eligible_descendant_classes": [], "unrepresented_descendant_classes": []}


def infer(saved, data, c, membership_cutoff=None, family_similarity_profile=None):
    _require_calibrated_model(saved,'Prediction')
    rule_path=c._BUNDLE/'config/annotation/family_identity_rules.json'
    expected_rule_sha=saved.get('provenance',{}).get('family_identity_rules')
    if expected_rule_sha and c.sha(rule_path)!=expected_rule_sha:
        raise ValueError('Family identity-rule source differs from the one frozen into this classifier')
    family_rules_path=c._BUNDLE/'config/annotation/family_rules.json'
    expected_family_rules_sha=saved.get('provenance',{}).get('family_rules')
    if expected_family_rules_sha and c.sha(family_rules_path)!=expected_family_rules_sha:
        raise ValueError('Family alias/rule source differs from the one frozen into this classifier')
    crosswalk_path=c._BUNDLE/c._SETTINGS['family_signature_crosswalk']['path']
    current_crosswalk_sha=c.sha(crosswalk_path)
    expected_crosswalk_sha=saved.get('provenance',{}).get('adhesome_family_rule_crosswalk_sha256')
    if not expected_crosswalk_sha:
        raise ValueError('The saved calibrated model lacks its crosswalk feature-library receipt; retrain with the current classifier')
    if expected_crosswalk_sha and current_crosswalk_sha!=expected_crosswalk_sha:
        raise ValueError('Adhesome family-rule crosswalk differs from the one frozen into this classifier')
    identity_rules,_,_=load_identity_rules(str(rule_path))
    vectors = [x[1] for x in data]
    binary = saved['binary']
    binary_raw_probabilities=binary.predict_proba([membership_model_vector(vector) for vector in vectors])
    binary_classes=list(binary.classes_)
    membership_probabilities=apply_binary_sigmoid(
        saved['probability_calibration']['membership'],binary_raw_probabilities,binary_classes)
    raw_membership_votes=binary_raw_probabilities[:,binary_classes.index('adhesome_associated')]
    raw_family_probabilities=saved['family'].predict_proba(vectors)
    raw_families=list(saved['family'].classes_)
    families,family_probabilities=apply_family_sigmoids(
        saved['probability_calibration']['family'],raw_family_probabilities,raw_families)
    policy = dict(saved['policy'])
    # Assignment policy is a decision-only layer over the frozen feature model;
    # use the current versioned policy so existing compatible models can adopt it.
    family_call_policy=c._SETTINGS.get('family_assignment_policy') or saved.get('family_assignment_policy',{})
    if family_call_policy.get('mode')!='raw_vote_evidence_gated_whitelist':
        raise ValueError('Saved model or current configuration has no supported family-assignment policy')
    similarity_settings=_family_similarity_config(c,saved)
    similarity_profile_name,similarity_selection=_family_similarity_profile(similarity_settings,family_similarity_profile)
    aliases=saved.get('aliases',{})
    if membership_cutoff is not None:
        policy['membership_probability_cutoff'] = membership_cutoff
    output = []
    for i, (row, vector) in enumerate(data):
        ranking = sorted(zip(families, map(float, family_probabilities[i])), key=lambda x: (-x[1], x[0]))
        raw_ranking=sorted(zip(raw_families,map(float,raw_family_probabilities[i])),key=lambda x:(-x[1],x[0]))
        family_probability_map=dict(ranking)
        proposals=[]
        for proposal in family_identity_rule_proposals(saved,vector,identity_rules):
            scoped_proposal=_resolve_family_call('specific',proposal['family'],
                'sequence_identity_rule_proposal',saved['parents'],family_call_policy)
            if scoped_proposal:
                proposals.append({**proposal,'source_family':proposal.get('source_family') or proposal['family'],
                    'family':scoped_proposal[1],'assignment_level':scoped_proposal[0]})
        proposal_families={x['family'] for x in proposals}
        proposal_source_families={x.get('source_family') or x['family'] for x in proposals}
        calibrated_category, calibrated_assignment, calibrated_reason = decide(float(membership_probabilities[i]), ranking, row, saved['profiles'], saved['counts'], saved['parents'], policy, c,
            vector=vector,identity_rules=identity_rules,aliases=aliases,rule_families=proposal_families,
            similarity_settings=similarity_settings,similarity_selection=similarity_selection)
        precalibration_category,precalibration_assignment,precalibration_reason=decide(
            float(raw_membership_votes[i]),raw_ranking,row,saved['profiles'],saved['counts'],
            saved['parents'],policy,c,vector=vector,identity_rules=identity_rules,
            aliases=aliases,rule_families=proposal_families,
            similarity_settings=similarity_settings,similarity_selection=similarity_selection)
        raw_policy=dict(policy)
        raw_thresholds=family_call_policy.get('raw_vote_thresholds',{})
        threshold_map={
            'membership_probability_cutoff':'membership_vote_cutoff',
            'specific_family_probability_cutoff':'specific_family_vote_fraction_cutoff',
            'specific_family_probability_margin':'specific_family_vote_margin',
            'broad_parent_probability_cutoff':'broad_parent_vote_fraction_cutoff',
            'minimum_family_examples':'minimum_family_examples'}
        for destination,source in threshold_map.items():
            if source in raw_thresholds:
                raw_policy[destination]=raw_thresholds[source]
        raw_vote_category,raw_vote_assignment,raw_vote_reason=decide(
            float(raw_membership_votes[i]),raw_ranking,row,saved['profiles'],saved['counts'],
            saved['parents'],raw_policy,c,vector=vector,identity_rules=identity_rules,
            aliases=aliases,rule_families=proposal_families,
            similarity_settings=similarity_settings,similarity_selection=similarity_selection,
            assignment_policy=family_call_policy,rule_source_families=proposal_source_families)
        groups=c.evidence_groups(row)
        membership_call=('unresolved' if not any(groups.values()) else
            'adhesome_associated' if float(membership_probabilities[i])>=policy['membership_probability_cutoff'] else 'non_adhesome')
        precalibration_membership_call=('unresolved' if not any(groups.values()) else
            'adhesome_associated' if float(raw_membership_votes[i])>=policy['membership_probability_cutoff'] else 'non_adhesome')
        selected_family_call=_resolve_family_call(raw_vote_category,raw_vote_assignment,raw_vote_reason,
            saved['parents'],family_call_policy)
        suppressed_outside_scope=bool(raw_vote_assignment and raw_vote_category in {'specific','broad'} and not selected_family_call)
        if membership_call!='adhesome_associated':
            category='non_adhesome' if membership_call=='non_adhesome' else 'unresolved'
            assignment=''
            reason='calibrated_membership_below_cutoff' if membership_call=='non_adhesome' else 'no_informative_features'
            selected_from_raw=False
            selected_tier=''
            selected_source_family=''
        elif selected_family_call:
            category,assignment,reason,selected_tier,selected_source_family=selected_family_call
            selected_from_raw=True
        else:
            category='domain_only' if domain_set(row,c) else 'unresolved'
            assignment=''
            reason=('raw_vote_family_outside_allowed_scope' if suppressed_outside_scope else
                    'no_supported_whitelisted_raw_vote_family')
            selected_from_raw=False
            selected_tier=''
            selected_source_family=''
        decision_ranking=raw_ranking if selected_from_raw else ranking
        candidate_family=assignment if category=='specific' else ''
        if category=='broad' and assignment:
            _,_,policy_aliases,policy_children,_=_family_policy_maps(family_call_policy)
            assignment_token=_family_policy_token(policy_aliases.get(_family_policy_token(assignment),assignment))
            related=[name for name,_ in decision_ranking
                     if assignment in ancestors(name,saved['parents']) or
                     _family_policy_token(name) in policy_children.get(assignment_token,set()) or
                     name==assignment]
            candidate_family=related[0] if related else assignment
        if not candidate_family and selected_from_raw:
            candidate_family=selected_source_family
        candidate_level='specific'
        matched_similarity=(_family_similarity_match(row,candidate_family,candidate_level,similarity_settings,
            similarity_selection,identity_rules,aliases) if candidate_family else None)
        if not matched_similarity:
            candidate_level='broad'
            matched_similarity=_family_similarity_match(row,candidate_family,candidate_level,similarity_settings,
                similarity_selection,identity_rules,aliases)
        if category=='broad' and assignment:
            related_candidates=[name for name,_ in decision_ranking if assignment in ancestors(name,saved['parents'])]
            _,_,policy_aliases,policy_children,_=_family_policy_maps(family_call_policy)
            assignment_token=_family_policy_token(policy_aliases.get(_family_policy_token(assignment),assignment))
            related_candidates.extend(name for name,_ in decision_ranking
                if _family_policy_token(name) in policy_children.get(assignment_token,set()))
            if assignment in decision_ranking and assignment not in related_candidates:
                related_candidates.append(assignment)
            family_hits=[]
            for name in related_candidates:
                hit=_family_similarity_match(row,name,'broad',similarity_settings,similarity_selection,
                    identity_rules,aliases)
                if hit:family_hits.append(hit)
            if family_hits:
                matched_similarity=min(family_hits,key=lambda hit:(float(hit['evalue']),-float(hit['bitscore'])))
                candidate_level='broad'
        family_probability_assignment=_family_probability_assignment_details(
            assignment,family_probability_map,saved['parents'],
            sorted(set(raw_families)|set(saved.get('counts',{}))|set(saved['parents'])),family_call_policy)
        assigned_family_probability=family_probability_assignment['probability']
        raw_top_vote=raw_ranking[0][1] if raw_ranking else 0.0
        best_family_probability=ranking[0][1] if ranking else None
        probability_status=('calibrated_over_eligible_family_classes' if ranking else
            'insufficient_independent_family_support')
        prediction={'protein_id': row['protein_id'], 'classifier_call': category,
            'predicted_family': assignment, 'best_family_candidate': ranking[0][0] if ranking else '',
            'best_family_candidate_raw': raw_ranking[0][0] if raw_ranking else '',
            'best_family_probability_conditional': best_family_probability if best_family_probability is not None else '',
            'best_family_raw_vote_fraction': raw_top_vote,
            'family_raw_vote_margin': raw_ranking[0][1] - (raw_ranking[1][1] if len(raw_ranking)>1 else 0) if raw_ranking else '',
            'family_assignment_raw_vote_fraction': _family_call_vote(assignment,raw_ranking,saved['parents'],family_call_policy) if assignment else '',
            'family_assignment_score_kind': 'random_forest_raw_vote_fraction' if assignment else '',
            'family_assignment_score_is_probability': False if assignment else '',
            'family_probability_margin': ranking[0][1] - (ranking[1][1] if len(ranking)>1 else 0) if ranking else '',
            'family_probability_conditional': assigned_family_probability if assigned_family_probability is not None else '',
            'family_probability_joint': (float(membership_probabilities[i])*assigned_family_probability
                if assigned_family_probability is not None else ''),
            'family_probability_assignment_basis':family_probability_assignment['basis'],
            'family_probability_eligible_descendant_classes':';'.join(family_probability_assignment['eligible_descendant_classes']),
            'family_probability_unrepresented_descendant_classes':';'.join(family_probability_assignment['unrepresented_descendant_classes']),
            'adhesome_probability': float(membership_probabilities[i]),
            'adhesome_raw_vote_fraction': float(raw_membership_votes[i]),
            'family_probabilities_json': json.dumps(family_probability_map,sort_keys=True),
            'family_ranking_json': json.dumps(ranking),
            'family_raw_ranking_json': json.dumps(raw_ranking),
            'family_probability_calibration_status': probability_status,
            'family_probability_class_set':'eligible_calibrated_family_classes' if ranking else '',
            'family_probability_eligible_class_count':len(families),
            'family_probability_ineligible_class_count':len(set(saved['probability_calibration']['family'].get('known_classes',raw_families))-set(families)),
            'family_probability_unavailable_classes':';'.join(sorted(saved['probability_calibration']['family'].get('unavailable_classes',{}))),
            'family_probability_unavailable_reasons_json':json.dumps(saved['probability_calibration']['family'].get('unavailable_classes',{}),sort_keys=True),
            'assignment_reason': reason, 'exact_sequence_match_used': 'no',
            'score_interpretation': ('Membership probability is empirically calibrated on the labelled control distribution. Named-family calls use the previous evidence-gated Random Forest raw-vote policy, restricted to the configured core and expanded-regulatory adhesome families; raw votes are decision scores, not probabilities. Any available calibrated family probability remains a separate conditional diagnostic.'),
            'membership_call':membership_call,
            'random_forest_call':category,'random_forest_family':assignment,
            'calibrated_classifier_call':calibrated_category,
            'calibrated_predicted_family':calibrated_assignment,
            'calibrated_assignment_reason':calibrated_reason,
            'family_call_policy_version':family_call_policy.get('version',''),
            'family_call_policy_sha256':_policy_digest(family_call_policy),
            'family_call_policy_call':category,
            'family_call_policy_family':assignment,
            'family_call_policy_tier':selected_tier,
            'family_call_policy_source_family':selected_source_family,
            'family_call_policy_suppressed_outside_scope':suppressed_outside_scope,
            'precalibration_policy_call':precalibration_category,
            'precalibration_predicted_family':precalibration_assignment,
            'precalibration_assignment_reason':precalibration_reason,
            'precalibration_membership_call':precalibration_membership_call,
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
            'family_signature_features_used':True,
            'family_signature_evidence_json':row.get('_family_signature_evidence_json','[]'),
            **c.prediction_evidence_summary(row)}
        families_proposed=sorted({x['family'] for x in proposals})
        rule_rf_conflict=(len(families_proposed)==1 and bool(assignment) and
            not (_family_call_compatible(assignment,families_proposed[0],saved['parents'],family_call_policy) or
                 _family_call_compatible(families_proposed[0],assignment,saved['parents'],family_call_policy)))
        prediction.update(rule_family_candidate=families_proposed[0] if len(families_proposed)==1 else '',
            rule_candidate_families=';'.join(families_proposed),
            rule_candidate_status=('unique_rule_candidate_provisional' if len(families_proposed)==1 else
                                   'ambiguous_multiple_family_rules' if families_proposed else 'no_strong_rule_match'),
            rule_candidate_rule_ids=';'.join(sorted({x['rule_id'] for x in proposals})),
            rule_candidate_evidence_json=json.dumps(proposals,sort_keys=True),
            rule_candidate_family_untrained=(len(families_proposed)==1 and families_proposed[0] not in raw_families),
            rule_conflicts_with_random_forest=rule_rf_conflict,
            prediction_method='calibrated_membership_with_raw_vote_evidence_gated_family_whitelist',
            prediction_policy_version=family_call_policy.get('version','raw_vote_evidence_gated_whitelist_v1'))
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
        proposals.append({'rule_id':rid,'family':family,
            'source_family':str(rule.get('family','')).strip(),
            'protein_name':rule.get('protein_name',''),
            'required_groups_satisfied':int(vector.get(prefix+'required_satisfied_count',0.0)),
            'required_group_count':len(groups),'specific_cdd_pattern':pattern,
            'specific_cdd_pattern_hit_count':int(vector.get(prefix+'specific_cdd_pattern_hit_count',0.0)),
            'ordered_groups_required':bool(rule.get('ordered_groups')),
            'ordered_groups_satisfied':bool(rule.get('ordered_groups') and vector.get(prefix+'ordered_groups_satisfied',0.0)>=0.5),
            'evidence_status':'sequence-rule-supported; provisional pending independent family-level resolution'})
    return proposals


def evaluate(rows, predictions, sequences, saved, c, family_similarity_profile=None):
    truth = [c.normalized_control_class(r) for r in rows]
    family_call_policy=c._SETTINGS.get('family_assignment_policy') or saved.get('family_assignment_policy',{})
    # Score membership from its own output; family-resolution profile changes must not
    # turn an unchanged membership call into a different membership metric.
    calls = ['adhesome_associated' if p.get('membership_call')=='adhesome_associated' else 'non_adhesome'
             for p in predictions]
    binary = c.class_metrics(truth, calls, ['adhesome_associated','non_adhesome'])
    binary['positive_precision'] = c.precision_score(truth, calls, pos_label='adhesome_associated', zero_division=0)
    binary['positive_recall'] = c.recall_score(truth, calls, pos_label='adhesome_associated', zero_division=0)
    binary['positive_f1'] = c.f1_score(truth, calls, pos_label='adhesome_associated', zero_division=0)
    precalibration_calls=['adhesome_associated' if p.get('precalibration_membership_call')=='adhesome_associated' else 'non_adhesome'
                          for p in predictions]
    precalibration_binary=c.class_metrics(truth,precalibration_calls,['adhesome_associated','non_adhesome'])
    precalibration_binary.update(
        positive_precision=c.precision_score(truth,precalibration_calls,pos_label='adhesome_associated',zero_division=0),
        positive_recall=c.recall_score(truth,precalibration_calls,pos_label='adhesome_associated',zero_division=0),
        positive_f1=c.f1_score(truth,precalibration_calls,pos_label='adhesome_associated',zero_division=0))
    binary['precalibration_policy_comparison']=precalibration_binary
    membership_truth=[int(value=='adhesome_associated') for value in truth]
    binary['probability_calibration']=probability_metrics_binary(membership_truth,
        [float(p['adhesome_probability']) for p in predictions],
        [float(p['adhesome_raw_vote_fraction']) for p in predictions],
        c._SETTINGS['probability_calibration']['reliability_bins'])
    labelled = [(i, trusted_family(r, saved['aliases'], saved['parents'])) for i,r in enumerate(rows) if truth[i]=='adhesome_associated']
    labelled = [(i,label) for i,label in labelled if label]
    expected = [label for _,label in labelled]
    assigned = [predictions[i]['predicted_family'] or '__abstained__' for i,_ in labelled]
    precalibration_assigned=[predictions[i].get('precalibration_predicted_family') or '__abstained__'
                             for i,_ in labelled]
    calibrated_assigned=[predictions[i].get('calibrated_predicted_family') or '__abstained__'
                         for i,_ in labelled]
    shared_family_labels=sorted(set(expected)|set(assigned)|set(precalibration_assigned)|set(calibrated_assigned))
    family = c.class_metrics(expected, assigned, shared_family_labels)
    family['precalibration_policy_comparison']=c.class_metrics(
        expected,precalibration_assigned,shared_family_labels)
    family['calibrated_policy_comparison']=c.class_metrics(
        expected,calibrated_assigned,shared_family_labels)
    family['labelled_positives'] = len(labelled)
    family['unmapped_positives'] = truth.count('adhesome_associated') - len(labelled)
    family['unseen_label_records'] = sum(label not in saved['counts'] for label in expected)
    family['correct_accepted'] = sum(a==b for a,b in zip(expected,assigned))
    family['accepted'] = sum(x!='__abstained__' for x in assigned)
    family['accepted_precision'] = family['correct_accepted'] / family['accepted'] if family['accepted'] else 0.0
    family['strict_recall'] = family['correct_accepted'] / len(labelled) if labelled else 0.0
    family['top_candidate_correct'] = sum(predictions[i]['best_family_candidate_raw']==label for i,label in labelled)
    family['calibrated_top_candidate_correct'] = sum(predictions[i]['best_family_candidate']==label for i,label in labelled)
    family['top_three_contains_truth'] = sum(label in [x[0] for x in json.loads(predictions[i]['family_raw_ranking_json'])[:3]] for i,label in labelled)
    _,_,_,broad_children,_=_family_policy_maps(family_call_policy)
    family['correct_broad_or_specific'] = sum(
        predictions[i]['predicted_family'] in [label]+ancestors(label,saved['parents']) or
        _family_policy_token(label) in broad_children.get(_family_policy_token(predictions[i]['predicted_family']),set())
        for i,label in labelled)
    eligible_family_classes=list(saved['probability_calibration']['family'].get('eligible_classes',[]))
    family_probability_rows=[]
    for i,label in labelled:
        if label not in eligible_family_classes:
            continue
        calibrated_map=json.loads(predictions[i].get('family_probabilities_json') or '{}')
        raw_map=dict(json.loads(predictions[i].get('family_raw_ranking_json') or '[]'))
        calibrated_values=[float(calibrated_map.get(name,0.0)) for name in eligible_family_classes]
        raw_values=[float(raw_map.get(name,0.0)) for name in eligible_family_classes]
        raw_total=sum(raw_values)
        if raw_total>0:
            raw_values=[value/raw_total for value in raw_values]
        family_probability_rows.append((label,calibrated_values,raw_values))
    family['probability_calibration']={
        'eligible_family_classes':eligible_family_classes,
        'family_labelled_records':len(labelled),
        'scored_records':len(family_probability_rows),
        'unavailable_true_family_records':len(labelled)-len(family_probability_rows),
        'coverage':len(family_probability_rows)/len(labelled) if labelled else 0.0,
        'calibration_availability_by_family':{
            label:{'calibrated':label in eligible_family_classes,
                   **saved['probability_calibration']['family'].get('class_support',{}).get(label,{}),
                   'family_training_positive_records':saved['probability_calibration']['family'].get('family_training_support',{}).get(label,{}).get('positive_records',0),
                   'family_training_positive_independent_groups':saved['probability_calibration']['family'].get('family_training_support',{}).get(label,{}).get('positive_independent_sequence_groups',0),
                   'unavailable_reason':saved['probability_calibration']['family'].get('unavailable_classes',{}).get(label,'')}
            for label in saved['probability_calibration']['family'].get('known_classes',saved['probability_calibration']['family'].get('classes',[]))},
        'metrics':probability_metrics_multiclass(
            [row[0] for row in family_probability_rows],
            [row[1] for row in family_probability_rows],eligible_family_classes,
            [row[2] for row in family_probability_rows],
            c._SETTINGS['probability_calibration']['reliability_bins'])}
    family_calibration=saved['probability_calibration']['family']
    min_positive=int(family_calibration.get('minimum_positive_groups',0) or 0)
    min_negative=int(family_calibration.get('minimum_negative_groups',0) or 0)
    min_eligible=int(family_calibration.get('minimum_eligible_classes',0) or 0)
    support=family_calibration.get('class_support',{})
    classes_meeting_floors=[]
    for name,counts in support.items():
        positive_groups=int(counts.get('positive_independent_groups',counts.get('independent_groups',0)) or 0)
        negative_groups=int(counts.get('negative_independent_groups',0) or 0)
        if positive_groups>=min_positive and negative_groups>=min_negative:
            classes_meeting_floors.append(name)
    family['probability_calibration'].update(
        classes_meeting_support_floors=sorted(classes_meeting_floors),
        classes_meeting_support_floor_count=len(classes_meeting_floors),
        minimum_positive_groups=min_positive,
        minimum_negative_groups=min_negative,
        minimum_eligible_classes=min_eligible,
        family_call_availability_reason=(
            f'{len(classes_meeting_floors)} {"class" if len(classes_meeting_floors)==1 else "classes"} met the per-class independent-group floors; '
            f'{min_eligible} eligible classes are required for family probability calibration'
            if len(eligible_family_classes)<min_eligible else 'Calibrated family classes are available'))
    family['per_family'] = {label: {'records': expected.count(label),
                                  'correct': sum(a==label and b==label for a,b in zip(expected,assigned)),
                                  'missed_or_wrong': sum(a==label and b!=label for a,b in zip(expected,assigned)),
                                  'calibrated_correct':sum(a==label and b==label for a,b in zip(expected,calibrated_assigned)),
                                  'calibrated_missed_or_wrong':sum(a==label and b!=label for a,b in zip(expected,calibrated_assigned)),
                                  'precalibration_correct':sum(a==label and b==label for a,b in zip(expected,precalibration_assigned)),
                                  'precalibration_missed_or_wrong':sum(a==label and b!=label for a,b in zip(expected,precalibration_assigned))}
                         for label in sorted(set(expected))}
    family['selected_family_call_policy']=_family_policy_metrics(
        rows,predictions,truth,[trusted_family(row,saved['aliases'],saved['parents']) for row in rows],
        saved,c,family_call_policy)
    family['family_compatible_accepted_precision']=family['selected_family_call_policy']['scope_precision_family_compatible']
    family['family_compatible_in_scope_recall']=family['selected_family_call_policy']['scope_recall_family_compatible']
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
            'family_signature_crosswalk':saved.get('family_signature_crosswalk',{}),
            'family_similarity_policy':{'profile':profile_name,'search_max_evalue':similarity_settings['search_max_evalue'],
                'thresholds':profile_values,'membership_model_uses_similarity_extension':False},
            'class_counts': dict(Counter(p['classifier_call'] for p in predictions)),
            'probability_calibration_policy':{
                'method':saved['probability_calibration'].get('method'),
                'grouping':saved['probability_calibration'].get('grouping'),
                'validation_used_for_calibration':saved['probability_calibration'].get('validation_used_for_calibration'),
                'membership':{key:value for key,value in saved['probability_calibration']['membership'].items()
                              if key!='estimator'},
                'family':{key:value for key,value in saved['probability_calibration']['family'].items()
                          if key!='estimators'}},
            'calibration_overlap_strata': strata,
            'family_assignment_policy':family_call_policy,
            'note': 'Validation was evaluation-only and did not fit calibrators or choose thresholds. Membership probabilities are conditional on the labelled control distribution. Family probabilities are conditional on adhesome association and are reported only for families meeting independent-group calibration support. Named family calls use the separately versioned, evidence-gated raw-vote shortlist policy; raw votes are not probabilities. Exact sequence overlap is reporting only; evolutionary homology overlap may remain.'}


def compare_validation(rows, baseline, current, saved, c):
    """Same-cohort pre/post comparison; baseline calls are outcomes, never features."""
    truth=[c.normalized_control_class(r) for r in rows]
    expected=[trusted_family(r,saved['aliases'],saved['parents']) for r in rows]
    def summarize(predictions):
        calls=[]
        for prediction in predictions:
            membership=prediction.get('membership_call')
            if membership is not None:
                calls.append('adhesome_associated' if membership=='adhesome_associated' else 'non_adhesome')
            else:
                calls.append('adhesome_associated' if prediction.get('classifier_call') in ('specific','broad','domain_only') else 'non_adhesome')
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
        ('Family-compatible correct calls, including configured broad RhoGEF matches',family['selected_family_call_policy']['family_compatible_correct_in_scope']),
        ('Calibrated family classes available',len(family['probability_calibration']['eligible_family_classes'])),
        ('Classes meeting family calibration support floors',family['probability_calibration']['classes_meeting_support_floor_count']),
        ('Calibrated family-call availability',family['probability_calibration']['family_call_availability_reason']),
        ('Strict family recall (abstentions count as missed)',f"{family['strict_recall']:.1%}"),
        ('Top raw-vote-ranked candidate correct (before acceptance)',family['top_candidate_correct']),
        ('Top calibrated-family candidate correct (before acceptance)',family['calibrated_top_candidate_correct']),
        ('Truth within top three candidates',family['top_three_contains_truth']),
        ('Positive controls without family labels',family['unmapped_positives'])]
    def table(items):
        return '<table>'+''.join('<tr><th>'+html.escape(str(k))+'</th><td>'+html.escape(str(v))+'</td></tr>' for k,v in items)+'</table>'
    def reliability_svg(title, calibrated, raw=None):
        width,height,left,top,plot_w,plot_h=480,330,55,28,390,240
        def point_string(summary):
            values=[]
            for item in (summary or {}).get('bins',[]):
                if item.get('count',0) and item.get('mean_probability') is not None:
                    x=left+float(item['mean_probability'])*plot_w
                    y=top+(1.0-float(item['observed_fraction']))*plot_h
                    values.append(f'{x:.1f},{y:.1f}')
            return ' '.join(values)
        diagonal=f'{left},{top+plot_h} {left+plot_w},{top}'
        svg=(f'<svg role="img" aria-label="{html.escape(title)}" viewBox="0 0 {width} {height}" '
             f'style="max-width:480px;width:100%;border:1px solid #ddd"><text x="{left}" y="18">'
             f'{html.escape(title)}</text><line x1="{left}" y1="{top+plot_h}" x2="{left+plot_w}" '
             f'y2="{top}" stroke="#888" stroke-dasharray="5,5"/><line x1="{left}" y1="{top}" '
             f'x2="{left}" y2="{top+plot_h}" stroke="#333"/><line x1="{left}" y1="{top+plot_h}" '
             f'x2="{left+plot_w}" y2="{top+plot_h}" stroke="#333"/>')
        raw_points=point_string(raw)
        cal_points=point_string(calibrated)
        if raw_points:
            svg+=f'<polyline points="{raw_points}" fill="none" stroke="#b55" stroke-width="2"/>'
        if cal_points:
            svg+=f'<polyline points="{cal_points}" fill="none" stroke="#1769aa" stroke-width="2"/>'
        svg+=(f'<text x="{left+90}" y="{top+plot_h+35}" font-size="12">Mean predicted probability</text>'
              f'<text x="{left+plot_w+5}" y="{top+plot_h+4}" font-size="11">1</text>'
              f'<text x="{left-15}" y="{top+5}" font-size="11">1</text>'
              f'<text x="{left+plot_w-115}" y="{height-12}" font-size="12" fill="#1769aa">Calibrated</text>')
        if raw_points:
            svg+=f'<text x="{left+plot_w-215}" y="{height-12}" font-size="12" fill="#b55">Raw votes</text>'
        return svg+'</svg>'
    document = '<!doctype html><meta charset="utf-8"><title>Family classifier validation</title><style>body{font:16px system-ui;max-width:1000px;margin:40px auto}td,th{padding:8px;border:1px solid #ccc;text-align:left}table{border-collapse:collapse}pre{white-space:pre-wrap}</style><h1>Independent panel evaluation</h1>'
    document += '<p>'+html.escape(report['note'])+'</p>'+table(rows)
    prior_membership=membership.get('precalibration_policy_comparison',{})
    if prior_membership:
        document+='<h2>Previous versus calibrated-policy calls</h2><p>Both policies use the same labelled records and the existing numeric cutoffs. The previous-policy column applies those cutoffs to raw forest votes; the calibrated-policy column applies them to calibrated probabilities. This is a fixed comparison, not validation tuning.</p>'
        document+='<table><tr><th>Membership metric</th><th>Previous raw-vote policy</th><th>Calibrated policy</th></tr>'
        for key,label in (('positive_precision','Precision'),('positive_recall','Recall'),('positive_f1','F1')):
            document+=f'<tr><th>{html.escape(label)}</th><td>{prior_membership.get(key,0.0):.3f}</td><td>{membership.get(key,0.0):.3f}</td></tr>'
        document+='</table><h3>Membership confusion matrices</h3><p>Previous raw-vote policy</p><pre>'+html.escape(json.dumps(prior_membership.get('confusion_matrix',[])))+'</pre><p>Calibrated policy</p><pre>'+html.escape(json.dumps(membership.get('confusion_matrix',[])))+'</pre>'
        prior_family=family.get('precalibration_policy_comparison',{})
        if prior_family:
            document+='<h3>Family assignment comparison on family-labelled positives</h3><table><tr><th>Metric</th><th>Previous raw-vote policy</th><th>Calibrated policy</th></tr>'
            for key,label in (('macro_precision','Macro precision'),('macro_recall','Macro recall'),('macro_f1','Macro F1')):
                document+=f'<tr><th>{html.escape(label)}</th><td>{prior_family.get(key,0.0):.3f}</td><td>{family.get(key,0.0):.3f}</td></tr>'
            document+='</table><p>Both family comparisons use the same labelled positives; abstentions remain explicit in the confusion matrices.</p><h4>Previous raw-vote family confusion matrix</h4><pre>'+html.escape(json.dumps(prior_family.get('confusion_matrix',[])))+'</pre><h4>Calibrated family confusion matrix</h4><pre>'+html.escape(json.dumps(family.get('confusion_matrix',[])))+'</pre>'
    membership_calibration=membership.get('probability_calibration',{})
    if membership_calibration:
        membership_rows=[('Scored records',membership_calibration.get('records',0)),
            ('Calibrated Brier score',membership_calibration.get('brier_score')),
            ('Raw-vote Brier score',membership_calibration.get('raw',{}).get('brier_score')),
            ('Calibrated log loss',membership_calibration.get('log_loss')),
            ('Raw-vote log loss',membership_calibration.get('raw',{}).get('log_loss')),
            ('Calibrated ECE',membership_calibration.get('reliability',{}).get('expected_calibration_error')),
            ('Raw-vote ECE',membership_calibration.get('raw',{}).get('reliability',{}).get('expected_calibration_error'))]
        document+='<h2>Membership probability calibration</h2>'+table(membership_rows)
        document+=reliability_svg('Membership reliability',membership_calibration.get('reliability'),
            membership_calibration.get('raw',{}).get('reliability'))
    family_calibration=family.get('probability_calibration',{})
    family_metrics=family_calibration.get('metrics',{})
    if family_calibration:
        family_rows=[('Family-labelled records',family_calibration.get('family_labelled_records',0)),
            ('Scored records',family_calibration.get('scored_records',0)),
            ('Unavailable true-family records',family_calibration.get('unavailable_true_family_records',0)),
            ('Calibrated family coverage',f"{family_calibration.get('coverage',0.0):.1%}"),
            ('Calibrated multiclass Brier score',family_metrics.get('brier_score')),
            ('Raw-vote multiclass Brier score',family_metrics.get('raw',{}).get('brier_score')),
            ('Calibrated multiclass log loss',family_metrics.get('log_loss')),
            ('Raw-vote multiclass log loss',family_metrics.get('raw',{}).get('log_loss')),
            ('Calibrated top-label ECE',family_metrics.get('top_label_reliability',{}).get('expected_calibration_error')),
            ('Raw-vote top-label ECE',family_metrics.get('raw',{}).get('top_label_reliability',{}).get('expected_calibration_error'))]
        document+='<h2>Family probability calibration</h2>'+table(family_rows)
        document+=reliability_svg('Top family reliability',family_metrics.get('top_label_reliability'),
            family_metrics.get('raw',{}).get('top_label_reliability'))
        availability=family_calibration.get('calibration_availability_by_family',{})
        if availability:
            document+='<h3>Per-family calibration support</h3><table><tr><th>Family</th><th>Available</th><th>Family-training positive records</th><th>Family-training positive groups</th><th>Calibrator positive groups</th><th>Calibrator negative groups</th><th>Reason unavailable</th></tr>'
            for label,details in sorted(availability.items()):
                document+=('<tr><td>'+html.escape(label)+'</td><td>'+str(bool(details.get('calibrated')))+
                    '</td><td>'+str(details.get('family_training_positive_records',0))+
                    '</td><td>'+str(details.get('family_training_positive_independent_groups',0))+
                    '</td><td>'+str(details.get('independent_groups',''))+'</td><td>'+str(details.get('negative_independent_groups',''))+
                    '</td><td>'+html.escape(str(details.get('unavailable_reason','')) )+'</td></tr>')
            document+='</table>'
        for label,calibration in family_metrics.get('per_class_reliability',{}).items():
            document+=reliability_svg('Family reliability: '+label,calibration)
    similarity=report.get('family_similarity_policy')
    if similarity:
        policy_rows=[('Family similarity profile',similarity['profile']),
            ('CDD search ceiling',similarity['search_max_evalue']),
            ('Specific-family threshold',json.dumps(similarity['thresholds']['specific'],sort_keys=True)),
            ('Broad-family threshold',json.dumps(similarity['thresholds']['broad'],sort_keys=True)),
            ('Used to alter membership model',similarity['membership_model_uses_similarity_extension'])]
        document += '<h2>Family-resolution similarity policy</h2><p>CDD-fragment sequence matches support family resolution by E-value and coverage. The membership model continues to use its existing feature set and decision cutoff.</p>'+table(policy_rows)
    crosswalk=report.get('family_signature_crosswalk')
    if crosswalk:
        crosswalk_rows=[('Crosswalk signatures',crosswalk.get('signature_combinations',0)),
            ('Mapped families',crosswalk.get('mapped_families',0)),
            ('Mapped family identity rules',crosswalk.get('mapped_identity_rules',0)),
            ('Alias proposals ignored',crosswalk.get('alias_proposal_rows_ignored',0)),
            ('Crosswalk SHA-256',crosswalk.get('sha256',''))]
        document += '<h2>Family-signature crosswalk feature source</h2><p>CDD hits contribute to family features only in combination with configured family identity-rule architecture predicates.</p>'+table(crosswalk_rows)
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
    selected_policy=family['selected_family_call_policy']
    document += '<h2>Selected family-call policy</h2><p>'+html.escape(selected_policy['score_interpretation'])+'</p>'
    policy_rows=[('Policy version',selected_policy['version']),
        ('Family-labelled controls within the allowed scope',selected_policy['in_scope_labelled_family_controls']),
        ('Selected family calls',selected_policy['selected_family_calls']),
        ('Exact correct family calls in scope',selected_policy['exact_correct_in_scope']),
        ('Compatible correct calls including broad family assignments',selected_policy['family_compatible_correct_in_scope']),
        ('Incorrect or unlabelled family calls',selected_policy['incorrect_or_unlabelled_family_calls']),
        ('Negative-control family calls',selected_policy['negative_control_family_calls']),
        ('Scope precision, counting configured broad matches',f"{selected_policy['scope_precision_family_compatible']:.1%}"),
        ('Scope recall, counting configured broad matches',f"{selected_policy['scope_recall_family_compatible']:.1%}"),
        ('Exact family recall across all labelled positives',f"{selected_policy['overall_labelled_exact_recall']:.1%}")]
    document += table(policy_rows)
    policy_tier_rows=[('Scope','Labelled controls','Calls','Family-compatible correct','Precision','Recall','Negative-control calls')]
    for tier,values in selected_policy['tiers'].items():
        policy_tier_rows.append((tier,values['labelled_family_controls'],values['family_calls'],
            values['family_compatible_correct'],f"{values['precision_family_compatible']:.1%}",
            f"{values['recall_family_compatible']:.1%}",values['negative_control_family_calls']))
    document += '<table>'+''.join('<tr>'+''.join(
        '<th>'+html.escape(str(value))+'</th>' if row_index==0 else '<td>'+html.escape(str(value))+'</td>'
        for value in row)+'</tr>' for row_index,row in enumerate(policy_tier_rows))+'</table>'
    document += '<h2>Per-family exact assignment counts</h2><table><tr><th>Family</th><th>Selected shortlist raw-vote policy correct / labelled</th><th>Calibrated policy correct / labelled</th><th>Previous unrestricted raw-vote policy correct / labelled</th></tr>'
    for label,values in family['per_family'].items():
        document+=('<tr><td>'+html.escape(label)+'</td><td>'+str(values['correct'])+' / '+str(values['records'])+
                   '</td><td>'+str(values['calibrated_correct'])+' / '+str(values['records'])+
                   '</td><td>'+str(values['precalibration_correct'])+' / '+str(values['records'])+'</td></tr>')
    document+='</table>'
    document += '<details><summary>Full metrics and overlap strata</summary><pre>'+html.escape(json.dumps(report,indent=2))+'</pre></details>'
    (directory/'internal_validation_report.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    (directory/'internal_validation_report.html').write_text(document,encoding='utf-8')


def train_family_gate(args, c):
    print('[classifier] Loading family rules, calibration labels, and verified Stage 1 evidence',flush=True)
    aliases, parents, motifs = load_resources(c)
    similarity_settings=c._SETTINGS['family_gate']['family_similarity']
    similarity_profile_name,similarity_selection=_family_similarity_profile(similarity_settings)
    _,crosswalk_metadata=_family_signature_library(c,aliases)
    print(f"[classifier] Family-signature crosswalk ready: {crosswalk_metadata['signature_combinations']} CDD/architecture combinations across {crosswalk_metadata['mapped_families']} families",flush=True)
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
    probability_settings=c._SETTINGS['probability_calibration']
    group_ids,grouping_metadata=_sequence_groups(rows,sequences,c,probability_settings['sequence_clustering'])
    print(f"[classifier] OOF leakage grouping: {grouping_metadata['method']} "
          f"(98% identity, 90% coverage of both sequences); "
          f"fallback={grouping_metadata.get('fallback_reason') or 'not used'}",flush=True)
    membership_matrix=[membership_model_vector(v) for _,v in data]
    print('[classifier] Generating group-disjoint OOF membership scores for probability calibration',flush=True)
    try:
        membership_classes,membership_oof,membership_folds,membership_fold_summary=\
            cross_fitted_probabilities(lambda:estimator(1729),membership_matrix,labels,group_ids,
                probability_settings['max_splits'],probability_settings['random_state'])
        membership_calibration=fit_binary_sigmoid(membership_oof,labels,membership_classes,group_ids,
            minimum_groups_per_class=probability_settings['membership_min_independent_groups_per_class'])
    except CalibrationUnavailable as exc:
        if grouping_metadata.get('method')!='cdhit_grouped_out_of_fold':
            raise ValueError(f'Membership probabilities cannot be calibrated from this control set: {exc}') from exc
        exact_settings=dict(probability_settings['sequence_clustering'],executable='')
        exact_group_ids,exact_metadata=_sequence_groups(rows,sequences,c,exact_settings)
        try:
            membership_classes,membership_oof,membership_folds,membership_fold_summary=\
                cross_fitted_probabilities(lambda:estimator(1729),membership_matrix,labels,exact_group_ids,
                    probability_settings['max_splits'],probability_settings['random_state'])
            membership_calibration=fit_binary_sigmoid(membership_oof,labels,membership_classes,exact_group_ids,
                minimum_groups_per_class=probability_settings['membership_min_independent_groups_per_class'])
        except CalibrationUnavailable as exact_exc:
            raise ValueError(f'Membership probabilities cannot be calibrated from this control set; '
                f'CD-HIT grouping failed ({exc}) and exact-sequence fallback failed ({exact_exc})') from exact_exc
        group_ids=exact_group_ids
        grouping_metadata={**exact_metadata,
            'cdhit_grouping_rejected_reason':str(exc),
            'fallback_reason':f'CD-HIT groups did not support membership calibration; exact-sequence grouping used: {exc}'}
        print(f"[classifier] OOF grouping fallback: {grouping_metadata['fallback_reason']}",flush=True)
    print('[classifier] Fitting final adhesome membership Random Forest',flush=True)
    binary = estimator(1729).fit(membership_matrix, labels)
    indices = [(i,trusted_family(row,aliases,parents)) for i,row in enumerate(rows) if labels[i]=='adhesome_associated']
    indices = [(i,label) for i,label in indices if label and any(c.evidence_groups(data[i][0]).values())]
    if len({label for _,label in indices}) < 2:
        raise ValueError('At least two trusted, feature-supported calibration families are required')
    counts = Counter(label for _,label in indices)
    model_indices=[(i,label) for i,label in indices if counts[label]>=2]
    print(f'[classifier] Fitting family model from {len(model_indices)} labelled, feature-supported controls',flush=True)
    family_matrix=[data[i][1] for i,_ in model_indices]
    family_labels=[label for _,label in model_indices]
    family_groups=[group_ids[i] for i,_ in model_indices]
    family = estimator(2718).fit(family_matrix,family_labels)
    family_classes=list(family.classes_)
    family_fold_assignments=[]
    print('[classifier] Generating group-disjoint OOF family scores for per-family calibration',flush=True)
    try:
        oof_classes,family_oof,family_folds,family_fold_summary=\
            cross_fitted_probabilities(lambda:estimator(2718),family_matrix,family_labels,family_groups,
                probability_settings['max_splits'],probability_settings['random_state'])
        family_calibration=fit_family_sigmoids(family_oof,family_labels,oof_classes,family_groups,
            minimum_positive_groups=probability_settings['family_min_positive_groups_per_class'],
            minimum_negative_groups=probability_settings['family_min_negative_groups_per_class'],
            minimum_eligible_classes=probability_settings['family_min_calibrated_classes'])
        for local_index,(global_index,_) in enumerate(model_indices):
            family_fold_assignments.append((global_index,int(family_folds[local_index]),family_groups[local_index]))
    except CalibrationUnavailable as exc:
        family_fold_summary=[]
        family_group_support={}
        for label in family_classes:
            positive_groups_for_label={family_groups[index] for index,value in enumerate(family_labels) if value==label}
            negative_groups_for_label={family_groups[index] for index,value in enumerate(family_labels) if value!=label}
            family_group_support[label]={
                'records':sum(value==label for value in family_labels),
                'independent_groups':len(positive_groups_for_label),
                'negative_records':sum(value!=label for value in family_labels),
                'negative_independent_groups':len(negative_groups_for_label)}
        family_calibration={'method':'unavailable_insufficient_grouped_oof_support','classes':family_classes,
            'eligible_classes':[],'estimators':{},'class_support':family_group_support,
            'unavailable_classes':{label:str(exc) for label in family_classes},
            'minimum_positive_groups':int(probability_settings['family_min_positive_groups_per_class']),
            'minimum_negative_groups':int(probability_settings['family_min_negative_groups_per_class']),
            'minimum_eligible_classes':int(probability_settings['family_min_calibrated_classes']),
            'parameters':{'estimator':'sklearn.linear_model.LogisticRegression','C':1e6,
                'solver':'lbfgs','max_iter':2000,
                'input':'per-family logit-transformed random-forest probability',
                'mapping':'one-vs-rest sigmoid, then normalize across eligible families'},
            'oof_records':0,'unavailable_reason':str(exc)}
        print(f'[classifier] Family probability calibration unavailable: {exc}',flush=True)
    family_calibration['known_classes']=sorted(counts)
    family_calibration['family_training_support']={
        label:{'positive_records':sum(family==label for _,family in indices),
               'positive_independent_sequence_groups':len({group_ids[index] for index,family in indices if family==label})}
        for label in sorted(counts)}
    for label in sorted(counts):
        if label not in family_classes:
            family_calibration.setdefault('unavailable_classes',{})[label]=(
                f"not included in the family Random Forest: only {counts[label]} feature-supported positive controls; requires at least 2")
    print(f"[classifier] Calibrated family classes available: {len(family_calibration['eligible_classes'])}/{len(family_classes)}",flush=True)
    membership_fold_assignments=[(index,int(fold),group_ids[index]) for index,fold in enumerate(membership_folds)]
    profiles = {label: sorted(set().union(*(domain_set(data[i][0],c) for i,f in indices if f==label))) for label in counts}
    calibration_provenance={'method':'sigmoid','grouping':probability_settings['grouping'],
        'grouping_details':grouping_metadata,
        'random_state':probability_settings['random_state'],'max_splits':probability_settings['max_splits'],
        'gate_threshold_policy':c._SETTINGS['probability_gate_threshold_policy'],
        'membership':membership_calibration,'family':family_calibration,
        'membership_fold_summary':membership_fold_summary,'family_fold_summary':family_fold_summary,
        'membership_calibration_records':len(labels),'family_calibration_records':len(family_labels),
        'validation_used_for_calibration':False}
    saved = dict(version='5.8',feature_schema=SCHEMA,binary=binary,family=family,config_vocab=sorted(vocab),
        aliases=aliases,parents=parents,motifs=motifs,counts=dict(counts),profiles=profiles,
        family_model_counts=dict(Counter(label for _,label in model_indices)),
        single_example_family_profiles=sorted(label for label,n in counts.items() if n==1),
        policy=c._SETTINGS['family_gate'],family_similarity_settings=similarity_settings,
        family_assignment_policy=c._SETTINGS['family_assignment_policy'],
        family_signature_crosswalk=crosswalk_metadata,
        probability_calibration=calibration_provenance,
        calibration_sequence_hashes={sequence_digest(sequences[c.row_id(r)]) for r in rows},
        provenance={**provenance,**receipts,**scan_hashes,'labels':c.sha(args.calibration/'reports/predictions.tsv'),
            'cdd_domain_fasta':c.sha(library/'cdd_domain_sequences.faa'),
            'cdd_domain_metadata':c.sha(library/'cdd_domain_metadata.tsv'),
            'family_identity_rules':c.sha(c._BUNDLE/'config/annotation/family_identity_rules.json'),
            'family_rules':c.sha(c._BUNDLE/'config/annotation/family_rules.json'),
            'posthoc_sequence_rules_sha256':scan_hashes['posthoc_sequence_rules_sha256'],
            'family_assignment_policy_sha256':_policy_digest(c._SETTINGS['family_assignment_policy'])})
    summary = {'classifier_version':'5.8','feature_schema':SCHEMA,'calibration_records':len(rows),
        'family_training_records':len(model_indices),'family_profile_records':len(indices),
        'family_model_counts':dict(Counter(label for _,label in model_indices)),
        'single_example_family_profiles':sorted(label for label,n in counts.items() if n==1),
        'family_profile_counts':dict(counts),
        'unmapped_or_ambiguous_positive_labels':sum(x=='adhesome_associated' for x in labels)-len(indices),
        'family_labels_excluded_from_multiclass_fit':sum(1 for n in counts.values() if n<2),
        'policy':saved['policy'],'family_assignment_policy':saved['family_assignment_policy'],
        'family_assignment_policy_sha256':_policy_digest(saved['family_assignment_policy']),
        'validation_used_for_fitting':False,
        'feature_contract':'Identical FASTA motif scans, verified Pfam hits, calibration CDD local similarities, domain order/copy/quality, directional homology, sequence-only boolean/token features parsed from post-hoc rule conditions, and family-rule crosswalk CDD-plus-architecture combinations; no post-hoc proposals, audit scores, crosswalk family/class text, or hypothesis labels.',
        'identity_rule_policy':'Only measurable sequence fields from family_identity_rules.json are added: required domain-group counts/satisfaction, specific CDD annotation pattern matches, ordered groups, maximum domain counts, repeat indices and protein length ranges. Topology, orthology, localization and signal peptide are excluded from this feature block.',
        'posthoc_rule_feature_policy':'Parse only the boolean conditions in specific() and broad() from schisto_pfam_cdd_reclassify.py. Include rule-condition hits, individual Pfam/CDD signature presence and counts, and concurrent-signal counts. Never read rule return values, proposed family labels, tiers, orthogroup arguments, post-hoc assignment tables, or validation outcomes as features.',
        'family_signature_crosswalk':crosswalk_metadata,
        'family_signature_crosswalk_policy':'For each reviewed crosswalk CDD accession, combine its calibration-bank HSP evidence with the linked family_identity_rules.json required architecture, order, pattern, copy and length predicates. Complete and partial combination features are learned by the family classifier. Use the same frozen feature builder on calibration, validation and target data. Proposed family/class text, same-protein family context, calibration counts and alias proposals are excluded.',
        'cdd_policy':'All SHA-256-verified CDD fragments from labelled calibration-positive controls are retained, including fragments without a mapped family label. Validation-origin fragments are excluded. Calibration queries exclude same-protein references. Query and subject coverage use inclusive BLAST HSP coordinate spans divided by full sequence lengths, capped at 1.0.',
        'optional_returned_policy':'Returned CDD hits are not required; candidate sequences are searched against the frozen calibration-only CDD fragment bank.'}
    summary['probability_calibration']={
        'method':'sigmoid; fitted from group-disjoint out-of-fold calibration predictions',
        'grouping':probability_settings['grouping'],
        'grouping_details':grouping_metadata,
        'membership':{'calibrated':True,'method':membership_calibration['method'],
            'parameters':membership_calibration['parameters'],
            'class_group_counts':membership_calibration['class_group_counts'],
            'minimum_groups_per_class':membership_calibration['minimum_groups_per_class'],
            'folds':membership_fold_summary},
        'family':{'method':family_calibration['method'],
            'parameters':family_calibration.get('parameters',{}),
            'eligible_classes':family_calibration['eligible_classes'],
            'unavailable_classes':family_calibration['unavailable_classes'],
            'class_support':family_calibration['class_support'],
            'minimum_positive_groups':family_calibration.get('minimum_positive_groups'),
            'minimum_negative_groups':family_calibration.get('minimum_negative_groups'),
            'minimum_eligible_classes':family_calibration.get('minimum_eligible_classes'),
            'folds':family_fold_summary},
        'validation_used_for_calibration':False,
        'gate_threshold_policy':c._SETTINGS['probability_gate_threshold_policy'],
        'interpretation':'Membership probability is conditional on the labelled control distribution. Family probability is conditional on adhesome association. Joint family probability is their product and is not a universal biological probability.'}
    summary['family_similarity_policy']={'selected_profile':similarity_profile_name,
        'search_max_evalue':similarity_settings['search_max_evalue'],
        'legacy_model_feature_max_evalue':similarity_settings['legacy_model_feature_max_evalue'],
        'profile_thresholds':similarity_selection,
        'membership_model_uses_family_similarity_extension':False,
        'membership_model_uses_family_signature_crosswalk':False,
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
    fold_rows=[]
    for index,fold,group_id in membership_fold_assignments:
        fold_rows.append({'model':'membership','protein_id':c.row_id(rows[index]),
                          'sequence_group_id':group_id,'heldout_fold':fold})
    for index,fold,group_id in family_fold_assignments:
        fold_rows.append({'model':'family','protein_id':c.row_id(rows[index]),
                          'sequence_group_id':group_id,'heldout_fold':fold})
    fold_path=args.out/'probability_calibration_folds.tsv'
    c.write_tsv(fold_path,fold_rows,fields=['model','protein_id','sequence_group_id','heldout_fold'])
    saved['provenance']['probability_calibration_folds']=c.sha(fold_path)
    saved['provenance']['frozen_cdd_reference_fasta']=c.sha(args.out/'cdd_calibration_references.faa')
    saved['provenance']['frozen_cdd_reference_metadata']=c.sha(args.out/'cdd_calibration_references.json')
    saved['provenance']['frozen_cdd_reference_manifest']=c.sha(args.out/'cdd_calibration_references.manifest.json')
    summary['cdd_reference_bank']=json.loads(frozen_cdd_manifest.decode('utf-8'))
    summary['probability_calibration']['fold_assignments_file']='probability_calibration_folds.tsv'
    summary['probability_calibration']['fold_assignments_sha256']=saved['provenance']['probability_calibration_folds']
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
    _require_calibrated_model(saved,'Evaluation')
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
    _require_calibrated_model(saved,'Prediction')
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
    receipt = {'classifier_version':'5.8','model_version':saved.get('version','unknown'),
               'rule_fallback_policy_version':'5.5_identity_rule_fallback_v1',
               'family_assignment_policy_version':(c._SETTINGS.get('family_assignment_policy') or saved.get('family_assignment_policy',{})).get('version',''),
               'family_assignment_policy_sha256':_policy_digest(c._SETTINGS.get('family_assignment_policy') or saved.get('family_assignment_policy',{})),
               'family_assignment_score_interpretation':'raw Random Forest vote fraction used as an evidence-gated decision score, not a probability',
               'family_similarity_profile':similarity_profile_name,
               'family_similarity_search_max_evalue':similarity_settings['search_max_evalue'],
               'family_similarity_assignment_thresholds':similarity_selection,
               'membership_model_uses_family_similarity_extension':False,
               'membership_model_uses_family_signature_crosswalk':False,
               'family_signature_crosswalk_sha256':saved.get('provenance',{}).get('adhesome_family_rule_crosswalk_sha256'),
               'family_signature_features_used':saved.get('feature_schema')==SCHEMA,
                'feature_schema':SCHEMA,'records':len(output),'class_counts':dict(Counter(r['classifier_call'] for r in output)),
                'family_assignment_counts':dict(Counter(r['family_call_policy_tier'] or 'no_family_call' for r in output)),
               'rule_candidate_status_counts':dict(Counter(r['rule_candidate_status'] for r in output)),
               'untrained_family_rule_candidates':sum(bool(r['rule_candidate_family_untrained']) for r in output),
               'rule_rf_conflicts':sum(bool(r['rule_conflicts_with_random_forest']) for r in output),
               'probability_calibration_policy':saved.get('probability_calibration',{}).get('method'),
               'family_probability_calibrated_classes':saved.get('probability_calibration',{}).get('family',{}).get('eligible_classes',[]),
               'model_sha256':c.sha(args.model/'adhesome_ml.joblib' if args.model.is_dir() else args.model),
               'input_sha256':c.sha(args.fasta),'threshold':args.threshold,'exact_sequence_match_used':False}
    args.out.with_suffix('.summary.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
    print(f'[classifier] COMPLETE: wrote {len(output)} predictions to {args.out}',flush=True)
    print(json.dumps(receipt,indent=2))
    return 0
