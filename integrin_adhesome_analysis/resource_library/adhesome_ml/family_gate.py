"""FASTA-derived family learning. Validation labels never enter fitting or thresholds."""
import hashlib
import html
import json
import re
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

SCHEMA = 'sequence_family_gate_v5_3_posthoc_sequence_rules'
from cdd_features import build_calibration_bank, similarity_features
from identity_rule_features import features as identity_rule_features
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
        if cdd and ident in cdd:
            vector.update(cdd[ident])
        vector.update(identity_rule_features(row,cdd.get(ident,{}) if cdd else {},
            c._BUNDLE/'config/annotation/family_identity_rules.json',c._domain_records))
        posthoc_features,source_sha=posthoc_sequence_features(row,cdd.get(ident,{}) if cdd else {},posthoc_script)
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


def decide(membership, ranking, row, profiles, counts, parents, policy, c):
    observed = domain_set(row, c)
    if not any(c.evidence_groups(row).values()):
        return 'unresolved', '', 'no_informative_features'
    if membership < policy['membership_cutoff']:
        return 'non_adhesome', '', 'membership_below_cutoff'
    top, score = ranking[0]
    margin = score - (ranking[1][1] if len(ranking) > 1 else 0)
    shared = observed & set(profiles[top])
    if shared and counts[top] >= policy['minimum_family_examples'] and score >= policy['specific_vote_cutoff'] and margin >= policy['specific_margin']:
        return ('broad' if any(top in ancestors(x, parents) for x in counts if x != top) else 'specific'), top, 'family_vote_margin_and_shared_domain_support'
    aggregated = defaultdict(float)
    for label, vote in ranking:
        for parent in ancestors(label, parents):
            aggregated[parent] += vote
    for parent, vote in sorted(aggregated.items(), key=lambda x: (-x[1], x[0])):
        related = [name for name, _ in ranking if parent in ancestors(name, parents)]
        if vote >= policy['broad_vote_cutoff'] and any(observed & set(profiles[name]) for name in related):
            return 'broad', parent, 'shared_parent_vote_and_domain_support'
    return ('domain_only' if observed else 'unresolved'), '', 'insufficient_family_resolution'


def infer(saved, data, c, membership_cutoff=None):
    vectors = [x[1] for x in data]
    binary = saved['binary']
    votes = binary.predict_proba(vectors)[:, list(binary.classes_).index('adhesome_associated')]
    family_votes = saved['family'].predict_proba(vectors)
    families = list(saved['family'].classes_)
    policy = dict(saved['policy'])
    if membership_cutoff is not None:
        policy['membership_cutoff'] = membership_cutoff
    output = []
    for i, (row, _) in enumerate(data):
        ranking = sorted(zip(families, map(float, family_votes[i])), key=lambda x: (-x[1], x[0]))
        category, assignment, reason = decide(float(votes[i]), ranking, row, saved['profiles'], saved['counts'], saved['parents'], policy, c)
        output.append({'protein_id': row['protein_id'], 'classifier_call': category,
            'predicted_family': assignment, 'best_family_candidate': ranking[0][0],
            'family_vote_fraction': ranking[0][1], 'family_vote_margin': ranking[0][1] - (ranking[1][1] if len(ranking)>1 else 0),
            'adhesome_vote_fraction': float(votes[i]), 'family_ranking_json': json.dumps(ranking),
            'assignment_reason': reason, 'exact_sequence_match_used': 'no',
            'score_interpretation': 'uncalibrated model vote; not a biological probability',
            **c.prediction_evidence_summary(row)})
    return output


def evaluate(rows, predictions, sequences, saved, c):
    truth = [c.normalized_control_class(r) for r in rows]
    calls = ['adhesome_associated' if p['classifier_call'] in ('specific','broad','domain_only') else 'non_adhesome' for p in predictions]
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
    overlap = [sequence_digest(sequences[c.row_id(row)]) in saved['calibration_sequence_hashes'] for row in rows]
    strata = {}
    for state in (True, False):
        idx = [i for i,x in enumerate(overlap) if x==state]
        strata['sequence_overlap' if state else 'no_exact_sequence_overlap'] = c.class_metrics([truth[i] for i in idx], [calls[i] for i in idx], ['adhesome_associated','non_adhesome'])
    return {'membership': binary, 'family': family, 'class_counts': dict(Counter(p['classifier_call'] for p in predictions)),
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
    return {'cohort_records':len(rows),'same_denominator':True,'baseline':old,'posthoc_sequence_features':new,
        'delta_posthoc_minus_baseline':{k:new[k]-v for k,v in old.items() if isinstance(v,(int,float)) and not isinstance(v,bool)},
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
    comparison=report.get('posthoc_feature_comparison')
    if comparison:
        metric_rows=[]
        old,new,delta=(comparison['baseline'],comparison['posthoc_sequence_features'],comparison['delta_posthoc_minus_baseline'])
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
        document += '<h2>Baseline vs post-hoc sequence features</h2><p>Same unique validation cohort and denominators; baseline predictions are comparison-only.</p>'+table([('Metric','Baseline | v5.3 | Change'),*metric_rows])
    document += '<p>Membership positives include specific, broad and domain-only calls. Unresolved calls count as negative for membership metrics. Family recall includes abstentions; membership performance does not establish family identity.</p><h2>Resolution tiers</h2>'+table(report['class_counts'].items())
    document += '<h2>Per-family counts</h2>'+table((label,f"{v['correct']} correct / {v['records']} labelled") for label,v in family['per_family'].items())
    document += '<details><summary>Full metrics and overlap strata</summary><pre>'+html.escape(json.dumps(report,indent=2))+'</pre></details>'
    (directory/'internal_validation_report.html').write_text(document,encoding='utf-8')


def train_family_gate(args, c):
    aliases, parents, motifs = load_resources(c)
    vocab, provenance = c.configured_vocab(c._BUNDLE)
    rows, sequences, receipts = read_cohort(args.calibration, c)
    library=c._BUNDLE/'resource_library/adhesome_ml/cdd_domain_library_v1'
    with tempfile.TemporaryDirectory(prefix='adhesome_cdd_train_') as temp:
        bank_fasta,bank_meta=build_calibration_bank(library,args.calibration,Path(temp))
        calibration_queries=[{'protein_id':c.row_id(r),'sequence':sequences[c.row_id(r)]} for r in rows]
        train_cdd=similarity_features(calibration_queries,bank_fasta,bank_meta,exclude_self=True)
        data, scan_hashes = extract(rows, sequences, args.calibration/'handoff/detection', motifs, vocab, c,cdd=train_cdd)
        frozen_cdd_fasta=bank_fasta.read_bytes()
        frozen_cdd_meta=(Path(temp)/'cdd_calibration_references.json').read_bytes()
        val_rows=val_sequences=val_receipts=val_hashes=val_data=None
        if args.validation:
            val_rows,val_sequences,val_receipts=read_cohort(args.validation,c)
            val_queries=[{'protein_id':c.row_id(r),'sequence':val_sequences[c.row_id(r)]} for r in val_rows]
            val_cdd=similarity_features(val_queries,bank_fasta,bank_meta)
            val_data,val_hashes=extract(val_rows,val_sequences,args.validation/'handoff/detection',motifs,vocab,c,cdd=val_cdd)
    labels = [c.normalized_control_class(r) for r in rows]
    if set(labels) != {'adhesome_associated','non_adhesome'}:
        raise ValueError('Both positive and negative calibration controls are required')
    def estimator(seed):
        settings = dict(c._SETTINGS['binary_model'], random_state=seed)
        return c.Pipeline([('vectorizer',c.DictVectorizer()),('model',c.random_forest(settings,args.jobs))])
    binary = estimator(1729).fit([v for _,v in data], labels)
    indices = [(i,trusted_family(row,aliases,parents)) for i,row in enumerate(rows) if labels[i]=='adhesome_associated']
    indices = [(i,label) for i,label in indices if label and any(c.evidence_groups(data[i][0]).values())]
    if len({label for _,label in indices}) < 2:
        raise ValueError('At least two trusted, feature-supported calibration families are required')
    counts = Counter(label for _,label in indices)
    model_indices=[(i,label) for i,label in indices if counts[label]>=2]
    family = estimator(2718).fit([data[i][1] for i,_ in model_indices], [label for _,label in model_indices])
    profiles = {label: sorted(set().union(*(domain_set(data[i][0],c) for i,f in indices if f==label))) for label in counts}
    saved = dict(version='5.3',feature_schema=SCHEMA,binary=binary,family=family,config_vocab=sorted(vocab),
        aliases=aliases,parents=parents,motifs=motifs,counts=dict(counts),profiles=profiles,
        family_model_counts=dict(Counter(label for _,label in model_indices)),
        single_example_family_profiles=sorted(label for label,n in counts.items() if n==1),
        policy=c._SETTINGS['family_gate'],calibration_sequence_hashes={sequence_digest(sequences[c.row_id(r)]) for r in rows},
        provenance={**provenance,**receipts,**scan_hashes,'labels':c.sha(args.calibration/'reports/predictions.tsv'),
            'cdd_domain_fasta':c.sha(library/'cdd_domain_sequences.faa'),
            'cdd_domain_metadata':c.sha(library/'cdd_domain_metadata.tsv'),
            'family_identity_rules':c.sha(c._BUNDLE/'config/annotation/family_identity_rules.json'),
            'posthoc_sequence_rules_sha256':scan_hashes['posthoc_sequence_rules_sha256']})
    summary = {'classifier_version':'5.3','feature_schema':SCHEMA,'calibration_records':len(rows),
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
        'cdd_policy':'CDD reference bank contains calibration-origin fragments only. Calibration sequences are scored with same-protein reference hits excluded; validation CDD fragments never enter the bank.',
        'optional_returned_policy':'Returned CDD hits are not required; candidate sequences are searched against the frozen calibration-only CDD fragment bank.'}
    val_predictions = []
    if args.validation:
        val_predictions = infer(saved,val_data,c)
        summary['validation'] = evaluate(val_rows,val_predictions,val_sequences,saved,c)
        summary['validation_provenance'] = {**val_receipts,**val_hashes,'labels':c.sha(args.validation/'reports/predictions.tsv')}
        if args.baseline_validation_predictions:
            baseline_predictions=read_baseline_predictions(args.baseline_validation_predictions,val_rows,c)
            summary['validation']['posthoc_feature_comparison']=compare_validation(
                val_rows,baseline_predictions,val_predictions,saved,c)
            summary['baseline_validation_predictions_sha256']=c.sha(args.baseline_validation_predictions)
        for row,prediction in zip(val_rows,val_predictions):
            prediction.update(expected_family=trusted_family(row,aliases,parents),control_class=row['control_class'],
                calibration_sequence_overlap=sequence_digest(val_sequences[c.row_id(row)]) in saved['calibration_sequence_hashes'])
    args.out.mkdir(parents=True,exist_ok=False)
    (args.out/'cdd_calibration_references.faa').write_bytes(frozen_cdd_fasta)
    (args.out/'cdd_calibration_references.json').write_bytes(frozen_cdd_meta)
    saved['provenance']['frozen_cdd_reference_fasta']=c.sha(args.out/'cdd_calibration_references.faa')
    saved['provenance']['frozen_cdd_reference_metadata']=c.sha(args.out/'cdd_calibration_references.json')
    c.joblib.dump(saved,args.out/'adhesome_ml.joblib')
    (args.out/'model_summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
    c.write_tsv(args.out/'family_training_labels.tsv',[{'protein_id':c.row_id(rows[i]),'family':label} for i,label in indices])
    if val_predictions:
        c.write_tsv(args.out/'validation_predictions.tsv',val_predictions)
        report = summary['validation']
        write_validation_report(args.out, report)
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
    if saved.get('feature_schema') != SCHEMA:
        raise ValueError('Evaluation requires a trained v5.3 post-hoc-sequence-rule family model')
    rows, sequences, receipts = read_cohort(args.validation, c)
    bank_path=args.model if args.model.is_dir() else path.parent
    bank_fasta=bank_path/'cdd_calibration_references.faa'; meta_path=bank_path/'cdd_calibration_references.json'
    if not bank_fasta.is_file() or not meta_path.is_file(): raise ValueError('Model lacks its frozen calibration-only CDD reference bank; retrain v5')
    bank_meta=json.loads(meta_path.read_text(encoding='utf-8'))
    queries=[{'protein_id':c.row_id(r),'sequence':sequences[c.row_id(r)]} for r in rows]
    cdd=similarity_features(queries,bank_fasta,bank_meta)
    data, hashes = extract(rows, sequences, args.validation/'handoff/detection', saved['motifs'], set(saved['config_vocab']), c,cdd=cdd)
    if hashes.get('posthoc_sequence_rules_sha256') != saved.get('provenance',{}).get('posthoc_sequence_rules_sha256'):
        raise ValueError('Post-hoc sequence-rule source differs from the one used to train this model')
    predictions = infer(saved, data, c)
    report = evaluate(rows, predictions, sequences, saved, c)
    for row, prediction in zip(rows, predictions):
        prediction.update(expected_family=trusted_family(row,saved['aliases'],saved['parents']),
            control_class=row['control_class'],
            calibration_sequence_overlap=sequence_digest(sequences[c.row_id(row)]) in saved['calibration_sequence_hashes'])
    args.out.mkdir(parents=True,exist_ok=False)
    c.write_tsv(args.out/'validation_predictions.tsv', predictions)
    write_validation_report(args.out, report)
    receipt = {'model_sha256':c.sha(path),'labels_sha256':c.sha(args.validation/'reports/predictions.tsv'),
        'records':len(rows),'validation_used_for_fitting':False,'sources':{**receipts,**hashes}}
    (args.out/'evaluation_manifest.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
    print(f"Validation report: {args.out/'internal_validation_report.html'}; records={len(rows)}; model unchanged")
    return 0


def predict_family_gate(args, c):
    if not args.fasta or not args.evidence_dir:
        raise ValueError('Family prediction requires --fasta and --evidence-dir for verified sequence-derived features')
    saved = c.joblib.load(args.model/'adhesome_ml.joblib' if args.model.is_dir() else args.model)
    if saved.get('feature_schema') != SCHEMA:
        raise ValueError('Retrain a v5.3 CDD, identity-rule and post-hoc-sequence-rule classifier; earlier models cannot use these features')
    records = c.read_fasta(args.fasta)
    if args.ids:
        requested = set(args.ids.read_text().splitlines())
        if requested - {r['protein_id'] for r in records}:
            raise ValueError('Requested IDs are absent from candidate FASTA')
        records = [r for r in records if r['protein_id'] in requested]
    if len(records)>250 and not args.batch:
        raise ValueError('Maximum 250 sequences; specify --batch for consecutive 250-record batches')
    sequences, _ = verified_sequences(args.evidence_dir,c)
    for row in records:
        if sequences.get(row['protein_id']) != row['sequence']:
            raise ValueError(f"Candidate sequence differs from its scanned FASTA: {row['protein_id']}")
    output = []
    cached = c.read_stage1_evidence(args.evidence_dir)
    posthoc_script=c._BUNDLE.parent/'schisto_pfam_cdd_reclassify/schisto_pfam_cdd_reclassify.py'
    current_posthoc_sha,_,_,_=load_posthoc_rules(str(posthoc_script))
    if current_posthoc_sha != saved.get('provenance',{}).get('posthoc_sequence_rules_sha256'):
        raise ValueError('Post-hoc sequence-rule source differs from the one used to train this model')
    bank_path=args.model if args.model.is_dir() else args.model.parent
    bank_fasta=bank_path/'cdd_calibration_references.faa'; meta_path=bank_path/'cdd_calibration_references.json'
    if not bank_fasta.is_file() or not meta_path.is_file(): raise ValueError('Model lacks its frozen calibration-only CDD reference bank; retrain v5')
    bank_meta=json.loads(meta_path.read_text(encoding='utf-8'))
    hypotheses = defaultdict(set)
    hypothesis_file = args.evidence_dir.parent/'hypotheses.tsv'
    if hypothesis_file.is_file():
        for hypothesis in c.read_tsv(hypothesis_file):
            hypotheses[hypothesis.get('protein_id','')].add(hypothesis.get('family',''))
    for start in range(0,len(records),250):
        chunk=records[start:start+250]
        cdd=similarity_features(chunk,bank_fasta,bank_meta)
        data,_ = extract(chunk,sequences,args.evidence_dir,saved['motifs'],set(saved['config_vocab']),c,cached,cdd)
        output.extend(infer(saved,data,c,args.threshold))
        print(f'Family classification: {len(output)}/{len(records)}',flush=True)
    if len(output)!=len(records):
        raise ValueError('Prediction cohort accounting failed')
    for prediction in output:
        prediction['stage1_hypotheses_context_only'] = ';'.join(sorted(hypotheses[prediction['protein_id']] - {''}))
    c.write_tsv(args.out,output)
    receipt = {'classifier_version':'5.3','feature_schema':SCHEMA,'records':len(output),'class_counts':dict(Counter(r['classifier_call'] for r in output)),
               'model_sha256':c.sha(args.model/'adhesome_ml.joblib' if args.model.is_dir() else args.model),
               'input_sha256':c.sha(args.fasta),'threshold':args.threshold,'exact_sequence_match_used':False}
    args.out.with_suffix('.summary.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(receipt,indent=2))
    return 0
