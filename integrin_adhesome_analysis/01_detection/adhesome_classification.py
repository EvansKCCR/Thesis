"""Comparative screening: retain evidence without equating recruitment with assignment."""
from collections import Counter, defaultdict
import json
from pathlib import Path
import re
from adhesome_core import sequence_dict, table, write_table, write_fasta, save_json, digest, pf, VERSION, detection_evidence
from adhesome_components import annotate_components, CANDIDATE_FIELDS
from reference_identity import identity_evidence, IDENTITY_FIELDS
from adhesome_motifs import load_library, scan_motifs, noncatalytic_alignment
from evidence_coverage import search_coverage, related_families, coverage_details, REPORT_FIELDS

EVIDENCE_KEYS = ('similarity', 'reciprocal', 'profile', 'motif')

def evidence_state(rows, key):
    values = {r.get(key, '').strip().lower() for r in rows} - {''}
    if values - {'true', 'false'}:
        raise ValueError(f'{key} evidence must be true, false, or empty')
    if values == {'true', 'false'}:
        return 'conflict'
    return {'true': 'supported', 'false': 'not_supported'}.get(next(iter(values), ''), 'not_tested')

def ordered_groups_present(hits,groups):
    position=0
    for hit in hits:
        if position<len(groups) and hit['pfam_accession'] in groups[position]:
            position+=1
    return position==len(groups)

def architecture_evidence(hits, rule):
    counts = Counter(r['pfam_accession'] for r in hits)
    alternatives = rule['domain_alternatives']
    def rank(group):
        return (sum(min(counts[p], n) for p, n in group.items()) / sum(group.values()),
                sum(counts[p] > 0 for p in group) / len(group))
    best = max(alternatives, key=rank)
    observed = {p for p in best if counts[p]}
    missing = ';'.join(f'{p}:{max(0,n-counts[p])}' for p,n in best.items() if counts[p] < n)
    complete = any(all(counts[p] >= n for p,n in group.items()) for group in alternatives)
    expected_order = rule.get('ordered_domains', [])
    position = 0
    for hit in hits:
        if position < len(expected_order) and hit['pfam_accession'] == expected_order[position]:
            position += 1
    order = 'not_configured' if not expected_order else ('consistent' if position == len(expected_order) else 'incomplete_or_reordered')
    if rule.get('ordered_domain_groups'):
        groups=rule['ordered_domain_groups']
        order='consistent' if ordered_groups_present(hits,groups) else 'incomplete_or_reordered'
    if order == 'not_configured' and all(len(group) == 1 for group in alternatives):
        order = 'not_applicable'
    order_reason = ('single_domain_type_or_repeat_counts_only' if order == 'not_applicable' else
                    'no_domain_order_rule' if order == 'not_configured' else 'configured_domain_order')
    support_order=any(ordered_groups_present(hits,g) for g in rule['support_order_alternatives']) if rule.get('support_order_alternatives') else order=='consistent'
    anchors = set(rule.get('screening_domains', []))
    fraction = len(observed) / len(best)
    # Repeat deficits and domain rearrangements affect confidence, not recruitment.
    partial = bool(anchors & counts.keys()) or (bool(observed) and fraction >= rule.get('min_partial_fraction', 0.5))
    diagnostic = bool(set(rule.get('diagnostic_domains', [])) & counts.keys())
    diagnostic = diagnostic or (rule.get('diagnostic', False) and complete)
    partial = partial or diagnostic
    specific = diagnostic or (complete and len(best) >= 2) or (rule.get('screening_specific', True) and bool(anchors & counts.keys()))
    return dict(complete=complete, partial=partial, diagnostic=diagnostic, specific=specific,support_order=support_order,
                fraction=fraction, copy_fraction=rank(best)[0], missing=missing, order=order, order_reason=order_reason,
                observed=';'.join(sorted(observed | (anchors & counts.keys()) | (set(rule.get('diagnostic_domains', [])) & counts.keys()))))

def classify_comparative(args):
    config = json.loads(Path(args.rules).read_text(encoding='utf-8'))
    rules = config['families']
    motif_path = (Path(args.rules).resolve().parent / config['motif_library']).resolve() if config.get('motif_library') else None
    motif_library = load_library(motif_path) if motif_path else {}
    sequences = sequence_dict(args.proteome)
    domains, evidence = defaultdict(list), defaultdict(list)
    for family, rule in rules.items():
        alternatives = rule['domain_alternatives']
        if not alternatives or any(not g or any(not isinstance(n,int) or n < 1 for n in g.values()) for g in alternatives):
            raise ValueError(f'Invalid domain alternatives: {family}')
        if not 0 < rule.get('min_partial_fraction', .5) <= 1 or rule['min_length'] < 1:
            raise ValueError(f'Invalid screening thresholds: {family}')
    for row in table(args.domains):
        identifier = row['query_name']
        if identifier not in sequences:
            raise ValueError(f'Unknown domain protein: {identifier}')
        if row.get('decision') not in (None, '', 'retained'):
            continue
        start, end = int(row['ali_from']), int(row['ali_to'])
        if not 1 <= start <= end <= len(sequences[identifier]):
            raise ValueError(f'Invalid domain coordinates for {identifier}')
        domains[identifier].append(row)
    evidence_paths = [p for p in [getattr(args,'homology',None), *(getattr(args,'evidence',None) or [])] if p]
    component_rows = table(args.crosswalk) if getattr(args,'crosswalk',None) else []
    component_symbols = {row['Official Symbol'] for row in component_rows}
    reference_counts, coverage_paths = search_coverage(getattr(args, 'homology', None), args.proteome)
    related = {family: related_families(family, rules, config.get('family_aliases', {})) for family in rules}
    for path in evidence_paths:
        for row in table(path):
            detection_evidence(row)
            component_only = row['family']=='component_reference' and row.get('component_symbol') in component_symbols
            if row['protein_id'] not in sequences or (row['family'] not in rules and not component_only) or not row.get('source'):
                raise ValueError('Evidence requires known protein/family and a source')
            for key in EVIDENCE_KEYS:
                evidence_state([row], key)
            evidence[(row['protein_id'],row['family'])].append(row)
    manifest = getattr(args, 'manifest', None)
    curated = {p for rule in rules.values() for group in rule['domain_alternatives'] for p in group}
    if manifest:
        curated.update(pf(r['Pfam accession']) for r in table(manifest))
    results, architecture_rows, motif_rows = [], [], []
    for identifier, sequence in sequences.items():
        first_result = len(results)
        hits = sorted(domains[identifier],key=lambda r:(int(r['ali_from']),int(r['ali_to'])))
        architecture = ';'.join(f'{r["pfam_accession"]}:{r["ali_from"]}-{r["ali_to"]}' for r in hits)
        if hits:
            architecture_rows.append(dict(protein_id=identifier,length=len(sequence),architecture=architecture))
        covered = set()
        for family, rule in rules.items():
            support = evidence[(identifier,family)]
            states = {key:evidence_state(support,key) for key in EVIDENCE_KEYS}
            context_support = [row for name in related[family] for row in evidence[(identifier, name)]]
            coverage = coverage_details(family, states, support, related[family], context_support, reference_counts, evidence_state)
            if coverage['similarity_reason'] == 'completed_search_no_qualifying_hit':
                states['similarity'] = 'no_qualifying_hit'
            positive = lambda key: states[key] == 'supported'
            resolved_identity,identity_details=identity_evidence(support,rule)
            screening=bool(rule.get('screening_only') or rule.get('classification_level')=='screen')
            arch = architecture_evidence(hits, rule)
            family_evidence = any(positive(key) for key in ('similarity','profile','reciprocal'))
            if rule.get('requires_family_evidence') and not family_evidence:
                continue
            if not arch['complete'] and not arch['partial'] and not family_evidence:
                continue
            covered.update(arch['observed'].split(';'))
            warnings = []
            if len(sequence) < rule['min_length']:
                warnings.append('below_expected_length_possible_fragment')
            if not arch['complete']:
                warnings.append('incomplete_expected_architecture')
            if arch['order'] == 'incomplete_or_reordered':
                warnings.append('domain_order_differs_or_incomplete')
            motif_seq = sequence[-rule['motif_c_terminal_window']:] if rule.get('motif_c_terminal_window') else sequence
            motif = ('detected' if all(re.search(p,motif_seq) for p in rule['motifs']) else 'not_detected') if rule.get('motifs') else 'not_configured'
            observations = scan_motifs(sequence,hits,rule,motif_library)
            motif_rows.extend(dict(protein_id=identifier,family=family,source=str(motif_path),**r) for r in observations)
            if observations:
                by_motif=defaultdict(set)
                for observation in observations:
                    by_motif[observation['motif']].add(observation['state'])
                detected=sum('detected' in states for states in by_motif.values())
                motif = ('detected' if detected==len(by_motif) else 'partially_detected' if detected else
                         'not_detected' if any('not_detected' in states for states in by_motif.values()) else 'not_tested_missing_domain')
                warnings.extend('motif_not_detected:'+name for name,states in by_motif.items() if 'not_detected' in states)
            if motif == 'not_detected':
                warnings.append('expected_motif_not_detected')
            conflicts = [key for key,state in states.items() if state == 'conflict']
            negatives = [key for key,state in states.items() if state == 'not_supported']
            warnings.extend('conflicting_' + key for key in conflicts)
            warnings.extend('negative_' + key for key in negatives)
            homology_support = [r for r in support if r.get('similarity','').lower() == 'true']
            anchor_only = bool(homology_support) and all(r.get('reference_tier') == 'computational_anchor' for r in homology_support)
            if anchor_only:
                warnings.append('computational_anchor_only_similarity')
            extended = (noncatalytic_alignment(support,hits,rule['noncatalytic_alignment'])
                        if rule.get('noncatalytic_alignment') else True)
            if not extended:
                warnings.append('noncatalytic_family_evidence_missing')
            motif_architecture = bool(rule.get('motif_corroborates_partial_architecture') and arch['partial'] and
                                      any(r['state']=='detected' and r['interpretation']=='family_corroboration' for r in observations))
            # Support is about family identity, not experimental adhesome membership.
            supported = (not screening and (not rule.get('requires_identity_evidence') or resolved_identity) and positive('reciprocal') and (positive('similarity') or positive('profile')) and
                         (arch['complete'] or arch['specific'] or motif_architecture) and not conflicts and not negatives and
                         extended and
                         (not rule.get('support_requires_order') or arch['support_order']) and
                         (not anchor_only or positive('profile')))
            if rule.get('requires_identity_evidence') and not resolved_identity:
                warnings.append('independent_reference_identity_missing')
            if screening:
                status, reason = 'unresolved', 'screening_bucket_not_a_family'
                warnings.append('not_an_orthology_class')
            elif supported:
                status, reason = 'supported_candidate', 'convergent_family_evidence'
            elif arch['specific'] or family_evidence or (rule.get('requires_identity_evidence') and arch['partial']):
                status, reason = 'provisional', 'family_hypothesis_requires_validation'
            else:
                status, reason = 'unresolved', 'generic_domains_retained_without_family_assignment'
            tier = ('diagnostic_domain' if arch['diagnostic'] else 'complete_architecture' if arch['complete'] and arch['specific']
                    else 'partial_architecture' if arch['specific'] else 'homology_or_profile' if family_evidence else 'generic_domain')
            results.append(dict(protein_id=identifier,family=family,assigned_family=family if supported else '',
                **identity_details,
                module=rule['module'],status=status,reason=reason,screening_tier=tier,length=len(sequence),
                minimum_length=rule['min_length'],length_check='within_expectation' if len(sequence)>=rule['min_length'] else 'shorter_than_expected',
                domain_match=arch['complete'],domain_type_fraction=arch['fraction'],domain_copy_fraction=arch['copy_fraction'],
                matched_domains=arch['observed'],missing_domains=arch['missing'],order_check=arch['order'],order_check_reason=arch['order_reason'],
                support_order_check='consistent' if arch['support_order'] else 'not_satisfied' if rule.get('support_requires_order') else 'not_required',
                motif_check=motif,reciprocal=states['reciprocal'],
                motif_details=';'.join(f"{r['motif']}:{r['state']}:{r['start']}-{r['end']}" for r in observations),
                similarity=states['similarity'],profile=states['profile'],**coverage,
                warnings=';'.join(warnings),architecture=architecture,
                evidence_sources=';'.join(sorted({r['source'] for r in support}))))
        unassigned = {r['pfam_accession'] for r in hits} & curated - covered
        if unassigned:
            results.append(dict(protein_id=identifier,family='unassigned_domain_candidate',assigned_family='',module='unassigned',
                status='unresolved',reason='curated_domains_without_matching_family_rule',screening_tier='generic_domain',
                length=len(sequence),matched_domains=';'.join(sorted(unassigned)),architecture=architecture,
                warnings='family_assignment_not_established'))
        component_support = evidence[(identifier,'component_reference')]
        if any(row.get('similarity','').lower()=='true' for row in component_support) and not any(
                row['status'] in ('supported_candidate','provisional') for row in results[first_result:]):
            results.append(dict(protein_id=identifier,family='component_reference_candidate',assigned_family='',
                module='crosswalk_components',status='provisional',reason='component_reference_hypothesis_requires_validation',
                screening_tier='component_reference_similarity',length=len(sequence),architecture=architecture,
                similarity=evidence_state(component_support,'similarity'),reciprocal=evidence_state(component_support,'reciprocal'),
                similarity_reason='component_evidence_only_not_family_assignment',reciprocal_reason='component_panel_relative_agreement',
                homology_check='completed' if reference_counts is not None else 'not_documented',
                order_check='not_applicable',order_check_reason='no_family_architecture_assigned',
                warnings='component_similarity_not_family_assignment',
                evidence_sources=';'.join(sorted({row['source'] for row in component_support}))))
    grouped = defaultdict(list)
    for row in results:
        grouped[row['protein_id']].append(row)
    for identifier, rows in grouped.items():
        supported = [r for r in rows if r['status']=='supported_candidate']
        if len(supported)>1:
            for row in supported:
                row.update(status='provisional',assigned_family='',reason='conflicting_family_assignments')
                row['warnings'] += ';conflicting_family_assignments'
    output = Path(args.out)
    output.mkdir(parents=True,exist_ok=True)
    fields = ['protein_id','family','assigned_family','module','status','reason','screening_tier','length','minimum_length',*IDENTITY_FIELDS,
              'length_check','domain_match','domain_type_fraction','domain_copy_fraction','matched_domains','missing_domains',
              'order_check','order_check_reason','support_order_check','motif_check','motif_details','reciprocal','similarity','profile',*REPORT_FIELDS,'warnings','architecture','evidence_sources']
    crosswalk=getattr(args,'crosswalk',None)
    if crosswalk:
        annotate_components(crosswalk,manifest,[r for hits in domains.values() for r in hits],evidence,results,output)
        fields+=CANDIDATE_FIELDS
    write_table(output/'candidates.csv',results,fields)
    coverage_rows = []
    for family in rules:
        family_rows = [row for row in results if row['family'] == family]
        coverage_rows.append(dict(family=family,
            search_completed=reference_counts is not None,
            reference_count=reference_counts.get(family, 0) if reference_counts is not None else '',
            related_reference_families=';'.join(sorted(name for name in related[family] if reference_counts is not None and reference_counts.get(name))),
            hypothesis_rows=len(family_rows),
            similarity_supported=sum(row['similarity']=='supported' for row in family_rows),
            no_qualifying_hit=sum(row['similarity']=='no_qualifying_hit' for row in family_rows),
            similarity_not_tested=sum(row['similarity']=='not_tested' for row in family_rows),
            calibration_action='add_family_specific_references' if reference_counts is not None and not reference_counts.get(family) else ''))
    write_table(output/'reference_coverage.csv',coverage_rows,list(coverage_rows[0]))
    write_table(output/'motif_observations.csv',motif_rows,['protein_id','family','motif','pattern','context','state','start','end','sequence','interpretation','source'])
    write_table(output/'architectures.csv',architecture_rows,['protein_id','length','architecture'])
    priorities = {'supported_candidate':0,'provisional':1,'unresolved':2}
    protein_rows = []
    for identifier, rows in sorted(grouped.items()):
        status = min((r['status'] for r in rows),key=priorities.get)
        protein_rows.append(dict(protein_id=identifier,status=status,
            assigned_family=';'.join(sorted({r['assigned_family'] for r in rows if r['assigned_family']})),
            family_hypotheses=';'.join(sorted({r['family'] for r in rows if r['family']!='unassigned_domain_candidate'})),
            warnings=';'.join(sorted({w for r in rows for w in r.get('warnings','').split(';') if w}))))
    write_table(output/'protein_summary.csv',protein_rows,['protein_id','status','assigned_family','family_hypotheses','warnings'])
    write_fasta(output/'candidates.faa',[(p,sequences[p]) for p in sorted(grouped)])
    prioritized = {r['protein_id'] for r in protein_rows if r['status']!='unresolved'}
    write_fasta(output/'prioritized_candidates.faa',[(p,sequences[p]) for p in sorted(prioritized)])
    save_json(output/'classification.json',dict(pipeline_version=VERSION,stage='detection',deferred_stages=['topology','localization','phylogeny'],mode='comparative',rule_status=config['status'],
        counts=dict(Counter(r['status'] for r in results)),protein_counts=dict(Counter(r['status'] for r in protein_rows)),
        input_proteins=len(sequences),retained_proteins=len(grouped),prioritized_proteins=len(prioritized),
        evidence_reporting=dict(similarity=dict(Counter(r.get('similarity','') for r in results)),
            similarity_reason=dict(Counter(r.get('similarity_reason','') for r in results)),
            order_check=dict(Counter(r.get('order_check','') for r in results))),
        checksums={str(p):digest(p) for p in [args.proteome,args.domains,args.rules,*evidence_paths,*coverage_paths,*([manifest] if manifest else []),*([crosswalk] if crosswalk else []),*([motif_path] if motif_path else [])]}))
