"""Version 1.0: isolated, paired Stage 2 reference-control regression.

Reuses sealed joined evidence; never runs recruitment, external tools, model
fitting, or the production audit. Gold labels do not enter the classifier.
"""

# Resolve source resources independently of the caller working directory.
import sys as _layout_sys
from pathlib import Path as _LayoutPath
_LAYOUT_BUNDLE_ROOT = next(p for p in _LayoutPath(__file__).resolve().parents
                          if (p/'config/layout_paths.json').is_file() and (p/'pipeline_paths.py').is_file())
if str(_LAYOUT_BUNDLE_ROOT) not in _layout_sys.path:
    _layout_sys.path.insert(0, str(_LAYOUT_BUNDLE_ROOT))
from pipeline_paths import activate as _layout_activate, resource_path as _layout_resource
_layout_activate()


import argparse
from collections import Counter, defaultdict
import csv
import importlib.util
import json
import os
from pathlib import Path
import re
from unittest.mock import patch

import evaluate_validation as evaluation
import family_identity
import identity_probability
from validation_reporting import reports, policy_expected_families, policy_prediction

ROOT = _LAYOUT_BUNDLE_ROOT
AUDIT = _layout_resource(ROOT / 'annotation')


def local_path(value):
    """Read the same immutable WSL receipts from Windows or WSL."""
    text = str(value)
    if os.name == 'nt' and re.match(r'^/mnt/[a-z]/', text):
        return _layout_resource(Path(text[5].upper() + ':/' + text[7:]))
    return _layout_resource(Path(text))


def validate_pair(families, identities):
    """Fail on unsupported rule fields instead of silently ignoring a gate."""
    allowed = {'id', 'parents', 'family', 'protein_name', 'required_groups',
               'specific_cdd_pattern', 'ordered_groups', 'length_range',
               'maximum_model_counts', 'allow_diagnostic_dispatch'}
    if str(families.get('version', families.get('pipeline_version', '1.0'))) != '1.0' or str(identities.get('pipeline_version')) != '1.0':
        raise ValueError('Both configurations must use workflow version 1.0')
    required = allowed - {'length_range', 'maximum_model_counts', 'allow_diagnostic_dispatch'}
    from family_catalogue import METADATA_FIELDS,GATE_FIELDS
    allowed |= METADATA_FIELDS | GATE_FIELDS | {'provisional_min_domain_fraction','confounders','evidence_source'}
    rules = identities['rules']
    if len({r['id'] for r in rules}) != len(rules):
        raise ValueError('Duplicate identity rule IDs')
    scope = set(families['families']) | {r['family'] for r in rules} | {'component_reference_candidate','unassigned_domain_candidate'}
    from stage2_resolution import BROAD_LABELS
    scope |= set(BROAD_LABELS)
    scope |= set(json.loads((_layout_resource(ROOT / 'motif_annotation/Motif_library/context_rules.json')).read_text())['families'])
    for r in rules:
        if set(r) - allowed or required - set(r):
            raise ValueError('Unsupported/missing fields: ' + r['id'])
        if not r['parents'] or set(r['parents']) - scope:
            raise ValueError('Unconfigured parent route: ' + r['id'])
        if not r['required_groups']:
            raise ValueError('An identity requires an observed domain core: ' + r['id'])
        for g in r['required_groups']:
            if set(g) != {'models', 'count'} or not g['models'] or type(g['count']) is not int or g['count'] < 1:
                raise ValueError('Invalid domain core: ' + r['id'])
        if any(type(i) is not int or not 0 <= i < len(r['required_groups']) for i in r['ordered_groups']):
            raise ValueError('Invalid order gate: ' + r['id'])
        re.compile(r['specific_cdd_pattern'])
    for pattern in {**identities['pfam_cdd_aliases'], **identities['additional_cdd_domains']}.values():
        re.compile(pattern)


def read_verified(path, expected, hashes):
    digest = evaluation.sha(path)
    if not expected or digest != expected:
        raise ValueError('Source differs from its recorded SHA256: ' + str(path))
    hashes[str(_layout_resource(path).resolve())] = digest


def source_context(handoff, evidence, target):
    hashes = {}
    hand = json.loads((_layout_resource(handoff / 'workflow_manifest.json')).read_text())
    sealed = json.loads((_layout_resource(evidence.parent / 'workflow_manifest.json')).read_text())
    for p in (_layout_resource(handoff / 'workflow_manifest.json'), _layout_resource(evidence.parent / 'workflow_manifest.json')):
        hashes[str(p)] = evaluation.sha(p)
    read_verified(evidence, sealed['files'].get(evidence.name), hashes)
    read_verified(_layout_resource(handoff / 'hypotheses.tsv'), hand['files'].get('hypotheses.tsv'), hashes)
    detection = local_path(hand['detection'])
    if not detection.is_absolute():
        detection = handoff / detection
    config = _layout_resource(detection / target / 'proteome.json')
    settings = json.loads(_layout_resource(config).read_text())
    raw = local_path(settings['proteomes'][0]['path'])
    if not raw.is_absolute():
        raw = config.parent / raw
    prep = _layout_resource(detection / target / 'prepared')
    receipt = json.loads((_layout_resource(prep / 'preparation.json')).read_text())['species'][target]
    read_verified(raw, receipt['input_sha256'], hashes)
    read_verified(prep / (target + '.faa'), receipt['output_sha256'], hashes)
    mapping_path = prep / (target + '.id_map.csv')
    source_hashes = {str(_layout_resource(local_path(p)).resolve()): sha for p, sha in hand.get('sources', {}).items()}
    for path in (config, _layout_resource(prep / 'preparation.json'), mapping_path,
                 _layout_resource(detection / target / 'candidates/candidates.csv')):
        # Preparation maps are linked by the detection receipt when present.
        expected = source_hashes.get(str(_layout_resource(path).resolve())) or hand['files'].get(path.relative_to(handoff).as_posix())
        if expected:
            read_verified(path, expected, hashes)
        else:
            hashes[str(_layout_resource(path).resolve())] = evaluation.sha(path)
    mapping = evaluation.read_table(mapping_path, ',')
    if len({r['original_id'] for r in mapping}) != len(mapping):
        raise ValueError('Duplicate original ID in preparation map')
    initial = defaultdict(set)
    for row in evaluation.read_table(_layout_resource(handoff / 'hypotheses.tsv')):
        if row.get('Specie') == target:
            initial[row['protein_id']].add(row['family'])
    csv.field_size_limit(100_000_000)
    with _layout_resource(evidence).open(encoding='utf-8-sig', newline='') as handle:
        joined = [r for r in csv.DictReader(handle, delimiter='\t') if r.get('Specie') == target]
    if {r['protein_id'] for r in joined} != set(initial):
        raise ValueError('Joined evidence changes the original exported target pool')
    sequences = {h.split()[0]: s for h, s in evaluation.records(prep / (target + '.faa'))}
    for row in joined:
        sequence = sequences[row['protein_id']]
        if evaluation.sequence_hash(sequence) != row['sequence_sha256']:
            raise ValueError('Joined evidence sequence differs from prepared protein')
        # Parse every full sidecar field before classification (no Excel truncation).
        hits = json.loads(row.get('NCBI_CDD_hits') or '[]')
        if not isinstance(hits, list):
            raise ValueError('Malformed structured CDD evidence')
    return raw, {r['original_id']: r for r in mapping}, sequences, initial, joined, hashes


def predictions(raw, mapping, sequences, initial, gold, trusted, refs, reverse, training, before, final, excluded, tiers):
    result = []; seen = set()
    for header, sequence in evaluation.records(raw):
        original = header.split()[0]; acc = evaluation.accession(original)
        digest = evaluation.sequence_hash(sequence)
        if acc in seen or original not in mapping:
            raise ValueError('Duplicate/unmapped validation original ID: ' + original)
        seen.add(acc); unit = mapping[original]; pid = unit['protein_id']
        if unit['decision'] == 'retained' and sequences.get(pid) != sequence:
            raise ValueError('Prepared protein differs from raw control: ' + pid)
        actual = gold.get(acc, {})
        prior = before.get(pid, {}); current = final.get(pid, {}); exclusion = excluded.get(pid, {})
        gate1 = pid in initial
        all_gates = str(exclusion.get('exclusion_gates', '')).split(';')
        without_exact = {} if 'diagnostic_negative_subfamily_and_complete_architecture' in all_gates else prior
        expected = actual.get('expected_families', set())
        row = dict(accession=acc, protein_id=pid, original_id=original, sequence_sha256=digest,
            evaluation_status='unlabelled_control' if not actual else 'labelled_sequence_mismatch' if acc in trusted and trusted[acc] != digest else 'scored',
            label_sequence_verification='not_labelled' if not actual else 'accession_label_only' if acc not in trusted else 'trusted_fasta_match' if trusted[acc] == digest else 'trusted_sequence_mismatch',
            control_class=actual.get('control_class', ''), trusted_protein_name=actual.get('protein_name', ''),
            trusted_genes=actual.get('genes', ''), expected_families=';'.join(sorted(expected)),
            expected_broad_families=';'.join(sorted(actual.get('expected_broad_families', set()))),
            expected_specific_families=';'.join(sorted(actual.get('expected_specific_families', set()))),
            preparation_decision=unit['decision'], coverage_tier=tiers.get(acc, 'previously_configured_or_deferred'),
            family_rule_coverage=actual.get('family_rule_coverage', ''),
            family_label_source=actual.get('family_label_source', ''),
            family_label_review_required=actual.get('family_label_review_required', ''),
            family_label_review_reason=actual.get('family_label_review_reason', ''),
            legacy_catalogue_families_ignored=actual.get('legacy_catalogue_families_ignored', ''),
            stage1_recruitment=gate1, stage2_candidate_retention=gate1 and bool(current),
            audit_before_negative_exclusion=gate1 and evaluation.supported(prior),
            audit_without_exact_negative_exclusion=gate1 and evaluation.supported(without_exact),
            exact_exclusion_ablation_status='complete_gate_provenance',
            without_exact_exclusion_family=without_exact.get('assigned_family', 'unassigned'),
            dual_gate_supported=gate1 and evaluation.supported(current),
            dual_gate_strong=gate1 and evaluation.supported(current, strong=True),
            audit_including_provisional=gate1 and current.get('grade') in ('A', 'B', 'C') and current.get('family_confidence') != 'ambiguous' and current.get('assigned_family') not in ('', 'unassigned', None),
            audited_family=current.get('assigned_family', 'unassigned'), pre_exclusion_family=prior.get('assigned_family', 'unassigned'),
            grade=current.get('grade', ''), assignment_confidence=current.get('family_confidence', ''),
            protein_name=current.get('protein_name', ''), protein_name_scope=current.get('protein_name_scope', ''),
            negative_exclusion=bool(exclusion), exclusion_gate=exclusion.get('exclusion_gate', ''),
            exclusion_gates=exclusion.get('exclusion_gates', ''),
            exact_forward_reference_overlap=digest in refs, exact_reverse_reference_overlap=digest in reverse,
            identity_calibration_training_overlap=digest in training,
            stage1_family_hypotheses=';'.join(sorted(initial.get(pid, set()))))
        row.update(evaluation.resolution_decisions(current))
        for field in ('identity_gate_trace', 'family_support_level', 'annotation_resolution'):
            row[field] = current.get(field) or prior.get(field) or ''
        for field in ('architecture_evidence', 'domain_score_0_4', 'orthology_score_0_3', 'identity_rule', 'identity_competitors', 'equally_ranked_hypotheses'):
            row['review_' + field] = prior.get(field, '')
        row['review_assignment_confidence'] = prior.get('family_confidence', '')
        result.append(row)
    if {r['original_id'] for r in result} != set(mapping):
        raise ValueError('Raw cohort differs from preparation map')
    return result


def write_rows(path, rows):
    fields = list(dict.fromkeys(k for row in rows for k in row)) or ['status']
    evaluation.write_table(path, rows, fields)


def run_arm(name, out, module, rows, family_path, identity_path, context, gold_context, tiers):
    out.mkdir()
    families = json.loads(_layout_resource(family_path).read_text()); identities = json.loads(_layout_resource(identity_path).read_text())
    validate_pair(families, identities)
    family_identity.LIBRARY = identity_path; family_identity.library.cache_clear()
    families['_audit_policy'] = json.loads((_layout_resource(AUDIT / 'audit_policy.json')).read_text())
    # A model fitted to different rules is inapplicable; confidence gates remain intact.
    with patch.object(identity_probability, 'annotate', return_value={'status': 'disabled_for_paired_configuration_regression'}):
        audited, retained = module.classify_rows(rows, families)
    before = {r['protein_id']: r for r in module.collapse_rows(audited)}
    final = {r['protein_id']: r for r in retained}
    excluded = {r['protein_id']: r for r in before.values() if r.get('stage2_excluded')}
    if set(before) != set(context[3]) or set(final) & set(excluded) or set(final) | set(excluded) != set(before):
        raise ValueError('Stage 2 protein accounting failed')
    pred = predictions(*context[:4], *gold_context, before, final, excluded, tiers)
    configured = set(families['families']) | {r['family'] for r in identities['rules']}
    for r in pred:
        expected = set(r['expected_families'].split(';')) - {''}
        r['family_rule_coverage'] = 'configured' if expected & configured else 'outside_configured_family_gates' if expected else 'needs_mapping_review'
    metrics, confusions, family_metrics, families_cells = evaluation.aggregate(pred)
    diagnostic = reports(pred, metrics)
    for filename, table in [('control_predictions.tsv', pred), ('metrics.tsv', metrics),
            ('confusion_matrices.tsv', confusions), ('family_metrics.tsv', family_metrics),
            ('family_confusion_matrix.tsv', families_cells), *diagnostic.items()]:
        write_rows(out / filename, table)
    # Full rule diagnostics are separate from the compact cohort decision report.
    traces = []
    for row in audited:
        for trace in row.get('_rule_gate_trace', []):
            if trace['reached'] or trace['diagnostic_cdd_matched']:
                traces.append(dict(protein_id=row['protein_id'], stage1_family=row['preaudit_family'],
                    **{k: json.dumps(v, sort_keys=True) if isinstance(v, (dict, list)) else v for k, v in trace.items()}))
    write_rows(_layout_resource(out / 'rule_gate_trace.tsv'), traces)
    counts = dict(Counter(r['control_class'] for r in pred if r['evaluation_status'] == 'scored'))
    summary = dict(arm=name, family_rules=len(families['families']), identity_rules=len(identities['rules']),
        family_rules_sha256=evaluation.sha(family_path), identity_rules_sha256=evaluation.sha(identity_path),
        cohort_counts=counts, input_sequences=len(pred), recruited_proteins=len(before),
        excluded_proteins=len(excluded), retained_proteins=len(final), probability_model_refit=False,
        metrics=[r for r in metrics if r['population'] == 'all_controls'], family_metrics=family_metrics)
    (_layout_resource(out / 'summary.json')).write_text(json.dumps(summary, indent=2) + '\n')
    return pred, summary


def compare(baseline, expanded):
    base = {r['accession']: r for r in baseline}; changes = []
    for row in expanded:
        old = base[row['accession']]
        expected = set(row['expected_families'].split(';')) - {''}
        changes.append(dict(accession=row['accession'], control_class=row['control_class'],
            trusted_genes=row['trusted_genes'], coverage_tier=row['coverage_tier'], expected_families=row['expected_families'],
            stage1_recruitment=row['stage1_recruitment'],
            baseline_family=old['audited_family'], expanded_family=row['audited_family'],
            baseline_supported=old['dual_gate_supported'], expanded_supported=row['dual_gate_supported'],
            baseline_correct_supported=old['dual_gate_supported'] and old['audited_family'] in expected,
            expanded_correct_supported=row['dual_gate_supported'] and row['audited_family'] in expected,
            baseline_before_exclusion=old['audit_before_negative_exclusion'],
            expanded_before_exclusion=row['audit_before_negative_exclusion'],
            baseline_without_exact=old['audit_without_exact_negative_exclusion'],
            expanded_without_exact=row['audit_without_exact_negative_exclusion'],
            expanded_confidence=row['assignment_confidence'], expanded_identity_rule=row['review_identity_rule'],
            expanded_gate_trace=row['identity_gate_trace']))
    return changes


def write_report(out, summary, changes):
    """Keep measured coverage distinct from membership, truth and deployment."""
    lines = ['# Family coverage development comparison — version 1.0', '',
        'Both supplied expanded JSON files were installed together under `installed_pair/` and evaluated with the production Stage 2 classifier against unchanged, sealed validation evidence. The baseline uses the current 96-family / 57-identity configuration and the same current classifier, labels and exclusions. This is not a comparison against the historical 70-positive output.', '',
        '**This is a reference-overlapping control self-test, not independent held-out validation or an estimate for novel lineages.** No Stage 1 run, external analysis, numerical model fitting, production installation change or backup was performed.', '',
        '| Measure | Baseline | Expanded |', '|---|---:|---:|']
    for label, field in [('Supported positive membership (of 235)', 'TP'), ('False-positive negative controls (of 348)', 'FP'), ('Membership recall', 'recall'), ('Membership PPV', 'PPV'), ('Membership F1', 'F1')]:
        items = [next(r for r in summary[arm]['metrics'] if r['policy'] == 'dual_gate_supported')[field] for arm in ('baseline', 'expanded')]
        lines.append(f'| {label} | {items[0]:.4g} | {items[1]:.4g} |')
    lines += ['', f"Of the **{summary['new_identity_positive_controls']} newly covered identity-reference positives**, {summary['new_identity_recruited']} were recruited and **{summary['new_identity_correct_supported']} received correct supported family calls**. There were {len(summary['new_identity_incorrect_supported'])} incorrect supported calls in that cohort. Missing support is not counted as successful recovery.", '',
        'The five controls without independently mapped expected families remain outside family-correctness scoring. All 235 positives remain in membership-recall denominators.', '', '| Regression | Result |', '|---|---|']
    for r in summary['regression_checks']:
        detail = r.get('observed_family') or (f"FP {r['baseline_FP']} → {r['expanded_FP']}" if 'baseline_FP' in r else r.get('lost_accessions') or '')
        lines.append(f"| {r['check']} | {'PASS' if r['passed'] else 'FAIL'} {detail} |")
    lines += ['', '## Remaining identity-proposal gaps', '',
        'See `unresolved_identity_proposals.tsv` for the required rule, missing observed models, anchor, order and coordinate gates for every missed control. Calibration workbook annotations were used to define gold labels and the proposal cohort; they were never injected into the classifier as observed target evidence. A model missing from the target evidence is not supplied from its expected label.', '',
        '| Accession | Genes | Expected family | Observed disposition |', '|---|---|---|---|']
    for r in changes:
        if r['coverage_tier'] == 'identity_proposal' and not r['expanded_correct_supported']:
            lines.append(f"| {r['accession']} | {r['trusted_genes']} | {r['expected_families']} | {r['expanded_family']} ({r['expanded_confidence']}) |")
    lines += ['', '## Screening-only controls and source integrity', '',
        'The 16 architecture-screen positives are reported separately from the 115 identity proposals. Defensin-6 (`Q01524`) is a screening-only proposal and was not recruited by Stage 1. Stage 2 did not add it or change recruitment.', '',
        'The expanded pair contains changes to four existing identity rules (kindlin, PTPN1, LAR and ERM), as well as the 79 additions. Their exact before/after definitions are listed in `configuration_diff.json`; the original source JSON files are unchanged.', '',
        f"The expanded family configuration refers to **{len(summary['expanded_models_not_in_bundled_scan_library'])} Pfam models absent from the bundled detection scan library**. This development replay does not install it into Stage 1. `missing_model_review.tsv` distinguishes missing target models, models outside the scan panel and any observed accession lost in normalization. Additional scanning or returned domain evidence requires a separate evidence-generation step; gates were not weakened to bypass it.", '',
        'All input, configuration and classifier hashes are recorded in `summary.json` and were checked again after replay. Both arms use the complete original validation cohort. Metrics, family macro-F1, exact-reference/calibration-overlap strata and exclusion ablations are available under each arm. Production outputs were read only.', '']
    (_layout_resource(out / 'report.md')).write_text('\n'.join(lines), encoding='utf-8')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--handoff', type=Path, required=True)
    p.add_argument('--evidence', type=Path, required=True, help='Sealed joined_evidence.tsv from a completed audit')
    p.add_argument('--out', type=Path, required=True, help='New development directory outside the source runs')
    p.add_argument('--target', default='validation')
    p.add_argument('--family-rules', type=Path, default=_layout_resource(AUDIT / 'family_rules.json'))
    p.add_argument('--identity-rules', type=Path, default=_layout_resource(AUDIT / 'family_identity_rules.json'))
    p.add_argument('--coverage-audit', type=Path, default=_layout_resource(AUDIT / 'family_coverage_expansion_audit.tsv'))
    args = p.parse_args(); out = _layout_resource(args.out).resolve(); handoff = _layout_resource(args.handoff).resolve(); evidence = _layout_resource(args.evidence).resolve()
    if _layout_resource(out).exists():
        raise FileExistsError('Choose a new development output directory: ' + str(out))
    if any(out.is_relative_to(path) for path in (handoff, evidence.parent)):
        raise ValueError('Development output must be outside the source run directories')
    pair = (_layout_resource(args.family_rules).resolve(), _layout_resource(args.identity_rules).resolve())
    validate_pair(*(json.loads(_layout_resource(path).read_text()) for path in pair))
    context = source_context(handoff, evidence, args.target)
    hashes = context[-1]
    bank_path = _layout_resource(ROOT / 'detection/databases/pfam/adhesome_v1.0.hmm.metadata.json')
    bundled_models = set(json.loads(_layout_resource(bank_path).read_text())['models'])
    for path in (*pair, args.coverage_audit, evaluation.LABELS, evaluation.CONTROL_FASTA,
                 _layout_resource(ROOT / 'detection/config/family_rules.json'), family_identity.LIBRARY,
                 _layout_resource(AUDIT / 'audit_policy.json'), bank_path):
        hashes[str(_layout_resource(path).resolve())] = evaluation.sha(path)
    gold, trusted, refs, reverse = evaluation.load_gold(evaluation.LABELS, evaluation.CONTROL_FASTA)
    training = {r['sequence_sha256'] for r in json.loads(_layout_resource(identity_probability.DEFAULT_MODEL).read_text()).get('samples', [])} if _layout_resource(identity_probability.DEFAULT_MODEL).exists() else set()
    gold_context = (gold, trusted, refs, reverse, training)
    tiers = {}; proposed = {}
    for row in evaluation.read_table(args.coverage_audit):
        if not row['decision'].startswith('ADD_'):
            continue
        tier = 'identity_proposal' if row['decision'] == 'ADD_IDENTITY_PROPOSAL' else 'screen_only'
        for acc in row['control_accessions'].split(';'):
            if acc in tiers or acc not in gold or gold[acc]['control_class'] != 'positive':
                raise ValueError('Invalid or overlapping proposed positive cohort: ' + acc)
            if row['expected_family'] not in gold[acc]['expected_families']:
                raise ValueError('Proposal family conflicts with frozen independent label: ' + acc)
            tiers[acc] = tier; proposed[acc] = row['expected_family']
    spec = importlib.util.spec_from_file_location('expansion_audit', _layout_resource(AUDIT / 'integrated_adhesome_pipeline.py'))
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    baseline_identity = family_identity.LIBRARY
    for path in (Path(__file__), _layout_resource(ROOT/'scripts/family_identity.py'), _layout_resource(ROOT/'scripts/stage2_resolution.py'),
                 _layout_resource(ROOT/'scripts/control_exclusion.py'), _layout_resource(ROOT/'scripts/evaluate_validation.py'), _layout_resource(ROOT/'scripts/validation_reporting.py'),
                 _layout_resource(ROOT/'scripts/label_ontology.py'), _layout_resource(ROOT/'scripts/calibrate_identity.py'), _layout_resource(Path(module.__file__))):
        hashes[str(_layout_resource(path).resolve())] = evaluation.sha(path)
    out.mkdir(parents=True)
    installed = _layout_resource(out / 'installed_pair'); installed.mkdir()
    for source, name in zip(pair, ('family_rules.json', 'family_identity_rules.json')):
        (installed / name).write_bytes(_layout_resource(source).read_bytes())
    original_id = json.loads(_layout_resource(baseline_identity).read_text())
    expanded_id = json.loads(_layout_resource(pair[1]).read_text())
    original_by_id = {r['id']: r for r in original_id['rules']}
    expanded_by_id = {r['id']: r for r in expanded_id['rules']}
    config_diff = dict(added_identity_rule_ids=sorted(set(expanded_by_id) - set(original_by_id)),
        removed_identity_rule_ids=sorted(set(original_by_id) - set(expanded_by_id)),
        changed_existing_rules={k: dict(baseline=v, expanded=expanded_by_id.get(k)) for k, v in original_by_id.items() if expanded_by_id.get(k) != v},
        changed_existing_families=[k for k, v in json.loads((_layout_resource(ROOT / 'detection/config/family_rules.json')).read_text())['families'].items()
                                  if json.loads(_layout_resource(pair[0]).read_text())['families'].get(k) != v])
    (_layout_resource(out / 'configuration_diff.json')).write_text(json.dumps(config_diff, indent=2) + '\n')
    try:
        print('Replaying baseline Stage 2 on the frozen validation target...', flush=True)
        baseline, bs = run_arm('baseline', _layout_resource(out / 'baseline'), module, context[4],
            _layout_resource(ROOT / 'detection/config/family_rules.json'), baseline_identity, context, gold_context, tiers)
        print('Replaying expanded Stage 2 with both development configurations...', flush=True)
        expanded, es = run_arm('expanded', _layout_resource(out / 'expanded'), module, context[4],
            _layout_resource(installed / 'family_rules.json'), _layout_resource(installed / 'family_identity_rules.json'), context, gold_context, tiers)
    finally:
        family_identity.LIBRARY = baseline_identity; family_identity.library.cache_clear()
    changes = compare(baseline, expanded)
    write_rows(_layout_resource(out / 'control_comparison.tsv'), changes)
    write_rows(_layout_resource(out / 'new_identity_positives.tsv'), [r for r in changes if r['coverage_tier'] == 'identity_proposal'])
    tier_metrics = []
    for arm, rows in [('baseline', baseline), ('expanded', expanded)]:
        for tier in sorted(set(tiers.values()) | {'previously_configured_or_deferred'}):
            cohort = [r for r in rows if r['control_class'] != 'positive' or r['coverage_tier'] == tier]
            metrics, _, fam, _ = evaluation.aggregate(cohort)
            for r in fam:
                tier_metrics.append(dict(arm=arm, coverage_tier=tier, **r))
    write_rows(_layout_resource(out / 'coverage_tier_family_metrics.tsv'), tier_metrics)
    regressions = []
    for gene, expected, forbidden in [('PTPRF', 'LAR_RPTP', ''), ('NF2', 'NF2_merlin_like', 'ERM_like'), ('TRIO', 'RhoGEF_like', 'spectrin_like')]:
        matches = [r for r in changes if gene in re.split(r'[;,\s]+', str(r['trusted_genes']))]
        passed = len(matches) == 1 and matches[0]['expanded_supported'] and matches[0]['expanded_family'] == expected
        regressions.append(dict(check=gene, passed=passed, expected_family=expected, forbidden_family=forbidden,
            accession=matches[0]['accession'] if len(matches) == 1 else '', observed_family=matches[0]['expanded_family'] if len(matches) == 1 else 'missing_or_multiple'))
    for policy in ('audit_before_negative_exclusion', 'audit_without_exact_negative_exclusion', 'dual_gate_supported'):
        before = {r['accession'] for r in baseline if r['evaluation_status'] == 'scored' and r['control_class'] != 'positive' and r[policy]}
        after = {r['accession'] for r in expanded if r['evaluation_status'] == 'scored' and r['control_class'] != 'positive' and r[policy]}
        regressions.append(dict(check='no_new_false_positives__' + policy, passed=not (after - before),
            baseline_FP=len(before), expanded_FP=len(after), added_false_positive_accessions=';'.join(sorted(after - before))))
    lost = [r['accession'] for r in changes if r['baseline_correct_supported'] and not r['expanded_correct_supported']]
    regressions.append(dict(check='previous_correct_supported_retained', passed=not lost, lost_accessions=';'.join(lost)))
    wrong = [r['accession'] for r in changes if r['control_class'] == 'positive' and r['expected_families'] and r['expanded_supported'] and not r['expanded_correct_supported']]
    regressions.append(dict(check='no_incorrect_supported_mapped_positive_families', passed=not wrong, incorrect_accessions=';'.join(wrong)))
    defensin = [r for r in changes if r['accession'] == 'Q01524']
    regressions.append(dict(check='Defensin6_remains_separate_Stage1_gap',
        passed=len(defensin) == 1 and defensin[0]['coverage_tier'] == 'screen_only' and not defensin[0]['stage1_recruitment'] and not defensin[0]['expanded_supported']))
    write_rows(_layout_resource(out / 'regression_checks.tsv'), regressions)
    unchanged = all(evaluation.sha(_layout_resource(Path(path))) == sha for path, sha in hashes.items())
    if not unchanged:
        raise ValueError('A read-only source changed during development replay')
    new = [r for r in changes if r['coverage_tier'] == 'identity_proposal']
    summary = dict(pipeline_version='1.0', status='regressions_passed' if all(r['passed'] for r in regressions) else 'needs_review',
        scope='Reference-overlapping development self-test; not held-out validation or novel-lineage accuracy',
        baseline=bs, expanded=es, new_identity_positive_controls=len(new),
        new_identity_recruited=sum(r['stage1_recruitment'] for r in new),
        new_identity_correct_supported=sum(r['expanded_correct_supported'] for r in new),
        new_identity_incorrect_supported=[r['accession'] for r in new if r['expanded_supported'] and not r['expanded_correct_supported']],
        new_identity_not_correct_supported=[r['accession'] for r in new if not r['expanded_correct_supported']],
        regression_checks=regressions, original_sources_unchanged=unchanged, input_sha256=hashes,
        expanded_models_not_in_bundled_scan_library=sorted(set(re.findall(r'PF\d{5}',json.dumps(json.loads(_layout_resource(pair[0]).read_text())['families']))) - bundled_models),
        stage1_rerun=False, orthology_rerun=False, phylogeny_run=False, probability_fit=False,
        production_installation_modified=False, backups_created=False)
    unresolved = {r['accession']: r for r in changes if r['coverage_tier'] == 'identity_proposal' and not r['expanded_correct_supported']}
    pids = {r['protein_id']: r['accession'] for r in expanded}
    gap_rows = {}
    for trace in evaluation.read_table(_layout_resource(out / 'expanded/rule_gate_trace.tsv')):
        acc = pids[trace['protein_id']]
        if acc in unresolved and trace['family'] in unresolved[acc]['expected_families'].split(';'):
            key = (acc, trace['rule_id'])
            gap_rows[key] = dict(accession=acc, trusted_genes=unresolved[acc]['trusted_genes'], **trace)
    write_rows(_layout_resource(out / 'unresolved_identity_proposals.tsv'), list(gap_rows.values()))
    raw_models = defaultdict(set)
    for row in context[4]:
        raw_models[row['protein_id']].update(re.findall(r'PF\d{5}', row.get('architecture', '')))
    missing_models = []
    for row in gap_rows.values():
        for group in json.loads(row['missing_domains']):
            for model in group['models']:
                missing_models.append(dict(accession=row['accession'], protein_id=row['protein_id'],
                    rule_id=row['rule_id'], required_model=model, required_count=group['required'], observed_count=group['observed'],
                    model_in_scan_library=model in bundled_models, raw_target_pfam_observed=model in raw_models[row['protein_id']],
                    diagnosis='observed_model_lost_in_normalization' if model in raw_models[row['protein_id']] else
                              'model_not_in_bundled_scan_library' if model not in bundled_models and model.startswith('PF') else
                              'required_model_not_observed_in_target'))
    write_rows(_layout_resource(out / 'missing_model_review.tsv'), missing_models)
    (_layout_resource(out / 'summary.json')).write_text(json.dumps(summary, indent=2) + '\n')
    write_report(out, summary, changes)
    print(json.dumps({k: summary[k] for k in ('status', 'new_identity_positive_controls', 'new_identity_recruited', 'new_identity_correct_supported', 'new_identity_incorrect_supported', 'new_identity_not_correct_supported', 'regression_checks')}, indent=2))


if __name__ == '__main__':
    main()
