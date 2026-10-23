"""Terminal evaluation only: trusted labels cannot affect analysis decisions."""
from collections import Counter, defaultdict
import csv
import html
import json
from pathlib import Path

from .datasets import dump, table
from evaluate_validation import binary_metrics, supported, rate
from calibrate_identity import expected_family_resolution


def read(path):
    with Path(path).open(encoding='utf-8-sig', newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def metrics(rows, policy):
    counts = Counter()
    for r in rows:
        actual = r['control_class'] == 'positive'; called = r[policy]
        counts['TP' if actual and called else 'FN' if actual else 'FP' if called else 'TN'] += 1
    result = binary_metrics(*(counts[k] for k in ('TP', 'FP', 'TN', 'FN')))
    result['FPR'] = rate(counts['FP'], counts['FP'] + counts['TN'])
    return result


def summarize_predictions(predictions, folder, kind):
    folder.mkdir(parents=True, exist_ok=True)
    policies = ('stage1_recruited', 'supported', 'provisional', 'including_provisional', 'before_exclusion')
    results = []; groups = [('all', 'all', predictions)]
    for field in ('species', 'control_class', 'length_bin', 'calibration_overlap', 'reference_overlap', 'training_overlap'):
        for value in sorted({str(r[field]) for r in predictions}):
            groups.append((field, value, [r for r in predictions if str(r[field]) == value]))
    for field, value, rows in groups:
        for policy in policies:
            results.append(dict(population=field, value=value, policy=policy, **metrics(rows, policy)))
    table(folder/'metrics.tsv', results); table(folder/'predictions.tsv', predictions)
    errors = [dict(r, error='false_positive' if r['control_class'] != 'positive' else
        'false_negative' if not r['supported'] else 'wrong_family') for r in predictions if
        (r['supported'] and r['control_class'] != 'positive') or
        (r['control_class'] == 'positive' and (not r['supported'] or
        (r['expected_families'] and r['predicted_family'] not in r['expected_families'].split(';'))))]
    table(folder/'errors.tsv', errors)
    # Family correctness: wrong accepted positive families contribute FP and FN.
    evaluable = [r for r in predictions if r['control_class'] != 'positive' or r['expected_families']]
    cells = Counter(); by_family = defaultdict(Counter)
    for r in evaluable:
        expected = set(r['expected_families'].split(';')) - {''}; predicted = r['predicted_family']
        correct = r['supported'] and predicted in expected
        if r['control_class'] == 'positive':
            cells['TP' if correct else 'FN'] += 1
            # Alternative gold labels are one admissible label set, not duplicated cases.
            gold_group = ';'.join(sorted(expected))
            by_family[gold_group]['TP' if correct else 'FN'] += 1
        if r['supported'] and not correct:
            cells['FP'] += 1; by_family[predicted]['FP'] += 1
    family_rows = []
    for name, c in sorted(by_family.items()):
        family_rows.append(dict(family=name, TP=c['TP'], FP=c['FP'], FN=c['FN'],
            precision=rate(c['TP'], c['TP']+c['FP']), recall=rate(c['TP'], c['TP']+c['FN']),
            F1=rate(2*c['TP'], 2*c['TP']+c['FP']+c['FN'])))
    table(folder/'family_metrics.tsv', family_rows)
    family = dict(TP=cells['TP'], FP=cells['FP'], FN=cells['FN'],
        micro_precision=rate(cells['TP'],cells['TP']+cells['FP']),
        micro_recall=rate(cells['TP'],cells['TP']+cells['FN']),
        micro_F1=rate(2*cells['TP'],2*cells['TP']+cells['FP']+cells['FN']),
        unmapped_positive_labels=sum(r['control_class']=='positive' and not r['expected_families'] for r in predictions))
    for key in ('precision', 'recall', 'F1'):
        values = [r[key] for r in family_rows if r[key] is not None]
        family['macro_'+key] = sum(values)/len(values) if values else None
        family['macro_'+key+'_denominator'] = len(values)
    table(folder/'family_summary.tsv', [family])
    grades = Counter(r['grade'] for r in predictions)
    summary = dict(kind=kind, records=len(predictions), grade_counts=dict(grades), family=family,
        supported=metrics(predictions,'supported'),
        interpretation='Internal self-test; no parameter fitting' if kind=='calibration' else
        'Final frozen-configuration evaluation. See overlap strata; exact overlaps are not independent.',
        previous_orthology_F1=.5411461870014839,
        comparison_note='Orthogroup pairwise F1 and candidate/family classification F1 measure different tasks; not a direct improvement comparison.')
    dump(folder/'summary.json', summary)
    render(folder/'report.html', summary, results)
    return summary


def render(path, summary, rows):
    esc = lambda v: html.escape(str(v))
    fields = list(rows[0]) if rows else []
    content = '<!doctype html><meta charset="utf-8"><title>Workflow 1.0 report</title><style>body{font:15px system-ui;margin:32px}td,th{padding:7px;border:1px solid #ccc}table{border-collapse:collapse}</style>'
    content += '<h1>'+esc(summary.get('kind','workflow'))+'</h1><pre>'+esc(json.dumps(summary,indent=2))+'</pre>'
    content += '<table><tr>'+''.join('<th>'+esc(f)+'</th>' for f in fields)+'</tr>'
    content += ''.join('<tr>'+''.join('<td>'+esc(r.get(f,''))+'</td>' for f in fields)+'</tr>' for r in rows)+'</table>'
    Path(path).write_text(content, encoding='utf-8')


def evaluate_branch(branch, manifest, audit, overlap_hashes, reference_hashes, training_hashes):
    """Join truth to completed predictions only. Missing recruitment remains FN."""
    from evaluate_validation import read_audit
    handoff = branch/'handoff'
    hypotheses = read(handoff/'hypotheses.tsv'); recruited = {r['protein_id'] for r in hypotheses}
    before = {}; final = {}; excluded = {}
    for code in sorted({r['Specie'] for r in manifest}):
        a,b,c = read_audit(audit,code); before.update(a); final.update(b); excluded.update(c)
    if set(before) != recruited:
        raise ValueError('Audit candidate pool differs from stage 1')
    output = []
    for source in manifest:
        pid = source['protein_id']; row = final.get(pid,{})
        expected = expected_family_resolution(source)['families'] if source['control_class']=='positive' else set()
        digest = source['sequence_sha256']; length = int(source['length'])
        name = row.get('putative_audited_family') or row.get('assigned_family') or 'unassigned'
        provisional = row.get('grade')=='C' and name not in ('','unassigned','unresolved')
        output.append(dict(source, expected_families=';'.join(sorted(expected)),
            predicted_family=name, stage1_recruited=pid in recruited, supported=supported(row),
            provisional=provisional, including_provisional=supported(row) or provisional,
            before_exclusion=supported(before.get(pid)), grade=row.get('grade','not_accepted'),
            exclusion_gate=excluded.get(pid,{}).get('exclusion_gate',''),
            domain_score=row.get('domain_score_0_4',''), architecture=row.get('architecture',''),
            orthogroup=row.get('Orthogroup') or row.get('orthogroup',''),
            orthology_score=row.get('orthology_score_0_3',''), ambiguity=row.get('assignment_confidence',''),
            length_bin='<100' if length<100 else '100-299' if length<300 else '300-999' if length<1000 else '>=1000',
            calibration_overlap=digest in overlap_hashes, reference_overlap=digest in reference_hashes,
            training_overlap=digest in training_hashes))
    return output


def target_report(audit, folder):
    folder.mkdir(parents=True,exist_ok=True)
    rows = read(audit/'audited_family_assignments.tsv')
    ids = [r['protein_id'] for r in rows]
    if len(ids) != len(set(ids)):
        raise ValueError('Atlas export contains duplicate candidate IDs')
    if not (audit/'Stage_2_output.xlsx').is_file():
        raise ValueError('Missing Atlas-compatible Stage 2 workbook')
    counts = Counter((r.get('Specie'),r.get('putative_audited_family') or r.get('assigned_family'),r.get('grade')) for r in rows)
    grouped = [dict(species=s,family=f,grade=g,candidates=n) for (s,f,g),n in sorted(counts.items())]
    table(folder/'candidate_counts.tsv',grouped)
    table(folder/'ambiguous_unresolved.tsv',[r for r in rows if r.get('grade') not in ('A','B') or r.get('assignment_confidence')=='ambiguous'])
    summary = dict(kind='target',candidates=len(rows),atlas_integrity='passed',
        atlas_workbook=str(audit/'Stage_2_output.xlsx'), atlas_tsv=str(audit/'audited_family_assignments.tsv'),
        grade_counts=dict(Counter(r.get('grade') for r in rows)))
    dump(folder/'summary.json',summary); render(folder/'report.html',summary,grouped)
    return summary
