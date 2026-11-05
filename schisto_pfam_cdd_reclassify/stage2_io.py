"""Stage 2 input contract for the independent post hoc domain-evidence review."""
from __future__ import annotations

import csv
import re
from pathlib import Path

ALIASES = {
    'protein_id': ('protein_id', 'Protein ID'),
    'Specie': ('Specie', 'Species', 'species'),
    'family_support_level': ('family_support_level',),
    'putative_audited_family': ('putative_audited_family', 'assigned_family'),
    'Pfam_architecture': ('Pfam_architecture', 'architecture', 'pfam_architecture'),
    'NCBI_CDD_architecture': ('NCBI_CDD_architecture', 'cdd_architecture'),
    'Orthology_Orthogroup': ('Orthology_Orthogroup', 'orthogroup', 'Orthogroup'),
    'stage1_family_hypothesis': ('stage1_family_hypothesis', 'preaudit_family', 'retained_stage1_family'),
    'length': ('length', 'protein_length'),
    'DeepLoc_2.1': ('DeepLoc_2.1',), 'DeepTMHMM': ('DeepTMHMM',),
    'SignalP.v6.0': ('SignalP.v6.0',),
}
SUPPORTED = {'specific_supported', 'broad_supported'}
REVIEW = {'unresolved', 'unassigned', 'provisional', 'ambiguous', 'conflicting',
          'unsupported', 'not_supported'}
NEGATIVE = {'negative', 'hard_negative', 'background_negative', 'excluded', 'rejected'}


def token(value):
    return re.sub(r'[\s-]+', '_', str(value or '').strip().casefold())


def text(value):
    return '' if value is None else str(value)


def key(row):
    return row['Specie'].strip(), row['protein_id'].strip()


def excluded(row):
    if any(token(row.get(k)) in {'true', 'yes', '1', '1.0'}
           for k in ('stage2_excluded', 'negative_control_excluded')):
        return True
    if token(row.get('control_class')) in NEGATIVE:
        return True
    return any(token(row.get(k)) in NEGATIVE for k in
               ('family_support_level', 'assignment_status', 'classification_state'))


def integer_length(value):
    if not str(value or '').strip():
        return 0
    number = float(value)
    if not number.is_integer() or number < 0:
        raise ValueError('Protein length must be a nonnegative integer: '+str(value))
    return int(number)


def load_table(path, sheet, xlsx_loader):
    path = Path(path).resolve()
    if path.suffix.lower() == '.xlsx':
        columns, values = xlsx_loader(path, sheet)
        actual_sheet = sheet
    elif path.suffix.lower() in ('.tsv', '.csv'):
        with path.open(encoding='utf-8-sig', newline='') as handle:
            rows = list(csv.reader(handle, delimiter='\t' if path.suffix.lower()=='.tsv' else ','))
        if not rows:
            raise ValueError('Input table has no header')
        columns, values = rows[0], rows[1:]
        actual_sheet = None
    else:
        raise ValueError('Input must be Stage_2_output.xlsx or its audit TSV/CSV export')
    columns = [str(c).strip() for c in columns]
    if not columns or any(not c for c in columns) or len(columns) != len(set(columns)):
        raise ValueError('Source header contains blank or duplicate column names')
    raw = []
    for index, row in enumerate(values, 2):
        if not any(str(v or '').strip() for v in row):
            continue
        if len(row) > len(columns):
            raise ValueError(f'Row {index} extends beyond the source header')
        raw.append(dict(zip(columns, [*map(text, row), *(['']*(len(columns)-len(row)))])))
    mapping = {name: next((alias for alias in aliases if alias in columns), None)
               for name, aliases in ALIASES.items()}
    required = ['protein_id', 'Specie'] + (['family_support_level'] if raw else [])
    missing = [name for name in required if not mapping[name]]
    if missing:
        raise ValueError('Missing required Stage 2 fields: '+', '.join(missing))
    if raw and not any(mapping[n] for n in ('Pfam_architecture','NCBI_CDD_architecture')):
        raise ValueError('Stage 2 input requires at least one Pfam/CDD architecture column')
    canonical, originals, duplicates, decisions = {}, {}, [], []
    for index, original in enumerate(raw, 2):
        row = dict(original)
        for name, source in mapping.items():
            row[name] = original.get(source, '') if source else ''
        species, protein = key(row)
        if not species or not protein:
            raise ValueError(f'Row {index}: species and protein_id must both be nonempty')
        pair = species, protein
        if pair in canonical:
            if row != canonical[pair]:
                raise ValueError(f'Conflicting duplicate species/protein pair: {species}/{protein}')
            duplicates.append({'Protein ID':protein, 'kept_species':species,
                               'duplicate_species':species, 'conflicting_data':'False'})
            continue
        if excluded(row):
            decision, reason = 'negative_excluded', 'Explicit source exclusion retained'
        else:
            level = token(row['family_support_level'])
            if level in SUPPORTED:
                decision, reason = 'supported_preserved', 'Supported audit assignment retained'
            elif level in REVIEW:
                decision, reason = 'reviewed', 'Unresolved, ambiguous or provisional audit assignment'
            else:
                raise ValueError(f'Unknown or blank family_support_level for {species}/{protein}: '
                                 +repr(row['family_support_level']))
        integer_length(row['length'])
        canonical[pair], originals[pair] = row, original
        decisions.append(dict(Specie=species, protein_id=protein,
                              original_support=row['family_support_level'],
                              posthoc_status=decision, reason=reason))
    return dict(columns=columns, rows=list(canonical.values()), originals=originals,
                source_rows=len(raw), duplicates=duplicates, decisions=decisions,
                column_mapping=mapping, sheet=actual_sheet)


def compare_prior(current, prior_path):
    """Compare by species/protein; permit old single-species-free ID keys only if unique."""
    with Path(prior_path).open(encoding='utf-8-sig', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        fields = reader.fieldnames or []
        previous = list(reader)
    if 'Protein ID' not in fields:
        raise ValueError('Prior comparison TSV requires Protein ID')
    use_species = 'Species' in fields
    def indexed(rows):
        result = {}
        for row in rows:
            pair = (row.get('Species',''), row['Protein ID']) if use_species else ('',row['Protein ID'])
            if pair in result:
                raise ValueError('Ambiguous duplicate comparison key: '+str(pair))
            result[pair] = row
        return result
    old, new = indexed(previous), indexed(current)
    differences = []
    for pair in sorted(old.keys() & new.keys()):
        for field in fields:
            a, b = text(old[pair].get(field)), text(new[pair].get(field))
            if a != b:
                differences.append({'Species':pair[0], 'Protein ID':pair[1],
                                    'Column':field, 'Prior':a, 'Reproduced':b})
    return dict(prior_rows=len(previous), new_rows=len(current),
                prior_only=[list(p) for p in sorted(old.keys()-new.keys())],
                new_only=[list(p) for p in sorted(new.keys()-old.keys())],
                difference_count=len(differences), differences=differences,
                all_rows_match=not differences and old.keys()==new.keys())
