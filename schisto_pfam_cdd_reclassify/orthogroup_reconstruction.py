"""Measure full-length sequence recovery from saved orthogroups and panel FASTAs."""
from collections import Counter, defaultdict
import csv
import hashlib
import re

DETAIL_FIELDS = ['species','protein_id','original_id','reference_column','orthogroup','control_in_panel',
                 'match_type','identical_sequence_available','identical_sequence_in_group',
                 'same_species_members','identical_sequence_matches']


def entries(path):
    with path.open(encoding='utf-8-sig', newline='') as handle:
        yield from csv.DictReader(handle, delimiter='\t')


def panel_inputs(orthogroups, panel, controls, target_columns, overrides):
    with orthogroups.open(encoding='utf-8-sig', newline='') as handle:
        columns = next(csv.reader(handle, delimiter='\t'), [])
    if not columns or columns[0] != 'Orthogroup' or len(columns) != len(set(columns)):
        raise ValueError('Expected a saved Orthogroups.tsv with unique species columns')
    if controls and not set(columns) & target_columns:
        raise ValueError('Orthogroups.tsv contains no target columns from this control manifest')
    species = {r['species'] for r in controls}
    mapping = {}
    for value in overrides:
        name, sep, column = value.partition('=')
        name, column = name.strip(), column.strip()
        if not sep or name not in species or name in mapping:
            raise ValueError('Use one --reference-column "species=orthogroup_column" per control species: '+value)
        if column and (column not in columns[1:] or column in target_columns):
            raise ValueError('Invalid same-species reference column: '+column)
        mapping[name] = column
    for name in sorted(species - mapping.keys()):
        slug = re.sub(r'[^A-Za-z0-9]+', '_', name).strip('_')
        found = [c for c in columns[1:] if c not in target_columns and
                 (c == slug or c == 'Bench_'+slug or c.startswith('Bench_'+slug+'_'))]
        if len(found)>1:
            raise ValueError('Multiple matching reference columns for '+name+'; use --reference-column')
        mapping[name] = found[0] if found else ''
    selected = ({r['Specie'] for r in controls} & set(columns)) | (set(mapping.values()) - {''})
    files = {}
    for column in sorted(selected):
        if not re.fullmatch(r'[A-Za-z0-9_.-]+', column) or column in {'.','..'}:
            raise ValueError('Unsafe panel column name: '+column)
        path = panel/(column+'.faa')
        if not path.is_file():
            raise ValueError('Missing saved input FASTA for '+column+': '+str(path))
        files[column] = path
    return mapping, files


def sequence_hashes(path):
    result = {}
    def save(identifier, parts):
        seq = ''.join(parts).upper().removesuffix('*')
        if not identifier or not seq or identifier in result:
            raise ValueError('Empty or duplicate FASTA record in '+str(path)+': '+identifier)
        result[identifier] = hashlib.sha256(seq.encode('ascii')).hexdigest()
    identifier = None
    parts = []
    with path.open(encoding='utf-8-sig') as handle:
        for raw in handle:
            line = raw.strip()
            if line.startswith('>'):
                if identifier is not None:save(identifier, parts)
                words = line[1:].split()
                identifier, parts = words[0] if words else '', []
            elif line:
                if identifier is None:raise ValueError('Sequence before first FASTA header: '+str(path))
                parts.append(line)
        if identifier is not None:save(identifier, parts)
    if not result:raise ValueError('Empty panel FASTA: '+str(path))
    return result


def reconstruction(orthogroups, controls, mapping, files):
    hashes = {column: sequence_hashes(path) for column,path in files.items()}
    targets = defaultdict(dict)
    for row in controls:
        column, protein, digest = row['Specie'], row['protein_id'], row['sequence_sha256']
        if protein in targets[column]:raise ValueError('Duplicate positive control: '+protein)
        targets[column][protein] = row
        panel_hash = hashes.get(column,{}).get(protein)
        if panel_hash is not None and panel_hash!=digest:
            raise ValueError('Saved target FASTA disagrees with control sequence hash: '+protein)
    references = set(mapping.values()) - {''}
    ref_sequences = defaultdict(lambda:defaultdict(set))
    for column in references:
        for protein,digest in hashes[column].items():
            ref_sequences[column][digest].add(protein)
    grouped = {}
    for row in entries(orthogroups):
        if None in row or any(v is None for v in row.values()) or not row['Orthogroup']:
            raise ValueError('Malformed row in saved Orthogroups.tsv')
        present = []
        for column,wanted in targets.items():
            present.extend((column,p) for p in members(row.get(column,'')) if p in wanted)
        if not present:continue
        reference_members = {c:members(row[c]) for c in references}
        for column,proteins in reference_members.items():
            if proteins - hashes[column].keys():
                raise ValueError('Orthogroup reference members absent from the saved panel: '+column)
        for pair in present:
            column,protein = pair
            if pair in grouped:raise ValueError('Positive control occurs in multiple orthogroups: '+protein)
            if protein not in hashes.get(column,{}):raise ValueError('Grouped control is absent from saved panel: '+protein)
            grouped[pair] = (row['Orthogroup'],reference_members)
    details = []
    for row in controls:
        species,column,protein = row['species'],row['Specie'],row['protein_id']
        reference = mapping[species]
        digest = row['sequence_sha256']
        group, refs = grouped.get((column,protein),('',{}))
        group_members = refs.get(reference,set())
        identical_ids = ref_sequences[reference].get(digest,set()) if reference else set()
        matched_sequences = identical_ids & group_members
        match_type = ('reference_not_available' if not reference else
                      'identical_sequence_recovered' if matched_sequences else
                      'identical_sequence_available_not_recovered' if identical_ids else
                      'other_same_species_member' if group_members else
                      'no_same_species_member')
        details.append(dict(species=species,protein_id=protein,original_id=row['original_id'],
            reference_column=reference,orthogroup=group,control_in_panel=protein in hashes.get(column,{}),
            match_type=match_type,identical_sequence_available=bool(identical_ids),
            identical_sequence_in_group=bool(matched_sequences),same_species_members=len(group_members),
            identical_sequence_matches=';'.join(sorted(matched_sequences))))
    by_species = defaultdict(list)
    for row in details:by_species[row['species']].append(row)
    summary = [summarize(name, rows, mapping[name]) for name,rows in sorted(by_species.items())]
    summary.append(summarize('All positive controls',details,''))
    return summary, details


def members(value):
    return {item.strip() for item in value.split(',') if item.strip()}


def summarize(species, rows, reference):
    classes = Counter(r['match_type'] for r in rows)
    assessable = len(rows)-classes['reference_not_available']
    available = sum(r['identical_sequence_available'] for r in rows)
    recovered = classes['identical_sequence_recovered']
    return dict(species=species,reference_column=reference,positive_controls=len(rows),
        controls_in_panel=sum(r['control_in_panel'] for r in rows),
        assigned_orthogroup=sum(bool(r['orthogroup']) for r in rows),
        controls_with_reference=assessable,any_same_species_member=sum(bool(r['same_species_members']) for r in rows),
        identical_sequence_available=available,identical_sequence_recovered=recovered,
        identical_sequence_available_not_recovered=classes['identical_sequence_available_not_recovered'],
        other_same_species_member=classes['other_same_species_member'],
        no_same_species_member=classes['no_same_species_member'],reference_not_available=classes['reference_not_available'],
        recovery_among_with_reference=recovered/assessable if assessable else None,
        recovery_among_available_matches=recovered/available if available else None)
