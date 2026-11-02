"""Read-only input contracts. Labels never enter sequence analysis inputs."""

# Resolve source resources independently of the caller working directory.
import sys as _layout_sys
from pathlib import Path as _LayoutPath
_LAYOUT_BUNDLE_ROOT = next(p for p in _LayoutPath(__file__).resolve().parents
                          if (p/'config/layout_paths.json').is_file() and (p/'pipeline_paths.py').is_file())
if str(_LAYOUT_BUNDLE_ROOT) not in _layout_sys.path:
    _layout_sys.path.insert(0, str(_LAYOUT_BUNDLE_ROOT))
from pipeline_paths import activate as _layout_activate, resource_path as _layout_resource
_layout_activate()


from collections import Counter, defaultdict
import csv
import hashlib
import json
from pathlib import Path
import re

ROOT = _LAYOUT_BUNDLE_ROOT
EXPECTED = {'calibration': {'positive': 234, 'hard_negative': 119, 'background_negative': 228},
            'validation': {'positive': 902, 'negative': 835}}
SPECIES = {'mus_musculus': 'Mus musculus', 'rattus_norvegicus': 'Rattus norvegicus',
           'danio_rerio': 'Danio rerio', 'gallus_gallus': 'Gallus gallus',
           'xenopus_tropicalis': 'Xenopus tropicalis'}


def sha(path):
    with _layout_resource(Path(path)).open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def sequence_hash(sequence):
    return hashlib.sha256(sequence.encode('ascii')).hexdigest()


def dump(path, value):
    _layout_resource(Path(path)).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n', encoding='utf-8')


def table(path, rows):
    fields = list(dict.fromkeys(k for r in rows for k in r)) or ['status']
    with _layout_resource(Path(path)).open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fields, delimiter='\t')
        writer.writeheader(); writer.writerows(rows)


def fasta(path):
    """Keep headers for species verification; never length-filter input proteins."""
    from evaluate_validation import records
    result = {}
    for header, sequence in records(_layout_resource(Path(path))):
        identifier = header.split()[0]
        if identifier in result or not sequence or set(sequence) - set('ACDEFGHIKLMNPQRSTVWYBXZJUO'):
            raise ValueError(f'Duplicate ID, empty or invalid protein: {identifier} in {path}')
        if len(sequence) >= 100 and set(sequence) <= set('ACGTN'):
            raise ValueError(f'Possible nucleotide input: {identifier}')
        result[identifier] = (header, sequence)
    if not result:
        raise ValueError('Empty protein FASTA: ' + str(path))
    return result


def label_rows(path):
    path = _layout_resource(Path(path))
    if path.suffix.lower() == '.xlsx':
        from openpyxl import load_workbook
        book = load_workbook(path, read_only=True, data_only=True)
        try:
            values = iter(book['control_identity'].values); fields = next(values)
            return [dict(zip(fields, row)) for row in values if any(v is not None for v in row)]
        finally:
            book.close()
    with _layout_resource(path).open(encoding='utf-8-sig', newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t' if path.suffix == '.tsv' else ','))


def category(value):
    key = re.sub('[^a-z]', '', str(value).lower())
    return {'positivecontrol': 'positive', 'hardnegativecontrol': 'hard_negative',
            'backgroundnegativecontrol': 'background_negative', 'negativecontrol': 'negative'}.get(key)


def load_dataset(kind, protein_file, labels_file, select_labelled=False):
    """Exact label/sequence match; selecting extras requires explicit policy."""
    seqs = fasta(protein_file); rows = label_rows(labels_file); result = []; seen = set()
    by_accession = defaultdict(list)
    for identifier in seqs:
        parts = identifier.split('|')
        by_accession[parts[1] if len(parts) > 2 else parts[0]].append(identifier)
    for label in rows:
        cat = category(label.get('class'))
        if cat not in EXPECTED[kind]:
            raise ValueError(f'{kind}: unknown/inappropriate class {label.get("class")}')
        if kind == 'calibration':
            matches = by_accession[str(label.get('accession') or '').strip()]
            if len(matches) != 1:
                raise ValueError('Calibration accession does not match one FASTA: ' + str(label.get('accession')))
            identifier = matches[0]
        else:
            identifier = str(label.get('fasta_id') or '').strip()
        if identifier not in seqs or identifier in seen:
            raise ValueError(f'{kind}: missing/duplicate labelled FASTA ID {identifier}')
        seen.add(identifier); header, sequence = seqs[identifier]
        length = label.get('length') if kind == 'calibration' else label.get('sequence_length')
        if length is None or int(length) != len(sequence):
            raise ValueError(f'{kind}: labelled length differs for {identifier}')
        species = 'Homo sapiens'
        if kind == 'validation' and cat == 'positive':
            code = identifier.split('|')[0]
            if code not in SPECIES or len(identifier.split('|')) != 3:
                raise ValueError('Unverified validation species/gene/protein ID: ' + identifier)
            species = SPECIES[code]
        elif not re.search(r'\btaxid=9606\b', header):
            raise ValueError('Missing verified human control taxid: ' + identifier)
        result.append(dict(original_id=identifier, species=species, control_class=cat,
            accession=str(label.get('accession') or label.get('mapped_accession') or ''),
            genes=label.get('genes') or '', protein_name=label.get('protein_name') or label.get('Protein_name') or '',
            trusted_family=label.get('trusted_family') or '', sequence=sequence,
            sequence_sha256=sequence_hash(sequence), length=len(sequence)))
    counts = dict(Counter(r['control_class'] for r in result))
    if counts != EXPECTED[kind]:
        raise ValueError(f'{kind}: expected {EXPECTED[kind]}, found {counts}')
    extras = sorted(set(seqs) - seen)
    if extras and not select_labelled:
        raise ValueError(f'{kind}: {len(extras)} unlabelled FASTA records; explicitly enable labelled-ID selection')
    classes = defaultdict(set)
    for r in result:
        classes[r['sequence_sha256']].add(r['control_class'] == 'positive')
    if any(len(c) > 1 for c in classes.values()):
        raise ValueError(kind + ': identical sequence has conflicting positive/negative labels')
    summary = dict(kind=kind, counts=counts, labelled_records=len(result), fasta_records=len(seqs),
        excluded_unlabelled_ids=extras, species=dict(Counter(r['species'] for r in result)),
        dataset_version=sha(labels_file), fasta_sha256=sha(protein_file), labels_sha256=sha(labels_file),
        protein_file=str(_layout_resource(Path(protein_file)).resolve()), labels_file=str(_layout_resource(Path(labels_file)).resolve()))
    return result, summary


def overlaps(left, right):
    by_hash = defaultdict(list)
    for r in left:
        by_hash[r['sequence_sha256']].append(r['original_id'])
    return [dict(sequence_sha256=r['sequence_sha256'], right_id=r['original_id'],
                 right_class=r.get('control_class', ''), left_ids=';'.join(by_hash[r['sequence_sha256']]))
            for r in right if r['sequence_sha256'] in by_hash]


def materialize(kind, rows, folder, qc):
    """Deterministic branch/species namespaces, no labels or names in FASTA headers."""
    folder.mkdir(parents=True, exist_ok=True)
    groups = defaultdict(list); mappings = []
    prefix = {'calibration': 'Cal', 'validation': 'Val', 'target': 'Tgt'}[kind]
    for r in rows:
        code = prefix + '_' + re.sub('[^A-Za-z0-9]', '_', r['species'])
        # Hash the original ID, not the sequence: preserve distinct biological units.
        unit = 'p' + hashlib.sha256(r['original_id'].encode()).hexdigest()[:24]
        pid = code + '__' + unit
        groups[code].append((unit, r['sequence']))
        mappings.append({k: v for k, v in dict(r, protein_id=pid, Specie=code).items() if k != 'sequence'})
    if len({r['protein_id'] for r in mappings}) != len(mappings):
        raise ValueError('Namespace collision')
    entries = []
    for code, proteins in sorted(groups.items()):
        path = folder/(code + '.faa')
        path.write_text(''.join(f'>{unit} gene={unit}\n{seq}\n' for unit, seq in proteins), encoding='ascii')
        entries.append(dict(code=code, species=next(r['species'] for r in mappings if r['Specie'] == code),
            assembly='labelled_control_protein_units', source_url='run input manifest',
            path=str(_layout_resource(path).resolve()), gene_regex=r'(?:^|\s)gene=([^\s]+)'))
    table(_layout_resource(folder/'sequence_manifest.tsv'), mappings)
    dump(_layout_resource(folder/'targets.json'), dict(qc=qc, proteomes=entries))
    return mappings


def reference_inventory(raw):
    """Freeze raw FASTAs, mappings, GFFs and ground truth by content, not path alone."""
    raw = _layout_resource(Path(raw)).resolve()
    files = {p.relative_to(raw).as_posix(): sha(p) for p in sorted(_layout_resource(raw).rglob('*')) if _layout_resource(p).is_file()}
    if not any(re.search(r'\.(fa|faa|fasta)(\.gz)?$', p, re.I) for p in files):
        raise ValueError('No reference FASTAs in ' + str(raw))
    return dict(root=str(raw), files=files, version=hashlib.sha256(json.dumps(files, sort_keys=True).encode()).hexdigest())


def reference_matches(rows, raw):
    """Exact overlap against the whole frozen orthology panel, not only search seeds."""
    import gzip
    wanted=defaultdict(list)
    for r in rows:wanted[r['sequence_sha256']].append(r['original_id'])
    result=[]
    for path in sorted(_layout_resource(Path(raw)).rglob('*')):
        if not _layout_resource(path).is_file() or not re.search(r'\.(fa|faa|fasta)(\.gz)?$',path.name,re.I):continue
        opener=gzip.open if path.suffix=='.gz' else open
        header=None;parts=[]
        def save():
            if header is None:return
            # Match the explicit original-benchmark stop-marker policy.
            sequence=''.join(parts).upper().removesuffix('*').replace('*','X')
            digest=sequence_hash(sequence)
            for identifier in wanted.get(digest,[]):
                result.append(dict(original_id=identifier,sequence_sha256=digest,
                    reference_file=str(_layout_resource(path).resolve()),reference_id=header.split()[0]))
        with opener(path,'rt',encoding='utf-8-sig') as handle:
            for line in handle:
                if line.startswith('>'):save();header=line[1:].strip();parts=[]
                elif line.strip():parts.append(''.join(line.split()))
        save()
    return result
