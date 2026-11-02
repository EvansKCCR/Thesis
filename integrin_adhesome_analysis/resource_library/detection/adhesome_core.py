"""Shared dependency-free utilities for Python >=3.9."""

# Resolve source resources independently of the caller working directory.
import sys as _layout_sys
from pathlib import Path as _LayoutPath
_LAYOUT_BUNDLE_ROOT = next(p for p in _LayoutPath(__file__).resolve().parents
                          if (p/'config/layout_paths.json').is_file() and (p/'pipeline_paths.py').is_file())
if str(_LAYOUT_BUNDLE_ROOT) not in _layout_sys.path:
    _layout_sys.path.insert(0, str(_LAYOUT_BUNDLE_ROOT))
from pipeline_paths import activate as _layout_activate, resource_path as _layout_resource
_layout_activate()

import csv
import hashlib
import json
import re
from pathlib import Path

VERSION = '1.0'
DEFERRED_EVIDENCE = ('phylogeny', 'phylogeny_concordant', 'topology', 'localization', 'subcellular_localization')


def detection_evidence(row):
    """Keep external context and phylogeny at the resume boundary."""
    present = [key for key in DEFERRED_EVIDENCE if str(row.get(key, '') or '').strip()]
    if present:
        raise ValueError('Detection cannot consume '+', '.join(present)+
                         ' evidence; return external context through adhesome_workflow.py resume')

def fasta(path):
    header, parts = None, []
    with open(_layout_resource(path), encoding='utf-8-sig') as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(parts).upper()
                header, parts = line[1:], []
            else:
                if header is None:
                    raise ValueError(f'{path}: sequence before FASTA header')
                parts.append(''.join(line.split()))
    if header is not None:
        yield header, ''.join(parts).upper()

def write_fasta(path, records):
    with open(_layout_resource(path), 'w', encoding='utf-8') as handle:
        for header, sequence in records:
            handle.write('>' + header + '\n')
            for start in range(0, len(sequence), 60):
                handle.write(sequence[start:start + 60] + '\n')

def sequence_dict(path):
    result = {}
    for header, sequence in fasta(path):
        identifier = header.split()[0]
        if identifier in result or not sequence:
            raise ValueError(f'Duplicate ID or empty sequence in {path}: {identifier}')
        result[identifier] = sequence
    if not result:
        raise ValueError(f'Empty FASTA: {path}')
    return result


def validate_protein_input(path,hmmer_limit=False):
    """Reject an overwhelmingly nucleotide FASTA before protein analysis starts.

    A/C/G/T are valid amino-acid symbols too. Require nucleotide alphabet,
    sufficient length, nucleotide diversity and overwhelming file-level support;
    do not reject ordinary proteins or short/simple low-complexity fixtures.
    """
    total=bases=suspect_count=suspect_bases=0;examples=[];oversized=[]
    for header,raw in fasta(path):
        sequence=raw.removesuffix('*');length=len(sequence);alphabet=set(sequence)
        total+=1;bases+=length
        if length>=100 and alphabet<=set('ACGTUN') and len(alphabet&set('ACGTU'))>=3:
            suspect_count+=1;suspect_bases+=length
            if len(examples)<3:examples.append(header.split()[0])
        if length>100000 and len(oversized)<3:oversized.append(f'{header.split()[0]} ({length} residues)')
    if not total:raise ValueError(f'No FASTA records in {path}')
    if suspect_count/total>=.9 and suspect_bases/max(bases,1)>=.9:
        raise ValueError(f'{path}: nucleotide/transcript FASTA supplied to a protein workflow '
            f'({suspect_count}/{total} nucleotide-only records; examples: {", ".join(examples)}). '
            'Use the matching protein/peptide FASTA, not mRNA_transcripts, CDS or genomic DNA. '
            'Do not rename or directly translate full mRNA to bypass this check.')
    if hmmer_limit and oversized:
        raise ValueError('HMMER comparison limit exceeded before scanning: '+', '.join(oversized)+
            '. Verify that this is a protein FASTA; genuine proteins above 100,000 residues need a dedicated long-protein analysis.')
    return dict(records=total,nucleotide_only_records=suspect_count)

def table(path):
    with open(_layout_resource(path), newline='', encoding='utf-8-sig') as handle:
        return list(csv.DictReader(handle))

def write_table(path, rows, fields):
    with open(_layout_resource(path), 'w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)

def digest(path):
    h = hashlib.sha256()
    with open(_layout_resource(path), 'rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

def save_json(path, value):
    _layout_resource(Path(path)).write_text(json.dumps(value, indent=2) + '\n', encoding='utf-8')

def pf(value):
    match = re.search(r'PF\d{5}', value or '')
    return match.group() if match else ''

def hmm_records(path):
    lines = []
    with open(_layout_resource(path), encoding='utf-8') as handle:
        for line in handle:
            lines.append(line)
            if line.strip() == '//':
                metadata = {}
                for item in lines:
                    fields = item.split()
                    if fields and fields[0] in ('NAME', 'ACC', 'LENG', 'GA'):
                        metadata[fields[0]] = fields[1:]
                if 'ACC' not in metadata or 'LENG' not in metadata:
                    raise ValueError('HMM record missing ACC or LENG')
                yield metadata, ''.join(lines)
                lines = []
    if any(line.strip() for line in lines):
        raise ValueError('Truncated HMM file: missing //')
