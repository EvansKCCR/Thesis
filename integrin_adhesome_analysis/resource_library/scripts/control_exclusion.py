"""Hard/background negatives are exclusion references, not identity classes."""

# Resolve source resources independently of the caller working directory.
import sys as _layout_sys
from pathlib import Path as _LayoutPath
_LAYOUT_BUNDLE_ROOT = next(p for p in _LayoutPath(__file__).resolve().parents
                          if (p/'config/layout_paths.json').is_file() and (p/'pipeline_paths.py').is_file())
if str(_LAYOUT_BUNDLE_ROOT) not in _layout_sys.path:
    _layout_sys.path.insert(0, str(_LAYOUT_BUNDLE_ROOT))
from pipeline_paths import activate as _layout_activate, resource_path as _layout_resource
_layout_activate()


from collections import defaultdict
from functools import lru_cache
import hashlib
import re
from pathlib import Path
from calibrate_identity import fasta,expected_families

ROOT=_LAYOUT_BUNDLE_ROOT


@lru_cache(maxsize=1)
def controls():
    from openpyxl import load_workbook
    folder=_layout_resource(ROOT/'validation');labels=_layout_resource(folder/'calibration_set.xlsx');proteins=_layout_resource(folder/'calibration.faa')
    if not _layout_resource(labels).is_file() or not _layout_resource(proteins).is_file():raise ValueError('Missing calibration exclusion inputs: '+str(folder))
    book=load_workbook(labels,read_only=True,data_only=True);values=iter(book['control_identity'].values)
    fields=next(values);rows=[dict(zip(fields,row)) for row in values if any(v is not None for v in row)];book.close()
    sequences=fasta(proteins);negative_by_hash=defaultdict(list);negative_families=defaultdict(list);positive_families=set();classes=defaultdict(set)
    positive_hashes=set();positive_cdd=set()
    for row in rows:
        if 'positive' in str(row['class']).lower():
            positive_cdd.update(re.findall(r'\bcd\d+\s+(\S+)\s*\[',str(row.get('NCBI_CDD') or ''),re.I))
    for row in rows:
        category=str(row['class']).lower();families=expected_families(row,{})
        if 'negative' in category:
            for family in families:negative_families[family].append(row)
        elif 'positive' in category:positive_families.update(families)
        else:raise ValueError('Unknown control class: '+str(row['class']))
        seq=sequences.get(row['accession'])
        if not seq:continue
        digest=hashlib.sha256(seq.encode()).hexdigest();classes[digest].add('negative' if 'negative' in category else 'positive')
        if 'negative' in category:negative_by_hash[digest].append(row)
        else:positive_hashes.add(digest)
    if any(len(group)>1 for group in classes.values()):raise ValueError('The same sequence has positive and negative control labels')
    for family,members in negative_families.items():
        for row in members:
            row['negative_specific_models']=set(re.findall(r'\bcd\d+\s+(\S+)\s*\[',str(row.get('NCBI_CDD') or ''),re.I))-positive_cdd
    return dict(negative_by_hash=dict(negative_by_hash),positive_hashes=positive_hashes,positive_families=positive_families,
        negative_families=dict(negative_families))


def apply(rows):
    bank=controls();retained=[];excluded=[]
    for row in rows:
        exact=bank['negative_by_hash'].get(str(row.get('sequence_sha256') or ''),[]);diagnostic=[]
        if str(row.get('sequence_sha256') or '') not in bank.get('positive_hashes',set()) and row.get('domain_score_0_4')==4 and row.get('family_specific_CDD_support'):
            from family_identity import cdd_hits
            from family_identity import library
            from cdd_evidence import interpretation
            observed={hit['short_name'] for hit in cdd_hits(row) if not library().get('precision_safeguards') or interpretation(hit)['specific']}
            diagnostic=[m for m in bank['negative_families'].get(row['assigned_family'],[]) if observed&m['negative_specific_models']]
        matches=exact or diagnostic
        gates=(['exact_negative_control_sequence'] if exact else [])+(['diagnostic_negative_subfamily_and_complete_architecture'] if diagnostic else [])
        if matches:
            row.update(stage2_excluded=True,exclusion_gate=gates[0],exclusion_gates=';'.join(gates),
                control_class='; '.join(sorted({m['class'] for m in matches})),
                control_accession='; '.join(sorted({m['accession'] for m in matches})),
                exclusion_reason='Hard/background-negative control evidence excludes this protein from final adhesome classification',
                excluded_putative_family=row['assigned_family'],excluded_protein_name=row['protein_name'])
            excluded.append(row)
        else:retained.append(row)
    return retained,excluded
