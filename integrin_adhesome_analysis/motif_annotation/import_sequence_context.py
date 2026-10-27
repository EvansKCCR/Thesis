"""Join DeepLoc, SignalP and CDD predictions to exact-topology-verified FASTA IDs."""
import csv
import re
from collections import defaultdict
from pathlib import Path


def load(context_dir, families, topology_audit, regions, include_all=False):
    root = Path(context_dir)
    exact = {r['sequence_id'] for r in topology_audit if r['status'] == 'exact_sequence_verified'}
    members = defaultdict(list)
    for family, records in families.items():
        for record in records:
            members[record['id']].append((family, record))
    for ident, entries in members.items():
        if len({r['sequence'] for _, r in entries}) != 1:
            raise ValueError('One ID has conflicting supplied sequences: ' + ident)
    with (root/'DeepLoc_adhesome_candidate_list.csv').open(encoding='utf-8-sig', newline='') as handle:
        deeploc = {}
        for row in csv.DictReader(handle):
            ident = row['Protein_ID']
            if ident in deeploc: raise ValueError('Duplicate DeepLoc ID: ' + ident)
            deeploc[ident] = row
    signalp = {}
    with (root/'SignalP_adhesome_candidate_list.txt').open(encoding='utf-8-sig') as handle:
        for line in handle:
            if not line.strip() or line.startswith('#'): continue
            cells = line.rstrip('\n').split('\t')
            if len(cells) < 5 or cells[1] not in ('OTHER', 'SP'):
                raise ValueError('Invalid SignalP row: ' + line[:100])
            if cells[0] in signalp: raise ValueError('Duplicate SignalP ID: ' + cells[0])
            match = re.search(r'CS pos:\s*(\d+)-(\d+)', cells[4])
            cleavage = int(match.group(1)) if match else None
            if cells[1] == 'SP' and cleavage is None: raise ValueError('SignalP SP without cleavage: ' + cells[0])
            signalp[cells[0]] = (cells[1], cleavage)
    context = {}
    for ident, entries in members.items():
        for family, record in entries:
            key = family, ident
            valid = ident in exact
            dl = deeploc.get(ident) if valid else None
            sp = signalp.get(ident) if valid else None
            cleavage = sp[1] if sp else None
            if cleavage is not None and not 1 <= cleavage < len(record['sequence']):
                raise ValueError('SignalP cleavage outside sequence: ' + ident)
            if sp and sp[0] == 'SP':
                regions[key].append(('signal_peptide', 1, cleavage))
                regions[key].append(('predicted_mature_chain', cleavage + 1, len(record['sequence'])))
            context[key] = dict(family=family, sequence_id=ident,
                verification='exact_topology_sequence' if valid else 'no_exact_topology_sequence',
                deeploc_localizations=dl['Localizations'] if dl else '',
                deeploc_cytoplasm=dl['Cytoplasm'] if dl else '',
                deeploc_extracellular=dl['Extracellular'] if dl else '',
                signalp_prediction=sp[0] if sp else '', signalp_cleavage_after=cleavage or '',
                cdd_specific_hits=0)
    hits = []
    seen = set()
    path = root/'NCBI_CDD_adhesome_hitdata.txt'
    with path.open(encoding='utf-8-sig', newline='') as handle:
        for line in handle:
            if not line.startswith('Q#'): continue
            cells = line.rstrip('\n').split('\t')
            if len(cells) < 10: continue
            category=cells[1].strip().lower()
            if category!='specific' and not include_all:continue
            match = re.search(r'>(\S+)', cells[0])
            if not match: raise ValueError('CDD hit without query ID')
            ident = match.group(1)
            if ident not in exact or ident not in members: continue
            a, b = int(cells[3]), int(cells[4])
            evalue = float(cells[5])
            if evalue<0:raise ValueError('Invalid CDD E-value: '+ident)
            if evalue>1e-3 and not include_all:continue
            accession, name = cells[7], cells[8]
            for family, record in members[ident]:
                if not 1 <= a <= b <= len(record['sequence']):
                    raise ValueError('CDD coordinates outside sequence: ' + ident)
                token = (family, ident, accession, a, b, category)
                if token in seen: continue
                seen.add(token)
                row = dict(family=family, sequence_id=ident, accession=accession,
                    short_name=name, start=a, end=b, evalue=evalue, bitscore=cells[6],
                    hit_type=category,model_identifier=cells[2],incomplete=cells[9],
                    superfamily=cells[10] if len(cells)>10 else '',
                    model_coverage='',specific_bitscore_cutoff='',specific_evalue_cutoff='',
                    category_source='NCBI_CD_Search_hitdata',
                    interpretation='CDD_specific_domain_prediction_not_family_or_motif_validation' if category=='specific' else 'CDD_'+category+'_domain_prediction')
                hits.append(row)
                if category=='specific' and evalue<=1e-3:
                    regions[family, ident].append(('CDD:'+accession, a, b))
                    context[family, ident]['cdd_specific_hits'] += 1
    return context, hits
