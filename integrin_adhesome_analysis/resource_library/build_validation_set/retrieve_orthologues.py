#!/usr/bin/env python3
"""Accession-first human -> nonhuman Ensembl one-to-one protein retrieval.
Python 3.9+; dependency: requests. No symbol-based rescue or best-hit fallback.
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
import csv
import hashlib
import json
import re
import sys
import time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import quote
import requests

DEFAULT_SPECIES = ['mus_musculus', 'rattus_norvegicus', 'danio_rerio',
                   'gallus_gallus', 'xenopus_tropicalis']

class API:
    def __init__(self, cache, refresh=False):
        self.cache = _layout_resource(Path(cache))
        self.cache.mkdir(parents=True, exist_ok=True)
        self.refresh = refresh
        self.session = requests.Session()
        self.session.headers.update({'User-Agent': 'strict-orthologues/1.0',
                                     'Accept': 'application/json'})
        self.last = 0

    def get(self, url, **params):
        key = hashlib.sha256(json.dumps([url, params], sort_keys=True).encode()).hexdigest()
        file = self.cache / (key + '.json')
        if _layout_resource(file).exists() and not self.refresh:
            return json.loads(_layout_resource(file).read_text())['body']
        for attempt in range(6):
            time.sleep(max(0, 0.25 - (time.monotonic() - self.last)))
            try:
                self.last = time.monotonic()
                r = self.session.get(url, params=params, timeout=(15, 120))
                if r.status_code == 429 or r.status_code >= 500:
                    try:
                        delay = float(r.headers.get('Retry-After', 2 ** attempt))
                    except ValueError:
                        delay = 2 ** attempt
                    time.sleep(min(max(delay, 1), 120))
                    continue
                r.raise_for_status()
                body = r.json()
                if isinstance(body, dict) and body.get('error'):
                    raise ValueError(str(body['error']))
                record = {'url': r.url, 'retrieved_utc': now(),
                          'response_headers': dict(r.headers), 'body': body}
                temp = file.with_suffix('.tmp')
                temp.write_text(json.dumps(record, indent=2))
                temp.replace(file)
                return body
            except (requests.ConnectionError, requests.Timeout):
                if attempt == 5:
                    raise
                time.sleep(2 ** attempt)
        raise RuntimeError('Retry limit exceeded: ' + url)

def now():
    return datetime.now(timezone.utc).isoformat()

def stable(identifier):
    return str(identifier).split('.')[0]

def read_input(path):
    with open(_layout_resource(path), encoding='utf-8-sig', newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        if not reader.fieldnames:
            raise ValueError('Input must be a headered TSV.')
        column = next((c for c in ('accession', 'input_accession') if c in reader.fieldnames), None)
        if not column:
            raise ValueError('Required column: accession or input_accession')
        rows = []
        seen = set()
        for row in reader:
            acc = row[column].strip()
            if not re.fullmatch(r'[A-Z0-9]+(?:-\d+)?', acc):
                raise ValueError('Invalid accession: ' + acc)
            if acc in seen:
                raise ValueError('Duplicate accession: ' + acc)
            seen.add(acc)
            rows.append({'accession': acc, 'input_genes': row.get('genes', row.get('selected_gene_label', ''))})
        if not rows:
            raise ValueError('Input contains no candidates.')
        return rows

def write_tsv(path, rows, fields):
    with open(_layout_resource(path), 'w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fields, delimiter='\t', extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)

def resolve(api, base, acc):
    # Query the accession, not potentially misleading aliases in the input.
    u = api.get('https://rest.uniprot.org/uniprotkb/' + quote(acc, safe='') + '.json')
    if u.get('organism', {}).get('taxonId') != 9606:
        raise ValueError('Source accession is not a verified human UniProt entry')
    if u.get('inactiveReason'):
        raise ValueError('Inactive UniProt entry; manual review required')
    primary = u.get('primaryAccession', '')
    names = [g.get('geneName', {}).get('value', '') for g in u.get('genes', [])]
    genes = set()
    for x in u.get('uniProtKBCrossReferences', []):
        if x.get('database') == 'Ensembl':
            for prop in x.get('properties', []):
                if prop.get('key') == 'GeneId' and str(prop.get('value', '')).startswith('ENSG'):
                    genes.add(stable(prop['value']))
    # Accession-linked Ensembl objects only. Never rescue using a gene symbol.
    if not genes:
        refs = api.get(base + '/xrefs/symbol/homo_sapiens/' + quote(primary, safe=''),
                       object_type='gene', **{'content-type': 'application/json'})
        genes = {stable(x['id']) for x in refs if x.get('type') == 'gene'
                 and str(x.get('id', '')).startswith('ENSG')}
    if len(genes) != 1:
        raise ValueError('Accession maps to %d Ensembl genes: %s' % (len(genes), ','.join(sorted(genes))))
    gene = next(iter(genes))
    info = api.get(base + '/lookup/id/' + gene, **{'content-type': 'application/json'})
    if info.get('species') != 'homo_sapiens' or info.get('object_type') != 'Gene':
        raise ValueError('Source gene lookup failed human-gene validation')
    return {'resolved_accession': primary, 'uniprot_gene': ';'.join(filter(None, names)),
            'uniprot_entry_type': u.get('entryType', ''), 'source_gene': gene,
            'source_gene_name': info.get('display_name', ''),
            'source_assembly': info.get('assembly_name', ''),
            'mapping_note': 'Gene-level mapping only; no exact isoform equivalence claimed.'}

def homologies(api, base, gene, compara):
    # Current species-qualified endpoint, with compatibility for older servers.
    params = {'type': 'orthologues', 'sequence': 'none', 'format': 'full',
              'compara': compara, 'content-type': 'application/json'}
    try:
        result = api.get(base + '/homology/id/homo_sapiens/' + gene, **params)
    except requests.HTTPError as e:
        if e.response is None or e.response.status_code not in (400, 404):
            raise
        result = api.get(base + '/homology/id/' + gene, **params)
    data = result.get('data')
    if not isinstance(data, list):
        raise ValueError('Malformed homology response: data missing')
    if not data:
        raise ValueError('Source gene absent from selected compara database')
    matching = [d for d in data if stable(d.get('id', '')) == stable(gene)]
    if len(matching) != 1 or not isinstance(matching[0].get('homologies'), list):
        raise ValueError('Malformed or mismatched source-gene homology response')
    return matching[0]['homologies']

def select(hits, species, high_confidence=False):
    found = [h for h in hits if h.get('target', {}).get('species') == species]
    if not found:
        return None, 'NO_ORTHOLOGUE_REPORTED'
    if any(h.get('type') != 'ortholog_one2one' for h in found):
        return None, 'EXCLUDED_NON_ONE2ONE'
    unique = {(h.get('target', {}).get('id'), h.get('target', {}).get('protein_id')) for h in found}
    if len(unique) != 1:
        return None, 'AMBIGUOUS_TARGET'
    h = found[0]
    if high_confidence and any(str(x.get('is_high_confidence', '')) != '1' for x in found):
        return None, 'CONFIDENCE_NOT_CONFIRMED'
    if not h.get('target', {}).get('id') or not h.get('target', {}).get('protein_id'):
        return None, 'MISSING_TARGET_IDENTIFIERS'
    return h, 'ONE2ONE_REPORTED'

def validate_sequence(seq):
    seq = ''.join(str(seq).split()).upper()
    if not seq or not re.fullmatch(r'[ACDEFGHIKLMNPQRSTVWYBXZJUO]+', seq):
        raise ValueError('Empty/invalid protein sequence, gap or stop character')
    return seq

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input', default='candidates.tsv')
    p.add_argument('--out', default='orthologue_results')
    p.add_argument('--species', nargs='+', default=DEFAULT_SPECIES,
                   help='Exact Ensembl species names, space-separated; humans forbidden')
    p.add_argument('--server', default='https://rest.ensembl.org')
    p.add_argument('--compara', default='vertebrates')
    p.add_argument('--require-high-confidence', action='store_true',
                   help='Also require is_high_confidence=1; absent flags are excluded')
    p.add_argument('--refresh', action='store_true', help='Ignore cached responses')
    args = p.parse_args()
    species = list(dict.fromkeys(args.species))
    if any(s in ('human', 'homo_sapiens') or not re.fullmatch('[a-z]+_[a-z0-9_]+', s) for s in species):
        p.error('Use exact nonhuman Ensembl species names.')
    out = _layout_resource(Path(args.out))
    out.mkdir(parents=True, exist_ok=True)
    # Avoid stale FASTA files being mistaken for this run's output.
    if (_layout_resource(out / 'run_metadata.json')).exists():
        p.error('Output already contains a run. Choose a new --out; reuse cache with a symlink if desired.')
    api = API(_layout_resource(out / 'raw_cache'), args.refresh)
    base = args.server.rstrip('/')
    inputs = read_input(args.input)
    meta = {'started_utc': now(), 'arguments': vars(args), 'input_sha256': hashlib.sha256(_layout_resource(Path(args.input)).read_bytes()).hexdigest(),
            'scope': 'Gene-level Ensembl ortholog_one2one; Compara representative target protein, not exact source isoform.',
            'input_count': len(inputs)}
    for endpoint in ('info/data', 'info/software', 'info/rest'):
        try:
            meta[endpoint] = api.get(base + '/' + endpoint, **{'content-type': 'application/json'})
        except Exception as e:
            meta[endpoint] = {'error': str(e)}
    (_layout_resource(out / 'run_metadata.json')).write_text(json.dumps(meta, indent=2))
    mapping, audit, records = [], [], {}
    for n, inp in enumerate(inputs, 1):
        print('[%d/%d] %s' % (n, len(inputs), inp['accession']), flush=True)
        m = dict(inp)
        try:
            m.update(resolve(api, base, inp['accession']))
            m['mapping_status'] = 'RESOLVED'
        except Exception as e:
            m.update(mapping_status='MAPPING_UNRESOLVED', error=str(e))
            mapping.append(m)
            for sp in species:
                audit.append(dict(inp, target_species=sp, status='MAPPING_UNRESOLVED', error=str(e)))
            continue
        mapping.append(m)
        try:
            hits = homologies(api, base, m['source_gene'], args.compara)
        except Exception as e:
            for sp in species:
                audit.append(dict(m, target_species=sp, status='HOMOLOGY_QUERY_ERROR', error=str(e)))
            continue
        for sp in species:
            a = dict(m, target_species=sp)
            h, status = select(hits, sp, args.require_high_confidence)
            a['status'] = status
            a['reported_types'] = ';'.join(sorted({x.get('type', '') for x in hits if x.get('target', {}).get('species') == sp}))
            if h:
                t = h['target']
                a.update(target_gene=t['id'], target_protein=t['protein_id'],
                         orthology_type=h['type'], is_high_confidence=h.get('is_high_confidence', ''),
                         source_percent_identity=h.get('source', {}).get('perc_id', ''),
                         target_percent_identity=t.get('perc_id', ''))
                try:
                    info = api.get(base + '/lookup/id/' + quote(t['protein_id'], safe=''), **{'content-type': 'application/json'})
                    if info.get('species') != sp or info.get('object_type') != 'Translation':
                        raise ValueError('Target protein species/type mismatch')
                    transcript = api.get(base + '/lookup/id/' + quote(info['Parent'], safe=''), **{'content-type': 'application/json'})
                    if transcript.get('species') != sp or stable(transcript.get('Parent', '')) != stable(t['id']):
                        raise ValueError('Target protein does not belong to reported target gene')
                    response = api.get(base + '/sequence/id/' + quote(t['protein_id'], safe=''), type='protein', **{'content-type': 'application/json'})
                    if stable(response.get('id', '')) != stable(t['protein_id']):
                        raise ValueError('Sequence identifier mismatch')
                    seq = validate_sequence(response.get('seq', ''))
                    if info.get('length') is not None and len(seq) != int(info['length']):
                        raise ValueError('Sequence length differs from target translation lookup')
                    fid = sp + '|' + t['id'] + '|' + t['protein_id']
                    if fid in records and records[fid][1] != seq:
                        raise ValueError('Conflicting sequences for same target protein')
                    records[fid] = (sp, seq)
                    a.update(status='ACCEPTED', fasta_id=fid, sequence_length=len(seq),
                             sequence_sha256=hashlib.sha256(seq.encode()).hexdigest(),
                             target_assembly=info.get('assembly_name', ''), target_transcript=info['Parent'])
                except Exception as e:
                    a.update(status='SEQUENCE_OR_VALIDATION_ERROR', error=str(e))
            audit.append(a)
        # Persist audit after each successfully queried candidate.
        save_tables(out, mapping, audit)
    save_tables(out, mapping, audit)
    by_species = _layout_resource(out / 'by_species')
    by_candidate = _layout_resource(out / 'by_candidate')
    by_species.mkdir(exist_ok=True)
    by_candidate.mkdir(exist_ok=True)
    def fasta(fid):
        seq = records[fid][1]
        return '>' + fid + '\n' + '\n'.join(seq[i:i+60] for i in range(0, len(seq), 60)) + '\n'
    with open(_layout_resource(out / 'all_one2one_proteins.fasta'), 'w') as f:
        for fid in sorted(records):
            f.write(fasta(fid))
    for sp in species:
        (by_species / (sp + '.fasta')).write_text(''.join(fasta(fid) for fid in sorted(records) if records[fid][0] == sp))
    for inp in inputs:
        ids = sorted({a['fasta_id'] for a in audit if a['accession'] == inp['accession'] and a['status'] == 'ACCEPTED'})
        if ids:
            (by_candidate / (inp['accession'] + '.fasta')).write_text(''.join(fasta(fid) for fid in ids))
    meta.update(completed_utc=now(), pair_status_counts=dict(Counter(a['status'] for a in audit)),
                unique_target_sequences=len(records), pair_count=len(audit))
    assert len(audit) == len(inputs) * len(species)
    (_layout_resource(out / 'run_metadata.json')).write_text(json.dumps(meta, indent=2))
    print(json.dumps({k: meta[k] for k in ('pair_status_counts', 'unique_target_sequences')}, indent=2))
    errors = {'MAPPING_UNRESOLVED', 'HOMOLOGY_QUERY_ERROR', 'SEQUENCE_OR_VALIDATION_ERROR'}
    return 2 if any(a['status'] in errors for a in audit) else 0

def save_tables(out, mapping, audit):
    mf = ['accession', 'input_genes', 'resolved_accession', 'uniprot_gene', 'uniprot_entry_type',
          'source_gene', 'source_gene_name', 'source_assembly', 'mapping_status', 'mapping_note', 'error']
    af = ['accession', 'input_genes', 'resolved_accession', 'uniprot_gene', 'source_gene',
          'target_species', 'status', 'reported_types', 'orthology_type', 'is_high_confidence',
          'target_gene', 'target_transcript', 'target_protein', 'target_assembly',
          'source_percent_identity', 'target_percent_identity', 'fasta_id',
          'sequence_length', 'sequence_sha256', 'error']
    write_tsv(_layout_resource(out / 'source_mapping.tsv'), mapping, mf)
    write_tsv(_layout_resource(out / 'pair_audit.tsv'), audit, af)
    write_tsv(_layout_resource(out / 'accepted_manifest.tsv'), [a for a in audit if a['status'] == 'ACCEPTED'], af)

if __name__ == '__main__':
    sys.exit(main())
