#!/usr/bin/env python3
"""Expand the version 1.0 paired panel to every crosswalk component using UniProt."""

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
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import shutil

from adhesome_core import digest, save_json, sequence_dict, table, write_fasta, write_table, VERSION
from adhesome_components import load_crosswalk
from build_reference_sets import UniProtClient, gene_names, resolve

ROOT = _layout_resource(Path(__file__)).resolve().parent
COMPONENT_FAMILY = 'component_reference'


def component_identity(entry, component):
    if str(entry.get('organism', {}).get('taxonId')) != '9606':
        return False
    names = {name.casefold() for name in gene_names(entry)}
    gene_ids = {r['id'] for r in entry.get('uniProtKBCrossReferences', []) if r['database'] == 'GeneID'}
    return component['Official Symbol'].casefold() in names or component['Gene ID'].strip() in gene_ids


def resolve_component(component, entries, client=None, repair=False):
    accessions = re.findall(r'\b[A-Z0-9]{6,10}(?:-\d+)?\b', component['reference_accession'])
    if not accessions:
        raise ValueError('No pinned accession: '+component['Official Symbol'])
    matches, rejected = {}, []
    for accession in accessions:
        entry = entries[accession]
        if component_identity(entry, component):
            matches[entry['primaryAccession']] = entry
        else:
            rejected.append(dict(accession=accession, resolved_accession=entry.get('primaryAccession', accession),
                genes=gene_names(entry), taxon_id=entry.get('organism', {}).get('taxonId'), reason='gene_or_taxon_mismatch'))
    if not matches and repair and client is not None:
        entry, method, alternatives = resolve(dict(accession='',gene=component['Official Symbol'],taxon_id='9606'),client)
        if not component_identity(entry,component):
            raise ValueError('Repaired component identity mismatch: '+component['Official Symbol'])
        matches[entry['primaryAccession']] = entry
    if len(matches) != 1:
        raise ValueError(f"{component['Official Symbol']}: expected one matching pinned reference, found {list(matches)}")
    entry = next(iter(matches.values()))
    if 'Swiss-Prot' not in entry.get('entryType', ''):
        raise ValueError('Component reference is not reviewed: '+component['Official Symbol'])
    sequence = entry.get('sequence', {}).get('value', '').upper()
    if not sequence or set(sequence)-set('ACDEFGHIKLMNPQRSTVWYXBZJUO'):
        raise ValueError('Invalid component sequence: '+component['Official Symbol'])
    return entry, sequence, rejected


def component_family_labels(entry, base_rows, rules):
    """Use explicit configured reference labels or an existing accession assignment.

    Domain overlap and functional category never select a family label.
    """
    names = {name.casefold() for name in gene_names(entry)}
    labels = {name for name, rule in rules.items() if names &
              {label.casefold() for label in rule.get('reference_labels', [])}}
    labels.update(row['family'] for row in base_rows if row['accession'] == entry['primaryAccession'])
    return sorted(name for name in labels if name in rules and not rules[name].get('screening_only')
                  and rules[name].get('classification_level') != 'screen')


def expand(args):
    base, out = _layout_resource(Path(args.base)).resolve(), _layout_resource(Path(args.out)).resolve()
    if _layout_resource(out).exists() or base == out or base in out.parents:
        raise ValueError('Choose a new output directory outside the base panel')
    lock = json.loads((_layout_resource(base/'reference_lock.json')).read_text(encoding='utf-8'))
    for name, expected in lock['files'].items():
        if digest(base/name) != expected:
            raise ValueError('Base reference lock mismatch: '+name)
    components = load_crosswalk(args.crosswalk)
    rules_path = getattr(args, 'rules', _layout_resource(ROOT/'config/family_rules.json'))
    rules = json.loads(_layout_resource(Path(rules_path)).read_text(encoding='utf-8'))['families']
    base_rows = [row for row in table(_layout_resource(base/'search_catalogue.csv')) if row.get('role') not in
                 ('crosswalk_component','crosswalk_family_reference')]
    client = UniProtClient(args.cache, getattr(args, 'offline', False))
    accessions = sorted({accession for row in components for accession in
        re.findall(r'\b[A-Z0-9]{6,10}(?:-\d+)?\b', row['reference_accession'])})
    entries = {}
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(client.get, 'https://rest.uniprot.org/uniprotkb/'+accession+'.json'):accession for accession in accessions}
        for future in as_completed(futures):
            accession = futures[future]
            entries[accession] = future.result()
            if len(entries) % 25 == 0 or len(entries) == len(accessions):
                print(f'Retrieved {len(entries)}/{len(accessions)} pinned accessions', flush=True)
    rows, records, coverage, provenance, family_rows, family_records = [], [], [], [], [], []
    for component in components:
        entry, sequence, rejected = resolve_component(component, entries, client, getattr(args, 'repair_identity', False))
        symbol, accession = component['Official Symbol'], entry['primaryAccession']
        rid = 'Hsa_component_'+symbol
        if not re.fullmatch(r'[A-Za-z0-9_.-]+', rid):
            raise ValueError('Unsafe component symbol: '+symbol)
        names = gene_names(entry)
        row = dict(reference_id=rid,family=COMPONENT_FAMILY,module='crosswalk_components',
            species=entry['organism']['scientificName'],taxon_id=9606,gene=';'.join(names),accession=accession,
            requested_accession=component['reference_accession'],
            protein_name=entry.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value',component['Protein name']),
            length=len(sequence),sequence_sha256=hashlib.sha256(sequence.encode()).hexdigest(),
            sequence_version=entry.get('entryAudit', {}).get('sequenceVersion',''),
            entry_version=entry.get('entryAudit', {}).get('entryVersion',''), reviewed=True,
            resolve_method='unique_gene_taxon_repair' if entry['primaryAccession'] not in re.findall(r'\b[A-Z0-9]{6,10}(?:-\d+)?\b',component['reference_accession']) else 'pinned_accession_gene_or_geneid_and_taxon',
            source_url='https://www.uniprot.org/uniprotkb/'+accession+'/entry',
            evidence_source='adhesome_component_crosswalk; UniProtKB:'+accession,
            experimental_evidence='Crosswalk membership assertion; no experimental adhesome verification inferred',
            evidence_status='crosswalk_component_identity_verified',adhesome_experiment_verified='false',
            experimental_annotation_count='',experimental_pubmed_ids='',role='crosswalk_component',
            reference_tier='component_reference',phylogenetic_group='Vertebrata',
            family_assignment_basis='No family assignment; component-specific human sequence reference',
            anchor_limitations='Similarity and panel-relative reciprocal agreement are hypotheses, not orthology assignments',
            component_symbol=symbol,component_gene_id=component['Gene ID'],
            component_functional_category=component['Functional Category'])
        rows.append(row)
        records.append((rid,sequence))
        labels = component_family_labels(entry, base_rows, rules)
        for family in labels:
            family_id = rid+'__'+family
            family_rows.append(dict(row,reference_id=family_id,family=family,module=rules[family]['module'],
                role='crosswalk_family_reference',reference_tier='panel_reference',
                family_assignment_basis='Explicit configured reference label or base-panel accession assignment'))
            family_records.append((family_id,sequence))
        coverage.append(dict(component_symbol=symbol,gene_id=component['Gene ID'],
            requested_accession=component['reference_accession'],resolved_accession=accession,
            resolved_gene=names[0] if names else '',reference_id=rid,status='covered',
            rejected_accessions=';'.join(r['accession'] for r in rejected),sequence_sha256=row['sequence_sha256'],
            family_labels=';'.join(labels)))
        provenance.append(dict(component_symbol=symbol,rejected=rejected,source_url=row['source_url'],
            entry_version=row['entry_version'],sequence_version=row['sequence_version'],sequence_sha256=row['sequence_sha256']))
    # Publish only after every component passes identity, taxon and sequence checks.
    shutil.copytree(base, out)
    write_table(_layout_resource(out/'component_catalogue.csv'), rows, list(rows[0]))
    write_fasta(_layout_resource(out/'component_references.faa'),records)
    write_table(_layout_resource(out/'component_family_catalogue.csv'), family_rows, list(rows[0]))
    write_fasta(_layout_resource(out/'component_family_references.faa'), family_records)
    write_table(_layout_resource(out/'component_coverage.csv'),coverage,list(coverage[0]))
    for catalogue_name, fasta_name in [('search_catalogue.csv','search_references.faa'),('metazoan_catalogue.csv','metazoan_references.faa')]:
        original_rows = [row for row in table(base/catalogue_name) if row.get('role') not in
                         ('crosswalk_component','crosswalk_family_reference')]
        sequences = sequence_dict(base/fasta_name)
        original_sequences = {row['reference_id']:sequences[row['reference_id']] for row in original_rows}
        if set(original_sequences) & {rid for rid,seq in records}:
            raise ValueError('Base already contains component references; select the pre-expansion base panel')
        fields = list(dict.fromkeys(key for row in original_rows+rows+family_rows for key in row))
        write_table(out/catalogue_name, original_rows+rows+family_rows, fields)
        write_fasta(out/fasta_name,list(original_sequences.items())+records+family_records)
    # Optional, explicitly targeted orthology panel; deduplicate by accession.
    search_rows = table(_layout_resource(out/'search_catalogue.csv'))
    search_sequences = sequence_dict(_layout_resource(out/'search_references.faa'))
    human = {}
    for row in search_rows:
        if str(row['taxon_id']) == '9606':
            human[row['accession']] = (row,search_sequences[row['reference_id']])
    write_fasta(_layout_resource(out/'HsaAdhesomeReferencePanel.faa'),[(accession+' gene='+row['gene'].split(';')[0],seq)
        for accession,(row,seq) in sorted(human.items())])
    save_json(_layout_resource(out/'HsaAdhesomeReferencePanel.reference.json'),dict(pipeline_version=VERSION,taxon_id=9606,
        species='Homo sapiens',complete_proteome=False,panel_scope='targeted_adhesome_crosswalk_and_retained_human_references',
        sequences=len(human),fasta_sha256=digest(_layout_resource(out/'HsaAdhesomeReferencePanel.faa')),
        crosswalk_sha256=digest(args.crosswalk),absence_interpretation='Unrepresented proteins are not evidence of biological absence'))
    save_json(_layout_resource(out/'component_expansion.json'), dict(pipeline_version=VERSION,created_at=datetime.now(timezone.utc).isoformat(),
        crosswalk_sha256=digest(args.crosswalk),rules_sha256=digest(rules_path),base_lock_sha256=digest(_layout_resource(base/'reference_lock.json')),
        builder_sha256=digest(__file__),components=len(components),covered=len(rows),
        forward_references=len(search_rows),family_labelled_component_references=len(family_rows),
        reciprocal_references=len(table(_layout_resource(out/'metazoan_catalogue.csv'))),
        orthology_panel=dict(file='HsaAdhesomeReferencePanel.faa',sequences=len(human),complete_proteome=False),
        references=provenance,requests=sorted(client.requests,key=lambda r:r['url'])))
    lock.update(pipeline_version=VERSION,created_at=datetime.now(timezone.utc).isoformat(),
        crosswalk_sha256=digest(args.crosswalk),component_expansion='component_expansion.json',
        interpretation='Family and component reference panel with retained competitors; component references do not establish family identity or experimentally verified adhesome membership.')
    for name in ('component_catalogue.csv','component_references.faa','component_family_catalogue.csv','component_family_references.faa','component_coverage.csv','component_expansion.json',
                 'search_catalogue.csv','search_references.faa','metazoan_catalogue.csv','metazoan_references.faa',
                 'HsaAdhesomeReferencePanel.faa','HsaAdhesomeReferencePanel.reference.json'):
        lock['files'][name] = digest(out/name)
    save_json(_layout_resource(out/'reference_lock.json'),lock)
    print(f'Published {len(rows)}/{len(components)} component references; {len(search_rows)} forward references: {out}',flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--version',action='version',version='adhesome-component-panel '+VERSION)
    parser.add_argument('--base',type=Path,default=_layout_resource(ROOT/'references/adhesome_v1.0'))
    parser.add_argument('--crosswalk',type=Path,default=_layout_resource(ROOT/'config/adhesome_component_crosswalk.csv'))
    parser.add_argument('--rules',type=Path,default=_layout_resource(ROOT/'config/family_rules.json'))
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--cache',type=Path,default=_layout_resource(ROOT/'references/cache/uniprot'))
    parser.add_argument('--offline',action='store_true')
    parser.add_argument('--repair-identity',action='store_true',help='Repair obsolete/mismatched pins only through a unique gene/taxon match, recording rejected accessions')
    parser.add_argument('--workers',type=int,default=4)
    args = parser.parse_args()
    if not 1 <= args.workers <= 8:
        parser.error('--workers must be between 1 and 8')
    expand(args)


if __name__ == '__main__':
    main()
