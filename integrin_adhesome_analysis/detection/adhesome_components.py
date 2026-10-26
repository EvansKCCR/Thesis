"""Crosswalk annotations retain possible components without claiming orthology."""
from collections import defaultdict
import re
from adhesome_core import table, write_table, digest, save_json, VERSION

FIELDS = ['protein_id', 'component_symbol', 'gene_id', 'protein_name', 'swiss_prot_id',
          'reference_accession', 'functional_category', 'adhesome_category', 'matched_pfams',
          'reference_pfams', 'signature_complete', 'reference_similarity', 'reference_reciprocal', 'interpretation']


def load_crosswalk(path, manifest=None):
    rows=table(path)
    required={'Official Symbol','Gene ID','Protein name','Swiss-Prot ID','Functional Category',
              'FA','reference_accession','reference_pfams'}
    if not rows or not required<=rows[0].keys():
        raise ValueError('Invalid adhesome component crosswalk columns')
    known={a for r in table(manifest) for a in re.findall(r'PF\d{5}',r['Pfam accession'])} if manifest else None
    seen=set()
    for row in rows:
        symbol=row['Official Symbol'].strip()
        if not symbol or symbol in seen:
            raise ValueError('Missing/duplicate component symbol: '+symbol)
        seen.add(symbol)
        text=row['reference_pfams'].strip()
        if text and not re.fullmatch(r'PF\d{5}(;PF\d{5})*',text):
            raise ValueError('Invalid reference_pfams for '+symbol)
        domains=set(text.split(';')) if text else set()
        if known is not None and domains-known:
            raise ValueError('Crosswalk domains absent from manifest: '+', '.join(sorted(domains-known)))
        row['_pfams']=domains
    return rows


def annotate_components(path, manifest, domains, evidence, candidates, output):
    components=load_crosswalk(path,manifest)
    observed=defaultdict(set)
    for row in domains:
        if row.get('decision') in (None,'','retained'):
            observed[row['query_name']].add(row['pfam_accession'])
    references=defaultdict(set)
    reference_components=defaultdict(set)
    reciprocal_components=defaultdict(set)
    for (protein,family),support in evidence.items():
        for row in support:
            if row.get('similarity','').lower()=='true' and row.get('reference_accession'):
                references[protein].add(row['reference_accession'])
                if row.get('component_symbol'):
                    reference_components[protein].add(row['component_symbol'])
                    if family=='component_reference' and row.get('reciprocal','').lower() in ('true','false'):
                        reciprocal_components[(protein,row['component_symbol'])].add(row['reciprocal'].lower())
    matches=[]
    by_protein=defaultdict(list)
    for protein in sorted(set(observed)|set(references)):
        for component in components:
            overlap=observed[protein]&component['_pfams']
            similar=bool(set(re.split(r'[,;\s]+',component['reference_accession'])) & references[protein]) or component['Official Symbol'] in reference_components[protein]
            reciprocal=reciprocal_components[(protein,component['Official Symbol'])]
            if not overlap and not similar:
                continue
            match=dict(protein_id=protein,component_symbol=component['Official Symbol'],
                gene_id=component['Gene ID'],protein_name=component['Protein name'],
                swiss_prot_id=component['Swiss-Prot ID'],reference_accession=component['reference_accession'],
                functional_category=component['Functional Category'],adhesome_category=component['FA'],
                matched_pfams=';'.join(sorted(overlap)),reference_pfams=component['reference_pfams'],
                signature_complete=bool(component['_pfams']) and component['_pfams']<=observed[protein],
                reference_similarity=similar,
                reference_reciprocal='conflict' if len(reciprocal)>1 else 'supported' if reciprocal=={'true'} else
                                     'not_supported' if reciprocal=={'false'} else 'not_tested',
                interpretation='Reference similarity supports a component hypothesis; orthology untested' if similar else
                               'Shared Pfam domains only; component identity and orthology untested')
            matches.append(match)
            by_protein[protein].append(match)
    for row in candidates:
        possible=by_protein[row['protein_id']]
        row['crosswalk_possible_components']=';'.join(sorted({r['component_symbol'] for r in possible}))
        row['crosswalk_reference_components']=';'.join(sorted({r['component_symbol'] for r in possible if r['reference_similarity']}))
        row['crosswalk_reciprocal_components']=';'.join(sorted({r['component_symbol'] for r in possible if r['reference_reciprocal']=='supported'}))
        row['crosswalk_functional_categories']=';'.join(sorted({r['functional_category'] for r in possible}))
    write_table(output/'component_hypotheses.csv',matches,FIELDS)
    save_json(output/'component_crosswalk.json',dict(pipeline_version=VERSION,
        crosswalk_sha256=digest(path),components=len(components),matches=len(matches),
        components_without_pfam=[r['Official Symbol'] for r in components if not r['_pfams']],
        interpretation='Crosswalk overlap does not promote detection status or establish family identity'))


CANDIDATE_FIELDS=['crosswalk_possible_components','crosswalk_reference_components','crosswalk_reciprocal_components','crosswalk_functional_categories']
