"""Offline ELM regex inspection. Reference associations are not candidate functions."""
import csv
import io
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

AA = set('ACDEFGHIKLMNPQRSTVWY')
RESOURCE_SCHEMAS = {
    'elm_classes.tsv': {'Accession', 'ELMIdentifier', 'Regex', 'Probability', 'Description'},
    'elm_interaction_domains.tsv': {'ELM identifier', 'Interaction Domain Id'},
    'elm_interactions.tsv': {'Elm', 'Domain', 'interactorElm', 'interactorDomain', 'PMID'},
    'elm_instances.tsv': {'Accession', 'ELMIdentifier', 'Primary_Acc', 'Start', 'End', 'InstanceLogic'},
    'elm_kegg.csv': {'uniprot_id', 'elm', 'kegg_id', 'kegg_name'},
    'elm_reactome.csv': {'uniprot_id', 'elm', 'pathway_id', 'display_name'},
}


def read_table(path, required):
    """ELM .csv exports may actually be TSV; metadata precedes quoted headers."""
    lines = Path(path).read_text(encoding='utf-8-sig').splitlines(keepends=True)
    metadata, body = [], []
    for line in lines:
        if line.startswith('#'):
            metadata.append(line.strip())
        elif line.strip():
            body.append(line)
    if not body:
        raise ValueError(f'{path}: empty ELM table')
    delimiter = '\t' if '\t' in body[0] else ','
    reader = csv.DictReader(io.StringIO(''.join(body)), delimiter=delimiter)
    if not required <= set(reader.fieldnames or []):
        raise ValueError(f'{path}: missing ELM columns {sorted(required-set(reader.fieldnames or []))}')
    rows = list(reader)
    if any(None in row or any(row[k] is None for k in required) for row in rows):
        raise ValueError(f'{path}: malformed ELM row')
    return rows, metadata


def matches(regex, sequence):
    """Scan original full sequence: anchors/lookarounds retain biological positions.

    One regex-selected match length per start is retained (Python alternation and
    greediness semantics), including overlapping starts. No artificial peptide termini.
    """
    for start in range(len(sequence)):
        match = regex.match(sequence, start)
        if match:
            peptide = match.group()
            if peptide and set(peptide) <= AA:
                yield start+1, match.end(), peptide


def load_resources(directory, read_fasta):
    directory = Path(directory)
    tables, audit, paths = {}, [], []
    for name, schema in RESOURCE_SCHEMAS.items():
        path = directory/name
        if not path.exists():
            if name == 'elm_classes.tsv':
                raise ValueError(f'Required ELM class file missing: {path}')
            tables[name] = []
            audit.append(dict(resource=name, status='not_supplied', rows=0, detail='Optional reference table unavailable'))
            continue
        rows, metadata = read_table(path, schema)
        tables[name] = rows
        paths.append(path)
        for tag in ('Num_Classes', 'NumInstances'):
            declared = [int(m.split(':', 1)[1]) for m in metadata if m.startswith('#'+tag+':')]
            if declared and declared[0] != len(rows):
                raise ValueError(f'{path}: declared {tag} does not match row count')
        audit.append(dict(resource=name, status='loaded', rows=len(rows), detail='; '.join(metadata)))
    classes = {}
    for row in tables['elm_classes.tsv']:
        ident = row['ELMIdentifier']
        if ident in classes:
            raise ValueError(f'Duplicate ELM class {ident}')
        try:
            regex = re.compile(row['Regex'])
            probability = float(row['Probability'])
        except (re.error, ValueError) as exc:
            raise ValueError(f'Invalid ELM class {ident}: {exc}') from exc
        if regex.match('') or not 0 <= probability <= 1:
            raise ValueError(f'Empty-matching regex or invalid probability: {ident}')
        classes[ident] = row | {'compiled': regex}
    for name, key in [('elm_interaction_domains.tsv', 'ELM identifier'), ('elm_interactions.tsv','Elm'),
                      ('elm_instances.tsv','ELMIdentifier'), ('elm_kegg.csv','elm'), ('elm_reactome.csv','elm')]:
        missing = sorted({r[key] for r in tables[name]} - classes.keys())
        if missing:
            audit.append(dict(resource=name, status='unresolved_class_references', rows=len(missing),
                              detail=';'.join(missing)))
    refs = {}
    fasta = directory/'elm_instances.fasta'
    if fasta.exists():
        for record in read_fasta(fasta):
            parts = record['id'].split('|')
            accession = parts[1] if len(parts) >= 3 and parts[0] in ('sp','tr') else parts[0]
            if accession in refs and refs[accession] != record['sequence']:
                raise ValueError(f'Conflicting ELM reference sequences: {accession}')
            refs[accession] = record['sequence']
        paths.append(fasta)
    audit.append(dict(resource='elm_instances.fasta', status='loaded' if fasta.exists() else 'not_supplied',
                      rows=len(refs), detail='Reference protein sequences; not a complete instance catalogue'))
    return dict(classes=classes, tables=tables, audit=audit, paths=paths, references=refs)


def reference_audit(resources):
    rows, exact = [], defaultdict(list)
    for instance in resources['tables']['elm_instances.tsv']:
        name, accession = instance['ELMIdentifier'], instance['Primary_Acc']
        sequence = resources['references'].get(accession)
        start, end = int(instance['Start']), int(instance['End'])
        if start < 1 or end < start:
            raise ValueError('Invalid ELM instance coordinates: '+instance['Accession'])
        peptide, recovered = '', False
        if name not in resources['classes']:
            state = 'class_not_in_supplied_catalogue'
        elif sequence is None:
            state = 'reference_sequence_not_supplied'
        elif end > len(sequence):
            raise ValueError('ELM instance exceeds reference sequence: '+instance['Accession'])
        else:
            peptide = sequence[start-1:end]
            recovered = any(a==start and b==end for a,b,_ in matches(resources['classes'][name]['compiled'], sequence))
            state = 'exact_regex_span_recovered' if recovered else 'curated_span_differs_from_regex_match'
            # Only same full sequence, exact span, and explicit true-positive annotation
            # can be linked to a candidate. This is still imported evidence, not a new experiment.
            if recovered and instance['InstanceLogic'].strip().lower() == 'true positive':
                exact[sequence, name, start, end].append(instance['Accession'])
        rows.append(dict(instance=instance['Accession'], elm_class=name, accession=accession,
                         start=start, end=end, peptide=peptide, state=state,
                         instance_logic=instance['InstanceLogic'], methods=instance.get('Methods',''),
                         organism=instance.get('Organism',''), references=instance.get('References','')))
    return rows, exact


def load_focus(path, classes, families):
    data = json.loads(Path(path).read_text(encoding='utf-8'))
    for family, rules in data['families'].items():
        for rule in rules:
            if rule['elm_class'] not in classes:
                raise ValueError(f'Focus class absent from ELM catalogue: {rule["elm_class"]}')
            if rule.get('context', 'whole_sequence') not in ('whole_sequence','domain','c_terminal_window','n_terminal_window'):
                raise ValueError('Unsupported ELM focus context')
    return data


def focus_context(rule, a, b, length, regions):
    if rule is None:
        return 'not_family_prioritized'
    context = rule.get('context','whole_sequence')
    relevant = [(x,y) for name,x,y in regions if name in rule.get('regions',[])]
    if relevant:
        return 'region_supported' if any(x <= a <= b <= y for x,y in relevant) else 'outside_expected_region'
    if context == 'domain':
        return 'required_region_unavailable'
    if context == 'c_terminal_window':
        return 'terminal_window_only' if a > max(0,length-rule['window']) else 'outside_expected_window'
    if context == 'n_terminal_window':
        return 'terminal_window_only' if b <= rule['window'] else 'outside_expected_window'
    return 'sequence_only'


def run(resources, families, regions, focus, discovered_sites, out, write_tsv):
    classes, tables = resources['classes'], resources['tables']
    instance_rows, exact = reference_audit(resources)
    domain_map = defaultdict(set)
    for row in tables['elm_interaction_domains.tsv']:
        domain_map[row['ELM identifier']].add(row['Interaction Domain Id'])
    annotations = []
    for ident, entry in classes.items():
        interactions = [r for r in tables['elm_interactions.tsv'] if r['Elm']==ident]
        kegg = [r for r in tables['elm_kegg.csv'] if r['elm']==ident]
        reactome = [r for r in tables['elm_reactome.csv'] if r['elm']==ident]
        annotations.append(dict(elm_class=ident, accession=entry['Accession'], regex=entry['Regex'],
            name=entry.get('FunctionalSiteName',''), description=entry['Description'],
            elm_pattern_probability=entry['Probability'], probability_interpretation='ELM_class_metadata_not_a_hit_pvalue',
            partner_domains=';'.join(sorted(domain_map[ident])),
            supplied_interaction_rows=len(interactions),
            reference_pmids=';'.join(sorted({p for r in interactions for p in re.split('[,; ]+',r['PMID']) if p and p!='None'})),
            reference_kegg_ids=';'.join(sorted({r['kegg_id'] for r in kegg})),
            reference_reactome_ids=';'.join(sorted({r['pathway_id'] for r in reactome})),
            association_scope='reference_proteins_only_not_candidate_pathway_assignments',
            url='https://elm.eu.org/elms/'+ident))
    # Retain reference protein/pathway provenance without projecting it onto candidates.
    pathway_rows = []
    for filename, database, id_key, name_key in [('elm_kegg.csv','KEGG','kegg_id','kegg_name'),
                                               ('elm_reactome.csv','Reactome','pathway_id','display_name')]:
        for row in tables[filename]:
            pathway_rows.append(dict(database=database, reference_uniprot=row['uniprot_id'], elm_class=row['elm'],
                                     pathway_id=row[id_key], pathway_name=row[name_key],
                                     class_resolved=row['elm'] in classes))
    all_sequences = {r['sequence'] for rs in families.values() for r in rs}
    # Cache unique-sequence scans; same proteins can occur in several candidate families.
    scan = {}
    positive = defaultdict(set)
    for sequence in sorted(all_sequences):
        sequence_hits = []
        for ident, entry in classes.items():
            found = list(matches(entry['compiled'], sequence))
            if found:
                positive[ident].add(sequence)
                sequence_hits.extend((ident,a,b,peptide) for a,b,peptide in found)
        scan[sequence] = sequence_hits
    focus_map = {f: {r['elm_class']:r for r in focus['families'].get(f,[])} for f in families}
    fields = 'family sequence_id elm_class start end peptide family_priority context_state overlapping_domains disorder_context exact_reference_instances reference_partner_domains evidence'.split()
    focused, hit_count, reference_match_count = [], 0, 0
    with (out/'elm_hits.tsv').open('w',encoding='utf-8',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=fields,delimiter='\t')
        writer.writeheader()
        for family, records in families.items():
            for record in records:
                seq, ident = record['sequence'], record['id']
                for elm,a,b,peptide in scan[seq]:
                    rs = regions[family,ident]
                    overlaps=sorted({n for n,x,y in rs if n.startswith('PF') and a <= y and b >= x})
                    disorder=[(x,y) for n,x,y in rs if n=='disordered_region']
                    disorder_state=('within_supplied_disorder' if any(x <= a <= b <= y for x,y in disorder) else
                                    'not_within_supplied_disorder' if disorder else 'not_assessed')
                    rule = focus_map[family].get(elm)
                    references = exact.get((seq,elm,a,b),[])
                    reference_match_count += bool(references)
                    row=dict(family=family,sequence_id=ident,elm_class=elm,start=a,end=b,peptide=peptide,
                             family_priority=bool(rule),context_state=focus_context(rule,a,b,len(seq),rs),
                             overlapping_domains=';'.join(overlaps),disorder_context=disorder_state,
                             exact_reference_instances=';'.join(references),reference_partner_domains=';'.join(sorted(domain_map[elm])),
                             evidence='same_sequence_curated_reference' if references else 'regex_candidate_only')
                    writer.writerow(row)
                    hit_count += 1
                    if rule:
                        focused.append(row | {'review_note':rule['note']})
    prevalence, state_rows = [], []
    for family, records in families.items():
        target = {r['sequence'] for r in records}
        other = all_sequences-target
        for elm in classes:
            a,b = len(target & positive[elm]),len(other & positive[elm])
            prevalence.append(dict(family=family,elm_class=elm,family_priority=elm in focus_map[family],
                target_positive=a,target_total=len(target),other_positive=b,other_total=len(other),
                target_fraction=round(a/len(target),6),other_fraction=round(b/len(other),6) if other else '',
                prevalence_difference=round(a/len(target)-b/len(other),6) if other else '',
                interpretation='whole_sequence_descriptive_no_significance'))
        for record in records:
            for elm in focus_map[family]:
                found=[h for h in scan[record['sequence']] if h[0]==elm]
                state_rows.append(dict(family=family,sequence_id=record['id'],elm_class=elm,
                                       hits=len(found),state='regex_candidates_detected' if found else 'not_detected_not_family_exclusion'))
    # Annotate discovery at original sites. Never rescan an isolated word: that
    # would incorrectly satisfy ^/$ anchors and discard required flanks.
    discovery_rows=[]
    by_id={(f,r['id']):r['sequence'] for f,rs in families.items() for r in rs}
    for site in discovered_sites:
        seq=by_id[site['family'],site['sequence_id']]
        for elm,a,b,peptide in scan[seq]:
            if a <= site['end'] and b >= site['start']:
                relation=('elm_within_discovered_word' if site['start'] <= a <= b <= site['end'] else
                          'discovered_word_within_elm' if a <= site['start'] <= site['end'] <= b else 'partial_overlap')
                discovery_rows.append(dict(family=site['family'],sequence_id=site['sequence_id'],word=site['motif'],
                    word_start=site['start'],word_end=site['end'],elm_class=elm,elm_start=a,elm_end=b,
                    relation=relation,interpretation='positional_overlap_not_functional_validation'))
    # Interaction-domain roles refer to a potential binder, not a required domain
    # on the motif-bearing protein. Report availability without inventing edges.
    binders=[]
    for family,records in families.items():
        for record in records:
            domains={n for n,a,b in regions[family,record['id']] if n.startswith('PF')}
            for elm,required in domain_map.items():
                shared=domains & required
                if shared and elm in classes:
                    binders.append(dict(family=family,sequence_id=record['id'],elm_class=elm,
                        observed_partner_domains=';'.join(sorted(shared)),
                        interpretation='domain_compatibility_only_no_interaction_inferred'))
    outputs = [
        ('elm_resource_audit',resources['audit'],'resource status rows detail'.split()),
        ('elm_class_annotations',annotations,'elm_class accession regex name description elm_pattern_probability probability_interpretation partner_domains supplied_interaction_rows reference_pmids reference_kegg_ids reference_reactome_ids association_scope url'.split()),
        ('elm_reference_instances',instance_rows,'instance elm_class accession start end peptide state instance_logic methods organism references'.split()),
        ('elm_reference_pathways',pathway_rows,'database reference_uniprot elm_class pathway_id pathway_name class_resolved'.split()),
        ('elm_family_focus_hits',focused,fields+['review_note']),
        ('elm_inspection_states',state_rows,'family sequence_id elm_class hits state'.split()),
        ('elm_family_prevalence',prevalence,'family elm_class family_priority target_positive target_total other_positive other_total target_fraction other_fraction prevalence_difference interpretation'.split()),
        ('elm_discovery_annotation',discovery_rows,'family sequence_id word word_start word_end elm_class elm_start elm_end relation interpretation'.split()),
        ('elm_partner_domain_candidates',binders,'family sequence_id elm_class observed_partner_domains interpretation'.split())]
    for filename,rows,columns in outputs:
        write_tsv(out/(filename+'.tsv'),rows,columns)
    (out/'elm_family_focus.snapshot.json').write_text(json.dumps(focus,indent=2),encoding='utf-8')
    summary=dict(classes=len(classes),unique_sequences=len(all_sequences),regex_hit_rows=hit_count,
                 family_focus_hit_rows=len(focused),reference_instances=len(instance_rows),
                 exact_reference_spans_recovered=sum(r['state']=='exact_regex_span_recovered' for r in instance_rows),
                 candidate_reference_matches=reference_match_count,
                 discovery_overlap_rows=len(discovery_rows),
                 supplied_disorder=any(n=='disordered_region' for rs in regions.values() for n,a,b in rs))
    (out/'elm_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    return summary
