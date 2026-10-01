"""Master-workbook adapter; supporting sheets never recruit excluded records."""
import csv
import hashlib
import json
import re
import warnings
from collections import Counter
from pathlib import Path
import curated_annotation as ca


def table(sheets,name):
    rows=sheets.get(name,[])
    if not rows: return []
    header=rows[0]['values']
    return [dict(zip(header,r['values']+['']*max(0,len(header)-len(r['values']))))|{'source_row':r['row']}
            for r in rows[1:] if any(v is not None and v!='' for v in r['values'])]


def load(path,rules):
    sheets=ca.read_workbook(path)
    if 'Master_motif_library' not in sheets: raise ValueError('Master_motif_library sheet required')
    rows=table(sheets,'Master_motif_library'); adhesion={r['motif_id']:r for r in table(sheets,'Adhesion_curated')}
    seen=set(); curated=[]; generic=[]; audit=[]
    exclusions=table(sheets,'QA_exclusions')
    excluded={str(r[k]) for r in exclusions for k in ('library_id','motif_id','source_identifier','elm_accession','elm_identifier') if r.get(k)}
    for raw in rows:
        r={k:('' if v is None else v) for k,v in raw.items()}
        required={'library_id','source_collection','family','scan_enabled','regex_valid','regex_python','screening_recommendation','motif_type'}
        if not required<=r.keys(): raise ValueError('Master schema missing required columns')
        mid=r['library_id']
        if not mid or mid in seen: raise ValueError('Duplicate/missing master ID: '+str(mid))
        seen.add(mid)
        if r['scan_enabled'] not in ('Yes','No') or r['regex_valid'] not in ('Yes','No'): raise ValueError('Invalid scan/regex flag: '+mid)
        reader=r['screening_recommendation']=='Annotation only' or r['motif_type']=='motif-reader annotation'
        forbidden=bool({str(r.get(k,'')) for k in ('library_id','source_identifier','elm_accession','elm_identifier')} & excluded)
        pattern=r['regex_python']; compiled=None; issue=''
        if pattern:
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter('error',FutureWarning)
                    compiled=re.compile(pattern)
                if compiled.match(''): compiled=None;issue='empty_matching_regex'
            except (re.error,FutureWarning): issue='invalid_or_ambiguous_regex';compiled=None
        else: issue='missing_pattern'
        enabled=r['scan_enabled']=='Yes' and r['regex_valid']=='Yes' and compiled is not None and not reader and not forbidden
        state='excluded_by_QA' if forbidden else 'reader_annotation_only' if reader else 'scan_enabled' if enabled else 'not_scanned'
        if r['regex_valid']=='Yes' and compiled is None and not reader: issue=issue or 'regex_flag_conflict'
        audit.append(dict(library_id=mid,source_collection=r['source_collection'],source_row=r['source_row'],
                          declared_scan_enabled=r['scan_enabled'],effective_state=state,issue=issue))
        r.update(compiled=compiled,enabled=enabled,reader_only=reader,effective_state=state)
        if r['source_collection']=='Curated adhesion SLiM library':
            if mid not in rules['motifs']: raise ValueError('No reviewed adhesion context rule for '+mid)
            rule=rules['motifs'][mid]
            if rule.get('source_family')!=r['family']: raise ValueError('Review changed family routing: '+mid)
            if reader!=(rule['role']=='reader_annotation'): raise ValueError('Reader routing mismatch: '+mid)
            if forbidden or (not enabled and not reader): continue
            source=adhesion.get(mid,{})
            # Citation status is joined only from the exact adhesion ID, never inferred from score.
            status=source.get('citation_status') or 'unresolved'
            pattern=r['regex_python']; normalization='none'
            if mid in ('SLIM012','SLIM013') and '[OP]' in pattern:
                pattern=pattern.replace('[OP]','P');normalization='encoded_proline_precursor_[OP]_to_P'
            score=r.get('confidence_score')
            if not isinstance(score,(int,float)) or not 0<=score<=100: raise ValueError('Invalid adhesion score: '+mid)
            curated.append(r|dict(motif_id=mid,source_sheet='Master_motif_library',scan_pattern=pattern,
                compiled=re.compile(pattern) if not reader else None,normalization=normalization,citation_status=status,
                primary_literature_citation=r.get('primary_citation',''),primary_literature_doi_or_pmid=r.get('primary_doi_or_pmid','')))
        else: generic.append(r)
    summary={r['metric']:r['value'] for r in table(sheets,'Library_summary')}
    checks={'Total non-redundant master records':len(rows),'Scan-enabled records':sum(r['scan_enabled']=='Yes' for r in rows),
            'Regex-valid records':sum(r['regex_valid']=='Yes' for r in rows)}
    discrepancies=[]
    for key,value in checks.items():
        if key in summary and summary[key]!=value:
            if (key=='Total non-redundant master records' and summary[key]==value+3
                    and summary.get('Curated adhesion records')==37
                    and summary.get('Current ELM class records')==352
                    and summary.get('Legacy candidates retained for manual review')==302):
                discrepancies.append(f'{key}: reported {summary[key]}, actual nonblank {value}; source summary retains three extra blank-row counts')
            else: raise ValueError('Workbook summary mismatch: '+key)
    metadata=dict(master_records=len(rows),collections=dict(Counter(r['source_collection'] for r in rows)),
        declared_scan_enabled=checks['Scan-enabled records'],summary_discrepancies=discrepancies,qa_sheet_present='QA_exclusions' in sheets,
        reported_exclusions=summary.get('Excluded/flagged source records',''),
        warning='' if 'QA_exclusions' in sheets else 'QA_exclusions sheet missing; excluded-record identities cannot be audited. Master remains authoritative.')
    return curated,sheets,[],generic,audit,metadata


def finish(generic,audit,metadata,sheets,families,regions,out,write_tsv,sequence_context=None):
    sequence_context=sequence_context or {}
    out=Path(out)
    # Retain detailed adhesion evidence, and provide one compact table for all sequence hits.
    (out/'motif_annotations.tsv').rename(out/'adhesion_motif_annotations.tsv')
    fields='family sequence_id library_id source_collection start end peptide site_group_id candidate_tier overlapping_domains predicted_topology_segment deeploc_localizations signalp_prediction context_requirements functional_hypothesis source_row'.split()
    def site_topology(family,ident,a,b):
        names={'predicted_inside','predicted_outside','transmembrane','signal_peptide','cytoplasmic_tail'}
        return ';'.join(sorted({n for n,x,y in regions[family,ident] if n in names and x<=a<=b<=y}))
    count=0; positive=Counter()
    with (out/'motif_annotations.tsv').open('w',encoding='utf-8',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=fields,delimiter='\t',extrasaction='ignore');writer.writeheader()
        with (out/'adhesion_motif_annotations.tsv').open(encoding='utf-8',newline='') as source:
            for row in csv.DictReader(source,delimiter='\t'):
                family,ident=row['family'],row['sequence_id']
                context=sequence_context.get((family,ident),{})
                writer.writerow(row|{'library_id':row['motif_id'],'source_collection':'Curated adhesion SLiM library',
                    'predicted_topology_segment':site_topology(family,ident,int(row['start']),int(row['end'])),
                    'deeploc_localizations':context.get('deeploc_localizations',''),
                    'signalp_prediction':context.get('signalp_prediction','')})
                count+=1
        seqs={r['sequence'] for rs in families.values() for r in rs}
        for spec in generic:
            if not spec['enabled']: continue
            cached={seq:list(ca.occurrences(spec['compiled'],seq)) for seq in seqs}
            tier='review_only' if spec['screening_recommendation']=='Review only' or spec['source_collection']=='Legacy ELM candidate library' else 'ELM_context_unverified'
            for family,records in families.items():
                for record in records:
                    for a,b,peptide in cached[record['sequence']]:
                        domains=sorted({n for n,x,y in regions[family,record['id']] if (n.startswith('PF') or n.startswith('CDD:')) and a<=y and b>=x})
                        group=hashlib.sha256(f'{family}|{record["id"]}|{a}|{b}|{peptide}'.encode()).hexdigest()[:16]
                        writer.writerow(dict(family=family,sequence_id=record['id'],library_id=spec['library_id'],
                            source_collection=spec['source_collection'],start=a,end=b,peptide=peptide,site_group_id=group,
                            candidate_tier=tier,overlapping_domains=';'.join(domains),context_requirements=spec.get('context_requirements',''),
                            predicted_topology_segment=site_topology(family,record['id'],a,b),
                            deeploc_localizations=sequence_context.get((family,record['id']),{}).get('deeploc_localizations',''),
                            signalp_prediction=sequence_context.get((family,record['id']),{}).get('signalp_prediction',''),
                            functional_hypothesis='Unvalidated sequence candidate; source-class function: '+spec.get('molecular_function',''),source_row=spec['source_row']))
                        count+=1;positive[spec['library_id']]+=1
    write_tsv(out/'master_library_audit.tsv',audit,'library_id source_collection source_row declared_scan_enabled effective_state issue'.split())
    write_tsv(out/'legacy_review_catalogue.tsv',[{k:v for k,v in r.items() if k!='compiled'} for r in generic if r['source_collection']=='Legacy ELM candidate library'],
              'library_id source_identifier motif_name regex_python scan_enabled regex_valid effective_state status context_requirements molecular_function source_row'.split())
    pathways=[]
    identifiers={r.get('elm_identifier') for r in table(sheets,'Master_motif_library')}
    for sheet in ('ELM_KEGG','ELM_Reactome'):
        for row in table(sheets,sheet):
            pathways.append(row|{'database':sheet,'class_in_master':row.get('ELMIdentifier') in identifiers,
                                  'interpretation':'reference_class_association_not_candidate_pathway_assignment'})
    write_tsv(out/'pathway_reference_associations.tsv',pathways,'database ELMIdentifier pathway_id pathway_name class_in_master interpretation'.split())
    write_tsv(out/'source_workbook_provenance.tsv',[r|{'verification':'recorded_source_checksum_not_independently_verified'} for r in table(sheets,'Provenance')],
              'file sha256 size_bytes verification'.split())
    write_tsv(out/'generic_class_counts.tsv',[dict(library_id=r['library_id'],state=r['effective_state'],hits=positive[r['library_id']]) for r in generic],
              'library_id state hits'.split())
    metadata.update(sequence_annotation_rows=count,generic_hit_rows=sum(positive.values()))
    (out/'master_summary.json').write_text(json.dumps(metadata,indent=2),encoding='utf-8')
    return metadata
