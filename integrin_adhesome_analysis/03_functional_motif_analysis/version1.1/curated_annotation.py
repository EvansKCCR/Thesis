"""Workbook-driven motif hypotheses with explicit, non-probabilistic evidence gates."""
import csv
import hashlib
import json
import posixpath
import re
import zipfile
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from pathlib import Path

AA = set('ACDEFGHIKLMNPQRSTVWY')
NS = {'m':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
REQUIRED = {'motif_id','family','motif_name','regex_python','motif_type','motif_location',
            'molecular_function','binding_partner_or_reader','evidence_tier','screening_recommendation',
            'context_requirements','confidence_score','confidence_class','citation_status','source_url',
            'primary_literature_citation','primary_literature_doi_or_pmid','notes'}


def read_workbook(path):
    """Read static XLSX values without changing the curated workbook or dependencies."""
    sheets = {}
    with zipfile.ZipFile(path) as z:
        shared=[]
        if 'xl/sharedStrings.xml' in z.namelist():
            for si in ET.fromstring(z.read('xl/sharedStrings.xml')).findall('m:si',NS):
                shared.append(''.join(t.text or '' for t in si.findall('.//m:t',NS)))
        rels={r.attrib['Id']:r.attrib['Target'] for r in ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))}
        for sheet in ET.fromstring(z.read('xl/workbook.xml')).findall('m:sheets/m:sheet',NS):
            target=rels[sheet.attrib['{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id']]
            member=target.lstrip('/') if target.startswith('/') else posixpath.normpath('xl/'+target)
            rows=[]
            for row in ET.fromstring(z.read(member)).findall('m:sheetData/m:row',NS):
                values={}
                for cell in row.findall('m:c',NS):
                    if cell.find('m:f',NS) is not None:
                        raise ValueError(f'Library formulas are not accepted: {sheet.attrib["name"]}!{cell.attrib["r"]}; supply evaluated static values')
                    col=0
                    for char in re.match('[A-Z]+',cell.attrib['r']).group():
                        col=col*26+ord(char)-64
                    kind=cell.attrib.get('t','n')
                    raw=cell.find('m:v',NS)
                    value=raw.text if raw is not None else None
                    if kind=='s' and value is not None:
                        value=shared[int(value)]
                    elif kind=='inlineStr':
                        value=''.join(t.text or '' for t in cell.findall('.//m:t',NS))
                    elif kind=='n' and value is not None:
                        value=float(value)
                        if value.is_integer(): value=int(value)
                    elif kind=='e':
                        raise ValueError(f'Excel error at {sheet.attrib["name"]}!{cell.attrib["r"]}')
                    values[col-1]=value
                rows.append({'row':int(row.attrib['r']),'values':[values.get(i) for i in range(max(values,default=-1)+1)]})
            sheets[sheet.attrib['name']]=rows
    return sheets


def load_library(path, rules):
    sheets=read_workbook(path)
    if 'Motif_library' not in sheets or not sheets['Motif_library']:
        raise ValueError('Workbook requires a Motif_library worksheet')
    header=sheets['Motif_library'][0]['values']
    if not REQUIRED <= set(header):
        raise ValueError('Missing workbook columns: '+str(sorted(REQUIRED-set(header))))
    if len(header)!=len(set(header)):
        raise ValueError('Duplicate workbook column')
    motifs, audit=[],[]
    seen=set()
    for source in sheets['Motif_library'][1:]:
        values=source['values']
        if not any(v is not None and v!='' for v in values): continue
        row={key:values[i] if i<len(values) and values[i] is not None else '' for i,key in enumerate(header)}
        ident=row['motif_id']
        if not ident or ident in seen: raise ValueError('Missing/duplicate motif ID: '+str(ident))
        seen.add(ident)
        if ident not in rules['motifs']: raise ValueError(f'No reviewed operational context rule for {ident}')
        rule=rules['motifs'][ident]
        if rule.get('source_family',row['family'])!=row['family']:
            raise ValueError('Workbook family changed; review context routing for '+ident)
        if row['screening_recommendation'] not in ('Include','Include with structural filter','Include with positional filter','Review only','Annotation only'):
            raise ValueError('Unknown screening recommendation: '+ident)
        if row['evidence_tier'] not in ('A','B','C') or row['citation_status'] not in ('direct','supporting','indirect','unresolved'):
            raise ValueError(f'Unknown evidence tier or citation status: {ident}')
        score=row['confidence_score']
        if not isinstance(score,(float,int)) or not 0<=score<=100: raise ValueError('Invalid library score: '+ident)
        reader=row['screening_recommendation']=='Annotation only' or row['motif_type']=='motif-reader annotation'
        if reader != (rule['role']=='reader_annotation'): raise ValueError('Reader/encoded role mismatch: '+ident)
        pattern=row['regex_python']
        transform='none'
        if ident in ('SLIM012','SLIM013') and '[OP]' in pattern:
            pattern=pattern.replace('[OP]','P')
            transform='encoded_proline_precursor_[OP]_to_P; O_is_not_hydroxyproline_in_FASTA'
        if not reader:
            if not pattern: raise ValueError('Scannable motif lacks regex: '+ident)
            try: compiled=re.compile(pattern)
            except re.error as exc: raise ValueError(f'{ident}: {exc}') from exc
            if compiled.match(''): raise ValueError('Empty matching motif: '+ident)
        else:
            compiled=None
        motifs.append(row | {'source_sheet':'Motif_library','source_row':source['row'],
                             'scan_pattern':pattern,'normalization':transform,'reader_only':reader,'compiled':compiled})
        issues=[]
        if row['citation_status'] in ('indirect','unresolved'): issues.append('citation_not_directly_supportive')
        if row['screening_recommendation']=='Review only': issues.append('cannot_be_high_confidence')
        if reader: issues.append('reader_property_not_scanned')
        if rule['role'] in ('partner_motif','partner_screen'): issues.append('routed_to_motif_bearing_partner')
        if ident=='SLIM022': issues.append('approximate_consensus_requires_structural_validation')
        if ident=='SLIM035': issues.append('source_notes_mention_alpha_actinin_in_filamin_entry; retained_verbatim_for_review')
        audit.append(dict(motif_id=ident,source_row=source['row'],source_family=row['family'],
                          target_families=';'.join(rule['target_families']),role=rule['role'],
                          original_regex=row['regex_python'],scan_regex=pattern,normalization=transform,
                          library_confidence=score,issues=';'.join(issues) or 'none'))
    return motifs,sheets,audit


def occurrences(regex,seq):
    # match(pos) preserves full-sequence terminal anchors and original group numbering.
    for i in range(len(seq)):
        m=regex.match(seq,i)
        if m and m.end()>i and set(m.group())<=AA:
            yield i+1,m.end(),m.group()


def load_evidence(path,families,rules):
    evidence=defaultdict(list)
    if not path: return evidence
    seqs={(f,r['id']):r['sequence'] for f,rs in families.items() for r in rs}
    with Path(path).open(encoding='utf-8-sig',newline='') as handle:
        reader=csv.DictReader(handle,delimiter='\t')
        required={'family','sequence_id','sequence_sha256','motif_id','start','end','criterion','status','source'}
        if not required<=set(reader.fieldnames or []): raise ValueError('Evidence TSV columns required: '+','.join(sorted(required)))
        for row in reader:
            key=row['family'],row['sequence_id']
            if key not in seqs or row['sequence_sha256']!=hashlib.sha256(seqs[key].encode()).hexdigest():
                raise ValueError('Unknown sequence or stale evidence hash: '+str(key))
            if row['motif_id'] not in rules['motifs'] or row['criterion'] not in set(rules['motifs'][row['motif_id']]['requirements'])|{'localization','domain_overlap_resolved','family_architecture'}:
                raise ValueError('Unknown evidence motif/criterion')
            if row['status'] not in ('supported','contradicted') or not row['source'].strip():
                raise ValueError('Evidence requires supported/contradicted and a source')
            a,b=int(row['start']),int(row['end'])
            if not 1<=a<=b<=len(seqs[key]): raise ValueError('Invalid evidence interval')
            evidence[key,row['motif_id']].append(row|{'start':a,'end':b})
    return evidence


def import_topology(paths,families,regions):
    seqs={(f,r['id']):r['sequence'] for f,rs in families.items() for r in rs}
    predictions={}
    for path in paths:
        if Path(path).suffix.lower()=='.zip':
            with zipfile.ZipFile(path) as archive:
                members=[n for n in archive.namelist() if n.endswith('predicted_topologies.3line')]
                if not members: raise ValueError(f'No topology records in {path}')
                sources=[(str(Path(path).resolve())+'!'+n,archive.read(n).decode('utf-8')) for n in members]
        else:
            sources=[(str(Path(path).resolve()),Path(path).read_text())]
        for source,text in sources:
            lines=[l.strip() for l in text.splitlines() if l.strip()]
            if len(lines)%3: raise ValueError(f'Invalid 3line topology file: {source}')
            for i in range(0,len(lines),3):
                header,seq,top=lines[i:i+3]
                if not header.startswith('>') or len(seq)!=len(top) or set(top)-set('IOMSBX'):
                    raise ValueError(f'Invalid topology record in {source}')
                key=header[1:].split()[0]
                value=(seq,top)
                if key in predictions and predictions[key][:2]!=value: raise ValueError('Conflicting topology: '+key)
                predictions[key]=(seq,top,source)
    audit=[]
    for (family,ident),seq in seqs.items():
        pred=predictions.get(ident)
        state='missing' if pred is None else 'sequence_mismatch' if pred[0]!=seq else 'exact_sequence_verified'
        audit.append(dict(family=family,sequence_id=ident,status=state,source=pred[2] if pred else ''))
        if state!='exact_sequence_verified': continue
        top=pred[1]
        for m in re.finditer(r'(I+|O+|M+|S+|B+|X+)',top):
            label={'I':'predicted_inside','O':'predicted_outside','M':'transmembrane','B':'transmembrane','S':'signal_peptide','X':'unknown_topology'}[m.group()[0]]
            regions[family,ident].append((label,m.start()+1,m.end()))
        # Tail only for a single membrane-spanning helix followed by terminal inside residues.
        tm=list(re.finditer('M+',top))
        tail=re.search('I+$',top)
        if len(tm)==1 and tail and tm[0].end()==tail.start():
            regions[family,ident].append(('cytoplasmic_tail',tail.start()+1,len(seq)))
    return audit


def localization(rule,a,b,length,rs):
    def inside(names): return any(n in names and x<=a<=b<=y for n,x,y in rs)
    mode=rule['location']
    if mode=='extreme_c_terminus': return 'supported' if b==length else 'contradicted'
    if mode=='c_terminal_window': return 'provisional' if a>max(0,length-50) else 'contradicted'
    if mode=='cytoplasmic_tail':
        if inside({'cytoplasmic_tail'}): return 'supported'
        if any(n=='cytoplasmic_tail' for n,x,y in rs): return 'contradicted'
        return 'provisional' if a>max(0,length-100) else 'unknown'
    if mode=='extracellular':
        # Outside prediction includes organelle lumen, so is never treated as
        # definitive extracellular localization without compartment evidence.
        if inside({'extracellular_region'}): return 'supported'
        if inside({'predicted_outside'}): return 'provisional'
        if inside({'predicted_inside','transmembrane','signal_peptide'}): return 'contradicted'
        return 'unknown'
    if mode=='intracellular':
        if inside({'cytoplasmic_region','cytoplasmic_tail'}): return 'supported'
        if inside({'predicted_inside'}): return 'provisional'
        if inside({'predicted_outside','transmembrane','signal_peptide'}): return 'contradicted'
        return 'unknown'
    raise ValueError('Unknown localization mode: '+mode)


def assess(motif,rule,family,record,a,b,rs,domain_state,evidence):
    criteria={key:'unknown' for key in rule['requirements']}
    criteria['family_architecture']='supported' if domain_state=='domain_count_rule_supported' else 'unknown'
    criteria['localization']=localization(rule,a,b,len(record['sequence']),rs)
    contained={n for n,x,y in rs if x<=a<=b<=y}
    if 'PF00041' in contained and 'fnIII_context' in criteria: criteria['fnIII_context']='supported'
    if 'PF01391' in contained and 'collagen_context' in criteria: criteria['collagen_context']='supported'
    if 'disordered_region' in contained and 'disordered_region' in criteria:
        criteria['disordered_region']='supported'
    overlapping={n for n,x,y in rs if n.startswith('PF') and a<=y and b>=x}
    if overlapping and not rule['allow_domain_overlap']: criteria['domain_overlap_resolved']='unknown'
    applicable=[r for r in evidence.get(((family,record['id']),motif['motif_id']),[]) if r['start']<=a<=b<=r['end']]
    for key in set(r['criterion'] for r in applicable):
        values={r['status'] for r in applicable if r['criterion']==key}
        if 'contradicted' in values: criteria[key]='contradicted'
        elif criteria.get(key)!='contradicted': criteria[key]='supported'
    missing=[k for k,v in criteria.items() if v!='supported']
    conflicts=[k for k,v in criteria.items() if v=='contradicted']
    permitted=(motif['confidence_score']>=75 and motif['evidence_tier'] in ('A','B') and
               motif['citation_status'] in ('direct','supporting') and motif['screening_recommendation'].startswith('Include'))
    if conflicts: tier='context_conflict'
    elif not permitted: tier='review_only'
    elif not missing: tier='high_confidence_hypothesis'
    elif criteria['family_architecture']=='supported' and criteria['localization'] in ('supported','provisional'):
        tier='context_supported_candidate'
    else: tier='sequence_only_candidate'
    return tier,criteria,missing,sorted(overlapping),applicable


def run(motifs,sheets,audit,rules,families,regions,domain_rows,evidence,sites,out,write_tsv):
    ds={(r['family'],r['sequence_id']):r['state'] for r in domain_rows}
    hits,states,readers,coverage,prevalence=[],[],[],[],[]
    allseq={r['sequence'] for rs in families.values() for r in rs}
    cache={}
    for motif in motifs:
        mid=motif['motif_id']; rule=rules['motifs'][mid]
        targets=[f for f in rule['target_families'] if f in families]
        coverage.append(dict(motif_id=mid,source_family=motif['family'],target_families=';'.join(targets),
                             state='reader_annotation_only' if motif['reader_only'] else 'scanned' if targets else 'target_family_not_supplied'))
        if motif['reader_only']:
            for family in targets:
                for record in families[family]:
                    readers.append(dict(family=family,sequence_id=record['id'],motif_id=mid,
                        domain_state=ds[family,record['id']],library_score=motif['confidence_score'],
                        hypothesis='If family identity and partner context are confirmed: '+motif['molecular_function'],
                        status='reader_hypothesis_not_a_sequence_motif',source_row=motif['source_row']))
            continue
        pattern=motif['scan_pattern']
        if pattern not in cache:
            cache[pattern]={s:list(occurrences(motif['compiled'],s)) for s in allseq}
        for family in targets:
            targetseq={r['sequence'] for r in families[family]}; bg=allseq-targetseq
            n=sum(bool(cache[pattern][s]) for s in targetseq); m=sum(bool(cache[pattern][s]) for s in bg)
            prevalence.append(dict(family=family,motif_id=mid,target_positive=n,target_total=len(targetseq),
                other_positive=m,other_total=len(bg),target_fraction=n/len(targetseq),other_fraction=m/len(bg) if bg else '',
                interpretation='descriptive_unmatched_controls_not_conservation_or_significance'))
            for record in families[family]:
                found=cache[pattern][record['sequence']]
                states.append(dict(family=family,sequence_id=record['id'],motif_id=mid,hits=len(found),
                                   state='detected' if found else 'not_detected_not_family_exclusion'))
                for a,b,peptide in found:
                    tier,criteria,missing,overlaps,ev=assess(motif,rule,family,record,a,b,regions[family,record['id']],ds[family,record['id']],evidence)
                    group=hashlib.sha256(f'{family}|{record["id"]}|{a}|{b}|{peptide}'.encode()).hexdigest()[:16]
                    hits.append(dict(family=family,sequence_id=record['id'],site_group_id=group,motif_id=mid,motif_name=motif['motif_name'],
                        source_family=motif['family'],role=rule['role'],start=a,end=b,peptide=peptide,
                        flank=record['sequence'][max(0,a-11):b+10],pattern=pattern,library_score=motif['confidence_score'],
                        library_confidence_class=motif['confidence_class'],evidence_tier=motif['evidence_tier'],
                        citation_status=motif['citation_status'],screening_recommendation=motif['screening_recommendation'],
                        candidate_tier=tier,domain_state=ds[family,record['id']],localization=criteria['localization'],
                        supported_criteria=';'.join(k for k,v in criteria.items() if v=='supported'),
                        missing_or_conflicting_criteria=';'.join(missing),criteria=json.dumps(criteria,sort_keys=True),
                        overlapping_domains=';'.join(overlaps),evidence_sources=';'.join(sorted({r['source'] for r in ev})),
                        functional_hypothesis=('Context conflicts with the proposed function: ' if tier=='context_conflict' else 'If the stated context and partner requirements hold, this site may support: ')+motif['molecular_function'],
                        partner=motif['binding_partner_or_reader'],context_requirements=motif['context_requirements'],
                        next_validation=';'.join(missing) or 'Test site-specific function; hypothesis is not experimentally validated',
                        source_sheet=motif['source_sheet'],source_row=motif['source_row'],source_url=motif['source_url'],
                        primary_citation=motif['primary_literature_citation'],primary_reference=motif['primary_literature_doi_or_pmid']))
    # Positional annotation of discovered sites; do not create terminal matches by
    # scanning an isolated word and do not count duplicate library entries as evidence.
    index=defaultdict(list)
    for hit in hits: index[hit['family'],hit['sequence_id']].append(hit)
    links=[]
    for site in sites:
        for h in index[site['family'],site['sequence_id']]:
            if h['start']<=site['end'] and h['end']>=site['start']:
                links.append(dict(family=site['family'],sequence_id=site['sequence_id'],word=site['motif'],word_start=site['start'],word_end=site['end'],
                                  motif_id=h['motif_id'],motif_start=h['start'],motif_end=h['end'],site_group_id=h['site_group_id'],
                                  interpretation='positional_overlap_not_independent_validation'))
    columns=('family sequence_id site_group_id motif_id motif_name source_family role start end peptide flank pattern library_score library_confidence_class evidence_tier citation_status screening_recommendation candidate_tier domain_state localization supported_criteria missing_or_conflicting_criteria criteria overlapping_domains evidence_sources functional_hypothesis partner context_requirements next_validation source_sheet source_row source_url primary_citation primary_reference').split()
    tables=[('motif_annotations',hits,columns),('high_confidence_hypotheses',[h for h in hits if h['candidate_tier']=='high_confidence_hypothesis'],columns),
            ('functional_hypotheses',[h for h in hits if h['candidate_tier'] in ('high_confidence_hypothesis','context_supported_candidate')],columns),
            ('inspection_states',states,'family sequence_id motif_id hits state'.split()),
            ('reader_hypotheses',readers,'family sequence_id motif_id domain_state library_score hypothesis status source_row'.split()),
            ('library_audit',audit,'motif_id source_row source_family target_families role original_regex scan_regex normalization library_confidence issues'.split()),
            ('library_coverage',coverage,'motif_id source_family target_families state'.split()),
            ('motif_prevalence',prevalence,'family motif_id target_positive target_total other_positive other_total target_fraction other_fraction interpretation'.split()),
            ('discovery_annotation',links,'family sequence_id word word_start word_end motif_id motif_start motif_end site_group_id interpretation'.split())]
    for name,rows,fields in tables: write_tsv(out/(name+'.tsv'),rows,fields)
    (out/'curated_library.snapshot.json').write_text(json.dumps(sheets,indent=2,ensure_ascii=False),encoding='utf-8')
    (out/'context_rules.snapshot.json').write_text(json.dumps(rules,indent=2),encoding='utf-8')
    counts=Counter(h['candidate_tier'] for h in hits)
    summary=dict(library_entries=len(motifs),reader_entries=sum(m['reader_only'] for m in motifs),
                 annotation_rows=len(hits),distinct_sites=len({h['site_group_id'] for h in hits}),
                 candidate_tiers=dict(counts),high_confidence_rows=counts['high_confidence_hypothesis'],
                 missing_evidence=dict(Counter(c for h in hits for c in h['missing_or_conflicting_criteria'].split(';') if c)))
    (out/'curated_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    report=['# Curated motif annotation','',f"{len(motifs)} library entries; {len(hits)} annotation rows representing {summary['distinct_sites']} distinct family/sequence sites.",'',
            'Workbook confidence scores describe library rules. Candidate tiers are transparent evidence gates, not calibrated probabilities or experimental validation.',
            'A high-confidence hypothesis requires a supported family architecture, every configured context criterion, independent orthology/positional conservation/accessibility/partner evidence, and an eligible library entry.',
            'Review-only and annotation-only entries cannot become high-confidence sequence annotations. Missing evidence remains unknown. Duplicate workbook entries share site-group identifiers.',
            'Predicted outside topology can include organelle lumen; it is only provisional extracellular support. Protein-family recurrence is not orthologous conservation.','',
            '| Candidate tier | Rows |','|---|---:|']
    report += [f'| {k} | {v} |' for k,v in sorted(counts.items())]
    report += ['',f"High-confidence hypotheses: {counts['high_confidence_hypothesis']}.",'',
               'Start with functional_hypotheses.tsv and its missing-evidence/next-validation columns. motif_annotations.tsv includes all matches and conflicts. reader_hypotheses.tsv contains conditional reader functions without motif scanning.',
               'The workbook is the only source of motif patterns. Phenylalanine variants or broader ELM patterns are not silently added. Collagen [OP] is normalized to encoded P as instructed by the workbook.',
               'See library_audit.tsv for routing, normalization and source issues; library_coverage.tsv includes unsupplied target families. All coordinates are 1-based inclusive.','',
               '## Main evidence gaps','']
    report += [f'- {k}: {v} annotation rows' for k,v in Counter(summary['missing_evidence']).most_common(10)]
    (out/'CURATED_FINDINGS.md').write_text('\n'.join(report)+'\n',encoding='utf-8')
    return summary
