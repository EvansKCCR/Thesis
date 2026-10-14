"""Current evidence-based, deduplicated motif snapshot adapter."""
import csv
import hashlib
import json
from pathlib import Path
import re
import warnings
import curated_annotation as ca


def load(path,scope='adhesion'):
    path=Path(path)
    if path.suffix.lower()=='.xlsx':
        from openpyxl import load_workbook
        book=load_workbook(path,read_only=True,data_only=True)
        try:
            sheet=book['Adhesion_core' if scope=='adhesion' else 'Master_library']
            rows=iter(sheet.values); header=next(rows)
            rules=[dict(zip(header,row)) for row in rows if row[0]]
        finally:
            book.close()
        snapshot=dict(schema_version='1.0',rules=rules,context_filters_implemented=False)
    else:
        snapshot=json.loads(path.read_text(encoding='utf-8'))
    required={'rule_id','motif_name','regex_python','scan_scope','context_requirements',
              'sequence_bearer','binding_reader','interpretation','regex_sha256'}
    rules=[]; ids=set(); patterns=set()
    for source in snapshot['rules']:
        rule={k:('' if v is None else v) for k,v in source.items()}
        if not required<=rule.keys():
            raise ValueError('Current motif rule schema is incomplete')
        if scope=='adhesion' and rule['scan_scope']!='adhesion_associated':
            continue
        mid,pattern=rule['rule_id'],rule['regex_python']
        if mid in ids or pattern in patterns:
            raise ValueError('Duplicate motif ID or regex in current snapshot: '+mid)
        ids.add(mid); patterns.add(pattern)
        if hashlib.sha256(pattern.encode()).hexdigest()!=rule['regex_sha256']:
            raise ValueError('Motif regex checksum mismatch: '+mid)
        with warnings.catch_warnings():
            warnings.simplefilter('error',FutureWarning)
            compiled=re.compile(pattern)
        if compiled.match(''):
            raise ValueError('Empty-matching motif: '+mid)
        rules.append(rule|dict(compiled=compiled))
    if not rules:
        raise ValueError('No current motif rules selected')
    return snapshot|dict(rules=rules)


def load_evidence(path,families,snapshot):
    if not path:return []
    sequences={(f,r['id']):r['sequence'] for f,rs in families.items() for r in rs}
    ids={r['rule_id'] for r in snapshot['rules']}
    result=[]
    with Path(path).open(encoding='utf-8-sig',newline='') as handle:
        reader=csv.DictReader(handle,delimiter='\t')
        required={'family','sequence_id','sequence_sha256','motif_id','start','end','criterion','status','source'}
        if not required<=set(reader.fieldnames or []):raise ValueError('Invalid site-evidence columns')
        for row in reader:
            key=row['family'],row['sequence_id']
            if key not in sequences or hashlib.sha256(sequences[key].encode()).hexdigest()!=row['sequence_sha256']:
                raise ValueError('Unknown sequence or stale evidence hash: '+str(key))
            if row['motif_id'] not in ids:raise ValueError('Use a current MOTIF_ rule ID: '+row['motif_id'])
            if row['status'] not in ('supported','contradicted') or not row['source'].strip() or not row['criterion'].strip():
                raise ValueError('Site evidence needs criterion, supported/contradicted status and source')
            a,b=int(row['start']),int(row['end'])
            if not 1<=a<=b<=len(sequences[key]):raise ValueError('Invalid site-evidence interval')
            result.append(row|dict(start=a,end=b))
    return result


def scan(snapshot,families,regions,evidence,out,write_tsv,sequence_context=None):
    sequence_context=sequence_context or {}
    out=Path(out)
    fields='family sequence_id motif_id library_id motif_name start end peptide site_group_id candidate_tier overlapping_domains predicted_topology_segment deeploc_localizations signalp_prediction context_requirements sequence_bearer binding_reader functional_hypothesis evidence_basis evidence_status source_class_ids source_row_ids source_url scan_scope site_evidence'.split()
    annotations=[]; cache={}
    for rule in snapshot['rules']:
        for family,records in families.items():
            for record in records:
                seq=record['sequence']
                key=(rule['rule_id'],seq)
                if key not in cache:cache[key]=list(ca.occurrences(rule['compiled'],seq))
                for a,b,peptide in cache[key]:
                    relevant=[r for r in evidence if r['family']==family and r['sequence_id']==record['id'] and
                              r['motif_id']==rule['rule_id'] and r['start']<=a<=b<=r['end']]
                    ctx=sequence_context.get((family,record['id']),{})
                    rs=regions[family,record['id']]
                    domains=';'.join(sorted({n for n,x,y in rs if (n.startswith('PF') or n.startswith('CDD:')) and a<=y and b>=x}))
                    top=';'.join(sorted({n for n,x,y in rs if n in ('predicted_inside','predicted_outside','transmembrane','signal_peptide','cytoplasmic_tail') and x<=a<=b<=y}))
                    tier='context_conflict' if any(r['status']=='contradicted' for r in relevant) else 'sequence_only_candidate'
                    annotations.append(dict(family=family,sequence_id=record['id'],motif_id=rule['rule_id'],library_id=rule['rule_id'],
                        motif_name=rule['motif_name'],start=a,end=b,peptide=peptide,
                        site_group_id=hashlib.sha256(f'{record["id"]}|{a}|{b}|{peptide}'.encode()).hexdigest()[:16],
                        candidate_tier=tier,overlapping_domains=domains,predicted_topology_segment=top,
                        deeploc_localizations=ctx.get('deeploc_localizations',''),signalp_prediction=ctx.get('signalp_prediction',''),
                        context_requirements=rule['context_requirements'],sequence_bearer=rule['sequence_bearer'],binding_reader=rule['binding_reader'],
                        functional_hypothesis=rule['interpretation'],evidence_basis=rule.get('evidence_basis',''),evidence_status=rule.get('evidence_status',''),
                        source_class_ids=rule.get('source_class_ids',''),source_row_ids=rule.get('source_row_ids',''),source_url=rule.get('source_url',''),
                        scan_scope=rule['scan_scope'],site_evidence=';'.join(f"{r['criterion']}:{r['status']}:{r['source']}" for r in relevant)))
    write_tsv(out/'motif_annotations.tsv',annotations,fields)
    write_tsv(out/'adhesion_motif_annotations.tsv',[r for r in annotations if r['scan_scope']=='adhesion_associated'],fields)
    write_tsv(out/'site_evidence.tsv',evidence,['family','sequence_id','motif_id','start','end','criterion','status','source'])
    serial=[{k:v for k,v in r.items() if k!='compiled'} for r in snapshot['rules']]
    write_tsv(out/'current_motif_rules.tsv',serial,list(serial[0]))
    return dict(master_records=len(serial),sequence_annotation_rows=len(annotations),
        generic_hit_rows=sum(r['scan_scope']!='adhesion_associated' for r in annotations),
        context_filters_implemented=False,warning='Hits are sequence candidates; context requirements remain to be verified.',summary_discrepancies=[])
