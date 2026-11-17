"""Current evidence-based, deduplicated motif snapshot adapter."""

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
from pathlib import Path
import re
import warnings
import curated_annotation as ca


def load(path,scope='adhesion'):
    path=_layout_resource(Path(path))
    if path.suffix.lower()=='.xlsx':
        from openpyxl import load_workbook
        book=load_workbook(path,read_only=True,data_only=True)
        try:
            sheet=book['Adhesion_core' if scope=='adhesion' else 'Master_library']
            rows=iter(sheet.values); header=next(rows)
            rules=[dict(zip(header,row)) for row in rows if row[0]]
            for name in ('Context_evidence','Evidence'):
                if name not in book.sheetnames:continue
                entries=iter(book[name].values); columns=next(entries)
                context={}
                for values in entries:
                    entry=dict(zip(columns,values)); rid=entry.get('rule_id')
                    if rid:context.setdefault(rid,[]).append(entry)
                for rule in rules:rule[name]=context.get(rule['rule_id'],[])
        finally:
            book.close()
        snapshot=dict(schema_version='1.0',rules=rules,context_filters_implemented=False)
    else:
        snapshot=json.loads(_layout_resource(path).read_text(encoding='utf-8'))
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
    with _layout_resource(Path(path)).open(encoding='utf-8-sig',newline='') as handle:
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


def _context_class(rule):
    structured=rule.get('Context_evidence',[])
    compartments=' '.join(str(x.get('compartment_context','')) for x in structured).casefold()
    if compartments:
        outside='extracellular' in compartments
        inside=any(x in compartments for x in ('cytoplas','cytosol','intracellular'))
        if outside != inside:return 'extracellular' if outside else 'intracellular'
    text=' '.join(str(rule.get(k,'')) for k in ('motif_name','context_requirements','sequence_bearer','binding_reader')).casefold()
    if any(x in text for x in ('extracellular','secreted ligand','laminin-derived','fibronectin')): return 'extracellular'
    if any(x in text for x in ('cytoplasmic','cytosolic','adaptor','sorting signal','sh2 reader','ptb reader','pdz reader','evh1 reader')): return 'intracellular'
    return 'unspecified'


def _context_compatible(rule, segment, domains, family):
    required=_context_class(rule)
    seg=str(segment).casefold()
    if required=='extracellular' and any(x in seg for x in ('predicted_inside','cytoplasmic_tail')): return False,'topology_conflict'
    if required=='intracellular' and any(x in seg for x in ('predicted_outside','signal_peptide')): return False,'topology_conflict'
    bearer=str(rule.get('sequence_bearer','')).casefold()
    fam=str(family).casefold().replace('_',' ')
    if required=='extracellular' and any(x in bearer for x in ('fibronectin','laminin','collagen')) and fam and not any(x in fam for x in ('fibronectin','laminin','collagen','extracellular','unassigned','component')):
        return False,'family_bearer_conflict'
    if required=='intracellular' and any(x in bearer for x in ('integrin alpha cytoplasmic','integrin beta cytoplasmic')) and 'integrin' not in fam:
        return False,'family_bearer_conflict'
    return True,'context_compatible' if required!='unspecified' else 'context_unresolved'


def scan(snapshot,families,regions,evidence,out,write_tsv,sequence_context=None,stage1_context=None):
    sequence_context=sequence_context or {}; stage1_context=stage1_context or {}
    out=_layout_resource(Path(out))
    fields='family sequence_id stage1_hypotheses stage1_modules stage1_status motif_id library_id motif_name start end peptide site_group_id candidate_tier overlapping_domains predicted_topology_segment deeploc_localizations signalp_prediction context_requirements sequence_bearer binding_reader functional_hypothesis context_filter_status evidence_basis evidence_status source_class_ids source_row_ids source_url scan_scope site_evidence'.split()
    adhesion=[]; regulatory=[]; seen_sites=set(); suppressed=0; redundant=0; feature_counts={}
    for rule in snapshot['rules']:
        for family,records in families.items():
            for record in records:
                seq=record['sequence']
                for a,b,peptide in ca.occurrences(rule['compiled'],seq):
                    relevant=[r for r in evidence if r['family']==family and r['sequence_id']==record['id'] and
                              r['motif_id']==rule['rule_id'] and r['start']<=a<=b<=r['end']]
                    ctx=sequence_context.get((family,record['id']),{})
                    rs=regions[family,record['id']]
                    domains=';'.join(sorted({n for n,x,y in rs if (n.startswith('PF') or n.startswith('CDD:')) and a<=y and b>=x}))
                    top=';'.join(sorted({n for n,x,y in rs if n in ('predicted_inside','predicted_outside','transmembrane','signal_peptide','cytoplasmic_tail') and x<=a<=b<=y}))
                    s1=stage1_context.get(record['id'],{})
                    hypotheses=s1.get('families',[]) or [family]
                    modules=s1.get('modules',[])
                    compatible_checks=[_context_compatible(rule,top,domains,h) for h in hypotheses]
                    compatible=any(x[0] for x in compatible_checks)
                    context_status=';'.join(sorted({x[1] for x in compatible_checks}))
                    if not compatible: suppressed+=1; continue
                    dedupe=(family,record['id'],a,b,peptide,rule.get('overlap_group') or rule['rule_id'])
                    if dedupe in seen_sites: redundant+=1; continue
                    seen_sites.add(dedupe)
                    feature_counts.setdefault((family,record['id'],rule['rule_id']),[]).append((a,b,peptide))
                    tier='context_conflict' if any(r['status']=='contradicted' for r in relevant) else 'context_supported_candidate' if any(r['status']=='supported' for r in relevant) else 'sequence_only_candidate'
                    target=adhesion if rule['scan_scope']=='adhesion_associated' else regulatory
                    target.append(dict(family=family,sequence_id=record['id'],stage1_hypotheses=';'.join(hypotheses),stage1_modules=';'.join(modules),stage1_status=s1.get('status',''),motif_id=rule['rule_id'],library_id=rule['rule_id'],
                        motif_name=rule['motif_name'],start=a,end=b,peptide=peptide,
                        site_group_id=hashlib.sha256(f'{record["id"]}|{a}|{b}|{peptide}'.encode()).hexdigest()[:16],
                        candidate_tier=tier,overlapping_domains=domains,predicted_topology_segment=top,
                        deeploc_localizations=ctx.get('deeploc_localizations',''),signalp_prediction=ctx.get('signalp_prediction',''),
                        context_requirements=rule['context_requirements'],sequence_bearer=rule['sequence_bearer'],binding_reader=rule['binding_reader'],
                        functional_hypothesis=rule['interpretation'],context_filter_status=context_status,evidence_basis=rule.get('evidence_basis',''),evidence_status=rule.get('evidence_status',''),
                        source_class_ids=rule.get('source_class_ids',''),source_row_ids=rule.get('source_row_ids',''),source_url=rule.get('source_url',''),
                        scan_scope=rule['scan_scope'],site_evidence=';'.join(f"{r['criterion']}:{r['status']}:{r['source']}" for r in relevant)))
    annotations=adhesion+regulatory
    write_tsv(_layout_resource(out/'motif_annotations.tsv'),adhesion,fields)
    write_tsv(_layout_resource(out/'adhesion_motif_annotations.tsv'),adhesion,fields)
    write_tsv(_layout_resource(out/'general_regulatory_motifs.tsv'),regulatory,fields)
    core_ids={r['rule_id'] for r in snapshot['rules'] if r['scan_scope']=='adhesion_associated'}
    feats=[dict(family=f,sequence_id=i,motif_id=m,occurrence_count=len(v),copy_number=len(v),positions_1based_inclusive=';'.join(f'{a}-{b}' for a,b,_ in v),matched_peptides=';'.join(x[2] for x in v)) for (f,i,m),v in sorted(feature_counts.items())]
    feature_fields=['family','sequence_id','motif_id','occurrence_count','copy_number','positions_1based_inclusive','matched_peptides']
    write_tsv(_layout_resource(out/'adhesion_core_motif_features.tsv'),[x for x in feats if x['motif_id'] in core_ids],feature_fields)
    write_tsv(_layout_resource(out/'regulatory_motif_features.tsv'),[x for x in feats if x['motif_id'] not in core_ids],feature_fields)
    write_tsv(_layout_resource(out/'site_evidence.tsv'),evidence,['family','sequence_id','motif_id','start','end','criterion','status','source'])
    serial=[{k:v for k,v in r.items() if k!='compiled'} for r in snapshot['rules']]
    write_tsv(_layout_resource(out/'current_motif_rules.tsv'),serial,list(serial[0]))
    return dict(master_records=len(serial),sequence_annotation_rows=len(annotations),adhesion_core_rows=len(adhesion),regulatory_rows=len(regulatory),
        generic_hit_rows=len(regulatory),
        context_filters_implemented=True,context_suppressed_rows=suppressed,redundant_rows_removed=redundant,supplied_site_context_evidence_used=bool(evidence),
        context_supported_rows=sum(r['candidate_tier']=='context_supported_candidate' for r in annotations),
        warning='Sequence-only matches retain partial support; supplied verified site-context evidence is reported separately.',summary_discrepancies=[])
