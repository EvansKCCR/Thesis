"""Join stage 2 evidence onto the preserved stage 1 species/protein hypotheses."""

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
from collections import defaultdict
import json
from pathlib import Path
import re


def read(path):
    if not _layout_resource(path).exists():return []
    with _layout_resource(path).open(encoding='utf-8-sig',newline='') as handle:return list(csv.DictReader(handle,delimiter='\t'))


def integer(value):
    return max([int(n) for n in re.findall(r'\d+',str(value or ''))],default=0)


def join(annotation,handoff,panel,orthology,rows,benchmark):
    """No reference accession, HOG or exact sequence can replace a family hypothesis."""
    motifs=_layout_resource(annotation/'motifs')
    if panel is not None:
        targets=json.loads((_layout_resource(panel/'targets.json')).read_text())
        panel_manifest=_layout_resource(panel/'workflow_manifest.json')
        panel_sha256=hashlib.sha256(_layout_resource(panel_manifest).read_bytes()).hexdigest() if _layout_resource(panel_manifest).exists() else ''
        species=json.loads((_layout_resource(panel/'species.json')).read_text());references=[s for s in species if s not in targets]
    else:
        manifest=_layout_resource(handoff/'workflow_manifest.json')
        panel_manifest_data=json.loads(_layout_resource(manifest).read_text(encoding='utf-8')) if _layout_resource(manifest).is_file() else {}
        targets=panel_manifest_data.get('target_species') or sorted({r.get('Specie','') for r in rows if r.get('Specie')})
        species=list(targets);references=[];panel_sha256=''
    report=read(_layout_resource(orthology/'candidates/candidate_report.tsv')) if orthology is not None else []
    reports={r['candidate_id']:r for r in report}
    hogs=read(_layout_resource(orthology/'candidates/candidate_HOGs.tsv')) if orthology is not None else []
    used=read(_layout_resource(motifs/'used_regions.tsv'));contexts=read(_layout_resource(motifs/'sequence_context.tsv'));cdds=read(_layout_resource(motifs/'CDD_specific_hits.tsv'))
    source=read(_layout_resource(annotation/'external_evidence.tsv'));external={r['protein_id']:r for r in source}
    context={r['sequence_id']:r for r in contexts}
    regions=defaultdict(set)
    for item in used:regions[item['sequence_id']].add((item['region'],item['start'],item['end']))
    cdd_by_id=defaultdict(dict)
    for item in cdds:
        cdd_by_id[item['sequence_id']][(item['accession'],item['start'],item['end'])]=item
    snapshot=read(_layout_resource(motifs/'current_motif_rules.tsv'))
    ids={r['rule_id'] for r in snapshot}
    if not ids:
        library=_layout_resource(Path(__file__).resolve().parents[1]/'motif_annotation/Motif_library/adhesion_core.json')
        ids={r['rule_id'] for r in json.loads(_layout_resource(library).read_text())['rules']}
    supplied=read(_layout_resource(motifs/'site_evidence.tsv'))
    evidence=defaultdict(list)
    for item in supplied:evidence[item['sequence_id'],item['motif_id']].append(item)
    annotations=read(_layout_resource(motifs/'motif_annotations.tsv')) or read(_layout_resource(motifs/'adhesion_motif_annotations.tsv'))
    sites=defaultdict(dict);ignored=defaultdict(int)
    priority={'context_conflict':3,'context_supported_candidate':2,'sequence_only_candidate':1}
    updated_annotations=[]
    for raw in annotations:
        item=dict(raw);pid=item['sequence_id'];mid=item['motif_id']
        if mid not in ids:ignored[pid]+=1;continue
        matching=[r for r in evidence[pid,mid] if int(r['start'])<=int(item['start'])<=int(item['end'])<=int(r['end']) and r.get('source') and r.get('criterion')]
        old=item.get('candidate_tier','sequence_only_candidate')
        tier='context_conflict' if old=='context_conflict' or any(r['status']=='contradicted' for r in matching) else 'context_supported_candidate' if 'context_supported' in old or any(r['status']=='supported' for r in matching) else 'sequence_only_candidate'
        item.update(reported_candidate_tier=old,candidate_tier=tier);updated_annotations.append(item)
        key=(mid,item['start'],item['end'],item.get('peptide',''))
        if key not in sites[pid] or priority[tier]>priority[sites[pid][key]['candidate_tier']]:sites[pid][key]=item
    joined=[];topology=[]
    for original in rows:
        row=dict(original);pid=row['protein_id'];rep=reports.get(pid,{})
        row.update(benchmark)
        row['orthology_panel_sha256']=panel_sha256
        row['domain_order_check']=row.get('order_check',row.get('domain_order_check',''))
        row['Orthology_mapping_status']=rep.get('mapping_status','not_run' if orthology is None else 'not_in_candidate_report')
        row['Orthology_Orthogroup']=rep.get('Orthogroup','');row['Orthology_HOG']=rep.get('HOG','')
        row['Orthology_HOG_status']=rep.get('HOG_status','not_run' if orthology is None else 'not_in_candidate_report')
        row['Orthology_analysis_completed']=bool(rep)
        row['Orthology_target_species']=';'.join(targets)
        reference_present=[];direct=[]
        linked={h.strip() for h in row['Orthology_HOG'].split(';') if h.strip()}
        source_hogs=[r for r in hogs if r.get('HOG',r.get('Orthogroup','')) in linked]
        for sp in species:
            copies=integer(rep.get(sp+'_HOG_copy_counts'))
            members=sorted({g.strip() for h in source_hogs for g in h.get(sp,'').split(',') if g.strip()})
            if not copies:copies=len(members)
            row['Orthology_'+sp+'_HOG_copy_count']=copies
            if sp in references:
                row['Orthology_'+sp+'_direct_orthologues']=rep.get(sp+'_orthologues','')
                row['Orthology_'+sp+'_relationship']=rep.get(sp+'_relationship','')
                row['Orthology_'+sp+'_HOG_members']='; '.join(members)
                if copies or members or rep.get(sp+'_orthologues'):reference_present.append(sp)
                if rep.get(sp+'_orthologues'):direct.append(sp)
        row['Orthology_HOG_reference_species_count']=len(reference_present)
        row['Orthology_direct_reference_count']=len(direct)
        row['Orthology_reference_support']='; '.join(reference_present)
        row['Orthology_interpretation']=('Reference-species representation in the sampled HOG panel; benchmark weight applies to orthology score only'
            if orthology is not None else 'Orthology not run; this preliminary audit uses available non-orthology evidence')
        ss=list(sites[pid].values());tiers={r['candidate_tier'] for r in ss}
        row.update(MotifPrediction_raw_match_count=len(ss),MotifPrediction_distinct_site_count=len({(r['start'],r['end'],r.get('peptide','')) for r in ss}),
            MotifPrediction_library_ids='; '.join(sorted({r['motif_id'] for r in ss})),MotifPrediction_candidate_tiers=';'.join(sorted(tiers)),
            MotifPrediction_context_supported_count=sum(r['candidate_tier']=='context_supported_candidate' for r in ss),
            MotifPrediction_context_conflict_count=sum(r['candidate_tier']=='context_conflict' for r in ss),
            MotifPrediction_sequence_only_count=sum(r['candidate_tier']=='sequence_only_candidate' for r in ss),
            MotifPrediction_priority_site_count=sum(r['candidate_tier']!='context_conflict' for r in ss),
            MotifPrediction_review_only_count=sum(r['candidate_tier']=='context_conflict' for r in ss),
            MotifPrediction_ignored_superseded_count=ignored[pid],
            MotifPrediction_priority_sites='; '.join(f"{r['motif_id']}:{r['start']}-{r['end']}" for r in ss[:20])+(' (additional sites in Motif_annotations)' if len(ss)>20 else ''),
            MotifPrediction_context_summary='; '.join(sorted({r.get('site_evidence','') for r in ss}-{''})))
        ctx=context.get(pid,{});ext=external.get(pid,{})
        row['DeepLoc_2.1']=ext.get('DeepLoc_2_1') or ctx.get('deeploc_localizations','')
        row['SignalP.v6.0']=ext.get('SignalP_v6_0') or ctx.get('signalp_prediction','')
        row['DeepTMHMM']=ext.get('DeepTMHMM') or '; '.join(('TM' if name=='transmembrane' else 'SP')+' '+a+'-'+b for name,a,b in sorted(regions[pid]) if name in ('transmembrane','signal_peptide'))
        cd=sorted(cdd_by_id[pid].values(),key=lambda hit:(int(hit['start']),int(hit['end']),hit['short_name']))
        row['NCBI_CDD_hits']=json.dumps(cd,sort_keys=True)
        row['NCBI_CDD_architecture']=';'.join(dict.fromkeys(text for text in
            [ext.get('NCBI_CDD_architecture',''),*(f"{r['short_name']}:{r['start']}-{r['end']}" for r in cd)] if text))
        row['NCBI_CDD_source']=ext.get('source','') or ('verified stage 2 CDD selected hits' if cd else '')
        row['NCBI_CDD_selected_hit_count']=len(cd)
        row['NCBI_CDD_best_evalue']=min((float(r['evalue']) for r in cd),default=None)
        row['External_evidence_source']=ext.get('source','')
        joined.append(row)
        topology.append(dict(protein_id=pid,DeepLoc_Localizations=row['DeepLoc_2.1'],DeepLoc_Cytoplasm=ctx.get('deeploc_cytoplasm',''),
            DeepLoc_Extracellular=ctx.get('deeploc_extracellular',''),DeepTMHMM_summary=row['DeepTMHMM'],SignalP_prediction=row['SignalP.v6.0']))
    topology=list({r['protein_id']:r for r in topology}.values())
    support=dict(reference_species=references,target_species=targets,benchmark=benchmark,tables=dict(
        Motif_annotations=updated_annotations,Orthology_candidate_report=report,Orthology_HOGs=hogs,
        Direct_orthologues=read(_layout_resource(orthology/'candidates/direct_orthologues.tsv')),NCBI_CDD_selected_hits=cdds,NCBI_CDD_all_hits=read(_layout_resource(motifs/'CDD_hits.tsv')) if (_layout_resource(motifs/'CDD_hits.tsv')).is_file() else cdds,Topology_predictions=topology))
    return joined,support
