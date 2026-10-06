#!/usr/bin/env python3
"""Auditable per-HOG MAFFT -> trimAl -> IQ-TREE workflow (Python >=3.11)."""
import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys

from pipeline import SPECIES, table_species, digest, members, read_tsv, records, write_tsv, configure_species


def dump(path, value):
    path.write_text(json.dumps(value, indent=2)+'\n', encoding='utf-8')


def fasta(path, sequences):
    with path.open('w', encoding='ascii') as f:
        for name, seq in sorted(sequences.items()):
            f.write(f'>{name}\n')
            f.write('\n'.join(seq[i:i+80] for i in range(0,len(seq),80))+'\n')


def load_sequences(path):
    result = {}
    for header, seq in records(path):
        name = header.split()[0]
        if name in result or not seq:
            raise ValueError(f'Duplicate/empty record in {path}: {name}')
        result[name] = seq
    return result


def family_plans(plans, family_map, family_column, active_species, candidate_report=None, selected=None,
                 audit_workbook=None, target_species=None, sequences=None):
    """Union whole source groups by candidate-family association, retaining provenance."""
    if audit_workbook:
        from openpyxl import load_workbook
        with audit_workbook.open('rb') as stream:
            book=load_workbook(stream,read_only=True,data_only=True)
            try:
                values=iter(book['Collapsed_candidates'].values)
                fields=next(values)
                annotations=[dict(zip(fields,row)) for row in values]
            finally:
                book.close()
        family_column='assigned_family'
        if 'grade' not in fields:
            raise ValueError('Collapsed_candidates requires an audit grade column')
        annotations=[r for r in annotations if r['grade'] in ('A','B','C')]
        if not target_species or not set(target_species)<=set(active_species):
            raise ValueError('Audit selection requires valid --target-species (panel/targets.json)')
        if sequences is None:
            raise ValueError('Audit selection requires prepared proteome sequences')
    else:
        fields, annotations=read_tsv(family_map)
    if 'protein_id' not in (fields or []) or family_column not in (fields or []):
        raise ValueError(f'Family map requires protein_id and {family_column} columns')
    aliases={}
    if candidate_report:
        _, report=read_tsv(candidate_report)
        for row in report:
            cid, prepared=row['candidate_id'],row['prepared_id']
            if cid in aliases and aliases[cid]!=prepared:
                raise ValueError(f'Conflicting candidate mapping: {cid}')
            aliases[cid]=prepared
    families=defaultdict(set)
    audit_details={}
    for row in annotations:
        family=str(row[family_column] or '').strip()
        if family in ('','unassigned','unresolved'):
            continue
        if not re.fullmatch(r'[A-Za-z0-9_.-]+',family) or family in ('.','..'):
            raise ValueError(f'Unsafe family name: {family}')
        ident=aliases.get(row['protein_id'],row['protein_id'])
        if audit_workbook and (ident not in sequences or ident.split('__')[0] not in target_species):
            raise ValueError(f'Audited candidate missing from target proteomes: {ident}')
        families[family].add(ident)
        if audit_workbook:
            if ident in audit_details:
                raise ValueError(f'Duplicate collapsed audited candidate: {ident}')
            audit_details[ident]=dict(protein_id=row['protein_id'],grade=row['grade'])
    if selected:
        absent=set(selected)-set(families)
        if absent:
            raise ValueError(f'Unknown families: {sorted(absent)}')
        families={f:ids for f,ids in families.items() if f in selected}
    merged,links,coverage=[],[],[]
    for family,wanted in sorted(families.items()):
        selected_groups=[(g,row,ids) for g,row,ids in plans if ids & wanted]
        found=set().union(*(ids & wanted for g,row,ids in selected_groups))
        for cid in sorted(wanted):
            coverage.append(dict(family=family,prepared_id=cid,status='included' if cid in found else
                                 ('included_without_source_group' if audit_workbook else 'not_in_selected_groups')))
            if audit_workbook:
                coverage[-1].update(audit_details[cid])
        if not selected_groups and not audit_workbook:
            continue
        ids=set().union(*(ids for g,row,ids in selected_groups))
        if audit_workbook:
            ids={n for n in ids if n.split('__')[0] not in target_species}|wanted
            found=wanted
        row={sp:', '.join(sorted(n for n in ids if n.startswith(sp+'__'))) for sp in active_species}
        row.update(Matched_candidates=', '.join(sorted(found)),Source_groups=';'.join(g for g,r,i in selected_groups),Family=family)
        source_by_tip=defaultdict(list)
        for group,original,group_ids in selected_groups:
            links.append(dict(family=family,source_group=group,trigger_candidates=';'.join(sorted(group_ids&wanted)),member_count=len(group_ids)))
            for ident in group_ids & ids:
                source_by_tip[ident].append(group)
        for ident in ids:
            source_by_tip.setdefault(ident,[])
        row['_source_by_tip']=source_by_tip
        merged.append(('family_'+family,row,ids))
    if not merged:
        raise ValueError('No source groups linked to the selected family hypotheses')
    return merged,links,coverage


def prepare(args):
    configure_species(getattr(args, 'species_list', None))
    fields, rows = read_tsv(args.hogs)
    active_species = table_species(fields)
    if args.duplications is None and not getattr(args,'allow_missing_duplications',False):
        raise ValueError('Supply duplication evidence or explicitly allow it to remain unavailable')
    files = [args.hogs,args.orthologues]+([args.duplications] if args.duplications else [])+[args.proteomes/f'{s}.faa' for s in active_species]
    for path in files:
        if not path.is_file():
            raise ValueError(f'Missing input: {path}')
    key = 'HOG' if 'HOG' in fields else 'Orthogroup'
    if key not in fields:
        raise ValueError('HOG table requires a HOG/Orthogroup identifier')
    if args.hog:
        absent = set(args.hog)-{r[key] for r in rows}
        if absent:
            raise ValueError(f'Unknown HOGs: {sorted(absent)}')
        rows = [r for r in rows if r[key] in args.hog]
    if (not rows and not getattr(args,'audit_workbook',None)) or len({r[key] for r in rows}) != len(rows):
        raise ValueError('Empty selection or duplicate group IDs')
    sequences = {}
    for sp in active_species:
        seqs = load_sequences(args.proteomes/f'{sp}.faa')
        if any(not n.startswith(sp+'__') for n in seqs):
            raise ValueError(f'Incorrect species prefix in {sp}.faa')
        sequences.update(seqs)
    _, orthologues = read_tsv(args.orthologues)
    dupfields, duplications = read_tsv(args.duplications) if args.duplications else (['Genes 1','Genes 2'],[])
    plans = []
    for row in rows:
        group = row[key]
        if not re.fullmatch(r'[A-Za-z0-9_.-]+', group) or group in ('.','..'):
            raise ValueError(f'Unsafe group ID: {group}')
        ids = set.union(*(members(row[s]) for s in active_species))
        if ids-set(sequences):
            raise ValueError(f'{group}: missing sequences {sorted(ids-set(sequences))}')
        for sp in active_species:
            if any(not n.startswith(sp+'__') for n in members(row[sp])):
                raise ValueError(f'{group}: incorrect species column {sp}')
        plans.append((group,row,ids))
    family_map=getattr(args,'family_map',None)
    audit_workbook=getattr(args,'audit_workbook',None)
    family_mode=bool(family_map or audit_workbook)
    target_species=None
    if audit_workbook:
        if not getattr(args,'target_species',None):
            raise ValueError('--audit-workbook requires --target-species panel/targets.json')
        target_species=json.loads(args.target_species.read_text())
        if not isinstance(target_species,list) or not all(isinstance(s,str) for s in target_species):
            raise ValueError('--target-species must contain a JSON list of species codes')
        files += [audit_workbook,args.target_species]
    links,coverage=[],[]
    if family_mode:
        plans,links,coverage=family_plans(plans,family_map,args.family_column,active_species,
            getattr(args,'candidate_report',None),getattr(args,'family',None),audit_workbook,target_species,sequences)
        files+=([family_map] if family_map else [])+([args.candidate_report] if getattr(args,'candidate_report',None) else [])
    elif getattr(args,'family',None):
        raise ValueError('--family requires --audit-workbook or --family-map')
    args.output.mkdir(parents=True, exist_ok=False)
    manifest = {'created_utc':datetime.now(timezone.utc).isoformat(),
        'inputs':{str(p.resolve()):digest(p) for p in files}, 'species':active_species,
        'duplication_evidence_status':'provided' if args.duplications else 'unavailable_not_evidence_of_no_duplication',
        'grouping':'family' if family_mode else 'orthogroup',
        'selection_policy':'audited_ABC_targets_plus_reference_group_members' if audit_workbook else 'full_group_union',
        'target_species':target_species,
        'family_column':'assigned_family' if audit_workbook else (args.family_column if family_map else None), 'groups':[]}
    if family_mode:
        write_tsv(args.output/'family_orthogroups.tsv',['family','source_group','trigger_candidates','member_count'],links)
        write_tsv(args.output/'family_candidate_coverage.tsv',['family','prepared_id','status']+
                  (['protein_id','grade'] if audit_workbook else []),coverage)
    if audit_workbook:
        manifest['selected_candidates']=len(coverage)
        manifest['selected_grades']=dict(Counter(r['grade'] for r in coverage))
        print(f"Audit selection: {len(coverage)} candidates; grades {manifest['selected_grades']}")
    summary = []
    for group,row,ids in plans:
        folder = args.output/group
        folder.mkdir()
        chosen = {n:sequences[n] for n in ids}
        fasta(folder/'members.faa',chosen)
        matched = members(row.get('Matched_candidates',''))
        tips = [dict(sequence_id=n,species=n.split('__')[0],length=len(sequences[n]),
                     candidate=int(n in matched),family=row.get('Family',''),
                     source_groups=';'.join(row.get('_source_by_tip',{}).get(n,[group]))) for n in sorted(ids)]
        write_tsv(folder/'tips.tsv',['sequence_id','species','length','candidate','family','source_groups'],tips)
        edges = []
        for r in orthologues:
            for other in members(r['orthologues']):
                if r['candidate_id'] in ids and other in ids:
                    edges.append(dict(candidate_id=r['candidate_id'],orthologue_id=other,
                                      relationship=r['relationship'],source_file=r['source_file']))
        write_tsv(folder/'orthology_edges.tsv',['candidate_id','orthologue_id','relationship','source_file'],edges)
        # Gene names, rather than old gene-tree node numbers, link external evidence to a new tree.
        events = [r for r in duplications if (members(r.get('Genes 1',''))|members(r.get('Genes 2',''))) & ids]
        write_tsv(folder/'duplication_evidence.tsv',dupfields,events)
        dump(folder/'duplication_status.json',dict(status=manifest['duplication_evidence_status']))
        status = 'ready' if len(ids)>=4 else 'too_few_taxa_for_supported_tree'
        entry = dict(group=group,family=row.get('Family',''),source_groups=row.get('Source_groups',group),sequences=len(ids),unique_sequences=len(set(chosen.values())),
                     candidates=len(matched & ids),orthology_edges=len(edges),duplication_rows=len(events),status=status)
        if audit_workbook:
            entry['reference_sequences']=len(ids-matched)
        summary.append(entry)
        manifest['groups'].append(dict(entry,files={p.name:digest(p) for p in folder.iterdir()}))
    write_tsv(args.output/'groups.tsv',list(summary[0]),summary)
    dump(args.output/'manifest.json',manifest)
    print(f'Prepared {len(plans)} separate groups in {args.output}; '+
          f'{sum(r["status"]=="ready" for r in summary)} have at least four taxa.')


def alignment_qc(path, expected, min_sites, max_gap):
    seqs = load_sequences(path)
    if set(seqs) != set(expected):
        dump(path.with_suffix('.membership.json'),dict(missing=sorted(set(expected)-set(seqs)),unexpected=sorted(set(seqs)-set(expected))))
        raise ValueError(f'{path}: alignment changed sequence identifiers; inspect {path.with_suffix(".membership.json")}')
    lengths = {len(s) for s in seqs.values()}
    if len(lengths)!=1:
        raise ValueError(f'{path}: unequal alignment lengths')
    length = lengths.pop()
    if length < min_sites:
        raise ValueError(f'{path}: only {length} sites; minimum is {min_sites}')
    rows = [dict(sequence_id=n,sites=length,gap_fraction=s.count('-')/length,
                 informative_residues=sum(c in 'ACDEFGHIKLMNPQRSTVWY' for c in s)) for n,s in seqs.items()]
    write_tsv(path.with_suffix('.sequence_qc.tsv'),list(rows[0]),rows)
    if any(r['gap_fraction']>max_gap or r['informative_residues']==0 for r in rows):
        raise ValueError(f'{path}: excessive gaps or no standard residues in a sequence; review alignment')
    informative = 0
    for column in zip(*seqs.values()):
        counts = Counter(c for c in column if c in 'ACDEFGHIKLMNPQRSTVWY')
        informative += sum(v>=2 for v in counts.values())>=2
    return seqs, rows, dict(sites=length,parsimony_informative_sites=informative,
                            distinct_aligned_sequences=len(set(seqs.values())))


def trim_input(aligned, work):
    """Use short, unique transport IDs so trimming cannot truncate biological IDs."""
    aliases = {f'S{i:07d}':name for i,name in enumerate(sorted(aligned),1)}
    fasta(work/'trimal_input.faa',{alias:aligned[name] for alias,name in aliases.items()})
    write_tsv(work/'trimal_identifiers.tsv',['alias','sequence_id'],
              [dict(alias=a,sequence_id=n) for a,n in aliases.items()])
    return aliases


def restore_trimmed(path, aliases, output):
    seqs = load_sequences(path)
    if set(seqs)!=set(aliases):
        dump(path.with_suffix('.membership.json'),dict(missing_aliases=sorted(set(aliases)-set(seqs)),unexpected_aliases=sorted(set(seqs)-set(aliases))))
        raise ValueError('trimAl changed sequence identifiers/membership; inspect trimmed alias membership report')
    fasta(output,{aliases[a]:seq for a,seq in seqs.items()})


def failure_category(exc, stage):
    if isinstance(exc,subprocess.CalledProcessError):
        return 'tool_failure'
    message=str(exc).lower()
    if 'identifier' in message or 'membership' in message or 'tree tips' in message:
        return 'identifier_mismatch'
    if any(word in message for word in ('gaps','residues','sites','aligned sequences','alignment lengths','column subsequence')):
        return 'alignment_qc_failure'
    return 'processing_failure'


def annotate(treefile, folder, outgroups):
    from Bio import Phylo
    tree = Phylo.read(treefile,'newick')
    ids = {t.name for t in tree.get_terminals()}
    if outgroups:
        if not set(outgroups)<ids:
            raise ValueError('Outgroup IDs must be present and leave at least one ingroup tip')
        tree.root_with_outgroup(*outgroups)
        if not tree.is_monophyletic([t for t in tree.get_terminals() if t.name in outgroups]):
            raise ValueError('Requested outgroup does not form an exclusive clade; rooting needs review')
        Phylo.write(tree,folder/'rooted.treefile','newick')
    _, edges = read_tsv(folder.parent/'orthology_edges.tsv')
    _, events = read_tsv(folder.parent/'duplication_evidence.tsv')
    rows = []
    for i,node in enumerate(tree.get_nonterminals(order='preorder'),1):
        descendants = {t.name for t in node.get_terminals()}
        if descendants == ids:
            continue
        taxa = sorted({n.split('__')[0] for n in descendants})
        label = node.name or (str(node.confidence) if node.confidence is not None else '')
        inside = [e for e in edges if {e['candidate_id'],e['orthologue_id']} <= descendants]
        exact = []
        for j,event in enumerate(events,1):
            a,b = members(event.get('Genes 1','')),members(event.get('Genes 2',''))
            childsets = [{t.name for t in c.get_terminals()} for c in node.clades]
            if outgroups and len(childsets)==2 and a and b and {frozenset(a),frozenset(b)}==set(map(frozenset,childsets)):
                exact.append(str(j))
        assessment = 'unrooted split; no ancestry or duplication timing assigned'
        if outgroups:
            assessment = 'mixed-species clade' if len(taxa)>1 else 'single-species clade; possible expansion, requires reconciliation'
            if exact:
                assessment += '; exact child-set match to supplied OrthoFinder duplication evidence'
        rows.append(dict(branch_id=f'B{i:04d}',support_label=label,tip_count=len(descendants),
            species=','.join(taxa),tip_ids=','.join(sorted(descendants)),
            direct_orthology_edges_within=len(inside),exact_duplication_rows=','.join(exact),interpretation=assessment))
    write_tsv(folder/'branch_evidence.tsv',['branch_id','support_label','tip_count','species','tip_ids',
        'direct_orthology_edges_within','exact_duplication_rows','interpretation'],rows)


def run(args):
    manifest = json.loads((args.project/'manifest.json').read_text())
    groups = [g for g in manifest['groups'] if not args.hog or g['group'] in args.hog]
    if args.hog and set(args.hog)-{g['group'] for g in groups}:
        raise ValueError('Unknown selected group')
    if args.outgroup and len(groups)!=1:
        raise ValueError('Specify exactly one --hog when providing --outgroup')
    if groups and all(g['status']!='ready' for g in groups):
        skipped=[dict(group=g['group'],status='skipped',category='too_few_taxa',stage='preparation',error=g['status'],stderr_log='') for g in groups]
        write_tsv(args.project/'run_summary.tsv',list(skipped[0]),skipped)
        dump(args.project/'run_summary.json',dict(status='complete',counts={'skipped':len(skipped)},groups=skipped))
        print(f'SKIP: {len(skipped)} groups have too few taxa; no tree inference attempted')
        return
    from Bio import Phylo  # Fail before starting expensive computation if unavailable.
    iq_threads=getattr(args,'iqtree_threads',2)
    if args.threads<1 or iq_threads<1 or args.seed<1 or args.min_sites<1 or not 0<=args.max_gap<1:
        raise ValueError('Invalid threads/seed/QC settings')
    tools = {name:shutil.which(value) for name,value in
             [('mafft',args.mafft),('trimal',args.trimal),('iqtree',args.iqtree)]}
    if any(v is None for v in tools.values()):
        raise ValueError(f'Missing executables: {tools}; install environment-phylogeny.yml')
    versions = {}
    for name, flag in [('mafft','--version'),('trimal','--version'),('iqtree','--version')]:
        result = subprocess.run([tools[name],flag],capture_output=True,text=True)
        versions[name] = result.stdout+result.stderr
    statuses = []
    for group in groups:
        folder = args.project/group['group']
        for name, sha in group['files'].items():
            if digest(folder/name)!=sha:
                raise ValueError(f'Prepared input changed: {folder/name}; prepare a new project')
        if group['status']!='ready':
            print(f'SKIP {group["group"]}: {group["status"]}')
            statuses.append(dict(group=group['group'],status='skipped',category='too_few_taxa',stage='preparation',error=group['status'],stderr_log=''))
            continue
        work = folder/'inference'
        if work.exists() and getattr(args,'retry_failed',False):
            previous=json.loads((work/'run.json').read_text())
            if previous['status']=='complete':
                statuses.append(dict(group=group['group'],status='complete',category='previously_complete',stage='complete',error='',stderr_log=''))
                continue
            if previous['status']!='failed':
                raise ValueError(f'{work}: run is not marked failed; inspect before retrying')
            attempt=1
            while (folder/f'inference_attempt_{attempt:03d}').exists():
                attempt+=1
            work.rename(folder/f'inference_attempt_{attempt:03d}')
        work.mkdir(exist_ok=False)
        info = {'versions':versions,'settings':vars(args).copy(),'commands':[], 'status':'running',
                'input_hashes':group['files'],'started_utc':datetime.now(timezone.utc).isoformat()}
        info['settings'] = {k:str(v) if isinstance(v,Path) else v for k,v in info['settings'].items() if k!='func'}
        def execute(command, stdout_name=None):
            info['commands'].append(command)
            dump(work/'run.json',info)
            index = len(info['commands'])
            info['last_stderr_log']=str(work/f'command_{index}.stderr.log')
            with (work/f'command_{index}.stderr.log').open('w') as err, (work/(stdout_name or f'command_{index}.stdout.log')).open('w') as out:
                subprocess.run(command,stdout=out,stderr=err,check=True)
        stage='input'
        try:
            raw = load_sequences(folder/'members.faa')
            if args.outgroup and not set(args.outgroup)<set(raw):
                raise ValueError('Outgroup must be a proper subset of this group')
            method = ['--localpair','--maxiterate','1000'] if args.aligner=='linsi' else ['--ep','0','--genafpair','--maxiterate','1000']
            stage='mafft'
            execute([tools['mafft'],'--amino','--thread',str(args.threads),'--threadit','0',*method,str((folder/'members.faa').resolve())],'aligned.faa')
            stage='alignment_qc'
            aligned, _, before = alignment_qc(work/'aligned.faa',raw,1,0.999999)
            if any(aligned[n].replace('-','')!=raw[n] for n in raw):
                raise ValueError('MAFFT changed residue sequences')
            aliases=trim_input(aligned,work)
            stage='trimal'
            execute([tools['trimal'],'-in',str((work/'trimal_input.faa').resolve()),'-out',str((work/'trimmed_aliases.faa').resolve()),'-automated1'])
            restore_trimmed(work/'trimmed_aliases.faa',aliases,work/'trimmed.faa')
            stage='trimmed_alignment_qc'
            trimmed, qc, after = alignment_qc(work/'trimmed.faa',raw,args.min_sites,args.max_gap)
            order = sorted(raw)
            original_columns = iter(zip(*(aligned[n] for n in order)))
            for column in zip(*(trimmed[n] for n in order)):
                if not any(candidate == column for candidate in original_columns):
                    raise ValueError('Trimmed alignment is not a column subsequence of the original alignment')
            write_tsv(work/'sequence_qc.tsv',list(qc[0]),qc)
            dump(work/'alignment_qc.json',dict(before=before,after=after,retained_fraction=after['sites']/before['sites']))
            if after['distinct_aligned_sequences']<4 or after['parsimony_informative_sites']==0:
                raise ValueError('Fewer than four distinct aligned sequences or no parsimony-informative sites')
            stage='iqtree'
            execute([tools['iqtree'],'-s',str((work/'trimmed.faa').resolve()),'-st','AA','-m','MFP',
                     '-B','1000','--alrt','1000','-bnni','-seed',str(args.seed),'-T',str(iq_threads),
                     '-keep-ident','--prefix',str((work/'gene_tree').resolve())])
            treefile = work/'gene_tree.treefile'
            names = [t.name for t in Phylo.read(treefile,'newick').get_terminals()]
            if len(names)!=len(raw) or set(names)!=set(raw):
                raise ValueError('Tree tips do not match extracted sequences')
            annotate(treefile,work,args.outgroup)
            info['status']='complete'
            info['outputs']={p.name:digest(p) for p in work.iterdir() if p.is_file() and p.name!='run.json'}
        except Exception as exc:
            info['status']='failed'
            info['error']=str(exc)
            info['failure_category']=failure_category(exc,stage)
            info['failed_stage']=stage
            if isinstance(exc,subprocess.CalledProcessError):
                info['exit_code']=exc.returncode
                log=Path(info['last_stderr_log'])
                info['stderr_tail']=log.read_text(errors='replace')[-6000:]
                info['diagnostic_note']='Exit code alone does not establish a resource/thread cause; inspect stderr and tool logs.'
            print(f'FAILED {group["group"]}: {exc}',file=sys.stderr)
        dump(work/'run.json',info)
        statuses.append(dict(group=group['group'],status=info['status'],category=info.get('failure_category',''),
            stage=info.get('failed_stage','complete'),error=info.get('error',''),stderr_log=info.get('last_stderr_log','')))
        print(f'{group["group"]}: {info["status"]}',flush=True)
    write_tsv(args.project/'run_summary.tsv',['group','status','category','stage','error','stderr_log'],statuses)
    dump(args.project/'run_summary.json',dict(status='partial_failure' if any(s['status']=='failed' for s in statuses) else 'complete',
        counts=dict(Counter(s['status'] for s in statuses)),groups=statuses))
    if any(s['status']=='failed' for s in statuses) and not getattr(args,'allow_partial',False):
        raise ValueError('One or more groups failed; inspect run_summary.tsv and per-group run.json/logs. Use --retry-failed after reviewing the cause.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(required=True)
    p = sub.add_parser('prepare')
    for name in ('hogs','orthologues','proteomes','output'):
        p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--duplications',type=Path)
    p.add_argument('--allow-missing-duplications',action='store_true',help='Keep unavailable duplication evidence explicitly unknown')
    p.add_argument('--hog',action='append',help='Repeat to select groups; default all rows')
    p.add_argument('--species-list',type=Path,help='JSON list of species codes for a custom panel')
    family_input=p.add_mutually_exclusive_group()
    family_input.add_argument('--family-map',type=Path,help='Legacy broad TSV mode: unions full orthogroups per family')
    family_input.add_argument('--audit-workbook',type=Path,help='Use only collapsed A/B/C targets by assigned_family, plus linked reference members')
    p.add_argument('--target-species',type=Path,help='panel/targets.json; required for --audit-workbook')
    p.add_argument('--family-column',default='family',help='Family-map column to use (default: family)')
    p.add_argument('--candidate-report',type=Path,help='candidate_report.tsv translating original candidate IDs to prepared IDs')
    p.add_argument('--family',action='append',help='Prepare only this family; repeat to select several')
    p.set_defaults(func=prepare)
    p = sub.add_parser('run')
    p.add_argument('--project',type=Path,required=True)
    p.add_argument('--hog',action='append')
    p.add_argument('--outgroup',action='append',default=[],help='Explicit justified tip ID; repeat for a multi-tip outgroup; one HOG only')
    p.add_argument('--threads',type=int,default=8)
    p.add_argument('--iqtree-threads',type=int,default=2,help='Independent IQ-TREE thread count; exit status 2 does not by itself imply resource exhaustion')
    p.add_argument('--retry-failed',action='store_true',help='Archive failed inference attempts and retry; preserve completed groups')
    p.add_argument('--allow-partial',action='store_true',help='Return success with an explicit partial_failure summary when some groups fail')
    p.add_argument('--seed',type=int,default=20260920)
    p.add_argument('--aligner',choices=['linsi','einsi'],default='linsi')
    p.add_argument('--min-sites',type=int,default=50)
    p.add_argument('--max-gap',type=float,default=0.9)
    p.add_argument('--mafft',default='mafft')
    p.add_argument('--trimal',default='trimal')
    p.add_argument('--iqtree',default='iqtree2')
    p.set_defaults(func=run)
    args=parser.parse_args()
    try:
        args.func(args)
    except (OSError,ValueError,ImportError) as exc:
        parser.exit(1,f'ERROR: {exc}\n')


if __name__=='__main__':
    main()
