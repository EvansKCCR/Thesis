#!/usr/bin/env python3
"""Auditable per-HOG MAFFT -> trimAl -> IQ-TREE workflow (Python >=3.11)."""
import argparse
from collections import Counter
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


def prepare(args):
    configure_species(getattr(args, 'species_list', None))
    fields, rows = read_tsv(args.hogs)
    active_species = table_species(fields)
    files = [args.hogs,args.orthologues,args.duplications]+[args.proteomes/f'{s}.faa' for s in active_species]
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
    if not rows or len({r[key] for r in rows}) != len(rows):
        raise ValueError('Empty selection or duplicate group IDs')
    sequences = {}
    for sp in active_species:
        seqs = load_sequences(args.proteomes/f'{sp}.faa')
        if any(not n.startswith(sp+'__') for n in seqs):
            raise ValueError(f'Incorrect species prefix in {sp}.faa')
        sequences.update(seqs)
    _, orthologues = read_tsv(args.orthologues)
    dupfields, duplications = read_tsv(args.duplications)
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
    args.output.mkdir(parents=True, exist_ok=False)
    manifest = {'created_utc':datetime.now(timezone.utc).isoformat(),
        'inputs':{str(p.resolve()):digest(p) for p in files}, 'species':active_species, 'groups':[]}
    summary = []
    for group,row,ids in plans:
        folder = args.output/group
        folder.mkdir()
        chosen = {n:sequences[n] for n in ids}
        fasta(folder/'members.faa',chosen)
        matched = members(row.get('Matched_candidates',''))
        tips = [dict(sequence_id=n,species=n.split('__')[0],length=len(sequences[n]),
                     candidate=int(n in matched)) for n in sorted(ids)]
        write_tsv(folder/'tips.tsv',['sequence_id','species','length','candidate'],tips)
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
        status = 'ready' if len(ids)>=4 else 'too_few_taxa_for_supported_tree'
        entry = dict(group=group,sequences=len(ids),unique_sequences=len(set(chosen.values())),
                     candidates=len(matched & ids),orthology_edges=len(edges),duplication_rows=len(events),status=status)
        summary.append(entry)
        manifest['groups'].append(dict(entry,files={p.name:digest(p) for p in folder.iterdir()}))
    write_tsv(args.output/'groups.tsv',list(summary[0]),summary)
    dump(args.output/'manifest.json',manifest)
    print(f'Prepared {len(plans)} separate groups in {args.output}; '+
          f'{sum(r["status"]=="ready" for r in summary)} have at least four taxa.')


def alignment_qc(path, expected, min_sites, max_gap):
    seqs = load_sequences(path)
    if set(seqs) != set(expected):
        raise ValueError(f'{path}: alignment changed sequence identifiers')
    lengths = {len(s) for s in seqs.values()}
    if len(lengths)!=1:
        raise ValueError(f'{path}: unequal alignment lengths')
    length = lengths.pop()
    if length < min_sites:
        raise ValueError(f'{path}: only {length} sites; minimum is {min_sites}')
    rows = [dict(sequence_id=n,sites=length,gap_fraction=s.count('-')/length,
                 informative_residues=sum(c in 'ACDEFGHIKLMNPQRSTVWY' for c in s)) for n,s in seqs.items()]
    if any(r['gap_fraction']>max_gap or r['informative_residues']==0 for r in rows):
        raise ValueError(f'{path}: excessive gaps or no standard residues in a sequence; review alignment')
    informative = 0
    for column in zip(*seqs.values()):
        counts = Counter(c for c in column if c in 'ACDEFGHIKLMNPQRSTVWY')
        informative += sum(v>=2 for v in counts.values())>=2
    return seqs, rows, dict(sites=length,parsimony_informative_sites=informative,
                            distinct_aligned_sequences=len(set(seqs.values())))


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
    from Bio import Phylo  # Fail before starting expensive computation if unavailable.
    manifest = json.loads((args.project/'manifest.json').read_text())
    groups = [g for g in manifest['groups'] if not args.hog or g['group'] in args.hog]
    if args.hog and set(args.hog)-{g['group'] for g in groups}:
        raise ValueError('Unknown selected group')
    if args.outgroup and len(groups)!=1:
        raise ValueError('Specify exactly one --hog when providing --outgroup')
    if args.threads<1 or args.seed<1 or args.min_sites<1 or not 0<=args.max_gap<1:
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
            continue
        work = folder/'inference'
        work.mkdir(exist_ok=False)
        info = {'versions':versions,'settings':vars(args).copy(),'commands':[], 'status':'running',
                'input_hashes':group['files'],'started_utc':datetime.now(timezone.utc).isoformat()}
        info['settings'] = {k:str(v) if isinstance(v,Path) else v for k,v in info['settings'].items() if k!='func'}
        def execute(command, stdout_name=None):
            info['commands'].append(command)
            dump(work/'run.json',info)
            index = len(info['commands'])
            with (work/f'command_{index}.stderr.log').open('w') as err, (work/(stdout_name or f'command_{index}.stdout.log')).open('w') as out:
                subprocess.run(command,stdout=out,stderr=err,check=True)
        try:
            raw = load_sequences(folder/'members.faa')
            if args.outgroup and not set(args.outgroup)<set(raw):
                raise ValueError('Outgroup must be a proper subset of this group')
            method = ['--localpair','--maxiterate','1000'] if args.aligner=='linsi' else ['--ep','0','--genafpair','--maxiterate','1000']
            execute([tools['mafft'],'--amino','--thread',str(args.threads),'--threadit','0',*method,str((folder/'members.faa').resolve())],'aligned.faa')
            aligned, _, before = alignment_qc(work/'aligned.faa',raw,1,0.999999)
            if any(aligned[n].replace('-','')!=raw[n] for n in raw):
                raise ValueError('MAFFT changed residue sequences')
            execute([tools['trimal'],'-in',str((work/'aligned.faa').resolve()),'-out',str((work/'trimmed.faa').resolve()),'-automated1'])
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
            execute([tools['iqtree'],'-s',str((work/'trimmed.faa').resolve()),'-st','AA','-m','MFP',
                     '-B','1000','--alrt','1000','-bnni','-seed',str(args.seed),'-T',str(args.threads),
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
            print(f'FAILED {group["group"]}: {exc}',file=sys.stderr)
        dump(work/'run.json',info)
        statuses.append(dict(group=group['group'],status=info['status']))
        print(f'{group["group"]}: {info["status"]}',flush=True)
    if any(s['status']=='failed' for s in statuses):
        raise ValueError('One or more groups failed; inspect per-group run.json and logs')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(required=True)
    p = sub.add_parser('prepare')
    for name in ('hogs','orthologues','duplications','proteomes','output'):
        p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--hog',action='append',help='Repeat to select groups; default all rows')
    p.add_argument('--species-list',type=Path,help='JSON list of species codes for a custom panel')
    p.set_defaults(func=prepare)
    p = sub.add_parser('run')
    p.add_argument('--project',type=Path,required=True)
    p.add_argument('--hog',action='append')
    p.add_argument('--outgroup',action='append',default=[],help='Explicit justified tip ID; repeat for a multi-tip outgroup; one HOG only')
    p.add_argument('--threads',type=int,default=8)
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
