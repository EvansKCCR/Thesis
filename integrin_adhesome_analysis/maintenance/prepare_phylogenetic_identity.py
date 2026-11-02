#!/usr/bin/env python3
"""Reuse existing trees; report pure supported reference clades, never infer new trees."""

# Resolve source resources independently of the caller working directory.
import sys as _layout_sys
from pathlib import Path as _LayoutPath
_LAYOUT_BUNDLE_ROOT = next(p for p in _LayoutPath(__file__).resolve().parents
                          if (p/'config/layout_paths.json').is_file() and (p/'pipeline_paths.py').is_file())
if str(_LAYOUT_BUNDLE_ROOT) not in _layout_sys.path:
    _layout_sys.path.insert(0, str(_LAYOUT_BUNDLE_ROOT))
from pipeline_paths import activate as _layout_activate, resource_path as _layout_resource
_layout_activate()

import sys
from pathlib import Path
ROOT=_LAYOUT_BUNDLE_ROOT;sys.path.insert(0,str(ROOT))
from pipeline_paths import activate
activate()
import argparse
import csv
import hashlib
import json
from phylogenetic_identity import project_inventory,sha


def prepare(project,out,reference_labels=None):
    inventory=project_inventory(project);out=_layout_resource(Path(out));out.mkdir(parents=True,exist_ok=False)
    (_layout_resource(out/'phylogenetic_integration.json')).write_text(json.dumps(inventory,indent=2)+'\n')
    fields='protein_id sequence_sha256 family clade_id tree_file tree_sha256 sh_alrt ufboot reference_count reference_taxa_count reference_label_source placement_status duplication_status'.split()
    evidence=[]
    if reference_labels:
        from Bio import Phylo
        with _layout_resource(Path(reference_labels)).open(encoding='utf-8-sig',newline='') as h:labels=list(csv.DictReader(h,delimiter='\t'))
        if any(not all(r.get(k) for k in ('sequence_id','family','species','source')) for r in labels):
            raise ValueError('Reference labels require sequence_id, family, species and source')
        if len({r['sequence_id'] for r in labels})!=len(labels):raise ValueError('Duplicate reference tip label')
        labels={r['sequence_id']:r for r in labels}
        manifest=json.loads((_layout_resource(project/'manifest.json')).read_text())
        for group in manifest['groups']:
            folder=project/group['group'];tree_path=_layout_resource(folder/'inference/gene_tree.treefile')
            if not _layout_resource(tree_path).is_file():continue
            for name in ('tips.tsv','members.faa'):
                if sha(folder/name)!=group['files'][name]:raise ValueError('Prepared phylogeny input changed: '+str(folder/name))
            with (_layout_resource(folder/'tips.tsv')).open(newline='') as h:tips={r['sequence_id']:r for r in csv.DictReader(h,delimiter='\t')}
            sequences={};key=None
            for line in (_layout_resource(folder/'members.faa')).read_text().splitlines():
                if line.startswith('>'):key=line[1:].split()[0];sequences[key]=''
                elif key:sequences[key]+=line.strip().upper()
            tree=Phylo.read(tree_path,'newick')
            if {t.name for t in tree.get_terminals()}!=set(tips):continue
            for pid,tip in tips.items():
                if tip['candidate']!='1':continue
                candidates=[]
                for index,clade in enumerate(tree.get_nonterminals()):
                    names={t.name for t in clade.get_terminals()}
                    if pid not in names or len(names)==len(tips):continue
                    refs=names-{i for i,t in tips.items() if t['candidate']=='1'}
                    if len(refs)<2 or not refs<=labels.keys():continue
                    families={labels[i]['family'] for i in refs};taxa={labels[i]['species'] for i in refs}
                    support=str(clade.name or clade.confidence or '').split('/')
                    if len(families)!=1 or len(taxa)<2 or len(support)!=2:continue
                    try:alrt,boot=map(float,support)
                    except ValueError:continue
                    if alrt<80 or boot<95:continue
                    candidates.append((len(names),dict(protein_id=pid,sequence_sha256=hashlib.sha256(sequences[pid].encode()).hexdigest(),
                        family=next(iter(families)),clade_id=group['group']+':'+str(index),tree_file=str(_layout_resource(tree_path).resolve()),
                        tree_sha256=sha(tree_path),sh_alrt=alrt,ufboot=boot,reference_count=len(refs),reference_taxa_count=len(taxa),
                        reference_label_source=sha(reference_labels),placement_status='pure_supported_clade',duplication_status='not_reconciled')))
                if candidates:evidence.append(min(candidates,key=lambda pair:pair[0])[1])
        # A candidate occurring in more than one tree must agree; conflicts remain unqualified.
        from collections import defaultdict
        by_id=defaultdict(list)
        for row in evidence:by_id[row['protein_id']].append(row)
        evidence=[rows[0] for rows in by_id.values() if len({r['family'] for r in rows})==1]
        inventory.update(reference_identity_labels_supplied=True,qualifying_candidates=len(evidence),
            orthology_interpretation='Supported unrooted clade concordance only; rooting/duplication reconciliation and paralog names remain unresolved')
    (_layout_resource(out/'phylogenetic_integration.json')).write_text(json.dumps(inventory,indent=2)+'\n')
    with (_layout_resource(out/'phylogenetic_identity_evidence.tsv')).open('w',newline='') as h:
        w=csv.DictWriter(h,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(evidence)
    return inventory


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--project',type=Path,required=True)
    p.add_argument('--reference-labels',type=Path);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();print(json.dumps(prepare(a.project,a.out,a.reference_labels),indent=2))


if __name__=='__main__':main()
