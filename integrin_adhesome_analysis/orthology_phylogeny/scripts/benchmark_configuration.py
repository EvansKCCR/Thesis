"""Version 1.0: explicit benchmark scope and user-defined acceptance criteria."""
import json
import math
from pathlib import Path
import re

DEFAULT_POLICY=Path(__file__).resolve().parents[1]/'config/benchmark_policy.json'


def quality_policy(path=None):
    path=Path(path) if path else DEFAULT_POLICY
    policy=json.loads(path.read_text())
    if not isinstance(policy,dict):raise ValueError('Acceptance policy must be a JSON object')
    if policy.get('metric')!='open_orthobench_weighted_pairwise':
        raise ValueError('Acceptance policy must use open_orthobench_weighted_pairwise')
    criteria={}
    for key,metric in [('minimum_precision','precision'),('minimum_f1','f1'),
                       ('minimum_recall','recall'),('minimum_exact_fraction','exact_reconstruction_fraction')]:
        threshold=policy.get(key)
        if threshold is not None:
            if isinstance(threshold,bool) or not isinstance(threshold,(float,int)) or not math.isfinite(threshold) or not 0<=threshold<=1:
                raise ValueError('Acceptance thresholds must be finite numbers between zero and one: '+key)
            criteria[metric]=float(threshold)
    if not criteria:
        raise ValueError('Acceptance policy must specify at least one minimum')
    return path,criteria


def acceptance(metrics, criteria):
    checks={name:dict(observed=metrics[name],minimum=minimum,passed=metrics[name]>=minimum)
            for name,minimum in criteria.items()}
    return dict(status='accepted' if all(check['passed'] for check in checks.values()) else 'rejected',checks=checks)


def labelled_clades(path):
    """Parse the safe species-code Newick produced by OrthoFinder, without a new dependency."""
    tokens=re.findall(r'[(),:;]|[^(),:;\s]+',Path(path).read_text())
    cursor=0;clades={};leaves=set()
    def label():
        nonlocal cursor
        if cursor>=len(tokens) or tokens[cursor] in '(),:;':
            raise ValueError('Malformed labelled species tree')
        result=tokens[cursor];cursor+=1
        if not re.fullmatch(r'[A-Za-z][A-Za-z0-9_]*',result):
            raise ValueError('Species tree requires safe explicit species/node labels')
        return result
    def node():
        nonlocal cursor
        if cursor>=len(tokens):raise ValueError('Incomplete species tree')
        if tokens[cursor]=='(':
            cursor+=1;members=node()
            while cursor<len(tokens) and tokens[cursor]==',':
                cursor+=1;members|=node()
            if cursor>=len(tokens) or tokens[cursor]!=')':raise ValueError('Malformed species tree')
            cursor+=1;name=label()
            if not re.fullmatch(r'N\d+',name) or name in clades:raise ValueError('Expected unique OrthoFinder N-node labels')
            clades[name]=set(members)
        else:
            name=label()
            if name in leaves:raise ValueError('Duplicate species-tree leaf: '+name)
            leaves.add(name);members={name}
        if cursor<len(tokens) and tokens[cursor]==':':
            cursor+=1
            if cursor>=len(tokens):raise ValueError('Missing species-tree branch length')
            length=float(tokens[cursor]);cursor+=1
            if not math.isfinite(length):raise ValueError('Invalid species-tree branch length')
        return members
    node()
    if cursor<len(tokens) and tokens[cursor]==';':cursor+=1
    if cursor!=len(tokens):raise ValueError('Unexpected species-tree content')
    return clades,leaves


def prediction_source(results,level='root',cohort_species=None):
    results=Path(results);log=results/'Log.txt';inputs=[]
    version=None
    if log.is_file():
        found=re.search(r'OrthoFinder version\s+(\d+)\.(\d+)\.(\d+)',log.read_text(errors='replace'))
        if found:version=tuple(map(int,found.groups()))
        inputs.append(log)
    root_options=[results/'Phylogenetic_Hierarchical_Orthogroups/N0.tsv',results/'Orthogroups/Orthogroups.tsv']
    if version and version>=(3,1,5):root_options.reverse()
    node='N0'
    if level=='cohort-lca':
        tree=results/'Species_Tree/SpeciesTree_rooted_node_labels.txt'
        if not tree.is_file() or not cohort_species:
            raise ValueError('cohort-lca requires the labelled OrthoFinder species tree and mapped benchmark species')
        clades,leaves=labelled_clades(tree);cohort=set(cohort_species)
        if not cohort<=leaves:raise ValueError('Benchmark species are missing from the labelled species tree')
        eligible=[name for name,members in clades.items() if cohort<=members]
        if not eligible:raise ValueError('No species-tree node contains the complete benchmark cohort')
        node=min(eligible,key=lambda name:len(clades[name]));inputs.append(tree)
        choices=root_options if node=='N0' else [results/'Phylogenetic_Hierarchical_Orthogroups'/(node+'.tsv')]
    elif level=='root':choices=root_options
    else:raise ValueError('Unknown prediction level: '+str(level))
    path=next((p for p in choices if p.is_file()),None)
    if path is None:
        raise ValueError('Complete orthogroup table unavailable for '+level+' ('+node+'); no candidate-filtered or alternate-clade fallback is permitted')
    return path,dict(requested_level=level,selected_node=node,
        orthofinder_version='.'.join(map(str,version)) if version else 'unknown',
        cohort_species=sorted(cohort_species or []),selection_basis='species-tree cohort membership; never benchmark score'),inputs
