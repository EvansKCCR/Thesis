#!/usr/bin/env python3
"""Score external gold/reference benchmark TSV after reference profiles and feature rows are populated.

Input must have architecture, length, expected_subfamily, split,
reference_exclusion_verified, and original gold provenance. Do not derive
reference_subfamily, orthology labels, or profile winners from expected_subfamily.
"""
import argparse,csv,json,sys
from collections import Counter,defaultdict
from pathlib import Path
sys.path.insert(0,str(Path(__file__).parent))
from integrated_adhesome_pipeline import audit_record,boolean


def ratio(n,d):return None if not d else round(n/d,4)

def evaluate(data,rules):
    assert data,'Empty benchmark subset'
    cm=Counter();g=defaultdict(lambda:Counter());faults=[];n_assigned=0;n_right=0;n_abstain=0
    for r in data:
        true=r['expected_subfamily'].strip();pred=audit_record(r,rules)
        called=pred['assigned_subfamily'] or pred['assigned_family']
        named=called not in ('unassigned','')
        correct=named and called==true
        n_assigned+=bool(named);n_right+=bool(correct);n_abstain+=not named
        cm[(true,called)]+=1
        group=r.get('expected_superfamily') or 'other';g[group]['n']+=1
        g[group]['correct']+=bool(correct);g[group]['assigned']+=bool(named)
        if not correct:faults.append({'species':r.get('species'),'symbol':r.get('reviewed_symbol'), 'expected':true,'called':called, 'state':pred['classification_state'],'notes':pred['rationale']})
    included={t for t,p in cm if any(tt==t and pp==t for tt,pp in cm)}
    allclasses={t for t,p in cm if t in rules['families'] and rules['families'][t].get('classification_level')=='subfamily'}
    def metrics_for(lbl):
        tp=cm[(lbl,lbl)]
        fn=sum(n for (t,p),n in cm.items() if t==lbl and p!=lbl)
        fp=sum(n for (t,p),n in cm.items() if t!=lbl and p==lbl)
        return {'n_true':tp+fn,'tp':tp,'fp':fp,'fn':fn,'precision':ratio(tp,tp+fp),'recall':ratio(tp,tp+fn),'f1':ratio(2*tp,2*tp+fp+fn)}
    classes={x:metrics_for(x) for x in sorted(allclasses)}
    nonempty_f1=[v['f1'] for v in classes.values() if v['f1'] is not None]
    other=[r for r in data if r.get('split')=='hard_negative_reference']
    other_err=sum(audit_record(r,rules)['assigned_subfamily'] not in ('',None) for r in other)
    return {'n':len(data),'coverage':ratio(n_assigned,len(data)),'abstention':ratio(n_abstain,len(data)),
        'exact_accuracy_with_abstention':ratio(n_right,len(data)),
        'selective_precision_assigned_only':ratio(n_right,n_assigned),
        'macro_subfamily_f1':round(sum(nonempty_f1)/len(nonempty_f1),4) if nonempty_f1 else None,
        'hard_negative_subfamily_false_assignment_rate':ratio(other_err,len(other)),
        'by_group':{x:dict(y) for x,y in sorted(g.items())},'by_class':classes,
        'errors':faults,'confusion':[{'true':x,'called':y,'count':z} for (x,y),z in sorted(cm.items())]}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--version',action='version',version='adhesome-benchmark 1.0')
    p.add_argument('filled_gold_tsv',type=Path)
    p.add_argument('--rules',type=Path,default=Path(__file__).resolve().parents[1]/'01_detection/config/family_rules.json')
    p.add_argument('--split',action='append',default=['phylogenetic_holdout','hard_negative_reference'])
    p.add_argument('--out',type=Path,default=Path('adhesome_benchmark_metrics.json'))
    args=p.parse_args()
    with args.filled_gold_tsv.open(newline='') as h:rows=list(csv.DictReader(h,delimiter='\t'))
    data=[r for r in rows if r.get('split') in set(args.split)]
    missing=[(r.get('reviewed_symbol'),r.get('split')) for r in data if not r.get('architecture') or not r.get('uniprot_accession') or not boolean(r.get('reference_exclusion_verified')) or not r.get('feature_generation_provenance')]
    if missing:raise ValueError(f'Benchmark invalid: {len(missing)} entries lack domain evidence, accession, excluded-test-set reference verification, or feature provenance. Example: {missing[:3]}')
    rules=json.loads(args.rules.read_text())
    result=evaluate(data,rules)
    result['pipeline_version']='1.0'
    args.out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('errors','confusion','by_class')},indent=2))

if __name__=='__main__':main()
