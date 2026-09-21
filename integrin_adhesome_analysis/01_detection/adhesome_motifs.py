"""Context-limited motif annotation; short patterns never recruit or assign families."""
import json
from pathlib import Path
import re

def load_library(path):
    data=json.loads(Path(path).read_text(encoding='utf-8'))
    modules={}
    for category in data['integrin_adhesome_network']:
        for entity in category['entities']:
            for module in entity['modules']:
                modules[(entity['name'],module['name'])]=module
    return modules

def scan_motifs(sequence, hits, rule, library):
    rows=[]
    for spec in rule.get('motif_checks',[]):
        key=(spec['entity'],spec['name'])
        if key not in library or library[key]['type']!='Motif':
            raise ValueError('Motif mapping is absent from the supplied library: '+str(key))
        # Human residue numbering and post-translational modifications are not regexes.
        pattern=spec['pattern']
        regex=re.compile(pattern)
        if regex.match(''):
            raise ValueError('Motif patterns must not match empty sequence')
        context=spec.get('context','whole_sequence')
        if context=='c_terminal_window':
            regions=[(max(1,len(sequence)-spec['window']+1),len(sequence))]
        elif context=='domain':
            regions=[(int(h['ali_from']),int(h['ali_to'])) for h in hits if h['pfam_accession'] in spec['domains']]
        elif context=='whole_sequence':
            regions=[(1,len(sequence))]
        else:
            raise ValueError('Unknown motif context: '+context)
        matches=set()
        for start,end in regions:
            for m in regex.finditer(sequence[start-1:end]):
                matches.add((start+m.start(),start+m.end()-1,m.group()))
        if matches:
            for start,end,matched in sorted(matches):
                rows.append(dict(motif=spec.get('label',spec['name']),pattern=pattern,context=context,state='detected',
                    start=start,end=end,sequence=matched,interpretation=spec.get('interpretation','supporting_only')))
        else:
            rows.append(dict(motif=spec.get('label',spec['name']),pattern=pattern,context=context,
                state='not_detected' if regions else 'not_tested_missing_domain',start='',end='',sequence='',
                interpretation=spec.get('interpretation','supporting_only')))
    return rows

def noncatalytic_alignment(support, hits, spec):
    """Conservative span outside a catalytic domain; not an exact gapped alignment count."""
    intervals=[(int(h['ali_from']),int(h['ali_to'])) for h in hits if h['pfam_accession'] in spec['domains']]
    if not intervals:
        return False
    for row in support:
        if row.get('similarity')!='true' or row.get('reciprocal')!='true':
            continue
        if not row.get('subject_start') or not row.get('subject_end'):
            continue
        start,end=sorted((int(row['subject_start']),int(row['subject_end'])))
        # Count only the span after every catalytic domain, avoiding receptor/PTP repeats.
        if end-max(start-1,max(b for a,b in intervals)) >= spec['minimum_residues']:
            return True
    return False
