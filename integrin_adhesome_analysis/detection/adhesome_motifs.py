"""Current motif definitions; functional hits remain sequence candidates."""
import json
from pathlib import Path
import re


def load_library(path):
    data=json.loads(Path(path).read_text(encoding="utf-8"))
    if "rules" not in data:raise ValueError("Use the current evidence-based motif snapshot")
    return {"current_rules":{r["rule_id"]:r for r in data["rules"]}}


def scan_motifs(sequence, hits, rule, library):
    rows=[]
    if not rule.get('motif_rule_ids'):return rows
    if 'current_rules' not in library:raise ValueError('Use the current evidence-based motif snapshot')
    if 'current_rules' in library:
        for mid in rule.get('motif_rule_ids',[]):
            if mid not in library['current_rules']:raise ValueError('Unknown current motif rule: '+mid)
            spec=library['current_rules'][mid];regex=re.compile(spec['regex_python'])
            if regex.match(''):raise ValueError('Motif patterns must not match empty sequence')
            for pos in range(len(sequence)):
                match=regex.match(sequence,pos)
                if match and match.end()>pos and set(match.group())<=set('ACDEFGHIKLMNPQRSTVWY'):
                    rows.append(dict(motif=mid,pattern=spec['regex_python'],context='full_sequence_context_unverified',
                        state='detected',start=pos+1,end=match.end(),sequence=match.group(),interpretation='sequence_candidate_only'))
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
