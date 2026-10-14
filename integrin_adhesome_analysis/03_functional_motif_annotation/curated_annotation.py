"""Current full-sequence motif matching and sequence-verified topology import."""
import re
import zipfile
from pathlib import Path
AA=set("ACDEFGHIKLMNPQRSTVWY")


def occurrences(regex,seq):
    # match(pos) preserves full-sequence terminal anchors and original group numbering.
    for i in range(len(seq)):
        m=regex.match(seq,i)
        if m and m.end()>i and set(m.group())<=AA:
            yield i+1,m.end(),m.group()


def import_topology(paths,families,regions):
    seqs={(f,r['id']):r['sequence'] for f,rs in families.items() for r in rs}
    predictions={}
    for path in paths:
        if Path(path).suffix.lower()=='.zip':
            with zipfile.ZipFile(path) as archive:
                members=[n for n in archive.namelist() if n.endswith('predicted_topologies.3line')]
                if not members: raise ValueError(f'No topology records in {path}')
                sources=[(str(Path(path).resolve())+'!'+n,archive.read(n).decode('utf-8')) for n in members]
        else:
            sources=[(str(Path(path).resolve()),Path(path).read_text())]
        for source,text in sources:
            lines=[l.strip() for l in text.splitlines() if l.strip()]
            if len(lines)%3: raise ValueError(f'Invalid 3line topology file: {source}')
            for i in range(0,len(lines),3):
                header,seq,top=lines[i:i+3]
                if not header.startswith('>') or len(seq)!=len(top) or set(top)-set('IOMSBX'):
                    raise ValueError(f'Invalid topology record in {source}')
                key=header[1:].split()[0]
                value=(seq,top)
                if key in predictions and predictions[key][:2]!=value: raise ValueError('Conflicting topology: '+key)
                predictions[key]=(seq,top,source)
    audit=[]
    for (family,ident),seq in seqs.items():
        pred=predictions.get(ident)
        state='missing' if pred is None else 'sequence_mismatch' if pred[0]!=seq else 'exact_sequence_verified'
        audit.append(dict(family=family,sequence_id=ident,status=state,source=pred[2] if pred else ''))
        if state!='exact_sequence_verified': continue
        top=pred[1]
        for m in re.finditer(r'(I+|O+|M+|S+|B+|X+)',top):
            label={'I':'predicted_inside','O':'predicted_outside','M':'transmembrane','B':'transmembrane','S':'signal_peptide','X':'unknown_topology'}[m.group()[0]]
            regions[family,ident].append((label,m.start()+1,m.end()))
        # Tail only for a single membrane-spanning helix followed by terminal inside residues.
        tm=list(re.finditer('M+',top))
        tail=re.search('I+$',top)
        if len(tm)==1 and tail and tm[0].end()==tail.start():
            regions[family,ident].append(('cytoplasmic_tail',tail.start()+1,len(seq)))
    return audit
