"""Discover reference FASTAs and explicit gene annotations; prepare an automatic panel."""
import argparse
import json
import re
from pathlib import Path
from pipeline import records, digest, write_tsv, VERSION
from refseq_gene_map import parse_gff
from prepare_reference import prepare


def build(raw, handoff, out):
    meta=json.loads((handoff/'workflow_manifest.json').read_text())
    targets=meta['target_species']
    files=sorted(p for p in raw.iterdir() if p.is_file() and
        re.search(r'\.(fa|faa|fasta)(\.gz)?$',p.name,re.I))
    if not files:
        raise ValueError(f'No reference FASTAs in {raw}')
    out.mkdir(parents=True,exist_ok=False)
    entries=[]
    for sp in targets:
        source=Path(meta['detection'])/sp/'prepared'/f'{sp}.faa'
        expected=meta['sources'].get(str(source.resolve()))
        if expected is None or digest(source)!=expected:
            raise ValueError(f'Target proteome changed or lacks handoff provenance: {source}')
        entries.append(dict(species=sp,fasta=str(source.resolve()),gene_level=True,source='Handoff-verified complete target proteome'))
    seen=set(targets)
    for source in files:
        code=re.sub(r'\.(fa|faa|fasta)(\.gz)?$','',source.name,flags=re.I)
        if code in seen or not re.fullmatch(r'[A-Za-z][A-Za-z0-9_]*',code) or '__' in code:
            raise ValueError(f'Use a unique safe reference filename/code: {source.name}')
        seen.add(code)
        folder=out/code
        folder.mkdir()
        maps=[p for p in (raw/f'{code}.protein_gene.tsv',raw/f'{code}.tsv') if p.exists()]
        gffs=[p for p in (raw/f'{code}.gff',raw/f'{code}.gff3') if p.exists()]
        if len(maps)>1 or len(gffs)>1:
            raise ValueError(f'Ambiguous mapping/GFF sidecars for {code}')
        rows=list(records(source))
        provenance=dict(fasta=str(source.resolve()),fasta_sha256=digest(source))
        if maps:
            mapping=maps[0]
            provenance.update(mapping=str(mapping.resolve()),mapping_sha256=digest(mapping),method='explicit_tsv')
        else:
            if gffs:
                gene_sets,metadata=parse_gff(gffs[0])
                if any(len(gene_sets.get(h.split()[0],set()))!=1 for h,s in rows):
                    raise ValueError(f'{code}: missing/conflicting CDS protein_id to GeneID mapping; provide an explicit TSV for other GFF schemas')
                mapping_rows=[dict(protein_id=h.split()[0],gene_id='GeneID_'+next(iter(gene_sets[h.split()[0]]))) for h,s in rows]
                provenance.update(gff=str(gffs[0].resolve()),gff_sha256=digest(gffs[0]),annotation_metadata=metadata,method='ncbi_gff')
            else:
                mapping_rows=[]
                for header,seq in rows:
                    pid=header.split()[0]
                    match=re.search(r'(?:^|\s)gene[=:](\S+)',header)
                    gene=match.group(1) if match else pid[len(code)+2:] if pid.startswith(code+'__') else None
                    if not gene:
                        raise ValueError(f'{code}: {pid} lacks an explicit gene ID; supply {code}.protein_gene.tsv or supported GFF')
                    mapping_rows.append(dict(protein_id=pid,gene_id=gene))
                provenance['method']='explicit_header_gene_or_prefixed_gene_level'
            mapping=folder/'protein_gene.tsv'
            write_tsv(mapping,['protein_id','gene_id'],mapping_rows)
        prepared=prepare(code,source,mapping,folder/'prepared')
        (folder/'source.json').write_text(json.dumps(provenance,indent=2)+'\n')
        entries.append(dict(species=code,fasta=str(prepared.resolve()),gene_level=True,source=json.dumps(provenance)))
    panel=out/'panel.auto.json'
    panel.write_text(json.dumps(dict(proteomes=entries),indent=2)+'\n')
    return panel


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--version',action='version',version='adhesome-reference-preparation '+VERSION)
    p.add_argument('--raw',type=Path,required=True)
    p.add_argument('--handoff',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    try:
        print(build(a.raw.resolve(),a.handoff.resolve(),a.out.resolve()))
    except (OSError,ValueError,KeyError) as exc:
        p.exit(1,f'ERROR: {exc}\n')
