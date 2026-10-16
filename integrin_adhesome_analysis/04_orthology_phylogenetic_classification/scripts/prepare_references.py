"""Discover reference FASTAs and explicit gene annotations; prepare an automatic panel."""
import argparse
import json
import re
from pathlib import Path
from pipeline import records, digest, write_tsv, VERSION
from refseq_gene_map import parse_gff
from prepare_reference import prepare, BENCHMARK_SEQUENCE_POLICY
from benchmark_workflow import bundle, read_groups, input_proteomes


def build(raw, handoff, out, benchmark_refogs=None, benchmark_input_dir=None):
    meta=json.loads((handoff/'workflow_manifest.json').read_text())
    targets=meta['target_species']
    files=sorted(p for p in raw.iterdir() if p.is_file() and
        re.search(r'\.(fa|faa|fasta)(\.gz)?$',p.name,re.I))
    refogs,input_dir=bundle(raw,benchmark_refogs,benchmark_input_dir)
    benchmark_files=set()
    if refogs is not None:
        truth,uncertain,truth_paths=read_groups(refogs)
        benchmark_files={p for p,seqs in input_proteomes(input_dir,truth)}
        benchmark_hashes={digest(p) for p in benchmark_files}
        files=[p for p in files if p in benchmark_files or digest(p) not in benchmark_hashes]
        files=sorted(set(files)|benchmark_files)
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
        is_benchmark=source in benchmark_files
        if is_benchmark:
            code='Bench_'+re.sub(r'[^A-Za-z0-9_]','_',code)
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
        panel_metadata=raw/f'{code}.reference.json'
        if panel_metadata.exists():
            scope=json.loads(panel_metadata.read_text(encoding='utf-8'))
            if scope.get('fasta_sha256')!=provenance['fasta_sha256']:
                raise ValueError(f'{code}: reference-scope metadata does not match the FASTA')
            provenance.update(reference_scope=scope,reference_scope_file=str(panel_metadata.resolve()),
                              reference_scope_sha256=digest(panel_metadata))
        if is_benchmark:
            # Preserve the supplied protein benchmark's units. Do not collapse
            # benchmark isoforms or claim that protein IDs are validated genes.
            mapping=folder/'protein_gene.tsv'
            write_tsv(mapping,['protein_id','gene_id'],[dict(protein_id=h.split()[0],gene_id=h.split()[0]) for h,s in rows])
            provenance.update(method='original_benchmark_protein_units',benchmark_refogs=str(refogs.resolve()),
                              benchmark_truth_sha256={str(p.relative_to(refogs)):digest(p) for p in truth_paths},
                              unit_level='benchmark_protein',sequence_policy=BENCHMARK_SEQUENCE_POLICY)
        elif maps:
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
        prepared=prepare(code,source,mapping,folder/'prepared',internal_stop_policy='mask' if is_benchmark else 'reject')
        preparation=folder/'prepared/preparation.json'
        provenance.update(preparation=str(preparation.resolve()),preparation_sha256=digest(preparation),
                          sequence_qc=json.loads(preparation.read_text())['sequence_qc'])
        if is_benchmark and provenance['sequence_qc']['internal_stop_proteins']:
            qc=provenance['sequence_qc']
            print(f"{code}: masked {qc['internal_stop_sites']} stop markers as X in {qc['internal_stop_proteins']} proteins; see {folder/'prepared/sequence_qc.tsv'}",flush=True)
        (folder/'source.json').write_text(json.dumps(provenance,indent=2)+'\n')
        entries.append(dict(species=code,fasta=str(prepared.resolve()),gene_level=not is_benchmark,
                            unit_level='benchmark_protein' if is_benchmark else 'gene',source=json.dumps(provenance)))
    panel=out/'panel.auto.json'
    panel.write_text(json.dumps(dict(proteomes=entries),indent=2)+'\n')
    return panel


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--version',action='version',version='adhesome-reference-preparation '+VERSION)
    p.add_argument('--raw',type=Path,required=True)
    p.add_argument('--handoff',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--benchmark-refogs',type=Path)
    p.add_argument('--benchmark-input-dir',type=Path)
    a=p.parse_args()
    try:
        print(build(a.raw.resolve(),a.handoff.resolve(),a.out.resolve(),a.benchmark_refogs,a.benchmark_input_dir))
    except (OSError,ValueError,KeyError) as exc:
        p.exit(1,f'ERROR: {exc}\n')
