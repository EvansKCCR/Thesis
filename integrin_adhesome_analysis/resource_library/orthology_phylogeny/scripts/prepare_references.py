"""Discover reference FASTAs and explicit gene annotations; prepare an automatic panel."""

# Resolve source resources independently of the caller working directory.
import sys as _layout_sys
from pathlib import Path as _LayoutPath
_LAYOUT_BUNDLE_ROOT = next(p for p in _LayoutPath(__file__).resolve().parents
                          if (p/'config/layout_paths.json').is_file() and (p/'pipeline_paths.py').is_file())
if str(_LAYOUT_BUNDLE_ROOT) not in _layout_sys.path:
    _layout_sys.path.insert(0, str(_LAYOUT_BUNDLE_ROOT))
from pipeline_paths import activate as _layout_activate, resource_path as _layout_resource
_layout_activate()

import argparse
import json
import re
from pathlib import Path
from pipeline import records, digest, write_tsv, VERSION
from refseq_gene_map import parse_gff
from prepare_reference import prepare, BENCHMARK_SEQUENCE_POLICY
from benchmark_workflow import bundle, read_groups, input_proteomes


def build(raw, handoff, out, benchmark_refogs=None, benchmark_input_dir=None):
    meta=json.loads((_layout_resource(handoff/'workflow_manifest.json')).read_text())
    targets=meta['target_species']
    files=sorted(p for p in raw.iterdir() if _layout_resource(p).is_file() and
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
        source=_layout_resource(Path(meta['detection'])/sp/'prepared'/f'{sp}.faa')
        expected=meta['sources'].get(str(_layout_resource(source).resolve()))
        if expected is None or digest(source)!=expected:
            raise ValueError(f'Target proteome changed or lacks handoff provenance: {source}')
        control_unit=meta.get('query_unit')=='control_protein'
        entries.append(dict(species=sp,fasta=str(_layout_resource(source).resolve()),gene_level=not control_unit,
            unit_level='control_protein' if control_unit else 'gene',
            source='Handoff-verified labelled protein units; not a complete proteome' if control_unit else 'Handoff-verified complete target proteome'))
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
        maps=[p for p in (_layout_resource(raw/f'{code}.protein_gene.tsv'),_layout_resource(raw/f'{code}.tsv')) if _layout_resource(p).exists()]
        gffs=[p for p in (_layout_resource(raw/f'{code}.gff'),_layout_resource(raw/f'{code}.gff3')) if _layout_resource(p).exists()]
        if len(maps)>1 or len(gffs)>1:
            raise ValueError(f'Ambiguous mapping/GFF sidecars for {code}')
        rows=list(records(source))
        provenance=dict(fasta=str(_layout_resource(source).resolve()),fasta_sha256=digest(source))
        panel_metadata=_layout_resource(raw/f'{code}.reference.json')
        if _layout_resource(panel_metadata).exists():
            scope=json.loads(_layout_resource(panel_metadata).read_text(encoding='utf-8'))
            if scope.get('fasta_sha256')!=provenance['fasta_sha256']:
                raise ValueError(f'{code}: reference-scope metadata does not match the FASTA')
            provenance.update(reference_scope=scope,reference_scope_file=str(_layout_resource(panel_metadata).resolve()),
                              reference_scope_sha256=digest(panel_metadata))
        if is_benchmark:
            # Preserve the supplied protein benchmark's units. Do not collapse
            # benchmark isoforms or claim that protein IDs are validated genes.
            mapping=_layout_resource(folder/'protein_gene.tsv')
            write_tsv(mapping,['protein_id','gene_id'],[dict(protein_id=h.split()[0],gene_id=h.split()[0]) for h,s in rows])
            provenance.update(method='original_benchmark_protein_units',benchmark_refogs=str(_layout_resource(refogs).resolve()),
                              benchmark_truth_sha256={str(p.relative_to(refogs)):digest(p) for p in truth_paths},
                              unit_level='benchmark_protein',sequence_policy=BENCHMARK_SEQUENCE_POLICY)
        elif maps:
            mapping=maps[0]
            provenance.update(mapping=str(_layout_resource(mapping).resolve()),mapping_sha256=digest(mapping),method='explicit_tsv')
        else:
            if gffs:
                gene_sets,metadata=parse_gff(gffs[0])
                if any(len(gene_sets.get(h.split()[0],set()))!=1 for h,s in rows):
                    raise ValueError(f'{code}: missing/conflicting CDS protein_id to GeneID mapping; provide an explicit TSV for other GFF schemas')
                mapping_rows=[dict(protein_id=h.split()[0],gene_id='GeneID_'+next(iter(gene_sets[h.split()[0]]))) for h,s in rows]
                provenance.update(gff=str(_layout_resource(gffs[0]).resolve()),gff_sha256=digest(gffs[0]),annotation_metadata=metadata,method='ncbi_gff')
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
            mapping=_layout_resource(folder/'protein_gene.tsv')
            write_tsv(mapping,['protein_id','gene_id'],mapping_rows)
        prepared=prepare(code,source,mapping,_layout_resource(folder/'prepared'),internal_stop_policy='mask' if is_benchmark else 'reject')
        preparation=_layout_resource(folder/'prepared/preparation.json')
        provenance.update(preparation=str(_layout_resource(preparation).resolve()),preparation_sha256=digest(preparation),
                          sequence_qc=json.loads(_layout_resource(preparation).read_text())['sequence_qc'])
        if is_benchmark and provenance['sequence_qc']['internal_stop_proteins']:
            qc=provenance['sequence_qc']
            print(f"{code}: masked {qc['internal_stop_sites']} stop markers as X in {qc['internal_stop_proteins']} proteins; see {_layout_resource(folder/'prepared/sequence_qc.tsv')}",flush=True)
        (_layout_resource(folder/'source.json')).write_text(json.dumps(provenance,indent=2)+'\n')
        entries.append(dict(species=code,fasta=str(_layout_resource(prepared).resolve()),gene_level=not is_benchmark,
                            unit_level='benchmark_protein' if is_benchmark else 'gene',source=json.dumps(provenance)))
    panel=_layout_resource(out/'panel.auto.json')
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
        print(build(_layout_resource(a.raw).resolve(),_layout_resource(a.handoff).resolve(),_layout_resource(a.out).resolve(),a.benchmark_refogs,a.benchmark_input_dir))
    except (OSError,ValueError,KeyError) as exc:
        p.exit(1,f'ERROR: {exc}\n')
