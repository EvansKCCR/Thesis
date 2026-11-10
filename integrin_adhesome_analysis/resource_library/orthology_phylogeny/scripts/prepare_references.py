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
from prepare_reference import prepare, DEFAULT_QC
from benchmark_workflow import bundle, read_groups, input_proteomes
from input_dependencies import InputRegistry


def build(raw, handoff, out, benchmark_refogs=None, benchmark_input_dir=None, reference_files=None, reference_panel=None):
    registry=InputRegistry(_LAYOUT_BUNDLE_ROOT)
    raw=registry.require('reference_raw','orthology.reference_preparation',raw)
    meta=json.loads((_layout_resource(handoff/'workflow_manifest.json')).read_text())
    targets=meta['target_species']
    selected=registry.select_references(panel=reference_panel,files=reference_files,raw=raw)
    files=selected if selected is not None else sorted(p for p in raw.iterdir() if _layout_resource(p).is_file() and
        re.search(r'\.(fa|faa|fasta)(\.gz)?$',p.name,re.I))
    if any(not re.search(r'\.(fa|faa|fasta)(\.gz)?$',p.name,re.I) for p in files):
        raise ValueError('Reference preparation supports .fa, .faa, .fasta (optionally .gz)')
    # An explicit panel never expands itself into unrelated benchmark proteomes.
    refogs,input_dir=bundle(raw,benchmark_refogs,benchmark_input_dir) if selected is None else (None,None)
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
            orthology_role='target',
            unit_level='control_protein' if control_unit else 'gene',
            source='Handoff-verified labelled protein units; not a complete proteome' if control_unit else 'Handoff-verified complete target proteome'))
    seen=set(targets)
    qc_config=json.loads(_layout_resource(_LAYOUT_BUNDLE_ROOT/'config/detection/proteomes.json').read_text(encoding='utf-8')).get('qc',DEFAULT_QC)
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
        # Benchmark FASTAs keep their upstream filename (often a long Ensembl
        # release name), while the prepared species code is prefixed with
        # ``Bench_``. Accept source-stem sidecars as well as species-code files
        # so the isoform filter always sees the real protein-to-gene mapping.
        sidecar_stems = [source.name, source.stem]
        maps=sorted({_layout_resource(source.parent/f'{stem}.protein_gene.tsv') for stem in sidecar_stems
                     if _layout_resource(source.parent/f'{stem}.protein_gene.tsv').is_file()})
        if len(maps)>1:
            raise ValueError(f'Ambiguous protein_gene.tsv sidecars for {code}; specify a single matching map')
        rows=list(records(source))
        provenance=dict(fasta=str(_layout_resource(source).resolve()),fasta_sha256=digest(source))
        metadata_candidates=sorted({_layout_resource(source.parent/f'{stem}.reference.json') for stem in [code,*sidecar_stems]
                                    if _layout_resource(source.parent/f'{stem}.reference.json').is_file()})
        if len(metadata_candidates)>1:
            raise ValueError(f'Ambiguous assembly/provenance metadata sidecars for {code}')
        panel_metadata=metadata_candidates[0] if metadata_candidates else _layout_resource(source.parent/f'{code}.reference.json')
        if _layout_resource(panel_metadata).is_file():
            scope=json.loads(_layout_resource(panel_metadata).read_text(encoding='utf-8'))
            if scope.get('fasta_sha256')!=provenance['fasta_sha256']:
                raise ValueError(f'{code}: reference-scope metadata does not match the FASTA')
            provenance.update(reference_scope=scope,reference_scope_file=str(_layout_resource(panel_metadata).resolve()),
                              reference_scope_sha256=digest(panel_metadata))
        mapping=maps[0] if maps else None
        if mapping is not None:
            provenance.update(mapping=str(_layout_resource(mapping).resolve()),mapping_sha256=digest(mapping),method='explicit_tsv')
        else:
            provenance['method']='header_suffix_and_cdhit_fallback'
        prepared=prepare(code,source,mapping,_layout_resource(folder/'prepared'),internal_stop_policy='mask',qc=qc_config)
        preparation=_layout_resource(folder/'prepared/preparation.json')
        provenance.update(preparation=str(_layout_resource(preparation).resolve()),preparation_sha256=digest(preparation),
                          sequence_qc=json.loads(_layout_resource(preparation).read_text())['sequence_qc'])
        if is_benchmark:
            provenance.update(benchmark_refogs=str(_layout_resource(refogs).resolve()),
                              benchmark_truth_sha256={str(p.relative_to(refogs)):digest(p) for p in truth_paths},
                              unit_level='gene',isoform_policy='longest QC-passing protein per mapped/derived gene; ties by protein ID')
        if is_benchmark and provenance['sequence_qc']['internal_stop_proteins']:
            qc=provenance['sequence_qc']
            print(f"{code}: masked {qc['internal_stop_sites']} stop markers as X in {qc['internal_stop_proteins']} proteins; see {_layout_resource(folder/'prepared/sequence_qc.tsv')}",flush=True)
        (_layout_resource(folder/'source.json')).write_text(json.dumps(provenance,indent=2)+'\n')
        entries.append(dict(species=code,fasta=str(_layout_resource(prepared).resolve()),gene_level=True,
                            unit_level='gene',orthology_role='benchmark_core' if is_benchmark else 'reference_outgroup',
                            source=json.dumps(provenance)))
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
    p.add_argument('--reference-file',action='append')
    p.add_argument('--reference-panel')
    a=p.parse_args()
    try:
        print(build(_layout_resource(a.raw).resolve(),_layout_resource(a.handoff).resolve(),_layout_resource(a.out).resolve(),a.benchmark_refogs,a.benchmark_input_dir,a.reference_file,a.reference_panel))
    except (OSError,ValueError,KeyError) as exc:
        p.exit(1,f'ERROR: {exc}\n')
