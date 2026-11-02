"""Read-only reuse audit and previous orthogroup benchmark investigation."""

# Resolve source resources independently of the caller working directory.
import sys as _layout_sys
from pathlib import Path as _LayoutPath
_LAYOUT_BUNDLE_ROOT = next(p for p in _LayoutPath(__file__).resolve().parents
                          if (p/'config/layout_paths.json').is_file() and (p/'pipeline_paths.py').is_file())
if str(_LAYOUT_BUNDLE_ROOT) not in _layout_sys.path:
    _layout_sys.path.insert(0, str(_LAYOUT_BUNDLE_ROOT))
from pipeline_paths import activate as _layout_activate, resource_path as _layout_resource
_layout_activate()

from collections import Counter
import json
from pathlib import Path
import statistics

from .datasets import fasta, sha, dump, table
from .reports import read


def length_summary(values):
    values=sorted(values)
    return dict(n=len(values),minimum=min(values) if values else None,
        median=statistics.median(values) if values else None,maximum=max(values) if values else None,
        below_100=sum(x<100 for x in values), below_300=sum(x<300 for x in values))


def audit_existing(root,out,corrected=None):
    """No old cluster is silently adopted into a new analytical context."""
    root=_layout_resource(Path(root)).resolve();out=_layout_resource(Path(out)).resolve();out.mkdir(parents=True,exist_ok=True)
    inventory=[]
    rules=[('**/prepared/*.id_map.csv','reusable_after_sequence_and_namespace_verification','Identifier mappings'),
           ('**/domains/scan.json','reusable_after_filtering_or_rescoring','Domain scans: verify sequence/library/threshold hashes'),
           ('**/motifs/*','reusable_after_filtering_or_rescoring','Verify sequence, regions and motif-library hashes'),
           ('**/joined_evidence.tsv','partial_recomputation','Independent features reusable; replace all orthogroup-dependent fields'),
           ('**/inference_complete.json','invalidated_for_isolated_inference','Old mixed-query inference cannot be fixed by filtering'),
           ('**/Stage_2_output.xlsx','partial_recomputation','Re-score after isolated orthology; retain as historical result'),
           ('**/benchmark/summary.json','directly_reusable_as_historical_comparison','Measured historical metrics; not a new inference')]
    for pattern,status,reason in rules:
        for p in sorted(_layout_resource(root).glob(pattern)):
            if _layout_resource(p).is_file():inventory.append(dict(path=str(p),sha256=sha(p),reuse_class=status,condition=reason))
    table(_layout_resource(out/'reuse_inventory.tsv'),inventory)
    lengths=[];fragment_headers=[];invalid=[]
    for p in sorted(_layout_resource(root).glob('handoff/detection/*/prepared/*.faa')):
        try:seqs=fasta(p)
        except ValueError as exc:
            invalid.append(dict(path=str(p),error=str(exc),interpretation='Historical failed input; not adopted into corrected inference'))
            continue
        lengths.append(dict(species=p.stem,**length_summary([len(s) for h,s in seqs.values()])))
        fragment_headers.extend(dict(file=str(p),identifier=i,reason='explicit_fragment_header') for i,(h,s) in seqs.items() if 'fragment' in h.lower())
    table(_layout_resource(out/'sequence_length_distributions.tsv'),lengths)
    table(_layout_resource(out/'explicit_fragment_annotations.tsv'),fragment_headers)
    table(_layout_resource(out/'invalid_historical_inputs.tsv'),invalid)
    benchmark=_layout_resource(root/'annotation/downstream/benchmark/summary.json')
    old=json.loads(_layout_resource(benchmark).read_text()) if _layout_resource(benchmark).exists() else {}
    group_paths=list(_layout_resource(root).glob('annotation/downstream/orthology/inference/Results_*/Orthogroups/Orthogroups.tsv'))
    groups=[];species=set()
    if len(group_paths)==1:
        for r in read(group_paths[0]):
            members={k:[v.strip() for v in str(value or '').split(',') if v.strip()] for k,value in r.items() if k!='Orthogroup'}
            species.update(members)
            control=sum(len(v) for k,v in members.items() if k.lower() in ('validation','calibration') or k.startswith(('Cal_','Val_')))
            if control:
                groups.append(dict(orthogroup=r.get('Orthogroup',''),evaluation_control_members=control,
                    total_members=sum(map(len,members.values())),species_with_multiple_copies=sum(len(v)>1 for v in members.values()),
                    interpretation='Multiple copies flag paralogy review; do not establish a mistaken merge alone'))
    table(_layout_resource(out/'control_containing_orthogroups.tsv'),groups)
    result=dict(status='read_only_audit_complete',historical_benchmark=str(benchmark),
        historical_benchmark_sha256=sha(benchmark) if _layout_resource(benchmark).exists() else None,
        historical_metrics=old.get('metrics'),old_inferred_species=sorted(species),
        mixed_evaluation_controls_detected=bool(groups),control_containing_groups=len(groups),
        short_sequence_causality='not_established',
        domain_only_similarity='Requires isolated inference and per-member similarity/domain-coverage inspection; length alone is insufficient.',
        reference_mapping='Existing benchmark sequence hashes and mapping diagnostics must be revalidated in each isolated context.',
        membership_comparison='pending_corrected_inference',
        no_valid_short_proteins_removed=True,artifacts_audited=len(inventory))
    if corrected:
        new=_layout_resource(Path(corrected)/'benchmark/summary.json')
        if not _layout_resource(new).exists():raise ValueError('Corrected inference has no benchmark summary')
        if not (_layout_resource(Path(corrected).parent.parent/'run_contract.json')).exists():
            raise ValueError('Corrected comparison requires an isolated-run contract')
        result.update(corrected_metrics=json.loads(_layout_resource(new).read_text()).get('metrics'),
            corrected_benchmark_sha256=sha(new),membership_comparison='Separate inference benchmark available; raw group IDs are not stable identifiers across runs')
    dump(_layout_resource(out/'diagnosis.json'),result)
    (_layout_resource(out/'README.md')).write_text('# Existing-results audit\n\nEvaluation controls were present in the old inference. '
        'This is evidence of dataset mixing, not proof that controls or short proteins caused the low F1. '
        'The causal comparison requires new isolated inference. Filtering old clusters is insufficient. '
        'No analysis outputs were changed. See reuse_inventory.tsv and diagnosis.json.\n',encoding='utf-8')
    return result
