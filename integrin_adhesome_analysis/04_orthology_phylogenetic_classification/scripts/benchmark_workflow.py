"""Workflow adapter for the supplied Emms/Kelly Open_Orthobench benchmark.

Ground truth is read only from supplied RefOGs, never inferred from predictions,
crosswalk labels or family rules. Original reference protein IDs are restored from
preparation maps before calling the supplied pairwise scoring implementation.
"""
import argparse
from collections import defaultdict
import contextlib
import csv
import io
import json
from pathlib import Path
import re
import sys

from pipeline import digest, records, write_tsv
from prepare_reference import normalize_sequence, BENCHMARK_SEQUENCE_POLICY
from benchmark_configuration import quality_policy, acceptance, prediction_source

VERSION = '1.0'
FASTA_PATTERN = re.compile(r'\.(?:fa|faa|fasta)(?:\.gz)?$', re.I)
SOURCE = 'https://github.com/davidemms/Open_Orthobench'
FAMILY_RULES = Path(__file__).resolve().parents[2]/'05_conservative_adhesome_family_assignment_audit/family_rules_expanded.json'
FIBRONECTIN_FAMILIES = {'fibronectin_like', 'fibronectin_III_ligand', 'fibronectin_like_candidate'}
SELECTION_NAME = 'adhesome_refog_selection.tsv'


def adhesome_selection(path, refogs):
    """Select independently curated whole RefOGs; never use predicted assignments."""
    families = set(json.loads(FAMILY_RULES.read_text())['families']) | FIBRONECTIN_FAMILIES
    if not path.is_file():
        raise ValueError('Adhesome scope requires a curated selection TSV: '+str(path)+
                         '. Columns: refog, family, curation_source, curation_status (curated).')
    with path.open(encoding='utf-8-sig', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        required = {'refog','family','curation_source','curation_status'}
        if not required <= set(reader.fieldnames or []):
            raise ValueError('Adhesome selection TSV requires: '+', '.join(sorted(required)))
        rows = [{key:(value or '').strip() for key,value in row.items() if key is not None}
                for row in reader if any(value for value in row.values())]
    selected = defaultdict(set)
    seen = set()
    for row in rows:
        name, family = row['refog'], row['family']
        if name not in refogs:
            raise ValueError('Selected RefOG is absent from curated ground truth: '+name)
        if family not in families:
            raise ValueError('Unknown adhesome family in selection: '+family+'; use a canonical expanded-rule family ID')
        if row['curation_status'] != 'curated' or not row['curation_source']:
            raise ValueError('Each selected RefOG/family link needs curated status and independent curation_source: '+name)
        if (name,family) in seen:
            raise ValueError('Duplicate RefOG/family selection: '+name+' / '+family)
        seen.add((name,family));selected[family].add(name)
    if not selected:
        raise ValueError('Adhesome selection contains no curated RefOG/family links; no accuracy score can be assigned')
    return dict(selected), rows, families


def family_reports(selected, catalogue, refogs, uncertain, predictions):
    """Keep separate RefOG truth partitions even when they share a family label."""
    scores, coverage = {}, []
    for family in sorted(catalogue):
        names = sorted(selected.get(family, set()))
        evaluable = [name for name in names if len(refogs[name]-uncertain[name]) >= 2]
        status = 'evaluated' if evaluable else 'too_few_certain_proteins' if names else 'no_curated_ground_truth'
        row = dict(family=family,status=status,selected_refogs=';'.join(names),
                   selected_refog_count=len(names),evaluated_refogs=len(evaluable),
                   certain_reference_proteins=len(set().union(*(refogs[n]-uncertain[n] for n in names))) if names else 0)
        if evaluable:
            metrics, _, _ = measure({n:refogs[n] for n in names}, uncertain, predictions)
            scores[family] = metrics
            row.update(precision=metrics['precision'],ppv=metrics['precision'],recall=metrics['recall'],f1=metrics['f1'])
        coverage.append(row)
    return scores, coverage


def bundle(raw, explicit_refogs=None, explicit_input=None):
    raw = Path(raw)
    if explicit_refogs:
        refogs = Path(explicit_refogs)
    else:
        choices = sorted({p for p in raw.glob('**/RefOGs') if p.is_dir()}) if raw.exists() else []
        if len(choices) > 1:
            raise ValueError('Multiple RefOGs directories; supply --benchmark-refogs explicitly')
        if not choices:
            return None, None
        refogs = choices[0]
    input_dir = Path(explicit_input) if explicit_input else refogs.parent/'Input'
    if not input_dir.is_dir():
        input_dir = refogs.parent
    return refogs, input_dir


def read_groups(folder):
    files = sorted(Path(folder).glob('RefOG*.txt')) or sorted(Path(folder).glob('*.txt'))
    if not files:
        raise ValueError('No curated reference orthogroup text files in '+str(folder))
    groups, uncertain, paths = {}, {}, []
    for path in files:
        genes = [line.strip() for line in path.read_text(encoding='utf-8-sig').splitlines()
                 if line.strip() and not line.lstrip().startswith('#')]
        if not genes or any(re.search(r'\s', gene) for gene in genes) or len(genes) != len(set(genes)):
            raise ValueError('RefOG files require unique, one-per-line protein identifiers: '+str(path))
        groups[path.stem] = set(genes)
        paths.append(path)
        uncertain_path = Path(folder)/'low_certainty_assignments'/path.name
        uncertain[path.stem] = set()
        if uncertain_path.exists():
            uncertain[path.stem] = {line.strip() for line in uncertain_path.read_text().splitlines()
                                    if line.strip() and not line.lstrip().startswith('#')}
            paths.append(uncertain_path)
    memberships = defaultdict(set)
    for name, genes in groups.items():
        for gene in genes-uncertain[name]:
            memberships[gene].add(name)
    duplicated = [gene for gene,names in memberships.items() if len(names)>1]
    if duplicated:
        raise ValueError('Certain reference proteins occur in multiple RefOGs: '+', '.join(duplicated[:10]))
    return groups, uncertain, paths


def input_proteomes(input_dir, refogs):
    """Select supplied benchmark proteomes by exact curated IDs, not gene symbols."""
    expected = set().union(*refogs.values())
    found = []
    for path in sorted(Path(input_dir).iterdir()):
        if not path.is_file() or not FASTA_PATTERN.search(path.name):
            continue
        seqs = {}
        for header, seq in records(path):
            identifier = header.split()[0]
            if identifier in seqs:
                raise ValueError('Duplicate benchmark input ID: '+identifier)
            seqs[identifier] = seq
        if set(seqs) & expected:
            found.append((path, {identifier:normalize_sequence(seq,identifier,internal_stop_policy='mask')
                                 for identifier,seq in seqs.items()}))
    if not found:
        raise ValueError('No benchmark proteome contains exact RefOG identifiers')
    return found


def predicted_groups(results, table=None):
    results = Path(results)
    options = [results/'Phylogenetic_Hierarchical_Orthogroups/N0.tsv', results/'Orthogroups/Orthogroups.tsv']
    path = Path(table) if table else next((p for p in options if p.is_file()), None)
    if path is None:
        raise ValueError('Benchmark needs the full root N0.tsv or Orthogroups.tsv; candidate-filtered groups cannot be evaluated')
    groups = []
    with path.open(encoding='utf-8-sig',newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        metadata = {'HOG','OG','Orthogroup','Gene Tree Parent Clade','Gene_Tree_Parent_Clade'}
        columns = [name for name in reader.fieldnames or [] if name not in metadata]
        if not columns:
            raise ValueError('No species columns in full orthogroup table')
        for row in reader:
            members = [gene.strip() for name in columns for gene in (row.get(name) or '').split(',') if gene.strip()]
            if len(members)!=len(set(members)):
                raise ValueError('Repeated protein inside a predicted orthogroup')
            if members:
                groups.append(set(members))
    seen = set()
    for group in groups:
        if seen & group:
            raise ValueError('Proteins occur in multiple predicted root orthogroups')
        seen.update(group)
    return groups, path


def reference_mapping(panel, benchmark_inputs):
    """Match original IDs and exact sequences under the declared stop-marker policy."""
    panel = Path(panel)
    snapshot = json.loads((panel/'panel.snapshot.json').read_text())
    targets = set(json.loads((panel/'targets.json').read_text()))
    universe, owners, prepared, input_paths, mapped_species = {}, {}, {}, [], set()
    for path, seqs in benchmark_inputs:
        input_paths.append(path)
        for identifier, sequence in seqs.items():
            if identifier in universe:
                raise ValueError('Benchmark protein IDs must be globally unique: '+identifier)
            universe[identifier] = sequence
            owners[identifier] = path.name
    for entry in snapshot['proteomes']:
        if entry['species'] in targets:
            continue
        source = (panel/entry['fasta']).resolve()
        mapping = source.parent/'identifier_map.tsv'
        for name in ('preparation.json','sequence_qc.tsv'):
            if (source.parent/name).is_file():
                input_paths.append(source.parent/name)
        seqs = {h.split()[0]:s.upper().removesuffix('*') for h,s in records(source)}
        rows = []
        if mapping.is_file():
            input_paths.append(mapping)
            with mapping.open(encoding='utf-8-sig',newline='') as handle:
                rows = list(csv.DictReader(handle,delimiter='\t'))
            for row in rows:
                if str(row.get('selected','1'))!='1':
                    continue
                identifier, pid = row['protein_id'], row['prepared_id']
                if identifier in universe and seqs.get(pid)==universe[identifier]:
                    if identifier in prepared and prepared[identifier]!=pid:
                        raise ValueError('Same benchmark input is represented more than once in panel: '+identifier)
                    prepared[identifier] = pid
                    mapped_species.add(entry['species'])
        else:
            # Explicitly gene-level benchmark IDs can already be prepared.
            for pid, sequence in seqs.items():
                candidates = {pid, pid.split('__',1)[-1]} & universe.keys()
                for identifier in candidates:
                    if sequence==universe[identifier]:
                        if identifier in prepared and prepared[identifier]!=pid:
                            raise ValueError('Ambiguous benchmark protein mapping: '+identifier)
                        prepared[identifier] = pid
                        mapped_species.add(entry['species'])
        input_paths.append(source)
    missing = set(universe)-prepared.keys()
    return universe, prepared, owners, missing, input_paths, mapped_species


def measure(refogs, uncertain, predictions):
    """Use the supplied Open_Orthobench equations, with explicit zero-case handling."""
    details, refs, uncertainties = [], [], []
    tp = fp = fn = 0.0
    for name, ref in refogs.items():
        excluded = uncertain.get(name,set())
        ref = ref-excluded
        if len(ref)<2:
            details.append(dict(refog=name,reference_proteins=len(ref),status='too_few_certain_proteins'))
            continue
        group_tp=group_fp=group_fn=0.0
        absent=set(ref)
        splits=0;foreign=set();merged_groups=0
        for predicted in predictions:
            predicted=predicted-excluded
            overlap=len(ref & predicted)
            if not overlap:
                continue
            splits+=1
            extra=predicted-ref
            foreign|=extra
            merged_groups+=bool(extra)
            absent-=predicted
            group_tp+=overlap*(overlap-1)/2
            group_fp+=overlap*(len(predicted)-overlap)
            group_fn+=(len(ref)-overlap)*overlap
        group_fn=(group_fn+len(absent)*(len(ref)-1))/2
        if group_tp+group_fn!=len(ref)*(len(ref)-1)/2:
            raise ValueError('Pairwise reference accounting failed: '+name)
        normalizer=len(ref)-1
        tp+=group_tp/normalizer;fp+=group_fp/normalizer;fn+=group_fn/normalizer
        refs.append(refogs[name]);uncertainties.append(excluded)
        group_precision=group_tp/(group_tp+group_fp) if group_tp+group_fp else 0.0
        group_recall=group_tp/(group_tp+group_fn) if group_tp+group_fn else 0.0
        group_f1=2*group_precision*group_recall/(group_precision+group_recall) if group_precision+group_recall else 0.0
        details.append(dict(refog=name,reference_proteins=len(ref),status='evaluated',
            true_positive_pairs=group_tp,false_positive_pairs=group_fp,false_negative_pairs=group_fn,
            fragments=splits,missing_proteins=len(absent),exactly_correct=group_fp==0 and group_fn==0,
            merged_predicted_groups=merged_groups,foreign_proteins=len(foreign),
            error_category='exact' if group_fp==0 and group_fn==0 else
                'merge_and_loss_or_split' if group_fp and group_fn else 'merge' if group_fp else
                'split_and_missing' if splits>1 and absent else 'split' if splits>1 else 'missing',
            precision=group_precision,ppv=group_precision,recall=group_recall,f1=group_f1))
    if not refs:
        raise ValueError('No RefOG has two certain proteins for pairwise evaluation')
    precision=tp/(tp+fp) if tp+fp else 0.0
    recall=tp/(tp+fn) if tp+fn else 0.0
    f1=2*precision*recall/(precision+recall) if precision+recall else 0.0
    # Exercise the user-supplied implementation as the authoritative scorer.
    from benchmark import calculate_benchmarks_pairwise
    capture=io.StringIO()
    try:
        with contextlib.redirect_stdout(capture):
            upstream=calculate_benchmarks_pairwise(refs,uncertainties,predictions)
        if any(abs(expected-observed/100)>1e-10 for expected,observed in zip((f1,precision,recall),upstream)):
            raise ValueError('Adapter metrics disagree with the supplied Open_Orthobench scorer')
    except ZeroDivisionError:
        if tp:
            raise
        capture.write('Zero true-positive pairs: precision, recall and F-score are reported as zero.\n')
    evaluated=[row for row in details if row['status']=='evaluated']
    exact=sum(row['exactly_correct'] for row in evaluated)
    return dict(f1=f1,precision=precision,ppv=precision,recall=recall,
        weighted_true_positive_pairs=tp,weighted_false_positive_pairs=fp,weighted_false_negative_pairs=fn,
        evaluated_refogs=len(refs),exact_reconstructed_refogs=exact,exact_reconstruction_fraction=exact/len(refs),
        merge_affected_refogs=sum(row['false_positive_pairs']>0 for row in evaluated),
        fragmented_refogs=sum(row['fragments']>1 for row in evaluated),
        macro_refog_precision=sum(row['precision'] for row in evaluated)/len(evaluated),
        macro_refog_recall=sum(row['recall'] for row in evaluated)/len(evaluated),
        macro_refog_f1=sum(row['f1'] for row in evaluated)/len(evaluated)), details, capture.getvalue()


def run(args):
    out = Path(args.out)
    out.mkdir(parents=True,exist_ok=False)
    summary=dict(pipeline_version=VERSION,source=SOURCE,status='unavailable',validated_weight=None,
        protocol=getattr(args,'protocol','open_orthobench'),
        scope=getattr(args,'scope','all'),
        configuration=dict(reference_raw=str(Path(args.reference_raw).resolve()),
            refogs=str(Path(args.refogs).resolve()) if getattr(args,'refogs',None) else None,
            input_dir=str(Path(args.input_dir).resolve()) if getattr(args,'input_dir',None) else None,
            scope=getattr(args,'scope','all'),
            selection=str(Path(args.selection).resolve()) if getattr(args,'selection',None) else None,
            level=getattr(args,'level','root'),
            quality_policy=str(Path(args.quality_policy).resolve()) if getattr(args,'quality_policy',None) else None),
        weight_policy='Reference-cohort pairwise F1 scales orthogroup-context evidence only; not a candidate identity probability',
        sequence_policy=BENCHMARK_SEQUENCE_POLICY,
        inputs={},diagnostics=[])
    entrypoint=Path(__file__).with_name('adhesome_benchmark_subset.py' if summary['scope']=='adhesome' else 'benchmark.py')
    summary['benchmark_entrypoint']=entrypoint.name
    try:
        policy_path,criteria=quality_policy(getattr(args,'quality_policy',None))
        summary['quality_acceptance']=dict(status='not_evaluated',criteria=criteria)
        if summary['scope'] not in ('all','adhesome'):
            raise ValueError('Unknown benchmark scope: '+str(summary['scope']))
        if getattr(args,'selection',None) and summary['scope'] != 'adhesome':
            raise ValueError('--selection requires --scope adhesome')
        refogs_path, inputs_path=bundle(args.reference_raw,getattr(args,'refogs',None),getattr(args,'input_dir',None))
        if refogs_path is None:
            summary['diagnostics'].append('No RefOGs directory found under reference_raw; proteomes alone cannot validate orthogroup accuracy')
        else:
            refogs,uncertain,truth_paths=read_groups(refogs_path)
            selected, selection_rows, catalogue = {}, [], set()
            if summary['scope']=='adhesome':
                selection_path=Path(args.selection) if getattr(args,'selection',None) else Path(args.reference_raw)/SELECTION_NAME
                selected,selection_rows,catalogue=adhesome_selection(selection_path,refogs)
                truth_paths += [selection_path,FAMILY_RULES]
                summary.update(selection_file=str(selection_path.resolve()),
                    selection_sha256=digest(selection_path),selected_families=sorted(selected),
                    weight_policy='Per-family curated RefOG F1 scales orthogroup-context evidence only for that exact family; untested families remain unvalidated; not a candidate identity probability')
                write_tsv(out/'selection.tsv',['refog','family','curation_source','curation_status','refog_sha256'],[
                    dict(**{key:row[key] for key in ('refog','family','curation_source','curation_status')},
                         refog_sha256=digest(refogs_path/(row['refog']+'.txt'))) for row in selection_rows])
            benchmark_inputs=input_proteomes(inputs_path,refogs)
            stop_proteins=stop_sites=terminal_stops=0
            affected_truth=set()
            truth_ids=set().union(*refogs.values())
            for path,seqs in benchmark_inputs:
                for header,sequence in records(path):
                    identifier=header.split()[0]
                    sites=sequence.removesuffix('*').count('*')
                    terminal_stops+=int(sequence.endswith('*'))
                    stop_proteins+=int(sites>0)
                    stop_sites+=sites
                    if sites and identifier in truth_ids:affected_truth.add(identifier)
            summary['sequence_qc']=dict(internal_stop_proteins=stop_proteins,internal_stop_sites=stop_sites,
                terminal_stops_removed=terminal_stops,affected_refog_proteins=sorted(affected_truth))
            if stop_proteins:
                summary['diagnostics'].append(f'{stop_sites} internal stop markers masked as X in {stop_proteins} benchmark proteins, including {len(affected_truth)} RefOG members; metrics describe this explicitly processed cohort, not unmodified-sequence performance')
            universe,mapping,owners,missing,mapping_paths,cohort_species=reference_mapping(args.panel,benchmark_inputs)
            predicted_path,prediction_metadata,prediction_inputs=prediction_source(args.results,getattr(args,'level','root'),cohort_species)
            predictions,predicted_path=predicted_groups(args.results,table=predicted_path)
            summary['prediction_configuration']=prediction_metadata
            panel_species={row['species'] for row in json.loads((Path(args.panel)/'panel.snapshot.json').read_text())['proteomes']}
            summary['inference_context']=dict(benchmark_species=len(cohort_species),panel_species=len(panel_species),
                additional_species=sorted(panel_species-cohort_species),
                scope='benchmark_only' if panel_species==cohort_species else 'embedded_reference_cohort')
            if summary['protocol']=='open_orthobench':
                expected=dict(reference_proteomes=12,reference_orthogroups=70,reference_members=1945,input_proteins=251378)
                observed=dict(reference_proteomes=len(benchmark_inputs),reference_orthogroups=len(refogs),
                              reference_members=sum(len(group) for group in refogs.values()),input_proteins=len(universe))
                summary['original_dataset_integrity']=dict(expected=expected,observed=observed)
                if observed!=expected:
                    raise ValueError('Original Open_Orthobench dataset is incomplete or replaced: '+str(observed)+'. Supply the complete original Input and RefOGs bundle; custom ground truth requires explicit --protocol custom.')
            inputs=[*truth_paths,*mapping_paths,predicted_path,Path(args.panel)/'workflow_manifest.json',
                    Path(args.panel)/'panel.snapshot.json',Path(__file__),Path(__file__).with_name('benchmark.py'),
                    Path(__file__).with_name('prepare_reference.py'),entrypoint,
                    Path(__file__).with_name('benchmark_configuration.py'),policy_path,*prediction_inputs]
            summary['inputs']={str(p.resolve()):digest(p) for p in inputs}
            summary.update(refogs_directory=str(refogs_path.resolve()),benchmark_input_directory=str(inputs_path.resolve()),
                panel_manifest_sha256=digest(Path(args.panel)/'workflow_manifest.json'),prediction_file=str(predicted_path.resolve()),
                prediction_sha256=digest(predicted_path),reference_proteomes=len(benchmark_inputs),
                expected_reference_proteins=len(universe),mapped_reference_proteins=len(mapping),
                reference_orthogroups=len(refogs),certain_reference_proteins=len(set().union(*(refogs[n]-uncertain[n] for n in refogs))))
            selected_names=set().union(*selected.values()) if selected else set(refogs)
            evaluation_refogs={name:group for name,group in refogs.items() if name in selected_names}
            summary.update(selected_reference_orthogroups=len(evaluation_refogs),
                selected_certain_reference_proteins=len(set().union(*(refogs[n]-uncertain[n] for n in evaluation_refogs))),
                evaluated_scope='adhesome_subset_reference_orthogroup_inference' if summary['scope']=='adhesome' else 'reference_cohort_orthogroup_inference',
                prediction_scope='Complete groups at the explicitly selected tree level projected to benchmark proteomes; no adhesome-only member filtering')
            write_tsv(out/'identifier_mapping.tsv',['benchmark_id','prepared_id','source_proteome','status'],[
                dict(benchmark_id=identifier,prepared_id=mapping.get(identifier,''),source_proteome=owners[identifier],
                     status='matched_id_and_sequence' if identifier in mapping else 'missing_or_sequence_mismatch') for identifier in sorted(universe)])
            missing_truth=set().union(*(refogs[n]-uncertain[n] for n in refogs))-universe.keys()
            if missing or missing_truth:
                summary['status']='incompatible_reference_cohort'
                summary['diagnostics'].append(f'{len(missing)} benchmark input proteins missing/changed in prepared panel; {len(missing_truth)} certain RefOG proteins missing from benchmark inputs')
            else:
                reverse={pid:identifier for identifier,pid in mapping.items()}
                projected=[{reverse[pid] for pid in group if pid in reverse} for group in predictions]
                projected=[group for group in projected if group]
                clustered=set().union(*projected) if projected else set()
                summary['clustered_reference_fraction']=len(clustered)/len(universe)
                # Only truth groups are selected. Retain ALL benchmark-universe proteins
                # in each prediction, including non-adhesome and unselected competitors.
                if summary['scope']=='adhesome':
                    scores,coverage=family_reports(selected,catalogue,refogs,uncertain,projected)
                    summary['family_metrics']=scores
                    summary['family_coverage']=dict(catalogue_families=len(catalogue),
                        selected_families=len(selected),evaluated_families=len(scores),
                        without_curated_ground_truth=sum(row['status']=='no_curated_ground_truth' for row in coverage))
                    write_tsv(out/'family_coverage.tsv',['family','status','selected_refogs','selected_refog_count',
                        'evaluated_refogs','certain_reference_proteins','precision','ppv','recall','f1'],coverage)
                metrics,details,log=measure(evaluation_refogs,uncertain,projected)
                summary['metrics']=metrics
                summary['quality_acceptance']=acceptance(metrics,criteria)
                write_tsv(out/'refog_scores.tsv',list(dict.fromkeys(k for row in details for k in row)),details)
                write_tsv(out/'orthogroup_error_matrix.tsv',['error_category','refog_count'],[
                    dict(error_category=category,refog_count=sum(row.get('error_category')==category for row in details))
                    for category in sorted({row['error_category'] for row in details if 'error_category' in row})])
                write_tsv(out/'acceptance_checks.tsv',['metric','observed','minimum','passed'],[
                    dict(metric=metric,**check) for metric,check in summary['quality_acceptance']['checks'].items()])
                (out/'upstream_scores.log').write_text(log,encoding='utf-8')
                if summary['clustered_reference_fraction']<0.5:
                    summary['status']='insufficient_prediction_coverage'
                    summary['diagnostics'].append('Fewer than half of benchmark input proteins are represented in root orthogroups; no audit weight applied')
                else:
                    if summary['quality_acceptance']['status']=='accepted':
                        summary.update(status='validated',validated_weight=metrics['f1'])
                    else:
                        summary.update(status='quality_rejected',validated_weight=None)
                        summary['diagnostics'].append('Valid benchmark inputs, but performance fails the configured acceptance criteria; measured scores are retained and no validated audit weight is applied')
    except (ValueError,OSError,KeyError,AssertionError) as exc:
        summary.update(status='invalid_benchmark_inputs',validated_weight=None)
        summary['diagnostics'].append(str(exc))
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
    (out/'README.md').write_text('# Evolutionary benchmark\n\nStatus: '+summary['status']+'\n\n'+
        '\n'.join(summary['diagnostics'])+'\n\n'+summary['weight_policy']+'\n',encoding='utf-8')
    return summary


def main(argv=None, *, subset_only=False):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--version',action='version',version='adhesome-orthobench '+VERSION)
    parser.add_argument('--results',type=Path,required=True)
    parser.add_argument('--panel',type=Path,required=True)
    parser.add_argument('--reference-raw',type=Path,required=True)
    parser.add_argument('--refogs',type=Path)
    parser.add_argument('--input-dir',type=Path)
    parser.add_argument('--protocol',choices=['open_orthobench','custom'],default='open_orthobench',help='Original dataset integrity checks are enabled by default; custom curated cohorts require an explicit opt-in')
    parser.add_argument('--scope',choices=['adhesome'] if subset_only else ['all','adhesome'],
        default='adhesome' if subset_only else 'all',help='Select truth RefOGs only; predicted benchmark group members are never filtered to adhesome proteins')
    parser.add_argument('--selection',type=Path,help='Curated RefOG/family TSV for adhesome scope; default: reference_raw/adhesome_refog_selection.tsv')
    parser.add_argument('--level',choices=['root','cohort-lca'],default='root',help='Root preserves the original benchmark; cohort-lca selects the smallest labelled species-tree node containing all benchmark species, independently of scores')
    parser.add_argument('--quality-policy',type=Path,help='Acceptance policy JSON; bundled default requires pairwise precision >=0.80 and F1 >=0.70')
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args(argv)
    summary=run(args)
    print('Benchmark report: '+str(args.out/'summary.json')+' ('+summary['status']+')')
    return 0


if __name__=='__main__':
    sys.exit(main())
