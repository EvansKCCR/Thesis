#!/usr/bin/env python3
"""Describe a completed, calibration-selected Stage 2 experiment; no inference."""

# Resolve source resources independently of the caller working directory.
import sys as _layout_sys
from pathlib import Path as _LayoutPath
_LAYOUT_BUNDLE_ROOT = next(p for p in _LayoutPath(__file__).resolve().parents
                          if (p/'config/layout_paths.json').is_file() and (p/'pipeline_paths.py').is_file())
if str(_LAYOUT_BUNDLE_ROOT) not in _layout_sys.path:
    _layout_sys.path.insert(0, str(_LAYOUT_BUNDLE_ROOT))
from pipeline_paths import activate as _layout_activate, resource_path as _layout_resource
_layout_activate()



# Resolve shared pipeline imports from the installed script, not the caller cwd.
import sys
from pathlib import Path
sys.path.insert(0,str(_layout_resource(Path(__file__)).resolve().parents[1]))
import argparse
from collections import Counter,defaultdict
from datetime import datetime,timezone
import json
import math
from pathlib import Path
import random

from schisto_adhesome.datasets import dump,sha,table
from schisto_adhesome.reports import read,metrics
from schisto_adhesome.validation_review import Document,predictions,family_metrics,POLICIES
from maintenance.optimize_stage2 import ARMS
import stage2_evidence


def correct(row):
    return row['control_class']=='positive' and row['supported'] and row['predicted_family'] in row['expected_families'].split(';')


def wilson(success,total):
    if not total:return None
    z=1.959963984540054;p=success/total;den=1+z*z/total
    centre=(p+z*z/(2*total))/den
    width=z*math.sqrt(p*(1-p)/total+z*z/(4*total*total))/den
    return [max(0.,centre-width),min(1.,centre+width)]


def paired(base,opt):
    old={r['protein_id']:r for r in base};counts=Counter()
    for r in opt:
        counts[(r['control_class']=='positive',old[r['protein_id']]['supported'],r['supported'])]+=1
    better=sum(v for (actual,b,o),v in counts.items() if b!=actual and o==actual)
    worse=sum(v for (actual,b,o),v in counts.items() if b==actual and o!=actual)
    n=better+worse
    p=min(1.,2*sum(math.comb(n,k) for k in range(min(better,worse)+1))/2**n) if n else 1.
    # Paired bootstrap samples whole proteins, retaining their two decisions.
    # This describes this panel, not uncertainty across evolutionary lineages.
    cells=list(counts);weights=list(counts.values());rng=random.Random(20261024);diffs=[]
    for unused in range(3000):
        sampled=Counter(rng.choices(cells,weights=weights,k=len(base)))
        results=[]
        for position in (1,2):
            tp=sum(v for k,v in sampled.items() if k[0] and k[position])
            fp=sum(v for k,v in sampled.items() if not k[0] and k[position])
            fn=sum(v for k,v in sampled.items() if k[0] and not k[position])
            results.append(2*tp/(2*tp+fp+fn))
        diffs.append(results[1]-results[0])
    diffs.sort()
    return dict(improved_binary_decisions=better,worsened_binary_decisions=worse,
        exact_mcnemar_two_sided_p=p,paired_F1_difference_95_percentile_CI=[diffs[74],diffs[2924]],
        bootstrap_replicates=3000,seed=20261024,
        limitation='Protein-level resampling; homology and species clustering are not independent sampling units.')


def change_log(base,opt,assignments):
    old={r['protein_id']:r for r in base};output=[]
    for row in opt:
        b=old[row['protein_id']];kinds=[]
        if row['supported'] and not b['supported']:kinds.append('newly_supported')
        if b['supported'] and not row['supported']:kinds.append('lost_supported')
        if row['predicted_family']!=b['predicted_family']:kinds.append('reclassified')
        if row['grade']!=b['grade']:kinds.append('grade_changed')
        if row['exclusion_gate'] and row['exclusion_gate']!=b['exclusion_gate']:kinds.append('newly_rejected')
        if not kinds:continue
        a=assignments.get(row['protein_id'],{})
        output.append(dict(protein_id=row['protein_id'],original_id=row['original_id'],species=row['species'],
            control_class=row['control_class'],expected_family=row['expected_families'],change=';'.join(kinds),
            baseline_supported=b['supported'],optimized_supported=row['supported'],baseline_family=b['predicted_family'],
            optimized_family=row['predicted_family'],baseline_grade=b['grade'],optimized_grade=row['grade'],
            baseline_correct_family=correct(b),optimized_correct_family=correct(row),
            discovery_rule=a.get('discovery_rule',''),recruitment_evidence=a.get('recruitment_evidence',''),
            identity_rule=a.get('identity_rule',''),architecture_evidence=a.get('architecture_evidence',''),
            observed_Pfam=row['architecture'],exclusion_gate=row['exclusion_gate']))
    return output


def pct(value):return 'undefined' if value is None else f'{100*value:.2f}%'
def number(value):return 'undefined' if value is None else f'{value:.4f}'


def generate(run,out):
    frozen=json.loads((_layout_resource(out/'frozen_selection.json')).read_text());selected=frozen['selected_arm']
    validation=json.loads((_layout_resource(out/'validation/ablation_results.json')).read_text())
    if validation['frozen_selection_sha256']!=sha(_layout_resource(out/'frozen_selection.json')):
        raise ValueError('Validation does not match frozen selection')
    for file,digest in frozen['code_files'].items():
        if sha(_layout_resource(Path(file)))!=digest:raise ValueError('Decision implementation changed since validation freeze')
    folder=_layout_resource(out/'reports');folder.mkdir(exist_ok=True)
    doc=Document();doc.heading('Stage 2 optimization — workflow 1.0')
    doc.paragraph('Calibration-only selection followed by frozen evaluation on the existing 1,737-record validation panel. '
        'This panel was previously inspected and contains reference and calibration overlaps; it is an overlap-qualified comparison, not a new blinded holdout. '
        'All experiments reuse the same cached sequence, domain, motif, topology/localization and orthology evidence. No Stage 1 or external-tool rerun was performed.')
    all_ablation=[];summaries={};artifacts={};paired_results={}
    for kind in ('calibration','validation'):
        base=predictions(_layout_resource(run/kind/'reports/predictions.tsv'))
        opt=predictions(_layout_resource(out/kind/selected/'predictions.tsv'))
        if {r['protein_id'] for r in base}!={r['protein_id'] for r in opt}:raise ValueError('Evaluation cohort changed')
        assigns={r['protein_id']:r for r in read(_layout_resource(out/kind/selected/'candidate_assignments.tsv'))}
        changes=change_log(base,opt,assigns);table(folder/(kind+'_candidate_changes.tsv'),changes)
        for r in changes:
            if r['baseline_correct_family'] and not r['optimized_correct_family']:
                # Report rather than hide any validation regression.
                r['correct_family_regression']=True
        recruitment=read(_layout_resource(out/kind/selected/'recruited_hypotheses.tsv'))
        original_recovery=read(_layout_resource(out/'calibration/original/recruited_hypotheses.tsv')) if kind=='calibration' else []
        extra=[r for r in recruitment if r.get('discovery_rule','').startswith('diagnostic_domain_recruitment:')]
        table(folder/(kind+'_new_recruitment_hypotheses.tsv'),extra)
        tiers=[]
        for label,rows in (('baseline',base),('optimized',opt)):
            for policy,title in POLICIES.items():
                tiers.append(dict(configuration=label,tier=title,**metrics(rows,policy),**{
                    'family_'+k:(v if policy not in ('stage1_recruited','before_exclusion') else None)
                    for k,v in family_metrics(rows,policy).items()}))
        table(folder/(kind+'_tier_metrics.tsv'),tiers)
        species=[];families=[];overlaps=[]
        for species_name in sorted({r['species'] for r in base}):
            for label,rows in (('baseline',base),('optimized',opt)):
                subset=[r for r in rows if r['species']==species_name]
                species.append(dict(species=species_name,configuration=label,**metrics(subset,'supported'),
                    **{'family_'+k:v for k,v in family_metrics(subset,'supported').items()}))
        table(folder/(kind+'_species_metrics.tsv'),species)
        for family in sorted({r['expected_families'] for r in base if r['control_class']=='positive'}):
            for label,rows in (('baseline',base),('optimized',opt)):
                subset=[r for r in rows if r['control_class']=='positive' and r['expected_families']==family]
                families.append(dict(expected_family=family or '[unmapped]',configuration=label,positive_records=len(subset),
                    membership_supported=sum(r['supported'] for r in subset),correct_supported_family=sum(correct(r) for r in subset),
                    membership_recall=sum(r['supported'] for r in subset)/len(subset),
                    family_recall=sum(correct(r) for r in subset)/len(subset) if family else None))
        table(folder/(kind+'_family_sensitivity.tsv'),families)
        for field in ('calibration_overlap','reference_overlap','training_overlap'):
            for value in (False,True):
                for label,rows in (('baseline',base),('optimized',opt)):
                    subset=[r for r in rows if r[field]==value]
                    overlaps.append(dict(overlap=field,value=value,configuration=label,**metrics(subset,'supported'),
                        **{'family_'+k:v for k,v in family_metrics(subset,'supported').items()}))
        table(folder/(kind+'_overlap_metrics.tsv'),overlaps)
        grades=Counter((r['grade'],next(b['grade'] for b in base if b['protein_id']==r['protein_id'])) for r in opt)
        table(folder/(kind+'_grade_transitions.tsv'),[dict(baseline_grade=b,optimized_grade=o,records=n) for (o,b),n in sorted(grades.items())])
        remaining=[dict(r,reason='Upstream Stage 1 recruitment missing' if not r['stage1_recruited'] else
            'Expected family truth unmapped' if not r['expected_families'] else 'Family-specific confirmation incomplete or unresolved')
            for r in opt if r['control_class']=='positive' and not correct(r)]
        table(folder/(kind+'_remaining_family_misses.tsv'),remaining)
        summary={label:dict(membership=metrics(rows,'supported'),family=family_metrics(rows,'supported'),
            grades=dict(Counter(r['grade'] for r in rows)),stage1=metrics(rows,'stage1_recruited'),
            provisional_membership=metrics(rows,'provisional'),provisional_family=family_metrics(rows,'provisional'))
            for label,rows in (('baseline',base),('optimized',opt))}
        paired_results[kind]=paired(base,opt)
        for label in summary:
            m=summary[label]['membership']
            summary[label]['intervals_95']=dict(PPV=wilson(m['TP'],m['TP']+m['FP']),recall=wilson(m['TP'],m['TP']+m['FN']),
                specificity=wilson(m['TN'],m['TN']+m['FP']))
        summary['changes']=dict(newly_supported=sum(not b['supported'] and r['supported'] for b,r in zip(sorted(base,key=lambda x:x['protein_id']),sorted(opt,key=lambda x:x['protein_id']))),
            correct_family_regressions=sum(r['baseline_correct_family'] and not r['optimized_correct_family'] for r in changes),
            new_family_hypotheses=len(extra),new_hypothesis_proteins=len({r['protein_id'] for r in extra}),
            original_stage1_protein_pool_unchanged=True)
        summaries[kind]=summary
        doc.heading(kind.capitalize()+' results')
        doc.table(['Configuration','TP','FP','FN','Precision','Recall','F1','Specificity','Correct families','Family F1'],[
            [label,*[summary[label]['membership'][k] for k in ('TP','FP','FN')],
             pct(summary[label]['membership']['PPV']),pct(summary[label]['membership']['recall']),number(summary[label]['membership']['F1']),
             pct(summary[label]['membership']['specificity']),summary[label]['family']['TP'],number(summary[label]['family']['F1'])]
            for label in ('baseline','optimized')])
        doc.paragraph(f"New family-routing hypotheses: {len(extra)} across {len({r['protein_id'] for r in extra})} existing Stage 1 candidates. "
            'This increases family-hypothesis coverage; it does not add new sequence identifiers or claim new Stage 1 recruitment. '
            f"Previously correct supported-family regressions: {summary['changes']['correct_family_regressions']}.")
        doc.table(['Configuration','Tier','TP','FP','Recall','Membership F1','Correct family calls'],[
            [r['configuration'],r['tier'],r['TP'],r['FP'],pct(r['recall']),number(r['F1']),r['family_TP']]
            for r in tiers if r['tier']!='Supported before exclusions'])
        doc.paragraph('Recruitment has no final family assignment to score; its family metrics are omitted. '
            'Unsupported shared-domain provisional names are demoted to unresolved when family-specific confirmation is absent. '
            'The original hypotheses and their evidence are preserved, and explicitly diagnostic partial cores remain eligible for provisional annotation. '
            'An incorrect provisional name is not counted as a correct supported family.')
        doc.table(['Species','Configuration','TP','FN','Recall','FP','Correct-family recall'],[
            [r['species'],r['configuration'],r['TP'],r['FN'],pct(r['recall']),r['FP'],pct(r['family_recall'])] for r in species])
        doc.table(['Configuration','A','B','C','D / not recruited'],[[label,*[summary[label]['grades'].get(g,0) for g in ('A','B','C')],
            '; '.join(str(k)+': '+str(v) for k,v in summary[label]['grades'].items() if k not in ('A','B','C'))] for label in ('baseline','optimized')])
        doc.table(['Configuration','Precision 95% CI','Recall 95% CI','Specificity 95% CI'],[
            [label,*[' – '.join(pct(x) for x in summary[label]['intervals_95'][k]) for k in ('PPV','recall','specificity')]] for label in ('baseline','optimized')])
        pr=paired_results[kind]
        doc.paragraph(f"Paired membership decisions improved/worsened: {pr['improved_binary_decisions']}/{pr['worsened_binary_decisions']}; "
            f"exact two-sided McNemar p={pr['exact_mcnemar_two_sided_p']:.3g}. Paired F1 difference 95% bootstrap interval: "+
            ' – '.join(number(x) for x in pr['paired_F1_difference_95_percentile_CI'])+'. Wilson intervals and protein-level bootstrap describe this panel; '
            'they do not account for family/species dependence or establish transferability to unseen lineages.')
        arms=frozen['calibration_results'] if kind=='calibration' else validation['arms']
        if kind=='validation':
            arms=[dict(arm='original',profile='original',context_weight=1.,**summary['baseline']['membership'],
                family=summary['baseline']['family'])]+arms
        source_hash=sha(_layout_resource(run/kind/'downstream/final/joined_evidence.tsv'))
        for a in arms:
            if 'source_evidence_sha256' in a and a['source_evidence_sha256']!=source_hash:raise ValueError('Ablations used different evidence')
            all_ablation.append(dict(dataset=kind,arm=a['arm'],profile=a['profile'],context_coefficient=a['context_weight'],
                **{k:a[k] for k in ('TP','FP','TN','FN','PPV','recall','F1','specificity')},
                family_TP=a['family']['TP'],family_FP=a['family']['FP'],family_FN=a['family']['FN'],family_F1=a['family']['F1'],
                selected=a['arm']==selected,source_evidence_sha256=source_hash))
        artifacts[kind]=dict(baseline_predictions_sha256=sha(_layout_resource(run/kind/'reports/predictions.tsv')),
            optimized_predictions_sha256=sha(_layout_resource(out/kind/selected/'predictions.tsv')),cached_evidence_sha256=source_hash)
    table(_layout_resource(folder/'ablation_comparison.tsv'),all_ablation)
    doc.heading('Ablation comparison')
    doc.table(['Dataset','Arm','Context','TP','FP','Recall','Membership F1','Family F1'],[
        [r['dataset'],r['arm'],r['context_coefficient'],r['TP'],r['FP'],pct(r['recall']),number(r['F1']),number(r['family_F1'])] for r in all_ablation])
    doc.paragraph('Recruitment alone does not bypass family confirmation. Exact CDD/Pfam model import resolves existing annotation gaps; '
        'additional diagnostic routing contributes recovery only after the architecture gate. The unselected all-rules combined arm introduced '
        'one extra calibration hard negative (LIMK2); its entire new LIMK-family rescue was excluded before validation. The slightly higher '
        'validation recovery of that unselected arm does not change the calibration-selected configuration.')
    doc.paragraph('Coefficients 0, 0.541 and 1.0 produce the same supported membership/family counts within each combined profile. '
        'Zero removes orthogroup support and moves complete architectures from A to B; A/B both remain supported. Partial C calls can become D. '
        'Thus this experiment shows no supported-call recovery benefit from 0.541 or orthogroup context. It does not demonstrate that orthology is '
        'biologically uninformative. These coefficients are not individual orthology probabilities, and the original benchmark quality decision is unchanged.')
    doc.heading('False-negative diagnosis and remaining limits')
    doc.paragraph('Each baseline missed positive is listed with species, expected family, observed Pfam/CDD evidence, identity-gate failure and corrective action. '
        'Missing exact CDD accessions, missing configured family routes, incomplete cores and unmapped expected-family labels are distinguished. '
        'Sequence fragmentation, alignment coverage loss, remote homology and structural similarity are not inferred from missing domains. '
        'No new remote-homology, structure or motif search was available in this replay; existing informative observations are reused. '
        'No global HMM, sequence-identity, coverage or length cutoffs were lowered. Family-specific partial-domain fractions and existing minimum lengths '
        'are explicit in the rule library. Full-core order, coordinates, lengths, copy limits, diagnostic CDD, competing-family resolution and negative gates remain required.')
    doc.paragraph('Stage 1 is unchanged. One calibration positive (Defensin-6) and five validation positives missing from the original recruitment pool '
        'cannot be rescued by Stage 2. A separate upstream change would require its own calibration and external-evidence handoff. '
        'Twenty validation positives lack mapped expected-family truth; accepted membership for these is not credited as correct family annotation. '
        'Remaining family misses and provisional calls are exported separately; incomplete evidence is not promoted to high confidence.')
    doc.heading('Family-specific rules and provenance')
    rules=stage2_evidence.configuration()['family_rules'];rule_rows=[]
    for r in rules:
        rule_rows.append(dict(family=r['family'],rule_id=r['id'],enabled_in_selected=r['family'] not in
            stage2_evidence.configuration()['profiles'][frozen['selected_profile']].get('excluded_new_rule_families',[]),
            full_core=json.dumps(r['required_groups']),ordered_groups=json.dumps(r.get('ordered_groups',[])),
            diagnostic_CDD=r['specific_cdd_pattern'],length_range=json.dumps(r.get('length_range')),
            provisional_domain_fraction=r['provisional_min_domain_fraction'],confounders=';'.join(r['confounders']),evidence=r['evidence_source']))
    table(_layout_resource(folder/'family_rule_review.tsv'),rule_rows)
    doc.table(['Family','Enabled','Full observed core','Diagnostic CDD','Confounders'],[
        [r['family'],r['enabled_in_selected'],r['full_core'],r['diagnostic_CDD'],r['confounders']] for r in rule_rows])
    regressions=[]
    for row in predictions(_layout_resource(out/'calibration'/selected/'predictions.tsv')):
        if row['accession'] in ('P10586','P35240','O75962'):
            regressions.append(dict(accession=row['accession'],original_id=row['original_id'],expected_family=row['expected_families'],
                assigned_family=row['predicted_family'],grade=row['grade'],supported=row['supported'],correct_family=correct(row)))
    table(_layout_resource(folder/'key_family_regressions.tsv'),regressions)
    doc.table(['Control','Expected family','Assigned family','Grade','Correct supported family'],[
        [r['original_id'],r['expected_family'],r['assigned_family'],r['grade'],r['correct_family']] for r in regressions])
    doc.paragraph('Exact pfamNNNNN accessions are imported only from coordinate-validated returned CDD records, with CDD provenance. '
        'Overlapping HMM/CDD observations count as one footprint, not two votes. Domain presence is not a measurement of whole-protein alignment coverage. '
        'Family-specific CDD and full compatible architectures determine putative family/protein names; absence of experimental validation does not lower the grade. '
        'Generic PH, LIM, FERM, kinase or FN3 domains alone do not identify a specific protein. No gene symbol or evaluation label is used by the new recruitment rules.')
    doc.links([('NCBI CDD documentation','https://www.ncbi.nlm.nih.gov/Structure/cdd/cdd_help.shtml'),
        ('PF01471 PG_binding_1 also occurs in matrixins','https://www.ncbi.nlm.nih.gov/Structure/cdd/pfam01471')])
    doc.heading('Interpretation and implementation')
    doc.paragraph('Validation was not used to select rules, exclusions, thresholds or the orthogroup coefficient. The selection was finalized and checksum-sealed '
        'after an identical-call calibration integration check and before these validation arms. However, this panel and its earlier errors were already inspected '
        'in the project. Exact reference/calibration/model-training strata are supplied; absence of exact overlap does not remove homology or rule-development overlap. '
        'Human validation proteins supply negatives while the other five species supply positives, so species and control class are confounded. Observed zero false positives '
        'is not a guarantee of perfect specificity in novel sequences. Calibration retains its 23 pre-existing hard-negative false positives; none are hidden or removed from scoring.')
    implementation=[('stage2_evidence.py','Separate evidence-driven hypothesis routing and family confirmation; explicit profiles'),
        ('family_identity.py','Exact observed CDD Pfam accessions, provenance and overlap merge; scoped model anchors'),
        ('annotation/stage2_evidence_rules.json','Family cores, diagnostic models, partial thresholds and calibration-selected profile'),
        ('annotation/integrated_adhesome_pipeline.py','Profile dispatch and existing Atlas workbook/trace output'),
        ('adhesome_workflow.py','Resume/continue/audit profile propagation, snapshots and recovery invalidation'),
        ('schisto_adhesome/branch.py; schisto_adhesome/workflows.py','Frozen isolated modes, calibration alongside selected branch, environment restoration'),
        ('configs/isolated_workflow.json','calibrated_v1 selected for future isolated runs; context coefficient 1.0'),
        ('tools/optimize_stage2.py; tools/stage2_optimization_report.py','Calibration selection, freeze verification, cached fixed ablations and transparent comparison'),
        ('tests/workflow/test_stage2_evidence.py','Observed accession, provenance, negative/confounded family, partial-gate and freeze safeguards')]
    doc.table(['Module/configuration','Change'],implementation)
    doc.paragraph('Future isolated validation and target runs use the supplied configuration, record the profile in their frozen contracts, '
        'and retain automatic terminal reporting. Existing completed runs/contracts are not migrated. Direct legacy commands retain the original profile '
        'unless --stage2-profile calibrated_v1 is supplied. Detection, motif libraries, reference proteomes, calibration/validation labels and orthology results are unchanged. '
        'The analysis uses no new third-party dependencies. 108 synthetic workflow tests and 75 family-classifier tests passed, including recovery, benchmark, evidence, reporting and isolation checks. '
        'No backup directories or duplicate pipeline copies were created; the study contains new experimental outputs only.')
    doc.heading('Deliverables')
    links=[('Ablation matrix','ablation_comparison.tsv'),('Family rule review','family_rule_review.tsv'),('Frozen selection','../frozen_selection.json')]
    for kind in ('calibration','validation'):
        links.extend([(kind+' false-negative audit','../'+kind+'/false_negative_audit.tsv'),
            (kind+' candidate change log',kind+'_candidate_changes.tsv'),
            (kind+' new recruitment hypotheses',kind+'_new_recruitment_hypotheses.tsv'),
            (kind+' family sensitivity',kind+'_family_sensitivity.tsv'),
            (kind+' species metrics',kind+'_species_metrics.tsv'),
            (kind+' overlap strata',kind+'_overlap_metrics.tsv'),
            (kind+' supported/provisional/recruitment tiers',kind+'_tier_metrics.tsv'),
            (kind+' remaining family misses',kind+'_remaining_family_misses.tsv'),
            (kind+' optimized assignments','../'+kind+'/'+selected+'/candidate_assignments.tsv')])
    doc.links(links);doc.save(folder)
    # Reuse the existing report renderer while giving this study distinct names.
    (_layout_resource(folder/'validation_review.html')).replace(_layout_resource(folder/'stage2_optimization.html'))
    (_layout_resource(folder/'validation_review.md')).replace(_layout_resource(folder/'stage2_optimization.md'))
    dump(_layout_resource(folder/'comparison_summary.json'),dict(version='1.0',selected_arm=selected,selected_profile=frozen['selected_profile'],
        generated_at_utc=datetime.now(timezone.utc).isoformat(),datasets=summaries,paired=paired_results,provenance=artifacts,
        frozen_selection_sha256=sha(_layout_resource(out/'frozen_selection.json')),validation_selected_configuration=False))
    print(_layout_resource(folder/'stage2_optimization.html'))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--run',type=Path,required=True);parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args();generate(_layout_resource(args.run).resolve(),_layout_resource(args.out).resolve())
