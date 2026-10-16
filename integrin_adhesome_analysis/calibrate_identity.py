"""Fit empirical identity correctness from labelled controls and scored audit rows.

Gold names/genes resolve positive-control evaluation labels, never classifier
inputs. Hard/background negatives are exclusion gates, not identity training data.
"""
import argparse
from collections import Counter,defaultdict
import csv
import hashlib
import json
from pathlib import Path
import re
from identity_probability import ROOT,sha,fingerprint,features,estimate,DEFAULT_MODEL
from classification_contract import general_family


def fasta(path):
    result={};key=None;parts=[]
    for line in path.read_text().splitlines():
        if line.startswith('>'):
            if key:result[key]=''.join(parts).upper().removesuffix('*')
            identifier=line[1:].split()[0];pieces=identifier.split('|');key=pieces[1] if len(pieces)>1 else identifier;parts=[]
            if key in result:raise ValueError('Duplicate calibration accession: '+key)
        else:parts.append(line.strip())
    if key:result[key]=''.join(parts).upper().removesuffix('*')
    return result


def expected_families(record,catalogue):
    """Explicit family ontology for independent gold names; no predicted fields."""
    genes=set(re.split(r'[;,\s]+',str(record.get('genes') or '').upper()))-{''}
    name=str(record.get('protein_name') or '');expected=set()
    for gene in genes:expected.update(catalogue.get(gene,set()))
    patterns={
        'talin':r'\btalin', 'kindlin':r'kindlin|FERMT', 'alpha_actinin':r'alpha.actinin',
        'filamin':r'filamin', 'spectrin_like':r'\bspectrin', 'actin':r'\bactin(?:\b|.related)',
        'cofilin':r'cofilin|actin.depolymeriz', 'profilin':r'profilin',
        'vinculin':r'vinculin|catenin alpha', 'paxillin_like':r'paxillin|leupaxin|paxillin.homolog',
        'zyxin_like':r'zyxin|ajuba|LIM domain.containing protein 1|Wilms tumor protein 1.interacting',
        'zyxin_ajuba_like':r'ajuba|LIM domain.containing protein 1|Wilms tumor protein 1.interacting',
        'PDLIM_family':r'PDZ and LIM|PDZ.LIM', 'PINCH_like':r'PINCH|LIM and senescent',
        'parvin_like':r'parvin', 'FAK':r'focal adhesion kinase|protein.tyrosine kinase 2',
        'ILK':r'integrin.linked', 'Syk':r'\bSyk\b|ZAP.70', 'Nck':r'\bNck',
        'Grb2':r'growth factor receptor.bound protein 2|\bGrb2', 'Crk':r'\bCrk', 'Shc':r'\bShc',
        'HIP1_HIP1R_like':r'Huntingtin.interacting protein 1',
        'integrin_alpha':r'integrin alpha', 'integrin_beta':r'integrin beta',
        'collagen_like':r'collagen', 'fibrillar_collagen_like':r'collagen alpha.*\((?:I|II|III|V|XI)\)',
        'collagen_IV_like':r'collagen.*\(IV\)', 'collagen_V_like':r'collagen.*\(V\)',
        'laminin_like':r'laminin', 'netrin_like':r'netrin', 'fibronectin_like':r'\bfibronectin\b',
        'vitronectin_like':r'vitronectin', 'PTPN1_PTP1B_like':r'tyrosine.protein phosphatase non.receptor type 1\b',
        'PTPRN_PTPRN2_like':r'tyrosine.protein phosphatase N(?:2)?\b',
        'LAR_RPTP':r'receptor.type tyrosine.protein phosphatase (?:F|D|S)\b',
        'RPTP_alpha_epsilon':r'receptor.type tyrosine.protein phosphatase (?:alpha|epsilon)\b',
        'ROBO_like':r'roundabout', 'DCC_neogenin_like':r'neogenin|deleted in colorectal',
        'UNC89_obscurin_SPEG_like':r'obscurin|\bSPEG\b|UNC.89', 'Eph_receptor_like':r'ephrin type.*receptor|\bEph.*receptor',
        'twinfilin_like':r'twinfilin', 'drebrin_like':r'drebrin', 'GMF_like':r'glia maturation',
        'PET_LIM_like':r'prickle|testin'}
    for family,pattern in patterns.items():
        if re.search(pattern,name,re.I):expected.add(family)
    for k,pattern in [('Src',r'SRC|FYN|LYN|YES1|HCK|LCK|BLK|FGR'),('Csk',r'CSK|MATK'),
        ('Abl_family',r'ABL[12]'),('Tec_family',r'TEC|BTK|BMX|ITK|TXK'),('PTP_PEST',r'PTPN12|PTPN18|PTPN22')]:
        if any(re.fullmatch(pattern,gene) for gene in genes):expected.add(k)
    if expected&{'Src','Csk'}:expected.add('Src_Csk_like')
    if any(g.startswith('PDLIM') for g in genes) or genes&{'LDB3','ZASP'}:expected.add('PDLIM_family')
    if genes&{'AJUBA','JUB','LIMD1','WTIP'}:expected.update(['zyxin_like','zyxin_ajuba_like'])
    if genes&{'LIMS1','LIMS2'}:expected.add('PINCH_like')
    if genes&{'CFL1','CFL2','DSTN'}:expected.add('cofilin')
    if expected&{'actin'}:
        subtype=re.search(r'actin.related protein\s*(1|2|3|4|5|6|8|10)(?:\b|[AB])',name,re.I)
        if subtype:expected.add('arp'+subtype[1]+'_like')
        elif 'related' not in name.lower():expected.add('actin_conventional')
    for chain in ('alpha','beta','gamma'):
        if re.search('laminin.*'+chain,name,re.I):expected.add('laminin_'+chain+'_like')
    if not expected and name.strip() and not re.search(r'^uncharacterized|^hypothetical|^unknown',name,re.I):
        expected.add('other_named_identity_outside_current_family_catalogue')
    return expected


def fit(calibration_set,feature_path,out,report_dir,minimum=20):
    from openpyxl import load_workbook
    books=list(calibration_set.glob('*.xlsx'));fastas=list(calibration_set.glob('*.faa'))
    if len(books)!=1 or len(fastas)!=1:raise ValueError('Calibration folder requires one XLSX label book and one .faa file')
    book=load_workbook(books[0],read_only=True,data_only=True);values=iter(book['control_identity'].values)
    headings=next(values);labels=[dict(zip(headings,row)) for row in values if any(x is not None for x in row)];book.close()
    if len({row['accession'] for row in labels})!=len(labels):raise ValueError('Calibration labels have duplicate accessions')
    seqs=fasta(fastas[0]);sequence_hash={key:hashlib.sha256(seq.encode()).hexdigest() for key,seq in seqs.items()}
    with feature_path.open(encoding='utf-8-sig',newline='') as handle:scored=list(csv.DictReader(handle,delimiter='\t'))
    by_hash=defaultdict(list)
    for row in scored:
        if row.get('sequence_sha256'):by_hash[row['sequence_sha256']].append(row)
    config=json.loads((ROOT/'01_detection/config/family_rules.json').read_text());catalogue=defaultdict(set);ref_hashes=set()
    with (ROOT/'01_detection/references/adhesome_v1.0/search_catalogue.csv').open(encoding='utf-8-sig',newline='') as handle:
        for row in csv.DictReader(handle):
            ref_hashes.add(row['sequence_sha256'])
            if row['family'] not in config['families']:continue
            for gene in re.split(r'[;,\s]+',row.get('gene','').upper()):
                catalogue[gene].add(general_family(row['family'],config['families']))
    samples=[];resolved=[];excluded=Counter();used_hashes=set()
    for gold in labels:
        accession=gold['accession'];digest=sequence_hash.get(accession)
        status='used';prediction={};truth=expected_families(gold,catalogue)
        if digest is None:status='missing_fasta_sequence'
        elif 'negative' in str(gold['class']).lower():status='negative_exclusion_control_not_identity_training'
        elif 'positive' not in str(gold['class']).lower():raise ValueError('Unknown control class: '+str(gold['class']))
        elif digest not in by_hash:status='no_scored_audit_evidence'
        elif digest in used_hashes:status='duplicate_sequence_group'
        else:
            rows=by_hash[digest]
            signatures={(r['assigned_family'],r['total_score_0_10'],r['family_confidence']) for r in rows}
            if len(signatures)!=1:raise ValueError('Conflicting scored features for '+accession)
            prediction=rows[0]
            if prediction['assigned_family']=='unassigned':status='classifier_abstained'
            elif not truth:status='truth_needs_family_mapping_review'
            else:
                sample=dict(features(prediction),correct=int(prediction['assigned_family'] in truth),accession=accession)
                samples.append(sample);used_hashes.add(digest)
        if status!='used':excluded[status]+=1
        resolved.append(dict(accession=accession,control_class=gold['class'],trusted_protein_name=gold['protein_name'],
            trusted_genes=gold['genes'],expected_families=';'.join(sorted(truth)),
            called_family=prediction.get('assigned_family',''),status=status,
            correct_identity_at_reported_scope=prediction.get('assigned_family') in truth if status=='used' else '',
            sequence_sha256=digest or '',exact_detection_reference_overlap=digest in ref_hashes))
    if len(samples)<minimum:raise ValueError(f'Only {len(samples)} scored, sequence-verified calibration controls; need {minimum}')
    cv=[]
    for query in samples:
        train=[sample for sample in samples if sample['sequence_sha256']!=query['sequence_sha256']]
        result=estimate(train,query,minimum)
        if 'probability' in result:cv.append((query['correct'],result['probability']))
    hashes={str(path.resolve().relative_to(ROOT)).replace('\\','/'):sha(path) for path in [books[0],fastas[0],Path(__file__)]}
    context=dict(label_rows=len(labels),fasta_sequences=len(seqs),scored_controls=len(samples),excluded=dict(excluded),
        exact_detection_reference_overlap=len(set(sequence_hash.values())&ref_hashes),
        exact_current_scored_sequence_overlap=len(set(sequence_hash.values())&set(by_hash)),
        reference_panels=sorted({sample.get('orthology_panel_sha256') for sample in samples}-{None,''}),
        positive_control_only=True,
        positive_family_scope=sorted({family for gold in labels if 'positive' in str(gold['class']).lower() for family in expected_families(gold,catalogue)}),
        interpretation='Empirical identity correctness conditional on positive adhesome controls and negative exclusion gates. Negative controls are excluded from identity calibration. Reference-overlapping controls; not independent validation for novel lineages.')
    model=dict(pipeline_version='1.0',method='isotonic_empirical_control_precision',minimum_controls=minimum,
        code_fingerprint=fingerprint(),input_hashes=hashes,scored_features_sha256=sha(feature_path),
        context=context,samples=samples,cross_validation=dict(method='leave_one_exact_sequence_group_out',
            evaluated=len(cv),brier_score=sum((y-p)**2 for y,p in cv)/len(cv) if cv else None,
            mean_predicted_probability=sum(p for y,p in cv)/len(cv) if cv else None,
            observed_correct_fraction=sum(y for y,p in cv)/len(cv) if cv else None))
    out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(model,indent=2)+'\n',encoding='utf-8')
    report_dir.mkdir(parents=True,exist_ok=True)
    with (report_dir/'control_identity_review.tsv').open('w',encoding='utf-8',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(resolved[0]),delimiter='\t');writer.writeheader();writer.writerows(resolved)
    (report_dir/'calibration_report.json').write_text(json.dumps({k:v for k,v in model.items() if k!='samples'},indent=2)+'\n',encoding='utf-8')
    return {k:v for k,v in model.items() if k in ('context','cross_validation')}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--calibration-set',type=Path,default=ROOT/'02_calibration_set');p.add_argument('--features',type=Path,required=True)
    p.add_argument('--out',type=Path,default=DEFAULT_MODEL);p.add_argument('--report-dir',type=Path,required=True)
    p.add_argument('--minimum-controls',type=int,default=20);args=p.parse_args()
    if args.minimum_controls<2:p.error('--minimum-controls must be at least two')
    print(json.dumps(fit(args.calibration_set.resolve(),args.features,args.out,args.report_dir,args.minimum_controls),indent=2))


if __name__=='__main__':main()
