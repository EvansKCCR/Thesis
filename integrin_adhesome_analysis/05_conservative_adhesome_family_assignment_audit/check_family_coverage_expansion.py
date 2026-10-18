#!/usr/bin/env python3
"""Static consistency checks of the reference-guided family expansion. Not a Stage 2 run."""
import json,csv,re,collections,hashlib
from pathlib import Path
from artifact_tool import Blob,SpreadsheetFile
p=Path(__file__).resolve().parent
f0=json.load(open(p/'family_rules.json'));i0=json.load(open(p/'family_identity_rules_recalibrated.json'))
f1=json.load(open(p/'family_rules_expanded.json'));i1=json.load(open(p/'family_identity_rules_expanded.json'))
old_f=f0['families'];new_f=f1['families'];old_i={r['id']:r for r in i0['rules']};new_i={r['id']:r for r in i1['rules']}
assert all(new_f[k]==v for k,v in old_f.items()), 'Modified an existing family'
assert all(new_i[k]==v for k,v in old_i.items()), 'Modified an existing identity rule'
assert len(new_i)==len(i1['rules'])
for k in new_f:
 assert new_f[k]['classification_level'] in ('family','screen','subfamily'),k
for r in i1['rules']:
 re.compile(r.get('specific_cdd_pattern','a'))
 assert r.get('required_groups') is not None, r['id']
if len(i1['rules'])!=len(old_i)+sum(x['decision']=='ADD_IDENTITY_PROPOSAL' for x in csv.DictReader(open(p/'family_coverage_expansion_audit.tsv'),delimiter='\t')): raise AssertionError('Unequal new identity rule counts')
audit=list(csv.DictReader(open(p/'family_coverage_expansion_audit.tsv'),delimiter='\t'))
accepted=[r for r in audit if r['decision'].startswith('ADD_')]
assert len(accepted)==len(new_f)-len(old_f)
assert all(int(r['n_arch_plus_cdd'])==int(r['n_positive']) and not int(r['n_other_positive_cross_hits']) and not int(r['n_negative_cross_hits']) for r in accepted)
assert len(new_f)==len(old_f)+len(accepted)
w=SpreadsheetFile.import_xlsx(Blob.load(str(p/'calibration_set.xlsx')))
raw=w.worksheets.get_item('control_identity').get_range('A1:K584').values
cal={str(r[1]):dict(zip(raw[0],r)) for r in raw[1:]}
assert len(cal)==583
for r in accepted:
 fam=r['expected_family'];idx='cal_'+re.sub(r'[^a-z0-9]+','_',fam.lower()).strip('_');
 if r['decision']=='ADD_ARCHITECTURE_SCREEN':
  assert new_f[fam]['screening_only'] is True and idx not in new_i
  continue
 rr=new_i[idx]
 assert rr['family']==fam
 assert rr['allow_diagnostic_dispatch']
 assert rr['required_groups'] and all(x['count']>=1 for x in rr['required_groups'])
 for ac in r['control_accessions'].split(';'):
  obs=cal[ac]; pf=set(re.findall('PF\\d{5}',str(obs['pfam'])))
  assert all(any(m in pf for m in g['models']) for g in rr['required_groups']), (fam,ac)
  assert re.search(rr['specific_cdd_pattern'],str(obs['NCBI_CDD']),re.I), (fam,ac)
print(json.dumps({'old_family_rules':len(old_f),'new_family_rules':len(new_f),'old_identity_rules':len(old_i),'new_identity_rules':len(new_i),'added_families':len(accepted),'identity_proposals':sum(z['decision']=='ADD_IDENTITY_PROPOSAL' for z in accepted),'screen_only':sum(z['decision']=='ADD_ARCHITECTURE_SCREEN' for z in accepted),'static_positive_reference_matches':sum(int(r['n_arch_plus_cdd']) for r in accepted),'new_rule_negative_cross_matches':sum(int(r['n_negative_cross_hits']) for r in accepted),'new_rule_other_family_cross_matches':sum(int(r['n_other_positive_cross_hits']) for r in accepted),'existing_family_rules_preserved':True,'existing_identity_rules_preserved':True},indent=2))
