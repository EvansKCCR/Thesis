#!/usr/bin/env python3
"""Evidence-gated rule expansion from the frozen calibration workbook.
Builds proposed JSON configurations and auditable static screens; does NOT run Stage 2.
Uses artifact_tool to inspect .xlsx (no openpyxl/pandas).
"""
from artifact_tool import Blob, SpreadsheetFile
import json,csv,re,copy,hashlib,collections,math
from pathlib import Path
P=Path('/mnt/data')
W=SpreadsheetFile.import_xlsx(Blob.load(str(P/'calibration_set.xlsx')))
tab=W.worksheets.get_item('control_identity').get_range('A1:K584').values
cal={str(r[1]):dict(zip(tab[0],r)) for r in tab[1:]}
maprows=list(csv.DictReader(open(P/'family_mapping_review.tsv',encoding='utf-8-sig'),delimiter='\t'))
missing=collections.defaultdict(list)
for r in maprows:
 if r['family_rule_coverage']=='outside_configured_family_gates':missing[r['expected_families']].append(r['accession'])
orig_f=json.load(open(P/'family_rules.json'))
orig_i=json.load(open(P/'family_identity_rules_recalibrated.json'))

def S(model_groups,cd,**kw):return dict(groups=model_groups,cdd=cd,**kw)
# Architecture is specified in terms of observed Pfam models, and MUST be supported by a
# second CDD observation. The named CDD signatures are corroboration, not experimental gold.
SPEC={
# multi-control families
'RhoGEF_like':S([['PF00621']],r'CH_(?:alphaPIX|betaPIX|VAV[123])|PH1_Kalirin_Trio_like|RGS_LARG|PH_ARHGEF2',module='downstream_signalling',caveat='Broad RhoGEF domain classification; TIAM/SOS refinements take precedence'),
'RhoGAP_like':S([['PF00620']],r'RhoGAP_(?:Graf|p190|ARHGAP|DLC|CdGAP)|RhoGAP\s*\[|RhoGAP_p',module='downstream_signalling'),
'ASAP_AGAP_like':S([['PF01412'],['PF12796']],r'ArfGap_(?:ASAP|AGAP)|PH_(?:ASAP|AGAP)',module='downstream_signalling'),
'PKC_like':S([['PF00069'],['PF00130']],r'(?:C1_[cn]PKC|C2_PKC|STKc_[cn]PKC)',module='downstream_signalling'),
'tubulin_like':S([['PF03953']],r'alpha_tubulin\b',module='actin_linking',caveat='Alpha-tubulin superfamily, not exact isoform'),
'SORBS_like':S([['PF02208']],r'SH3_(?:Sorbs[123]|Vinexin)|\bSorb\s*\[',module='focal_adhesion_scaffolds'),
'ABI_like':S([['PF07815']],r'Abi_HHR|SH3_Abi',module='focal_adhesion_scaffolds'),
'GIT_like':S([['PF01412'],['PF08518']],r'ArfGap_GIT[12]',module='downstream_signalling'),
'INPP5_SHIP_like':S([['PF00017']],r'SH2_SHIP|INPP5c_SHIP',module='downstream_signalling'),
'MAPK_like':S([['PF00069']],r'STKc_(?:ERK1_2_like|JNK)',module='downstream_signalling'),
'PPM_like':S([['PF00481']],r'PP2Cc',module='downstream_signalling',caveat='Broad PP2C phosphatase family; generic model requires competitors'),
'RNF5_185_like':S([['PF13920','PF00097']],r'RING-HC_RNF(?:5|185)\b',module='downstream_signalling'),
'Rho_GTPase_like':S([['PF00071']],r'(?:Wrch_1|RhoA_like)',module='downstream_signalling',caveat='Broad RHO identity; other Ras subfamilies must compete'),
'calpain_like':S([['PF00648'],['PF01067']],r'Calpain_III|CysPc',module='actin_linking'),
'neuropilin_like':S([['PF00431'],['PF00629'],['PF00754']],r'MAM|CUB',module='extracellular_matrix'),
'parvin_like':S([['PF00307']],r'CH_PARV[AB]',module='membrane_proximal_adaptors'),
'syndecan_like':S([['PF01034']],r'Syndecan',module='integrin_receptors'),
 'tensin_like':S([['PF10409']],r'PTB_tensin|SH2_Tensin_like|PTP_tensin',module='focal_adhesion_scaffolds'),
 'tetraspanin_like':S([['PF00335']],r'Tetraspanin|CD151_like_LEL|uroplakin_I_like_LEL',module='integrin_receptors'),
# singleton with distinctive architecture or curated CDD model
'ADAM_like':S([['PF01421'],['PF00200']],r'ZnMc_adamalysin|\bADAM_CR\b',module='extracellular_matrix'),
'AKT_like':S([['PF00069'],['PF00169']],r'STKc_PKB|PH_PKB',module='downstream_signalling'),
'ARF_GTPase_like':S([['PF00025']],r'Arf1_5_like',module='downstream_signalling'),
'ARPC2_like':S([['PF04045']],r'P34-Arc',module='actin_linking'),
'BCAR3_like':S([['PF00017'],['PF00617']],r'SH2_BCAR3',module='focal_adhesion_scaffolds'),
'CBL_like':S([['PF02262'],['PF02761'],['PF02762']],r'RING-HC_Cbl|SH2_Cbl',module='downstream_signalling'),
'CD47_like':S([['PF04549'],['PF08204']],r'IgV_CD47',module='integrin_receptors'),
'CEACAM_like':S([['PF07686'],['PF13927']],r'IgV_CEACAM',module='integrin_receptors'),
'DOCK_like':S([['PF16172'],['PF14429']],r'SH3_DOCK|C2_Dock',module='downstream_signalling'),
'ELMO_like':S([['PF04727'],['PF11841']],r'PH_ELMO|ELMO_CED12',module='downstream_signalling'),
'FABP_like':S([['PF00061']],r'FABP3',module='downstream_signalling',caveat='FABP family evidence specific to reference, not lineage generalization'),
'GAB_like':S([['PF00169']],r'PH_Gab1_Gab2',module='downstream_signalling',caveat='Cannot use generic PH domain alone'),
'GRB7_family':S([['PF00017'],['PF00169']],r'RA_GRB7',module='downstream_signalling'),
'HSP70_like':S([['PF00012']],r'ASKHA_NBD_HSP70|\bHSP70\b',module='downstream_signalling'),
'ICAP1_like':S([['PF10480']],r'PTB_ICAP1|ICAP-1_inte_bdg',module='membrane_proximal_adaptors'),
'INSR_like':S([['PF01030'],['PF07714']],r'PTKc_INSR|IR_IRR|InsR|Recep_L_domain',module='integrin_receptors'),
'IRS_like':S([['PF00169'],['PF02174']],r'PH_IRS|PTB_IRS',module='downstream_signalling'),
'ITGB3BP_like':S([['PF06729']],r'CENP-R',module='membrane_proximal_adaptors',caveat='CENPR own identity, not integrin-beta identity'),
'KCNH_like':S([['PF00520'],['PF13426']],r'PAS|Ion_trans',module='integrin_receptors'),
'LDLR_related_like':S([['PF00057'],['PF00058']],r'LDLa',module='extracellular_matrix',caveat='LDLR superfamily screening only unless independently corroborated'),
'LIMK_like':S([['PF00412'],['PF00595']],r'LIM1_LIMK|PDZ_LIMK',module='downstream_signalling'),
'MAPK8IP_like':S([['PF09744'],['PF16471']],r'Jnk-SapK_ap_N|JIP_LZII',module='downstream_signalling'),
'MARCKS_like':S([['PF02063']],r'MARCKS',module='membrane_proximal_adaptors'),
'MMP_like':S([['PF00413'],['PF01471']],r'ZnMc_MMP',module='extracellular_matrix'),
'NDEL1_like':S([['PF04880']],r'NUDE_C',module='actin_linking'),
'NF2_merlin_like':S([['PF00373'],['PF00769']],r'FERM_F1_Merlin',module='membrane_proximal_adaptors'),
'OSTF1_like':S([['PF12796'],['PF00018','PF14604']],r'SH3_OSTF1',module='focal_adhesion_scaffolds'),
'PABP_like':S([['PF00658']],r'RRM1_I_PABPs|PABP-1234',module='downstream_signalling'),
'PAK_like':S([['PF00786'],['PF00069']],r'STKc_PAK|CRIB_PAK',module='downstream_signalling'),
'PDE4_like':S([['PF00233'],['PF18100']],r'PDE4_UCR',module='downstream_signalling'),
'PDPK1_like':S([['PF00069'],['PF14593']],r'STKc_PDK1|PH_PDK1',module='downstream_signalling'),
'PI3K_like':S([['PF00454'],['PF00794']],r'PI3K',module='downstream_signalling'),
'PIP5K_like':S([['PF01504']],r'PIPKc_PIP5K1C|\bPIP5K\b',module='downstream_signalling'),
'PKA_like':S([['PF00069']],r'STKc_PKA',module='downstream_signalling'),
'PLAUR_like':S([['PF00021']],r'uPAR_rpt|TFP_LU_ECD_uPAR',module='integrin_receptors'),
'PLCG_like':S([['PF00387'],['PF00388']],r'PI-PLCc_gamma|PH_PLC_gamma',module='downstream_signalling'),
'PLD_like':S([['PF00787'],['PF00614']],r'PLDc_vPLD1_2',module='downstream_signalling'),
'PPP2_like':S([['PF00149']],r'MPP_PP2A_PP4_PP6',module='downstream_signalling',caveat='Broad PP2A/4/6, not PP2A-specific without extra evidence'),
'PRKD_like':S([['PF00069'],['PF00130']],r'C1_PKD1|PH_PKD|STKc_PKD',module='downstream_signalling'),
'PTEN_like':S([['PF10409']],r'PTP_PTEN',module='downstream_signalling'),
'PTPN2_like':S([['PF00102']],r'PTPc-N2',module='downstream_signalling'),
'PTPRH_like':S([['PF00041'],['PF00102']],r'R-PTPc-H|PTPRH',module='downstream_signalling'),
'PTPRO_like':S([['PF00041'],['PF00102']],r'R-PTPc-O|PTPRO',module='downstream_signalling'),
'RACK1_like':S([['PF00400']],r'RACK1',module='membrane_proximal_adaptors',caveat='Generic WD40 alone cannot determine RACK1'),
'RAPGEF_like':S([['PF00617'],['PF00618']],r'REM|RasGEF',module='downstream_signalling',caveat='Could overlap SOS; must resolve with subtype CDD'),
'ROCK_like':S([['PF00069'],['PF08912']],r'STKc_ROCK|HR1_ROCK',module='downstream_signalling'),
'RasGAP_like':S([['PF00616']],r'SH2_Nterm_RasGAP|SH3_RasGAP',module='downstream_signalling'),
'Ras_GTPase_like':S([['PF00071']],r'H_N_K_Ras_like',module='downstream_signalling'),
'SH2B_like':S([['PF08916'],['PF00017']],r'SH2_SH2B_family|PH_SH2B_family',module='downstream_signalling'),
'SH3KBP1_like':S([['PF14604']],r'SH3_CIN85_',module='focal_adhesion_scaffolds'),
'SHARPIN_like':S([['PF16764']],r'PH_SHARPIN|Sharpin_PH',module='focal_adhesion_scaffolds'),
'SIRP_like':S([['PF07654','PF07686']],r'IgV_SIRP|IgC1_SIRP',module='integrin_receptors'),
'SLC16_like':S([['PF07690']],r'MFS_MCT3_4',module='integrin_receptors'),
'SLC3_like':S([['PF16028']],r'AmyAc_SLC3A2|SLC3A2_N',module='integrin_receptors'),
'SLC9_like':S([['PF00999']],r'Na_H_Exchanger',module='integrin_receptors',caveat='Transporter superfamily; exact lineage unresolved'),
'SMPX_like':S([['PF15355']],r'Chisel',module='actin_linking'),
'SOS_like':S([['PF00617'],['PF00618']],r'HFD_SOS1|PH_SOS',module='downstream_signalling'),
'STAT_like':S([['PF01017'],['PF02864']],r'STAT3_(?:CCD|DBD)|STAT_alpha',module='downstream_signalling'),
'TIAM_like':S([['PF00621'],['PF02196']],r'PH1_Tiam|PH_Tiam|RhoGEF_Tiam',module='downstream_signalling'),
'TRPM_like':S([['PF16519'],['PF18139']],r'Alpha_kinase_ChaK1_TRMP7|TRPM_tetra',module='integrin_receptors'),
'ZFYVE21_like':S([['PF01363'],['PF16696']],r'FYVE_ZF21|ZFYVE21_C',module='membrane_proximal_adaptors'),
'caspase_like':S([['PF00656'],['PF01335']],r'DED_Caspase_8|CASc',module='downstream_signalling'),
'caveolin_like':S([['PF01146']],r'Caveolin',module='membrane_proximal_adaptors'),
'cytohesin_like':S([['PF01369'],['PF00169']],r'PH_GRP1-like|\bSec7\b',module='downstream_signalling'),
'defensin_like':S([['PF00323'],['PF00879']],r'Defensin_propep|Defensin_1',module='extracellular_matrix'),
'dynamin_like':S([['PF00350'],['PF01031']],r'DLP_1|Dynamin_N',module='downstream_signalling'),
'intermediate_filament_like':S([['PF00038'],['PF04732']],r'Filament_head',module='actin_linking'),
'kinectin_like':S([['PF05104']],r'Rib_recp_KP_reg',module='actin_linking'),
'liprin_like':S([['PF00536']],r'SAM_liprin',module='focal_adhesion_scaffolds'),
'migfilin_like':S([['PF00412']],r'LIM2_FBLP',module='focal_adhesion_scaffolds'),
'myosin_like':S([['PF00063'],['PF01576']],r'Myosin_head|MYSc',module='actin_linking'),
'nectin_PVR_like':S([['PF08205']],r'IgV_1_Nectin|IgC1_2_Nectin',module='integrin_receptors'),
'nischarin_like':S([['PF00787']],r'PX_IRAS',module='focal_adhesion_scaffolds'),
'palladin_like':S([['PF07679']],r'IgI_1_Palladin_C',module='actin_linking'),
'plakin_like':S([['PF03501'],['PF00307']],r'CH_PLEC|S10_plectin',module='actin_linking'),
'prion_like':S([['PF11587']],r'Prion_bPrPp|\bPRP\s*\[',module='extracellular_matrix'),
'slingshot_like':S([['PF00782'],['PF08766']],r'DSP_slingshot|SSH-N',module='downstream_signalling'),
'small_HSP_like':S([['PF00011']],r'ACD_HspB1_like',module='actin_linking'),
'syntenin_like':S([['PF00595']],r'PDZ[12]_syntenin',module='membrane_proximal_adaptors'),
}
# Safety constraint: predicates are measured on the reference panel and must not be misrepresented as independent thresholds.
assert not set(SPEC).intersection(orig_f['families'])
for fam in SPEC:assert fam in missing,fam

def pfams(entry):return set(re.findall(r'PF\d{5}',str(entry['pfam'])))
def match_arch(spec,entry):
 models=pfams(entry)
 return all(any(m in models for m in g) for g in spec['groups'])
def match_cdd(spec,entry):return bool(re.search(spec['cdd'],str(entry['NCBI_CDD']),flags=re.IGNORECASE))
def targetmatch(spec,entry):return match_arch(spec,entry) and match_cdd(spec,entry)
labels={r['accession']:r for r in maprows}
positives=set(labels)
records=[]
accepted={}
for family,sp in SPEC.items():
 targets=missing[family];targetpos=[a for a in targets if targetmatch(sp,cal[a])]
 nonpos=[a for a in cal if a not in targets and a in positives and targetmatch(sp,cal[a])]
 neg=[a for a in cal if a not in positives and targetmatch(sp,cal[a])]
 coverage=len(targetpos)/len(targets)
 # Require every reference example match both architecture and independent CDD signature. Allow 1 unmatched only for 4+ controls with missing CDD evidence; but those get deferred.
 is_supported=(len(targetpos)==len(targets) and not neg and not nonpos)
 named_cd_all=all(any(re.match(r'^\s*cd\d+',seg,re.I) and re.search(sp['cdd'],seg,re.I) for seg in str(cal[a]['NCBI_CDD']).split(';')) for a in targets)
 decision=('ADD_IDENTITY_PROPOSAL' if named_cd_all else 'ADD_ARCHITECTURE_SCREEN') if is_supported else ('DEFER_INCOMPLETE_REFERENCE_CDD' if len(targetpos)!=len(targets) else 'DEFER_CROSS_FAMILY_OR_NEGATIVE_MATCH')
 records.append(dict(expected_family=family,n_positive=len(targets),n_arch=sum(match_arch(sp,cal[a]) for a in targets),n_arch_plus_cdd=len(targetpos),n_other_positive_cross_hits=len(nonpos),n_negative_cross_hits=len(neg),control_accessions=';'.join(targets),positive_matched=';'.join(targetpos),other_positive_cross_hits=';'.join(nonpos),negative_cross_hits=';'.join(neg),configured_pfams=';'.join('/'.join(g) for g in sp['groups']),diagnostic_cdd_pattern=sp['cdd'],named_curated_cdd_for_all=str(named_cd_all),decision=decision,caveat=sp.get('caveat','')))
 if is_supported:accepted[family]=sp
# Compile paired reference-based family screens and identity refinement rules.
f=copy.deepcopy(orig_f);i=copy.deepcopy(orig_i)
BASE_PARENTS=['component_reference_candidate','unassigned_domain_candidate']
for fam,sp in accepted.items():
 accs=missing[fam]; lengths=[int(cal[a]['length']) for a in accs]
 min_length=max(45, int(min(lengths)*0.6))
 # The family configuration is non-diagnostic by itself for generic PF models; specific CDD + family competition is essential.
 famdef={
   'module':sp['module'],'min_length':min_length,'diagnostic':False,
   'domain_alternatives':[{m:1 for g in sp['groups'] for m in g[:1]}],
   'min_partial_fraction':0.5,
   'screening_domains':list(dict.fromkeys(m for g in sp['groups'] for m in g)),
   'screening_specific':False,
   'classification_level':('family' if next(x for x in records if x['expected_family']==fam)['decision']=='ADD_IDENTITY_PROPOSAL' else 'screen'),'family_group':fam,
   'requires_identity_evidence':True,'negative_class_competition':True,
   'benchmark_status':'reference_required_not_calibrated',
   'classification_note':'Observed reference-panel Pfam/C DD architecture; requires independent family-specific CDD evidence and competing-family assessment. Not a calibrated novelty threshold. '+sp.get('caveat',''),
   'label':fam
 }
 # represent alternatives properly as Cartesian combinations
 import itertools
 alts=[dict((g[j],1) for g,j in zip(sp['groups'],ix)) for ix in itertools.product(*[range(len(g)) for g in sp['groups']])]
 famdef['domain_alternatives']=alts
 if len(sp['groups'])>1:
  famdef['ordered_domain_groups']=sp['groups']
  famdef['support_requires_order']=False # CDD reports may not retain coordinates for all models
 if next(x for x in records if x['expected_family']==fam)['decision']=='ADD_ARCHITECTURE_SCREEN':
  famdef['screening_only']=True
  famdef['classification_note']+=' Architecture-screen only: all characterized positives lack a matching named curated cd model; do not treat as supported identity.'
 f['families'][fam]=famdef
 if famdef['classification_level']=='screen':continue
 rule_id='cal_'+re.sub(r'[^a-z0-9]+','_',fam.lower()).strip('_')
 assert rule_id not in [z['id'] for z in i['rules']]
 i['rules'].append({'id':rule_id,'parents':BASE_PARENTS + [fam],
   'family':fam, 'protein_name':fam.replace('_',' ').replace('like','-like')+' protein',
   'required_groups':[{'models':g,'count':1} for g in sp['groups']],
   'specific_cdd_pattern':sp['cdd'],'ordered_groups':[],
   'allow_diagnostic_dispatch':True})
i['rule_revision_context']=str(orig_i.get('rule_revision_context',''))+' Expansion adds evidence-gated identity rules only where all family positives match Pfam+specific CDD and no other labelled control does in the frozen calibration panel; retain conservative thresholds and validate runtime.'
(P/'family_rules_expanded.json').write_text(json.dumps(f,indent=2,ensure_ascii=False)+'\n')
(P/'family_identity_rules_expanded.json').write_text(json.dumps(i,indent=2,ensure_ascii=False)+'\n')
with open(P/'family_coverage_expansion_audit.tsv','w',newline='') as h:
 w=csv.DictWriter(h,fieldnames=list(records[0]),delimiter='\t');w.writeheader();w.writerows(records)
# all outside target not covered by accepted proposals
uncovered={fam:accs for fam,accs in missing.items() if fam not in accepted}
print(json.dumps({'candidate_specs':len(SPEC),'accepted':len(accepted),'added_positive_controls':sum(len(missing[k]) for k in accepted),'identity_proposals':sum(x['decision']=='ADD_IDENTITY_PROPOSAL' for x in records),'screen_only':sum(x['decision']=='ADD_ARCHITECTURE_SCREEN' for x in records),'total_outside_controls':sum(len(z) for z in missing.values()),'uncovered_family_sets':len(uncovered),'uncovered_controls':sum(len(z) for z in uncovered.values()),'deferred_requested_specs':[{k:r[k] for k in ('expected_family','n_positive','n_arch_plus_cdd','n_other_positive_cross_hits','n_negative_cross_hits')} for r in records if not r['decision'].startswith('ADD_')],'rules_orig':len(orig_i['rules']),'rules_new':len(i['rules'])},indent=2))
