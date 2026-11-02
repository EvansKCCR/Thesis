"""Trusted-name evaluation ontology, version 1.0. Never a classifier input.

Names outside the configured architecture catalogue remain evaluable labels;
having a label does not create a classification rule or establish acceptance.
Patterns resolve the supplied trusted protein names, not predicted annotations.
"""

# Resolve source resources independently of the caller working directory.
import sys as _layout_sys
from pathlib import Path as _LayoutPath
_LAYOUT_BUNDLE_ROOT = next(p for p in _LayoutPath(__file__).resolve().parents
                          if (p/'config/layout_paths.json').is_file() and (p/'pipeline_paths.py').is_file())
if str(_LAYOUT_BUNDLE_ROOT) not in _layout_sys.path:
    _layout_sys.path.insert(0, str(_LAYOUT_BUNDLE_ROOT))
from pipeline_paths import activate as _layout_activate, resource_path as _layout_resource
_layout_activate()


import re

# More specific names precede broad classes. Labels have family scope, not an
# exact paralogue or an assertion of experimental adhesome membership.
NAME_FAMILIES=[
    # Binding partners are separate identities. These precede the broad actin,
    # integrin, GRB2 and LIM names used by the legacy recruitment catalogue.
    ('GAB_like',r'^GRB2.associated.binding protein [123]\b'),
    ('ICAP1_like',r'^integrin beta.1.binding protein 1\b|^integrin cytoplasmic domain.associated protein 1\b'),
    ('ARPC2_like',r'^actin.related protein 2/3 complex subunit 2\b'),
    ('migfilin_like',r'^filamin.binding LIM protein 1\b|^migfilin\b'),
    ('Testin_like',r'^testin\b'),('Prickle_like',r'^prickle(?:.like)?(?: protein)?(?: [1-4])?\b'),
    ('MACF1_spectraplakin_like',r'microtubule.actin cross.linking factor'),
    ('plakin_like',r'\bplectin\b|\bdystonin\b|\bdesmoplakin\b'),
    ('alpha_actinin',r'alpha.actinin'),
    ('actin',r'^actin(?:$|,|\s+(?:cytoplasmic|aortic|skeletal|cardiac|alpha|beta|gamma|muscle))|^actin.related protein\b(?!.*\bcomplex subunit\b)'),
    ('coronin_like',r'^coronin'),('cortactin_HS1_like',r'\bcortactin\b|hematopoietic lineage cell.specific protein'),
    ('LASP_like',r'^LIM and SH3 domain protein'),('ENA_VASP_like',r'vasodilator.stimulated phosphoprotein|protein enabled homolog|Ena/VASP'),
    ('nexilin_like',r'^nexilin'),('supervillin_like',r'^supervillin'),('SORBS_like',r'^sorbin and SH3|^vinexin'),
    ('BCAR_CAS_like',r'breast cancer anti.estrogen resistance protein 1\b|enhancer of filamentation|Cas scaffolding protein'),
    ('caveolin_like',r'^caveolin'),('SMPX_like',r'^small muscular protein'),
    ('SH3KBP1_like',r'^SH3 domain.containing kinase.binding protein'),
    ('ERM_like',r'^ezrin\b|^radixin\b|^moesin\b'),('NF2_merlin_like',r'^merlin\b'),
    ('FHL_like',r'^four and a half LIM'),('GRB7_family',r'growth factor receptor.bound protein (?:7|10|14)\b'),
    ('liprin_like',r'^liprin'),('palladin_like',r'^palladin'),('RACK1_like',r'\bRACK1\b'),
    ('OSTF1_like',r'^osteoclast.stimulating factor'),('synemin_like',r'^synemin'),('syntenin_like',r'^syntenin'),
    ('tensin_like',r'^tensin'),('NDEL1_like',r'nuclear distribution protein nudE.like'),
    ('SH2B_like',r'^SH2B adapter'),('ZFYVE21_like',r'^zinc finger FYVE domain.containing protein 21\b'),
    ('SLC3_like',r'amino acid transporter heavy chain'),('kinectin_like',r'^kinectin'),
    ('LDLR_related_like',r'low.density lipoprotein receptor.related'),('nectin_PVR_like',r'^poliovirus receptor\b|^nectin'),
    ('syndecan_like',r'^syndecan'),('neuropilin_like',r'^neuropilin'),('tetraspanin_like',r'^tetraspanin|^CD151 antigen'),
    ('PDE4_like',r'phosphodiesterase 4[A-D]\b'),('polycystin_like',r'^polycystin'),('TRPM_like',r'cation channel subfamily M member'),
    ('calreticulin_like',r'^calreticulin'),('ASAP_AGAP_like',r'^Arf.GAP with (?:SH3|GTPase)'),
    ('GIT_like',r'^ARF GTPase.activating protein GIT'),('RhoGAP_like',r'^Rho GTPase.activating|^StAR.related lipid transfer protein 13\b'),
    ('RasGAP_like',r'^Ras GTPase.activating'),('defensin_like',r'^defensin'),
    ('DOCK_like',r'^dedicator of cytokinesis'),('ELMO_like',r'^engulfment and cell motility'),
    ('RhoGEF_like',r'^Rho guanine nucleotide exchange|^triple functional domain|^proto.oncogene vav|^guanine nucleotide exchange factor VAV'),
    ('dynamin_like',r'^dynamin'),('Rho_GTPase_like',r'^Rho.related GTP.binding|^transforming protein Rho'),
    ('PLCG_like',r'phosphodiesterase gamma'),('INPP5_SHIP_like',r'inositol.*5.phosphatase (?:1|2)\b'),
    ('ITGB3BP_like',r'^centromere protein R$'),('RAVER_like',r'^ribonucleoprotein PTB.binding'),
    ('STAT_like',r'^signal transducer and activator of transcription'),('SPTLC_like',r'^serine palmitoyltransferase'),
    ('PAK_like',r'^serine/threonine.protein kinase PAK'),('PDPK1_like',r'phosphoinositide.dependent protein kinase 1'),
    ('PKC_like',r'^protein kinase C '),('PRKD_like',r'^serine/threonine.protein kinase D'),
    ('PPM_like',r'^protein phosphatase 1[EFM]\b'),('PPP2_like',r'^serine/threonine.protein phosphatase 2A catalytic'),
    ('PEAK_like',r'^inactive tyrosine.protein kinase PEAK'),('SHP_PTPN6_11_like',r'phosphatase non.receptor type (?:6|11)\b'),
    ('PTPN2_like',r'phosphatase non.receptor type 2\b'),('prion_like',r'^major prion protein'),
    ('myosin_like',r'^myosin'),('MARCKS_like',r'^myristoylated alanine.rich C.kinase substrate'),
    ('ABI_like',r'^Abl interactor|^ABI gene family'),('PPP6_regulatory_like',r'phosphatase 6 regulatory ankyrin'),
    ('CSRP_like',r'^cysteine and glycine.rich protein'),('IRS_like',r'^insulin receptor substrate'),
    ('MAPK8IP_like',r'^C.Jun.amino.terminal kinase.interacting'),('myomesin_like',r'^myomesin'),
    ('tubulin_like',r'^tubulin'),('intermediate_filament_like',r'^vimentin'),('SHARPIN_like',r'^sharpin'),
    ('FABP_like',r'^fatty acid.binding protein'),('nischarin_like',r'^nischarin'),('CIB_like',r'^calcium and integrin.binding'),
    ('SRCIN1_like',r'^SRC kinase signaling inhibitor'),('ADAM_like',r'^disintegrin and metalloproteinase'),
    ('CEACAM_like',r'cell adhesion molecule CEACAM'),('CD47_like',r'^leukocyte surface antigen CD47'),('layilin_like',r'^layilin'),
    ('SIRP_like',r'^tyrosine.protein phosphatase non.receptor type substrate'),('THY1_like',r'^Thy.1 '),
    ('PLAUR_like',r'^urokinase plasminogen activator surface receptor'),('KCNH_like',r'potassium channel KCNH'),
    ('SLC16_like',r'^monocarboxylate transporter'),('SLC9_like',r'^sodium/hydrogen exchanger'),
    ('small_HSP_like',r'^heat shock protein beta'),('HSP70_like',r'^heat shock.related 70'),
    ('CBL_like',r'^E3 ubiquitin.protein ligase CBL'),('RNF5_185_like',r'^E3 ubiquitin.protein ligase RNF(?:5|185)\b'),
    ('BCAR3_like',r'breast cancer anti.estrogen resistance protein 3\b'),('RAPGEF_like',r'^Rap guanine nucleotide exchange'),
    ('SOS_like',r'^son of sevenless'),('TIAM_like',r'^Rho guanine nucleotide exchange factor TIAM'),
    ('cytohesin_like',r'^cytohesin'),('ARF_GTPase_like',r'^ADP.ribosylation factor \d'),('Ras_GTPase_like',r'^GTPase [HKN]Ras'),
    ('PLD_like',r'^phospholipase D'),('calpain_like',r'^calpain'),('caspase_like',r'^caspase'),('MMP_like',r'^matrix metalloproteinase'),
    ('PI3K_like',r'phosphatidylinositol.*3.kinase catalytic'),('PIP5K_like',r'phosphatidylinositol 4.phosphate 5.kinase'),
    ('PTEN_like',r'protein phosphatase PTEN'),('PABP_like',r'^polyadenylate.binding'),('slingshot_like',r'^protein phosphatase Slingshot'),
    ('AKT_like',r'^RAC.alpha serine/threonine.protein kinase'),('MAPK_like',r'^mitogen.activated protein kinase'),
    ('LIMK_like',r'^LIM domain kinase'),('PKA_like',r'^cAMP.dependent protein kinase catalytic'),
    ('ROCK_like',r'^Rho.associated protein kinase'),('INSR_like',r'^insulin receptor$'),('TESK_like',r'^dual specificity testis.specific'),
    ('PTPRO_like',r'^receptor.type tyrosine.protein phosphatase O\b'),('PTPRH_like',r'^receptor.type tyrosine.protein phosphatase H\b')]


def named_family(record):
    name=str(record.get('protein_name') or '')
    # TIAM is narrower than the generic RhoGEF description.
    if re.search(r'^Rho guanine nucleotide exchange factor TIAM',name,re.I):return 'TIAM_like'
    return next((family for family,pattern in NAME_FAMILIES if re.search(pattern,name,re.I)),None)


# Exact symbols are own-identity aliases in the supplied trusted label table.
# They are not mined from a recruitment crosswalk, interacting proteins, Pfam,
# CDD, the proposed benchmark panel or classifier predictions.
GENE_FAMILIES={
    'GAB1':{'GAB_like'},'GAB2':{'GAB_like'},'GAB3':{'GAB_like'},
    'ITGB1BP1':{'ICAP1_like'},'ICAP1':{'ICAP1_like'},
    'ARPC2':{'ARPC2_like'},'ARC34':{'ARPC2_like'},
    'ACTR2':{'arp2_like'},'ACTR3':{'arp3_like'},
    'PPM1F':{'PPM_like'},'POPX2':{'PPM_like'},
    'FBLIM1':{'migfilin_like'},'FBLP1':{'migfilin_like'},
    'MACF1':{'MACF1_spectraplakin_like'},'ACF7':{'MACF1_spectraplakin_like'},
    'TES':{'Testin_like'},
    **{'PRICKLE'+str(n):{'Prickle_like'} for n in range(1,5)},
    **{gene:{'actin_conventional'} for gene in ('ACTA1','ACTA2','ACTB','ACTG1','ACTG2','ACTC1')},
    **{gene:{'talin'} for gene in ('TLN1','TLN2')},
    **{gene:{'kindlin'} for gene in ('FERMT1','FERMT2','FERMT3')},
    **{gene:{'alpha_actinin'} for gene in ('ACTN1','ACTN2','ACTN3','ACTN4')},
    **{gene:{'filamin'} for gene in ('FLNA','FLNB','FLNC')},
    'VCL':{'vinculin'},'PXN':{'paxillin_like'},'LPXN':{'paxillin_like'},
    'GRB2':{'Grb2'},'ILK':{'ILK'},
}


def trusted_genes(record):
    return set(re.split(r'[;,\s]+',str(record.get('genes') or '').upper()))-{''}


def gene_families(record):
    return set().union(*(GENE_FAMILIES.get(gene,set()) for gene in trusted_genes(record)))
