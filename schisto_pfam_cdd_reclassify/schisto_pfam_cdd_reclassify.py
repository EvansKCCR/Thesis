#!/usr/bin/env python3
"""Independent post hoc refinement of unresolved Stage 2 family annotations.

Reads saved Pfam/CDD evidence; no external searches or audit reruns are launched.
"""
from __future__ import annotations
import argparse
import collections
import csv
import hashlib
import json
import re
import statistics
import sys
import zipfile
from pathlib import Path
import xml.etree.ElementTree as ET

VERSION = "1.0"
RULESET_VERSION = "2026-11-04-frozen-r1"
MAIN_NS = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
OFFICE_REL_NS = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
PKG_REL_NS = "http://schemas.openxmlformats.org/package/2006/relationships"


def _excel_col_index(ref: str) -> int:
    match = re.match(r"([A-Z]+)", ref)
    if not match:
        raise ValueError(f"Invalid Excel cell reference: {ref!r}")
    index = 0
    for character in match.group(1):
        index = index * 26 + (ord(character) - 64)
    return index - 1


def load_xlsx(path: Path, sheet_name: str) -> tuple[list[str], list[list[str]]]:
    """Read values from an Excel worksheet using only standard-library XLSX XML.

    Intentionally reads saved cell values, not computed formulas. This input's
    annotation fields are ordinary cells. No modifications are made to XLSX.
    """
    ns = "{" + MAIN_NS + "}"
    with zipfile.ZipFile(path) as archive:
        workbook = ET.fromstring(archive.read("xl/workbook.xml"))
        relationships = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
        mapping = {node.get("Id"): node.get("Target")
                   for node in relationships.findall(f"{{{PKG_REL_NS}}}Relationship")}
        available = {}
        for sh in workbook.find(ns + "sheets"):
            rel_id = sh.get(f"{{{OFFICE_REL_NS}}}id")
            target = mapping[rel_id].lstrip("/")
            if target.startswith("xl/"):
                sheet_path = target
            else:
                sheet_path = "xl/" + target
            available[sh.get("name")] = sheet_path
        if sheet_name not in available:
            raise ValueError(f"Worksheet {sheet_name!r} not found. Available: {', '.join(available)}")
        strings = []
        if "xl/sharedStrings.xml" in archive.namelist():
            with archive.open("xl/sharedStrings.xml") as stream:
                for _, element in ET.iterparse(stream, events=("end",)):
                    if element.tag == ns + "si":
                        strings.append("".join((t.text or "") for t in element.iter(ns + "t")))
                        element.clear()
        rows = []
        with archive.open(available[sheet_name]) as stream:
            for _, element in ET.iterparse(stream, events=("end",)):
                if element.tag != ns + "row":
                    continue
                cells = {}
                for cell in element.findall(ns + "c"):
                    ci = _excel_col_index(cell.get("r", ""))
                    value_node = cell.find(ns + "v")
                    datatype = cell.get("t", "")
                    if datatype == "inlineStr":
                        value = "".join((t.text or "") for t in cell.iter(ns + "t"))
                    else:
                        value = (value_node.text or "") if value_node is not None else ""
                        if datatype == "s" and value:
                            value = strings[int(value)]
                    cells[ci] = value
                rows.append(cells)
                element.clear()
        if not rows:
            raise ValueError("Worksheet contains no rows")
        if not rows[0]:
            raise ValueError("Worksheet has an empty header row")
        columns = [rows[0].get(i, "") for i in range(max(rows[0]) + 1)]
        records = [[row.get(i, "") for i in range(len(columns))] for row in rows[1:]]
        return columns, records


def write_tsv(path: Path, rows: list[dict], fields: list[str] | None = None) -> None:
    if fields is None:
        fields = list(rows[0]) if rows else []
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


from stage2_io import load_table, compare_prior, resolve_audit_input, key, token, integer_length

def contains_any(keys,terms):return any(k in terms for k in keys)
def part(keys,terms):return any(any(k.lower() in t.lower() for t in terms) for k in keys)
# Note: each classification is a new hypothesis, never a reassignment of the original audit.
CAT={'Talin-like':'Membrane–actin clutch','Parvin-like':'Membrane–actin clutch'}
def specific(p,cc,og):
 f=lambda *xs: all(x in p for x in xs)
 a=lambda *xs: any(x in cc for x in xs)
 b=lambda *xs: any(x in p for x in xs)
 # distinct diagnostic models and multidomain agreement
 if a('FERM_C_Talin','FERM_F1_TLN','FERM_F0_TLN','talin-RS') and b('PF00373','PF01608'): return ('Talin-like','A: concordant diagnostic domains','talin-specific FERM CDD and Pfam FERM/LWEQ','Membrane–actin clutch')
 if a('CH_PARV_rpt1','CH_PARV_rpt2') and b('PF00307'):return ('Parvin-like','A: concordant diagnostic domains','parvin-specific CH repeat model and PF00307','Membrane–actin clutch')
 if a('LIM1_PINCH','LIM2_PINCH','LIM3_PINCH','LIM4_PINCH','LIM5_PINCH') and b('PF00412'):return ('PINCH-like','A: concordant diagnostic domains','PINCH-specific LIM CDD and PF00412','Membrane–actin clutch')
 if b('PF18124') and b('PF00373'):return ('Kindlin-like','A: concordant diagnostic domains','kindlin-related Pfam FERM/PH architecture','Membrane–actin clutch')
 if a('FERM_F0_kindlins') and b('PF00373'):return ('Kindlin-like','A: concordant diagnostic domains','kindlin-specific FERM CDD and Pfam','Membrane–actin clutch')
 if a('CH_ACTN_rpt1','CH_ACTN_rpt2') and b('PF00307'):return ('α-actinin-like','A: concordant diagnostic domains','ACTN-specific CH CDD and Pfam CH','Membrane–actin clutch')
 if a('CH_dFLNA-like_rpt1','CH_FLN_rpt2','IG_FLMN') and b('PF00307','PF00630'):return ('Filamin-like','A: concordant diagnostic domains','filamin-specific CDD and Pfam domains','Membrane–actin clutch')
 if b('PF01044'):return ('Vinculin-like','B: family-characteristic Pfam','vinculin Pfam family domain PF01044','Membrane–actin clutch')
 if a('FERM_C_ERM') and b('PF00373'):return ('Ezrin/radixin/moesin-like FERM protein','A: concordant diagnostic domains','ERM-specific FERM CDD and PF00373','Membrane–actin clutch')
 if a('Integrin_beta','Integrin_b_cyt') and b('PF00362','PF08725'):return ('Integrin β-like','A: concordant diagnostic domains','integrin-β-specific Pfam and CDD','Membrane and receptors')
 if a('Int_alpha','FG-GAP') and b('PF01839','PF20806'):return ('Integrin α-like','A: concordant diagnostic domains','integrin-α-specific CDD and Pfam domains','Membrane and receptors')
 if b('PF20806'):return ('Integrin α-like (Ig3 segment)','B: family-characteristic Pfam','PF20806 integrin-α Ig-like region; whole receptor architecture not recovered','Membrane and receptors')
 if b('PF00335'):return ('Tetraspanin-like','B: family-characteristic Pfam','tetraspanin Pfam PF00335','Membrane and receptors')
 if a('CH_PLS_rpt1','CH_PLS_FIM_rpt2','CH_PLS_rpt3') and pf_count(p,'PF00307')>=2:return ('Plastin/fimbrin-like','A: concordant diagnostic domains','plastin/fimbrin-specific CDD CH repeats','Actin and regulator')
 if a('PTKc_Syk_like','SH2_N-SH2_Zap70_Syk_like') and f('PF00017','PF00069'):return ('Syk/ZAP70-like kinase','A: concordant diagnostic domains','Syk-specific PTKc/SH2 models and kinase Pfam','Kinase and phosphatase')
 if a('PTKc_Abl','SH2_ABL','SH3_Abl') and f('PF00017','PF00069'):return ('Abl-family kinase','A: concordant diagnostic domains','ABL-specific SH2/kinase CDD and Pfam','Kinase and phosphatase')
 if a('PTKc_Fes_like','SH2_Fps_family') and f('PF00017','PF00069'):return ('Fes/Fps-family kinase','A: concordant diagnostic domains','Fes/Fps-specific CDD and SH2 kinase architecture','Kinase and phosphatase')
 if a('PTKc_Src_like','SH2_Src_family','SH3_Src_like') and f('PF00017','PF00069'):return ('Src-family kinase','A: concordant diagnostic domains','Src-specific CDD and SH2/kinase Pfam','Kinase and phosphatase')
 if a('STKc_PAK','CRIB_PAK_like') and f('PF00069','PF00786'):return ('PAK-family kinase','A: concordant diagnostic domains','PAK-specific CDD and p21-binding/kinase Pfam','Kinase and phosphatase')
 if f('PF00069','PF00786') and a('PBD'):return ('PAK-family kinase','B: family-associated architecture','p21-binding domain and catalytic kinase domain','Kinase and phosphatase')
 if a('STKc_PKA','STKc_PKA_like') and b('PF00069'):return ('PKA-like kinase','A: concordant diagnostic domains','PKA-specific kinase CDD and PF00069','Kinase and phosphatase')
 if a('STKc_ERK1_2_like','STKc_JNK') and b('PF00069'):return ('MAPK-like kinase','A: concordant diagnostic domains','MAPK-specific CDD and catalytic kinase Pfam','Kinase and phosphatase')
 if a('STKc_PDK1') and b('PF00069'):return ('PDPK1-like kinase','A: concordant diagnostic domains','PDK1-specific CDD and catalytic Pfam','Kinase and phosphatase')
 if a('STKc_PKB') and b('PF00069'):return ('AKT/PKB-like kinase','A: concordant diagnostic domains','PKB-specific CDD and catalytic Pfam','Kinase and phosphatase')
 if a('STKc_PKD','C1_PKD_rpt1') and b('PF00069'):return ('PKD-family kinase','A: concordant diagnostic domains','PKD-specific CDD and kinase Pfam','Kinase and phosphatase')
 if a('STKc_cPKC','STKc_PKC','C1_cPKC_rpt2') and b('PF00069'):return ('PKC-family kinase','A: concordant diagnostic domains','PKC-specific CDD and kinase Pfam','Kinase and phosphatase')
 if a('STKc_GSK3') and b('PF00069'):return ('GSK3-family kinase','A: concordant diagnostic domains','GSK3-specific CDD and kinase Pfam','Kinase and phosphatase')
 if a('C1_aPKC') and f('PF00069','PF00130'):return ('Atypical PKC-like kinase','B: family-associated architecture','aPKC-type C1 and kinase domain','Kinase and phosphatase')
 if a('PTPc_motif','PTP_PTEN') and b('PF10409'):return ('PTEN-like phosphatase','A: concordant diagnostic domains','PTEN-specific CDD and phosphatase Pfam','Kinase and phosphatase')
 if a('CH_PARV_rpt1','CH_PARV_rpt2') and b('PF00307'):return ('Parvin-like','A: concordant diagnostic domains','Parvin-specific CH repeat','Membrane–actin clutch')
 if a('Cofilin_ADF') and b('PF00241'):
  return ('Twinfilin-like (two ADF-H domains)' if p.count('PF00241')>=2 else 'ADF/cofilin-domain protein', 'B: family-associated architecture','ADF-H Pfam and cofilin/ADF CDD; repeat count resolves subtype','Actin and regulator')
 if a('PTZ00335','alpha_tubulin') and b('PF00091'):return ('α-tubulin-like','A: concordant diagnostic domains','alpha-tubulin CDD and tubulin Pfam','Actin and regulator')
 if a('beta_tubulin') and b('PF00091','PF03953'):return ('β-tubulin-like','A: concordant diagnostic domains','beta-tubulin CDD and tubulin Pfam','Actin and regulator')
 if a('gamma_tubulin') and b('PF00091','PF03953'):return ('γ-tubulin-like','A: concordant diagnostic domains','gamma-tubulin CDD and tubulin Pfam','Actin and regulator')
 if a('delta_zeta_tubulin-like') and b('PF00091'):return ('δ/ζ-tubulin-related','A: concordant diagnostic domains','delta/zeta-specific CDD and tubulin Pfam','Actin and regulator')
 if a('Calreticulin') and b('PF00262'):return ('Calreticulin-like','A: concordant diagnostic domains','calreticulin CDD and Pfam','Actin and regulator')
 if a('ASKHA_NBD_HSP70_HSPA1','HSP70','DnaK') and b('PF00012'):return ('HSP70-family chaperone','A: concordant diagnostic domains','HSP70/DnaK CDD and Pfam','Actin and regulator')
 if a('PIP5K','PIPKc_PIP5KI') and b('PF01504'):return ('PIP5K-like','A: concordant diagnostic domains','PIP5K CDD and Pfam','Kinase and phosphatase')
 if a('P34-Arc') and b('PF04045'):return ('ARPC2-like','A: concordant diagnostic domains','ARPC2 CDD/Pfam','Actin and regulator')
 if a('Myosin_head','MYSc') and b('PF00063'):return ('Myosin-like motor protein','A: concordant diagnostic domains','myosin-specific head CDD and Pfam','Actin and regulator')
 if a('Dynamin_N','DYNc','DLP_1') and b('PF00169'):return ('Dynamin-family GTPase','B: family-associated architecture','dynamin-specific CDD and PH domain','Adaptors and small GTPases')
 if a('RasGEF','RasGEF_N') and b('PF00617','PF00618'):return ('RasGEF-like','A: concordant diagnostic domains','RasGEF CDD and Pfam','Adaptors and small GTPases')
 if a('RhoGEF') and b('PF00621'):return ('RhoGEF-like','A: concordant diagnostic domains','RhoGEF CDD and Pfam','Adaptors and small GTPases')
 if a('ArfGap','ArfGap_AGAP') and b('PF01412'):return ('ARF-GAP-like','A: concordant diagnostic domains','ArfGAP CDD and Pfam','Adaptors and small GTPases')
 if a('Sec7') and b('PF01369'):return ('Cytohesin/Sec7-family ARF-GEF','A: concordant diagnostic domains','Sec7 CDD/Pfam','Adaptors and small GTPases')
 if a('RhoA_like','RHO') and b('PF00071'):return ('Rho-family small GTPase','A: concordant diagnostic domains','Rho-specific CDD and Ras Pfam','Adaptors and small GTPases')
 if a('RAB','Rab11_like','Rab') and b('PF00071'):return ('Rab-family small GTPase','A: concordant diagnostic domains','Rab CDD and Ras Pfam','Adaptors and small GTPases')
 if a('Arf1_5_like','Arf') and b('PF00025','PF00071'):return ('ARF-family small GTPase','A: concordant diagnostic domains','Arf-specific CDD and GTPase Pfam','Adaptors and small GTPases')
 if a('H_N_K_Ras_like','RAS') and b('PF00071'):return ('Ras-family small GTPase','A: concordant diagnostic domains','RAS-specific CDD and GTPase Pfam','Adaptors and small GTPases')
 if a('LIM2_Ajuba_like') and b('PF00412'):return ('Ajuba/zyxin-like LIM adaptor','A: concordant diagnostic domains','Ajuba-specific LIM CDD and Pfam','Adaptors and small GTPases')
 if a('LIM4_Leupaxin') and b('PF00412'):return ('Paxillin/leupaxin-like','A: concordant diagnostic domains','Leupaxin-specific LIM CDD and Pfam','Adaptors and small GTPases')
 if a('SH2_Grb2_like') and f('PF00017','PF00018'):return ('GRB2-like adaptor','A: concordant diagnostic domains','Grb2-specific SH2 and SH3 Pfam','Adaptors and small GTPases')
 if a('SH3_Nck_1','SH3_Nck_2') and f('PF00017','PF00018'):return ('NCK-like adaptor','A: concordant diagnostic domains','NCK-specific SH3 and SH2/SH3 Pfam','Adaptors and small GTPases')
 if a('PTB_Shc') and b('PF00640'):return ('SHC-like adaptor','A: concordant diagnostic domains','SHC-specific PTB CDD and Pfam','Adaptors and small GTPases')
 if a('STKc_PAK','PBD') and f('PF00069','PF00786'):return ('PAK-family kinase','A: concordant diagnostic domains','PAK p21 binding + kinase','Kinase and phosphatase')
 if a('ZnMc_MMP') and b('PF00413'):return ('MMP-family metalloprotease','A: concordant diagnostic domains','MMP CDD and Pfam','Extracellular')
 if a('ZnMc_adamalysin_II_like','Disintegrin') and b('PF00200','PF01421'):return ('ADAM-family metalloprotease','A: concordant diagnostic domains','ADAM CDD and Pfam','Membrane and receptors')
 if a('Collagen') and b('PF01391'):return ('Collagen-like','A: concordant diagnostic domains','collagen CDD and PF01391','Extracellular')
 if a('Laminin_EGF','Laminin_N') and b('PF00055','PF00053'):return ('Laminin-like','A: concordant diagnostic domains','laminin-specific N/EGF CDD and Pfam','Extracellular')
 if a('STKc_GSK3') and b('PF00069'):return ('GSK3-like','A: concordant diagnostic domains','GSK3 CDD','Kinase and phosphatase')
 if a('PTPc') and b('PF00102') and b('PF00041','PF00047','PF07679','PF13927'):return ('Receptor-like tyrosine phosphatase (FN3/Ig–PTP)','B: family-associated architecture','PTP catalytic domain with extracellular Ig/FN3 domains','Kinase and phosphatase')
 if a('PTPc') and b('PF00102'):return ('Protein tyrosine phosphatase (PTP-domain)','B: family-associated architecture','PTP CDD and PF00102','Kinase and phosphatase')
 if a('PP2Ac') and b('PF00149'):return ('PPP/PP2A-like phosphatase','B: family-associated architecture','PP2A catalytic CDD and PF00149','Kinase and phosphatase')
 if a('PP2Cc','PP2C') and b('PF00481'):return ('PP2C/PPM-type phosphatase','A: concordant diagnostic domains','PP2C CDD and Pfam','Kinase and phosphatase')
 if a('PDEase_I','HDc') and b('PF00233'):return ('Cyclic nucleotide phosphodiesterase (PDE)','A: concordant diagnostic domains','PDE catalytic CDD and PF00233','Kinase and phosphatase')
 if a('Peptidase_C14','CASc') and b('PF00656'):return ('Caspase-family cysteine protease','A: concordant diagnostic domains','caspase CDD and Pfam','Actin and regulator')
 if a('Na_H_Exchanger','NhaP') and b('PF00999'):return ('Na+/H+ exchanger-family transporter','A: concordant diagnostic domains','Na+/H+ exchanger CDD and Pfam','Membrane and receptors')
 if a('Ion_trans') and f('PF00027','PF00520') and b('PF13426'):return ('KCNH-like channel','B: family-associated architecture','ion channel + cNMP + PAS signatures','Membrane and receptors')
 if a('Ion_trans') and f('PF00027','PF00520'):return ('Cyclic nucleotide-binding ion channel-like','B: family-associated architecture','ion-transport and cyclic nucleotide-binding domains','Membrane and receptors')
 if a('FERM_C_FARP1-like') and b('PF00373'):return ('FARP-like FERM–PH adaptor','A: concordant diagnostic domains','FARP-specific FERM CDD','Adaptors and small GTPases')
 if b('PF00046') and b('PF00412') and a('HOX','homeodomain'):return ('LIM-homeobox transcription factor','A: concordant diagnostic domains','homeodomain and LIM Pfam + HOX CDD','Adaptors and small GTPases')
 if a('PDZ_ZASP52-like') and b('PF00412','PF00595'):return ('PDZ–LIM/ZASP-like protein','A: concordant diagnostic domains','PDZ–LIM CDD and Pfam','Actin and regulator')
 if b('PF00307','PF00435') and a('Spectrin','SPEC'):return ('Spectrin-family cytoskeletal protein','B: family-associated architecture','spectrin repeats and calponin-homology Pfam/CDD','Membrane–actin clutch')
 if a('BTB') and b('PF00651') and b('PF01344'):return ('BTB–Kelch adaptor','B: family-associated architecture','BTB/BACK/Kelch multi-domain architecture','Adaptors and small GTPases')
 if a('RRM1_I_PABPs','PABP') and b('PF00076'):return ('Poly(A)-binding protein-like','A: concordant diagnostic domains','PABP-specific RNA-binding CDD','Actin and regulator')
 return None

def pf_count(p,token):return p.count(token)

def broad(p,cc):
 b=lambda *s: any(x in p for x in s)
 a=lambda *s: any(x in cc for x in s)
 if b('PF00069','PF07714') and a('S_TKc','PTKc','PK_Tyr_Ser-Thr','SPS1','Pkinase','STYKc','TyrKc'):
  if b('PF00017') and b('PF00018'):return ('SH2/SH3 tyrosine-kinase domain protein','C: catalytic/domain class','SH2–SH3 and tyrosine-kinase domains','Kinase and phosphatase')
  if b('PF00017'):return ('SH2-containing tyrosine kinase','C: catalytic/domain class','SH2 domain and kinase core','Kinase and phosphatase')
  if b('PF00130'):return ('C1-domain-containing kinase','C: catalytic/domain class','C1 and kinase domains; PKC/PKD subtype open','Kinase and phosphatase')
  return ('Protein kinase superfamily (subtype unresolved)','C: catalytic/domain class','PF00069/PF07714; general kinase CDD','Kinase and phosphatase')
 if b('PF00071','PF00025') and a('Ras','Gem1','small_GTP','Rab','Arf','RHO','RAS'):
  return ('Ras-superfamily GTPase (subtype unresolved)','C: catalytic/domain class','Ras-related Pfam and broad CDD GTPase models','Adaptors and small GTPases')
 if b('PF00400') and a('WD40'):return ('WD40-repeat protein (family unresolved)','C: domain superfamily','WD40 repeats in Pfam/CDD','Adaptors and small GTPases')
 if b('PF00307') and a('CH','Calponin','CH_SF'):
  return ('Calponin-homology (CH)-domain protein','C: domain superfamily','CH Pfam and CDD; multiple possible actin-binding lineages','Actin and regulator')
 if b('PF00373') and a('B41','FERM_M','FERM_C'):
  return ('FERM-domain protein (subtype unresolved)','C: domain superfamily','FERM Pfam and CDD','Membrane–actin clutch')
 if b('PF00412') and a('LIM'):
  return ('LIM-domain protein (subtype unresolved)','C: domain superfamily','LIM domain Pfam/CDD','Adaptors and small GTPases')
 if b('PF12796','PF00023') and a('Ank_2','ANKYR','ANK'):return ('Ankyrin-repeat protein','C: domain superfamily','ANK repeat Pfam and CDD','Adaptors and small GTPases')
 if b('PF00169') and a('PH'):return ('PH-domain protein (subtype unresolved)','C: domain superfamily','PH domain Pfam and CDD','Adaptors and small GTPases')
 if b('PF00036') and a('EFh','EF-hand_7','FRQ1'):return ('EF-hand calcium-binding protein','C: domain superfamily','EF-hand Pfam/CDD','Actin and regulator')
 if b('PF01344','PF07707') and a('Kelch','BACK'):return ('Kelch/BACK-repeat protein','C: domain superfamily','Kelch/BACK CDD and Pfam','Adaptors and small GTPases')
 if b('PF00041') and a('FN3','fn3'):return ('Fibronectin type III-domain protein','C: domain superfamily','FN3 Pfam and CDD (not fibronectin protein identity)','Membrane and receptors')
 if b('PF00149') and a('Metallophos'):return ('Metallophosphoesterase/phosphatase','C: catalytic/domain class','metallophosphoesterase Pfam/CDD','Kinase and phosphatase')
 if b('PF00091','PF03953') and a('Tubulin'):return ('Tubulin-family protein (isotype unresolved)','C: domain superfamily','tubulin Pfam and CDD','Actin and regulator')
 if b('PF00027','PF00520') and a('cNMP','Ion_trans'):return ('Cyclic nucleotide-binding channel/regulator','C: domain superfamily','cNMP/ion transport domain','Membrane and receptors')
 if b('PF00054','PF02210') and a('LamG','Laminin_G_2'):return ('Laminin-G-domain protein (ligand unresolved)','C: domain superfamily','laminin G Pfam/CDD','Extracellular')
 if b('PF00076') and a('RRM'):return ('RNA recognition motif protein','C: domain superfamily','RRM Pfam/CDD','Unassigned/other')
 if b('PF00046') and a('homeodomain','HOX'):return ('Homeodomain transcription factor','C: domain superfamily','homeodomain Pfam/CDD','Unassigned/other')
 if b('PF00038'):return ('Intermediate-filament-related protein','C: domain superfamily','PF00038','Actin and regulator')
 if b('PF00651') and a('BTB'):return ('BTB-domain protein','C: domain superfamily','BTB Pfam/CDD','Adaptors and small GTPases')
 if b('PF00041','PF00047','PF07679') and a('IG_like','Ig','I-set'):return ('Immunoglobulin-fold protein','C: domain superfamily','Ig domains Pfam/CDD','Membrane and receptors')
 if b('PF00018') and a('SH3','SH3_1'):return ('SH3-containing adaptor protein','C: domain superfamily','SH3 Pfam/CDD','Adaptors and small GTPases')
 if b('PF00017') and a('SH2'):return ('SH2-containing signalling protein','C: domain superfamily','SH2 Pfam/CDD','Adaptors and small GTPases')
 if b('PF00595') and a('PDZ'):return ('PDZ-domain protein','C: domain superfamily','PDZ Pfam/CDD','Adaptors and small GTPases')
 if b('PF00097') and a('RING'):return ('RING-finger ubiquitin ligase-related protein','C: domain superfamily','RING Pfam/CDD','Adaptors and small GTPases')
 if b('PF07690') and a('MFS_1'):return ('Major facilitator superfamily transporter','C: domain superfamily','MFS Pfam/CDD','Membrane and receptors')
 if b('PF00620') and a('RhoGAP'):return ('RhoGAP-domain protein','C: catalytic/domain class','RhoGAP Pfam/CDD','Adaptors and small GTPases')
 if b('PF00241') and a('Cofilin_ADF'):return ('ADF-H-domain protein','C: domain superfamily','ADF/cofilin homology domain','Actin and regulator')
 if b('PF00435') and a('SPEC','Spectrin'):return ('Spectrin-repeat protein','C: domain superfamily','Spectrin Pfam and CDD','Actin and regulator')
 if b('PF00168') and a('C2'):return ('C2-domain-containing protein','C: domain superfamily','C2 Pfam/CDD','Membrane and receptors')
 if b('PF00155') and a('BioF','Aminotran_1_2','KBL_like'):return ('Aminotransferase / serine palmitoyltransferase-related','C: catalytic/domain class','PLP-dependent aminotransferase Pfam/CDD','Unassigned/other')
 if b('PF00262'):return ('Calreticulin-family protein','C: domain superfamily','calreticulin Pfam','Actin and regulator')
 if b('PF01576'):return ('Myosin N-terminal motor-associated domain protein','C: domain superfamily','Pfam myosin-associated domain','Actin and regulator')
 if b('PF01608') and a('I_LWEQ'):return ('I/LWEQ actin-binding tail protein','C: domain superfamily','I/LWEQ Pfam/CDD','Membrane–actin clutch')
 if b('PF00307'):return ('CH-domain protein (family unconfirmed)','C: domain superfamily','calponin-homology Pfam','Actin and regulator')
 if b('PF00069','PF07714'):return ('Protein kinase superfamily (subtype unresolved)','C: catalytic/domain class','kinase Pfam, no discriminating CDD','Kinase and phosphatase')
 if b('PF00025','PF00071'):return ('Ras-related GTPase-domain protein','C: catalytic/domain class','GTPase Pfam, subtype unresolved','Adaptors and small GTPases')
 if b('PF00412'):return ('LIM-domain protein (family unresolved)','C: domain superfamily','LIM Pfam','Adaptors and small GTPases')
 if b('PF00400'):return ('WD40-repeat protein (family unresolved)','C: domain superfamily','WD40 Pfam','Adaptors and small GTPases')
 if b('PF00041'):return ('FN3 domain-containing protein','C: domain superfamily','FN3 Pfam, receptor identity unresolved','Membrane and receptors')
 if b('PF00373'):return ('FERM-domain protein','C: domain superfamily','FERM Pfam','Membrane–actin clutch')
 if b('PF00169'):return ('PH-domain protein','C: domain superfamily','PH Pfam','Adaptors and small GTPases')
 return ('Insufficient discriminating Pfam/CDD evidence','D: unclassified','No concordant or distinctive family-defining domain signature','Unassigned/other')


REVIEW_FIELDS = ['Species','Protein ID','Orthogroup','Original support','Original audited family',
    'Suggested family / class','Evidence tier','Functional module (proposed)','Priority','Pfam IDs',
    'CDD models (selected)','Pfam architecture (full)','CDD architecture (full)','Evidence rationale',
    'Supported-family orthogroup context','Stage 1 hypothesis','DeepLoc 2.1','DeepTMHMM',
    'SignalP 6.0','Protein length (aa)']
POSTHOC_FIELDS = ['posthoc_status','posthoc_proposed_family_or_class','posthoc_evidence_tier',
    'posthoc_functional_module','posthoc_priority','posthoc_evidence_rationale',
    'posthoc_supported_orthogroup_context','posthoc_ruleset_version',
    'posthoc_supported_crosscheck_status','posthoc_supported_domain_family',
    'posthoc_supported_domain_tier','posthoc_supported_domain_evidence']


def pf(row):
    return re.findall(r'\b(PF\d{5})(?:\.\d+)?\b', row['Pfam_architecture'])


def cdd(row):
    return [re.sub(r':\s*\d+\s*-\s*\d+\s*$', '', x).strip()
            for x in row['NCBI_CDD_architecture'].split(';') if x.strip()]


def orthogroup(row):
    value = row['Orthology_Orthogroup'].strip()
    # A collapsed multi-group value is not a single orthogroup membership.
    return value if value and not re.search(r'[;,|\s]', value) else ''


def build_context(supported):
    groups = collections.defaultdict(list)
    for row in supported:
        if orthogroup(row) and row['putative_audited_family'].strip():
            groups[orthogroup(row)].append(row)
    result = {}
    for group, rows in groups.items():
        labels = {}
        for row in rows:
            family = row['putative_audited_family'].strip()
            labels.setdefault(family.casefold(), family)
        counts = collections.Counter(r['putative_audited_family'].strip().casefold() for r in rows)
        # Explicitly show ties; they cannot supply unanimous anchoring below.
        maximum = max(counts.values())
        family = ' / '.join(labels[k] for k in sorted(counts) if counts[k] == maximum)
        common = collections.Counter(x for row in rows for x in set(pf(row)))
        result[group] = (family, maximum, len(rows),
            {x for x,n in common.items() if n/len(rows)>=0.4},
            int(statistics.median(integer_length(r['length']) for r in rows)))
    return result


def classify(rows, context):
    output = []
    for row in rows:
        p, cc, group = pf(row), cdd(row), orthogroup(row)
        hit = specific(p, cc, group) or broad(p, cc)
        proposal, tier, rule, category = hit
        og_context = ''
        if group in context:
            family, n, total, common_pf, median = context[group]
            og_context = f'{family} ({n}/{total} supported)'
            shared = set(p) & common_pf
            if not tier.startswith('A') and n>=2 and n==total and shared and p and (
                    len(shared)/len(set(p))>=0.5 or len(set(p))==1):
                family_token = family.casefold().replace(' ','')
                proposal_token = proposal.casefold().replace(' ','')
                if tier.startswith('B') and (family_token in proposal_token or
                        proposal.startswith('Integrin α') or proposal.startswith('Talin')):
                    rule += '; anchored by supported orthogroup members'
                if tier.startswith('C') and family_token not in proposal_token:
                    proposal += f' [OG suggests {family}]'
                    rule += '; orthogroup contains supported '+family+' (not sufficient for family assignment)'
            if proposal.startswith('Insufficient') and n==total and n>=2:
                proposal = 'No diagnostic model; orthogroup context: '+family
                rule = 'Supported members of this orthogroup exist but candidate lacks discriminating domain evidence'
        og_prior = group in context and context[group][1]>=3 and context[group][1]==context[group][2]
        score = {'A':4,'B':3,'C':2,'D':1}[tier[0]]
        priority = ('P1 — targeted review' if score>=3 and (category in
                    ('Membrane–actin clutch','Membrane and receptors','Kinase and phosphatase') or og_prior)
                    else 'P2 — annotate domain family' if score>=3
                    else 'P3 — domain-class cataloguing' if score==2 else 'P4 — further profiling')
        if tier.startswith('C') and og_prior and category=='Membrane–actin clutch':
            priority = 'P1 — targeted review'
        cc_short = '; '.join(dict.fromkeys(cc))
        if len(cc_short)>280:
            cc_short = cc_short[:277]+'...'
        output.append(dict(zip(REVIEW_FIELDS, [row['Specie'],row['protein_id'],
            row['Orthology_Orthogroup'],row['family_support_level'],row['putative_audited_family'],
            proposal,tier,category,priority,'; '.join(dict.fromkeys(p)),cc_short,
            row['Pfam_architecture'],row['NCBI_CDD_architecture'],rule,og_context,
            row['stage1_family_hypothesis'],row['DeepLoc_2.1'],row['DeepTMHMM'],
            row['SignalP.v6.0'],integer_length(row['length']) if row['length'].strip() else None])))
    return output


def family_label(row):
    return row['Suggested family / class'].split(' [OG suggests')[0]


_FAMILY_GENERIC={'family','protein','domain','related','class','subfamily','superfamily','canonical','like'}
_FAMILY_SYNONYMS={'parva':'parvin','parvb':'parvin','actopaxin':'parvin',
    'fermt1':'kindlin','fermt2':'kindlin','fermt3':'kindlin',
    'actn1':'actinin','actn2':'actinin','actn3':'actinin','actn4':'actinin',
    'ezrin':'erm','radixin':'erm','moesin':'erm',
    'fps':'fesfps','fer':'fesfps','fes':'fesfps'}


def family_tokens(value):
    text=str(value or '').casefold().translate(str.maketrans({'α':' alpha ','β':' beta ','γ':' gamma ','δ':' delta ','ζ':' zeta '}))
    tokens={t for t in re.findall(r'[a-z0-9]+',text) if t not in _FAMILY_GENERIC and len(t)>1}
    return {_FAMILY_SYNONYMS.get(t,t) for t in tokens}


def supported_domain_crosschecks(supported):
    """Flag only A-tier domain evidence that names a different audit family."""
    crosschecks=[]
    for original,proposal in zip(supported,classify(supported,{})):
        audited=original.get('putative_audited_family','').strip()
        proposed=family_label(proposal)
        tier=proposal['Evidence tier']
        audited_tokens=family_tokens(audited);proposed_tokens=family_tokens(proposed)
        if not audited or audited.casefold() in {'unassigned','unresolved','unknown'}:
            status='not_assessable'
            reason='No specific audited family is available for comparison'
        elif token(original.get('family_support_level'))=='broad_supported':
            status='broad_scope_not_comparable'
            reason='The audit call is intentionally broad; a specific domain proposal alone cannot establish a false positive'
        elif tier.startswith('A') and proposed_tokens and not (audited_tokens & proposed_tokens):
            status='potential_false_positive_candidate'
            reason='A-tier Pfam/CDD evidence proposes a different family label; review the full architecture and audit rationale'
        elif tier.startswith('A') and audited_tokens & proposed_tokens:
            status='domain_concordant'
            reason='A-tier domain proposal shares a normalized family label with the audit assignment'
        else:
            status='no_high_specificity_domain_conflict'
            reason='No conflicting A-tier Pfam/CDD family proposal was produced'
        crosschecks.append(dict(Specie=original['Specie'],protein_id=original['protein_id'],
            family_support_level=original.get('family_support_level',''),
            audited_family=audited,domain_proposed_family=proposed,evidence_tier=tier,
            domain_evidence=proposal['Evidence rationale'],status=status,reason=reason,
            Pfam_architecture=original.get('Pfam_architecture',''),
            NCBI_CDD_architecture=original.get('NCBI_CDD_architecture','')))
    return crosschecks


def execute(args):
    source = resolve_audit_input(args.input,args.kind)
    initial_hash = sha256(source)
    data = load_table(source, args.sheet, load_xlsx)
    scopes = collections.Counter(r['posthoc_status'] for r in data['decisions'])
    summary = dict(version=VERSION, input=str(source), input_sha256=initial_hash,
        kind=args.kind, sheet=data['sheet'], input_rows=data['source_rows'], unique_species_protein_pairs=len(data['rows']),
        duplicate_rows=len(data['duplicates']), scope_counts=dict(scopes), column_mapping=data['column_mapping'])
    statuses = {(r['Specie'],r['protein_id']):r['posthoc_status'] for r in data['decisions']}
    supported = [r for r in data['rows'] if statuses[key(r)]=='supported_preserved']
    crosschecks = supported_domain_crosschecks(supported)
    flagged_crosschecks = [r for r in crosschecks if r['status']=='potential_false_positive_candidate']
    summary.update(supported_assignment_domain_crosschecks=len(crosschecks),
                   potential_false_positive_candidates=len(flagged_crosschecks))
    if args.check_inputs:
        summary.update(status='inputs_valid', analysis_launched=False, files_written=False)
        print(json.dumps(summary, indent=2))
        return 0
    if args.outdir is None:
        raise ValueError('--outdir is required unless --check-inputs is used')
    destination = args.outdir.resolve()
    if destination.exists():
        raise ValueError('Output directory already exists; choose a new directory: '+str(destination))
    if destination == source.parent or destination in source.parents:
        raise ValueError('Keep post hoc outputs separate from the source audit')
    if (source.parent/'workflow_manifest.json').is_file() and source.parent in destination.parents:
        raise ValueError('Keep post hoc outputs outside the sealed source audit directory')
    crosscheck_by_key={(r['Specie'],r['protein_id']):r for r in crosschecks}
    unresolved = [r for r in data['rows'] if statuses[key(r)]=='reviewed']
    output = classify(unresolved, build_context(supported))
    tiers = collections.Counter(r['Evidence tier'][0] for r in output)
    expected = {'A':609, 'B':270, 'C':2239, 'D':163}
    reproduces = len(output)==3281 and dict(tiers)==expected and len(data['rows'])==4359 and len(supported)==1078
    comparison = compare_prior(output,args.compare) if args.compare else None
    # Retain every source column and value in the all-record reconciliation table.
    if set(data['columns']) & set(POSTHOC_FIELDS):
        raise ValueError('Input already contains posthoc fields; supply the original Stage 2 audit')
    proposals = {(r['Species'].strip(),r['Protein ID'].strip()):r for r in output}
    reconciled = []
    for row in data['rows']:
        proposal = proposals.get(key(row),{})
        reconciled.append(dict(data['originals'][key(row)],
            posthoc_status=statuses[key(row)],
            posthoc_proposed_family_or_class=proposal.get('Suggested family / class',''),
            posthoc_evidence_tier=proposal.get('Evidence tier',''),
            posthoc_functional_module=proposal.get('Functional module (proposed)',''),
            posthoc_priority=proposal.get('Priority',''),
            posthoc_evidence_rationale=proposal.get('Evidence rationale',''),
            posthoc_supported_orthogroup_context=proposal.get('Supported-family orthogroup context',''),
            posthoc_ruleset_version=RULESET_VERSION,
            posthoc_supported_crosscheck_status=crosscheck_by_key.get(key(row),{}).get('status','not_applicable'),
            posthoc_supported_domain_family=crosscheck_by_key.get(key(row),{}).get('domain_proposed_family',''),
            posthoc_supported_domain_tier=crosscheck_by_key.get(key(row),{}).get('evidence_tier',''),
            posthoc_supported_domain_evidence=crosscheck_by_key.get(key(row),{}).get('domain_evidence','')))
    if sha256(source)!=initial_hash:
        raise ValueError('Source audit changed while preparing the refinement; no outputs written')
    destination.mkdir(parents=True, exist_ok=False)
    write_tsv(destination/'Schisto_unassigned_Pfam_CDD_classification.tsv',output,REVIEW_FIELDS)
    write_tsv(destination/'Stage_2_posthoc_refinement.tsv',reconciled,data['columns']+POSTHOC_FIELDS)
    write_tsv(destination/'scope_accounting.tsv',data['decisions'],
        ['Specie','protein_id','original_support','posthoc_status','reason'])
    write_tsv(destination/'supported_assignment_crosschecks.tsv',crosschecks,
        ['Specie','protein_id','family_support_level','audited_family','domain_proposed_family','evidence_tier','domain_evidence','status','reason','Pfam_architecture','NCBI_CDD_architecture'])
    write_tsv(destination/'potential_false_positive_candidates.tsv',flagged_crosschecks,
        ['Specie','protein_id','family_support_level','audited_family','domain_proposed_family','evidence_tier','domain_evidence','status','reason','Pfam_architecture','NCBI_CDD_architecture'])
    write_tsv(destination/'duplicate_protein_ids.tsv',data['duplicates'],
        ['Protein ID','kept_species','duplicate_species','conflicting_data'])
    write_tsv(destination/'tier_summary.tsv',[
        dict(Tier=t,Candidates=tiers[t],Meaning={
            'A':'Concordant diagnostic domain heuristics; inspect the stated rule',
            'B':'Family-characteristic Pfam or multidomain pattern',
            'C':'Broad domain or superfamily only',
            'D':'No discriminating domain profile'}[t]) for t in 'ABCD'],['Tier','Candidates','Meaning'])
    families = collections.defaultdict(list)
    for row in output:
        families[family_label(row)].append(row)
    write_tsv(destination/'family_summary.tsv',[
        {'Proposed family / class':family, 'Candidates':len(rows),
         **{t:sum(r['Evidence tier'].startswith(t) for r in rows) for t in 'ABCD'},
         'Species':len({r['Species'] for r in rows}),
         'Orthogroups':len({r['Orthogroup'] for r in rows if r['Orthogroup']})}
        for family,rows in sorted(families.items(),key=lambda item:(-len(item[1]),item[0]))],
        ['Proposed family / class','Candidates','A','B','C','D','Species','Orthogroups'])
    modules = collections.Counter(r['Functional module (proposed)'] for r in output)
    write_tsv(destination/'module_summary.tsv',[
        {'Functional module':m,'Candidates':n} for m,n in sorted(modules.items(),key=lambda x:(-x[1],x[0]))],
        ['Functional module','Candidates'])
    write_tsv(destination/'priority_review.tsv',[r for r in output if r['Priority'].startswith('P1')],REVIEW_FIELDS)
    if comparison is not None:
        write_tsv(destination/'comparison_differences.tsv',comparison.pop('differences'),
            ['Species','Protein ID','Column','Prior','Reproduced'])
    failed = args.strict and not reproduces or comparison is not None and not comparison['all_rows_match']
    manifest = dict(summary, status='validation_failed' if failed else 'complete',
        algorithm_version=VERSION, ruleset_version=RULESET_VERSION,
        code_sha256=sha256(Path(__file__)), adapter_sha256=sha256(Path(__file__).with_name('stage2_io.py')),
        supported_reference_records=len(supported), unassigned_review_records=len(unresolved),
        supported_assignment_domain_crosschecks=len(crosschecks),
        potential_false_positive_candidates=len(flagged_crosschecks),
        proposed_rows=len(output), tier_counts=dict(sorted(tiers.items())),
        source_record_accounting_passed=(len(reconciled)==sum(scopes.values()) and
            len(reconciled)+len(data['duplicates'])==data['source_rows']),
        source_values_preserved=True, original_family_support_was_preserved=True,
        classification_is_only_a_review_hypothesis=True,
        orthogroup_context_is_not_a_primary_diagnostic_criterion=True,
        external_searches_recomputed=False, audit_rerun=False,
        expected_from_reference_run=expected, reference_counts_match=reproduces,
        strict_baseline_requested=args.strict, comparison=comparison,
        output_sha256={p.name:sha256(p) for p in sorted(destination.iterdir()) if p.is_file()})
    (destination/'run_manifest.json').write_text(json.dumps(manifest,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print(f'Input rows: {data["source_rows"]}; unique species/protein pairs: {len(reconciled)}')
    print(f'Supported preserved: {len(supported)}; reviewed: {len(unresolved)}; excluded: {scopes["negative_excluded"]}')
    print(f'Potential false-positive audit assignments flagged for review: {len(flagged_crosschecks)}')
    print('Post hoc review tiers: '+', '.join(f'{t}={tiers[t]}' for t in 'ABCD'))
    print('Output: '+str(destination/'Stage_2_posthoc_refinement.tsv'))
    if failed:
        print('ERROR: Requested frozen-baseline/comparison check failed; inspect run_manifest.json',file=sys.stderr)
    return 1 if failed else 0


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--version',action='version',version='schisto-pfam-cdd-posthoc '+VERSION)
    parser.add_argument('--input',type=Path,required=True,help='Stage_2_output.xlsx/TSV/CSV, final-audit directory, or calibration/validation run directory')
    parser.add_argument('--kind',choices=['calibration','validation'],help='Select one branch when --input is a run directory containing both calibration and validation audits')
    parser.add_argument('--sheet',default='Audited_family_assignment',help='Worksheet for XLSX input')
    parser.add_argument('--outdir','--out',type=Path,help='New separate output directory; existing outputs are refused')
    parser.add_argument('--check-inputs',action='store_true',help='Validate schema/scope/keys only; write nothing')
    parser.add_argument('--compare',type=Path,help='Prior review TSV for field-level comparison')
    parser.add_argument('--strict',action='store_true',help='Historical 2026-11-04 reproduction counts only; omit for new audit outputs')
    args = parser.parse_args(argv)
    try:
        return execute(args)
    except (OSError,ValueError,KeyError,zipfile.BadZipFile,ET.ParseError) as exc:
        print('ERROR: '+str(exc),file=sys.stderr)
        return 1


if __name__=='__main__':
    sys.exit(main())

