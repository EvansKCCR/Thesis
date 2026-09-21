# Requested-domain recalibration (v4)

The active library contains 96 models: 92 previously selected models plus PF00039, PF00040, PF01033 and PF12796. Every requested domain is included. Versioned HMMs, source metadata and checksums for the four additions are in `databases/pfam/targeted_v4/`. This is a documented composite library, not a claim that all models came from one Pfam release.

## Family rules

| Target | Architecture and interpretation |
|---|---|
| Fibronectin-like | FN type I (PF00039), type II (PF00040) and type III (PF00041). The reference-like architecture uses all three types in order, without requiring vertebrate copy counts. Partial hits remain visible but are not automatically fibronectin orthologues. |
| Vitronectin-like | Somatomedin-B (PF01033) before hemopexin-like (PF00045), together with reciprocal family evidence. Either domain alone is generic. RGD is annotated as supporting context only. |
| Alpha integrin | PF01839, PF08441, PF20805 and PF20806 are all used; PF00357 remains an additional previously validated signature. Variable repeat counts and partial-domain rescue routes remain. The head-to-stalk order is reported, without excluding incomplete candidates. GFFxR variants remain allowed and an alpha-I domain is not mandatory. |
| FAK family | PF21477, PF00373, PF18038, PF03623 and PF07714 are all used. PF18038 now explicitly participates in architecture alternatives and recruitment, rather than only the ordering check. Generic FERM-N alone is not specific. The diagnostic ordered PF21477-to-kinase partial-core rescue remains, with missing FAT reported. PF00069 remains an alternative kinase model. |
| ILK | PF12796 followed by PF07714 (or PF00069), with reciprocal family evidence. The older PF00023-repeat architecture remains an alternative. PF12796 is an ankyrin-repeat model covering three copies; one model match is sufficient for the alternative. Hits from overlapping ankyrin models are not added together to manufacture repeat counts. A catalytic-activity motif is not required for this pseudokinase family. |

The legacy `fibronectin_III_ligand` rule is now explicitly **screening-only**. Its single PF00041 domain is not sufficient to identify fibronectin. The actual FN1 reference was explicitly relabelled `fibronectin_like` in the v4 seeds and paired catalogues; old homology results are not silently relabelled. `generic_ECM_domain_screen` preserves isolated FN1/FN2/SMB/hemopexin hits without assigning a family.

Supported assignments still require reciprocal sequence evidence and compatible architecture/independent evidence. An architecture match or short motif alone is insufficient. Topology that has not been tested remains a warning; none was inferred from an extracellular-sounding domain name. Confidence categories are qualitative evidence tiers, **not calibrated probabilities**.

These retention statements refer to the default comparative mode. Optional strict mode retains its harder architecture, length and topology requirements; it should not be used to infer gene-family absence from divergent or incomplete sequences.

## Reference controls and limitations

The new paired set has **55 forward** and **126 reciprocal** records, including the existing 20 computational spiralian anchors. Human vitronectin (P04004) was added as an annotated supplemental forward reference, not automatically promoted to an experimentally verified adhesome seed.

Seven added reciprocal controls expose the ambiguity of the new domains:

| Reference | Accession | Relevant shared domains |
|---|---|---|
| ENPP1 | P22413 | Somatomedin-B |
| ENPP2 | Q13822 | Somatomedin-B |
| Hemopexin | P02790 | Hemopexin-like |
| MMP2 | P08253 | FNII and hemopexin-like |
| MMP9 | P14780 | FNII and hemopexin-like |
| Tissue-type plasminogen activator (PLAT) | P00750 | FNI |
| PLA2R1 | Q13018 | FNII |

The legacy competitor label `non_adhesome` means outside the currently assigned target families; it is not proof that a protein never participates in adhesion biology. The panel is targeted and human-weighted, especially for ECM and ILK. No nearer-taxon ECM or ILK references were newly claimed in this update. Failure to support a schistosome family must not be interpreted as evolutionary absence.

## Calibration and actual searches

The completed [v4 rerun](../results/species_reference_v4/DOMAIN_RESULTS.md) supports two alpha-integrin, one FAK-family and one ILK candidate per species. No fibronectin-like or vitronectin-like candidate is currently supported; shared-domain hits remain in the screening output. Retained pools contain 1131 Sman, 1087 Shae and 1087 Sjap proteins: increases of 3, 4 and 5 respectively, with zero losses from v3. All 61 regression tests pass and the reference set rebuilds offline. Exact IDs and competing evidence are in [candidate review](../results/species_reference_v4/candidate_review.csv).

`results/reference_calibration_v4/architecture_calibration.csv` records every relevant positive/competitor architecture. All 10 alpha-integrin, 7 FAK, 1 ILK, 1 fibronectin and 1 vitronectin reference architectures were recovered. No other panel member met the complete architecture for these target rules. This is calibration on the same small panel used to develop the rules; it does not establish independent sensitivity/specificity or phylogenetic orthology.

## Reproduction after cleanup

Run from `01_detection` in the WSL environment with HMMER and BLAST+:

```bash
python scripts/run_domain_v4.py --cpu 8
python scripts/calibrate_domain_v4.py
python build_reference_sets.py --offline --strict-reviewed-curated --out results/reference_rebuild_check
python -B -m unittest discover -s tests -v
```

The full-run wrapper scans all 96 models and defaults to a fresh `results/species_reference_v4_full` directory. It does not depend on historical scans. The retained latest results were originally generated by verified incremental GA scanning; their original logs and comparative audit tables remain unchanged as provenance. Old baseline result trees have been removed. Prepared proteomes now reside within each latest species result directory. The one-time migration and comparison scripts have been removed.

## Primary sources

The following annotations were checked against cached UniProt records and/or the retrieved InterPro API metadata:

- [PF00039: fibronectin type I](https://www.ebi.ac.uk/interpro/entry/pfam/PF00039/) and [PF00040: fibronectin type II](https://www.ebi.ac.uk/interpro/entry/pfam/PF00040/); [human FN1](https://www.uniprot.org/uniprotkb/P02751/entry).
- [PF01033: somatomedin-B](https://www.ebi.ac.uk/interpro/entry/pfam/PF01033/) and [PF00045: hemopexin-like](https://www.ebi.ac.uk/interpro/entry/pfam/PF00045/); [human vitronectin](https://www.uniprot.org/uniprotkb/P04004/entry).
- [PF12796: ankyrin repeats (3 copies)](https://www.ebi.ac.uk/interpro/entry/pfam/PF12796/); [human ILK](https://www.uniprot.org/uniprotkb/Q13418/entry).
Additional motif, kindlin and PTP-PEST methods are retained below.

## Motif and targeted-family methods


The user-provided `references/motifs/adhesome_motif.json` is now an input to classification. The original JSON and helper script are preserved. Explicit mappings in `config/family_rules.json` distinguish sequence patterns from domain descriptions and residue positions. Each comparative run writes `motif_observations.csv`, with one-based inclusive coordinates, context, detected/not-detected/not-tested states, and the source library path. Both the rule file and motif library are checksummed in the classification report.

Motifs annotate candidates already recruited by domains or homology. They never recruit a protein by themselves, establish experimental adhesome membership, or identify a vertebrate integrin subtype. Missing motifs do not discard potentially divergent or truncated proteins.

| Library entry | Implementation | Interpretation limit |
|---|---|---|
| Alpha-tail GFFKR | Evolutionary `GFF.R` pattern in final 100 residues; exact matched sequence reported | Includes GFFKR, GFFNR and GFFHR; supports an existing alpha-integrin hypothesis with reciprocal homology; terminal location does not prove cytoplasmic topology |
| Beta-tail NPxY / NxxY | `N(?:P.|..)[YF]` in final 100 residues | Explicit Tyr/Phe extension of the supplied pattern; short motif alone is nonspecific |
| MIDAS DxSxS | Within a detected beta-integrin PF00362 domain | Optional context annotation; absence of the domain means not tested, not negative |
| Collagen GFOGER | Search encoded `GFPGER` within PF01391 | Hydroxyproline is a modification, not a genomically encoded `O`; detection does not establish hydroxylation or ligand binding |
| RGD | Within a detected FNIII domain | Supports sequence inspection, not fibronectin identity or demonstrated integrin binding |
| FAK Y397 | No position-397 test | Requires an alignment to transfer the human site; no phosphosite claim is made |
| PxxP, YxxP, LD, PHSRN, LDV | Not automatic family-assignment criteria | Short/nonspecific, insufficiently defined, or unsupported by the current family catalogue; retained in the source library |
| Named domains | Assessed by Pfam architecture | Descriptive domain labels are not sequence regular expressions |

## Architecture corrections

**Alpha integrin:** retain variable FG-GAP repeat counts; add PF08441, PF20805, PF20806 and PF00357 to screening. Two detected FG-GAP repeats, FG-GAP plus an alpha-specific stalk/tail model, or the configured stalk pair satisfy alternative screening architectures. A single FG-GAP repeat is neither diagnostic nor sufficient for a supported assignment. Alpha-specific domain evidence or a terminal GFFxR can corroborate an incomplete architecture, but supported assignments still require reciprocal family evidence. No alpha-I/VWA requirement is imposed: alpha-I domains are not universal. These are pragmatic screening alternatives, not a claim that a partial sequence reconstructs a complete receptor.

**FAK:** PF03623 (`Focal_AT`) is the relevant focal-adhesion targeting model. The previous PF02259 (`FAT`) model describes the different PIKK-associated FAT region. PF02259 remains a generic domain in the broader screen but contributes no FAK evidence. Added PF21477 (FAK FERM-C) and PF18038 (FERM-N). With reciprocal family evidence, comparative mode accepts either the full N-terminal FERM → kinase → focal-adhesion targeting order, or the diagnostic PF21477 → kinase core. The second route preserves a supported family identity while flagging the missing FAT domain; generic FERM+kinase alone does not qualify. `order_check` describes the full expected architecture and `support_order_check` records whether a permitted diagnostic core is ordered. Independent phylogeny is an additional rescue route. Strict mode retains its full-architecture requirement. The family label includes PTK2/PTK2B/kin-32 relatives and is not a one-to-one PTK2 claim; catalytic inactivity is not an exclusion.

**Kindlin:** added PF18124 (kindlin N-terminal domain). A kindlin-specific model plus FERM or PH supports the expected architecture, in that order. Generic FERM+PH alone no longer establishes a kindlin architecture. A missing PH call is allowed because the PH insertion can interrupt FERM detection.

**PTP-PEST:** forward references now include PTPN12, PTPN18 and PTPN22. Reciprocal competitors include other cytoplasmic/receptor PTPs. A shared PF00102 catalytic domain or PEST/proline-rich composition cannot resolve the family. In addition to reciprocal evidence, a supported call requires a BLAST alignment span reaching at least 60 residues beyond the final catalytic domain, or independent phylogenetic support. This threshold is an explicit conservative heuristic, not an empirically validated boundary; span endpoints do not count exact ungapped aligned residues. Older homology tables lacking alignment coordinates remain provisional for this family until rerun.




- [Pfam PF03623: focal adhesion targeting region](https://www.ebi.ac.uk/interpro/entry/pfam/PF03623/), including the structural study of the paxillin-binding FAT bundle.
- [UniProt human FAK/PTK2 Q05397](https://www.uniprot.org/uniprotkb/Q05397/entry): FERM and Focal_AT cross-references.
- [UniProt kindlin-1 Q9BQL6](https://www.uniprot.org/uniprotkb/Q9BQL6/entry): interrupted FERM-domain detection and kindlin-specific annotation.
- [UniProt integrin alpha-5 P08648](https://www.uniprot.org/uniprotkb/P08648/entry): FG-GAP plus alpha-integrin Ig-like domains.
- [Experimental analysis of UNC-112/kindlin](https://pmc.ncbi.nlm.nih.gov/articles/PMC7747009/): invertebrate kindlin context.
- [Drosophila alpha-PS2 tail experiments](https://pmc.ncbi.nlm.nih.gov/articles/PMC5155636/) and [oyster integrin-family study](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2020.01420/full): evolutionary GFFxR patterns; these do not validate the function of a particular schistosome motif.
