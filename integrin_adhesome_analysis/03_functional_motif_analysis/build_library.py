"""Rebuild the reviewed local motif panel; no network access required."""
import json
from pathlib import Path

LOCAL = '../01_detection/config/family_rules.json'
motifs = {}


def add(name, pattern, context, regions, interpretation, source=LOCAL):
    motifs[name] = dict(pattern=pattern, context=context, regions=regions,
                        interpretation=interpretation, source=source)
    if context == 'c_terminal_window':
        motifs[name]['window'] = 100


add('alpha_signature', '[FYWS][RK].GFF.R', 'c_terminal_window', ['cytoplasmic_tail'],
    'PROSITE alpha signature; independent architecture and topology required.', 'https://prosite.expasy.org/PS00242')
add('alpha_GFFxR', 'GFF.R', 'c_terminal_window', ['cytoplasmic_tail'],
    'Relaxed alpha core, lower specificity than full signature; tail window is not topology.')
add('beta_NPxY', 'NP.Y', 'c_terminal_window', ['cytoplasmic_tail'],
    'Candidate adaptor site; membrane proximal/distal assignment requires topology and alignment.')
add('beta_NxxYF', 'N..[YF]', 'c_terminal_window', ['cytoplasmic_tail'],
    'Broad variant screen including NPxY; nested hits are not independent evidence or binding proof.')
add('MIDAS_core', 'D.S.S', 'domain', ['PF00362', 'PF00092'],
    'Partial metal-binding pattern; domain and other coordinating residues required. Alpha-I is optional.')
add('RGD', 'RGD', 'domain', ['PF00041'],
    'Short ligand candidate; FNIII context and structural accessibility required; not fibronectin or binding proof.')
add('PHSRN', 'PHSRN', 'domain', ['PF00041'],
    'Canonical fibronectin synergy screen; absence in divergent candidates does not exclude family.',
    '../01_detection/references/motifs/adhesome_motif.json')
add('GFPGER_precursor', 'GFPGER', 'domain', ['PF01391'],
    'Encoded GFOGER precursor; hydroxylation and triple-helical context unknown. O is not hydroxyproline.')
add('collagen_GXY_run', '(?:G..){5,}', 'whole_sequence', [],
    'Heuristic: at least five GXY units; repeat composition alone does not establish collagen identity.',
    'Local exploratory repeat rule, not a curated diagnostic signature')
add('PxxP', 'P..P', 'whole_sequence', [],
    'Generic proline-rich candidate site; not family-specific or sufficient for SH3/EVH1 binding.',
    '../01_detection/references/motifs/adhesome_motif.json')
add('EVH1_FPPPP', '[FW]PPPP', 'whole_sequence', [],
    'Conservative EVH1 ligand core screen; flanks and accessibility matter; not zyxin-specific.',
    'https://pmc.ncbi.nlm.nih.gov/articles/PMC5724578/')
add('PTP_active_site', '[LIVMF]HC..G..R[STC][STAGP]', 'domain', ['PF00102'],
    'PROSITE phosphatase active-site pattern; does not distinguish PTP-PEST or prove activity.',
    'https://prosite.expasy.org/PS00383')

families = {}
data = {
    'actin': ([], 'Actin domain/profile and homology required.'),
    'alpha_actinin': ([], 'CH, spectrin-repeat and EF-hand architecture/homology required.'),
    'cofilin': ([], 'ADF/cofilin domain and homology required; no residue-number transfer.'),
    'collagen_like': (['GFPGER_precursor', 'collagen_GXY_run'], 'Repeat/precursor inspection; hydroxylation unknown.'),
    'FAK': (['PxxP'], 'FERM-kinase-FAT architecture required. Human Y397 cannot be transferred without alignment.'),
    'fibronectin_like_candidate': (['RGD', 'PHSRN'], 'Candidate label retained; FNIII alone is not fibronectin orthology.'),
    'filamin': ([], 'CH and Ig-repeat architecture/homology required; no arbitrary repeat-21 assignment.'),
    'ILK': ([], 'Ankyrin and pseudokinase architecture required; catalytic consensus absence is not rejection.'),
    'integrin_alpha': (['alpha_signature', 'alpha_GFFxR', 'MIDAS_core'], 'Alpha architecture and TM topology required; alpha-I/MIDAS is optional.'),
    'integrin_beta': (['beta_NPxY', 'beta_NxxYF', 'MIDAS_core'], 'Tail and beta-I context corroborate but do not assign family.'),
    'kindlin': ([], 'FERM/PH organization required; NxxY ligand belongs on integrin beta, not kindlin.'),
    'laminin_like': ([], 'Laminin domain architecture required; no universal short integrin-binding signature.'),
    'paxillin_like': ([], 'LIM architecture and alignment-based LD-helix inspection required; LD dipeptides are not diagnostic.'),
    'profilin': ([], 'Profilin profile/homology required; singleton discovery skipped.'),
    'PTP_PEST': (['PTP_active_site', 'PxxP'], 'Generic PTP patterns; PEST-family support requires noncatalytic homology.'),
    'Src': ([], 'SH3-SH2-kinase architecture required; generic kinase patterns cannot distinguish Src.'),
    'talin': ([], 'FERM/rod architecture required; mechanically exposed VBS not reliably inferred by short regex.'),
    'vinculin': (['EVH1_FPPPP'], 'Optional EVH1-ligand inspection; vinculin head/tail architecture required.'),
    'zyxin_like': (['EVH1_FPPPP', 'PxxP'], 'Generic interaction patterns; LIM architecture and homology required.')
}
for name, (checks, note) in data.items():
    families[name] = dict(motifs=checks, note=note)

if __name__ == '__main__':
    rules_path = Path(__file__).resolve().parent.parent/'01_detection/config/family_rules.json'
    rules = json.loads(rules_path.read_text(encoding='utf-8'))['families']
    for name, config in families.items():
        rule_name = 'fibronectin_like' if name == 'fibronectin_like_candidate' else name
        config['domain_alternatives'] = rules[rule_name].get('domain_alternatives', [])
        config['domain_rule_source'] = LOCAL
    Path(__file__).with_name('motif_library.json').write_text(json.dumps(dict(
        version='1.0.0', reviewed='2026-09-19',
        scope='Conservative inspection panel, not exhaustive signatures or a classifier.',
        motifs=motifs, families=families), indent=2)+'\n', encoding='utf-8')
