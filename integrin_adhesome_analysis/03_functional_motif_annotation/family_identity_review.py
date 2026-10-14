"""Conservative architecture triage; never changes supplied family membership."""
import hashlib


def interval_count(regions, accession):
    intervals = sorted(set((a, b) for name, a, b in regions if name == accession))
    merged = []
    for a, b in intervals:
        if merged and a <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(b, merged[-1][1]))
        else:
            merged.append((a, b))
    return len(merged)


def review(families, regions):
    rows = []
    for family in ('PINCH_like_candidates','PINCH_like','actin','actin_conventional','alpha_actinin','parvin_like'):
        for record in families.get(family, []):
            rs = regions[family, record['id']]
            lim, ch, actin, spectrin = [interval_count(rs, p) for p in ('PF00412', 'PF00307', 'PF00022', 'PF00435')]
            pinch = family in ('PINCH_like_candidates','PINCH_like')
            priority = (lim == 5) if pinch else (ch >= 2 and not actin and not spectrin)
            rows.append(dict(family=family, sequence_id=record['id'], length=len(record['sequence']),
                sequence_sha256=hashlib.sha256(record['sequence'].encode()).hexdigest(),
                LIM_domains=lim, CH_domains=ch, actin_domains=actin, spectrin_domains=spectrin,
                hypothesis='PINCH' if pinch else 'possible_Parvin',
                review_priority='architecture_followup' if priority else 'unresolved_review',
                identity_status='unresolved_original_family_retained',
                missing_evidence='PINCH-specific LIM1/CDD;orthology;IPP partner context' if pinch else 'Parvin-specific CDD/homology;orthology;CH domain order;IPP partner context',
                interpretation='Domain counts are a triage screen only. Missing domains do not exclude identity; tandem CH also occurs in other proteins. Parvin-binding motifs in paxillin cannot identify Parvin.'))
    return rows
