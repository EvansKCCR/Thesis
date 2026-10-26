"""Describe search coverage without turning missing evidence into a negative call."""
from collections import Counter
import json
from pathlib import Path

from adhesome_core import digest, table


REPORT_FIELDS = ['homology_check', 'reference_count', 'similarity_reason', 'reciprocal_reason',
                 'related_reference_count', 'related_reference_families', 'related_similarity',
                 'related_reciprocal', 'related_evidence_sources']


def search_coverage(homology, proteome):
    """Only a completed, input-matched paired search can document a missing hit.

    Read the run-local snapshot, so moved runs do not need their original absolute
    reference paths. Standalone imported evidence has unknown search coverage.
    """
    if not homology:
        return None, []
    folder = Path(homology).resolve().parent
    manifest = folder/'homology.json'
    if not manifest.exists():
        return None, []
    meta = json.loads(manifest.read_text(encoding='utf-8'))
    commands = meta.get('commands', [])
    if len(commands) != 2:
        raise ValueError('Homology completion record requires both forward and reverse searches')
    def argument(command, flag):
        return command[command.index(flag)+1]
    try:
        target = argument(commands[0], '-subject')
        if argument(commands[1], '-query') != target:
            raise ValueError('Paired homology searches have different targets')
    except (ValueError, IndexError) as exc:
        raise ValueError('Invalid paired homology completion record') from exc
    if meta.get('checksums', {}).get(target) != digest(proteome):
        raise ValueError('Homology search target checksum does not match the classified proteome')
    snapshot = folder/'reference_snapshot'
    paths = [manifest]
    for name, expected in meta.get('artifact_checksums', {}).items():
        if name not in ('forward.tsv', 'reverse.tsv', 'homology.csv'):
            raise ValueError('Unknown homology artifact: '+name)
        path = folder/name
        if not path.is_file() or digest(path) != expected:
            raise ValueError('Homology search artifact missing or changed: '+str(path))
        paths.append(path)
    for name in ('search_catalogue.csv', 'search_references.faa', 'metazoan_catalogue.csv', 'metazoan_references.faa'):
        path = snapshot/name
        if not path.is_file() or meta.get('snapshot_checksums', {}).get(name) != digest(path):
            raise ValueError('Homology reference snapshot missing or changed: '+str(path))
        paths.append(path)
    return Counter(row['family'] for row in table(snapshot/'search_catalogue.csv')), paths


def related_families(family, rules, aliases):
    """Ancestors and legacy aliases provide context, never subfamily identity."""
    lineage = {family}
    parent = rules[family].get('parent_family')
    while parent and parent not in lineage:
        lineage.add(parent)
        parent = rules.get(parent, {}).get('parent_family')
    return (lineage | {name for name, destination in aliases.items() if destination in lineage}) - {family}


def coverage_details(family, states, support, related, related_support, counts, state_function):
    direct = counts.get(family, 0) if counts is not None else ''
    def reason(key):
        if states[key] != 'not_tested':
            return 'family_specific_evidence_rows'
        if counts is None:
            return 'search_completion_not_documented'
        if not direct:
            return 'no_family_reference'
        if key == 'similarity':
            return 'completed_search_no_qualifying_hit'
        if not any(row.get('similarity', '').lower() == 'true' for row in support):
            return 'no_qualifying_forward_hit_for_reciprocal_assessment'
        return 'reciprocal_result_not_supplied'
    return dict(homology_check='completed' if counts is not None else 'not_documented',
                reference_count=direct, similarity_reason=reason('similarity'), reciprocal_reason=reason('reciprocal'),
                related_reference_count=sum(counts.get(name, 0) for name in related) if counts is not None else '',
                related_reference_families=';'.join(sorted(name for name in related if
                    (counts is not None and counts.get(name)) or any(r['family'] == name for r in related_support))),
                related_similarity=state_function(related_support, 'similarity'),
                related_reciprocal=state_function(related_support, 'reciprocal'),
                related_evidence_sources=';'.join(sorted({r['source'] for r in related_support})))
