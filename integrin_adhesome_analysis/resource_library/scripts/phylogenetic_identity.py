"""Optional, traceable tree concordance; never infer one-to-one orthology from a HOG."""

# Resolve source resources independently of the caller working directory.
import sys as _layout_sys
from pathlib import Path as _LayoutPath
_LAYOUT_BUNDLE_ROOT = next(p for p in _LayoutPath(__file__).resolve().parents
                          if (p/'config/layout_paths.json').is_file() and (p/'pipeline_paths.py').is_file())
if str(_LAYOUT_BUNDLE_ROOT) not in _layout_sys.path:
    _layout_sys.path.insert(0, str(_LAYOUT_BUNDLE_ROOT))
from pipeline_paths import activate as _layout_activate, resource_path as _layout_resource
_layout_activate()

import csv
import hashlib
import math
from pathlib import Path


def sha(path):
    with _layout_resource(Path(path)).open('rb') as h:return hashlib.file_digest(h,'sha256').hexdigest()


def apply(rows,evidence=None,enabled=False):
    """Missing/failed optional evidence cannot prevent the completed family audit."""
    bank={};error='';path=_layout_resource(Path(evidence)) if evidence else None
    if enabled:
        try:
            if path is None:raise ValueError('No phylogenetic identity evidence supplied')
            with _layout_resource(path).open(encoding='utf-8-sig',newline='') as h:
                for item in csv.DictReader(h,delimiter='\t'):
                    key=item['protein_id']
                    if key in bank:raise ValueError('Duplicate phylogenetic candidate: '+key)
                    bank[key]=item
        except (OSError,ValueError,KeyError) as exc:bank={};error=str(exc)
    for row in rows:
        row['phylogenetic_identity_status']='not_requested' if not enabled else 'unavailable' if error else 'placement_unresolved'
        row['phylogenetic_identity_reason']=error or ('Optional resolution disabled' if not enabled else 'No qualifying placement')
        if not enabled or error:continue
        hit=bank.get(row['protein_id'])
        if not hit:continue
        try:
            tree=_layout_resource(Path(hit['tree_file']))
            if not tree.is_absolute():tree=path.parent/tree
            if sha(tree)!=hit['tree_sha256']:raise ValueError('Tree checksum differs')
            if not row.get('sequence_sha256') or hit['sequence_sha256']!=row['sequence_sha256']:
                raise ValueError('Candidate sequence hash differs')
            boot,alrt=float(hit['ufboot']),float(hit['sh_alrt'])
            if not math.isfinite(boot) or not math.isfinite(alrt) or not 95<=boot<=100 or not 80<=alrt<=100:
                raise ValueError('Weak or absent branch support')
            if int(hit['reference_count'])<2 or int(hit['reference_taxa_count'])<2:
                raise ValueError('Insufficient labelled reference or taxon coverage')
            if not hit.get('reference_label_source') or not hit.get('clade_id'):raise ValueError('Reference/clade provenance missing')
            if hit.get('placement_status')!='pure_supported_clade':raise ValueError('Mixed or unresolved clade')
            family=hit['family']
            # No automatic relabelling or vertebrate paralog transfer. A tree is
            # corroboration/review evidence at the stated level, not a new family gate.
            agrees=family==row.get('assigned_family')
            row.update(phylogenetic_identity_status='family_concordant' if agrees else 'conflicting_placement_review',
                phylogenetic_identity_reason='Supported labelled clade; original family call retained',
                phylogenetic_family=family,phylogenetic_clade_id=hit['clade_id'],
                phylogenetic_tree_sha256=hit['tree_sha256'],phylogenetic_branch_support=hit['sh_alrt']+'/'+hit['ufboot'],
                phylogenetic_duplication_status=hit.get('duplication_status','not_assessed'),
                phylogenetic_validation_scope='evidence_concordance_not_independent_accuracy')
            if not agrees:row['identity_outcome']='conflicting_identity'
            if agrees and hit.get('subfamily') and hit.get('taxon_scope')=='target_lineage' and hit.get('rooting_status')=='outgroup_supported' and hit.get('duplication_status')=='reconciled' and hit.get('relationship_evidence_source'):
                row.update(phylogenetic_subfamily=hit['subfamily'],evolutionary_identity_level='supported_subfamily')
                relationship=hit.get('orthology_relationship')
                if relationship in ('ortholog','co_ortholog','paralog'):
                    row['orthology_paralogy_status']=relationship
        except (OSError,KeyError,ValueError,TypeError) as exc:
            row.update(phylogenetic_identity_status='evidence_not_qualifying',phylogenetic_identity_reason=str(exc))
    return dict(requested=bool(enabled),status='unavailable' if error else 'complete' if enabled else 'not_requested',
        error=error,source=str(path) if enabled and path else None,source_sha256=sha(path) if enabled and path and _layout_resource(path).is_file() else None,
        identity_transfer='No automatic relabelling; no one-to-one orthology or paralog name inferred from topology alone')


def project_inventory(project):
    import json
    project=_layout_resource(Path(project));manifest=json.loads((_layout_resource(project/'manifest.json')).read_text())
    summary=json.loads((_layout_resource(project/'run_summary.json')).read_text()) if (_layout_resource(project/'run_summary.json')).is_file() else {}
    counts={}
    for row in summary.get('groups',[]):counts[row.get('status','unknown')]=counts.get(row.get('status','unknown'),0)+1
    trees=list(_layout_resource(project).glob('*/inference/gene_tree.treefile'))
    return dict(project=str(_layout_resource(project).resolve()),groups_prepared=len(manifest.get('groups',[])),
        inference_status=summary.get('status','unknown'),group_status_counts=counts,tree_files=len(trees),
        grouping=manifest.get('grouping'),duplication_evidence_status=manifest.get('duplication_evidence_status'),
        rooted_reconciled_orthology=False,reference_identity_labels_supplied=False,
        interpretation='Trees/tip maps can be reused. Curated reference identity labels and supported placements are required before identity adjudication.')
