#!/usr/bin/env python3
"""Extend the Pfam resource without invoking hmmpress, scans or any analysis.

Existing model records/thresholds and analytical code/configuration are preserved.
Temporary build files are validated then replaced; no backups are created.
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

import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import csv
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import re
import sys
import tempfile
import time
from urllib.request import Request, urlopen
from uuid import uuid4

ROOT=_LAYOUT_BUNDLE_ROOT
sys.path.insert(0,str(_layout_resource(ROOT/'detection')))
from adhesome_core import digest,hmm_records,pf
from build_adhesome_pfam_library import build


def records(path):
    result={}
    for meta,block in hmm_records(path):
        accession=pf(meta['ACC'][0])
        if accession in result:raise ValueError('Duplicate HMM: '+accession)
        if not block.startswith('HMMER3/') or not re.search(r'^ALPH\s+amino\s*$',block,re.M):
            raise ValueError('Not a protein HMMER3 model: '+accession)
        length=int(meta['LENG'][0])
        states=[int(n) for n in re.findall(r'^\s*(\d+)\s+',block,re.M)]
        if states!=list(range(1,length+1)):raise ValueError('Incomplete HMM states: '+accession)
        if 'GA' in meta and (len(meta['GA'])!=2 or not all(math.isfinite(float(x.rstrip(';'))) for x in meta['GA'])):
            raise ValueError('Invalid GA thresholds: '+accession)
        result[accession]=(meta,block)
    return result


def expand(root=ROOT,check_only=False):
    root=_layout_resource(Path(root)).resolve();base=_layout_resource(root/'detection/databases/pfam');cache=_layout_resource(base/'source_v1.0')
    library=_layout_resource(base/'adhesome_v1.0.hmm');metadata=_layout_resource(Path(str(library)+'.metadata.json'))
    manifest=_layout_resource(root/'detection/config/adhesome_pfam_manifest.csv');rules_path=_layout_resource(root/'detection/config/family_rules.json')
    crosswalk=_layout_resource(root/'detection/config/adhesome_component_crosswalk.csv')
    source=_layout_resource(cache/'combined_source.hmm');provenance=_layout_resource(cache/'download_provenance.json')
    old_meta=json.loads(_layout_resource(metadata).read_text());old_records=records(library)
    if digest(library)!=old_meta['library_sha256']:raise ValueError('Existing library hash mismatch')
    if digest(manifest)!=old_meta['manifest_sha256']:raise ValueError('Existing manifest hash mismatch')
    if set(old_records)!=set(old_meta['models']):raise ValueError('Existing model metadata mismatch')
    with _layout_resource(manifest).open(encoding='utf-8-sig',newline='') as h:
        reader=csv.DictReader(h);fields=reader.fieldnames;rows=list(reader)
    listed={a for row in rows for a in re.findall(r'PF\d{5}',row['Pfam accession'])}
    if listed!=set(old_records):raise ValueError('Manifest/model membership mismatch')
    rules=json.loads(_layout_resource(rules_path).read_text());required=set(re.findall(r'PF\d{5}',json.dumps(rules['families'])))
    requested_path=_layout_resource(root/'maintenance/pfam_requested_additions.json')
    requested_data=json.loads(requested_path.read_text(encoding='utf-8'))
    if requested_data.get('pipeline_version')!='1.0' or not isinstance(requested_data.get('accessions'),list):
        raise ValueError('Invalid requested Pfam addition list')
    requested=requested_data['accessions']
    if any(not isinstance(a,str) or not re.fullmatch(r'PF\d{5}',a) for a in requested) or len(requested)!=len(set(requested)):
        raise ValueError('Requested Pfam accessions must be unique PFxxxxx identifiers')
    missing=sorted((required|set(requested))-set(old_records))
    plan=dict(pipeline_version='1.0',existing_models=len(old_records),additional_models=len(missing),
        requested_additions=len(set(requested)-set(old_records)),final_models=len(old_records)+len(missing),added_accessions=missing)
    if check_only or not missing:return plan
    # Changes to any analytical input during retrieval must abort before publication.
    protected=[p for p in _layout_resource(root).rglob('*.py') if not {'inputs','runs','returned','results','__pycache__','maintenance'}&set(p.relative_to(root).parts)]
    protected += [rules_path,_layout_resource(root/'annotation/family_rules.json'),_layout_resource(root/'annotation/family_identity_rules.json'),
        _layout_resource(root/'annotation/stage2_evidence_rules.json'),_layout_resource(root/'annotation/audit_policy.json'),_layout_resource(root/'configs/isolated_workflow.json'),crosswalk]
    protected.append(requested_path)
    fingerprints={p:digest(p) for p in protected}
    additions_library=_layout_resource(base/'adhesome_v1.0.additions.hmm')
    additions_metadata=_layout_resource(base/'adhesome_v1.0.additions.metadata.json')
    if additions_library.exists()!=additions_metadata.exists():
        raise ValueError('Incomplete existing additive Pfam library and metadata pair')
    previous_addition_records=records(additions_library) if additions_library.exists() else {}
    previous_addition_metadata=json.loads(additions_metadata.read_text(encoding='utf-8')) if additions_metadata.exists() else {}
    if previous_addition_records and (digest(additions_library)!=previous_addition_metadata.get('hmm_sha256') or
            set(previous_addition_records)!=set(previous_addition_metadata.get('accessions',[])) or
            set(previous_addition_records)-set(old_records)):
        raise ValueError('Existing additive Pfam library checksum or membership mismatch')
    initial={p:digest(p) for p in (library,metadata,manifest,source,provenance)}
    for p in (additions_library,additions_metadata):
        if p.exists():initial[p]=digest(p)
    old_provenance=json.loads(_layout_resource(provenance).read_text());prior_entries={r['accession']:r for r in old_provenance['models']}
    if set(prior_entries)!=set(old_records):raise ValueError('Incomplete existing download provenance')
    def fetch(accession):
        path=cache/(accession+'.hmm');url=f'https://www.ebi.ac.uk/interpro/api/entry/pfam/{accession}/?annotation=hmm'
        reused=_layout_resource(path).is_file()
        if not reused:
            for attempt in range(3):
                try:
                    with urlopen(Request(url,headers={'User-Agent':'adhesome-workflow/1.0'}),timeout=60) as response:data=response.read()
                    if data[:2]==b'\x1f\x8b':data=gzip.decompress(data)
                    if not data.startswith(b'HMMER3/') or not data.rstrip().endswith(b'//'):raise ValueError('Invalid response for '+accession)
                    path.write_bytes(data);break
                except Exception:
                    if attempt==2:raise
                    time.sleep(1+attempt)
        found=records(path)
        if set(found)!={accession}:raise ValueError('Downloaded accession differs: '+accession)
        meta,block=found[accession]
        print('Validated '+meta['ACC'][0],flush=True)
        return accession,block,dict(accession=accession,model_version=meta['ACC'][0],url=url,
            sha256=digest(path),retrieved_utc=datetime.now(timezone.utc).isoformat(),
            retrieval='existing_download_cache' if reused else 'EMBL_EBI_InterPro_Pfam_API')
    with ThreadPoolExecutor(max_workers=4) as pool:additions=dict((a,(b,p)) for a,b,p in pool.map(fetch,missing))
    # Preserve every original manifest byte; append only the newly required entries.
    text=io.StringIO(newline='');writer=csv.DictWriter(text,fieldnames=fields,lineterminator='\r\n' if b'\r\n' in _layout_resource(manifest).read_bytes() else '\n')
    for accession in missing:
        users=[name for name,rule in rules['families'].items() if accession in json.dumps(rule)]
        writer.writerow({'Pfam accession':accession,'Domain or family':records(cache/(accession+'.hmm'))[accession][0]['NAME'][0],
            'Adhesome role':('Current family-rule coverage: '+'; '.join(sorted(users))) if users else
                             'User-requested additive Pfam profile; no family-rule assignment added'})
    prefix=_layout_resource(manifest).read_bytes()
    if not prefix.endswith(b'\n'):prefix+=b'\n'
    expanded_manifest=prefix+text.getvalue().encode('utf-8')
    combined=''.join(block for meta,block in old_records.values())+''.join(additions[a][0] for a in missing)
    with tempfile.TemporaryDirectory(prefix='.pfam_build_',dir=base) as temp:
        temp=_layout_resource(Path(temp));staged_manifest=_layout_resource(temp/'manifest.csv');staged_source=_layout_resource(temp/'source.hmm');staged_library=temp/library.name
        staged_manifest.write_bytes(expanded_manifest);staged_source.write_bytes(combined.encode('utf-8'))
        found=build(staged_manifest,staged_source,staged_library,crosswalk=crosswalk)
        new_records=records(staged_library)
        if set(new_records)!=set(old_records)|set(missing):raise ValueError('Expanded model membership differs')
        if any(new_records[a]!=old_records[a] or found[a]!=old_meta['models'][a] for a in old_records):
            raise ValueError('An existing HMM or threshold changed')
        split={kind:records(_layout_resource(Path(str(staged_library)+'.'+kind+'.hmm'))) for kind in ('ga','fallback')}
        if set(split['ga'])&set(split['fallback']) or set(split['ga'])|set(split['fallback'])!=set(new_records):
            raise ValueError('GA/fallback split is inconsistent')
        if any(digest(p)!=expected for p,expected in {**fingerprints,**initial}.items()):
            raise ValueError('An analytical input or library changed during retrieval')
        staged_meta=_layout_resource(Path(str(staged_library)+'.metadata.json'));new_meta=json.loads(_layout_resource(staged_meta).read_text())
        staged_additions=_layout_resource(temp/'adhesome_v1.0.additions.hmm')
        prior_addition_bytes=additions_library.read_bytes() if additions_library.exists() else b''
        staged_additions.write_bytes(prior_addition_bytes+''.join(additions[a][0] for a in missing).encode('utf-8'))
        staged_additions_metadata=_layout_resource(temp/'adhesome_v1.0.additions.metadata.json')
        cumulative_additions=sorted(set(previous_addition_records)|set(missing))
        staged_additions_metadata.write_text(json.dumps(dict(pipeline_version='1.0',accessions=cumulative_additions,
            previous_additions=sorted(previous_addition_records),new_accessions=missing,
            model_count=len(cumulative_additions),hmm_sha256=digest(staged_additions),source_manifest_sha256=digest(staged_manifest),
            note='Additive source for this expansion; the main adhesome_v1.0.hmm contains these models merged with the preserved previous library.'),indent=2)+'\n')
        receipt=dict(**plan,date=datetime.now(timezone.utc).date().isoformat(),completed_utc=datetime.now(timezone.utc).isoformat(),
            original_model_records_and_ga_unchanged=True,analytical_files_unchanged=True,main_profile='calibrated_v1',
            previous_library_sha256=initial[library],previous_manifest_sha256=initial[manifest],
            library_sha256=digest(staged_library),manifest_sha256=digest(staged_manifest),source_sha256=digest(staged_source),
            ga_models=len(split['ga']),fallback_models=len(split['fallback']),
            additions_library=str(additions_library.relative_to(root).as_posix()),
            additions_library_sha256=digest(staged_additions),
            additions=[additions[a][1] for a in missing],
            preserved_model_records_sha256={a:hashlib.sha256(old_records[a][1].encode()).hexdigest() for a in old_records},
            protected_files_sha256={p.relative_to(root).as_posix():v for p,v in fingerprints.items()},
            analysis_run=False,hmmpress_run=False,backups_created=False)
        new_meta.update(source='source_v1.0/combined_source.hmm',source_path_base='metadata_directory',
            manifest='../../config/adhesome_pfam_manifest.csv',manifest_path_base='metadata_directory',
            family_rules_sha256=digest(rules_path),model_count=len(found),
            expansion={k:receipt[k] for k in ('date','existing_models','additional_models','final_models','added_accessions',
                'original_model_records_and_ga_unchanged','main_profile','previous_library_sha256','analysis_run','backups_created')})
        staged_meta.write_text(json.dumps(new_meta,indent=2)+'\n')
        new_provenance=dict(old_provenance)
        new_provenance.update(updated_utc=receipt['completed_utc'],manifest_sha256=digest(staged_manifest),
            models=[*old_provenance['models'],*(additions[a][1] for a in missing)],
            expansion_note='Only missing models retrieved; existing model versions, payloads and prior retrieval provenance preserved')
        staged_provenance=_layout_resource(temp/'download_provenance.json');staged_provenance.write_text(json.dumps(new_provenance,indent=2)+'\n')
        for staged,target in [(staged_manifest,manifest),(staged_source,source),(staged_library,library),
            *( (_layout_resource(Path(str(staged_library)+'.'+k+'.hmm')),_layout_resource(Path(str(library)+'.'+k+'.hmm'))) for k in ('ga','fallback')),
            (staged_meta,metadata),(staged_provenance,provenance),
            (staged_additions,additions_library),(staged_additions_metadata,additions_metadata)]:
            # A Windows TemporaryDirectory has a private ACL. Publish through a
            # normal adjacent file so the final resource inherits its own folder's
            # read permissions instead of the private build directory's ACL.
            publication=_layout_resource(target.with_name(target.name+'.publishing-'+uuid4().hex))
            try:
                with _layout_resource(publication).open('xb') as h:h.write(_layout_resource(staged).read_bytes())
                publication.replace(target)
            finally:
                if _layout_resource(publication).exists():publication.unlink()
    for p,expected in fingerprints.items():
        if digest(p)!=expected:raise ValueError('Protected analytical file changed: '+str(p))
    (_layout_resource(root/'maintenance/pfam_expansion_latest.json')).write_text(json.dumps(receipt,indent=2)+'\n')
    return {k:v for k,v in receipt.items() if k not in ('additions','protected_files_sha256','preserved_model_records_sha256')}


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--check-only',action='store_true')
    a=p.parse_args();print(json.dumps(expand(check_only=a.check_only),indent=2))


if __name__=='__main__':main()
