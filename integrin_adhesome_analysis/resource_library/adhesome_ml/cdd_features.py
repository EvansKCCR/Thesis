"""CDD-fragment local similarity features, with calibration-only reference bank."""
import csv
import hashlib
import json
import re
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

BLAST_FIELDS = "qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen"
CLASSES = {"specific", "non-specific", "broad"}


def _token(text):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(text or "")).strip("_") or "unknown"


def _tsv(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def build_calibration_bank(library_dir, calibration_branch, out_dir):
    """Make a labelled bank from calibration hits only; validation rows are never read."""
    source_fasta = library_dir / "cdd_domain_sequences.faa"
    metadata = _tsv(library_dir / "cdd_domain_metadata.tsv")
    labels = {r["protein_id"]: r for r in _tsv(calibration_branch / "reports/predictions.tsv")
              if str(r.get("control_class", "")).casefold() in {"positive", "positive_control", "adhesome_associated"}}
    seqs = {}
    current = None
    with source_fasta.open(encoding="utf-8-sig") as fh:
        for line in fh:
            line=line.strip()
            if not line: continue
            if line.startswith(">"):
                current=line[1:].split()[0]; seqs[current]=[]
            elif current is None: raise ValueError("Malformed CDD sequence FASTA")
            else: seqs[current].append(line)
    seqs={k:"".join(v) for k,v in seqs.items()}
    bank_rows={}
    for row in metadata:
        if row.get("origin") != "calibration": continue
        ident=row.get("source_protein_id","")
        label=labels.get(ident,{}).get("expected_families","") or labels.get(ident,{}).get("trusted_family","")
        family=_token(label.split(";")[0].split("|")[0].strip()) if label else ""
        if not family: continue
        if row.get("hit_classification") not in CLASSES: continue
        sid=row["library_id"]
        item=bank_rows.setdefault(sid,{"library_id":sid,"families":set(),"classes":set(),"accessions":set(),"annotations":set(),"source_ids":set()})
        item["families"].add(family); item["classes"].add(row["hit_classification"])
        if row.get("cdd_accession"): item["accessions"].add(row["cdd_accession"])
        if row.get("domain_annotation"): item["annotations"].add(row["domain_annotation"])
        item["source_ids"].add(ident)
    if not bank_rows: raise ValueError("No calibration-labelled CDD domain references found")
    out_dir.mkdir(parents=True,exist_ok=True)
    bank_fasta=out_dir/"cdd_calibration_references.faa"
    with bank_fasta.open("w",encoding="utf-8",newline="\n") as fh:
        for sid,item in sorted(bank_rows.items()):
            if sid not in seqs: raise ValueError(f"CDD sequence missing for {sid}")
            fh.write(f">{sid}\n")
            sequence=seqs[sid]
            for i in range(0,len(sequence),60): fh.write(sequence[i:i+60]+"\n")
    bank_meta={k:{"families":sorted(v["families"]),"classes":sorted(v["classes"]),
        "accessions":sorted(v["accessions"]),"annotations":sorted(v["annotations"]),
        "source_ids":sorted(v["source_ids"])} for k,v in bank_rows.items()}
    (out_dir/"cdd_calibration_references.json").write_text(json.dumps(bank_meta,indent=2)+"\n",encoding="utf-8")
    return bank_fasta, bank_meta


def similarity_features(records, bank_fasta, bank_meta, blastp="blastp", exclude_self=False, workdir=None):
    """BLAST full protein queries against calibration domain fragments and encode hits."""
    if not records: return {}
    with tempfile.TemporaryDirectory(prefix="cdd_features_",dir=str(workdir) if workdir else None) as temp:
        root=Path(temp); query=root/"queries.faa"; output=root/"hits.tsv"
        with query.open("w",encoding="utf-8",newline="\n") as fh:
            for row in records:
                fh.write(f">{row['protein_id']}\n")
                seq=row["sequence"]
                for i in range(0,len(seq),60): fh.write(seq[i:i+60]+"\n")
        cmd=[blastp,"-query",str(query),"-subject",str(bank_fasta),"-seg","yes","-evalue","1e-3",
            "-max_target_seqs","500","-num_threads","1","-outfmt",f"6 {BLAST_FIELDS}","-out",str(output)]
        run=subprocess.run(cmd,capture_output=True,text=True)
        if run.returncode: raise RuntimeError(f"CDD reference blastp failed ({run.returncode}): {run.stderr[-2000:]}")
        hits=defaultdict(list)
        if output.exists():
            with output.open(encoding="utf-8") as fh:
                for line in fh:
                    f=line.rstrip("\n").split("\t")
                    if len(f)!=12: raise ValueError("Unexpected CDD BLAST output format")
                    q,sid=f[0],f[1]
                    if sid not in bank_meta: continue
                    pident,alen,qstart,qend,sstart,send,evalue,bitscore,qlen,slen=map(float,f[2:])
                    hits[q].append((sid,pident,alen,qstart,qend,sstart,send,evalue,bitscore,qlen,slen))
    result={}
    for row in records:
        ident=row["protein_id"]; seqhits=hits.get(ident,[]); features={"cdd:hit_count":0.0}
        retained=[]
        for h in seqhits:
            sid=h[0]; meta=bank_meta[sid]
            if exclude_self and ident in meta["source_ids"]: continue
            retained.append(h)
        features["cdd:hit_count"]=float(len(retained))
        per_class=defaultdict(list); per_family=defaultdict(list); per_accession=defaultdict(list); per_family_class=defaultdict(list); per_annotation=defaultdict(list)
        for h in retained:
            sid,pident,alen,qs,qe,ss,se,evalue,bitscore,qlen,slen=h; meta=bank_meta[sid]
            qcov=min(1.0,abs(qe-qs)+1.0)/qlen if qlen else 0.0
            scov=min(1.0,abs(se-ss)+1.0)/slen if slen else 0.0
            for cl in meta["classes"]:
                per_class[cl].append((bitscore,pident,qcov,scov,qs,qe))
                for family in meta["families"]: per_family_class[(cl,family)].append((bitscore,pident,qcov,scov,qs,qe))
            for family in meta["families"]: per_family[family].append((bitscore,pident,qcov,scov,qs,qe))
            for acc in meta["accessions"]: per_accession[acc].append((bitscore,pident,qcov,scov,qs,qe))
            for annotation in meta.get("annotations",[]): per_annotation[annotation].append((bitscore,pident,qcov,scov,qs,qe))
        for cl,values in per_class.items():
            best=max(values,key=lambda x:x[0]); features[f"cdd:class:{cl}:hit_count"]=float(len(values))
            features[f"cdd:class:{cl}:best_bitscore"]=best[0]; features[f"cdd:class:{cl}:best_identity"]=best[1]
            features[f"cdd:class:{cl}:query_coverage"]=best[2]; features[f"cdd:class:{cl}:subject_coverage"]=best[3]
        ordered=[]
        for family,values in per_family.items():
            values.sort(key=lambda x:x[4]); best=max(values,key=lambda x:x[0])
            token=_token(family); features[f"cdd:family:{token}:hit_count"]=float(len(values))
            features[f"cdd:family:{token}:best_bitscore"]=best[0]; features[f"cdd:family:{token}:best_identity"]=best[1]
            features[f"cdd:family:{token}:query_coverage"]=best[2]
            features[f"cdd:family:{token}:domain_copy_number"]=float(len(values))
            ordered.append((min(x[4] for x in values),token))
        for (cl,family),values in per_family_class.items():
            best=max(values,key=lambda x:x[0]); key=f"cdd:evidence:{cl}:family:{_token(family)}"
            features[key+":hit_count"]=float(len(values)); features[key+":best_bitscore"]=best[0]
            features[key+":best_identity"]=best[1]; features[key+":query_coverage"]=best[2]
        for acc,values in per_accession.items():
            best=max(values,key=lambda x:x[0]); token=_token(acc)
            features[f"cdd:accession:{token}:hit_count"]=float(len(values))
            features[f"cdd:accession:{token}:best_bitscore"]=best[0]
        for annotation,values in per_annotation.items():
            best=max(values,key=lambda x:x[0]); token=_token(annotation)
            features[f"cdd:annotation:{token}:hit_count"]=float(len(values))
            features[f"cdd:annotation:{token}:best_bitscore"]=best[0]
            features[f"cdd:annotation:{token}:best_identity"]=best[1]
            features[f"cdd:annotation:{token}:query_coverage"]=best[2]
            features[f"cdd:annotation:{token}:n_start"]=min(x[4] for x in values)
            features[f"cdd:annotation:{token}:c_end"]=max(x[5] for x in values)
        ordered.sort()
        for (leftpos,left),(rightpos,right) in zip(ordered,ordered[1:]):
            if left!=right: features[f"cdd:family_order:{left}>{right}"]=features.get(f"cdd:family_order:{left}>{right}",0.0)+1.0
        result[ident]=features
    return result

