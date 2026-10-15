#!/usr/bin/env python3
"""Archive unreferenced task-owned DE CAS shards after a final catalog closes.

No scientific/base data file is modified. Ownership requires our exact JS
wrapper and content-address digest; every current catalog/base shard is kept.
"""
import argparse,base64,gzip,hashlib,json,re,time
from pathlib import Path

KINDS={"download","de_rows","core_reliability","review","ca1_program_review","fits"}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def owned(p):
    text=p.read_text();a=re.search(r'^;\(function\(w\)\{var id=("[^"]+"),entry=(\{"encoding":"gzip-base64","data":"[A-Za-z0-9+/=]+"\});',text)
    if not a:return None
    ident=json.loads(a.group(1));entry=json.loads(a.group(2))
    if ident!=p.stem or ident.split("__",1)[0] not in KINDS:return None
    raw=gzip.decompress(base64.b64decode(entry["data"],validate=True));digest=hashlib.sha256(raw).hexdigest()
    if ident.startswith("fits__"):
        if ident.rsplit("__",1)[-1]!=digest[:16]:return None
    elif ident.split("__",1)[1]!=digest:return None
    expected=";(function(w){var id="+json.dumps(ident)+",entry="+json.dumps(entry,separators=(",",":"))+";(w.SATB2_SHARDS=w.SATB2_SHARDS||{})[id]=entry;if(typeof w.registerShard==='function'){w.registerShard(id,entry.data);}})(window);\n"
    if text!=expected:return None
    return digest
def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument("--bundle-root",required=True,type=Path);args=ap.parse_args()
    root=args.bundle_root.resolve(strict=True);catp=root/"DE_CATALOG.json";rp=root/"DE_CATALOG_COMPLETE.json";mp=root/"data/manifest.json"
    cat=json.loads(catp.read_text());receipt=json.loads(rp.read_text());manifest=json.loads(mp.read_text())
    if cat.get("status")!="COMPLETE" or receipt.get("status")!="COMPLETE":raise ValueError("Final closed catalog is required")
    keep={r["path"] for r in cat["shards"].values()}|{r["path"] for r in manifest["shards"].values()}
    for relative,digest in receipt["output_sha256"].items():
        if sha(root/relative)!=digest:raise ValueError("Current catalog output changed")
    archive=root.parent/"offline_dashboard_catalog_history"/(time.strftime("%Y%m%dT%H%M%SZ",time.gmtime())+"__verified_aborted_orphans__"+sha(catp)[:12])
    moved=[];skipped=[]
    for p in sorted((root/"data/shards").glob("*.js")):
        relative=str(p.relative_to(root))
        if relative in keep:continue
        if p.is_symlink():raise ValueError("Refusing orphan symlink")
        digest=owned(p)
        if digest is None:skipped.append(relative);continue
        file_digest=sha(p);size=p.stat().st_size;dst=archive/relative;dst.parent.mkdir(parents=True,exist_ok=True)
        if dst.exists():raise ValueError("Archive destination exists")
        p.rename(dst)
        if sha(dst)!=file_digest:raise ValueError("Archive bytes changed")
        moved.append({"path":relative,"sha256":file_digest,"payload_sha256":digest,"bytes":size})
    # Verify that cleanup preserved every current scientific display/download.
    for relative,digest in receipt["output_sha256"].items():
        if sha(root/relative)!=digest:raise ValueError("Current output changed during orphan archival")
    archive.mkdir(parents=True,exist_ok=True)
    result={"status":"PASS","stage":"archive_task_owned_unreferenced_DE_shards","catalog_sha256":sha(catp),"native_outputs_modified":False,"base_assets_modified":False,"scientific_fitting_performed":False,"archived_files":moved,"archived_bytes":sum(r["bytes"] for r in moved),"unrecognized_unreferenced_files_preserved":skipped,"current_catalog_outputs_reverified":True}
    (archive/"ORPHAN_ARCHIVE_RECEIPT.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({"status":"PASS","archived_files":len(moved),"archived_MiB":round(result["archived_bytes"]/1024**2,2),"receipt":str(archive/"ORPHAN_ARCHIVE_RECEIPT.json")},indent=2))
if __name__=="__main__":main()
