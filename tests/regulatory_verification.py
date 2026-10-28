#!/usr/bin/env python3
"""Pure synthetic transport/cohort verification; no scientific code or data."""
import base64
import copy
import gzip
import hashlib
import importlib.util
import json
import tempfile
from pathlib import Path
from unittest.mock import patch

REPO = Path(__file__).resolve().parent.parent
spec = importlib.util.spec_from_file_location("verify_regulatory", REPO / "scripts/verify_regulatory_atlas.py")
v = importlib.util.module_from_spec(spec); spec.loader.exec_module(v)


def reject(call, phrase):
    try:
        call()
    except ValueError as error:
        assert phrase in str(error), (phrase, str(error))
    else:
        raise AssertionError("Expected rejection: " + phrase)


with tempfile.TemporaryDirectory(prefix="SYNTHETIC_REGULATORY_VERIFIER_ONLY_") as tmp:
    site = Path(tmp) / "site/regulatory"; site.mkdir(parents=True)
    payload = {"schema":"satb2-regulatory-atlas/v1","analysis_id":"SYNTHETIC_VERIFIER_TEST_ONLY","source_fingerprint":"a"*64,"status":"partial",
               "endpoints":[{"id":ep,"family":a,"source_endpoint_id":b,"n_nuclei":c,"membership_sha256":d} for ep,(a,b,c,d) in v.EXPECTED.items()],
               "samples":[{"id":m,"genotype":g} for m,g in v.MICE.items()],"features":[],"values":[],"embeddings":[],"downloads":[],
               "provenance":{"synthetic":True,"TF_activity":"COMPLETE_CURATED_COLLECTRI_ULM","motif_supported_regulons":"BLOCKED_SYNTHETIC_TEST_ONLY","no_scientific_fits_in_publication_builder":True,"components":{"hdWGCNA":{},"curated_tf_activity":{}},"component_identities":{"hdWGCNA":{"analysis_id":"SYNTHETIC_HD_TEST_ONLY","source_fingerprint":"b"*64},"curated_tf_activity":{"analysis_id":"SYNTHETIC_TF_TEST_ONLY","source_fingerprint":"c"*64}}}}
    payload["methods"]=[{"id":"hdwgcna_"+lineage+"_"+branch,"family":"hdWGCNA","status":"complete" if branch=="pooled" else "not_available"} for lineage in ("CA1","IN") for branch in ("pooled","consensus")]
    payload["methods"] += [{"id":"collectri_ulm","family":"curated_tf_activity","status":"complete"},{"id":"pyscenic_motif","family":"motif_regulon","status":"not_available"}]
    for ep in v.EXPECTED:
        for method in ("hdwgcna_"+v.EXPECTED[ep][0]+"_pooled","collectri_ulm"):
            f=method+"::"+ep
            payload["features"].append({"id":f,"method_id":method,"endpoint_id":ep})
            payload["values"] += [{"method_id":method,"endpoint_id":ep,"feature_id":f,"sample_id":m,"value":float(i),"eligible":True} for i,m in enumerate(v.MICE)]
    manifest={"schema":"satb2-regulatory-combined-transport/v1","analysis_id":payload["analysis_id"],"source_fingerprint":payload["source_fingerprint"],"TF_activity_run":True,"motif_regulon_run":False,
              "components":{"hdWGCNA":{"branches":{l+"_pooled":{"validated_complete":True} for l in ("CA1","IN")}},"curated_tf_activity":{"completed_primary":"curated_mouse_CollecTRI_ULM_and_independent_mouse_associations","blocked_optional":"pySCENIC_motif_cache_missing"}}}
    for component,identity in payload["provenance"]["component_identities"].items():manifest["components"][component].update(identity)
    raw=json.dumps(payload).encode();script=b"window.SATB2_REGULATORY_DATA_GZIP="+json.dumps(base64.b64encode(gzip.compress(raw)).decode()).encode()+b";\n"
    (site/"data.js").write_bytes(script)
    manifest["transport"]={"data.js":{"bytes":len(script),"sha256":hashlib.sha256(script).hexdigest()}}
    (site/"manifest.json").write_text(json.dumps(manifest))
    fig=site/"figures/SYNTHETIC_ONLY.png";fig.parent.mkdir();fig.write_bytes(base64.b64decode("iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mP8/x8AAusB9Wl6bAAAAABJRU5ErkJggg=="))
    gallery={"source_fingerprint":payload["source_fingerprint"],"native_bytes_preserved":True,"raster_resizing_performed":False,"figures":[{"component":"hdWGCNA","files":{"png":{"bytes":fig.stat().st_size,"sha256":v.sha(fig),"href":"figures/SYNTHETIC_ONLY.png","location":"Pages"}}}]}
    (site/"figure_manifest.json").write_text(json.dumps(gallery))
    for name in ("index.html","gallery.html","downloads.html"):(site/name).write_text('<a href="index.html">Synthetic test</a>')
    reject(lambda:v.verify(site),"Synthetic")
    print("PASS Production verifier rejects explicitly synthetic payloads")
    # Only the synthetic-exclusion helper is mocked for the positive transport
    # fixture. Production CLI has no bypass for it or for any scientific gate.
    with patch.object(v,"reject_synthetic",lambda _value:None):
        assert v.verify(site)["status"]=="PASS"
        assert v.verify(site)["network_calls"]==0
        with (site/"data.js").open("ab") as f:f.write(b"\n")
        reject(lambda:v.verify(site),"Bytes/SHA256 changed")
        (site/"data.js").write_bytes(script)
        fig.write_bytes(b"changed synthetic test bytes")
        reject(lambda:v.verify(site),"Bytes/SHA256 changed")
    print("PASS Data.js and local native figure bytes are bound by SHA256 receipts")
    changed=copy.deepcopy(payload);changed["endpoints"][4]["n_nuclei"]=568
    reject(lambda:v.validate_payload(changed,manifest),"Frozen endpoint")
    changed=copy.deepcopy(payload);changed["endpoints"][4]["membership_sha256"]="0"*64
    reject(lambda:v.validate_payload(changed,manifest),"Frozen endpoint")
    changed=copy.deepcopy(payload);changed["samples"][0]["genotype"]="cKO"
    reject(lambda:v.validate_payload(changed,manifest),"genotype roster")
    changed=copy.deepcopy(payload);changed["values"].append(changed["values"][0])
    reject(lambda:v.validate_payload(changed,manifest),"Duplicate independent-mouse")
    changed=copy.deepcopy(payload);changed["values"].pop()
    reject(lambda:v.validate_payload(changed,manifest),"all twelve")
    print("PASS Frozen counts / memberships / mouse roster / duplicate and missing mouse rows rejected")
    changed=copy.deepcopy(payload);changed["methods"][-1]["status"]="complete"
    reject(lambda:v.validate_payload(changed,manifest),"Blocked motif")
    changed=copy.deepcopy(manifest);changed["components"]["hdWGCNA"]["branches"]["IN_pooled"]["validated_complete"]=False
    reject(lambda:v.validate_payload(payload,changed),"Native HD component")
    changed=copy.deepcopy(payload);changed["downloads"]=[{"path":"../tf_activity/missing.tsv"}]
    reject(lambda:v.validate_payload(changed,manifest),"broken local scientific path")
    print("PASS Uncompleted components / false motif claims / active broken payload downloads rejected")

print(json.dumps({"status":"PASS","checks":4,"synthetic_only":True,"scientific_fits":0,"network_calls":0,"published":False},indent=2))
