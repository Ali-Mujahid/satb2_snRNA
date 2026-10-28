#!/usr/bin/env python3
"""Synthetic completed-fit copy contract tests; no scientific models or data."""
from pathlib import Path
import copy
import importlib.util
import json
import shutil
import tempfile

ROOT = Path(__file__).resolve().parent.parent
spec = importlib.util.spec_from_file_location("publication", ROOT/"scripts/build_regulatory_publication.py")
m = importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value,sort_keys=True)+"\n")


def reject(call, phrase):
    try:call()
    except ValueError as error:assert phrase in str(error),(phrase,str(error))
    else:raise AssertionError("Expected refusal: "+phrase)


with tempfile.TemporaryDirectory(prefix="SYNTHETIC_RECOVERY_PUBLICATION_ONLY_") as temp:
    t=Path(temp).resolve();old_run=t/"original";old=old_run/"results";new=t/"recovery/results";source=t/"recovery_source"
    original_cfg=t/"original_config.json";cfg_path=t/"recovery_config.json"
    original={"run_root":str(old_run),"workflow_version":"1.0.1","scientific_parameters":{"SYNTHETIC":"UNCHANGED"}}
    write(original_cfg,original);native_record=m.receipt(original_cfg)
    r=source/"scripts/run_hdwgcna.R";r.parent.mkdir(parents=True);r.write_text("# SYNTHETIC TEST ONLY; never executed\n")
    old_r=old_run/"workflow/scripts/run_hdwgcna.R";old_r.parent.mkdir(parents=True);shutil.copy2(r,old_r)
    write(old/"prepared/source_hashes.json",{"configuration":{k:native_record["sha256"] for k in ("expected_sha256","before_sha256","after_sha256")}})
    write(old/"prepared/COMPLETE.json",{"status":"COMPLETE","synthetic":True})
    write(old/"hdwgcna/run_parameters.json",{"script_md5":m.digest(r,"md5")})
    for lineage in ("CA1","IN"):
        for network in ("pooled","consensus"):
            d=old/"hdwgcna"/lineage/network
            write(d/"status.json",{"status":"NETWORK_COMPLETE","synthetic":True})
            write(d/"run_identity.json",{"fingerprint":"SYNTHETIC_"+lineage+network})
            write(d/"completed.json",{"fingerprint":"SYNTHETIC_"+lineage+network,"files":{p.name:m.digest(p,"md5") for p in d.iterdir()}})
    rows=[m.receipt(p,old) for p in sorted(old.rglob("*")) if p.is_file()]
    inv=t/"inventory.json";write(inv,{"schema":"satb2-hdwgcna-copy-inventory/v1","status":"FROZEN_COMPLETE_FIT_COPY_SOURCE","source_results":str(old),"native_config":native_record,"files":rows})
    cfg={**original,"run_root":str(new.parent),"workflow_version":"1.0.2","reuse_completed_native":{"schema":"satb2-hdwgcna-copy-recovery/v1","source_results":str(old),"native_config":native_record,"inventory":m.receipt(inv),"no_network_refits":True}}
    write(cfg_path,cfg)
    new.mkdir(parents=True)
    for name in ("prepared","hdwgcna"):shutil.copytree(old/name,new/name)
    proof={"schema":"satb2-hdwgcna-reused-native-receipt/v1","status":"VERIFIED_EXACT_COPY","scientific_fits_run":False,"source_results":str(old),"native_config":native_record,"inventory":m.receipt(inv),"recovery_config":m.receipt(cfg_path),"copied_files":rows}
    flags=("all_four_completed_checkpoints_MD5_verified","prepared_complete_SHA256_verified","all_source_and_destination_SHA256_verified","scientific_configuration_unchanged","native_fitting_script_MD5_unchanged")
    proof.update({k:True for k in flags});proof_path=new/"REUSED_NATIVE_RECEIPT.json";write(proof_path,proof)
    meta={"schema":"satb2-hdwgcna-native-fit-reuse-validation/v1","source_results":str(old),"native_configuration":native_record,"native_configuration_sha256":native_record["sha256"],"native_fit_configuration":original,"recovery_configuration":m.receipt(cfg_path),"source_inventory":m.receipt(inv),"reuse_receipt":m.receipt(proof_path,new),"copied_file_count":len(rows),"scientific_config_keys_identical":True,"allowed_configuration_differences":["reuse_completed_native","run_root","workflow_version"],"network_refits_run":False,"native_fit_R_script_md5":m.digest(r,"md5"),"preserved_R_script_sha256":m.digest(r)}
    final={"native_fit_reuse":meta,"native_fit_configuration_sha256":native_record["sha256"],"completed_native_fits_reused_without_refit":True}
    payload={"provenance":{"native_fit_reuse":meta}}
    verify=lambda:m.validate_native_reuse(new,source,cfg_path,cfg,final,payload)
    bindings=verify();assert "recovery:native_config" in bindings and len([k for k in bindings if k.startswith("reused_native:")])==len(rows)
    print("PASS Explicit no-refit copy verifies every original/recovered byte, source/config and all four checkpoint seals")
    changed=copy.deepcopy(cfg);changed["scientific_parameters"]["SYNTHETIC"]="CHANGED"
    reject(lambda:m.validate_native_reuse(new,source,cfg_path,changed,final,payload),"scientific configuration")
    print("PASS Recovery cannot change scientific parameters")
    target=new/rows[0]["path"];previous=target.read_bytes();target.write_bytes(b"SYNTHETIC CHANGED COPY\n")
    reject(verify,"bytes / SHA256")
    target.write_bytes(previous)
    print("PASS Changed copied native bytes are rejected")
    for relative in ("hdwgcna/CA1/pooled/SYNTHETIC_UNRECEIPTED.tsv", "prepared/.snakemake/SYNTHETIC_UNRECEIPTED"):
        extra=new/relative;extra.parent.mkdir(parents=True,exist_ok=True);extra.write_text("SYNTHETIC UNRECEIPTED NATIVE BYTES\n")
        reject(verify,"not exhaustive")
        extra.unlink()
    assert verify()
    print("PASS Preserved scientific-tree inventory rejects every unknown file, including reserved names; no global enumeration exemption")
    previous=r.read_bytes();r.write_bytes(previous+b"# CHANGED\n");reject(verify,"fitting source differs");r.write_bytes(previous)
    print("PASS Recovery cannot change original executed fitting source")
    bad=copy.deepcopy(proof);bad["scientific_fits_run"]=True;write(proof_path,bad)
    reject(verify,"verified no-fit receipt");write(proof_path,proof)
    changed=copy.deepcopy(final);changed["native_fit_reuse"]["network_refits_run"]=True
    reject(lambda:m.validate_native_reuse(new,source,cfg_path,cfg,changed,payload),"provenance differs")
    print("PASS New fit claims and inconsistent final/payload provenance are rejected")

print(json.dumps({"status":"PASS","checks":6,"synthetic_only":True,"real_outputs_opened":False,"scientific_fits":0,"published":False},indent=2))
