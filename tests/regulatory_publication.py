#!/usr/bin/env python3
"""Synthetic transport tests only. No real run is opened or published.

The full-build transport fixture explicitly mocks scientific validation; the
production validation gate is tested separately and has no CLI bypass.
"""
import argparse
import base64
import gzip
import hashlib
import importlib.util
import json
import tempfile
import copy
import subprocess
import sys
import zipfile
from pathlib import Path
from unittest.mock import patch

REPO = Path(__file__).resolve().parent.parent
spec = importlib.util.spec_from_file_location("regulatory_publication", REPO / "scripts/build_regulatory_publication.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def expect_error(function, phrase):
    try:
        function()
    except ValueError as e:
        assert phrase in str(e), (phrase, str(e))
    else:
        raise AssertionError("Expected rejection: " + phrase)


with tempfile.TemporaryDirectory(prefix="SATB2_SYNTHETIC_TRANSPORT_ONLY_") as temp:
    root = Path(temp)
    results = root / "SYNTHETIC_NATIVE_RUN/results"
    results.mkdir(parents=True)
    source = root / "SYNTHETIC_IMMUTABLE_SOURCE"
    (source / "scripts").mkdir(parents=True)
    config = root / "SYNTHETIC_CONFIG.json"
    config.write_text("{}")
    expect_error(lambda: module.validate_complete(results, source, config), "COMPLETE-only")
    (results / "FINAL_VALIDATION.json").write_text(json.dumps({"status": "RUNNING"}))
    expect_error(lambda: module.validate_complete(results, source, config), "has not completed")
    print("PASS Complete-only gate rejects missing and running results")
    expect_error(lambda: module.validate_tf_complete(results, source, config), "has not completed")
    print("PASS Independent TF completion gate rejects running / unvalidated results")

    fixture = {"schema": module.SCHEMA, "analysis_id": "SYNTHETIC_TRANSPORT_TEST_ONLY", "source_fingerprint": "SYNTHETIC_NOT_A_SCIENTIFIC_RESULT", "status": "partial", "methods": [], "values": [], "provenance": {"synthetic": True}}
    raw = json.dumps(fixture).encode()
    atlas = results / "atlas"
    atlas.mkdir()
    data = b"window.SATB2_REGULATORY_DATA_GZIP=" + json.dumps(base64.b64encode(gzip.compress(raw)).decode()).encode() + b";\n"
    (atlas / "data.js").write_bytes(data)
    expect_error(lambda: module.decode_data(atlas / "data.js"), "Synthetic fixtures")
    print("PASS Synthetic result payload cannot pass the production publication gate")

    setup = root / "SYNTHETIC_SETUP"
    setup.mkdir()
    wrapper = setup / "run_R.sh"
    wrapper.write_text("#!/bin/sh\n# SYNTHETIC TEST ONLY; no program is run\n")
    (setup / "installed_packages.tsv").write_text("SYNTHETIC\tNO_PACKAGES\n")
    (source / "scripts/package_results.py").write_text("# SYNTHETIC SOURCE; never executed\n")
    (source / "scripts/run_hdwgcna.R").write_text("# SYNTHETIC SOURCE; never executed\n")
    native = results / "hdwgcna/CA1/pooled/figures"
    native.mkdir(parents=True)
    png = base64.b64decode("iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mP8/x8AAusB9Wl6bAAAAABJRU5ErkJggg==")
    (native / "SYNTHETIC_FIGURE.png").write_bytes(png)
    (native / "SYNTHETIC_FIGURE.pdf").write_bytes(b"%PDF-1.4\n% SYNTHETIC TEST ONLY\n" + b"x" * 5000)
    (native / "SYNTHETIC_FIGURE.svg").write_text('<svg xmlns="http://www.w3.org/2000/svg"><title>SYNTHETIC TEST ONLY</title></svg>')
    (native.parent / "SYNTHETIC_TABLE.tsv").write_text("gene\tvalue\nSYNTHETIC\tNA\n")
    (native.parent / "analysis_state.rds").write_bytes(b"SYNTHETIC RDS TEST BYTES; not an R object")
    for lineage in ("CA1", "IN"):
        folder = results / "prepared/lineages" / lineage
        folder.mkdir(parents=True)
        (folder / "counts.mtx.gz").write_bytes(gzip.compress(b"SYNTHETIC MTX TEST BYTES; not a scientific matrix"))
        (folder / "cells.tsv").write_text("SYNTHETIC_CELLS\n")
        (folder / "genes.tsv").write_text("SYNTHETIC_GENES\n")
    (atlas / "manifest.json").write_text(json.dumps({"source_fingerprint": fixture["source_fingerprint"]}))
    final = {"status": "COMPLETE", "analysis_id": fixture["analysis_id"], "source_fingerprint": fixture["source_fingerprint"]}
    (results / "FINAL_VALIDATION.json").write_text(json.dumps(final))
    cfg = {"run_root": str(results.parent), "R_wrapper": {"path": str(wrapper)}}
    config.write_text(json.dumps(cfg))
    files = {str(p.relative_to(results)): module.receipt(p, results) for p in results.rglob("*") if p.is_file()}
    bindings = {"R_wrapper": {"path": str(wrapper), "bytes": wrapper.stat().st_size, "sha256": module.digest(wrapper)}}
    mocked = {"final": final, "config": cfg, "payload": fixture, "manifest": {}, "files": files, "bindings": bindings}
    out = root / "SYNTHETIC_PUBLICATION_STAGE"
    site_target = root / "SYNTHETIC_SITE/regulatory"
    args = argparse.Namespace(results=results, source=source, config=config, outdir=out, site_dir=site_target, setup_dir=setup, tag="atlas-v3.7.0", stage_only=True)
    # Force the unchanged test PDF to overflow the Pages figure allowance. The
    # test PNG still fits and must keep its exact bytes as the thumbnail source.
    base = sum((REPO / "site/regulatory" / name).stat().st_size for name in ("styles.css", "app.js", "README.md")) + (atlas / "data.js").stat().st_size + (atlas / "manifest.json").stat().st_size
    with patch.object(module, "validate_complete", return_value=mocked), patch.object(module, "MAX_NEW_PAGES", base + 1_000_000 + 150):
        report = module.build(args)
        assert module.verify_outdir(out)["status"] == "PASS"
    assert not site_target.exists(), "Stage-only transport fixture must not copy into a site"
    assert report["scientific_fits_performed"] is False and report["GitHub_actions_performed"] is False
    assert report["native_figure_groups"] == 1 and report["archive_members_SHA256_and_CRC_verified"] is True
    staged = out / "site/regulatory"
    gallery = (staged / "gallery.html").read_text()
    assert "SYNTHETIC FIGURE" in gallery
    figure_manifest = json.loads((staged / "figure_manifest.json").read_text())
    f = figure_manifest["figures"][0]["files"]
    assert f["png"]["location"] == "Pages" and (staged / f["png"]["href"]).read_bytes() == png
    assert f["pdf"]["location"] == "required_release_asset"
    assert f["pdf"]["href"].startswith(module.RELEASE + "atlas-v3.7.0/")
    assert not list(staged.rglob("*.rds")) and not list(staged.rglob("*.mtx.gz")) and not list(staged.rglob("*.svg"))
    assets = json.loads((out / "RELEASE_ASSETS.json").read_text())["assets"]
    assert len(assets) == 3 and all(Path(a["local_path"]).is_file() for a in assets)
    standalone = next(Path(a["local_path"]) for a in assets if a["kind"] == "interactive_offline_HTML")
    text = standalone.read_text()
    assert '<script src="' not in text and '<link rel="stylesheet"' not in text
    assert 'href="https://ali-mujahid.github.io/satb2_snRNA/' in text
    assert 'src="data.js' not in text and "SATB2_REGULATORY_DATA_GZIP" in text
    assert (staged / "index.html").read_text().index('src="data.js?') < (staged / "index.html").read_text().index('src="app.js?')
    print("PASS Synthetic transport fixture: ZIP SHA256/CRC, all native formats, exact PNG bytes, overflow release PDF, offline inlining and no site installation")

    large = root / "SYNTHETIC_SPLIT.zip"
    original = bytes(range(256)) * 3
    large.write_bytes(original)
    with patch.object(module, "MAX_RELEASE_ASSET", 200):
        parts = module.split_release_asset(large, out / "release_assets")
    assert len(parts) == 4 and b"".join(p.read_bytes() for p in parts) == original
    assert all(p.stat().st_size <= 200 for p in parts)
    print("PASS Release byte-splitting reconstructs the exact original ZIP bytes")

    # Pure transport merge: synthetic independent components, never model fits.
    endpoint_rows = [{"id": ep, "family": "CA1" if ep.startswith("CA1") else "IN", "source_endpoint_id": "SYNTHETIC_"+ep,
                      "membership_sha256": hashlib.sha256(ep.encode()).hexdigest(), "n_nuclei": n, "label": "SYNTHETIC "+ep}
                     for ep, n in [("CA1_DL",8750),("CA1_SL",11706),("IN_C1",1718),("IN_C6",501),("IN_C8",468)]]
    mice = [{"id": g+str(i), "genotype": g} for g in ("WT", "cKO") for i in range(6)]
    hp = {**copy.deepcopy(fixture), "endpoints": endpoint_rows, "samples": mice,
          "methods": [{"id":"hd_primary","family":"hdWGCNA","status":"complete"}],
          "features": [{"id":"hd_feature","method_id":"hd_primary","label":"Synthetic module"}],
          "effects":[],"nodes":[],"edges":[],"hubs":[],"embeddings":[],"coverage":[],
          "provenance":{"synthetic":True,"TF_activity":"NOT_RUN_AWAITING_METHOD_SELECTION"}}
    hp["values"] = [{"method_id":"hd_primary","endpoint_id":ep["id"],"feature_id":"hd_feature","sample_id":m["id"],"value":float(i),"eligible":True}
                    for ep in endpoint_rows for i,m in enumerate(mice)]
    tp = {**copy.deepcopy(fixture), "analysis_id":"SYNTHETIC_TF_TRANSPORT_TEST_ONLY", "source_fingerprint":"SYNTHETIC_TF_SOURCE", "endpoints":copy.deepcopy(endpoint_rows), "samples":copy.deepcopy(mice),
          "methods":[{"id":"collectri_ulm","family":"curated_tf_activity","status":"complete"},{"id":"pyscenic_motif","family":"motif_regulon","status":"not_available"}],
          "features":[{"id":"tf_model","method_id":"collectri_ulm","label":"SYNTHETIC TF"},{"id":"tf_unavailable","method_id":"collectri_ulm","label":"SYNTHETIC unavailable TF"}],
          "effects":[],"nodes":[],"edges":[],"hubs":[],"embeddings":[],"coverage":[],
          "provenance":{"synthetic":True,"TF_activity":"NATIVE_COLLECTRI_ULM_COMPLETE"},
          "wald_enrichment":[{"TF":"SYNTHETIC","escore":1.2,"gene_set_fdr":.3}],"native_supplements":{"example":"SYNTHETIC_ONLY"}}
    tp["values"] = [{"method_id":"collectri_ulm","endpoint_id":ep["id"],"feature_id":f,"sample_id":m["id"],"value":None if f=="tf_unavailable" else float(i),"eligible":f!="tf_unavailable","missing_reason":"FEWER_THAN_5_TARGETS" if f=="tf_unavailable" else ""}
                    for ep in endpoint_rows for f in ("tf_model","tf_unavailable") for i,m in enumerate(mice)]
    hp_before,tp_before = copy.deepcopy(hp),copy.deepcopy(tp)
    merged = module.merge_payloads(hp,tp)
    assert hp==hp_before and tp==tp_before, "Transport merge must not mutate source components"
    assert merged["values"] == hp["values"]+tp["values"]
    assert merged["wald_enrichment"]==tp["wald_enrichment"]
    assert merged["provenance"]["TF_activity"]=="COMPLETE_CURATED_COLLECTRI_ULM"
    assert merged["provenance"]["components"]["hdWGCNA"]["TF_activity"]=="NOT_RUN_AWAITING_METHOD_SELECTION"
    assert merged["provenance"]["motif_supported_regulons"].startswith("BLOCKED_")
    assert sum(r["value"] is None for r in merged["values"])==60
    changed=copy.deepcopy(tp); changed["endpoints"][0]["membership_sha256"]="0"*64
    expect_error(lambda:module.merge_payloads(hp,changed),"memberships")
    changed=copy.deepcopy(tp); changed["samples"][0]["genotype"]="cKO";changed["samples"][6]["genotype"]="WT"
    expect_error(lambda:module.merge_payloads(hp,changed),"rosters")
    changed=copy.deepcopy(tp);changed["methods"][0]["id"]="hd_primary"
    expect_error(lambda:module.merge_payloads(hp,changed),"methods IDs overlap")
    changed=copy.deepcopy(tp);changed["features"][0]["id"]="hd_feature"
    expect_error(lambda:module.merge_payloads(hp,changed),"features IDs overlap")
    changed=copy.deepcopy(tp);changed["values"].append(changed["values"][0])
    expect_error(lambda:module.merge_payloads(hp,changed),"Duplicate combined")
    changed=copy.deepcopy(tp);changed["values"][0]["value"]=None
    expect_error(lambda:module.merge_payloads(hp,changed),"must be finite or explicitly unavailable")
    changed=copy.deepcopy(tp);changed["values"].pop()
    expect_error(lambda:module.merge_payloads(hp,changed),"all twelve mouse values")
    print("PASS Exact memberships / rosters / distinct source keys; native rows, null eligibility and component-scoped provenance preserved")

    tf_results=root/"SYNTHETIC_TF_RUN/results";tf_source=root/"SYNTHETIC_TF_SOURCE";tf_config=root/"SYNTHETIC_TF_CONFIG.json"
    (tf_results/"tf_activity/CA1_Deep").mkdir(parents=True);(tf_source/"scripts").mkdir(parents=True)
    tf_config.write_text("{}")
    for name in ("run_tf_activity.py","analyze_tf_mouse.R","package_tf_results.py"):
        (tf_source/"scripts"/name).write_text("# SYNTHETIC TEST ONLY; never executed\n")
    (tf_results/"tf_activity/CA1_Deep/SYNTHETIC_TF_FIGURE.png").write_bytes(png)
    (tf_results/"tf_activity/CA1_Deep/SYNTHETIC_TF_FIGURE.pdf").write_bytes(b"%PDF-1.4\n% SYNTHETIC TF FIGURE\n")
    (tf_results/"tf_activity/CA1_Deep/SYNTHETIC_TF_TABLE.tsv").write_text("sample\tvalue\nSYNTHETIC\tNA\n")
    (tf_results/"FINAL_VALIDATION.json").write_text(json.dumps({"status":"SYNTHETIC_MOCKED_COMPLETE_ONLY"}))
    tf_files={str(p.relative_to(tf_results)):module.receipt(p,tf_results) for p in tf_results.rglob("*") if p.is_file()}
    shared_inputs={group:{ep:bindings["R_wrapper"] for ep in module.TF_ENDPOINTS} for group in ("source_pseudobulks","existing_deseq2")}
    tv={"final":{"status":"SYNTHETIC_MOCKED_COMPLETE_ONLY","analysis_id":tp["analysis_id"],"source_fingerprint":tp["source_fingerprint"]},"config":shared_inputs,"payload":tp,"manifest":{},"files":tf_files,"bindings":{"prior":bindings["R_wrapper"]}}
    hv={**mocked,"payload":hp,"config":{**cfg,**shared_inputs}}
    combined_args=argparse.Namespace(**vars(args));combined_args.outdir=root/"SYNTHETIC_COMBINED_STAGE"
    combined_args.tf_results=tf_results;combined_args.tf_source=tf_source;combined_args.tf_config=tf_config
    native_before={str(p):p.read_bytes() for p in (tf_results/"FINAL_VALIDATION.json",tf_source/"scripts/run_tf_activity.py",tf_config)}
    with patch.object(module,"validate_complete",return_value=hv),patch.object(module,"validate_tf_complete",return_value=tv):
        combined_report=module.build(combined_args)
    assert module.verify_outdir(combined_args.outdir)["status"]=="PASS"
    combined_site=combined_args.outdir/"site/regulatory"
    combined_json=json.loads(gzip.decompress(base64.b64decode(json.loads((combined_site/"data.js").read_text().split("=",1)[1].strip().rstrip(";")))))
    assert combined_json["values"]==merged["values"] and combined_json["wald_enrichment"]==tp["wald_enrichment"]
    assert combined_report["native_figure_groups"]==2 and combined_report["TF_activity_run"] is True and combined_report["motif_regulon_run"] is False
    assert "Curated CollecTRI TF activity is complete" in (combined_site/"gallery.html").read_text()
    contents=json.loads((combined_args.outdir/"ARCHIVE_CONTENTS.json").read_text())
    assert any(r["archive_path"]=="tf_component/source/scripts/run_tf_activity.py" for r in contents["files"])
    assert any(r["archive_path"]=="tf_component/results/tf_activity/CA1_Deep/SYNTHETIC_TF_TABLE.tsv" for r in contents["files"])
    assert contents["generated_combined_transport"]["sha256"]==hashlib.sha256(module.encoded(merged)).hexdigest()
    extracted=root/"SYNTHETIC_EXTRACTED"
    with zipfile.ZipFile(combined_report["complete_native_zip"]["path"]) as z:z.extractall(extracted)
    validator=next(extracted.rglob("VALIDATE_EXTRACTED.py"))
    proc=subprocess.run([sys.executable,"-B",str(validator)],capture_output=True,text=True,check=True)
    assert proc.stdout.startswith("PASS ")
    assert all(Path(p).read_bytes()==b for p,b in native_before.items()) and not site_target.exists()
    print("PASS Combined synthetic archive / figure gallery / standalone transport; both native sources immutable and extraction validator passes")

    incomplete_args=argparse.Namespace(**vars(args));incomplete_args.outdir=root/"SYNTHETIC_REJECTED_STAGE";incomplete_args.tf_results=tf_results
    with patch.object(module,"validate_complete",return_value=hv):
        expect_error(lambda:module.build(incomplete_args),"together")
    assert not incomplete_args.outdir.exists()
    print("PASS Incomplete TF argument set is rejected before staging")

    different_inputs=copy.deepcopy(shared_inputs)
    different_inputs["source_pseudobulks"]["IN_C8"]["sha256"]="0"*64
    expect_error(lambda:module.validate_shared_inputs(shared_inputs,different_inputs),"original input bytes disagree")
    print("PASS Components cannot merge different original raw-count or DESeq2 input bytes")

    # Exercise the real TF checksum gate independently of the mocked build.
    gate=root/"SYNTHETIC_TF_RECEIPT_GATE";gate_results=gate/"results";gate_results.mkdir(parents=True)
    marker=gate_results/"tf_activity/COMPLETE.json";marker.parent.mkdir();marker.write_text('{"SYNTHETIC":"TEST_ONLY"}\n')
    gate_cfg=root/"SYNTHETIC_TF_GATE_CONFIG.json";gate_cfg.write_text(json.dumps({"schema":"satb2-tf-activity/v1","run_root":str(gate)}))
    gate_final={"schema":"satb2-tf-regulatory-final-validation/v1","status":"COMPLETE_WITH_BLOCKED_MOTIF",
                "TF_activity_run":True,"motif_regulon_run":False,"activity_UMAP_packaged":False,
                "configuration_sha256":module.digest(gate_cfg),"executed_source_sha256":{str(p.relative_to(tf_source)):module.digest(p) for p in (tf_source/"scripts").iterdir()},
                "native_file_receipts":[module.receipt(marker,gate_results)],"package_file_receipts":[]}
    for k in ("twelve_mice_per_endpoint_TF","unique_mouse_keys","raw_mouse_values_finite","frozen_membership_hashes_verified","native_score_receipts_verified","all_configured_inputs_SHA256_before_after_verified","original_analysis_files_unchanged","no_scientific_fits_in_packager"):
        gate_final[k]=True
    (gate_results/"FINAL_VALIDATION.json").write_text(json.dumps(gate_final))
    marker.write_text('{"SYNTHETIC":"CHANGED_TEST_BYTES"}\n')
    expect_error(lambda:module.validate_tf_complete(gate_results,tf_source,gate_cfg),"bytes / SHA256 differ")
    print("PASS Real TF completion gate rejects changed native receipt bytes before transport merging")

    # Isolate the receipt filesystem contract with a synthetic decoded payload.
    # Only decode_data is mocked: source/config/input/result hashes and the
    # exhaustive result-file comparison execute normally. No study data/models.
    marker.write_text('{"SYNTHETIC":"TEST_ONLY"}\n')
    mouse_marker=gate_results/"tf_mouse/COMPLETE.json";mouse_marker.parent.mkdir()
    mouse_marker.write_text('{"SYNTHETIC":"TEST_ONLY"}\n')
    gate_input=root/"SYNTHETIC_TF_BOUND_INPUT.tsv";gate_input.write_text("SYNTHETIC\tNO_STUDY_DATA\n")
    gate_binding=module.receipt(gate_input)
    gate_config={"schema":"satb2-tf-activity/v1","run_root":str(gate),"samples":gate_binding,"references_mouse_collectri":gate_binding,
                 **{name:{ep:gate_binding for ep in module.TF_ENDPOINTS} for name in module.TF_INPUT_MAPS}}
    gate_cfg.write_text(json.dumps(gate_config))
    gate_payload={**copy.deepcopy(tp),"analysis_id":"SYNTHETIC_FILESYSTEM_GATE_ONLY","source_fingerprint":"SYNTHETIC_FILESYSTEM_SOURCE"}
    gate_raw=module.encoded(gate_payload)
    gate_atlas=gate_results/"atlas";gate_atlas.mkdir()
    (gate_atlas/"data.json").write_bytes(gate_raw)
    (gate_atlas/"data.js").write_text("/* SYNTHETIC DECODE MOCK ONLY; cannot pass decode_data */\n")
    gate_manifest={"analysis_id":gate_payload["analysis_id"],"source_fingerprint":gate_payload["source_fingerprint"],
                   "transport":{name:module.receipt(gate_atlas/name,gate_atlas) for name in ("data.js","data.json")}}
    (gate_atlas/"manifest.json").write_text(json.dumps(gate_manifest))
    gate_final.update({"analysis_id":gate_payload["analysis_id"],"source_fingerprint":gate_payload["source_fingerprint"],
                       "configuration_sha256":module.digest(gate_cfg),
                       "native_file_receipts":[module.receipt(p,gate_results) for p in (marker,mouse_marker)],
                       "package_file_receipts":[module.receipt(p,gate_results) for p in sorted(gate_atlas.iterdir())]})
    (gate_results/"FINAL_VALIDATION.json").write_text(json.dumps(gate_final))
    with patch.object(module,"decode_data",return_value=(gate_payload,gate_raw)):
        baseline=module.validate_tf_complete(gate_results,tf_source,gate_cfg)
        controller=gate_results/".snakemake/log/SYNTHETIC_CONTROLLER.log";controller.parent.mkdir(parents=True)
        controller.write_text("SYNTHETIC CONTROLLER LOG; mutable after packaging\n")
        metadata=gate_results/".snakemake/metadata/SYNTHETIC_RECORD";metadata.parent.mkdir(parents=True)
        metadata.write_text("SYNTHETIC CONTROLLER METADATA\n")
        verified=module.validate_tf_complete(gate_results,tf_source,gate_cfg)
        assert verified["files"]==baseline["files"] and not any(".snakemake" in Path(p).parts for p in verified["files"])
        for relative in ("tf_activity/UNRECEIPTED_SCIENTIFIC.tsv","UNRECEIPTED_ROOT.tsv",".snakemake_backup/UNRECEIPTED.tsv"):
            unexpected=gate_results/relative;unexpected.parent.mkdir(parents=True,exist_ok=True);unexpected.write_text("SYNTHETIC UNRECEIPTED BYTES\n")
            expect_error(lambda:module.validate_tf_complete(gate_results,tf_source,gate_cfg),"not exhaustively covered")
            unexpected.unlink()
        assert module.validate_tf_complete(gate_results,tf_source,gate_cfg)["files"]==baseline["files"]
    # The generic enumerator retains .snakemake for scientific/source callers.
    source_metadata=tf_source/".snakemake/SYNTHETIC_SOURCE_FILE";source_metadata.parent.mkdir()
    source_metadata.write_text("SYNTHETIC SOURCE BYTES\n")
    assert source_metadata.resolve() in module.all_files(tf_source)
    source_metadata.unlink()
    source_link=tf_source/"SYNTHETIC_LINK";source_link.symlink_to(gate_input)
    expect_error(lambda:module.all_files(tf_source),"Do not follow a symlink")
    source_link.unlink()
    print("PASS TF receipt comparison ignores only .snakemake controller metadata; unknown scientific/root/near-name files fail and source enumeration remains strict")

    # A changed staged file must invalidate the output receipt.
    with (staged / "data.js").open("ab") as handle:
        handle.write(b"\n")
    expect_error(lambda: module.verify_outdir(out), "bytes / SHA256 differ")
    print("PASS Changed staged content is rejected by the publication receipt")

print(json.dumps({"status": "PASS", "checks": 12, "synthetic_only": True, "scientific_validation_mocked_only_for_transport_fixture": True, "real_outputs_opened": False, "published": False}, indent=2))
