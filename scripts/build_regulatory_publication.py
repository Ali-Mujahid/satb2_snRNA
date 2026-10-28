#!/usr/bin/env python3
"""Stage COMPLETE native regulatory results for Pages/offline/release transport.

Standard-library packaging only: no scientific subprocess, fitting, upload,
Git command, tag, release or deployment is performed. Do not invoke on real
outputs before FINAL_VALIDATION.json records complete validated primary networks.
"""
from __future__ import annotations

import argparse
import base64
import gzip
import hashlib
import html
import json
import math
import os
import re
import shutil
import zipfile
from pathlib import Path
from urllib.parse import quote, urljoin

ROOT = Path(__file__).resolve().parent.parent
TAG = "atlas-v3.8.0"
RELEASE = "https://github.com/Ali-Mujahid/satb2_snRNA/releases/download/"
LIVE = "https://ali-mujahid.github.io/satb2_snRNA/"
MAX_NEW_PAGES = 38_000_000
MAX_SITE = 950_000_000
MAX_RELEASE_ASSET = 1_900_000_000
SCHEMA = "satb2-regulatory-atlas/v1"
ENDPOINT_IDS = {"CA1_DL", "CA1_SL", "IN_C1", "IN_C6", "IN_C8"}
TF_INPUT_MAPS = ("source_pseudobulks", "fitted_size_factors", "existing_deseq2", "existing_deseq2_stable_id", "endpoint_mouse_metadata")
TF_ENDPOINTS = {"CA1_Deep", "CA1_Superficial", "IN_C1", "IN_C6", "IN_C8"}
PYSCENIC_METHODS = {"pyscenic_motif_CA1", "pyscenic_motif_IN"}
PYSCENIC_SOURCE = ("scripts/run_pyscenic_native.py", "scripts/analyze_regulon_mouse.R", "scripts/package_pyscenic_results.py")
PYSCENIC_NATIVE_CHECKS = ("positive_motif_regulons_only", "full_genotype_neutral_RNA_background", "no_DEG_selection", "native_nucleus_AUCell", "exact_nucleus_axes_verified", "all12_mouse_roster", "min_measured_targets_10", "Satb2_self_target_sensitivity_separate", "inputs_config_source_SHA256_unchanged")
PYSCENIC_MOUSE_CHECKS = ("all12_mouse_roster", "native_AUCell_means_preserved", "global_primary_BH_both_lineages", "separate_sensitivity_BH", "no_activity_threshold_filter", "no_nucleus_pseudoreplication", "no_new_gene_DESeq2", "positive_motif_regulons_only", "inputs_config_source_SHA256_unchanged")
_cache = {}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def digest(path, algorithm="sha256", force=False):
    path = Path(path).resolve()
    st = path.stat()
    key = (str(path), algorithm, st.st_size, st.st_mtime_ns)
    if not force and key in _cache:
        return _cache[key]
    h = hashlib.new(algorithm)
    with path.open("rb") as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(block)
    _cache[key] = h.hexdigest()
    return _cache[key]


def encoded(value):
    return (json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False, sort_keys=True) + "\n").encode()


def load(path):
    return json.loads(Path(path).read_text())


def receipt(path, relative_to=None):
    path = Path(path)
    return {"path": str(path.relative_to(relative_to)) if relative_to else str(path.resolve()), "bytes": path.stat().st_size, "sha256": digest(path)}


def child(root, relative):
    root = Path(root).resolve()
    path = (root / relative).resolve()
    require(path.is_relative_to(root), f"Path escapes source root: {relative}")
    return path


def matching(path, expected):
    path = Path(path)
    require(path.is_file(), f"Required validated file is missing: {path}")
    require(path.stat().st_size == expected["bytes"] and digest(path) == expected["sha256"], f"Validated bytes / SHA256 differ: {path}")


def decode_data(path):
    text = Path(path).read_text()
    match = re.fullmatch(r"\s*(?:/\*[\s\S]*?\*/\s*)?window\.SATB2_REGULATORY_DATA_GZIP\s*=\s*(\"[A-Za-z0-9+/=]*\")\s*;\s*", text)
    require(match, "Expected the generated compressed regulatory data assignment")
    raw = gzip.decompress(base64.b64decode(json.loads(match.group(1)), validate=True))
    value = json.loads(raw)
    require(value.get("schema") == SCHEMA, "Unsupported completed regulatory payload schema")
    require(not value.get("provenance", {}).get("synthetic") and "SYNTHETIC" not in value.get("analysis_id", ""), "Synthetic fixtures must never be published")
    return value, raw


def bound_inputs(config):
    values = {}
    for collection in ("inputs", "existing_deseq2", "source_pseudobulks"):
        require(isinstance(config.get(collection), dict) and config[collection], f"Missing bound input collection: {collection}")
        values.update({collection + ":" + k: v for k, v in config[collection].items()})
    values.update({key: config[key] for key in ("strict_lr_pooled", "runtime_manifest", "R_wrapper")})
    return values


def validate_complete(results, source, config_path):
    results, source, config_path = Path(results).resolve(), Path(source).resolve(), Path(config_path).resolve()
    final_path = results / "FINAL_VALIDATION.json"
    require(final_path.is_file(), "Publication is COMPLETE-only: FINAL_VALIDATION.json is not present")
    final = load(final_path)
    require(final.get("schema") == "satb2-regulatory-final-validation/v1" and final.get("status") in ("COMPLETE", "COMPLETE_PRIMARY_WITH_OPTIONAL_LIMITATIONS"), "Primary scientific validation has not completed")
    for lineage in ("CA1", "IN"):
        require(final.get("primary_networks", {}).get(lineage, {}).get("validated_complete") is True, f"Primary {lineage} network is not validated complete")
    for key in ("twelve_mice_per_endpoint_module", "unique_mouse_keys", "raw_mouse_values_finite", "frozen_membership_hashes_verified", "native_checkpoint_MD5_verified", "prepared_SHA256_receipts_verified", "all_configured_inputs_SHA256_before_after_verified", "original_analysis_files_unchanged", "no_scientific_fits_in_packager"):
        require(final.get(key) is True, f"Required native validation did not pass: {key}")
    require(final.get("TF_activity_run") is False and final.get("motif_regulon_run") is False, "This transport contract covers hdWGCNA; TF/regulon inference requires its own reviewed publication contract")
    config = load(config_path)
    require(results == Path(config["run_root"]).resolve() / "results", "Results are not the immutable configured run")
    require(digest(config_path) == final.get("configuration_sha256"), "The immutable source configuration has changed")
    expected_files = {}
    for r in final.get("native_file_receipts", []) + final.get("package_file_receipts", []):
        p = child(results, r["path"])
        matching(p, r)
        expected_files[r["path"]] = r
    require(expected_files, "Complete validation has no native/package file receipts")
    for lineage in ("CA1", "IN"):
        for name in ("counts.mtx.gz", "cells.tsv", "genes.tsv"):
            require(f"prepared/lineages/{lineage}/{name}" in expected_files, f"Missing selected raw matrix/axis receipt: {lineage}/{name}")
        require(f"hdwgcna/{lineage}/pooled/analysis_state.rds" in expected_files, f"Missing recoverable native {lineage} scientific state")
    atlas = results / "atlas"
    manifest = load(atlas / "manifest.json")
    for name, r in manifest.get("transport", {}).items():
        matching(child(atlas, name), r)
    payload, raw = decode_data(atlas / "data.js")
    require(raw == (atlas / "data.json").read_bytes(), "Compressed UI transport differs from exact saved JSON")
    require(payload.get("analysis_id") == final.get("analysis_id") == manifest.get("analysis_id"), "Analysis identities disagree")
    require(payload.get("source_fingerprint") == final.get("source_fingerprint") == manifest.get("source_fingerprint"), "Source fingerprints disagree")
    require(payload.get("status") in ("complete", "partial") and payload.get("values"), "No completed primary mouse values are available")
    require(not payload.get("embeddings"), "Activity UMAP transport is unavailable in this hdWGCNA-only publication contract")
    require(all(m.get("family") == "hdWGCNA" or m.get("status") not in ("complete", "completed", "partial") for m in payload.get("methods", [])), "An undeclared TF/regulon result cannot be relabeled as a completed method")
    packager = source / "scripts/package_results.py"
    require(packager.is_file() and digest(packager) == payload.get("provenance", {}).get("packager_sha256"), "Supplied scientific source differs from the executed immutable packager")
    r_params = load(results / "hdwgcna/run_parameters.json")
    require(digest(source / "scripts/run_hdwgcna.R", "md5") == r_params.get("script_md5"), "Supplied R source differs from the native executed R script")
    bindings = bound_inputs(config)
    for name, item in bindings.items():
        path = Path(item["path"])
        require(path.is_absolute() and path.is_file() and digest(path) == item["sha256"], f"Original bound source changed or is unavailable: {name}")
        require(item.get("bytes") is None or path.stat().st_size == item["bytes"], f"Original bound source size changed: {name}")
    if config.get("reuse_completed_native") is not None:
        require("REUSED_NATIVE_RECEIPT.json" in expected_files, "Recovered fit-copy proof is absent from final native receipts")
        bindings.update(validate_native_reuse(results, source, config_path, config, final, payload))
    else:
        require(not final.get("native_fit_reuse") and not payload.get("provenance", {}).get("native_fit_reuse"), "Native reuse cannot be claimed without an explicit configuration/proof")
    # The final receipt itself is also included byte-identically in the archive.
    expected_files["FINAL_VALIDATION.json"] = receipt(final_path, results)
    return {"final": final, "config": config, "payload": payload, "manifest": manifest, "files": expected_files, "bindings": bindings}


def validate_native_reuse(results, source, config_path, cfg, final, payload):
    """Independently verify an explicit copy-only export recovery, never refit."""
    results, source, config_path = map(lambda p: Path(p).resolve(), (results, source, config_path))
    reuse = cfg.get("reuse_completed_native", {})
    require(reuse.get("schema") == "satb2-hdwgcna-copy-recovery/v1" and reuse.get("no_network_refits") is True, "Unsupported or non-copy native recovery")
    extra = {}
    for name in ("native_config", "inventory"):
        row = reuse.get(name, {})
        p = Path(row.get("path", ""))
        require(p.is_absolute() and not p.is_symlink() and re.fullmatch(r"[0-9a-f]{64}", row.get("sha256", "")), f"Unbound recovery source: {name}")
        matching(p, row)
        extra["recovery:" + name] = row
    native_record, inv_record = reuse["native_config"], reuse["inventory"]
    original = load(native_record["path"])
    allowed = {"run_root", "workflow_version", "reuse_completed_native"}
    require({k: v for k, v in cfg.items() if k not in allowed} == {k: v for k, v in original.items() if k not in allowed}, "Export recovery changed scientific configuration keys")
    old = Path(reuse.get("source_results", ""))
    require(old.is_absolute() and old.is_dir() and not old.is_symlink(), "Original native recovery source is unavailable")
    old = old.resolve()
    require(old == Path(original["run_root"]).resolve() / "results" and old != results and not old.is_relative_to(results) and not results.is_relative_to(old), "Original and recovered native results must be distinct nonoverlapping trees")
    inv = load(inv_record["path"])
    require(inv.get("schema") == "satb2-hdwgcna-copy-inventory/v1" and inv.get("status") == "FROZEN_COMPLETE_FIT_COPY_SOURCE", "Recovery inventory is not a frozen completed-fit source")
    require(Path(inv.get("source_results", "")).resolve() == old and inv.get("native_config") == native_record, "Recovery inventory source/config identity differs")
    def records(rows):
        require(isinstance(rows, list) and rows, "Empty completed-native copy inventory")
        out = {}
        for row in rows:
            rel = row.get("path", "")
            parts = Path(rel).parts
            require(isinstance(rel, str) and parts and parts[0] in ("prepared", "hdwgcna") and not Path(rel).is_absolute() and ".." not in parts and "\\" not in rel, "Copy inventory escapes prepared/hdwgcna")
            require(rel not in out and isinstance(row.get("bytes"), int) and row["bytes"] >= 0 and re.fullmatch(r"[0-9a-f]{64}", row.get("sha256", "")), "Duplicate or invalid native copy receipt")
            out[rel] = {k: row[k] for k in ("path", "bytes", "sha256")}
        return out
    expected = records(inv.get("files"))
    proof_path = results / "REUSED_NATIVE_RECEIPT.json"
    proof = load(proof_path)
    require(proof.get("schema") == "satb2-hdwgcna-reused-native-receipt/v1" and proof.get("status") == "VERIFIED_EXACT_COPY" and proof.get("scientific_fits_run") is False, "Recovered fit copy lacks an explicit verified no-fit receipt")
    require(proof.get("native_config") == native_record and proof.get("inventory") == inv_record and Path(proof.get("source_results", "")).resolve() == old and proof.get("recovery_config") == receipt(config_path), "Native reuse proof source/config binding differs")
    require(records(proof.get("copied_files")) == expected, "Native reuse proof does not cover the full frozen inventory")
    for flag in ("all_four_completed_checkpoints_MD5_verified", "prepared_complete_SHA256_verified", "all_source_and_destination_SHA256_verified", "scientific_configuration_unchanged", "native_fitting_script_MD5_unchanged"):
        require(proof.get(flag) is True, f"Native recovery verification did not pass: {flag}")
    for root in (old, results):
        actual = set()
        for name in ("prepared", "hdwgcna"):
            tree = root / name
            require(tree.is_dir() and not tree.is_symlink(), f"Missing regular preserved native tree: {tree}")
            actual.update(str(p.relative_to(root)) for p in all_files(tree))
        require(actual == set(expected), f"Preserved native inventory is not exhaustive: {root}")
        for rel, row in expected.items():
            matching(child(root, rel), row)
    for rel, row in expected.items():
        extra["reused_native:" + rel] = {**row, "path": str(child(old, rel))}
    cfg_hashes = load(results / "prepared/source_hashes.json").get("configuration", {})
    require(all(cfg_hashes.get(k) == native_record["sha256"] for k in ("expected_sha256", "before_sha256", "after_sha256")), "Prepared counts were not bound to the preserved original fitting configuration")
    for lineage in ("CA1", "IN"):
        for network in ("pooled", "consensus"):
            folder = results / "hdwgcna" / lineage / network
            c = load(folder / "completed.json")
            require(load(folder / "status.json").get("status") == "NETWORK_COMPLETE" and c.get("fingerprint") == load(folder / "run_identity.json").get("fingerprint"), "Copied network is not a completed original checkpoint")
            for rel, expected_md5 in c.get("files", {}).items():
                require(digest(child(folder, rel), "md5") == expected_md5, "Copied native checkpoint MD5 differs")
    r_script = source / "scripts/run_hdwgcna.R"
    old_r_script = Path(original["run_root"]) / "workflow/scripts/run_hdwgcna.R"
    require(old_r_script.is_file() and not old_r_script.is_symlink() and digest(old_r_script) == digest(r_script), "Recovery fitting source differs from the preserved original executed source")
    extra["recovery:native_fitting_source"] = receipt(old_r_script)
    md5 = load(results / "hdwgcna/run_parameters.json").get("script_md5")
    require(digest(r_script, "md5") == md5, "Recovery fitting R source differs from the native executed MD5")
    reuse_meta = final.get("native_fit_reuse", {})
    require(reuse_meta == payload.get("provenance", {}).get("native_fit_reuse") and reuse_meta.get("schema") == "satb2-hdwgcna-native-fit-reuse-validation/v1", "Native reuse provenance differs between final receipt and payload")
    checks = {"source_results": str(old), "native_configuration": native_record, "native_configuration_sha256": native_record["sha256"], "native_fit_configuration": original, "recovery_configuration": receipt(config_path), "source_inventory": inv_record, "reuse_receipt": receipt(proof_path, results), "copied_file_count": len(expected), "scientific_config_keys_identical": True, "allowed_configuration_differences": sorted(allowed), "network_refits_run": False, "native_fit_R_script_md5": md5, "preserved_R_script_sha256": digest(r_script)}
    require(all(reuse_meta.get(k) == value for k, value in checks.items()) and final.get("native_fit_configuration_sha256") == native_record["sha256"] and final.get("completed_native_fits_reused_without_refit") is True, "Native reuse scientific/source identity metadata is inconsistent")
    return extra


def tf_bound_inputs(config):
    """Find explicit file receipts, never guess scientific dependencies by path."""
    for name in ("samples", "references_mouse_collectri"):
        require(isinstance(config.get(name), dict) and {"path", "bytes", "sha256"} <= set(config[name]), f"Missing TF bound input: {name}")
    for name in TF_INPUT_MAPS:
        require(isinstance(config.get(name), dict) and set(config[name]) == TF_ENDPOINTS, f"TF inputs must bind exactly five endpoints: {name}")
    values = {}
    def visit(obj, prefix=""):
        if isinstance(obj, dict):
            if "path" in obj and "sha256" in obj:
                require("bytes" in obj, f"TF input lacks a byte receipt: {prefix}")
                values[prefix] = obj
            else:
                for key, value in obj.items():
                    visit(value, prefix + ":" + key if prefix else key)
    visit(config)
    require(values, "No explicit TF bound input receipts")
    return values


def validate_tf_complete(results, source, config_path):
    """Verify the independent TF native completion contract; perform no analysis."""
    results, source, config_path = map(lambda p: Path(p).resolve(), (results, source, config_path))
    final_path = results / "FINAL_VALIDATION.json"
    require(final_path.is_file(), "TF publication is COMPLETE-only: FINAL_VALIDATION.json is not present")
    final = load(final_path)
    require(final.get("schema") == "satb2-tf-regulatory-final-validation/v1" and final.get("status") == "COMPLETE_WITH_BLOCKED_MOTIF", "TF scientific validation has not completed under the reviewed contract")
    for key in ("twelve_mice_per_endpoint_TF", "unique_mouse_keys", "raw_mouse_values_finite", "frozen_membership_hashes_verified", "native_score_receipts_verified", "all_configured_inputs_SHA256_before_after_verified", "original_analysis_files_unchanged", "no_scientific_fits_in_packager"):
        require(final.get(key) is True, f"Required TF validation did not pass: {key}")
    require(final.get("TF_activity_run") is True and final.get("motif_regulon_run") is False and final.get("activity_UMAP_packaged") is False, "TF curated activity and blocked motif/UMAP statuses are not explicit")
    config = load(config_path)
    require(config.get("schema") == "satb2-tf-activity/v1" and results == Path(config["run_root"]).resolve() / "results", "TF results are not the immutable configured run")
    require(digest(config_path) == final.get("configuration_sha256"), "The immutable TF configuration has changed")
    source_hashes = final.get("executed_source_sha256", {})
    for name in ("scripts/run_tf_activity.py", "scripts/analyze_tf_mouse.R", "scripts/package_tf_results.py"):
        require(name in source_hashes, f"Missing executed TF scientific source receipt: {name}")
    for rel, sha in source_hashes.items():
        require(re.fullmatch(r"[0-9a-f]{64}", sha) and digest(child(source, rel)) == sha, f"Supplied TF scientific source differs from execution: {rel}")
    files = {}
    for row in final.get("native_file_receipts", []) + final.get("package_file_receipts", []):
        require(row["path"] not in files, f"Repeated TF native/package receipt: {row['path']}")
        matching(child(results, row["path"]), row)
        files[row["path"]] = row
    require(files and "tf_activity/COMPLETE.json" in files and "tf_mouse/COMPLETE.json" in files, "TF native completion markers are not covered by receipts")
    # Snakemake writes mutable controller logs/metadata after native packaging.
    # Exempt only that reserved directory from this result-receipt comparison;
    # leave all_files strict for executed source and preserved scientific trees.
    actual = {str(p.relative_to(results)) for p in all_files(results)
              if ".snakemake" not in p.relative_to(results).parts} - {"FINAL_VALIDATION.json"}
    require(actual == set(files), "TF result files are not exhaustively covered by final receipts")
    atlas = results / "atlas"
    manifest = load(atlas / "manifest.json")
    require(manifest.get("transport"), "TF atlas transport has no checksums")
    for name, row in manifest["transport"].items():
        matching(child(atlas, name), row)
    payload, raw = decode_data(atlas / "data.js")
    require(raw == (atlas / "data.json").read_bytes(), "TF compressed transport differs from the native saved JSON")
    require(payload.get("analysis_id") == final.get("analysis_id") == manifest.get("analysis_id") and payload.get("source_fingerprint") == final.get("source_fingerprint") == manifest.get("source_fingerprint"), "TF component identities disagree")
    require(payload.get("status") in ("complete", "partial") and payload.get("values") and not payload.get("embeddings"), "TF component lacks completed mouse activity or declares an unsupported activity UMAP")
    completed = [m for m in payload.get("methods", []) if m.get("status") in ("complete", "completed", "partial")]
    require(completed and all(m.get("family") == "curated_tf_activity" for m in completed), "Only curated-prior TF activity may be declared complete by this TF contract")
    motif = [m for m in payload.get("methods", []) if m.get("family") == "motif_regulon"]
    require(motif and all(m.get("status") == "not_available" for m in motif), "Blocked motif-supported inference must be declared unavailable")
    bindings = tf_bound_inputs(config)
    for name, row in bindings.items():
        p = Path(row["path"])
        require(p.is_absolute() and not p.is_symlink(), f"TF bound input is not an absolute regular source: {name}")
        matching(p, row)
    files["FINAL_VALIDATION.json"] = receipt(final_path, results)
    return {"final": final, "config": config, "payload": payload, "manifest": manifest, "files": files, "bindings": bindings}


def pyscenic_bound_inputs(config):
    """Explicit recursive native input records, including ranking lists."""
    require(config.get("schema") == "satb2-pyscenic/v1", "Unexpected pySCENIC configuration schema")
    require({"samples", "prepared_complete", "lineages", "references", "additional_hash_inputs"} <= set(config), "Missing pySCENIC bound inputs")
    require(set(config["lineages"]) == {"CA1", "IN"}, "pySCENIC needs exact CA1/IN lineage inputs")
    for lineage in ("CA1", "IN"):
        require({"counts", "genes", "cells"} <= set(config["lineages"][lineage]), "pySCENIC lineage inputs are incomplete")
    require({"tf_names", "motif_annotations", "rankings"} <= set(config["references"]) and isinstance(config["references"]["rankings"], list) and len(config["references"]["rankings"]) == 2, "pySCENIC requires the two matched ranking references")
    require({"frozen_membership_audit", "mouse_coverage"} <= set(config["additional_hash_inputs"]), "pySCENIC frozen cohort evidence is missing")
    found = {}
    def walk(value, key=""):
        if isinstance(value, dict):
            if "path" in value or "sha256" in value:
                require({"path", "bytes", "sha256"} <= set(value), f"Incomplete pySCENIC input receipt: {key}")
                found[key] = value
            else:
                for name, child_value in value.items():
                    walk(child_value, key + ":" + name if key else name)
        elif isinstance(value, list):
            for i, child_value in enumerate(value):
                walk(child_value, key + ":" + str(i))
    walk(config)
    require(found, "No explicit pySCENIC input bindings")
    return found


def validate_pyscenic_complete(results, source, config_path):
    """Validate a separately completed positive native motif/AUCell component."""
    results, source, config_path = map(lambda p: Path(p).resolve(), (results, source, config_path))
    final_path = results / "FINAL_VALIDATION.json"
    require(final_path.is_file(), "pySCENIC publication is COMPLETE-only: FINAL_VALIDATION.json is not present")
    final = load(final_path)
    require(final.get("schema") == "satb2-pyscenic-regulatory-final-validation/v1" and final.get("status") == "COMPLETE", "pySCENIC native scientific validation has not completed")
    for flag in ("TF_activity_run", "motif_regulon_run", "positive_motif_regulons_only", "twelve_mice_per_endpoint_regulon", "unique_mouse_keys", "raw_mouse_values_finite", "frozen_membership_hashes_verified", "native_score_receipts_verified", "exact_nucleus_axes_verified", "Satb2_self_target_sensitivity_separate", "all_configured_inputs_SHA256_before_after_verified", "original_analysis_files_unchanged", "no_scientific_fits_in_packager"):
        require(final.get(flag) is True, f"Required pySCENIC validation did not pass: {flag}")
    require(final.get("minimum_measured_targets") == 10 and final.get("activity_UMAP_packaged") is False, "pySCENIC target eligibility / unavailable UMAP contract differs")
    require(isinstance(final.get("Satb2_official_TF_list_membership"), bool), "pySCENIC official Satb2 TF-list membership must be explicitly audited")
    config = load(config_path)
    require(results == Path(config["run_root"]).resolve() / "results" and digest(config_path) == final.get("configuration_sha256"), "Immutable pySCENIC run configuration differs")
    source_hashes = final.get("executed_source_sha256", {})
    require(set(source_hashes) == set(PYSCENIC_SOURCE), "pySCENIC executed scientific source receipts differ")
    for rel, sha in source_hashes.items():
        require(re.fullmatch(r"[0-9a-f]{64}", sha) and digest(child(source, rel)) == sha, f"Supplied pySCENIC source differs from execution: {rel}")
    files = {}
    for row in final.get("native_file_receipts", []) + final.get("package_file_receipts", []):
        require(row["path"] not in files, "Repeated pySCENIC native/package receipt")
        matching(child(results, row["path"]), row)
        files[row["path"]] = row
    require({"pyscenic_native/COMPLETE.json", "regulon_mouse/COMPLETE.json"} <= set(files), "Native pySCENIC and mouse completion markers are not sealed")
    for relative, schema, script in (("pyscenic_native/COMPLETE.json", "satb2-pyscenic-native-complete/v1", PYSCENIC_SOURCE[0]), ("regulon_mouse/COMPLETE.json", "satb2-regulon-mouse-complete/v1", PYSCENIC_SOURCE[1])):
        marker = load(child(results, relative))
        require(marker.get("schema") == schema and marker.get("status") == "COMPLETE" and marker.get("config", {}).get("sha256") == digest(config_path) and marker.get("source", {}).get("sha256") == source_hashes[script], "Native pySCENIC step seal/configuration/source identity differs")
        required = PYSCENIC_NATIVE_CHECKS if relative.startswith("pyscenic_native/") else PYSCENIC_MOUSE_CHECKS
        require(all(marker.get("checks", {}).get(flag) is True for flag in required), "Native pySCENIC step scientific checks did not all pass")
        folder = child(results, str(Path(relative).parent))
        sealed = marker.get("output_files", marker.get("outputs", {}))
        require(sealed and {str(p.relative_to(folder)) for p in all_files(folder)} == set(sealed) | {"COMPLETE.json"}, "Native pySCENIC step output seal is not exhaustive")
        for name, row in sealed.items():
            matching(child(folder, name), row)
    actual = {str(p.relative_to(results)) for p in all_files(results) if ".snakemake" not in p.relative_to(results).parts} - {"FINAL_VALIDATION.json"}
    require(actual == set(files), "pySCENIC result files are not exhaustively covered by final receipts")
    for lineage, count in (("CA1", 20456), ("IN", 2687)):
        axis = final.get("native_nucleus_axis_audit", {}).get(lineage, {})
        require(axis.get("n_nuclei") == count and axis.get("exact_source_nucleus_axis_verified") is True and axis.get("no_score_recomputation_or_averaging") is True, f"Exact native pySCENIC nucleus axis was not validated: {lineage}")
        for name in ("auc_table", "metadata_table"):
            require("pyscenic_native/" + axis.get(name, "") in files, f"Unsealed native pySCENIC nucleus-axis file: {lineage}/{name}")
    atlas = results / "atlas"
    manifest = load(atlas / "manifest.json")
    require(manifest.get("completed_primary") == "positive_motif_supported_regulons_native_nucleus_AUCell_and_independent_mouse_associations" and manifest.get("transport"), "Native pySCENIC transport completion differs")
    for name, row in manifest["transport"].items():
        matching(child(atlas, name), row)
    payload, raw = decode_data(atlas / "data.js")
    require(raw == (atlas / "data.json").read_bytes(), "pySCENIC compressed transport differs from native JSON")
    require(payload.get("analysis_id") == final.get("analysis_id") == manifest.get("analysis_id") and payload.get("source_fingerprint") == final.get("source_fingerprint") == manifest.get("source_fingerprint"), "pySCENIC component source identities disagree")
    require(payload.get("status") == "complete" and not payload.get("embeddings"), "pySCENIC workflow is incomplete or declares an unvalidated activity UMAP")
    methods = {m.get("id"): m for m in payload.get("methods", [])}
    require(len(payload.get("methods", [])) == 2 and set(methods) == PYSCENIC_METHODS and all(m.get("family") == "motif_regulon" and m.get("status") in ("complete", "not_available") for m in methods.values()), "pySCENIC lineage-specific method identities/statuses differ")
    require(digest(source / "scripts/package_pyscenic_results.py") == payload.get("provenance", {}).get("packager_sha256"), "pySCENIC transport does not match its executed packager")
    bindings = pyscenic_bound_inputs(config)
    for name, row in bindings.items():
        p = Path(row["path"])
        require(p.is_absolute() and not p.is_symlink(), f"pySCENIC bound input is not an absolute regular source: {name}")
        matching(p, row)
    files["FINAL_VALIDATION.json"] = receipt(final_path, results)
    return {"final": final, "config": config, "payload": payload, "manifest": manifest, "files": files, "bindings": bindings}


def validate_pyscenic_shared_inputs(hd, ps):
    """No scientific input substitution while adding the independent component."""
    for lineage in ("CA1", "IN"):
        for key, filename in (("counts", "counts.mtx.gz"), ("genes", "genes.tsv"), ("cells", "cells.tsv")):
            expected = hd["files"].get("prepared/lineages/" + lineage + "/" + filename, {})
            actual = ps["config"]["lineages"][lineage][key]
            require(all(actual.get(k) == expected.get(k) for k in ("bytes", "sha256")), f"pySCENIC selected raw lineage input differs from completed HD preparation: {lineage}/{key}")
    for ps_record, hd_record, label in ((ps["config"]["samples"], hd["config"]["inputs"]["samples"], "samples"), (ps["config"]["prepared_complete"], hd["files"].get("prepared/COMPLETE.json", {}), "prepared_complete")):
        require(all(ps_record.get(k) == hd_record.get(k) for k in ("bytes", "sha256")), f"pySCENIC original shared input bytes disagree: {label}")


def merge_payloads(hd, tf, pyscenic=None):
    """New transport only: immutable native JSON and all scientific rows survive."""
    endpoint_fields = ("id", "family", "source_endpoint_id", "membership_sha256", "n_nuclei")
    def endpoints(payload):
        rows = payload.get("endpoints", [])
        require(len(rows) == 5 and {r.get("id") for r in rows} == ENDPOINT_IDS, "Combined publication requires exactly the same five endpoints")
        require(all(re.fullmatch(r"[0-9a-f]{64}", r.get("membership_sha256", "")) for r in rows), "Endpoint membership hashes are missing or invalid")
        return {r["id"]: {k: r.get(k) for k in endpoint_fields} for r in rows}
    require(endpoints(hd) == endpoints(tf), "HD/TF frozen endpoint memberships, identities or nucleus counts disagree")
    def samples(payload):
        rows = payload.get("samples", [])
        require(len(rows) == 12 and len({r.get("id") for r in rows}) == 12 and all(r.get("id") for r in rows), "Combined publication requires twelve unique independent mice")
        require(sum(r.get("genotype") == "WT" for r in rows) == 6 and sum(r.get("genotype") == "cKO" for r in rows) == 6, "Combined publication requires six WT and six cKO mice")
        return {r["id"]: r["genotype"] for r in rows}
    hd_samples, tf_samples = samples(hd), samples(tf)
    require(hd_samples == tf_samples, "HD/TF mouse rosters or genotypes disagree")
    def unique(rows, label):
        keys = [r.get("id") for r in rows]
        require(all(isinstance(k, str) and k for k in keys) and len(set(keys)) == len(keys), f"Invalid/duplicate {label} IDs")
        return set(keys)
    for field in ("methods", "features", "nodes", "embeddings"):
        require(not unique(hd.get(field, []), "HD " + field) & unique(tf.get(field, []), "TF " + field), f"HD/TF {field} IDs overlap")
    seen = set()
    for component in (hd, tf):
        method_ids = {r["id"] for r in component.get("methods", [])}
        features = {r["id"]: r for r in component.get("features", [])}
        for row in component.get("values", []):
            key = tuple(row.get(k) for k in ("method_id", "endpoint_id", "feature_id", "sample_id"))
            require(key not in seen, "Duplicate combined independent-mouse value key")
            seen.add(key)
            require(key[0] in method_ids and key[1] in ENDPOINT_IDS and key[2] in features and key[3] in hd_samples, "Combined mouse row has an absent component source key")
            feature = features[key[2]]
            require(feature["method_id"] == key[0] and (not feature.get("endpoint_id") or feature["endpoint_id"] == key[1]), "Combined mouse row crosses feature source memberships")
            require(isinstance(row.get("eligible"), bool), "Combined mouse eligibility must be explicit")
            finite = isinstance(row.get("value"), (int, float)) and not isinstance(row["value"], bool) and math.isfinite(row["value"])
            missing = row.get("value") is None and row["eligible"] is False and bool(row.get("missing_reason"))
            require(finite or missing, "Combined raw mouse scores must be finite or explicitly unavailable; missing is not zero")
    value_rosters = {}
    for row in tf.get("values", []):
        value_rosters.setdefault(tuple(row[k] for k in ("method_id", "endpoint_id", "feature_id")), set()).add(row["sample_id"])
    require(value_rosters and all(ids == set(tf_samples) for ids in value_rosters.values()), "A TF feature/endpoint does not retain all twelve mouse values")
    # JSON copying avoids aliasing even for nested provenance/settings dictionaries.
    merged = json.loads(json.dumps(hd, allow_nan=False))
    merged["analysis_id"] = hd["analysis_id"] + "__" + tf["analysis_id"]
    material = {"schema": "satb2-regulatory-combined-source/v1", "components": {"hdWGCNA": {"analysis_id": hd["analysis_id"], "source_fingerprint": hd["source_fingerprint"]}, "curated_tf_activity": {"analysis_id": tf["analysis_id"], "source_fingerprint": tf["source_fingerprint"]}}, "endpoint_memberships": endpoints(hd), "samples": hd_samples, "publication_builder_sha256": digest(__file__)}
    merged["source_fingerprint"] = hashlib.sha256(encoded(material)).hexdigest()
    merged["status"] = "partial" if any(p.get("status") == "partial" for p in (hd, tf)) else "complete"
    for field in ("methods", "features", "values", "effects", "hubs", "nodes", "edges", "embeddings", "coverage"):
        merged[field] = json.loads(json.dumps(hd.get(field, []) + tf.get(field, []), allow_nan=False))
    merged["wald_enrichment"] = json.loads(json.dumps(hd.get("wald_enrichment", []) + tf.get("wald_enrichment", []), allow_nan=False))
    merged["native_supplements"] = json.loads(json.dumps(tf.get("native_supplements", {}), allow_nan=False))
    merged["component_supplements"] = {"hdWGCNA": hd.get("native_supplements", {}), "curated_tf_activity": tf.get("native_supplements", {})}
    merged["downloads"] = [{"component": name, "native_component_record": row} for name, component in (("hdWGCNA", hd), ("curated_tf_activity", tf)) for row in component.get("downloads", [])]
    merged["provenance"] = {"publication_transport": "verified_complete_components_merged_without_scientific_refit", "publication_builder_sha256": material["publication_builder_sha256"], "component_identities": material["components"], "historical_provenance_scope": "Each preserved component provenance describes only its own executed workflow; historical hdWGCNA TF_NOT_RUN flags do not describe the combined release.", "components": {"hdWGCNA": hd.get("provenance", {}), "curated_tf_activity": tf.get("provenance", {})}, "TF_activity": "COMPLETE_CURATED_COLLECTRI_ULM", "motif_supported_regulons": "BLOCKED_OFFICIAL_MOUSE_CISTARGET_RESOURCES_UNAVAILABLE_NOT_INFERRED", "activity_UMAP": "NOT_PACKAGED_NO_VERIFIED_COORDINATE_PLUS_SCORE_AXIS", "no_scientific_fits_in_publication_builder": True, "original_component_bytes_unchanged": True, "interpretation": "Coexpression modules, curated-prior TF activity, and their independent-mouse associations are distinct exploratory outputs. Neither curated activity nor signed prior edges establish motif-supported regulons, direct binding or causality."}
    return merge_pyscenic_payload(merged, pyscenic) if pyscenic is not None else merged


def merge_pyscenic_payload(combined, ps):
    """Add actual motif results; historical component payloads stay untouched."""
    fields = ("id", "family", "source_endpoint_id", "membership_sha256", "n_nuclei")
    def endpoints(p):
        rows = p.get("endpoints", [])
        require(len(rows) == 5 and {r.get("id") for r in rows} == ENDPOINT_IDS, "pySCENIC must retain the same five frozen endpoint identities")
        return {r["id"]: {k: r.get(k) for k in fields} for r in rows}
    require(endpoints(combined) == endpoints(ps), "pySCENIC frozen endpoint memberships or counts disagree")
    roster = {r["id"]: r.get("genotype") for r in combined["samples"]}
    require(len(ps.get("samples", [])) == 12 and len({r.get("id") for r in ps["samples"]}) == 12 and {r["id"]: r.get("genotype") for r in ps["samples"]} == roster, "pySCENIC frozen mouse roster disagrees")
    for name in ("methods", "features", "nodes", "embeddings"):
        ids = [r.get("id") for r in ps.get(name, [])]
        require(all(isinstance(i, str) and i for i in ids) and len(set(ids)) == len(ids) and not set(ids) & {r["id"] for r in combined.get(name, [])}, f"Duplicate/overlapping pySCENIC {name} IDs")
    require({m.get("id") for m in ps.get("methods", [])} == PYSCENIC_METHODS and all(m.get("family") == "motif_regulon" and m.get("status") in ("complete", "not_available") for m in ps["methods"]), "Unverified motif method cannot replace historical unavailable placeholder")
    features = {r["id"]: r for r in ps.get("features", [])}
    methods = {r["id"]: r for r in ps["methods"]}
    for feature in features.values():
        method = feature.get("method_id", "");lineage = method.rsplit("_", 1)[-1]
        require(method in methods and feature.get("native_name", "").endswith("(+)") and feature.get("source_id", "").startswith(lineage + "::"), "Motif features must retain positive native lineage-specific definitions")
    groups = {}
    seen = set()
    for row in ps.get("values", []):
        k = tuple(row.get(n) for n in ("method_id", "endpoint_id", "feature_id", "sample_id"))
        require(k not in seen and k[0] in methods and k[1] in ENDPOINT_IDS and k[2] in features and k[3] in roster, "Duplicate or absent pySCENIC independent-mouse source key")
        seen.add(k);feature = features[k[2]]
        require(feature["method_id"] == k[0] and k[1] in methods[k[0]].get("endpoint_ids", []), "pySCENIC mouse score crosses native lineage definitions")
        value = row.get("value");finite = isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)
        missing = value is None and row.get("eligible") is False and bool(row.get("missing_reason"))
        require(isinstance(row.get("eligible"), bool) and (finite and 0 <= value <= 1 or missing), "pySCENIC AUCell must retain finite [0,1] AUC or explicit unavailable status")
        require(methods[k[0]]["status"] != "not_available" or missing, "Unavailable pySCENIC lineage cannot carry an active score")
        groups.setdefault(k[:3], set()).add(k[3])
    require(all(ids == set(roster) for ids in groups.values()), "A pySCENIC regulon/endpoint lacks all twelve mouse rows")
    expected_groups = {(f["method_id"], ep, f["id"]) for f in features.values() for ep in methods[f["method_id"]].get("endpoint_ids", [])}
    require(set(groups) == expected_groups, "A pySCENIC feature is missing an entire matching endpoint/mouse group")
    for method, meta in methods.items():
        has_activity = any(r["method_id"] == method and r.get("eligible") is True and r.get("value") is not None for r in ps.get("values", []))
        require((meta["status"] == "complete") == has_activity, "pySCENIC lineage activity status does not match eligible native values")
    result = json.loads(json.dumps(combined, allow_nan=False))
    placeholder = [m for m in result["methods"] if m.get("id") == "pyscenic_motif"]
    require(len(placeholder) == 1 and placeholder[0].get("family") == "motif_regulon" and placeholder[0].get("status") == "not_available", "Historical motif placeholder is not safely replaceable")
    require(not any(r.get("method_id") == "pyscenic_motif" for name in ("features", "values", "effects", "hubs", "nodes", "edges", "coverage", "embeddings") for r in result.get(name, [])), "Historical motif placeholder has scientific rows; cannot discard them")
    result["methods"] = [m for m in result["methods"] if m.get("id") != "pyscenic_motif"]
    for name in ("methods", "features", "values", "effects", "hubs", "nodes", "edges", "coverage", "embeddings"):
        result[name].extend(json.loads(json.dumps(ps.get(name, []), allow_nan=False)))
    result["analysis_id"] += "__" + ps["analysis_id"]
    provenance = result["provenance"]
    identity = {"analysis_id": ps["analysis_id"], "source_fingerprint": ps["source_fingerprint"]}
    provenance["component_identities"]["motif_regulon"] = identity
    provenance["components"]["motif_regulon"] = json.loads(json.dumps(ps.get("provenance", {})))
    provenance["motif_supported_regulons"] = "COMPLETE_NATIVE_POSITIVE_PYSCENIC"
    provenance["historical_provenance_scope"] = "Each preserved native component records its own executed workflow and historical resource state. Earlier motif-unavailable flags remain unchanged inside the archived HD/curated-TF components. Only the top-level empty placeholder is superseded by the independently validated motif component."
    provenance["interpretation"] = "Coexpression, curated-prior TF activity and native positive motif-supported regulons are distinct exploratory outputs. Motif enrichment and expression-derived candidate orientation do not establish direct binding, physical interaction or causality."
    result["component_supplements"]["motif_regulon"] = json.loads(json.dumps(ps.get("native_supplements", {})))
    if "Satb2_status" in ps:
        result["Satb2_status"] = json.loads(json.dumps(ps["Satb2_status"], allow_nan=False))
    result["downloads"].extend({"component": "motif_regulon", "native_component_record": row} for row in ps.get("downloads", []))
    result["status"] = "complete" if all(m.get("status") in ("complete", "completed", "partial") for m in result["methods"]) else "partial"
    result["source_fingerprint"] = hashlib.sha256(encoded({"schema": "satb2-regulatory-three-component-source/v1", "components": provenance["component_identities"], "memberships": endpoints(result), "samples": roster, "publication_builder_sha256": digest(__file__)})).hexdigest()
    return result


def validate_shared_inputs(hd_config, tf_config):
    """Bind both methods to the same original endpoint counts and gene fits."""
    for collection in ("source_pseudobulks", "existing_deseq2"):
        hd = hd_config.get(collection, {})
        tf = tf_config.get(collection, {})
        require(set(hd) == set(tf) == TF_ENDPOINTS, f"HD/TF original source endpoint collections differ: {collection}")
        for endpoint in TF_ENDPOINTS:
            require(all(hd[endpoint].get(k) == tf[endpoint].get(k) for k in ("bytes", "sha256")), f"HD/TF original input bytes disagree: {collection}/{endpoint}")


def write_combined_transport(folder, payload, hd_validation, tf_validation, pyscenic_validation=None):
    raw = encoded(payload)
    packed = gzip.compress(raw, compresslevel=9, mtime=0)
    script = b"/* Verified completed components; display transport merge only. */\nwindow.SATB2_REGULATORY_DATA_GZIP=" + json.dumps(base64.b64encode(packed).decode()).encode() + b";\n"
    data = Path(folder) / "data.js"
    require(not data.exists(), "Refusing to replace a staged transport")
    data.write_bytes(script)
    manifest = {"schema": "satb2-regulatory-combined-transport/v1", "analysis_id": payload["analysis_id"], "source_fingerprint": payload["source_fingerprint"], "status": payload["status"], "components": {"hdWGCNA": hd_validation["manifest"], "curated_tf_activity": tf_validation["manifest"]}, "record_counts": {k: len(payload[k]) for k in ("methods", "endpoints", "samples", "features", "values", "effects", "hubs", "nodes", "edges", "embeddings", "coverage")}, "transport": {"data.js": {k: v for k, v in receipt(data).items() if k != "path"}}, "scientific_fits_in_publication_builder": False, "native_component_bytes_unchanged": True, "TF_activity_run": True, "motif_regulon_run": False}
    if pyscenic_validation is not None:
        manifest["components"]["motif_regulon"] = pyscenic_validation["manifest"]
        manifest["motif_regulon_run"] = True
        final = pyscenic_validation["final"]
        manifest["component_validation"] = {"motif_regulon": {k: final[k] for k in ("schema", "status", "positive_motif_regulons_only", "twelve_mice_per_endpoint_regulon", "exact_nucleus_axes_verified", "Satb2_self_target_sensitivity_separate", "Satb2_official_TF_list_membership", "minimum_measured_targets", "activity_UMAP_packaged")}}
        manifest["component_validation"]["motif_regulon"]["FINAL_VALIDATION_sha256"] = pyscenic_validation["files"]["FINAL_VALIDATION.json"]["sha256"]
    (Path(folder) / "manifest.json").write_bytes(encoded(manifest))
    return manifest


def copy_exact(source, target):
    target = Path(target)
    target.parent.mkdir(parents=True, exist_ok=True)
    require(not target.exists(), f"Refusing to overwrite staged content: {target}")
    shutil.copy2(source, target)
    require(digest(source) == digest(target), f"Copy checksum mismatch: {target}")


def release_name(relative):
    relative = str(relative)
    stem = re.sub(r"[^A-Za-z0-9_.-]+", "__", relative)
    return "regulatory__" + hashlib.sha256(relative.encode()).hexdigest()[:10] + "__" + stem


def all_files(folder):
    folder = Path(folder).resolve()
    result = []
    for path in sorted(folder.rglob("*")):
        require(not path.is_symlink(), f"Do not follow a symlink while publishing: {path}")
        if path.is_file() and not any(p in (".git", "__pycache__") for p in path.relative_to(folder).parts) and path.name not in (".DS_Store",) and path.suffix != ".pyc":
            result.append(path)
    return result


def make_archive(path, results, source, config_path, validation, setup_dir, tf_component=None, combined_payload=None, pyscenic_component=None):
    top = "SATB2_regulatory_" + validation["final"]["analysis_id"]
    members = []
    for rel, r in sorted(validation["files"].items()):
        members.append({"source": str(child(results, rel)), "archive_path": "results/" + rel, "role": "validated_native_result", **{k: r[k] for k in ("bytes", "sha256")}})
    for p in all_files(source):
        members.append({"source": str(p), "archive_path": "source/" + str(p.relative_to(source)), "role": "exact_scientific_source", **receipt(p)})
    # receipt(path) has a local path field; archive_path remains the portable key.
    members.append({"source": str(config_path), "archive_path": "config/" + config_path.name, "role": "immutable_configuration", **receipt(config_path)})
    setup_files = []
    if setup_dir.is_dir():
        for p in sorted(setup_dir.iterdir()):
            if p.is_file() and not p.is_symlink() and p.suffix.lower() in (".r", ".py", ".sh", ".md", ".json", ".tsv", ".txt", ".lock", ".toml", ".yaml", ".yml"):
                setup_files.append(p)
                members.append({"source": str(p), "archive_path": "setup/" + p.name, "role": "runtime_setup_or_manifest", **receipt(p)})
    for name, item in validation["bindings"].items():
        if name.startswith(("existing_deseq2:", "source_pseudobulks:", "recovery:")) or name in ("strict_lr_pooled", "runtime_manifest", "R_wrapper"):
            p = Path(item["path"])
            members.append({"source": str(p), "archive_path": "bound_reference_inputs/" + name.replace(":", "/") + "/" + p.name, "role": "unchanged_bound_reference_input", **receipt(p)})
    if tf_component:
        tv = tf_component["validation"]
        p = Path(__file__).resolve()
        members.append({"source": str(p), "archive_path": "publication_source/" + p.name, "role": "exact_transport_merge_source_no_analysis", **receipt(p)})
        for rel, row in sorted(tv["files"].items()):
            p = child(tf_component["results"], rel)
            members.append({"source": str(p), "archive_path": "tf_component/results/" + rel, "role": "validated_native_TF_result", **{k: row[k] for k in ("bytes", "sha256")}})
        for p in all_files(tf_component["source"]):
            members.append({"source": str(p), "archive_path": "tf_component/source/" + str(p.relative_to(tf_component["source"])), "role": "exact_TF_scientific_source", **receipt(p)})
        p = tf_component["config"]
        members.append({"source": str(p), "archive_path": "tf_component/config/" + p.name, "role": "immutable_TF_configuration", **receipt(p)})
        for name, item in tv["bindings"].items():
            p = Path(item["path"])
            members.append({"source": str(p), "archive_path": "tf_component/bound_reference_inputs/" + name.replace(":", "/") + "/" + p.name, "role": "unchanged_TF_bound_reference_input", **receipt(p)})
    ps_input_locations = []
    if pyscenic_component:
        pv = pyscenic_component["validation"]
        for rel, row in sorted(pv["files"].items()):
            p = child(pyscenic_component["results"], rel)
            members.append({"source": str(p), "archive_path": "pyscenic_component/results/" + rel, "role": "validated_native_motif_regulon_result", **{k: row[k] for k in ("bytes", "sha256")}})
        for p in all_files(pyscenic_component["source"]):
            members.append({"source": str(p), "archive_path": "pyscenic_component/source/" + str(p.relative_to(pyscenic_component["source"])), "role": "exact_motif_regulon_scientific_source", **receipt(p)})
        p = pyscenic_component["config"]
        members.append({"source": str(p), "archive_path": "pyscenic_component/config/" + p.name, "role": "immutable_motif_regulon_configuration", **receipt(p)})
        # Keep a single byte-identical copy of shared selected RNA inputs. Large
        # public ranking/annotation files stay external, hash/URL-bound, rather
        # than duplicating ~589 MB on a space-limited research workstation.
        known = {(r["bytes"], r["sha256"]): r["archive_path"] for r in members}
        reference_metadata = {}
        for item in pv["bindings"].values():
            if Path(item["path"]).name == "READY_MANIFEST.json":
                record = load(item["path"])
                for row in record.get("resources", []):
                    reference_metadata[row.get("sha256")] = {k: row[k] for k in ("requested_url", "final_url", "head_content_length_bytes", "official_checksums_verified", "tls_validation") if k in row}
        for name, item in sorted(pv["bindings"].items()):
            p = Path(item["path"]);row = {"input_key": name, **receipt(p)}
            ident = (row["bytes"], row["sha256"])
            if name.startswith("references:rankings:") or name == "references:motif_annotations":
                details = reference_metadata.get(row["sha256"], {})
                urls = [item.get(k) for k in ("source_url", "url", "requested_url") if item.get(k)]
                if not details and urls:
                    details = {"requested_url": urls[0]}
                require(str(details.get("requested_url", "")).startswith("https://resources.aertslab.org/cistarget/"), "Large external pySCENIC reference lacks hash-bound official download provenance")
                row.update(location="external_official_reference", included_in_archive=False, download_provenance=details, reason="Public reference retained as exact URL/length/SHA256 binding; native results, source, configuration and reference metadata are included")
            elif ident in known:
                row.update(location="shared_byte_identical_archive_member", archive_path=known[ident], included_in_archive=True)
            else:
                archive_path = "pyscenic_component/bound_reference_inputs/" + name.replace(":", "/") + "/" + p.name
                members.append({"source": str(p), "archive_path": archive_path, "role": "unchanged_motif_regulon_bound_input_or_reference_metadata", **receipt(p)})
                known[ident] = archive_path
                row.update(location="native_bound_input_archive_member", archive_path=archive_path, included_in_archive=True)
            ps_input_locations.append(row)
    names = [r["archive_path"] for r in members]
    require(len(names) == len(set(names)), "Archive path collision")
    portable_members = [{k: v for k, v in r.items() if k not in ("source", "path")} for r in members]
    identity = combined_payload or validation["payload"]
    contents = {"schema": "satb2-regulatory-native-archive/v1", "analysis_id": identity["analysis_id"], "source_fingerprint": identity["source_fingerprint"], "files": portable_members, "scientific_bytes_preserved": True, "scientific_fits_performed": False, "installed_runtime_packages_included": False, "selected_raw_MTX_and_axes_included": True, "native_RDS_and_all_figures_included": True, "external_original_inputs": validation["bindings"], "TF_component_included": bool(tf_component), "motif_regulon_run": False, "note": "The selected raw lineage matrices/axes reproduce the downstream analysis. Original all-nuclei H5AD and annotation inputs remain external hash-bound inputs; installed package libraries are excluded."}
    if tf_component:
        contents["components"] = {"hdWGCNA": {"analysis_id": validation["final"]["analysis_id"], "source_fingerprint": validation["final"]["source_fingerprint"], "results_prefix": "results/", "source_prefix": "source/"}, "curated_tf_activity": {"analysis_id": tf_component["validation"]["final"]["analysis_id"], "source_fingerprint": tf_component["validation"]["final"]["source_fingerprint"], "results_prefix": "tf_component/results/", "source_prefix": "tf_component/source/"}}
        contents["TF_bound_inputs"] = tf_component["validation"]["bindings"]
    if pyscenic_component:
        pv = pyscenic_component["validation"]
        contents["motif_regulon_run"] = True
        contents["pySCENIC_component_included"] = True
        contents["components"]["motif_regulon"] = {"analysis_id": pv["final"]["analysis_id"], "source_fingerprint": pv["final"]["source_fingerprint"], "results_prefix": "pyscenic_component/results/", "source_prefix": "pyscenic_component/source/"}
        contents["pySCENIC_bound_input_locations"] = ps_input_locations
        contents["large_external_reference_files_included"] = False
        contents["large_external_reference_policy"] = "Two official mm10 ranking files and matching MGI motif annotations remain external, exact HTTPS URL/length/SHA256-bound; included READY/schema metadata records upstream integrity evidence. Shared selected raw input matrices appear once in existing exact native archive members."
    readme = ("SATB2 completed hdWGCNA analysis\n\nThis archive preserves exact validated result bytes, selected raw matrices/axes, recoverable native RDS state, all PDF/PNG/SVG figures, exact scientific source, immutable configuration, runtime setup records, and unchanged original DESeq2/LR reference tables.\n\nValidate extracted files with: python3 VALIDATE_EXTRACTED.py\n\nAnalysis scope: exploratory undirected RNA coexpression and independent-mouse eigengene associations. No TF activity or motif-supported regulon analysis is claimed. The configuration retains original absolute input paths; review source/README.md and native provenance before reproducing the workflow. No model is run by the validation helper.\n").encode()
    if tf_component:
        readme = ("SATB2 completed hdWGCNA and curated CollecTRI TF activity\n\nThe hdWGCNA component is preserved under results/, source/ and config/. The independently validated TF component is preserved under tf_component/results/, tf_component/source/ and tf_component/config/. Every receipt-covered native file, figure, exact scientific source/configuration and TF bound reference input is included byte-identically. Historical TF_NOT_RUN statements inside the hdWGCNA component describe that workflow only.\n\nValidate extracted files with: python3 VALIDATE_EXTRACTED.py\n\nCurated signed prior edges are not de novo inferred motif-supported regulons or demonstrated direct binding. Full pySCENIC remains blocked on unavailable official mouse cisTarget resources. Mouse-score associations are exploratory and technical partition partly confounds genotype. This transport build performs no scientific model. Installed runtime libraries and the original all-nuclei H5AD remain external dependencies.\n").encode()
    if pyscenic_component:
        readme = ("SATB2 completed hdWGCNA, curated CollecTRI TF activity and native positive pySCENIC motif regulons\n\nThree independently completed native components are preserved under results/, tf_component/results/ and pyscenic_component/results/, with their exact source and immutable configurations alongside. Every receipt-covered native result, per-nucleus AUCell table, mouse score/association table, target definition, figure, source and configuration is included byte-identically. Historical resource-blocked statements inside earlier native components remain unchanged and describe those components at execution time. The combined transport supersedes only the empty top-level unavailable motif placeholder.\n\nValidate extracted files with: python3 VALIDATE_EXTRACTED.py\n\nARCHIVE_CONTENTS.json lists all native member hashes and pySCENIC_bound_input_locations. Shared selected raw matrices are included once and referenced by their exact archive paths. The two large official mm10 ranking files and MGI annotation table remain external with exact official HTTPS download URLs, expected bytes and SHA256, plus included hash-bound reference metadata. Installed software libraries and original full all-nuclei H5AD/annotation inputs also remain external dependencies. Download external references and verify their hashes before reproducing the native workflow.\n\nMotif enrichment, native positive target membership and candidate TF-to-target arrows do not establish direct binding, physical contacts or causality. Independent-mouse associations remain exploratory because genotype partly confounds technical partition. No scientific model is run by publication or validation helpers.\n").encode()
        lines = ["\nVerified reference availability / compatibility metadata included in this archive:"]
        for row in ps_input_locations:
            if Path(row["path"]).name in {"READY_MANIFEST.json", "SCHEMA_INSPECTION.json"}:
                lines.extend(["  " + row["archive_path"], "  Bytes: " + str(row["bytes"]), "  SHA256: " + row["sha256"]])
        lines.append("\nExternal public motif reference files (not included in this ZIP; needed to reproduce inference, not to use the saved offline explorer):")
        for row in ps_input_locations:
            if row["location"] == "external_official_reference":
                lines.extend(["\n  " + Path(row["path"]).name, "  Official HTTPS download: " + row["download_provenance"]["requested_url"], "  Expected exact bytes: " + str(row["bytes"]), "  Expected SHA256: " + row["sha256"]])
        readme += ("\n".join(lines) + "\n").encode()
    helper = b'''#!/usr/bin/env python3
"""Validate extracted native result bytes only; never execute analysis."""
import hashlib,json
from pathlib import Path
root=Path(__file__).resolve().parent
manifest=json.loads((root/'ARCHIVE_CONTENTS.json').read_text())
for r in manifest['files']:
 p=(root/r['archive_path']).resolve()
 assert p.is_relative_to(root) and p.is_file(),r['archive_path']
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 assert p.stat().st_size==r['bytes'] and h.hexdigest()==r['sha256'],r['archive_path']
if manifest.get('generated_combined_transport'):
 r=manifest['generated_combined_transport']; p=(root/r['archive_path']).resolve()
 assert p.is_relative_to(root) and p.is_file(),r['archive_path']
 assert p.stat().st_size==r['bytes'] and hashlib.sha256(p.read_bytes()).hexdigest()==r['sha256'],r['archive_path']
print('PASS',len(manifest['files']),'native/source/config/setup files; no analysis executed')
'''
    generated = {"ARCHIVE_CONTENTS.json": encoded(contents), "README.txt": readme, "VALIDATE_EXTRACTED.py": helper}
    if combined_payload:
        generated["combined_transport/data.json"] = encoded(combined_payload)
        contents["generated_combined_transport"] = {"archive_path": "combined_transport/data.json", "bytes": len(generated["combined_transport/data.json"]), "sha256": hashlib.sha256(generated["combined_transport/data.json"]).hexdigest(), "operation": "merge_verified_transport_no_scientific_fit"}
        generated["ARCHIVE_CONTENTS.json"] = encoded(contents)
    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6, allowZip64=True) as archive:
        for r in members:
            archive.write(r["source"], top + "/" + r["archive_path"])
        for name, content in generated.items():
            archive.writestr(top + "/" + name, content)
    # Stream every member back out: this verifies both ZIP CRC and native SHA256.
    expected = {top + "/" + r["archive_path"]: r for r in members}
    expected.update({top + "/" + n: {"bytes": len(v), "sha256": hashlib.sha256(v).hexdigest()} for n, v in generated.items()})
    with zipfile.ZipFile(path) as archive:
        require(set(archive.namelist()) == set(expected) and len(archive.namelist()) == len(expected), "Archive membership mismatch")
        for name, r in expected.items():
            h, n = hashlib.sha256(), 0
            with archive.open(name) as f:
                for block in iter(lambda: f.read(8 * 1024 * 1024), b""):
                    h.update(block)
                    n += len(block)
            require(n == r["bytes"] and h.hexdigest() == r["sha256"], f"Archive native-byte validation failed: {name}")
    return contents, members


def split_release_asset(path, asset_dir):
    if path.stat().st_size < MAX_RELEASE_ASSET:
        return [path]
    parts = []
    with path.open("rb") as inp:
        i = 1
        while True:
            first = inp.read(min(8 * 1024 * 1024, MAX_RELEASE_ASSET))
            if not first:
                break
            part = asset_dir / (path.name + f".part{i:03d}")
            with part.open("wb") as out:
                out.write(first)
                remaining = MAX_RELEASE_ASSET - len(first)
                while remaining:
                    block = inp.read(min(8 * 1024 * 1024, remaining))
                    if not block:
                        break
                    out.write(block)
                    remaining -= len(block)
            parts.append(part)
            i += 1
    return parts


def portable(index, target, origin):
    sources = []
    def js(match):
        p = child(index.parent, match.group(1).split("?")[0]) if not match.group(1).startswith("../") else (index.parent / match.group(1).split("?")[0]).resolve()
        require(p.is_file(), f"Portable script is unavailable: {p}")
        sources.append(receipt(p))
        return "<script>\n" + p.read_text().replace("</script", "<\\/script") + "\n</script>"
    def css(match):
        p = child(index.parent, match.group(1).split("?")[0])
        sources.append(receipt(p))
        return "<style>\n" + p.read_text() + "\n</style>"
    text = re.sub(r'<script src="([^"]+)"></script>', js, index.read_text())
    text = re.sub(r'<link rel="stylesheet" href="([^"]+)">', css, text)
    text = re.sub(r'href="([^"]+)"', lambda m: m.group(0) if m.group(1).startswith(("https:", "http:", "#", "data:", "mailto:")) else 'href="' + urljoin(origin, m.group(1)) + '"', text)
    require(not re.search(r'<script\s+src=|<link\b[^>]*rel="stylesheet"', text), "Portable document still has external rendering dependencies")
    require("SATB2_REGULATORY_DATA_GZIP" in text and "gunzipSync" in text and "jsPDF" in text, "Portable renderer is incomplete")
    target.write_text(text)
    return {"file": target.name, "bytes": target.stat().st_size, "sha256": digest(target), "embedded_sources": sources, "external_render_dependencies": [], "scientific_fits_performed": False}


def gallery(records, archive_url, offline_url, status, tf_complete=False, pyscenic_complete=False):
    cards = []
    for item in records:
        files = item["files"]
        png = files.get("png")
        thumb = '<a href="' + html.escape(png["href"], quote=True) + '"><img loading="lazy" src="' + html.escape(png["href"], quote=True) + '" alt="' + html.escape(item["label"], quote=True) + '"></a>' if png else '<p class="no-thumb">Native PDF / SVG figure; no PNG thumbnail was produced.</p>'
        links = " ".join('<a href="' + html.escape(r["href"], quote=True) + '" download>' + ext.upper() + '</a>' for ext, r in files.items() if ext in ("png", "pdf"))
        if "svg" in files:
            links += ' <a href="' + archive_url + '">SVG in complete ZIP</a>'
        cards.append('<article data-lineage="' + item["lineage"] + '" data-network="' + item["network"] + '"><h2>' + html.escape(item["label"]) + '</h2>' + thumb + '<div class="links">' + links + '</div><details><summary>Exact source and hashes</summary><pre>' + html.escape(json.dumps(item, indent=2)) + '</pre></details></article>')
    script = """const f=()=>{const q=document.querySelector('#q').value.toLowerCase(),l=document.querySelector('#lineage').value,n=document.querySelector('#network').value;let shown=0;document.querySelectorAll('article').forEach(a=>{const yes=(!q||a.textContent.toLowerCase().includes(q))&&(!l||a.dataset.lineage===l)&&(!n||a.dataset.network===n);a.hidden=!yes;if(yes)shown++;});document.querySelector('#count').textContent=shown+' native figure groups shown';};['q','lineage','network'].forEach(id=>document.querySelector('#'+id).addEventListener('input',f));f();"""
    page = ('<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>SATB2 · Native regulatory figures</title><style>body{font:15px/1.5 system-ui,sans-serif;background:#f3f7f9;color:#263e4d;margin:0}header,main{padding:24px;max-width:1600px;margin:auto}nav,.filters,.links{display:flex;flex-wrap:wrap;gap:12px}a{color:#216a89}input,select{font:inherit;padding:9px;border:1px solid #b8ccd8;border-radius:6px}.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(330px,1fr));gap:18px}article{border:1px solid #d3e0e8;border-radius:10px;background:white;padding:15px}h2{font-size:16px}img{width:100%;height:250px;object-fit:contain}pre{white-space:pre-wrap;overflow-wrap:anywhere;font-size:11px}summary{cursor:pointer;margin-top:12px}.no-thumb{padding:50px 15px;color:#617683}[hidden]{display:none!important}</style></head><body><header><h1>Native hdWGCNA figures</h1><p>Completed scientific outputs · ' + html.escape(status) + '. Every figure is copied without changing its native bytes. Technical partition partly confounds genotype; mouse-level module associations are exploratory.</p><nav><a href="index.html">Interactive regulatory explorer</a><a href="' + archive_url + '">Complete results / source ZIP</a><a href="' + offline_url + '">Standalone offline explorer</a><a href="../communication/index.html">Existing DEG-only LR views</a></nav><p>All figures, SVG originals, native RDS state, selected raw matrices and tables are in the complete ZIP. Some figure downloads use release assets to keep Pages within its size limit. CSV exports are available in the interactive explorer. No TF activity or motif-supported regulon inference is claimed.</p><div class="filters"><input id="q" type="search" placeholder="Find a module, plot or source path" aria-label="Search figures"><select id="lineage" aria-label="Lineage"><option value="">All lineages</option><option>CA1</option><option>IN</option><option value="other">Other / combined</option></select><select id="network" aria-label="Network"><option value="">All networks</option><option>pooled</option><option>consensus</option><option value="other">Other / combined</option></select></div><p id="count"></p></header><main><div class="grid">' + "".join(cards) + '</div></main><script>' + script + '</script></body></html>')
    if tf_complete:
        page = page.replace('Native hdWGCNA figures', 'Native hdWGCNA and curated TF activity figures').replace('No TF activity or motif-supported regulon inference is claimed.', 'Curated CollecTRI TF activity is complete. Full pySCENIC motif-supported inference remains blocked on unavailable official mouse cisTarget resources; curated prior edges are not inferred motif regulons.').replace('<option>consensus</option>', '<option>consensus</option><option value="curated_tf_activity">Curated TF activity</option>')
    if pyscenic_complete:
        page = page.replace('Native hdWGCNA and curated TF activity figures', 'Native coexpression, TF activity and motif regulon figures').replace('Curated CollecTRI TF activity is complete. Full pySCENIC motif-supported inference remains blocked on unavailable official mouse cisTarget resources; curated prior edges are not inferred motif regulons.', 'Curated CollecTRI TF activity and native positive pySCENIC motif-supported regulon discovery/AUCell are complete as separate methods. Motif enrichment and candidate target arrows do not establish direct binding or causality. Historical resource limitations remain preserved inside earlier native components.').replace('<option value="curated_tf_activity">Curated TF activity</option>', '<option value="curated_tf_activity">Curated TF activity</option><option value="motif_regulon">Native positive motif regulons</option>')
    return page.encode()


def downloads_page(zip_assets, release, offline_url, zip_sha, tables, tf_complete=False, pyscenic_complete=False):
    parts = len(zip_assets) > 1
    archive_links = "".join('<li><a href="' + release + quote(p.name) + '">' + html.escape(p.name) + '</a></li>' for p in zip_assets)
    explanation = '<p><b>Download every archive part before extracting.</b> The parts are byte segments of one ZIP. Concatenate them in numeric order:</p><pre>cat SATB2_regulatory_complete_native_results_source.zip.part* &gt; SATB2_regulatory_complete_native_results_source.zip</pre>' if parts else '<p>The complete ZIP contains every validated result, selected raw lineage matrices, native RDS state, PDF/PNG/SVG figure, scientific source, configuration and setup records.</p>'
    table_rows = "".join('<tr><td>' + html.escape(r["path"]) + '</td><td>' + str(r["bytes"]) + '</td><td><code>' + r["sha256"] + '</code></td></tr>' for r in tables)
    page = ('<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>SATB2 · Completed regulatory downloads</title><style>body{font:15px/1.5 system-ui,sans-serif;color:#263e4d;max-width:1300px;margin:25px auto;padding:0 20px}table{border-collapse:collapse;width:100%;font-size:12px}td,th{border-bottom:1px solid #d5e0e7;padding:8px;text-align:left}code,pre{overflow-wrap:anywhere;white-space:pre-wrap}a{color:#216a89}.wrap{overflow:auto}</style></head><body><h1>Completed regulatory results and source</h1><p><a href="index.html">Interactive explorer</a> · <a href="gallery.html">Native figure gallery</a> · <a href="' + offline_url + '">Standalone offline HTML</a></p><h2 id="native-archive">Exact native results / source archive</h2>' + explanation + '<ul>' + archive_links + '</ul><p>Complete reconstructed ZIP SHA256: <code>' + zip_sha + '</code>.</p><p>After extraction, run <code>python3 VALIDATE_EXTRACTED.py</code> inside the extracted top folder to verify native/source bytes. The helper performs no scientific analysis. Original full all-nuclei H5AD/annotation inputs and installed runtime packages remain external hash-bound dependencies; the selected raw lineage matrices and metadata are included.</p><h2>Native tables included in the archive</h2><p>These exact TSV/TSV.gz files are in the complete archive, preserving their original schemas and missing values. Use the explorer CSV buttons for filtered current-view values.</p><div class="wrap"><table><thead><tr><th>Path inside results/</th><th>Bytes</th><th>SHA256</th></tr></thead><tbody>' + table_rows + '</tbody></table></div><p>Mouse-level module associations are exploratory continuous-score effects, not gene DEGs. TF activity and motif-supported regulons were not inferred by this workflow.</p></body></html>')
    if tf_complete:
        page = page.replace('Path inside results/', 'Path inside complete ZIP').replace('TF activity and motif-supported regulons were not inferred by this workflow.', 'Curated CollecTRI ULM activity and separate independent-mouse score associations are complete under tf_component/. Full pySCENIC motif-supported regulon inference remains resource-blocked. Historical hdWGCNA TF_NOT_RUN statements are component-scoped and preserved unchanged.')
    if pyscenic_complete:
        page = page.replace('Curated CollecTRI ULM activity and separate independent-mouse score associations are complete under tf_component/. Full pySCENIC motif-supported regulon inference remains resource-blocked. Historical hdWGCNA TF_NOT_RUN statements are component-scoped and preserved unchanged.', 'Curated CollecTRI ULM is preserved under tf_component/. Native positive motif-supported pySCENIC regulons, exact nucleus AUCell axes, independent-mouse means/associations and the separate Satb2 selfgene sensitivity are preserved under pyscenic_component/. Large public cisTarget rankings/MGI annotations remain external official-URL/length/SHA256-bound references; their verified provenance metadata and all scientific code/configuration are included. ARCHIVE_CONTENTS.json lists their download details and shared raw-matrix archive locations. Historical blocked/TF_NOT_RUN component statements remain unchanged. Motif enrichment does not prove direct binding or causality.')
    return page.encode()


def verify_outdir(outdir):
    outdir = Path(outdir).resolve()
    r = load(outdir / "PUBLICATION_RECEIPT.json")
    require(r.get("status") == "PREPARED_COMPLETE_ONLY", "Not a completed publication staging receipt")
    for item in r["staged_files"]:
        matching(child(outdir, item["path"]), item)
    require(r["new_pages_bytes"] <= MAX_NEW_PAGES and r["prospective_site_bytes"] <= MAX_SITE, "Recorded Pages size exceeds hard limits")
    return {"status": "PASS", "files": len(r["staged_files"]), "source_fingerprint": r["source_fingerprint"], "new_pages_bytes": r["new_pages_bytes"], "scientific_fits_performed": False}


def build(args):
    results, source, config_path, outdir, site_dir = [Path(p).resolve() for p in (args.results, args.source, args.config, args.outdir, args.site_dir)]
    require(not outdir.exists(), "Use a fresh, nonexistent deliverable output directory")
    require(not outdir.is_relative_to(results) and not results.is_relative_to(outdir) and not outdir.is_relative_to(source) and not source.is_relative_to(outdir), "Publication staging must not overlap native results or immutable scientific source")
    require(not outdir.is_relative_to(site_dir.parent), "Publication staging must stay outside the deployed site")
    v = validate_complete(results, source, config_path)
    tf_args = [getattr(args, k, None) for k in ("tf_results", "tf_source", "tf_config")]
    require(not any(tf_args) or all(tf_args), "Supply --tf-results, --tf-source and --tf-config together")
    tf_component = None
    payload = v["payload"]
    if all(tf_args):
        require(all(Path(p).is_absolute() for p in tf_args), "Use absolute TF bound paths")
        tf_results, tf_source, tf_config = map(lambda p: Path(p).resolve(), tf_args)
        for p in (tf_results, tf_source):
            require(not outdir.is_relative_to(p) and not p.is_relative_to(outdir), "Publication staging must not overlap native TF outputs/source")
        tv = validate_tf_complete(tf_results, tf_source, tf_config)
        validate_shared_inputs(v["config"], tv["config"])
        payload = merge_payloads(v["payload"], tv["payload"])
        tf_component = {"results": tf_results, "source": tf_source, "config": tf_config, "validation": tv}
    ps_args = [getattr(args, k, None) for k in ("pyscenic_results", "pyscenic_source", "pyscenic_config")]
    require(not any(ps_args) or all(ps_args), "Supply --pyscenic-results, --pyscenic-source and --pyscenic-config together")
    pyscenic_component = None
    if all(ps_args):
        require(tf_component is not None, "A third pySCENIC component requires the separately completed TF component")
        require(all(Path(p).is_absolute() for p in ps_args), "Use absolute pySCENIC bound paths")
        ps_results, ps_source, ps_config = map(lambda p: Path(p).resolve(), ps_args)
        for p in (ps_results, ps_source):
            require(not outdir.is_relative_to(p) and not p.is_relative_to(outdir), "Publication staging must not overlap native pySCENIC outputs/source")
        pv = validate_pyscenic_complete(ps_results, ps_source, ps_config)
        validate_pyscenic_shared_inputs(v, pv)
        payload = merge_payloads(v["payload"], tv["payload"], pv["payload"])
        pyscenic_component = {"results": ps_results, "source": ps_source, "config": ps_config, "validation": pv}
    require(re.fullmatch(r"atlas-v[0-9]+(?:\.[0-9]+)+", args.tag), "Use an explicit versioned atlas release tag")
    if pyscenic_component:
        require(tuple(int(i) for i in args.tag.removeprefix("atlas-v").split(".")) >= (3, 8, 0), "The completed pySCENIC component requires a new additive release (atlas-v3.8.0 or later); preserve v3.7")
    release = RELEASE + args.tag + "/"
    outdir.mkdir(parents=True)
    assets = outdir / "release_assets"
    assets.mkdir()
    staged = outdir / "site/regulatory"
    staged.mkdir(parents=True)
    ui = ROOT / "site/regulatory"
    for name in ("styles.css", "app.js", "README.md"):
        copy_exact(ui / name, staged / name)
    if tf_component:
        write_combined_transport(staged, payload, v, tv, pv if pyscenic_component else None)
    else:
        copy_exact(results / "atlas/data.js", staged / "data.js")
        copy_exact(results / "atlas/manifest.json", staged / "manifest.json")
    native_manifest = outdir / "NATIVE_FILE_RECEIPTS.json"
    native_receipts = {"files": list(v["files"].values()), "analysis_id": payload["analysis_id"], "source_fingerprint": payload["source_fingerprint"]}
    if tf_component:
        native_receipts["component_receipts"] = {"hdWGCNA": list(v["files"].values()), "curated_tf_activity": list(tv["files"].values())}
        native_receipts["native_prefixes"] = {"hdWGCNA": "results/", "curated_tf_activity": "tf_component/results/"}
    if pyscenic_component:
        native_receipts["component_receipts"]["motif_regulon"] = list(pv["files"].values())
        native_receipts["native_prefixes"]["motif_regulon"] = "pyscenic_component/results/"
    native_manifest.write_bytes(encoded(native_receipts))
    zip_path = assets / "SATB2_regulatory_complete_native_results_source.zip"
    setup = Path(args.setup_dir).resolve() if args.setup_dir else Path(v["config"]["R_wrapper"]["path"]).parent
    archive_contents, archive_members = make_archive(zip_path, results, source, config_path, v, setup, tf_component, payload if tf_component else None, pyscenic_component)
    zip_assets = split_release_asset(zip_path, assets)
    zip_url = release + quote(zip_assets[0].name) if len(zip_assets) == 1 else "downloads.html#native-archive"
    offline_name = "SATB2_regulatory_interactive_offline.html"
    offline_url = release + offline_name
    # Reserve space for complete gallery metadata, frontend and a compact manifest.
    base_bytes = sum(p.stat().st_size for p in all_files(staged))
    current_site_bytes = sum(p.stat().st_size for p in all_files(site_dir.parent)) if site_dir.parent.is_dir() else 0
    replaced_base = sum((site_dir / name).stat().st_size for name in ("styles.css", "app.js", "README.md", "data.js", "manifest.json", "index.html", "gallery.html", "downloads.html", "figure_manifest.json") if (site_dir / name).is_file())
    figure_budget = min(MAX_NEW_PAGES - base_bytes - 1_000_000,
                        MAX_SITE - current_site_bytes + replaced_base - base_bytes - 1_000_000)
    require(figure_budget > 0, "The compact analysis page and metadata reserve exceed the available Pages budget")
    figures = {}
    components = [("hdWGCNA", results, v, "results/")]
    if tf_component:
        components.append(("curated_tf_activity", tf_results, tv, "tf_component/results/"))
    if pyscenic_component:
        components.append(("motif_regulon", ps_results, pv, "pyscenic_component/results/"))
    component_roots = {name: root for name, root, _, _ in components}
    for component, component_root, validation, archive_prefix in components:
        for rel, r in validation["files"].items():
            p = Path(rel)
            if p.suffix.lower() not in (".png", ".pdf", ".svg"):
                continue
            stem = ({"curated_tf_activity": "tf_component/", "motif_regulon": "pyscenic_component/"}.get(component, "")) + str(p.with_suffix(""))
            parts = p.parts
            lineage = next((x for x in ("CA1", "IN") if any(z == x or z.startswith(x + "_") for z in parts)), "other")
            label = " · ".join(part.replace("_", " ") for part in p.with_suffix("").parts if part not in ("figures", "hdwgcna"))
            item = figures.setdefault(stem, {"source_stem": stem, "label": {"curated_tf_activity": "TF activity · ", "motif_regulon": "Motif regulons · "}.get(component, "") + label, "component": component, "lineage": lineage, "network": component if component != "hdWGCNA" else next((x for x in parts if x in ("pooled", "consensus")), "other"), "files": {}})
            item["files"][p.suffix[1:].lower()] = {"native_path": archive_prefix + rel if tf_component else rel, "component_relative_path": rel, "component": component, "bytes": r["bytes"], "sha256": r["sha256"]}
    require(figures, "Completed outputs contain no verified scientific figures")
    remote_assets = []
    # Native PNGs receive first priority, so thumbnails use Pages whenever possible.
    for ext in ("png", "pdf", "svg"):
        for item in figures.values():
            if ext not in item["files"]:
                continue
            record = item["files"][ext]
            p = child(component_roots[record["component"]], record["component_relative_path"])
            with p.open("rb") as handle:
                header = handle.read(1024)
            if ext == "png":
                require(header.startswith(b"\x89PNG\r\n\x1a\n"), f"Invalid native PNG signature: {p}")
            elif ext == "pdf":
                require(header.startswith(b"%PDF-"), f"Invalid native PDF signature: {p}")
            if ext == "svg":
                record["href"] = zip_url
                record["location"] = "complete_native_ZIP"
                continue
            name = release_name(record["native_path"])
            if record["bytes"] <= figure_budget:
                target = staged / "figures" / name
                copy_exact(p, target)
                record["href"] = "figures/" + name
                record["location"] = "Pages"
                figure_budget -= record["bytes"]
            else:
                target = assets / name
                copy_exact(p, target)
                record["href"] = release + quote(name)
                record["location"] = "required_release_asset"
                remote_assets.append(target)
    records = [figures[k] for k in sorted(figures)]
    (staged / "gallery.html").write_bytes(gallery(records, zip_url, offline_url, payload["status"], bool(tf_component), bool(pyscenic_component)))
    tables = [{**r, "path": prefix + r["path"] if tf_component else r["path"]} for _, _, validation, prefix in components for r in validation["files"].values() if r["path"].endswith((".tsv", ".tsv.gz", ".txt"))]
    (staged / "downloads.html").write_bytes(downloads_page(zip_assets, release, offline_url, digest(zip_path), sorted(tables, key=lambda r: r["path"]), bool(tf_component), bool(pyscenic_component)))
    (staged / "figure_manifest.json").write_bytes(encoded({"schema": "satb2-regulatory-native-gallery/v1", "source_fingerprint": payload["source_fingerprint"], "figures": records, "native_bytes_preserved": True, "raster_resizing_performed": False}))
    index = ui.joinpath("index.html").read_text()
    index = re.sub(r'<nav id="publication-links"[\s\S]*?</nav>', "", index)
    index = re.sub(r'<script src="data\.js(?:\?[^"]*)?"></script>', "", index)
    index = re.sub(r'(["\'])(app\.js|styles\.css)(?:\?v=[a-zA-Z0-9]+)?\1', lambda m: m.group(1) + m.group(2) + m.group(1), index)
    require('<script src="app.js"></script>' in index, "Cannot locate the declared regulatory frontend entry point")
    index = index.replace('<script src="app.js"></script>', '<script src="data.js"></script><script src="app.js"></script>')
    index = index.replace('<section id="summary"', '<nav id="publication-links" class="header-links" aria-label="Native completed result downloads"><a href="gallery.html">Native figures and individual PDF / PNG</a><a href="downloads.html">Native tables and complete source ZIP</a><a href="' + offline_url + '">Standalone offline explorer</a></nav><section id="summary"')
    for name in ("styles.css", "app.js", "data.js"):
        p = staged / name
        index = index.replace('"' + name + '"', '"' + name + '?v=' + digest(p)[:12] + '"')
    (staged / "index.html").write_text(index)
    # Vendor libraries are staged only outside Pages to resolve ../vendor inlining.
    vendor = outdir / "site/vendor"
    vendor.mkdir()
    for name in ("fflate.min.js", "jspdf.umd.min.js"):
        copy_exact(ROOT / "site/vendor" / name, vendor / name)
    portable_receipt = portable(staged / "index.html", assets / offline_name, LIVE + "regulatory/")
    (outdir / "PORTABLE_RECEIPT.json").write_bytes(encoded(portable_receipt))
    if len(zip_assets) > 1:
        (assets / "REASSEMBLE_NATIVE_ZIP.txt").write_text("Download every .partNNN asset, then concatenate in numeric order.\n\nmacOS/Linux: cat SATB2_regulatory_complete_native_results_source.zip.part* > SATB2_regulatory_complete_native_results_source.zip\n\nExpected reconstructed ZIP SHA256: " + digest(zip_path) + "\nThe first-part link is not a complete ZIP. Download all parts before extracting.\n")
        remote_assets.append(assets / "REASSEMBLE_NATIVE_ZIP.txt")
    required_assets = list(dict.fromkeys([*zip_assets, assets / offline_name, *remote_assets]))
    asset_rows = [{"asset_name": p.name, "local_path": str(p), "bytes": p.stat().st_size, "sha256": digest(p), "kind": "complete_native_archive" if p in zip_assets else "interactive_offline_HTML" if p.name == offline_name else "native_figure_or_instructions"} for p in required_assets]
    require(all(r["bytes"] < 2 * 1024**3 for r in asset_rows), "A staged upload asset exceeds the GitHub per-asset limit")
    (outdir / "RELEASE_ASSETS.json").write_bytes(encoded({"tag": args.tag, "status": "PREPARED_NOT_UPLOADED", "assets": asset_rows, "every_listed_asset_required_before_publication": True}))
    (outdir / "ARCHIVE_CONTENTS.json").write_bytes(encoded(archive_contents))
    new_bytes = sum(p.stat().st_size for p in all_files(staged))
    require(new_bytes <= MAX_NEW_PAGES, f"New Pages subtree exceeds {MAX_NEW_PAGES} bytes: {new_bytes}")
    current_bytes = sum(p.stat().st_size for p in all_files(site_dir.parent)) if site_dir.parent.is_dir() else 0
    replaced = sum((site_dir / p.relative_to(staged)).stat().st_size for p in all_files(staged) if (site_dir / p.relative_to(staged)).is_file())
    prospective = current_bytes - replaced + new_bytes
    require(prospective <= MAX_SITE, f"Prospective complete Pages site exceeds 950 MB: {prospective}")
    # Confirm all scientific, source and configuration files remain unchanged.
    for row in archive_members:
        require(digest(row["source"], force=True) == row["sha256"], f"Native/source file changed during publication staging: {row['source']}")
    for name, item in v["bindings"].items():
        require(digest(item["path"], force=True) == item["sha256"], f"Original bound input changed during staging: {name}")
    if tf_component:
        for name, item in tv["bindings"].items():
            require(digest(item["path"], force=True) == item["sha256"], f"TF bound input changed during staging: {name}")
    installed = []
    if not args.stage_only:
        for p in all_files(staged):
            target = site_dir / p.relative_to(staged)
            target.parent.mkdir(parents=True, exist_ok=True)
            tmp = target.with_name(target.name + ".publication-writing")
            require(not tmp.exists(), f"An incomplete prior installation exists: {tmp}")
            shutil.copy2(p, tmp)
            require(digest(tmp) == digest(p), f"Installed Pages bytes mismatch: {target}")
            os.replace(tmp, target)
            installed.append(str(target))
    staged_files = [receipt(p, outdir) for p in all_files(outdir)]
    report = {"schema": "satb2-regulatory-publication/v1", "status": "PREPARED_COMPLETE_ONLY", "scientific_fits_performed": False, "GitHub_actions_performed": False, "analysis_id": payload["analysis_id"], "source_fingerprint": payload["source_fingerprint"], "native_final_validation_sha256": digest(results / "FINAL_VALIDATION.json"), "scientific_source": str(source), "configuration": receipt(config_path), "new_pages_bytes": new_bytes, "prospective_site_bytes": prospective, "new_pages_limit": MAX_NEW_PAGES, "total_site_limit": MAX_SITE, "native_figure_groups": len(records), "required_release_assets": len(asset_rows), "complete_native_zip": receipt(zip_path), "archive_members_SHA256_and_CRC_verified": True, "native_source_inputs_unchanged": True, "standalone_HTML": portable_receipt, "stage_only": args.stage_only, "installed_site_files": installed, "staged_files": staged_files, "publication_condition": "Upload and verify every RELEASE_ASSETS.json asset before deploying the site. This script does not publish anything."}
    if tf_component:
        report["components"] = {"hdWGCNA": {"analysis_id": v["final"]["analysis_id"], "source_fingerprint": v["final"]["source_fingerprint"], "native_final_validation_sha256": digest(results / "FINAL_VALIDATION.json"), "source": str(source), "configuration": receipt(config_path)}, "curated_tf_activity": {"analysis_id": tv["final"]["analysis_id"], "source_fingerprint": tv["final"]["source_fingerprint"], "native_final_validation_sha256": digest(tf_results / "FINAL_VALIDATION.json"), "source": str(tf_source), "configuration": receipt(tf_config)}}
        report["TF_activity_run"] = True
        report["motif_regulon_run"] = False
        report["frozen_memberships_and_mouse_rosters_identical"] = True
    if pyscenic_component:
        report["components"]["motif_regulon"] = {"analysis_id": pv["final"]["analysis_id"], "source_fingerprint": pv["final"]["source_fingerprint"], "native_final_validation_sha256": digest(ps_results / "FINAL_VALIDATION.json"), "source": str(ps_source), "configuration": receipt(ps_config)}
        report["motif_regulon_run"] = True
        report["large_public_motif_references_external_hash_bound"] = True
    (outdir / "PUBLICATION_RECEIPT.json").write_bytes(encoded(report))
    verify_outdir(outdir)
    return report


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--results", type=Path)
    p.add_argument("--source", type=Path)
    p.add_argument("--config", type=Path)
    p.add_argument("--outdir", type=Path)
    p.add_argument("--site-dir", type=Path, default=ROOT / "site/regulatory")
    p.add_argument("--setup-dir", type=Path)
    p.add_argument("--tf-results", type=Path, help="Optional separately completed native TF results")
    p.add_argument("--tf-source", type=Path, help="Executed immutable TF scientific source; requires all three TF arguments")
    p.add_argument("--tf-config", type=Path, help="Immutable TF configuration; requires all three TF arguments")
    p.add_argument("--pyscenic-results", type=Path, help="Optional third independently completed native positive pySCENIC results")
    p.add_argument("--pyscenic-source", type=Path, help="Executed immutable pySCENIC source; requires all three pySCENIC arguments plus TF component")
    p.add_argument("--pyscenic-config", type=Path, help="Immutable pySCENIC configuration; requires all three pySCENIC arguments plus TF component")
    p.add_argument("--tag", default=TAG)
    p.add_argument("--stage-only", action="store_true", help="Prepare complete publication artifacts without copying into the repository site")
    p.add_argument("--verify-outdir", type=Path, help="Validate an already prepared deliverable directory; no model or upload is executed")
    args = p.parse_args()
    if args.verify_outdir:
        result = verify_outdir(args.verify_outdir)
    else:
        require(all(getattr(args, k) is not None for k in ("results", "source", "config", "outdir")), "--results, --source, --config and a fresh --outdir are required")
        require(all(getattr(args, k).is_absolute() for k in ("results", "source", "config", "outdir")), "Use absolute bound paths")
        result = build(args)
        result = {k: result[k] for k in ("status", "analysis_id", "source_fingerprint", "new_pages_bytes", "prospective_site_bytes", "native_figure_groups", "required_release_assets", "stage_only")}
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
