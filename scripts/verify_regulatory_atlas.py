#!/usr/bin/env python3
"""Verify published regulatory transport and local figures; no network or analysis.

The fixed cohort is the frozen Seurat five-endpoint follow-up. This validator
checks presentation/source contracts, not scientific model adequacy. Remote
release assets are checked for declared immutable URLs/hashes only; upload-byte
verification remains a separate publication prerequisite.
"""
from __future__ import annotations

import argparse
import base64
import gzip
import hashlib
import io
import json
import math
import re
from pathlib import Path
from urllib.parse import unquote, urlsplit


ROOT = Path(__file__).resolve().parent.parent
EXPECTED = {
    "CA1_DL": ("CA1", "P0244_seurat_CA1_ProS_LAYER_PROGRAM", 8750, "b92fba8fbfdbef55982732b72d50f56fe81c46f4e6c38e0865103ee4ac4211c4"),
    "CA1_SL": ("CA1", "P0247_seurat_CA1_ProS_LAYER_PROGRAM", 11706, "71cde18ce1da9b2d2439bdc687e69ffeced31af52a962343e44f944f68dca6f6"),
    "IN_C1": ("IN", "P0140_seurat_Inhibitory_CLUSTER", 1718, "bff0ec716a0d4b3d76c38828465fb59ab50184690cf609f3d184fa4ed1db6189"),
    "IN_C6": ("IN", "P0145_seurat_Inhibitory_CLUSTER", 501, "ad7adf54e947a9e7c50da527fe06026801f8a9ec682d59f25ef79440c3ce0b64"),
    "IN_C8": ("IN", "P0147_seurat_Inhibitory_CLUSTER", 468, "1252e05b37f591b850bbc4c43fda7b3ace78be52fde723e65382a8e9c371c7fd"),
}
MICE = {**{s: "WT" for s in ("wt_1122", "wt_1123", "wt_1178", "wt_1220", "wt_1226", "wt_1228")},
        **{s: "cKO" for s in ("cKO_933", "cKO_991", "cKO_1138", "cKO_1139", "cKO_1160", "cKO_1165")}}
RELEASE = "https://github.com/Ali-Mujahid/satb2_snRNA/releases/download/atlas-v"
READY = {"complete", "completed", "partial"}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def local(root, relative):
    parts = urlsplit(relative)
    require(not parts.scheme and not parts.netloc and not Path(unquote(parts.path)).is_absolute(), f"Expected a local relative asset: {relative}")
    p = (Path(root) / unquote(parts.path)).resolve()
    require(p.is_relative_to(Path(root).resolve()) and p.is_file() and not p.is_symlink(), f"Missing/escaping local regulatory asset: {relative}")
    return p


def matching(path, row):
    require(isinstance(row.get("bytes"), int) and row["bytes"] >= 0 and re.fullmatch(r"[0-9a-f]{64}", row.get("sha256", "")), "Invalid transport/figure receipt")
    require(Path(path).stat().st_size == row["bytes"] and sha(path) == row["sha256"], f"Bytes/SHA256 changed: {path}")


def reject_synthetic(value):
    if isinstance(value, dict):
        require(value.get("synthetic") is not True, "Synthetic fixtures cannot be production regulatory results")
        for key, item in value.items():
            if key == "analysis_id":
                require("SYNTHETIC" not in str(item).upper(), "Synthetic analysis identity is not publishable")
            reject_synthetic(item)
    elif isinstance(value, list):
        for item in value:
            reject_synthetic(item)


def payload_from(path):
    match = re.fullmatch(r"\s*(?:/\*[\s\S]*?\*/\s*)?window\.SATB2_REGULATORY_DATA_GZIP\s*=\s*(\"[A-Za-z0-9+/=]*\")\s*;\s*", Path(path).read_text())
    require(match, "Regulatory data must be the declared compressed JSON assignment")
    packed = base64.b64decode(json.loads(match.group(1)), validate=True)
    with gzip.GzipFile(fileobj=io.BytesIO(packed)) as stream:
        raw = stream.read(512_000_001)
    require(len(raw) <= 512_000_000, "Regulatory transport exceeds the bounded JSON size")
    payload = json.loads(raw)
    reject_synthetic(payload)
    return payload


def unique(rows, label):
    ids = [r.get("id") for r in rows]
    require(all(isinstance(i, str) and i for i in ids) and len(set(ids)) == len(ids), f"Duplicate/invalid {label} IDs")
    return {r["id"]: r for r in rows}


def validate_payload(payload, manifest, require_tf=True):
    require(payload.get("schema") == "satb2-regulatory-atlas/v1" and payload.get("status") in ("complete", "partial"), "No completed production regulatory payload")
    require(payload.get("analysis_id") == manifest.get("analysis_id") and payload.get("source_fingerprint") == manifest.get("source_fingerprint") and re.fullmatch(r"[0-9a-f]{64}", payload.get("source_fingerprint", "")), "Transport source identities disagree")
    eps = unique(payload.get("endpoints", []), "endpoint")
    require(set(eps) == set(EXPECTED), "Wrong frozen regulatory endpoint cohort")
    for key, (family, source, count, membership) in EXPECTED.items():
        require(tuple(eps[key].get(k) for k in ("family", "source_endpoint_id", "n_nuclei", "membership_sha256")) == (family, source, count, membership), f"Frozen endpoint count/source/membership changed: {key}")
    samples = unique(payload.get("samples", []), "mouse")
    require({k: r.get("genotype") for k, r in samples.items()} == MICE, "Frozen twelve-mouse genotype roster differs")
    methods = unique(payload.get("methods", []), "method")
    hd_ids = {"hdwgcna_" + lineage + "_" + branch for lineage in ("CA1", "IN") for branch in ("pooled", "consensus")}
    for name in ("hdwgcna_CA1_pooled", "hdwgcna_IN_pooled"):
        require(methods.get(name, {}).get("family") == "hdWGCNA" and methods[name].get("status") == "complete", f"Primary network is not complete: {name}")
    combined = manifest.get("schema") == "satb2-regulatory-combined-transport/v1"
    require(combined or not require_tf, "This publication requires independently completed HD and curated TF components")
    if combined:
        require(set(methods) == hd_ids | {"collectri_ulm", "pyscenic_motif"}, "Unexpected combined source method IDs")
        require(manifest.get("TF_activity_run") is True and manifest.get("motif_regulon_run") is False, "Combined method statuses are not explicit")
        require(methods.get("collectri_ulm", {}).get("family") == "curated_tf_activity" and methods["collectri_ulm"].get("status") == "complete", "Curated TF primary method is not complete")
        require(methods.get("pyscenic_motif", {}).get("family") == "motif_regulon" and methods["pyscenic_motif"].get("status") == "not_available", "Blocked motif inference must remain unavailable")
        components = manifest.get("components", {})
        for lineage in ("CA1", "IN"):
            require(components.get("hdWGCNA", {}).get("branches", {}).get(lineage + "_pooled", {}).get("validated_complete") is True, f"Native HD component is not complete: {lineage}")
        tf = components.get("curated_tf_activity", {})
        require(tf.get("completed_primary") == "curated_mouse_CollecTRI_ULM_and_independent_mouse_associations" and tf.get("blocked_optional") == "pySCENIC_motif_cache_missing", "Native TF component completion/block status differs")
        p = payload.get("provenance", {})
        require(p.get("TF_activity") == "COMPLETE_CURATED_COLLECTRI_ULM" and p.get("motif_supported_regulons", "").startswith("BLOCKED_") and p.get("no_scientific_fits_in_publication_builder") is True, "Combined provenance misstates curated activity or motif inference")
        require(set(p.get("components", {})) == {"hdWGCNA", "curated_tf_activity"}, "Historical component provenance is missing")
        for name in ("hdWGCNA", "curated_tf_activity"):
            ident = p.get("component_identities", {}).get(name, {})
            require(ident.get("analysis_id") == components.get(name, {}).get("analysis_id") and ident.get("source_fingerprint") == components.get(name, {}).get("source_fingerprint") and re.fullmatch(r"[0-9a-f]{64}", ident.get("source_fingerprint", "")), f"Native component source identity differs: {name}")
    else:
        require(set(methods) == hd_ids and all(m.get("family") == "hdWGCNA" for m in methods.values()), "An HD-only transport cannot declare completed TF results")
    require(all(m.get("family") in {"hdWGCNA", "curated_tf_activity", "motif_regulon"} for m in methods.values()), "Undeclared regulatory method family")
    require(all(m.get("status") == "not_available" for m in methods.values() if m.get("family") == "motif_regulon"), "Unavailable motif inference cannot declare an active result")
    require(not payload.get("embeddings"), "No native nucleus-activity embedding was validated for this publication")
    features = unique(payload.get("features", []), "feature")
    for f in features.values():
        require(f.get("method_id") in methods and (not f.get("endpoint_id") or f["endpoint_id"] in eps), "Feature source is unavailable")
    seen, groups = set(), {}
    for r in payload.get("values", []):
        key = tuple(r.get(k) for k in ("method_id", "endpoint_id", "feature_id", "sample_id"))
        method, ep, feature, mouse = key
        require(key not in seen, "Duplicate independent-mouse regulatory value")
        seen.add(key)
        require(method in methods and ep in eps and feature in features and mouse in samples and methods[method].get("status") in READY, "Mouse activity has an absent or unavailable source key")
        f = features[feature]
        require(f["method_id"] == method and (not f.get("endpoint_id") or f["endpoint_id"] == ep), "Mouse value crosses incompatible feature source")
        finite = isinstance(r.get("value"), (float, int)) and not isinstance(r["value"], bool) and math.isfinite(r["value"])
        missing = r.get("value") is None and r.get("eligible") is False and bool(r.get("missing_reason"))
        require(isinstance(r.get("eligible"), bool) and (finite or missing), "Mouse score must be finite or explicitly unavailable")
        groups.setdefault(key[:3], set()).add(mouse)
    require(groups and all(ids == set(MICE) for ids in groups.values()), "A recorded feature/endpoint lacks all twelve independent mice")
    for download in payload.get("downloads", []):
        for key in ("path", "href", "url"):
            if key in download:
                require(str(download[key]).startswith("https://"), "An active payload download points to a broken local scientific path")
    return {"combined": combined, "methods": len(methods), "features": len(features), "mouse_values": len(seen)}


def verify(site_dir, require_tf=True):
    site = Path(site_dir).resolve()
    manifest = json.loads(local(site, "manifest.json").read_text())
    reject_synthetic(manifest)
    require("data.js" in manifest.get("transport", {}), "No production data.js checksum")
    for name, row in manifest["transport"].items():
        matching(local(site, name), row)
    payload = payload_from(local(site, "data.js"))
    summary = validate_payload(payload, manifest, require_tf)
    for field, count in manifest.get("record_counts", {}).items():
        require(isinstance(payload.get(field), list) and len(payload[field]) == count, f"Declared transport row count differs: {field}")
    gallery = json.loads(local(site, "figure_manifest.json").read_text())
    require(gallery.get("source_fingerprint") == payload["source_fingerprint"] and gallery.get("native_bytes_preserved") is True and gallery.get("raster_resizing_performed") is False, "Figure gallery source contract differs")
    require(gallery.get("figures"), "No verified native figures")
    checked, remote = 0, 0
    for group in gallery["figures"]:
        require(group.get("component") in {"hdWGCNA", "curated_tf_activity"}, "Figure has no declared source component")
        for ext, row in group.get("files", {}).items():
            require(ext in {"png", "pdf", "svg"} and re.fullmatch(r"[0-9a-f]{64}", row.get("sha256", "")) and isinstance(row.get("bytes"), int), "Figure receipt is invalid")
            href = row.get("href", "")
            if row.get("location") == "Pages":
                p = local(site, href); matching(p, row)
                with p.open("rb") as handle:
                    magic = handle.read(8)
                require(ext == "png" and magic == b"\x89PNG\r\n\x1a\n" or ext == "pdf" and magic.startswith(b"%PDF-"), "Local native figure has an invalid signature")
                checked += 1
            elif row.get("location") == "required_release_asset":
                require(href.startswith(RELEASE) and bool(re.fullmatch(r"https://github\.com/Ali-Mujahid/satb2_snRNA/releases/download/atlas-v[0-9.]+/[^/]+", href)), "Native remote figure is not a versioned declared release asset")
                remote += 1
            else:
                require(ext == "svg" and row.get("location") == "complete_native_ZIP" and (href.startswith(RELEASE) or href == "downloads.html#native-archive"), "Native SVG/archive link contract differs")
    for page in ("index.html", "gallery.html", "downloads.html"):
        text = local(site, page).read_text()
        for attr, href in re.findall(r'\b(src|href)="([^"]+)"', text):
            if href.startswith(("https://", "http://", "#", "data:", "mailto:")):
                continue
            # Parent atlas/vendor links are allowed only within the site root.
            p = (site / unquote(urlsplit(href).path)).resolve()
            parent_ok = p.is_relative_to(site.parent) and p.is_file()
            if not parent_ok and href.startswith("../"):
                # A fresh deliverable stages only the new subtree/vendor files.
                # Existing parent navigation is checked against its destination.
                q = (ROOT / "site/regulatory" / unquote(urlsplit(href).path)).resolve()
                parent_ok = q.is_relative_to(ROOT / "site") and q.is_file()
            require(parent_ok, f"Broken local {page} {attr}: {href}")
    return {"status": "PASS", **summary, "exact_endpoints": 5, "independent_mice": 12, "local_native_figures_checked": checked, "remote_release_figures_declared_not_network_checked": remote, "network_calls": 0, "scientific_fits": 0}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--site-dir", type=Path, default=ROOT / "site/regulatory")
    parser.add_argument("--allow-hd-only", action="store_true", help="Validate an intentionally HD-only completed release; never infer TF completion")
    args = parser.parse_args()
    print(json.dumps(verify(args.site_dir, require_tf=not args.allow_hd_only), indent=2))


if __name__ == "__main__":
    main()
