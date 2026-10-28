#!/usr/bin/env python3
"""Package validated completed hdWGCNA outputs; never fit a scientific model.

Called only by the registered, approval-bound native workflow. Source authoring
and synthetic tests do not authorize execution against real analysis outputs.
"""
from __future__ import annotations

import argparse
import base64
import csv
import gzip
import hashlib
import json
import math
import re
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

SCHEMA = "satb2-regulatory-atlas/v1"
ENDPOINT_MAP = {"CA1_Deep": "CA1_DL", "CA1_Superficial": "CA1_SL", "IN_C1": "IN_C1", "IN_C6": "IN_C6", "IN_C8": "IN_C8"}
LINEAGES = {"CA1": ["CA1_Deep", "CA1_Superficial"], "IN": ["IN_C1", "IN_C6", "IN_C8"]}
LABELS = {"CA1_Deep": "CA1 deep-associated expression program", "CA1_Superficial": "CA1 superficial-associated expression program", "IN_C1": "Inhibitory C1", "IN_C6": "Inhibitory C6", "IN_C8": "Inhibitory C8 · mixed identity review"}
NATIVE_COUNTS = {"CA1_Deep": 8750, "CA1_Superficial": 11706, "IN_C1": 1718, "IN_C6": 501, "IN_C8": 468}
_hash_cache = {}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def digest(path, algorithm="sha256", force=False):
    path = Path(path).resolve()
    st = path.stat()
    key = (str(path), algorithm, st.st_size, st.st_mtime_ns)
    if not force and key in _hash_cache:
        return _hash_cache[key]
    h = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(block)
    value = h.hexdigest()
    _hash_cache[key] = value
    return value


def json_bytes(value, compact=False):
    return (json.dumps(value, ensure_ascii=False, allow_nan=False, sort_keys=True,
                       separators=(",", ":") if compact else None,
                       indent=None if compact else 2) + "\n").encode("utf-8")


def load_json(path):
    return json.loads(Path(path).read_text())


def table(path, required=(), optional=False):
    path = Path(path)
    if not path.exists():
        require(optional, f"Required native table is missing: {path}")
        return []
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        # R's empty data.frame() may produce a zero-column/zero-row file.
        if not reader.fieldnames or reader.fieldnames == [""]:
            require(optional, f"Required native table is empty: {path}")
            return []
        require(set(required) <= set(reader.fieldnames), f"Native columns missing from {path}: {set(required) - set(reader.fieldnames)}")
        rows = list(reader)
        require(all(None not in r for r in rows), f"Malformed TSV row: {path}")
        return rows


def number(value, required=False, field="number"):
    if value is None or str(value).strip() in ("", "NA", "NaN", "nan", "null", "None"):
        require(not required, f"Required finite native {field} is missing")
        return None
    out = float(value)
    require(math.isfinite(out), f"Non-finite native {field}: {value}")
    return out


def integer(value, required=True, field="count"):
    out = number(value, required, field)
    if out is None:
        return None
    require(out >= 0 and out.is_integer(), f"Expected nonnegative integer {field}: {value}")
    return int(out)


def safe_child(root, relative):
    root = Path(root).resolve()
    path = (root / relative).resolve()
    require(path.is_relative_to(root), f"Native receipt escapes its output root: {relative}")
    return path


def verify_input_bindings(cfg, phase):
    records = []
    collections = [("inputs", cfg.get("inputs", {})), ("existing_deseq2", cfg.get("existing_deseq2", {})), ("source_pseudobulks", cfg.get("source_pseudobulks", {}))]
    require(all(items for _, items in collections), "Missing hash-bound scientific input collection")
    for group, items in collections:
        records.extend((group + ":" + name, item) for name, item in items.items())
    records.extend((key, cfg[key]) for key in ("strict_lr_pooled", "runtime_manifest", "R_wrapper"))
    result = {}
    for name, item in records:
        path = Path(item["path"])
        require(path.is_absolute() and path.is_file(), f"Bound input is unavailable: {name}")
        require(re.fullmatch(r"[a-fA-F0-9]{64}", item.get("sha256", "")) is not None, f"Input lacks SHA256: {name}")
        observed = digest(path, force=phase == "after")
        require(observed == item["sha256"].lower(), f"Source input changed ({phase}): {name}")
        require(item.get("bytes") is None or path.stat().st_size == int(item["bytes"]), f"Source input size changed: {name}")
        result[name] = {"path": str(path.resolve()), "bytes": path.stat().st_size, "sha256": observed}
    return result


def verify_prepared(root, cfg):
    prep = root / "prepared"
    complete = load_json(prep / "COMPLETE.json")
    require(complete.get("status") == "COMPLETE" and complete.get("source_files_unchanged") is True and complete.get("source_all_gene_pseudobulk_matches") is True, "Prepared counts / memberships did not pass source verification")
    require(complete.get("group_counts") == NATIVE_COUNTS == cfg.get("expected_group_counts"), "Prepared frozen endpoint counts changed")
    for rel, expected in complete.get("files", {}).items():
        p = safe_child(prep, rel)
        require(p.is_file() and digest(p) == expected["sha256"] and p.stat().st_size == expected["bytes"], f"Prepared output receipt mismatch: {rel}")
    audit = table(prep / "frozen_membership_audit.tsv", ["endpoint", "population_key", "retained_nuclei", "membership_sha256", "native_population_membership_sha256"])
    by_endpoint = {r["endpoint"]: r for r in audit}
    require(len(by_endpoint) == len(audit) == 5 and set(by_endpoint) == set(ENDPOINT_MAP), "Frozen membership audit must contain five unique requested endpoints")
    for ep, r in by_endpoint.items():
        require(integer(r["retained_nuclei"]) == NATIVE_COUNTS[ep], f"Frozen nucleus count differs: {ep}")
        require(re.fullmatch(r"[a-f0-9]{64}", r["membership_sha256"]) and r["membership_sha256"] == r["native_population_membership_sha256"], f"Frozen/native membership hashes disagree: {ep}")
        require(r["population_key"] == cfg["existing_deseq2"][ep]["source_endpoint_id"], f"Original DESeq2 endpoint differs: {ep}")
    coverage = table(prep / "mouse_coverage.tsv", ["endpoint", "sample", "condition", "n_nuclei"])
    sample_conditions = {}
    for r in coverage:
        require(r["endpoint"] in ENDPOINT_MAP and r["condition"] in ("WT", "cKO"), "Unknown prepared endpoint / genotype")
        require(r["sample"] not in sample_conditions or sample_conditions[r["sample"]] == r["condition"], "Mouse genotype changed across endpoints")
        sample_conditions[r["sample"]] = r["condition"]
    require(len(sample_conditions) == 12 and Counter(sample_conditions.values()) == {"WT": 6, "cKO": 6}, "Expected six independent WT and six cKO mice")
    by_key = {(r["endpoint"], r["sample"]): r for r in coverage}
    require(len(by_key) == len(coverage) == 60, "Prepared coverage must contain 60 unique endpoint × mouse rows")
    for ep in ENDPOINT_MAP:
        require({sample for endpoint, sample in by_key if endpoint == ep} == set(sample_conditions), f"Incomplete twelve-mouse coverage: {ep}")
        require(sum(integer(r["n_nuclei"]) for r in coverage if r["endpoint"] == ep) == NATIVE_COUNTS[ep], f"Coverage nucleus total differs: {ep}")
    return complete, by_endpoint, sample_conditions, by_key


def branch_status(folder):
    if (folder / "failure.json").is_file():
        return load_json(folder / "failure.json")
    if (folder / "status.json").is_file():
        return load_json(folder / "status.json")
    return {"status": "NOT_AVAILABLE", "reason": "No completed native branch status"}


def verify_checkpoint(folder):
    receipt = load_json(folder / "completed.json")
    require(receipt.get("fingerprint") and isinstance(receipt.get("files"), dict), f"Native completed branch lacks a checksum receipt: {folder}")
    for rel, expected in receipt["files"].items():
        p = safe_child(folder, rel)
        require(p.is_file() and digest(p, "md5") == expected, f"Native completed checkpoint checksum mismatch: {p}")
    return receipt


def feature_id(method, module):
    return method + "::" + str(module)


def node_id(method, module, symbol):
    return json.dumps([method, str(module), symbol], separators=(",", ":"), ensure_ascii=False)


def append_branch(payload, cfg, root, lineage, network, sample_conditions, source_coverage):
    folder = root / "hdwgcna" / lineage / network
    status = branch_status(folder)
    method = "hdwgcna_" + lineage + "_" + network
    eps = LINEAGES[lineage]
    meta = {"id": method, "label": f"{lineage} {network} hdWGCNA · exploratory", "family": "hdWGCNA", "status": "complete" if status.get("status") == "NETWORK_COMPLETE" else "not_available", "score_label": "Native unharmonized positive-RNA-oriented module PC1", "units": "Native module eigengene PC1 units", "activity_unit": "mouse", "endpoint_ids": [ENDPOINT_MAP[e] for e in eps], "network_edge_basis": "Undirected top-hub metacell RNA bicor coexpression", "interpretation": "Exploratory fixed lineage modules. PC1 is oriented to correlate positively with mean module RNA; eigengenes are not harmonized by mouse or genotype. Associations are cKO minus WT continuous-score effects, not gene DEGs. Technical partition partly confounds genotype. Shared lineage edges are not refitted by endpoint selection. No TF binding or causal regulation is inferred.", "eligibility_policy": "One finite mean score per independent mouse × exact frozen endpoint × native module. Metacell-discovery eligibility is recorded separately; all frozen nuclei are projected and scored.", "resources": {"hdWGCNA": cfg["hdwgcna"], "native_status": status, "network": network, "lineage": lineage}, "source_files": []}
    payload["methods"].append(meta)
    if meta["status"] != "complete":
        require(network != "pooled", f"Primary {lineage} network unavailable: {status}")
        return {"status": status.get("status"), "validated_complete": False, "method_id": method}
    checkpoint = verify_checkpoint(folder)
    native = {name: folder / filename for name, filename in {"means": "module_eigengenes_mouse.tsv", "membership": "module_membership.tsv", "hubs": "hub_genes.tsv", "edges": "network_edges.tsv"}.items()}
    require(all(p.is_file() for p in native.values()), f"Completed branch is missing a stable native alias: {method}")
    means = table(native["means"], ["sample_id", "condition", "endpoint", "module", "n_nuclei", "eigengene_mean"])
    membership = table(native["membership"], ["module", "gene_symbol", "gene_ids", "kME", "assigned"])
    hubs = table(native["hubs"], ["module", "gene_symbol", "gene_ids", "kME", "hub_rank"], optional=True)
    edges = table(native["edges"], ["module", "gene1", "gene2", "bicor", "edge_type"], optional=True)
    modules = sorted({r["module"] for r in means})
    require(modules and set(modules) <= {r["module"] for r in membership}, f"Module means / native membership disagree: {method}")
    counts = Counter((r["endpoint"], r["module"]) for r in means)
    require(set(counts) == {(ep, m) for ep in eps for m in modules} and all(n == 12 for n in counts.values()), f"Every endpoint × module needs twelve mouse scores: {method}")
    seen = set()
    for r in means:
        sample, ep, module = r["sample_id"], r["endpoint"], r["module"]
        require(ep in eps and sample in sample_conditions and sample_conditions[sample] == r["condition"], f"Mouse / genotype / endpoint mismatch: {method}")
        k = (ep, module, sample)
        require(k not in seen, f"Duplicate module × mouse result: {method} {k}")
        seen.add(k)
        require(integer(r["n_nuclei"]) == integer(source_coverage[(ep, sample)]["n_nuclei"]), f"Scored nucleus coverage differs from exact source: {k}")
        value = number(r["eigengene_mean"], True, "eigengene_mean")
        if r.get("value") not in (None, "", "NA"):
            require(number(r["value"], True) == value, f"Native value alias differs: {method} {k}")
        payload["values"].append({**r, "method_id": method, "endpoint_id": ENDPOINT_MAP[ep], "feature_id": feature_id(method, module), "sample_id": sample, "value": value, "eligible": True, "missing_reason": "", "n_nuclei": integer(r["n_nuclei"]), "eigengene_mean": value, "module_RNA_mean": number(r.get("module_RNA_mean")), "source_fit_id": method, "normalization": "Native fixed PC1 projected to original normalized nuclei, then arithmetic mean within each independent mouse and endpoint"})
    for m in modules:
        rows = [r for r in membership if r["module"] == m]
        colors = {r.get("color", "") for r in rows} - {""}
        # WGCNA color names are retained as provenance; only hex values style the UI.
        native_hex = rows[0].get("color_hex", "") if rows else ""
        require(not native_hex or re.fullmatch(r"#[0-9a-fA-F]{6}", native_hex), f"Invalid native module color_hex: {method} {m}")
        payload["features"].append({"id": feature_id(method, m), "method_id": method, "label": f"{lineage} {network} · {m}", "kind": "module", "source_id": m, "module_size": len(rows), "color": native_hex or None, "native_module_colors": sorted(colors), "assigned_module": m not in ("grey", "gray", "0", "unassigned")})
    nodes = {}
    for r in hubs:
        require(r["module"] in modules and r["gene_symbol"], f"Hub without a native module / symbol: {method}")
        n = node_id(method, r["module"], r["gene_symbol"])
        require(n not in nodes, f"Duplicate native hub symbol within module: {n}")
        row = {**r, "method_id": method, "endpoint_id": None, "feature_id": feature_id(method, r["module"]), "gene_id": r["gene_ids"], "kME": number(r["kME"], True, "kME"), "hub_rank": integer(r["hub_rank"])}
        payload["hubs"].append(row)
        nodes[n] = {"id": n, "method_id": method, "endpoint_id": None, "feature_id": row["feature_id"], "gene_id": r["gene_ids"], "label": r["gene_symbol"], "kind": "module_gene", "color": r.get("color_hex") or None, "kME": row["kME"], "hub_rank": row["hub_rank"]}
    payload["nodes"].extend(nodes.values())
    skipped = []
    for i, r in enumerate(edges):
        a, b = node_id(method, r["module"], r["gene1"]), node_id(method, r["module"], r["gene2"])
        require(a in nodes and b in nodes and a != b, f"Native weighted edge has absent or identical hub nodes: {method} row {i}")
        weight = number(r["bicor"])
        if weight is None:
            skipped.append({"native_row_index": i, "reason": "Native bicor unavailable; not converted to zero", **r})
            continue
        require(abs(weight) <= 1 + 1e-9, f"Invalid native bicor: {weight}")
        payload["edges"].append({**r, "method_id": method, "endpoint_id": None, "feature_id": feature_id(method, r["module"]), "source": a, "target": b, "weight": weight, "sign": "positive" if weight > 0 else "negative" if weight < 0 else "unknown", "directed": False, "edge_basis": r["edge_type"]})
    mc_path = folder / "metacell_coverage.tsv"
    if not mc_path.exists() and network == "consensus":
        mc_path = root / "hdwgcna" / lineage / "pooled/metacell_coverage.tsv"
    mc = table(mc_path, ["sample_id", "endpoint", "n_nuclei", "n_metacells", "discovery_status"])
    require(len(mc) == len(eps) * 12 and len({(r["endpoint"], r["sample_id"]) for r in mc}) == len(mc), f"Metacell coverage axis differs: {method}")
    for r in mc:
        require(r["endpoint"] in eps and r["sample_id"] in sample_conditions, f"Unknown metacell coverage key: {method}")
        payload["coverage"].append({**r, "method_id": method, "endpoint_id": ENDPOINT_MAP[r["endpoint"]], "eligible": r["discovery_status"] == "included", "reason": r["discovery_status"], "n_nuclei": integer(r["n_nuclei"]), "n_metacells": integer(r["n_metacells"]), "coverage_kind": "metacell_network_discovery; separate from finite mouse-score eligibility", "all_nuclei_scored": True})
    effect_path = root / "hdwgcna" / f"{network}_module_genotype_associations.tsv"
    effects = table(effect_path, ["lineage", "network", "endpoint", "module", "effect_cKO_minus_WT", "statistic", "pvalue", "padj", "model_status", "fdr_family"])
    chosen = [r for r in effects if r["lineage"] == lineage and r["network"] == network]
    require(len(chosen) == len(eps) * len(modules), f"Globally BH-adjusted effect coverage differs: {method}")
    effect_keys = set()
    for r in chosen:
        k = (r["endpoint"], r["module"])
        require(k not in effect_keys and k[0] in eps and k[1] in modules, f"Duplicate or unknown adjusted effect: {method} {k}")
        effect_keys.add(k)
        p, q = number(r["pvalue"]), number(r["padj"])
        require(all(v is None or 0 <= v <= 1 for v in (p, q)), "Invalid saved effect P/FDR")
        require("all tested CA1+IN" in r["fdr_family"], "Use the output-root CA1+IN global-BH table, never branch unadjusted aliases")
        payload["effects"].append({**r, "method_id": method, "endpoint_id": ENDPOINT_MAP[k[0]], "feature_id": feature_id(method, k[1]), "contrast": "cKO_minus_WT", "effect_name": "Independent-mouse mean native module eigengene difference", "effect": number(r["effect_cKO_minus_WT"]), "statistic": number(r["statistic"]), "pvalue": p, "padj": q, "tested": r["model_status"] == "mouse_level_lm", "n_WT": integer(r.get("n_WT"), False), "n_cKO": integer(r.get("n_cKO"), False), "source_fit_id": method, "source_adjusted_table": str(effect_path.relative_to(root))})
    meta["resources"]["skipped_unavailable_edge_weights"] = skipped
    meta["resources"]["completed_checkpoint_fingerprint"] = checkpoint["fingerprint"]
    for path in [*native.values(), effect_path, folder / "status.json", folder / "completed.json", mc_path]:
        meta["source_files"].append({"path": str(path.relative_to(root)), "sha256": digest(path), "bytes": path.stat().st_size})
    return {"status": "NETWORK_COMPLETE", "validated_complete": True, "method_id": method, "n_modules": len(modules), "n_mouse_values": len(means), "n_hubs": len(hubs), "n_native_edges": len(edges), "n_unavailable_edge_weights": len(skipped)}


def write_identical(path, content):
    path = Path(path)
    if path.exists():
        require(path.read_bytes() == content, f"Refusing to overwrite different existing package content: {path}")
    else:
        path.parent.mkdir(parents=True, exist_ok=True)
        tmp = path.with_name(path.name + ".writing")
        require(not tmp.exists(), f"Incomplete write already exists: {tmp}")
        tmp.write_bytes(content)
        tmp.replace(path)


def bound_record(item, label, force=False):
    require(isinstance(item, dict) and {"path", "bytes", "sha256"} <= set(item), f"Missing complete recovery binding: {label}")
    path = Path(item["path"])
    require(path.is_absolute() and path.is_file() and not path.is_symlink(), f"Recovery binding is not a regular absolute file: {label}")
    require(re.fullmatch(r"[a-f0-9]{64}", item["sha256"]) is not None, f"Invalid recovery SHA256: {label}")
    require(path.stat().st_size == item["bytes"] and digest(path, force=force) == item["sha256"], f"Recovery bound file changed: {label}")
    return {"path": str(path.resolve()), "bytes": path.stat().st_size, "sha256": item["sha256"]}


def verify_native_reuse(cfg, config_path, results, phase="before"):
    """Validate copy-only recovery without relabeling its original network fits."""
    reuse = cfg.get("reuse_completed_native")
    if reuse is None:
        return None
    require(isinstance(reuse, dict) and reuse.get("schema") == "satb2-hdwgcna-copy-recovery/v1", "Unsupported explicit native-fit reuse schema")
    require(reuse.get("no_network_refits") is True, "Copy-only recovery must prohibit network refits")
    force = phase == "after"
    native_record = bound_record(reuse["native_config"], "original native configuration", force)
    inventory_record = bound_record(reuse["inventory"], "frozen native-copy inventory", force)
    native_cfg = load_json(native_record["path"])
    allowed = {"run_root", "workflow_version", "reuse_completed_native"}
    require({k: v for k, v in native_cfg.items() if k not in allowed} == {k: v for k, v in cfg.items() if k not in allowed},
            "Recovery changed a scientific input/parameter/configuration key")
    source = Path(reuse["source_results"])
    require(source.is_absolute() and source.is_dir() and not source.is_symlink(), "Original completed native results are unavailable")
    source = source.resolve()
    require(source == Path(native_cfg["run_root"]).resolve() / "results" and source != results,
            "Original native configuration/result identity differs, or recovery would overwrite it")
    require(not results.is_relative_to(source) and not source.is_relative_to(results), "Recovery and original results must be separate trees")
    inventory = load_json(inventory_record["path"])
    require(inventory.get("schema") == "satb2-hdwgcna-copy-inventory/v1" and inventory.get("status") == "FROZEN_COMPLETE_FIT_COPY_SOURCE", "Native-copy source inventory is not frozen complete")
    require(Path(inventory["source_results"]).resolve() == source and inventory.get("native_config") == reuse["native_config"],
            "Frozen native-copy inventory refers to a different original run/configuration")
    def records(items):
        require(isinstance(items, list) and items, "Native-copy inventory has no scientific files")
        by_path = {}
        for item in items:
            rel = item["path"]
            require(isinstance(rel, str) and not Path(rel).is_absolute() and ".." not in Path(rel).parts and "\\" not in rel
                    and Path(rel).parts[0] in ("prepared", "hdwgcna"), "Native-copy receipt includes files outside prepared/hdwgcna")
            require(rel not in by_path and re.fullmatch(r"[a-f0-9]{64}", item.get("sha256", "")), "Duplicate or unhashed native-copy file")
            require(isinstance(item.get("bytes"), int) and item["bytes"] >= 0, "Invalid native-copy size")
            by_path[rel] = {k: item[k] for k in ("path", "bytes", "sha256")}
        return by_path
    expected = records(inventory["files"])
    proof_path = results / "REUSED_NATIVE_RECEIPT.json"
    proof = load_json(proof_path)
    require(proof.get("schema") == "satb2-hdwgcna-reused-native-receipt/v1" and proof.get("status") == "VERIFIED_EXACT_COPY"
            and proof.get("scientific_fits_run") is False, "Recovery lacks a verified no-fit copy receipt")
    require(Path(proof["source_results"]).resolve() == source and proof.get("inventory") == reuse["inventory"]
            and proof.get("native_config") == reuse["native_config"], "Native-copy proof binding differs from configuration")
    require(records(proof["copied_files"]) == expected, "Native-copy proof does not cover the complete frozen inventory")
    recovery_record = {"path": str(config_path), "bytes": config_path.stat().st_size, "sha256": digest(config_path, force=force)}
    require(proof.get("recovery_config") == recovery_record, "Copy proof is bound to a different recovery configuration")
    for key in ("all_four_completed_checkpoints_MD5_verified", "prepared_complete_SHA256_verified", "all_source_and_destination_SHA256_verified",
                "scientific_configuration_unchanged", "native_fitting_script_MD5_unchanged"):
        require(proof.get(key) is True, f"Copy proof validation did not pass: {key}")
    for root in (source, results):
        actual = set()
        for name in ("prepared", "hdwgcna"):
            tree = root / name
            require(tree.is_dir() and not tree.is_symlink(), f"Missing regular native tree: {tree}")
            for path in tree.rglob("*"):
                require(not path.is_symlink(), f"Symlink in preserved native tree: {path}")
                if path.is_file():
                    actual.add(str(path.relative_to(root)))
        require(actual == set(expected), f"Copied native tree inventory is not exhaustive: {root}")
        for rel, item in expected.items():
            path = safe_child(root, rel)
            require(path.is_file() and path.stat().st_size == item["bytes"] and digest(path, force=force) == item["sha256"],
                    f"Preserved native scientific file changed: {root}/{rel}")
    for lineage in LINEAGES:
        for network in ("pooled", "consensus"):
            folder = results / "hdwgcna" / lineage / network
            checkpoint = verify_checkpoint(folder)
            require(load_json(folder / "status.json").get("status") == "NETWORK_COMPLETE"
                    and checkpoint["fingerprint"] == load_json(folder / "run_identity.json").get("fingerprint"),
                    "Native copy has unfinished or mismatched checkpoint identity")
    script = Path(__file__).resolve().parent / "run_hdwgcna.R"
    parameters = load_json(results / "hdwgcna/run_parameters.json")
    require(digest(script, "md5", force) == parameters.get("script_md5"), "Recovered native fit R source differs from executed original MD5")
    return {"schema": "satb2-hdwgcna-native-fit-reuse-validation/v1", "source_results": str(source),
            "native_configuration": native_record, "native_configuration_sha256": native_record["sha256"],
            "native_fit_configuration": native_cfg, "recovery_configuration": recovery_record,
            "source_inventory": inventory_record, "reuse_receipt": {"path": "REUSED_NATIVE_RECEIPT.json", "bytes": proof_path.stat().st_size, "sha256": digest(proof_path, force=force)},
            "copied_file_count": len(expected), "scientific_config_keys_identical": True,
            "allowed_configuration_differences": sorted(allowed), "network_refits_run": False,
            "native_fit_R_script_md5": parameters["script_md5"], "preserved_R_script_sha256": digest(script, force=force),
            "interpretation": "Original completed native network fits and scientific tables were copied byte-identically. Only enrichment/export/transport recovery ran under the new recovery configuration."}


def package(config_path, results):
    config_path, results = Path(config_path).resolve(), Path(results).resolve()
    cfg = load_json(config_path)
    require(cfg.get("schema") == "satb2-regulatory-preparation/v1", "Wrong immutable analysis configuration schema")
    require(results == Path(cfg["run_root"]).resolve() / "results", "Results must be exactly the configured new run root / results")
    config_sha = digest(config_path)
    inputs_before = verify_input_bindings(cfg, "before")
    native_reuse = verify_native_reuse(cfg, config_path, results)
    native_config_sha = native_reuse["native_configuration_sha256"] if native_reuse else config_sha
    prepared, frozen, sample_conditions, source_coverage = verify_prepared(results, cfg)
    prepared_sources = load_json(results / "prepared/source_hashes.json")
    require(prepared_sources.get("configuration", {}).get("expected_sha256") == native_config_sha, "Prepared inputs belong to a different original native configuration")
    if native_reuse:
        require(all(prepared_sources.get("configuration", {}).get(k) == native_config_sha for k in ("expected_sha256", "before_sha256", "after_sha256")),
                "Copied prepared source configuration receipts differ from the original native fit")
    for name, item in cfg["inputs"].items():
        record = prepared_sources.get(name, {})
        require(all(record.get(k) == item["sha256"] for k in ("expected_sha256", "before_sha256", "after_sha256")), f"Prepared/source binding differs: {name}")
    workflow = load_json(results / "hdwgcna/workflow_status.json")
    require(not workflow.get("primary_failures"), "Native workflow recorded primary network failures")
    enrichment = load_json(results / "enrichment/COMPLETE.json")
    require(enrichment.get("status") == "COMPLETE" and not enrichment.get("skipped_lineages"), "Both primary lineage enrichment outputs must complete without a skipped lineage")
    payload = {"schema": SCHEMA, "analysis_id": results.parent.name, "source_fingerprint": "", "created_at": prepared.get("finished_utc"), "status": "complete", "methods": [], "endpoints": [], "samples": [{"id": s, "genotype": sample_conditions[s]} for s in sorted(sample_conditions, key=lambda s: (sample_conditions[s] != "WT", s))], "features": [], "values": [], "effects": [], "hubs": [], "nodes": [], "edges": [], "embeddings": [], "coverage": [], "downloads": [], "provenance": {"configuration_sha256": config_sha, "configuration": cfg, "bound_inputs": inputs_before, "prepared_complete_sha256": digest(results / "prepared/COMPLETE.json"), "frozen_membership_audit_sha256": digest(results / "prepared/frozen_membership_audit.tsv"), "native_workflow_status": workflow, "enrichment_status": enrichment, "no_scientific_fits_in_packager": True, "original_DESeq2_fits_unchanged": True, "strict_LR_results_unchanged": True, "TF_activity": "NOT_RUN_AWAITING_METHOD_SELECTION", "motif_supported_regulons": "NOT_RUN", "activity_UMAP": "NOT_PACKAGED_NO_VERIFIED_COORDINATE_PLUS_SCORE_AXIS", "interpretation": "Exploratory coexpression and independent-mouse continuous-score associations. They are not gene DEGs, causal TF regulation, formal module preservation or physical cell-cell communication."}}
    payload["provenance"]["native_fit_configuration_sha256"] = native_config_sha
    payload["provenance"]["native_fit_reuse"] = native_reuse
    for ep, ui_id in ENDPOINT_MAP.items():
        payload["endpoints"].append({"id": ui_id, "label": LABELS[ep], "family": "CA1" if ep.startswith("CA1_") else "IN", "source_endpoint_id": frozen[ep]["population_key"], "membership_sha256": frozen[ep]["membership_sha256"], "n_nuclei": NATIVE_COUNTS[ep], "identity_warning": "Molecular expression program; physical layer undetermined" if ep.startswith("CA1_") else "Mixed atlas identity remains under review" if ep == "IN_C8" else "Molecular reference resemblance; not connectivity or morphology"})
    branches = {}
    for network in ("pooled", "consensus"):
        for lineage in LINEAGES:
            branches[lineage + "_" + network] = append_branch(payload, cfg, results, lineage, network, sample_conditions, source_coverage)
    require(all(branches[lin + "_pooled"]["validated_complete"] for lin in LINEAGES), "Both primary lineage networks must validate")
    if native_reuse:
        for method in payload["methods"]:
            method["resources"]["native_fit_configuration_sha256"] = native_config_sha
            method["resources"]["recovery_configuration_sha256"] = config_sha
            method["resources"]["completed_native_fit_reused_without_refit"] = True
    if any(not branches[lin + "_consensus"]["validated_complete"] for lin in LINEAGES):
        payload["status"] = "partial"
    # Exact native tables remain separate downloadable artifacts, not new fits.
    # Controller logs/metadata continue changing after this task exits; they are
    # execution evidence, not immutable scientific artifacts or transport inputs.
    scientific_files = [p for p in results.rglob("*") if p.is_file() and ".snakemake" not in p.relative_to(results).parts and not p.is_relative_to(results / "atlas") and p.name not in ("FINAL_VALIDATION.json", "FINAL_VALIDATION_FAILED.json") and not p.name.endswith(".writing")]
    payload["provenance"]["mutable_controller_metadata_excluded_from_scientific_receipts"] = True
    receipts = []
    for p in sorted(scientific_files):
        record = {"path": str(p.relative_to(results)), "bytes": p.stat().st_size, "sha256": digest(p)}
        receipts.append(record)
        if p.suffix in (".tsv", ".json") or p.name.endswith(".tsv.gz"):
            payload["downloads"].append({"label": str(p.relative_to(results)), "path": "../" + str(p.relative_to(results)), "bytes": record["bytes"], "sha256": record["sha256"], "native_bytes_preserved": True})
    fingerprint_material = {"configuration_sha256": config_sha, "inputs": inputs_before, "memberships": {ep: frozen[ep]["membership_sha256"] for ep in ENDPOINT_MAP}, "native_receipts": receipts, "packager_sha256": digest(Path(__file__))}
    payload["source_fingerprint"] = hashlib.sha256(json_bytes(fingerprint_material, True)).hexdigest()
    payload["provenance"]["native_output_receipts"] = receipts
    payload["provenance"]["packager_sha256"] = fingerprint_material["packager_sha256"]
    # Membership, means and reference keys are checked before any output is written.
    require(len({(r["method_id"], r["endpoint_id"], r["feature_id"], r["sample_id"]) for r in payload["values"]}) == len(payload["values"]), "Duplicate packaged independent-mouse key")
    require(all(isinstance(r["value"], (int, float)) and math.isfinite(r["value"]) for r in payload["values"]), "Non-finite raw native mouse values")
    inputs_after = verify_input_bindings(cfg, "after")
    require(inputs_before == inputs_after and digest(config_path, force=True) == config_sha, "An original source or immutable configuration changed while packaging")
    require(verify_native_reuse(cfg, config_path, results, "after") == native_reuse, "Preserved native-fit copy/configuration proof changed while packaging")
    out = results / "atlas"
    raw = json_bytes(payload, True)
    packed = gzip.compress(raw, compresslevel=9, mtime=0)
    script = b"/* Completed native saved results; no scientific fitting in this transport. */\nwindow.SATB2_REGULATORY_DATA_GZIP=" + json.dumps(base64.b64encode(packed).decode()).encode() + b";\n"
    outputs = {"data.json": raw, "data.json.gz": packed, "data.js": script, "downloads.json": json_bytes(payload["downloads"])}
    transport = {name: {"bytes": len(content), "sha256": hashlib.sha256(content).hexdigest()} for name, content in outputs.items()}
    manifest = {"schema": "satb2-regulatory-transport/v1", "analysis_id": payload["analysis_id"], "source_fingerprint": payload["source_fingerprint"], "status": payload["status"], "branches": branches, "record_counts": {k: len(payload[k]) for k in ("methods", "endpoints", "samples", "features", "values", "effects", "hubs", "nodes", "edges", "embeddings", "coverage")}, "transport": transport, "scientific_fits_in_packager": False, "original_sources_unchanged": True}
    for name, content in outputs.items():
        write_identical(out / name, content)
    write_identical(out / "manifest.json", json_bytes(manifest))
    final = {"schema": "satb2-regulatory-final-validation/v1", "status": "COMPLETE_PRIMARY_WITH_OPTIONAL_LIMITATIONS" if payload["status"] == "partial" else "COMPLETE", "configuration_sha256": config_sha, "analysis_id": payload["analysis_id"], "source_fingerprint": payload["source_fingerprint"], "primary_networks": {lin: branches[lin + "_pooled"] for lin in LINEAGES}, "optional_consensus": {lin: branches[lin + "_consensus"] for lin in LINEAGES}, "frozen_endpoint_counts": NATIVE_COUNTS, "twelve_mice_per_endpoint_module": True, "unique_mouse_keys": True, "raw_mouse_values_finite": True, "frozen_membership_hashes_verified": True, "native_checkpoint_MD5_verified": True, "prepared_SHA256_receipts_verified": True, "all_configured_inputs_SHA256_before_after_verified": True, "original_analysis_files_unchanged": True, "TF_activity_run": False, "motif_regulon_run": False, "activity_UMAP_packaged": False, "no_scientific_fits_in_packager": True, "enrichment": enrichment, "atlas": manifest, "native_file_receipts": receipts, "package_file_receipts": [{"path": str(p.relative_to(results)), "bytes": p.stat().st_size, "sha256": digest(p)} for p in sorted(out.iterdir()) if p.is_file()]}
    final["native_fit_configuration_sha256"] = native_config_sha
    final["native_fit_reuse"] = native_reuse
    final["completed_native_fits_reused_without_refit"] = native_reuse is not None
    write_identical(results / "FINAL_VALIDATION.json", json_bytes(final))
    return final


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--results", required=True, type=Path)
    args = parser.parse_args()
    require(args.config.is_absolute() and args.results.is_absolute(), "Use absolute immutable config and results paths")
    result = package(args.config, args.results)
    print(json.dumps({"status": result["status"], "analysis_id": result["analysis_id"], "source_fingerprint": result["source_fingerprint"], "atlas": result["atlas"]["record_counts"]}, indent=2))


if __name__ == "__main__":
    main()
