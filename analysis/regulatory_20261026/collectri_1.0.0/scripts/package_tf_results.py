#!/usr/bin/env python3
"""Validate and transport completed native CollecTRI/independent-mouse results.

CLI --config ABS --results ABS_RUN/results. No scientific model, normalization,
FDR correction, source edit, motif inference or scientific subprocess is run.
Requires results/tf_activity/COMPLETE.json and results/tf_mouse/COMPLETE.json.
Config is the hash-bound satb2-tf-activity/v1 config used by both native steps;
additional_hash_inputs must bind frozen_membership_audit and mouse_coverage.
Writes only a new/strictly empty results/atlas and FINAL_VALIDATION.json last.
Original ULM, Wald enrichment, saved primary association, same-partition/LOO/
bootstrap and Satb2 no-self sensitivities remain distinct. Missing activity is
explicit null, never zero. All native tables/plots retain independent receipts.
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
METHOD = "collectri_ulm"
ENDPOINTS = {
    "CA1_Deep": ("CA1_DL", "P0244_seurat_CA1_ProS_LAYER_PROGRAM", 8750, "CA1 deep-associated expression program"),
    "CA1_Superficial": ("CA1_SL", "P0247_seurat_CA1_ProS_LAYER_PROGRAM", 11706, "CA1 superficial-associated expression program"),
    "IN_C1": ("IN_C1", "P0140_seurat_Inhibitory_CLUSTER", 1718, "Inhibitory C1"),
    "IN_C6": ("IN_C6", "P0145_seurat_Inhibitory_CLUSTER", 501, "Inhibitory C6"),
    "IN_C8": ("IN_C8", "P0147_seurat_Inhibitory_CLUSTER", 468, "Inhibitory C8 · mixed identity review"),
}
EXPRESSION = "mouse_log1p_original_sizefactor_expression"
WALD = "original_DESeq2_Wald_supplement"
REQUIRED_NATIVE = ["TF_ACTIVITY_PER_MOUSE.tsv", "TF_WALD_CONTRAST_ENRICHMENT.tsv", "TF_TARGET_COVERAGE.tsv",
                   "SATB2_NO_SELF_TARGET_ACTIVITY_PER_MOUSE.tsv", "SATB2_NO_SELF_TARGET_WALD_CONTRAST_ENRICHMENT.tsv",
                   "SATB2_NO_SELF_TARGET_COVERAGE.tsv", "PROVENANCE.json", "STARTED.json"]
REQUIRED_MOUSE = ["TF_ASSOCIATIONS.tsv", "TF_SAME_PARTITION_ASSOCIATIONS.tsv", "TF_LOO_STABILITY.tsv",
                  "TF_LOO_MODELS.tsv", "TF_MOUSE_BOOTSTRAP.tsv", "TF_MODEL_ELIGIBILITY.tsv", "TF_MOUSE_SCORES.tsv",
                  "TF_PLOT_SELECTION.tsv", "SATB2_NO_SELF_TARGET_ASSOCIATIONS.tsv"]


def require(ok, message):
    if not ok:
        raise ValueError(message)


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8*1024*1024), b""):
            h.update(block)
    return h.hexdigest()


def receipt(path, root=None):
    path = Path(path).resolve()
    return {"path": str(path.relative_to(Path(root).resolve())) if root else str(path),
            "bytes": path.stat().st_size, "sha256": digest(path)}


def raw_json(value, compact=False):
    return (json.dumps(value, ensure_ascii=False, allow_nan=False, sort_keys=True,
                       separators=(",", ":") if compact else None, indent=None if compact else 2)+"\n").encode()


def load(path):
    return json.loads(Path(path).read_text())


def table(path, columns=()):
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require(reader.fieldnames and len(reader.fieldnames) == len(set(reader.fieldnames)), f"Missing/duplicate TSV header: {path}")
        require(set(columns) <= set(reader.fieldnames), f"Missing native columns {set(columns)-set(reader.fieldnames)}: {path}")
        rows = list(reader)
        require(all(None not in row for row in rows), f"Malformed native TSV: {path}")
    return rows


def num(value, required=False):
    if value is None or str(value) in ("", "NA", "NaN", "nan", "null", "None"):
        require(not required, "Required native numeric value is missing")
        return None
    value = float(value)
    require(math.isfinite(value), "Native numeric value is nonfinite")
    return value


def integer(value):
    value = num(value, True)
    require(value >= 0 and value.is_integer(), "Native count is not a nonnegative integer")
    return int(value)


def boolean(value):
    if value in (True, "TRUE", "True", "true", "1", 1):
        return True
    require(value in (False, "FALSE", "False", "false", "0", 0), f"Native Boolean missing/invalid: {value}")
    return False


def same(a, b):
    return a == b or a is not None and b is not None and math.isclose(a, b, rel_tol=1e-13, abs_tol=1e-14)


def child(root, relative):
    root = Path(root).resolve()
    path = (root/relative).resolve()
    require(path.is_relative_to(root), f"Native receipt escapes its result root: {relative}")
    require(not (root/relative).is_symlink(), f"Symlink is not a native result: {relative}")
    return path


def check_record(record):
    require(isinstance(record, dict) and {"path", "bytes", "sha256"} <= set(record), "Input lacks a complete hash binding")
    path = Path(record["path"])
    require(path.is_absolute() and path.is_file() and not path.is_symlink(), f"Unavailable bound input: {path}")
    require(re.fullmatch(r"[a-f0-9]{64}", record["sha256"]) is not None, "Invalid input SHA256")
    actual = receipt(path)
    require(actual["bytes"] == record["bytes"] and actual["sha256"] == record["sha256"], f"Bound input changed: {path}")
    return actual


def verify_inputs(cfg):
    require(cfg.get("schema") == "satb2-tf-activity/v1", "Unexpected TF configuration schema")
    output = {}
    for name in ["samples", "references_mouse_collectri"]:
        output[name] = check_record(cfg[name])
    for name in ["source_pseudobulks", "fitted_size_factors", "existing_deseq2", "existing_deseq2_stable_id", "endpoint_mouse_metadata"]:
        require(set(cfg[name]) == set(ENDPOINTS), f"Wrong five-endpoint map: {name}")
        for endpoint, record in cfg[name].items():
            output[name+"/"+endpoint] = check_record(record)
    extra = cfg.get("additional_hash_inputs", {})
    require({"frozen_membership_audit", "mouse_coverage"} <= set(extra), "TF transport requires exact frozen membership/coverage bindings")
    for name, record in extra.items():
        output["additional/"+name] = check_record(record)
    return output


def sealed_step(folder, schema, required, cfg_record, source_path):
    require(not (folder/"FAILED.json").exists(), f"Native failure marker exists: {folder}")
    marker_path = folder/"COMPLETE.json"
    marker = load(marker_path)
    require(marker.get("schema") == schema and marker.get("status") == "COMPLETE", f"Native step is not complete: {folder}")
    require(check_record(marker["config"])["sha256"] == cfg_record["sha256"], "Native step configuration differs")
    source = check_record(marker["source"])
    require(source["sha256"] == digest(source_path), f"Executed source differs: {source_path}")
    bound = marker.get("output_files", marker.get("outputs", {}))
    require(bound and set(required) <= set(bound), f"Native step receipt does not cover required tables: {folder}")
    for relative, record in bound.items():
        path = child(folder, relative)
        require(path.is_file() and path.stat().st_size == record["bytes"] and digest(path) == record["sha256"], f"Native output changed: {path}")
    for record in marker.get("inputs", {}).values():
        check_record(record)
    # No extra unsealed output may enter the archive or UI silently.
    found = {str(p.relative_to(folder)) for p in folder.rglob("*") if p.is_file()}
    require(found == set(bound)|{"COMPLETE.json"}, f"Unsealed native output files: {found-set(bound)-{'COMPLETE.json'}}")
    return marker


def frozen_cohort(cfg):
    extra = cfg["additional_hash_inputs"]
    audit = table(extra["frozen_membership_audit"]["path"], ["endpoint", "population_key", "retained_nuclei", "membership_sha256", "native_population_membership_sha256"])
    frozen = {r["endpoint"]: r for r in audit}
    require(len(frozen) == len(audit) == 5 and set(frozen) == set(ENDPOINTS), "Frozen cohort does not have exactly five endpoints")
    for ep, r in frozen.items():
        require(r["population_key"] == ENDPOINTS[ep][1] and integer(r["retained_nuclei"]) == ENDPOINTS[ep][2], f"Frozen endpoint identity/count differs: {ep}")
        require(re.fullmatch(r"[a-f0-9]{64}", r["membership_sha256"]) and r["membership_sha256"] == r["native_population_membership_sha256"], f"Membership hashes disagree: {ep}")
    sample_rows = table(cfg["samples"]["path"], ["sample", "condition"])
    samples = {r["sample"]: r["condition"] for r in sample_rows}
    require(len(samples) == len(sample_rows) == 12 and Counter(samples.values()) == {"WT": 6, "cKO": 6}, "Mouse roster is not six WT + six cKO")
    rows = table(extra["mouse_coverage"]["path"], ["endpoint", "sample", "condition", "n_nuclei"])
    coverage = {(r["endpoint"], r["sample"]): r for r in rows}
    require(len(coverage) == len(rows) == 60 and set(coverage) == {(ep, s) for ep in ENDPOINTS for s in samples}, "Frozen coverage does not contain 60 exact endpoint-mouse rows")
    for (ep, sample), r in coverage.items():
        require(r["condition"] == samples[sample] and integer(r["n_nuclei"]) > 0, "Frozen mouse genotype/nucleus coverage differs")
    for ep in ENDPOINTS:
        require(sum(integer(coverage[ep, s]["n_nuclei"]) for s in samples) == ENDPOINTS[ep][2], "Frozen mouse counts do not sum to endpoint")
    return frozen, samples, coverage


def feature(tf):
    return METHOD+"::"+tf


def node_id(ep, tf, symbol):
    return "ct_"+hashlib.sha256(raw_json([ep, tf, symbol], True)).hexdigest()[:24]


def probability(value):
    out = num(value)
    require(out is None or 0 <= out <= 1, "Native P/FDR outside [0,1]")
    return out


def package(config_path, results):
    config_path, results = Path(config_path).resolve(), Path(results).resolve()
    cfg_record, cfg = receipt(config_path), load(config_path)
    require(results == Path(cfg["run_root"]).resolve()/"results", "Results differ from immutable configured run_root")
    require(not (results/"FINAL_VALIDATION.json").exists(), "Refusing to replace completed TF transport")
    out = results/"atlas"
    require(not out.is_symlink() and (not out.exists() or out.is_dir() and not any(out.iterdir())), "Refusing nonempty or nonregular atlas output")
    inputs_before = verify_inputs(cfg)
    source_dir = Path(__file__).resolve().parent
    executed = {"scripts/"+name: digest(source_dir/name) for name in ["run_tf_activity.py", "analyze_tf_mouse.R", "package_tf_results.py"]}
    native = results/"tf_activity"
    mouse = results/"tf_mouse"
    native_marker = sealed_step(native, "satb2-tf-activity-complete/v1", REQUIRED_NATIVE, cfg_record, source_dir/"run_tf_activity.py")
    mouse_marker = sealed_step(mouse, "satb2-tf-mouse-complete/v1", REQUIRED_MOUSE, cfg_record, source_dir/"analyze_tf_mouse.R")
    for key in ("all12_mouse_roster", "original_native_scores_preserved", "global_primary_BH", "separate_sensitivity_BH",
                "native_gene_set_FDR_not_association_FDR", "no_nucleus_pseudoreplication", "no_new_gene_DESeq2", "inputs_config_source_SHA256_unchanged"):
        require(mouse_marker.get("checks", {}).get(key) is True, f"Native mouse analysis check did not pass: {key}")
    provenance = load(native/"PROVENANCE.json")
    require(provenance["source"]["sha256"] == executed["scripts/run_tf_activity.py"] and provenance["config"]["sha256"] == cfg_record["sha256"], "Native provenance differs from executed source/config")
    require(provenance.get("genotype_association_run") is False and provenance.get("motif_regulon_run") is False and provenance.get("nucleus_activity_UMAP_run") is False, "Native ULM scope differs from contract")
    require(provenance.get("Satb2_self_target_sensitivity", {}).get("main_full_prior_unchanged") is True, "Prespecified no-self-target sensitivity is not independently recorded")
    frozen, samples, source_coverage = frozen_cohort(cfg)
    cov_rows = table(native/"TF_TARGET_COVERAGE.tsv", ["endpoint", "mode", "TF", "n_prior_targets", "n_measured_targets", "n_used_targets", "n_background_genes", "eligible_tmin5", "reason"])
    primary_cov = {(r["endpoint"], r["TF"]): r for r in cov_rows if r["mode"] == EXPRESSION}
    require(len(primary_cov) == sum(r["mode"] == EXPRESSION for r in cov_rows) and {ep for ep, _ in primary_cov} == set(ENDPOINTS), "Duplicate/incomplete native TF coverage")
    tfs = sorted({tf for _, tf in primary_cov})
    require(set(primary_cov) == {(ep, tf) for ep in ENDPOINTS for tf in tfs}, "Each frozen prior TF needs every endpoint coverage record")
    prior_rows = table(cfg["references_mouse_collectri"]["path"], ["source", "target", "weight"])
    prior = {(r["source"], r["target"]): num(r["weight"], True) for r in prior_rows}
    require(len(prior) == len(prior_rows) and set(tfs) == {s for s, _ in prior}, "Frozen prior and native TF universe differ")
    require(("Satb2", "Satb2") in prior, "The full-prior Satb2 self-target caveat is missing")
    raw_values = table(native/"TF_ACTIVITY_PER_MOUSE.tsv", ["endpoint", "endpoint_id", "source_endpoint_id", "TF", "sample", "condition", "escore", "gene_set_fdr", "n_nuclei", "n_used_targets", "n_background_genes"])
    values = {(r["endpoint"], r["TF"], r["sample"]): r for r in raw_values}
    expected_keys = {(ep, tf, s) for (ep, tf), c in primary_cov.items() if boolean(c["eligible_tmin5"]) for s in samples}
    require(len(values) == len(raw_values) and set(values) == expected_keys, "Native per-mouse TF score keys are incomplete/duplicated")
    r_scores = table(mouse/"TF_MOUSE_SCORES.tsv", ["endpoint", "TF", "sample", "escore", "gene_set_fdr", "condition", "n_nuclei"])
    r_keys = {(r["endpoint"], r["TF"], r["sample"]): r for r in r_scores}
    require(len(r_keys) == len(r_scores) and set(r_keys) == set(values), "R native score axes differ from ULM producer")
    for key, r in r_keys.items():
        original = values[key]
        require(r["condition"] == original["condition"] and integer(r["n_nuclei"]) == integer(original["n_nuclei"]), "R native mouse membership differs")
        require(same(num(r["escore"], True), num(original["escore"], True)) and same(probability(r["gene_set_fdr"]), probability(original["gene_set_fdr"])), "R changed original native score/gene-set FDR")
    payload = {"schema": SCHEMA, "analysis_id": cfg["analysis_id"], "source_fingerprint": "", "created_at": datetime.now(timezone.utc).isoformat(),
               "status": "partial", "methods": [], "endpoints": [], "samples": [{"id": s, "genotype": samples[s]} for s in sorted(samples, key=lambda s: (samples[s] != "WT", s))],
               "features": [], "values": [], "effects": [], "hubs": [], "nodes": [], "edges": [], "embeddings": [], "coverage": [], "downloads": []}
    payload["methods"] = [{"id": METHOD, "label": "Mouse CollecTRI ULM · curated TF activity", "family": "curated_tf_activity", "status": "complete",
        "score_label": "Native CollecTRI ULM t-value", "units": "ULM t-value", "activity_unit": "mouse", "endpoint_ids": [v[0] for v in ENDPOINTS.values()],
        "network_edge_basis": "Directed signed frozen mouse CollecTRI prior; widths reflect prior weights, not measured coexpression or inferred interaction strength",
        "interpretation": "Activity is inferred from signed target expression in independent mouse pseudobulks, not TF RNA itself or motif/binding evidence. Full genotype-neutral background; original saved DESeq2 size factors. Native gene-set FDR is within each mouse, distinct from original continuous-score genotype-association FDR. Technical partition partly confounds genotype. Full-prior Satb2 includes a Satb2→Satb2 self-target edge; own RNA may contribute. Prespecified Satb2-only no-self sensitivity is a separate native supplement, never the primary fit.",
        "eligibility_policy": "Exactly 12 independent mice per modelable TF × frozen endpoint; >=5 measured eligible targets. >=5 raw counts in >=3 all12 mice; ambiguous symbols excluded. Unmodelable activities stay explicit null.",
        "resources": {"mouse_CollecTRI": cfg["references_mouse_collectri"], "software": provenance["software"], "native_ULM_source": provenance["source"],
                      "same_partition_LOO_bootstrap": "Separate native sensitivity tables; primary activity/association source unchanged", "Satb2_self_target_sensitivity": provenance["Satb2_self_target_sensitivity"]}},
        {"id": "pyscenic_motif", "label": "pySCENIC motif regulons · unavailable", "family": "motif_regulon", "status": "not_available",
         "score_label": "No motif-supported regulon activity available", "units": "Unavailable", "activity_unit": "mouse",
         "endpoint_ids": [v[0] for v in ENDPOINTS.values()], "interpretation": "Not run: required frozen cisTarget motif rankings/annotations are unavailable and no cache was supplied. Curated CollecTRI activity is not substituted for motif-supported regulons.",
         "resources": {"status": "BLOCKED_MISSING_MOTIF_CACHE", "evidence": cfg.get("pyscenic", cfg.get("motif_regulon", {})),
                       "blocked_resource_manifest": cfg.get("additional_hash_inputs", {}).get("cisTarget_blocked_manifest")}}]
    for ep, (ui, source_id, count, label) in ENDPOINTS.items():
        payload["endpoints"].append({"id": ui, "label": label, "family": "CA1" if ep.startswith("CA1_") else "IN", "source_endpoint_id": source_id,
                                     "membership_sha256": frozen[ep]["membership_sha256"], "n_nuclei": count,
                                     "identity_warning": "Molecular expression program; physical layer undetermined" if ep.startswith("CA1_") else "Mixed atlas identity remains under review" if ep == "IN_C8" else "Molecular reference resemblance; not connectivity or morphology"})
    for tf in tfs:
        payload["features"].append({"id": feature(tf), "method_id": METHOD, "label": tf, "kind": "tf", "source_id": tf, "units": "Native ULM t-value", "color": "#9c65a9",
                                    "self_target_caveat": "Own RNA contributes through a curated self-target prior edge; see separate no-self sensitivity" if (tf, tf) in prior else ""})
    for (ep, tf), cov in sorted(primary_cov.items()):
        eligible = boolean(cov["eligible_tmin5"])
        cov_payload = {k: v for k, v in cov.items() if not k.endswith("targets") or k.startswith("n_")}
        for k in list(cov_payload):
            if k.startswith("n_"):
                cov_payload[k] = integer(cov_payload[k])
            elif k in ("eligible_tmin5", "self_target"):
                cov_payload[k] = boolean(cov_payload[k])
            elif k == "coverage_fraction":
                cov_payload[k] = num(cov_payload[k], True)
        payload["coverage"].append({**cov_payload, "method_id": METHOD, "endpoint_id": ENDPOINTS[ep][0], "feature_id": feature(tf),
                                     "eligible": eligible, "measured_targets": integer(cov["n_measured_targets"]), "required_targets": 5})
        for sample, genotype in samples.items():
            r = values.get((ep, tf, sample))
            nucleus_count = integer(source_coverage[ep, sample]["n_nuclei"])
            row = {"method_id": METHOD, "endpoint_id": ENDPOINTS[ep][0], "feature_id": feature(tf), "sample_id": sample,
                   "value": None, "eligible": eligible, "n_nuclei": nucleus_count, "coverage_n": integer(cov["n_used_targets"]),
                   "coverage_total": integer(cov["n_prior_targets"]), "source_fit_id": METHOD+"_expression", "missing_reason": "" if eligible else cov["reason"],
                   "normalization": provenance["normalization"], "gene_set_fdr": None}
            if r:
                require(r["endpoint_id"] == ENDPOINTS[ep][0] and r["source_endpoint_id"] == ENDPOINTS[ep][1] and r["condition"] == genotype, "Native score membership/genotype differs")
                require(integer(r["n_nuclei"]) == nucleus_count and integer(r["n_used_targets"]) == integer(cov["n_used_targets"]), "Native score coverage differs")
                row.update({**r, "value": num(r["escore"], True), "escore": num(r["escore"], True), "eligible": True,
                            "gene_set_fdr": probability(r["gene_set_fdr"]), "n_nuclei": nucleus_count,
                            "n_used_targets": integer(r["n_used_targets"]), "n_background_genes": integer(r["n_background_genes"]),
                            "genotype_association_tested": False})
            payload["values"].append(row)
    assoc = table(mouse/"TF_ASSOCIATIONS.tsv", ["endpoint", "TF", "tested", "effect_cKO_minus_WT", "statistic", "pvalue", "association_fdr", "n_WT", "n_cKO", "fdr_family"])
    seen_effects = set()
    for r in assoc:
        ep, tf = r["endpoint"], r["TF"]
        require((ep, tf) in primary_cov and (ep, tf) not in seen_effects, "Unknown/duplicate native TF association")
        seen_effects.add((ep, tf))
        tested = boolean(r["tested"])
        require(not tested or boolean(primary_cov[ep, tf]["eligible_tmin5"]), "Association tested an unavailable activity")
        require(not tested or (integer(r["n_WT"]), integer(r["n_cKO"])) == (6, 6), "Primary TF association changed six-mouse/genotype cohort")
        effect = {**r, "method_id": METHOD, "endpoint_id": ENDPOINTS[ep][0], "feature_id": feature(tf), "contrast": "cKO_vs_WT",
                  "effect_name": "Mean native ULM activity difference (cKO − WT)", "effect": num(r["effect_cKO_minus_WT"], tested),
                  "statistic": num(r["statistic"], tested), "pvalue": probability(r["pvalue"]), "padj": probability(r["association_fdr"]),
                  "association_fdr": probability(r["association_fdr"]), "tested": tested, "n_WT": integer(r["n_WT"]), "n_cKO": integer(r["n_cKO"]),
                  "source_fit_id": "native_mouse_activity_lm_primary", "source_table": "tf_mouse/TF_ASSOCIATIONS.tsv"}
        for col in ["mean_WT", "mean_cKO", "effect_cKO_minus_WT", "SE", "CI_low", "CI_high"]:
            if col in r:
                effect[col] = num(r[col])
        require(not tested or effect["pvalue"] is not None and effect["padj"] is not None and r["fdr_family"], "Tested association lacks native P/FDR family")
        payload["effects"].append(effect)
    require({(ep, tf) for ep, tf, _ in expected_keys} <= seen_effects, "Modelable native activities lack primary association status")
    nodes = {}
    for ep in ENDPOINTS:
        edges = table(native/ep/"mouse_activity_used_signed_prior_edges.tsv", ["source", "target", "weight", "self_target", "edge_basis"])
        edge_seen, counts = set(), Counter()
        for r in edges:
            tf, target = r["source"], r["target"]
            require((tf, target) not in edge_seen and (tf, target) in prior and (ep, tf) in primary_cov, "Unexpected/duplicate native prior edge")
            edge_seen.add((tf, target)); counts[tf] += 1
            weight = num(r["weight"], True)
            require(weight == prior[tf, target] and boolean(r["self_target"]) == (tf == target), "Curated prior weight/self flag changed")
            require(boolean(primary_cov[ep, tf]["eligible_tmin5"]), "Used edge belongs to unavailable TF")
            for symbol, kind in [(tf, "tf"), (target, "target")]:
                ident = node_id(ep, tf, symbol)
                n = {"id": ident, "method_id": METHOD, "endpoint_id": ENDPOINTS[ep][0], "feature_id": feature(tf),
                     "label": symbol, "source_id": symbol, "kind": "tf" if symbol == tf else kind, "self_target": symbol == tf and (tf, tf) in edge_seen}
                require(ident not in nodes or all(nodes[ident][k] == n[k] for k in ["method_id", "endpoint_id", "feature_id", "label"]), "Network node ID collision")
                nodes[ident] = n
            payload["edges"].append({"method_id": METHOD, "endpoint_id": ENDPOINTS[ep][0], "feature_id": feature(tf),
                                      "source": node_id(ep, tf, tf), "target": node_id(ep, tf, target), "weight": weight,
                                      "sign": "positive" if weight > 0 else "negative", "directed": True, "self_target": tf == target,
                                      "edge_basis": "curated_mouse_CollecTRI_signed_prior_not_measured_strength_not_motif_inference"})
        require(all(counts[tf] == integer(cov["n_used_targets"]) for (endpoint, tf), cov in primary_cov.items() if endpoint == ep), "Used-prior edge totals disagree with target coverage")
    payload["nodes"] = list(nodes.values())
    wald = table(native/"TF_WALD_CONTRAST_ENRICHMENT.tsv", ["endpoint", "TF", "escore", "gene_set_fdr"])
    wald_cov = {(r["endpoint"], r["TF"]): r for r in cov_rows if r["mode"] == WALD}
    require(len({(r["endpoint"], r["TF"]) for r in wald}) == len(wald), "Duplicate native Wald supplement rows")
    require({(r["endpoint"], r["TF"]) for r in wald} == {k for k, r in wald_cov.items() if boolean(r["eligible_tmin5"])}, "Native finite-Wald supplement coverage differs")
    payload["wald_enrichment"] = [{**r, "escore": num(r["escore"], True), "gene_set_fdr": probability(r["gene_set_fdr"]),
                                  "genotype_association_tested": False} for r in wald]
    # Supplement tables stay native downloads; no averaging/retesting or fake
    # motif/UMAP values are introduced by the adapter.
    native_files = []
    for folder in (native, mouse):
        for path in sorted(folder.rglob("*")):
            require(not path.is_symlink(), "Symlink in native result tree")
            if path.is_file():
                native_files.append(receipt(path, results))
    require(native_files, "No sealed native outputs")
    for r in native_files:
        if r["path"].endswith((".tsv", ".tsv.gz", ".json")):
            payload["downloads"].append({**r, "native_bytes_preserved": True, "label": r["path"]})
    payload["native_supplements"] = {"wald_enrichment": {"table": "tf_activity/TF_WALD_CONTRAST_ENRICHMENT.tsv", "interpretation": "Native TF enrichment of original gene Wald statistics; not independent-mouse genotype testing"},
        "Satb2_no_self_target": {"tables": [r for r in native_files if "SATB2_NO_SELF_TARGET" in r["path"] or "/self_excluded/" in r["path"]], "primary_full_prior_unchanged": True,
                                 "caveat": "Own RNA may contribute through the primary Satb2→Satb2 prior link; sensitivity excludes only that link, with unchanged background genes and separate association BH family"},
        "mouse_sensitivities": [r for r in native_files if any(s in r["path"] for s in ("SAME_PARTITION", "LOO_", "BOOTSTRAP"))]}
    inputs_after = verify_inputs(cfg)
    require(inputs_before == inputs_after and receipt(config_path) == cfg_record, "Original source/config changed while packaging")
    for r in native_files:
        require(receipt(child(results, r["path"]), results) == r, "Native output changed while packaging")
    require(executed == {name: digest(source_dir/Path(name).name) for name in executed}, "TF workflow source changed while packaging")
    fingerprint = {"configuration_sha256": cfg_record["sha256"], "inputs": inputs_before, "native_outputs": native_files,
                   "memberships": {ep: frozen[ep]["membership_sha256"] for ep in ENDPOINTS}, "executed_source_sha256": executed}
    payload["source_fingerprint"] = hashlib.sha256(raw_json(fingerprint, True)).hexdigest()
    payload["provenance"] = {"configuration_sha256": cfg_record["sha256"], "configuration": cfg, "bound_inputs": inputs_before,
                              "executed_source_sha256": executed, "packager_sha256": digest(__file__), "native_output_receipts": native_files,
                              "native_ULM_provenance": provenance, "native_mouse_association_complete": mouse_marker,
                              "no_scientific_fits_in_packager": True, "TF_activity": "NATIVE_COLLECTRI_ULM_COMPLETE",
                              "motif_supported_regulons": "BLOCKED_MISSING_MOTIF_CACHE_NOT_RUN", "activity_UMAP": "NOT_PACKAGED_NO_NATIVE_NUCLEUS_ACTIVITY",
                              "native_gene_set_FDR_is_not_genotype_association_FDR": True,
                              "original_DESeq2_and_LR_results_unchanged": True}
    raw = raw_json(payload, True)
    compressed = gzip.compress(raw, compresslevel=9, mtime=0)
    outputs = {"data.json": raw, "data.json.gz": compressed,
               "data.js": b"/* Completed native curated TF activity; no fits in this transport. */\nwindow.SATB2_REGULATORY_DATA_GZIP="+json.dumps(base64.b64encode(compressed).decode()).encode()+b";\n",
               "downloads.json": raw_json(payload["downloads"])}
    transport = {name: {"bytes": len(content), "sha256": hashlib.sha256(content).hexdigest()} for name, content in outputs.items()}
    manifest = {"schema": "satb2-regulatory-transport/v1", "analysis_id": cfg["analysis_id"], "source_fingerprint": payload["source_fingerprint"],
                "status": "partial", "completed_primary": "curated_mouse_CollecTRI_ULM_and_independent_mouse_associations", "blocked_optional": "pySCENIC_motif_cache_missing",
                "record_counts": {k: len(payload[k]) for k in ("methods", "endpoints", "samples", "features", "values", "effects", "nodes", "edges", "coverage", "embeddings")},
                "transport": transport, "scientific_fits_in_packager": False, "original_sources_unchanged": True}
    out.mkdir(parents=True, exist_ok=True)
    for name, content in {**outputs, "manifest.json": raw_json(manifest)}.items():
        with (out/name).open("xb") as handle:
            handle.write(content)
    final = {"schema": "satb2-tf-regulatory-final-validation/v1", "status": "COMPLETE_WITH_BLOCKED_MOTIF",
             "analysis_id": cfg["analysis_id"], "configuration_sha256": cfg_record["sha256"], "source_fingerprint": payload["source_fingerprint"],
             "executed_source_sha256": executed, "TF_activity_run": True, "motif_regulon_run": False, "activity_UMAP_packaged": False,
             "twelve_mice_per_endpoint_TF": True, "unique_mouse_keys": True, "raw_mouse_values_finite": True,
             "frozen_membership_hashes_verified": True, "native_score_receipts_verified": True,
             "all_configured_inputs_SHA256_before_after_verified": True, "original_analysis_files_unchanged": True, "no_scientific_fits_in_packager": True,
             "frozen_endpoint_counts": {ep: info[2] for ep, info in ENDPOINTS.items()}, "native_file_receipts": native_files,
             "package_file_receipts": [receipt(p, results) for p in sorted(out.iterdir()) if p.is_file()], "atlas": manifest}
    with (results/"FINAL_VALIDATION.json").open("xb") as handle:
        handle.write(raw_json(final))
    return final


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--results", required=True, type=Path)
    args = parser.parse_args()
    require(args.config.is_absolute() and args.results.is_absolute(), "Use absolute immutable config and results paths")
    final = package(args.config, args.results)
    print(json.dumps({"status": final["status"], "analysis_id": final["analysis_id"], "atlas": final["atlas"]["record_counts"],
                      "gzip_bytes": final["atlas"]["transport"]["data.json.gz"]["bytes"]}, indent=2))


if __name__ == "__main__":
    main()
