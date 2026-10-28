#!/usr/bin/env python3
"""Native mouse CollecTRI ULM activity; no gene-DE or genotype-association fit.

Run only through the reviewed workflow: python run_tf_activity.py --config ABS
--outdir ABS_FRESH_DIRECTORY. A strictly empty real directory made by the
workflow controller is accepted; files, symlinks and nonempty dirs are refused.

JSON contract, schema ``satb2-tf-activity/v1``:
  analysis_id: nonempty string; decoupler_version: "2.3.0";
  samples: FILE; references_mouse_collectri: FILE;
  source_pseudobulks: {ENDPOINT: FILE}; fitted_size_factors: {ENDPOINT: FILE};
  existing_deseq2: {ENDPOINT: FILE};
  existing_deseq2_stable_id: {ENDPOINT: FILE};
  endpoint_mouse_metadata: {ENDPOINT: FILE};
  expected_gene_rows: 32285; expected_group_counts: {ENDPOINT: positive int};
  tmin: 5; min_gene_count: 5; min_gene_mice: 3;
  satb2_no_self_target_sensitivity: true;
  plot_top_tfs: 30; boxplot_top_tfs: 6; plots: true;
  optional additional_hash_inputs: {name: FILE}.
FILE = {path: absolute path, bytes: integer, sha256: lowercase hex64}.
ENDPOINT keys must be exactly CA1_Deep, CA1_Superficial, IN_C1, IN_C6, IN_C8.
Samples: unique 12 rows, sample/condition, six WT and six cKO. Counts: unique
gene_id/gene_symbol followed by exactly those 12 sample columns. Size factors:
sample/fitted_size_factor from the ORIGINAL matching DESeq2 fit. Metadata:
sample/condition/n_nuclei/total_UMIs (matching original MICE_ALL12_AUDIT.tsv).
ALL_RESULTS.txt: native seven columns gene/baseMean/log2FoldChange/lfcSE/stat/
pvalue/padj. Its original ALL_RESULTS_STABLE_ID_AUDIT.tsv must match every
numeric result and exported display symbol. The audit maps stable gene IDs to
source_gene_symbol. Prior TSV: source/target/weight, exact mouse symbols and
unique pairs, signed finite nonzero weights. It is not downloaded by this script.

Expression is log1p(symbol raw count sum / saved original size factor), using
all genotype-neutral eligible symbols (>=5 raw counts in >=3 of all12 mice).
All32k raw features remain represented in normalization-source checks; factors
are NEVER re-estimated after symbol/target/gene filtering. Missing symbols and
symbols shared by multiple stable IDs are excluded (no arbitrary ID choice).
No alias inference, case conversion, complex splitting, DEG selection, SCT,
Harmony, per-gene scaling, genotype-specific filter or normalization is applied.

Supplementary contrast ULM uses every finite unambiguous original DESeq2 Wald
statistic, including zeros and non-DEGs, without expression transformation.
ULM escore is the native t-value. Native gene_set_fdr is BH across modeled TFs
WITHIN each mouse or the single contrast; it is NOT genotype-association FDR.
No WT-cKO activity test is performed here. CollecTRI edges are signed curated
prior associations, NOT inferred motif regulons or demonstrated direct binding.

Output: native score and FDR matrices per endpoint/mode; all-TF target coverage,
used signed prior edges, all-gene symbol/filter audit, mouse/sample audit,
long mouse activities and Wald enrichments, genotype-neutral descriptive plots,
software/input/output SHA256 receipts and COMPLETE.json (written last).
No nucleus UMAP or motif/regulon inference is performed.
The main prior is unchanged, including any self-target edge. A prespecified
separate Satb2-only sensitivity excludes only Satb2 -> Satb2, retaining all
background genes and tmin5. Native sensitivity gene-set FDR is within the single
Satb2 set, not the primary across-TF FDR or a genotype P-value.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import inspect
import json
import math
import os
import platform
import re
import sys
import traceback
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


ENDPOINTS = {
    "CA1_Deep": ("CA1_DL", "P0244_seurat_CA1_ProS_LAYER_PROGRAM", 8750),
    "CA1_Superficial": ("CA1_SL", "P0247_seurat_CA1_ProS_LAYER_PROGRAM", 11706),
    "IN_C1": ("IN_C1", "P0140_seurat_Inhibitory_CLUSTER", 1718),
    "IN_C6": ("IN_C6", "P0145_seurat_Inhibitory_CLUSTER", 501),
    "IN_C8": ("IN_C8", "P0147_seurat_Inhibitory_CLUSTER", 468),
}
NUMERIC_RESULTS = ["baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]
INPUT_MAPS = ["source_pseudobulks", "fitted_size_factors", "existing_deseq2",
              "existing_deseq2_stable_id", "endpoint_mouse_metadata"]
NA_SYMBOLS = {"", "NA", "N/A", "nan", "None", "null"}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def receipt(path):
    path = Path(path).resolve()
    return {"path": str(path), "bytes": path.stat().st_size, "sha256": sha256(path)}


def json_write(value, path):
    with Path(path).open("x", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, allow_nan=False)
        handle.write("\n")


def table_write(frame, path):
    require(not Path(path).exists(), f"Refusing to replace output: {path}")
    frame.to_csv(path, sep="\t", index=False, na_rep="NA", float_format="%.17g")


def table_read(path, required):
    frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    require(set(required) <= set(frame), f"Missing columns {set(required)-set(frame)} in {path}")
    require(not frame.columns.duplicated().any(), f"Duplicate header in {path}")
    return frame


def numbers(values, label, finite=True):
    out = pd.to_numeric(values.replace({"": np.nan, "NA": np.nan, "NaN": np.nan,
                                         "Inf": np.inf, "-Inf": -np.inf}), errors="raise")
    if finite:
        require(np.isfinite(out.to_numpy(dtype=float)).all(), f"Nonfinite {label}")
    return out.astype(float)


def check_record(value, name):
    require(isinstance(value, dict) and {"path", "bytes", "sha256"} <= set(value), f"Unbound input {name}")
    path = Path(value["path"])
    require(path.is_absolute() and path.is_file() and not path.is_symlink(), f"Invalid input path {name}")
    require(re.fullmatch(r"[0-9a-f]{64}", value["sha256"]) is not None, f"Invalid SHA256 {name}")
    observed = receipt(path)
    require(observed["bytes"] == value["bytes"] and observed["sha256"] == value["sha256"], f"Changed input {name}")
    return observed


def fresh_output(path):
    path = Path(path)
    require(not path.is_symlink(), "Refusing a symlink output directory")
    require(not path.exists() or path.is_dir() and not any(path.iterdir()),
            "Refusing to reuse/overwrite a file or nonempty output directory")


def verify_config_inputs(cfg):
    require(cfg.get("schema") == "satb2-tf-activity/v1", "Unsupported TF activity config schema")
    require(isinstance(cfg.get("analysis_id"), str) and cfg["analysis_id"].strip(), "Missing analysis_id")
    require(cfg.get("decoupler_version", "2.3.0") == "2.3.0", "This API contract requires decoupler 2.3.0")
    require(cfg.get("tmin", 5) == 5 and cfg.get("min_gene_count", 5) == 5 and cfg.get("min_gene_mice", 3) == 3,
            "This protocol fixes tmin=5 and raw count>=5 in>=3 all12 mice")
    require(cfg.get("satb2_no_self_target_sensitivity", cfg.get("satb2_self_excluded_sensitivity", True)) is True,
            "Prespecified Satb2 self-target exclusion sensitivity must remain separate and enabled")
    records = {}
    for name in ["samples", "references_mouse_collectri"]:
        records[name] = check_record(cfg[name], name)
    for group in INPUT_MAPS:
        require(set(cfg[group]) == set(ENDPOINTS), f"Exactly five endpoints required for {group}")
        for endpoint, value in cfg[group].items():
            if "source_endpoint_id" in value:
                require(value["source_endpoint_id"] == ENDPOINTS[endpoint][1], f"Wrong source endpoint {group}/{endpoint}")
            records[f"{group}/{endpoint}"] = check_record(value, f"{group}/{endpoint}")
    for name, value in cfg.get("additional_hash_inputs", {}).items():
        records[f"additional/{name}"] = check_record(value, name)
    return records


def sample_roster(cfg):
    frame = table_read(cfg["samples"]["path"], ["sample", "condition"])
    require(len(frame) == 12 and not frame["sample"].duplicated().any() and frame["sample"].str.len().gt(0).all(),
            "Samples must be 12 unique independent mice")
    require(frame["condition"].value_counts().to_dict() == {"WT": 6, "cKO": 6}, "Expected six WT and six cKO mice")
    return frame


def prior_network(cfg):
    frame = table_read(cfg["references_mouse_collectri"]["path"], ["source", "target", "weight"])
    require(len(frame) > 0 and not frame.duplicated(["source", "target"]).any(), "Prior is empty or has repeated source-target pairs")
    for col in ["source", "target"]:
        require(frame[col].map(lambda x: x not in NA_SYMBOLS and x == x.strip()).all(), f"Missing/whitespace prior {col}")
    frame["weight"] = numbers(frame["weight"], "prior weight")
    require(frame["weight"].ne(0).all(), "Zero prior weights are unsupported")
    return frame


def endpoint_inputs(cfg, endpoint, roster):
    samples = roster["sample"].tolist()
    raw = table_read(cfg["source_pseudobulks"][endpoint]["path"], ["gene_id", "gene_symbol"] + samples)
    require(set(raw) == {"gene_id", "gene_symbol", *samples}, f"Raw counts contain an unexpected sample/column: {endpoint}")
    require(len(raw) == int(cfg.get("expected_gene_rows", 32285)), f"Wrong full gene count: {endpoint}")
    require(raw["gene_id"].str.len().gt(0).all() and not raw["gene_id"].duplicated().any(), f"Stable IDs not unique: {endpoint}")
    counts = raw[samples].apply(lambda s: numbers(s, f"raw counts {endpoint}"))
    values = counts.to_numpy()
    require((values >= 0).all() and np.equal(values, np.floor(values)).all() and (values <= np.iinfo(np.int32).max).all(), f"Raw counts are not exact native nonnegative R integers: {endpoint}")
    counts = counts.astype(np.int64)
    require(counts.sum(axis=0).gt(0).all(), f"Zero full-gene mouse library: {endpoint}")
    sf = table_read(cfg["fitted_size_factors"][endpoint]["path"], ["sample", "fitted_size_factor"])
    require(len(sf) == 12 and not sf["sample"].duplicated().any() and set(sf["sample"]) == set(samples), f"Size factors not exactly all12 mice: {endpoint}")
    sf = sf.set_index("sample").loc[samples]
    sf["fitted_size_factor"] = numbers(sf["fitted_size_factor"], f"saved size factors {endpoint}")
    require(sf["fitted_size_factor"].gt(0).all(), f"Nonpositive original size factor: {endpoint}")
    md = table_read(cfg["endpoint_mouse_metadata"][endpoint]["path"], ["sample", "condition", "n_nuclei", "total_UMIs"])
    require(len(md) == 12 and not md["sample"].duplicated().any() and set(md["sample"]) == set(samples), f"Mouse endpoint coverage differs: {endpoint}")
    md = md.set_index("sample").loc[samples].reset_index()
    require(md["condition"].tolist() == roster["condition"].tolist(), f"Genotype mismatch: {endpoint}")
    for col in ["n_nuclei", "total_UMIs"]:
        md[col] = numbers(md[col], col)
        require(md[col].gt(0).all() and np.equal(md[col], np.floor(md[col])).all(), f"Nonpositive/noninteger {col}: {endpoint}")
        md[col] = md[col].astype(np.int64)
    require(np.array_equal(md["total_UMIs"].to_numpy(), counts.sum(axis=0).to_numpy()), f"Original complete count-library sum differs: {endpoint}")
    expected = cfg.get("expected_group_counts", {}).get(endpoint, ENDPOINTS[endpoint][2])
    require(int(md["n_nuclei"].sum()) == int(expected), f"Frozen endpoint nucleus count differs: {endpoint}")
    if "included_in_model" in md:
        require(md["included_in_model"].isin(["TRUE", "True", "true", "1"]).all(), f"Original PRIMARY model does not include all12 mice: {endpoint}")
    readable = table_read(cfg["existing_deseq2"][endpoint]["path"], ["gene"] + NUMERIC_RESULTS)
    require(readable.columns.tolist() == ["gene"] + NUMERIC_RESULTS, f"Original seven-column ALL_RESULTS schema differs: {endpoint}")
    audit = table_read(cfg["existing_deseq2_stable_id"][endpoint]["path"], ["gene_id", "source_gene_symbol", "display_symbol"] + NUMERIC_RESULTS)
    require(len(audit) == len(readable) and not audit["gene_id"].duplicated().any(), f"DESeq2 stable audit row mismatch: {endpoint}")
    require(audit["display_symbol"].tolist() == readable["gene"].tolist(), f"DESeq2 export symbol mismatch: {endpoint}")
    raw_symbols = raw.set_index("gene_id")["gene_symbol"]
    require(audit["gene_id"].isin(raw_symbols.index).all(), f"Wald stable ID absent from raw counts: {endpoint}")
    require(audit["source_gene_symbol"].tolist() == raw_symbols.loc[audit["gene_id"]].tolist(), f"DESeq2/raw symbol mapping differs: {endpoint}")
    for col in NUMERIC_RESULTS:
        a, b = numbers(audit[col], col, finite=False), numbers(readable[col], col, finite=False)
        require(np.array_equal(a.to_numpy(), b.to_numpy(), equal_nan=True), f"ALL_RESULTS numeric audit differs for {col}: {endpoint}")
        audit[col] = a
    return raw[["gene_id", "gene_symbol"]], counts, sf["fitted_size_factor"], md, audit


def prepare_expression(mapping, counts, sf):
    symbols = mapping["gene_symbol"]
    valid = symbols.map(lambda s: s not in NA_SYMBOLS and s == s.strip())
    ambiguous = valid & symbols.duplicated(keep=False)
    symbol_ok = valid & ~ambiguous
    gene_audit = mapping.copy()
    gene_audit["symbol_missing_or_invalid"] = ~valid
    gene_audit["symbol_ambiguous_multiple_stable_ids"] = ambiguous
    gene_audit["n_mice_raw_count_ge5"] = counts.ge(5).sum(axis=1)
    gene_audit["all_zero_all12"] = counts.sum(axis=1).eq(0)
    eligible = symbol_ok & gene_audit["n_mice_raw_count_ge5"].ge(3) & ~gene_audit["all_zero_all12"]
    gene_audit["eligible_activity_background"] = eligible
    # No duplicated symbol survives; explicit groupby documents symbol aggregation.
    summed = counts.loc[symbol_ok].groupby(symbols.loc[symbol_ok], sort=True).sum()
    keep = summed.ge(5).sum(axis=1).ge(3) & summed.sum(axis=1).gt(0)
    summed = summed.loc[keep]
    require(set(summed.index) == set(symbols.loc[eligible]), "Symbol aggregation changed eligible background")
    normalized = summed.div(sf, axis=1)
    expression = np.log1p(normalized).T
    require(expression.shape[1] >= 6 and np.isfinite(expression.to_numpy()).all(), "Too few or nonfinite background genes")
    require(expression.std(axis=1, ddof=1).gt(0).all(), "Constant mouse expression background")
    return expression, gene_audit


def prepare_wald(mapping, audit):
    valid = mapping["gene_symbol"].map(lambda s: s not in NA_SYMBOLS and s == s.strip())
    unambiguous = set(mapping.loc[valid & ~mapping["gene_symbol"].duplicated(keep=False), "gene_symbol"])
    good = audit["source_gene_symbol"].isin(unambiguous) & np.isfinite(audit["stat"])
    selected = audit.loc[good, ["gene_id", "source_gene_symbol", "stat"]].copy()
    require(not selected["source_gene_symbol"].duplicated().any(), "Repeated Wald-statistic symbol")
    require(len(selected) >= 6 and selected["stat"].std(ddof=1) > 0, "Too few or constant finite Wald statistics")
    stats = pd.DataFrame([selected.set_index("source_gene_symbol")["stat"]], index=["cKO_vs_WT_Wald"])
    wald_audit = audit[["gene_id", "source_gene_symbol", "stat"]].copy()
    wald_audit["eligible_finite_unambiguous_wald_background"] = good
    return stats, wald_audit


def coverage_and_edges(net, background, endpoint, mode, assayed_genes=None, measured_genes=None):
    present = set(background)
    assayed = set(assayed_genes) if assayed_genes is not None else present
    measured = set(measured_genes) if measured_genes is not None else present
    rows, used = [], []
    for tf, part in net.groupby("source", sort=True):
        overlap = part.loc[part["target"].isin(present)]
        eligible = len(overlap) >= 5
        rows.append({"endpoint": endpoint, "endpoint_id": ENDPOINTS[endpoint][0],
                     "source_endpoint_id": ENDPOINTS[endpoint][1], "mode": mode, "TF": tf,
                     "n_prior_targets": len(part), "n_assayed_unambiguous_targets": int(part["target"].isin(assayed).sum()),
                     "n_prior_self_targets": int(part["source"].eq(part["target"]).sum()),
                     "n_measured_targets": int(part["target"].isin(measured).sum()), "n_eligible_background_targets": len(overlap),
                     "n_used_targets": len(overlap) if eligible else 0, "n_background_genes": len(present),
                     "n_used_self_targets": int(overlap["source"].eq(overlap["target"]).sum()) if eligible else 0,
                     "self_target": bool(part["source"].eq(part["target"]).any()),
                     "eligible_tmin5": eligible, "coverage_fraction": len(overlap)/len(part),
                     "prior_targets": ";".join(part["target"]),
                     "measured_targets": ";".join(part.loc[part["target"].isin(measured), "target"]),
                     "eligible_background_targets": ";".join(overlap["target"]),
                     "used_targets": ";".join(overlap["target"]) if eligible else "",
                     "reason": "MODELED" if eligible else "FEWER_THAN_5_ELIGIBLE_TARGETS"})
        if eligible:
            used.append(overlap)
    edges = pd.concat(used, ignore_index=True) if used else net.iloc[:0].copy()
    edges["self_target"] = edges["source"].eq(edges["target"])
    return pd.DataFrame(rows), edges


def run_native_ulm(matrix, network, dc):
    """Use full matrix background; never preselect only the prior target columns."""
    require(matrix.columns.is_unique and matrix.index.is_unique, "ULM matrix axes not unique")
    require(len(network) and matrix.shape[1] > 2, "No modelable prior/background")
    result = dc.mt.ulm(data=matrix, net=network[["source", "target", "weight"]], tmin=5,
                       raw=False, empty=False, verbose=True, tval=True)
    require(isinstance(result, tuple) and len(result) == 2, "Unexpected native decoupler ULM return")
    es, fdr = result
    require(es.index.equals(matrix.index) and fdr.index.equals(es.index) and fdr.columns.equals(es.columns), "ULM output axes changed")
    require(set(es.columns) == set(network["source"]), "ULM silently changed modelable TF coverage")
    require(np.isfinite(es.to_numpy()).all() and np.isfinite(fdr.to_numpy()).all(), "Native ULM returned nonfinite estimate/FDR")
    require(((fdr.to_numpy() >= 0) & (fdr.to_numpy() <= 1)).all(), "Native gene-set FDR outside [0,1]")
    return es, fdr


def long_scores(es, fdr, endpoint, mode, coverage, md=None):
    cov = coverage.set_index("TF")
    mouse = md.set_index("sample") if md is not None else None
    rows = []
    for sample in es.index:
        for tf in es.columns:
            row = {"endpoint": endpoint, "endpoint_id": ENDPOINTS[endpoint][0],
                   "source_endpoint_id": ENDPOINTS[endpoint][1], "mode": mode, "TF": tf,
                   "sample": sample, "escore": float(es.loc[sample, tf]),
                   "gene_set_fdr": float(fdr.loc[sample, tf]),
                   "gene_set_fdr_family": "BH_across_TFs_within_one_input_mouse_or_contrast",
                   "n_used_targets": int(cov.loc[tf, "n_used_targets"]),
                   "n_background_genes": int(cov.loc[tf, "n_background_genes"]),
                   "score_unit": "native_ULM_t_value", "genotype_association_tested": False}
            if mouse is not None:
                row.update(condition=str(mouse.loc[sample, "condition"]), n_nuclei=int(mouse.loc[sample, "n_nuclei"]),
                           complete_gene_UMIs=int(mouse.loc[sample, "total_UMIs"]), eligible=True)
            else:
                row.update(contrast="cKO_vs_WT", interpretation="target_set_enrichment_of_original_gene_Wald_statistics")
            rows.append(row)
    return pd.DataFrame(rows)


def plots(es, md, endpoint, folder, cfg):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    n_heat = min(max(int(cfg.get("plot_top_tfs", 30)), 1), 60, len(es.columns))
    n_box = min(max(int(cfg.get("boxplot_top_tfs", 6)), 1), 24, len(es.columns))
    # Selection does not use genotype effects, genotype P-values or DEGs.
    order = es.var(axis=0, ddof=1).sort_values(ascending=False, kind="stable").index
    selected = order[:n_heat]
    sd = es[selected].std(axis=0, ddof=1)
    z = es[selected].sub(es[selected].mean()).div(sd.replace(0, np.nan)).T
    fig, ax = plt.subplots(figsize=(max(8, len(es)*.55), max(4, n_heat*.24)))
    masked = np.ma.masked_invalid(z.to_numpy())
    image = ax.imshow(masked, cmap="RdBu_r", vmin=-2.5, vmax=2.5, aspect="auto", interpolation="none")
    ax.set_xticks(range(len(es)), [f"{s}\n{md.set_index('sample').loc[s,'condition']}" for s in es.index], rotation=60, ha="right", fontsize=8)
    ax.set_yticks(range(n_heat), selected, fontsize=8)
    ax.set_title(f"{endpoint}: CollecTRI ULM activity\nDescriptive within-TF z across all12 mice; top by pooled variance")
    fig.colorbar(image, ax=ax, label="Display z (sample SD); native score tables unchanged")
    fig.tight_layout()
    for ext in ["png", "pdf", "svg"]:
        fig.savefig(folder / f"activity_heatmap_display_z.{ext}", dpi=160)
    plt.close(fig)
    table_write(pd.DataFrame({"TF": selected, "selection_rank_by_genotype_neutral_variance": np.arange(1, n_heat+1),
                              "pooled_mouse_mean": es[selected].mean().to_numpy(), "pooled_mouse_sample_sd": sd.to_numpy()}), folder / "plot_TF_selection.tsv")
    nc, nr = min(n_box, 3), math.ceil(n_box/3)
    fig, axes = plt.subplots(nr, nc, figsize=(nc*3.5, nr*3.2), squeeze=False)
    conditions = md.set_index("sample")["condition"]
    for ax, tf in zip(axes.flat, order[:n_box]):
        for position, condition, color in [(0,"WT","#397eb5"),(1,"cKO","#cf8334")]:
            subset = es.loc[conditions.eq(condition), tf]
            box = ax.boxplot([subset.to_numpy()], positions=[position], widths=.45, patch_artist=True, showfliers=False)
            box["boxes"][0].set_facecolor(color+"55")
            ax.scatter(position+np.linspace(-.11,.11,len(subset)), subset, color=color, s=27, edgecolors="white", linewidths=.4, zorder=3)
        ax.set_xticks([0,1],["WT (6 mice)","cKO (6 mice)"])
        ax.set_title(tf); ax.set_ylabel("Native ULM t-value")
    for ax in list(axes.flat)[n_box:]:
        ax.set_visible(False)
    fig.suptitle(f"{endpoint}: independent mouse activities\nTop TFs by pooled variance; no genotype-association P-values")
    fig.tight_layout()
    for ext in ["png", "pdf", "svg"]:
        fig.savefig(folder / f"mouse_box_dot_native_units.{ext}", dpi=160)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    require(args.config.is_absolute() and args.outdir.is_absolute(), "Use absolute config and output paths")
    require(args.config.is_file() and not args.config.is_symlink(), "Config must be a regular file")
    fresh_output(args.outdir)
    cfg_receipt = receipt(args.config)
    cfg = json.loads(args.config.read_text())
    before = verify_config_inputs(cfg)
    source_receipt = receipt(__file__)
    require(importlib.metadata.version("decoupler") == "2.3.0", "Unexpected decoupler version; re-review API")
    import decoupler as dc
    roster, net = sample_roster(cfg), prior_network(cfg)
    fresh_output(args.outdir)
    args.outdir.mkdir(parents=True, exist_ok=True)
    out = args.outdir
    json_write({"status": "RUNNING", "started_UTC": datetime.now(timezone.utc).isoformat(), "config": cfg_receipt,
                "source": source_receipt, "inputs": before}, out / "STARTED.json")
    all_activity, all_contrast, all_coverage, summaries = [], [], [], []
    self_activity, self_contrast, self_coverage, self_status = [], [], [], {}
    satb2 = net.loc[net["source"].eq("Satb2")].copy()
    require(len(satb2) > 0, "Frozen mouse prior must contain Satb2 for prespecified sensitivity")
    removed_self = satb2.loc[satb2["target"].eq("Satb2")].copy()
    require(len(removed_self) == 1, "Expected exactly one frozen Satb2 -> Satb2 link for prespecified sensitivity")
    no_self_prior = satb2.loc[~satb2["target"].eq("Satb2")].copy()
    try:
        for endpoint in ENDPOINTS:
            folder = out / endpoint
            folder.mkdir()
            mapping, counts, factors, md, audit = endpoint_inputs(cfg, endpoint, roster)
            expression, gene_audit = prepare_expression(mapping, counts, factors)
            wald, wald_audit = prepare_wald(mapping, audit)
            table_write(gene_audit, folder / "all_gene_activity_background_audit.tsv")
            table_write(wald_audit, folder / "all_results_wald_background_audit.tsv")
            md["saved_original_DESeq2_size_factor"] = factors.loc[md["sample"]].to_numpy()
            table_write(md, folder / "mouse_normalization_coverage_audit.tsv")
            table_write(expression.rename_axis("sample").reset_index(), folder / "mouse_activity_input_log1p_original_sizefactor.tsv.gz")
            table_write(wald.rename_axis("contrast").reset_index(), folder / "wald_contrast_input_original_finite_statistics.tsv.gz")
            symbol_ok = ~gene_audit["symbol_missing_or_invalid"] & ~gene_audit["symbol_ambiguous_multiple_stable_ids"]
            assayed = mapping.loc[symbol_ok, "gene_symbol"]
            measured = mapping.loc[symbol_ok & counts.sum(axis=1).gt(0), "gene_symbol"]
            for mode, matrix in [("mouse_log1p_original_sizefactor_expression", expression), ("original_DESeq2_Wald_supplement", wald)]:
                cov, used = coverage_and_edges(net, matrix.columns, endpoint, mode, assayed,
                                               measured if mode.startswith("mouse_") else matrix.columns)
                require(len(used) > 0, f"No tmin5 eligible TFs for {endpoint}/{mode}")
                es, fdr = run_native_ulm(matrix, used, dc)
                stem = "mouse_activity" if mode.startswith("mouse_") else "wald_contrast_enrichment"
                table_write(es.rename_axis("sample").reset_index(), folder / f"{stem}_native_escore.tsv")
                table_write(fdr.rename_axis("sample").reset_index(), folder / f"{stem}_native_gene_set_fdr.tsv")
                table_write(cov, folder / f"{stem}_all_TF_target_coverage.tsv")
                edges = used.copy()
                edges["edge_basis"] = "curated_mouse_CollecTRI_signed_prior_not_motif_inference"
                edges["endpoint"] = endpoint; edges["mode"] = mode
                table_write(edges, folder / f"{stem}_used_signed_prior_edges.tsv")
                long = long_scores(es, fdr, endpoint, mode, cov, md if stem == "mouse_activity" else None)
                table_write(long, folder / f"{stem}_long.tsv")
                all_coverage.append(cov)
                if stem == "mouse_activity":
                    require(len(long) == 12*es.shape[1] and not long.duplicated(["endpoint", "sample", "TF"]).any(), "Mouse activity completeness failed")
                    all_activity.append(long)
                    if cfg.get("plots", True):
                        plots(es, md, endpoint, folder, cfg)
                else:
                    all_contrast.append(long)
                summaries.append({"endpoint": endpoint, "mode": mode, "n_input_rows": len(matrix),
                                  "n_background_genes": matrix.shape[1], "n_modeled_TFs": es.shape[1],
                                  "n_score_rows": len(long), "all_native_estimates_finite": True})
            self_folder = folder / "self_excluded"
            self_folder.mkdir()
            statuses = {}
            for mode, matrix in [("mouse_log1p_original_sizefactor_expression_no_Satb2_self_target", expression),
                                 ("original_DESeq2_Wald_supplement_no_Satb2_self_target", wald)]:
                cov, used = coverage_and_edges(no_self_prior, matrix.columns, endpoint, mode, assayed,
                                               measured if mode.startswith("mouse_") else matrix.columns)
                cov["network_variant"] = "SATB2_only_self_target_excluded"
                cov["original_Satb2_prior_links"] = len(satb2)
                cov["excluded_self_target_links"] = len(removed_self)
                self_coverage.append(cov)
                stem = "mouse_activity" if mode.startswith("mouse_") else "wald_contrast_enrichment"
                table_write(cov, self_folder / f"{stem}_target_coverage.tsv")
                table_write(used, self_folder / f"{stem}_used_signed_prior_edges.tsv")
                if not len(used):
                    statuses[mode] = {"status": "UNAVAILABLE", "reason": "FEWER_THAN_5_ELIGIBLE_NONSELF_TARGETS",
                                      "n_eligible_targets": int(cov.iloc[0]["n_eligible_background_targets"])}
                    continue
                es, fdr = run_native_ulm(matrix, used, dc)
                require(es.columns.tolist() == ["Satb2"], "Self-excluded supplement must contain only Satb2")
                table_write(es.rename_axis("sample").reset_index(), self_folder / f"{stem}_native_escore.tsv")
                table_write(fdr.rename_axis("sample").reset_index(), self_folder / f"{stem}_native_gene_set_fdr.tsv")
                long = long_scores(es, fdr, endpoint, mode, cov, md if stem == "mouse_activity" else None)
                long["network_variant"] = "SATB2_only_self_target_excluded"
                long["gene_set_fdr_family"] = "native_BH_within_single_Satb2_set_NOT_primary_across_TFs_NOT_genotype_association"
                table_write(long, self_folder / f"{stem}_long.tsv")
                (self_activity if stem == "mouse_activity" else self_contrast).append(long)
                statuses[mode] = {"status": "COMPLETE", "n_input_rows": len(matrix),
                                  "n_background_genes": matrix.shape[1], "n_used_targets": len(used)}
            self_status[endpoint] = statuses
            json_write({"endpoint": endpoint, "variant": "SATB2_only_self_target_excluded",
                        "main_prior_unchanged": True, "removed_prior_link": removed_self.to_dict(orient="records"),
                        "background_genes_unchanged": True, "modes": statuses}, self_folder / "status.json")
        table_write(pd.concat(all_activity, ignore_index=True), out / "TF_ACTIVITY_PER_MOUSE.tsv")
        table_write(pd.concat(all_contrast, ignore_index=True), out / "TF_WALD_CONTRAST_ENRICHMENT.tsv")
        table_write(pd.concat(all_coverage, ignore_index=True), out / "TF_TARGET_COVERAGE.tsv")
        sensitivity_cols = list(all_activity[0].columns) + ["network_variant"]
        table_write(pd.concat(self_activity, ignore_index=True) if self_activity else pd.DataFrame(columns=sensitivity_cols),
                    out / "SATB2_NO_SELF_TARGET_ACTIVITY_PER_MOUSE.tsv")
        sensitivity_cols = list(all_contrast[0].columns) + ["network_variant"]
        table_write(pd.concat(self_contrast, ignore_index=True) if self_contrast else pd.DataFrame(columns=sensitivity_cols),
                    out / "SATB2_NO_SELF_TARGET_WALD_CONTRAST_ENRICHMENT.tsv")
        table_write(pd.concat(self_coverage, ignore_index=True), out / "SATB2_NO_SELF_TARGET_COVERAGE.tsv")
        after = verify_config_inputs(cfg)
        require(before == after and receipt(args.config) == cfg_receipt and receipt(__file__) == source_receipt, "Source/config/input changed during analysis")
        versions = {name: importlib.metadata.version(name) for name in ["decoupler", "numpy", "pandas", "scipy", "numba", "anndata", "matplotlib"]}
        implementation = {name: receipt(inspect.getsourcefile(module)) for name, module in [("ULM", dc.mt.ulm.func), ("Method", dc.mt.ulm.__class__)]}
        json_write({"schema": "satb2-tf-activity-provenance/v1", "analysis_id": cfg["analysis_id"],
                    "config": cfg_receipt, "source": source_receipt, "inputs_before": before, "inputs_after": after,
                    "software": versions, "python": sys.version, "executable": sys.executable, "platform": platform.platform(),
                    "native_api_implementation": implementation, "endpoint_summaries": summaries,
                    "Satb2_self_target_sensitivity": {"main_full_prior_unchanged": True,
                        "single_removed_link": removed_self.to_dict(orient="records"), "background_genes_unchanged": True,
                        "statuses": self_status, "native_gene_set_fdr_family": "within_single_Satb2_set_not_primary_across_TFs",
                        "caveat": "Full-prior Satb2 activity may include its own RNA via curated Satb2->Satb2; separate no-self-target sensitivity is not the primary fit."},
                    "normalization": "log1p(symbol_raw_counts / saved_original_matching_DESeq2_size_factor)",
                    "gene_filter": "raw_symbol_count>=5_in>=3_of_all12_genotype_neutral_mice; allzero_removed",
                    "ambiguous_symbol_policy": "exclude_all_symbols_mapping_to_multiple_stable_IDs; no_alias_guessing",
                    "background": "all_eligible_genes_not_just_network_targets_not_DEG_filtered",
                    "contrast_background": "all_finite_unambiguous_original_DESeq2_Wald_statistics_without_expression_transform_or_FDR_filter",
                    "gene_set_fdr": "native_decoupler_BH_across_TFs_within_each_input_mouse_or_contrast_NOT_genotype_association",
                    "genotype_association_run": False, "new_gene_DESeq2_run": False, "edgeR_run": False,
                    "TF_RNA_expression_substituted_for_activity": False, "motif_regulon_run": False,
                    "nucleus_activity_UMAP_run": False}, out / "PROVENANCE.json")
        outputs = {str(p.relative_to(out)): receipt(p) for p in sorted(out.rglob("*")) if p.is_file()}
        json_write({"schema": "satb2-tf-activity-complete/v1", "status": "COMPLETE", "analysis_id": cfg["analysis_id"],
                    "finished_UTC": datetime.now(timezone.utc).isoformat(), "source": source_receipt, "config": cfg_receipt,
                    "checks": {"all_five_endpoints": True, "twelve_mice_per_TF_endpoint": True,
                               "six_mice_per_genotype": True, "native_score_and_gene_set_fdr_finite": True,
                               "no_genotype_specific_or_DEG_filter": True, "original_size_factors_reused": True,
                               "full_raw_counts_library_sums_match": True, "ALL_RESULTS_stable_ID_audit_matches": True,
                               "input_source_config_SHA256_unchanged": True, "no_new_gene_DE_or_genotype_score_test": True},
                    "output_files": outputs}, out / "COMPLETE.json")
        print(json.dumps({"status": "COMPLETE", "outdir": str(out), "endpoint_summaries": summaries}))
    except Exception as error:
        json_write({"status": "FAILED", "error": str(error), "traceback": traceback.format_exc(),
                    "complete": False}, out / "FAILED.json")
        raise


if __name__ == "__main__":
    main()
