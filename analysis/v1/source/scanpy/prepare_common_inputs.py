#!/usr/bin/env python3
"""Approved preparation stage: original counts, frozen roster, common HVGs."""
from __future__ import annotations
import argparse
import json
import os
import shutil
from pathlib import Path

for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
            "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(key, "4")

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import psutil
import scanpy as sc

from common import (boolean, join_mapmycells, load_original_counts, population_mask, sha256,
                    sparse_digest, versions, write_h5ad, write_json, write_table)


def resource_guard(cfg):
    available = psutil.virtual_memory().available
    minimum = float(cfg.get("minimum_available_gib", 24)) * 1024 ** 3
    if available < minimum:
        raise RuntimeError(f"Available RAM {available / 1024**3:.2f} GiB below {minimum / 1024**3:g} GiB guard")


def select_hvg(data, mask, cfg):
    counts = data.X if mask.all() else data.X[mask]
    detectable = np.asarray((counts > 0).sum(axis=0)).ravel() > 0
    symbols = data.var.gene_symbol.astype(str)
    # Exclusions apply to discovery features only, never to saved original counts.
    eligible = detectable.copy()
    for pattern in cfg.get("feature_discovery_exclude_patterns", []):
        eligible &= ~symbols.str.contains(pattern, regex=True).to_numpy()
    candidate = ad.AnnData(X=counts[:, eligible],
                          obs=data.obs.loc[mask, ["sample"]].copy(),
                          var=data.var.loc[eligible].copy())
    top = int(cfg.get("n_hvg", 3000))
    if candidate.n_vars < top:
        raise ValueError("Too few eligible genes for the agreed common HVG panel")
    sc.pp.highly_variable_genes(candidate, flavor=cfg.get("hvg_flavor", "seurat_v3_paper"),
                               n_top_genes=top, batch_key="sample", span=.3, check_values=True)
    selected = candidate.var.loc[candidate.var.highly_variable].copy()
    selected = selected.sort_values(["highly_variable_nbatches", "highly_variable_rank"],
                                    ascending=[False, True], kind="stable")
    if len(selected) != top:
        raise ValueError(f"HVG selection returned {len(selected)} rather than {top}")
    selected.insert(0, "gene_id", selected.index.astype(str))
    selected.insert(1, "common_hvg_order", np.arange(1, len(selected) + 1))
    return selected.reset_index(drop=True)


def reuse_prepared_shared(cfg, config_path, out):
    """Copy an explicitly verified completed count/HVG stage; do not refit it."""
    old = Path(cfg["reuse_prepared_shared_dir"]).resolve()
    branch = cfg.get("count_branch", "original_full")
    if branch not in ("original_full", "paired_cellbender"):
        raise ValueError("Shared cache reuse is authorized only for original_full or paired_cellbender")
    if old == out:
        raise ValueError("Shared cache source and destination must differ")
    receipt_path = old / "COMPLETE.json"
    previous = json.loads(receipt_path.read_text())
    if previous.get("status") != "COMPLETE" or previous.get("count_branch") != branch:
        raise ValueError("Shared cache lacks a complete matching count-branch receipt")
    expected_nuclei = int(cfg.get("expected_nuclei", 108062))
    if int(previous.get("n_nuclei", -1)) != expected_nuclei:
        raise ValueError("Shared cache cohort size differs from the current declared branch")
    populations = cfg.get("populations", ["allcells", "Excitatory", "Inhibitory", "CA1_ProS"])
    previous_populations = {row["population"]: row for row in previous["populations"]}
    if len(previous_populations) != len(populations) or set(previous_populations) != set(populations):
        raise ValueError("Shared cache populations differ from the current declared panels")
    old_config_path = old.parent.parent / "inputs" / f"scanpy_{branch}.json"
    if sha256(old_config_path) != previous["config_sha256"]:
        raise ValueError("Shared cache's original configuration identity changed")
    old_config = json.loads(old_config_path.read_text())
    for key, default in (("n_hvg", 3000), ("hvg_flavor", "seurat_v3_paper"),
                         ("feature_discovery_exclude_patterns", []),
                         ("keep_column", "passes_primary_QC")):
        if cfg.get(key, default) != old_config.get(key, default):
            raise ValueError(f"Shared cache scientific preparation setting changed: {key}")
    old_files = [path for path in old.iterdir() if path.name != "COMPLETE.json"]
    if not old_files or any(not path.is_file() or path.is_symlink() for path in old_files):
        raise ValueError("Shared cache must contain regular prepared artifact files only")
    required = {"branch_counts.h5ad", "nucleus_metadata.tsv", "gene_features.tsv"}
    required.update(f"{prefix}_{population}.tsv" for population in populations for prefix in ("common_hvg", "membership"))
    if not required.issubset({path.name for path in old_files}):
        raise ValueError("Shared cache is missing a declared prepared artifact")
    if sha256(old / "branch_counts.h5ad") != previous["branch_counts_h5ad_sha256"]:
        raise ValueError("Shared cached count artifact checksum changed")
    with h5py.File(old / "branch_counts.h5ad", "r") as handle:
        if tuple(handle["X"].attrs["shape"]) != (expected_nuclei, int(previous["n_genes"])):
            raise ValueError("Shared cached count shape differs from its completed receipt")
    current = pd.read_csv(cfg["qc_roster"], sep="\t", dtype={"nucleus_id": str, "sample": str,
                                                            "condition": str, "barcode": str}, low_memory=False)
    keep_column = cfg.get("keep_column", "passes_primary_QC")
    if keep_column in current:
        current = current.loc[boolean(current[keep_column], keep_column)].copy()
    stored = pd.read_csv(old / "nucleus_metadata.tsv", sep="\t", dtype={"nucleus_id": str,
                      "sample": str, "condition": str, "barcode": str}, low_memory=False)
    if current.nucleus_id.duplicated().any() or stored.nucleus_id.duplicated().any():
        raise ValueError("Shared-cache/current roster has duplicate nucleus IDs")
    if len(current) != expected_nuclei or set(current.nucleus_id) != set(stored.nucleus_id):
        raise ValueError("Shared cached roster does not match the current authoritative nuclei")
    stored = stored.set_index("nucleus_id").loc[current.nucleus_id]
    for column in ("sample", "condition", "barcode"):
        if not np.array_equal(current[column].astype(str).to_numpy(), stored[column].astype(str).to_numpy()):
            raise ValueError(f"Shared cached roster conflicts with current {column} identity")
    if current["sample"].nunique() != 12 or current.groupby("sample", observed=True).condition.nunique().max() != 1:
        raise ValueError("Shared reuse requires the same twelve independent mouse identities")
    sample_sheet = pd.read_csv(cfg["sample_sheet"], sep="\t", dtype={"sample": str, "condition": str})
    source_by_sample = {row["sample"]: row for row in previous["sources"]}
    if sample_sheet["sample"].duplicated().any() or set(sample_sheet["sample"]) != set(source_by_sample):
        raise ValueError("Shared cache capture inventory differs from current sample sheet")
    for row in sample_sheet.itertuples(index=False):
        source = source_by_sample[row.sample]
        if row.source_h5_sha256 != source["sha256"] or int(row.n_nuclei) != source["n_nuclei"]:
            raise ValueError(f"Shared cached count source changed for {row.sample}")
        if hasattr(row, "branch_UMIs") and int(row.branch_UMIs) != source["UMIs"]:
            raise ValueError(f"Shared cached library totals changed for {row.sample}")
        if branch == "original_full":
            evidence = source.get("mapmycells_query_verification", {})
            if evidence.get("status") != "VERIFIED_EXACT_ORIGINAL_INTEGER_COUNTS" or evidence.get("normalization") != "raw" or evidence.get("unequal_sparse_count_entries") != 0:
                raise ValueError(f"Shared cache lacks accepted exact-original MapMyCells-query evidence for {row.sample}")
            if evidence.get("shape") != [source["n_nuclei"], int(previous["n_genes"])]:
                raise ValueError("Cached MapMyCells query evidence has inconsistent gene/nucleus axes")
    if sha256(cfg["mapmycells_predictions"]) != previous["mapmycells_prediction_sha256"]:
        raise ValueError("Shared cached MapMyCells annotation input changed")
    flags = {"Excitatory": "REF_EX_HS", "Inhibitory": "REF_IN_HS", "CA1_ProS": "REF_CA1_PROS_HS"}
    for population in populations:
        panel_path = old / f"common_hvg_{population}.tsv"
        if sha256(panel_path) != previous_populations[population]["hvg_sha256"]:
            raise ValueError(f"Shared cached HVG panel changed: {population}")
        panel = pd.read_csv(panel_path, sep="\t", dtype={"gene_id": str})
        if len(panel) != int(cfg.get("n_hvg", 3000)) or panel.gene_id.duplicated().any():
            raise ValueError("Shared cached common-HVG panel is malformed")
        if cfg.get("hvg_source_dir") and sha256(Path(cfg["hvg_source_dir"]) / panel_path.name) != sha256(panel_path):
            raise ValueError("Paired shared cache does not use the current original-full fixed HVGs")
        if population != "allcells":
            flag = flags[population]
            if not np.array_equal(boolean(current[flag], flag), boolean(stored[flag], flag)):
                raise ValueError("Shared cached atlas population flags differ from the current roster")
            expected_ids = set(current.loc[boolean(current[flag], flag), "nucleus_id"])
        else:
            expected_ids = set(current.nucleus_id)
        members = pd.read_csv(old / f"membership_{population}.tsv", sep="\t", dtype={"nucleus_id": str})
        if members.nucleus_id.duplicated().any() or set(members.nucleus_id) != expected_ids or len(members) != previous_populations[population]["n_nuclei"]:
            raise ValueError(f"Shared cached population membership differs: {population}")
    copied_hashes = {}
    for path in old_files:
        expected_hash = sha256(path)
        destination = out / path.name
        shutil.copyfile(path, destination)
        if sha256(destination) != expected_hash:
            raise ValueError(f"Shared artifact copy checksum mismatch: {path.name}")
        copied_hashes[path.name] = expected_hash
    fresh = dict(previous)
    fresh.update(config_sha256=sha256(config_path), roster_sha256=sha256(cfg["qc_roster"]),
                 mapmycells_prediction_sha256=sha256(cfg["mapmycells_predictions"]),
                 common_hvg_source_dir=cfg.get("hvg_source_dir", "selected_from_this_branch"),
                 reused_prepared_artifacts_from=str(old), reused_prepared_receipt_sha256=sha256(receipt_path),
                 copied_artifact_sha256=copied_hashes, scientific_preparation_refit_performed=False,
                 raw_count_sources_reloaded=False)
    write_json(out / "COMPLETE.json", fresh)
    print(f"Reused verified shared preparation: {branch}, {expected_nuclei:,} nuclei", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    cfg = json.loads(args.config.read_text())
    resource_guard(cfg)
    out = args.outdir.resolve()
    if out.exists() and (not out.is_dir() or any(out.iterdir())):
        raise RuntimeError("Output destination contains data; use workflow caching or a new run directory")
    # Snakemake creates output parent directories before launching this stage.
    out.mkdir(parents=True, exist_ok=True)
    if cfg.get("reuse_prepared_shared_dir"):
        reuse_prepared_shared(cfg, args.config, out)
        return
    sc.settings.n_jobs = int(cfg.get("threads", 4))
    data, sources = load_original_counts(cfg)
    before = sparse_digest(data.X)
    join_mapmycells(data, cfg["mapmycells_predictions"], cfg.get("mapmycells_allow_superset", False))
    data.uns["MapMyCells_annotation_evidence"] = cfg.get("annotation_provenance", data.uns["MapMyCells_annotation_evidence"])
    write_table(out / "nucleus_metadata.tsv", data.obs.reset_index(drop=True))
    write_table(out / "gene_features.tsv", data.var.rename_axis("gene_id").reset_index())
    populations = cfg.get("populations", ["allcells", "Excitatory", "Inhibitory", "CA1_ProS"])
    summaries = []
    for population in populations:
        mask = population_mask(data, population)
        cohort = data.obs.loc[mask, ["nucleus_id", "sample", "condition"]].copy()
        if cohort["sample"].nunique() != 12:
            raise ValueError(f"{population} does not have nuclei from all twelve mice")
        write_table(out / f"membership_{population}.tsv", cohort)
        hvg_source_dir = cfg.get("hvg_source_dir")
        if hvg_source_dir:
            source_hvg = Path(hvg_source_dir) / f"common_hvg_{population}.tsv"
            hvg = pd.read_csv(source_hvg, sep="\t", dtype={"gene_id": str})
            if len(hvg) != int(cfg.get("n_hvg", 3000)) or hvg.gene_id.duplicated().any():
                raise ValueError("Inherited original-full HVG file has invalid size or duplicates")
            if not set(hvg.gene_id).issubset(data.var_names):
                raise ValueError("Paired count branch lacks original-full common HVG IDs")
        else:
            hvg = select_hvg(data, mask, cfg)
        if hvg_source_dir:
            shutil.copyfile(source_hvg, out / f"common_hvg_{population}.tsv")
        else:
            write_table(out / f"common_hvg_{population}.tsv", hvg)
        summaries.append({"population": population, "n_nuclei": int(mask.sum()),
                          "n_hvg": len(hvg), "hvg_sha256": sha256(out / f"common_hvg_{population}.tsv")})
        print(f"Common HVGs prepared: {population}, {int(mask.sum()):,} nuclei", flush=True)
    if sparse_digest(data.X) != before:
        raise RuntimeError("Original counts changed during HVG preparation")
    data.uns["branch_counts_digest"] = before
    data.uns["common_config_sha256"] = sha256(args.config)
    write_h5ad(out / "branch_counts.h5ad", data)
    receipt = {"status": "COMPLETE", "config_sha256": sha256(args.config),
               "versions": versions(), "n_nuclei": data.n_obs, "n_genes": data.n_vars,
               "branch_counts_digest": before, "sources": sources, "populations": summaries,
               "roster_sha256": sha256(cfg["qc_roster"]),
               "mapmycells_prediction_sha256": sha256(cfg["mapmycells_predictions"]),
               "branch_counts_h5ad_sha256": sha256(out / "branch_counts.h5ad"),
               "count_branch": cfg.get("count_branch", "original_full"),
               "common_hvg_source_dir": cfg.get("hvg_source_dir", "selected_from_this_branch"),
               "normalization_for_future_matched_branches": "counts/total_UMIs * 10000, natural log1p"}
    write_json(out / "COMPLETE.json", receipt)


if __name__ == "__main__":
    main()
