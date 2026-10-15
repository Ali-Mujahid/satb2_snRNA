#!/usr/bin/env python3
"""Fresh matched Scanpy representation/clustering; no biological condition DE."""
from __future__ import annotations
import argparse
import gc
import hashlib
import json
import os
from pathlib import Path
import shutil
import time

for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
            "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(key, "4")

import anndata as ad
import harmonypy
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd
import psutil
import scanpy as sc
from scipy import sparse

from common import (population_mask, read_table, sha256, sparse_digest, versions,
                    write_h5ad, write_json, write_table)


def resource_guard(cfg):
    minimum = float(cfg.get("minimum_available_gib", 24))
    available = psutil.virtual_memory().available / 1024 ** 3
    if available < minimum:
        raise RuntimeError(f"Available RAM {available:.2f} GiB below {minimum:g} GiB guard")


def representation_settings(cfg):
    """Effective scientific settings, independent of run-specific input paths."""
    defaults = {
        "target_sum": 10000, "scale_max": 10, "n_pcs": 30,
        "n_neighbors": 20, "umap_metric": "cosine", "umap_min_dist": .3,
        "resolutions": [.2, .4, .6], "main_resolution": .4,
        "harmony_key": "sample", "harmony_theta": 2, "harmony_sigma": .1,
        "harmony_lambda": 1.0, "harmony_nclust": 100, "harmony_max_iter": 10,
        "harmony_max_iter_kmeans": 4, "harmony_epsilon": .01,
        "seed": 20261013, "threads": 4, "minimum_marker_cluster_nuclei": 20,
        "marker_top_n": 100,
    }
    settings = {key: cfg.get(key, value) for key, value in defaults.items()}
    settings["umap_n_neighbors"] = cfg.get("umap_n_neighbors", cfg.get("umap_neighbors", 30))
    settings["marker_nuclei_per_mouse_cluster"] = cfg.get(
        "marker_nuclei_per_mouse_cluster", cfg.get("marker_max_nuclei_per_sample_cluster", 100))
    return settings


def saved_count_layer_digest(path):
    """Stream the saved integer CSR count layer; never load or fit the matrix."""
    with h5py.File(path, "r") as handle:
        group = handle["layers/counts"]
        if group.attrs.get("encoding-type") != "csr_matrix":
            raise ValueError("Cached representation count layer must be CSR")
        digest = hashlib.sha256()
        digest.update(str(tuple(map(int, group.attrs["shape"]))).encode())
        for key in ("data", "indices", "indptr"):
            values = group[key]
            if key == "data" and not np.issubdtype(values.dtype, np.integer):
                raise ValueError("Cached representation counts must remain integer")
            digest.update(str(values.dtype).encode())
            for start in range(0, len(values), 1_000_000):
                part = np.ascontiguousarray(values[start:start + 1_000_000])
                digest.update(memoryview(part).cast("B"))
        return digest.hexdigest()


def reuse_completed_population(cfg, config_path, shared_dir, population, out):
    """Copy an explicitly selected completed population with verified fit inputs."""
    old = Path(cfg["reuse_completed_populations"][population]).resolve()
    branch = cfg.get("count_branch", "original_full")
    if old == out:
        raise ValueError("Population cache source and destination must differ")
    receipt_path = old / "COMPLETE.json"
    previous = json.loads(receipt_path.read_text())
    if previous.get("status") != "COMPLETE" or previous.get("count_branch") != branch or previous.get("population") != population:
        raise ValueError("Population cache lacks a complete matching branch/population receipt")
    old_results = old.parents[2]
    old_config_path = old_results / "inputs" / f"scanpy_{branch}.json"
    if sha256(old_config_path) != previous["config_sha256"]:
        raise ValueError("Population cache's original configuration identity changed")
    old_config = json.loads(old_config_path.read_text())
    current_settings = representation_settings(cfg)
    if current_settings != representation_settings(old_config):
        raise ValueError("Population cache scientific representation/marker parameters changed")
    for key in ("markers_tsv",):
        current_path, old_path = cfg.get(key), old_config.get(key)
        if bool(current_path) != bool(old_path) or (current_path and sha256(current_path) != sha256(old_path)):
            raise ValueError("Population cache marker-panel input changed")
    shared_receipt = json.loads((shared_dir / "COMPLETE.json").read_text())
    old_shared_receipt = json.loads((old_results / "shared" / branch / "COMPLETE.json").read_text())
    for receipt in (shared_receipt, old_shared_receipt):
        if receipt.get("status") != "COMPLETE" or receipt.get("count_branch") != branch:
            raise ValueError("Population cache requires completed matching shared count inputs")
    for key in ("branch_counts_h5ad_sha256", "branch_counts_digest", "mapmycells_prediction_sha256"):
        if shared_receipt[key] != old_shared_receipt[key]:
            raise ValueError(f"Population cache count/annotation input changed: {key}")
    if sha256(shared_dir / "branch_counts.h5ad") != shared_receipt["branch_counts_h5ad_sha256"]:
        raise ValueError("Current shared count artifact checksum changed")
    hvg_hash = sha256(shared_dir / f"common_hvg_{population}.tsv")
    if hvg_hash != previous["common_hvg_sha256"]:
        raise ValueError("Population cache does not use the current exact common HVGs")
    old_files = [path for path in old.iterdir() if path.name != "COMPLETE.json"]
    if not old_files or any(not path.is_file() or path.is_symlink() for path in old_files):
        raise ValueError("Population cache must contain regular completed artifact files only")
    required = {"scanpy.h5ad", "nucleus_metadata.tsv", "cluster_counts_per_mouse.tsv", "representation_review.pdf"}
    if not required.issubset({path.name for path in old_files}):
        raise ValueError("Population cache is missing a completed representation artifact")
    saved_path = old / "scanpy.h5ad"
    if sha256(saved_path) != previous["scanpy_h5ad_sha256"]:
        raise ValueError("Cached representation H5AD checksum changed")
    current = read_table(shared_dir / f"membership_{population}.tsv")
    stored = read_table(old / "nucleus_metadata.tsv")
    if current.nucleus_id.duplicated().any() or stored.nucleus_id.duplicated().any():
        raise ValueError("Population cache/current membership has duplicate nucleus IDs")
    if len(current) != int(previous["n_nuclei"]) or set(current.nucleus_id) != set(stored.nucleus_id):
        raise ValueError("Population cache nucleus roster differs from current shared membership")
    stored = stored.set_index("nucleus_id").loc[current.nucleus_id]
    for column in ("sample", "condition"):
        if not np.array_equal(current[column].astype(str).to_numpy(), stored[column].astype(str).to_numpy()):
            raise ValueError(f"Population cache conflicts with current {column} identity")
    backed = ad.read_h5ad(saved_path, backed="r")
    try:
        saved_obs = backed.obs
        if backed.shape != (len(current), int(previous["n_genes"])) or set(saved_obs.nucleus_id.astype(str)) != set(current.nucleus_id):
            raise ValueError("Population cache H5AD gene/nucleus axes differ from its receipt")
        saved_obs = saved_obs.set_index("nucleus_id").loc[current.nucleus_id]
        for column in ("sample", "condition"):
            if not np.array_equal(current[column].astype(str).to_numpy(), saved_obs[column].astype(str).to_numpy()):
                raise ValueError("Population cache H5AD mouse/genotype metadata conflicts")
        for basis in ("unintegrated", "harmony"):
            key = "cluster_" + basis
            if saved_obs[key].astype(str).tolist() != stored[key].astype(str).tolist() or saved_obs[key].nunique() != previous["clusters_" + basis]:
                raise ValueError("Population cache H5AD and TSV clustering metadata disagree")
    finally:
        backed.file.close()
    if not previous.get("count_preservation_verified") or saved_count_layer_digest(saved_path) != previous["branch_counts_digest"]:
        raise ValueError("Population cache saved integer count layer digest changed")
    copied_hashes = {}
    for path in old_files:
        expected_hash = sha256(path)
        destination = out / path.name
        shutil.copyfile(path, destination)
        if sha256(destination) != expected_hash:
            raise ValueError(f"Population artifact copy checksum mismatch: {path.name}")
        copied_hashes[path.name] = expected_hash
    fresh = dict(previous)
    fresh.update(config_sha256=sha256(config_path), reused_completed_population_from=str(old),
                 reused_completed_receipt_sha256=sha256(receipt_path),
                 original_config_sha256=previous["config_sha256"], copied_artifact_sha256=copied_hashes,
                 representation_refit_performed=False, validated_representation_settings=current_settings,
                 current_shared_counts_h5ad_sha256=shared_receipt["branch_counts_h5ad_sha256"])
    write_json(out / "COMPLETE.json", fresh)
    print(f"Reused verified completed Scanpy representation: {branch}/{population}, {len(current):,} nuclei", flush=True)


def scatter_panel(ax, coordinates, labels, title):
    labels = pd.Series(labels).astype(object).fillna("NA").astype(str)
    categories = list(pd.unique(labels))
    colors = plt.colormaps["tab20"].resampled(max(20, len(categories)))
    for i, category in enumerate(categories):
        mask = labels.eq(category).to_numpy()
        ax.scatter(coordinates[mask, 0], coordinates[mask, 1], s=1.2, linewidths=0,
                   c=[colors(i)], alpha=.7, rasterized=True, label=category)
    ax.set_title(title, fontsize=10)
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_xticks([])
    ax.set_yticks([])
    if len(categories) <= 40:
        ax.legend(loc="center left", bbox_to_anchor=(1, .5), frameon=False,
                  markerscale=4, fontsize=6)


def figures(data, out, population, cfg):
    with PdfPages(out / "representation_review.pdf") as pdf:
        for field in ["sample", "condition", "sequencing_partition", "mmc_reference_broad_prediction",
                      "cluster_unintegrated", "cluster_harmony"]:
            if field not in data.obs:
                continue
            fig, axes = plt.subplots(1, 2, figsize=(13, 5), constrained_layout=True)
            scatter_panel(axes[0], data.obsm["X_umap_unintegrated"], data.obs[field], f"{population}: unintegrated / {field}")
            scatter_panel(axes[1], data.obsm["X_umap_harmony"], data.obs[field], f"{population}: Harmony / {field}")
            pdf.savefig(fig)
            plt.close(fig)
        ratio = data.uns["pca"]["variance_ratio"]
        fig, ax = plt.subplots(figsize=(8, 4), constrained_layout=True)
        ax.plot(np.arange(1, len(ratio)+1), ratio, marker=".")
        ax.set(xlabel="Principal component", ylabel="Explained variance ratio",
               title=f"{population}: fixed common HVGs; log-normalized scaled PCA")
        pdf.savefig(fig)
        plt.close(fig)
        for field in ["cluster_unintegrated", "cluster_harmony"]:
            proportions = pd.crosstab(data.obs["sample"], data.obs[field], normalize="index")
            fig, ax = plt.subplots(figsize=(10, 4.5), constrained_layout=True)
            proportions.plot.bar(stacked=True, ax=ax, width=.85, legend=False)
            ax.set(ylabel="Fraction of this population in each mouse", title=f"{population}: {field}")
            ax.tick_params(axis="x", rotation=45)
            if proportions.shape[1] <= 40:
                ax.legend(loc="center left", bbox_to_anchor=(1, .5), frameon=False, fontsize=6)
            pdf.savefig(fig)
            plt.close(fig)


def descriptive_markers(data, out, population, cfg):
    # Annotation markers compare nuclei between clusters; genotype tests belong to
    # the separate mouse-level pseudobulk stage, never to this cell-level test.
    for basis in ["unintegrated", "harmony"]:
        groupby = "cluster_" + basis
        # Deterministic equal contribution from each represented mouse within a
        # cluster, capped at 100 nuclei/mouse; biological DE uses all nuclei.
        rng = np.random.default_rng(int(cfg.get("seed", 20261013)))
        selected = []
        for cluster in data.obs[groupby].cat.categories:
            cluster_md = data.obs.loc[data.obs[groupby].eq(cluster)]
            per_mouse = cluster_md.groupby("sample", observed=True).size()
            quota = min(int(cfg.get("marker_nuclei_per_mouse_cluster",
                                    cfg.get("marker_max_nuclei_per_sample_cluster", 100))), int(per_mouse.min()))
            for sample in per_mouse.index:
                ids = cluster_md.loc[cluster_md["sample"].eq(sample)].index.to_numpy()
                selected.extend(rng.choice(ids, size=quota, replace=False).tolist())
        marker_md = data.obs.loc[selected, ["nucleus_id", "sample", "condition", groupby]].copy()
        marker_md["marker_selection_seed"] = int(cfg.get("seed", 20261013))
        marker_cap = int(cfg.get("marker_nuclei_per_mouse_cluster", cfg.get("marker_max_nuclei_per_sample_cluster", 100)))
        marker_md["marker_selection_rule"] = f"equal_per_represented_mouse_within_cluster_capped_at_{marker_cap}"
        write_table(out / f"marker_selected_roster_{basis}.tsv", marker_md.reset_index(drop=True))
        marker_data = ad.AnnData(X=data.X[data.obs_names.get_indexer(selected)].copy(),
                                obs=marker_md, var=data.var.copy())
        marker_data.obs[groupby] = pd.Categorical(marker_data.obs[groupby])
        sizes = marker_data.obs[groupby].value_counts()
        eligible_groups = sizes.loc[sizes >= int(cfg.get("minimum_marker_cluster_nuclei", 20))].index.tolist()
        if len(eligible_groups) < 1 or marker_data.obs[groupby].nunique() < 2:
            write_json(out / f"markers_{basis}_status.json", {"status": "SKIPPED", "reason": "Insufficient cluster groups/nuclei"})
            continue
        key = "markers_" + basis
        sc.tl.rank_genes_groups(marker_data, groupby=groupby, groups=eligible_groups,
                               reference="rest", method="wilcoxon", use_raw=False,
                               n_genes=int(cfg.get("marker_top_n", 100)), pts=True,
                               tie_correct=True, corr_method="benjamini-hochberg", key_added=key)
        frame = sc.get.rank_genes_groups_df(marker_data, group=None, key=key)
        frame = frame.rename(columns={"names": "gene_id", "group": "cluster"})
        frame.insert(0, "population", population)
        frame.insert(1, "representation", basis)
        frame["gene_symbol"] = frame.gene_id.map(data.var.gene_symbol)
        frame["inference_scope"] = "descriptive_cell_level_cluster_markers_not_genotype_DE"
        write_table(out / f"cluster_markers_{basis}.tsv", frame)
        del marker_data
        gc.collect()


def marker_panel_review(data, out, cfg):
    path = cfg.get("markers_tsv")
    if not path:
        return
    panel = read_table(path)
    if not {"group", "gene_symbol"}.issubset(panel):
        raise ValueError("Marker panel requires group and gene_symbol columns")
    if "eligible" in panel:
        from common import boolean
        panel = panel.loc[boolean(panel.eligible, "eligible")]
    lookup = data.var.reset_index().rename(columns={data.var.index.name or "index": "gene_id"})
    audit = panel.merge(lookup[["gene_id", "gene_symbol"]], on="gene_symbol", how="left", validate="many_to_many")
    write_table(out / "marker_panel_gene_id_audit.tsv", audit)
    rows = []
    for cluster in data.obs.cluster_harmony.cat.categories:
        mask = data.obs.cluster_harmony.eq(cluster).to_numpy()
        x = data.X[mask]
        for row in audit.dropna(subset=["gene_id"]).itertuples(index=False):
            position = data.var_names.get_loc(row.gene_id)
            values = x[:, position]
            rows.append({"cluster": str(cluster), "marker_group": row.group,
                         "gene_id": row.gene_id, "gene_symbol": row.gene_symbol,
                         "mean_log_normalized_expression": float(values.mean()),
                         "fraction_detected": float((values > 0).mean()), "n_nuclei": int(mask.sum())})
    if rows:
        write_table(out / "marker_panel_cluster_expression.tsv", pd.DataFrame(rows))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--config", required=True, type=Path)
    p.add_argument("--shared-dir", required=True, type=Path)
    p.add_argument("--population", required=True, choices=["allcells", "Excitatory", "Inhibitory", "CA1_ProS"])
    p.add_argument("--outdir", required=True, type=Path)
    args = p.parse_args()
    cfg = json.loads(args.config.read_text())
    resource_guard(cfg)
    out = args.outdir.resolve()
    if out.exists() and (not out.is_dir() or any(out.iterdir())):
        raise RuntimeError("Output destination contains data; use workflow caching or a new run directory")
    # Snakemake creates output parent directories before launching this stage.
    out.mkdir(parents=True, exist_ok=True)
    if args.population in cfg.get("reuse_completed_populations", {}):
        reuse_completed_population(cfg, args.config, args.shared_dir, args.population, out)
        return
    started = time.monotonic()
    seed = int(cfg.get("seed", 20261013))
    sc.settings.n_jobs = int(cfg.get("threads", 4))
    shared_path = args.shared_dir / "branch_counts.h5ad"
    shared_receipt = json.loads((args.shared_dir / "COMPLETE.json").read_text())
    if sha256(shared_path) != shared_receipt["branch_counts_h5ad_sha256"]:
        raise ValueError("Shared count artifact checksum changed")
    # Backed subsetting loads only this population into memory; allcells loads once.
    backed = ad.read_h5ad(shared_path, backed="r")
    mask = population_mask(backed, args.population)
    data = backed[mask].to_memory()
    backed.file.close()
    del backed
    data.obs["analysis_population"] = args.population
    for c in ["sample", "condition", "sequencing_partition", "mmc_reference_broad_prediction"]:
        if c in data.obs:
            data.obs[c] = pd.Categorical(data.obs[c])
    original_digest = sparse_digest(data.X)
    # Counts is the same sparse object until X is replaced with a float32 copy.
    data.layers["counts"] = data.X
    data.X = data.X.astype(np.float32)
    sc.pp.normalize_total(data, target_sum=float(cfg.get("target_sum", 10000)))
    sc.pp.log1p(data)
    hvg_path = args.shared_dir / f"common_hvg_{args.population}.tsv"
    hvg = read_table(hvg_path)
    positions = data.var_names.get_indexer(hvg.gene_id)
    if (positions < 0).any() or hvg.gene_id.duplicated().any():
        raise ValueError("Common HVGs absent or duplicated")
    data.var["highly_variable"] = data.var_names.isin(hvg.gene_id)
    work = ad.AnnData(X=data.X[:, positions].copy(), obs=data.obs[["sample"]].copy(),
                      var=data.var.iloc[positions].copy())
    sc.pp.scale(work, zero_center=True, max_value=float(cfg.get("scale_max", 10)))
    n_pcs = min(int(cfg.get("n_pcs", 30)), work.n_obs-1, work.n_vars-1)
    sc.pp.pca(work, n_comps=n_pcs, svd_solver="arpack", random_state=seed)
    data.obsm["X_pca"] = np.asarray(work.obsm["X_pca"], dtype=np.float32)
    data.uns["pca"] = work.uns["pca"]
    loadings = np.zeros((data.n_vars, n_pcs), dtype=np.float32)
    loadings[positions] = work.varm["PCs"]
    data.varm["PCs"] = loadings
    del work, loadings
    gc.collect()
    # Capture-level integration is a representation sensitivity only; retain the
    # unintegrated primary control and use counts for mouse-level inference.
    harmony_key = cfg.get("harmony_key", "sample")
    if harmony_key not in data.obs or data.obs[harmony_key].isna().any():
        raise ValueError("Harmony batch metadata missing")
    if harmony_key == "condition":
        raise ValueError("Never integrate away the tested WT/cKO condition")
    if data.obs[harmony_key].nunique() < 2:
        raise ValueError("Harmony requires at least two technical batches")
    result = harmonypy.run_harmony(data.obsm["X_pca"], data.obs, harmony_key,
                                   max_iter_harmony=int(cfg.get("harmony_max_iter", 10)),
                                   theta=float(cfg.get("harmony_theta", 2)),
                                   sigma=float(cfg.get("harmony_sigma", .1)),
                                   lamb=cfg.get("harmony_lambda", 1.0),
                                   nclust=int(cfg.get("harmony_nclust", 100)),
                                   max_iter_kmeans=int(cfg.get("harmony_max_iter_kmeans", 4)),
                                   epsilon_harmony=float(cfg.get("harmony_epsilon", .01)),
                                   ncores=int(cfg.get("threads", 4)),
                                   random_state=seed, verbose=False)
    corrected = np.asarray(result.Z_corr)
    # harmonypy 2.1.0 preserves the input (cells x PCs) orientation.
    if corrected.shape != (data.n_obs, n_pcs) or not np.isfinite(corrected).all():
        raise ValueError("Harmony returned malformed corrected PCs")
    data.obsm["X_pca_harmony"] = corrected.astype(np.float32)
    del result, corrected
    gc.collect()
    resolutions = cfg.get("resolutions", [.2, .4, .6])
    main_resolution = float(cfg.get("main_resolution", .4))
    if main_resolution not in [float(r) for r in resolutions]:
        raise ValueError("Main clustering resolution absent from sensitivity grid")
    for basis, representation in [("unintegrated", "X_pca"), ("harmony", "X_pca_harmony")]:
        neighbor_key = "neighbors_" + basis
        sc.pp.neighbors(data, n_neighbors=int(cfg.get("n_neighbors", 20)), use_rep=representation,
                        n_pcs=n_pcs, random_state=seed, key_added=neighbor_key)
        umap_neighbor_key = "umap_neighbors_" + basis
        sc.pp.neighbors(data, n_neighbors=int(cfg.get("umap_n_neighbors", cfg.get("umap_neighbors", 30))), use_rep=representation,
                        n_pcs=n_pcs, random_state=seed, metric=cfg.get("umap_metric", "cosine"),
                        key_added=umap_neighbor_key)
        sc.tl.umap(data, neighbors_key=umap_neighbor_key, random_state=seed,
                   min_dist=float(cfg.get("umap_min_dist", .3)))
        data.obsm["X_umap_" + basis] = data.obsm["X_umap"].copy()
        for resolution in resolutions:
            key = f"leiden_{basis}_r{float(resolution):g}"
            sc.tl.leiden(data, resolution=float(resolution), neighbors_key=neighbor_key,
                         key_added=key, random_state=seed, flavor="igraph", directed=False,
                         n_iterations=2)
        data.obs["cluster_" + basis] = data.obs[f"leiden_{basis}_r{main_resolution:g}"].copy()
    data.obsm["X_umap"] = data.obsm["X_umap_harmony"].copy()
    descriptive_markers(data, out, args.population, cfg)
    marker_panel_review(data, out, cfg)
    figures(data, out, args.population, cfg)
    if sparse_digest(data.layers["counts"]) != original_digest:
        raise RuntimeError("Saved original counts changed during representation analysis")
    metadata = data.obs.copy()
    for basis in ["unintegrated", "harmony"]:
        coordinates = data.obsm["X_umap_" + basis]
        metadata[f"UMAP_{basis}_1"] = coordinates[:, 0]
        metadata[f"UMAP_{basis}_2"] = coordinates[:, 1]
    write_table(out / "nucleus_metadata.tsv", metadata.reset_index(drop=True))
    cluster_counts = data.obs.groupby(["sample", "condition", "cluster_harmony"], observed=True).size().rename("n_nuclei").reset_index()
    write_table(out / "cluster_counts_per_mouse.tsv", cluster_counts)
    data.uns["analysis_provenance"] = {"population": args.population,
        "branch_counts_digest": original_digest, "common_hvg_sha256": sha256(hvg_path),
        "config_sha256": sha256(args.config), "harmony_key": harmony_key,
        "marker_inference_scope": "descriptive nuclei-vs-rest cluster markers; not WT/cKO inference"}
    write_h5ad(out / "scanpy.h5ad", data)
    write_json(out / "COMPLETE.json", {"status": "COMPLETE", "population": args.population,
        "n_nuclei": data.n_obs, "n_genes": data.n_vars, "n_pcs": n_pcs,
        "clusters_unintegrated": int(data.obs.cluster_unintegrated.nunique()),
        "clusters_harmony": int(data.obs.cluster_harmony.nunique()),
        "branch_counts_digest": original_digest, "count_preservation_verified": True,
        "count_branch": cfg.get("count_branch", "original_full"),
        "versions": versions(), "config_sha256": sha256(args.config),
        "common_hvg_sha256": sha256(hvg_path), "scanpy_h5ad_sha256": sha256(out / "scanpy.h5ad"),
        "seconds": time.monotonic()-started})


if __name__ == "__main__":
    main()
