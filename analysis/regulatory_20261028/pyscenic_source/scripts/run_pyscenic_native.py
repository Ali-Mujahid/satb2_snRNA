#!/usr/bin/env python3
"""Additive, checkpointed native pySCENIC inference for the frozen SATB2 cohort.

No scientific action runs on import. Execute only from the registered native
workflow. This script never downloads references, patches a package, overwrites
a checkpoint, or treats nuclei as independent genotype replicates.
"""
import argparse
import gc
import gzip
import hashlib
import importlib.metadata
import json
import logging
import math
import os
import pickle
import random
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import io, sparse


ENDPOINTS = {
    "CA1_Deep": ("CA1_DL", "P0244_seurat_CA1_ProS_LAYER_PROGRAM", 8750, "CA1"),
    "CA1_Superficial": ("CA1_SL", "P0247_seurat_CA1_ProS_LAYER_PROGRAM", 11706, "CA1"),
    "IN_C1": ("IN_C1", "P0140_seurat_Inhibitory_CLUSTER", 1718, "IN"),
    "IN_C6": ("IN_C6", "P0145_seurat_Inhibitory_CLUSTER", 501, "IN"),
    "IN_C8": ("IN_C8", "P0147_seurat_Inhibitory_CLUSTER", 468, "IN"),
}
ROSTER = {"wt_" + s: "WT" for s in ("1122", "1123", "1178", "1220", "1226", "1228")}
ROSTER.update({"cKO_" + s: "cKO" for s in ("1138", "1139", "1160", "1165", "933", "991")})
FIXED = {
    "seed": 1337, "num_workers": 4, "min_detection_fraction": 0.05,
    "normalization_cp10k": 10000, "module_min_genes": 20,
    "rho_threshold": 0.03, "rank_threshold": 1500,
    "ctx_auc_threshold": 0.05, "nes_threshold": 3.0,
    "motif_similarity_fdr": 0.001, "orthologous_identity_threshold": 0.0,
    "min_regulon_targets": 10, "aucell_threshold": 0.05,
    "aucell_batch_size": 1000, "aucell_num_workers": 1,
}
BASE_COLS = ["lineage", "endpoint", "endpoint_id", "source_endpoint_id",
             "regulon_id", "native_name", "TF", "n_targets", "n_background_genes"]
MEAN_COLS = BASE_COLS + ["sample", "condition", "n_nuclei", "aucell_mean", "score_unit", "eligible", "network_variant"]
COV_COLS = BASE_COLS + ["n_native_targets", "n_self_targets", "eligible", "reason", "network_variant"]
TARGET_COLS = ["lineage", "regulon_id", "native_name", "TF", "target", "weight", "self_target", "edge_origin", "edge_basis"]


def need(value, message):
    if not value:
        raise ValueError(message)


def utc():
    return datetime.now(timezone.utc).isoformat()


def receipt(path):
    raw = Path(path)
    need(not raw.is_symlink(), "Symlink file rejected: " + str(raw))
    p = raw.resolve()
    need(p.is_file(), "Expected regular file: " + str(p))
    h = hashlib.sha256()
    with p.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            h.update(block)
    return {"path": str(p), "bytes": p.stat().st_size, "sha256": h.hexdigest()}


def verify(r):
    actual = receipt(r["path"])
    need(actual["bytes"] == r["bytes"] and actual["sha256"] == r["sha256"],
         "Input bytes or SHA256 differ: " + r["path"])
    return actual


def json_write(path, data):
    p = Path(path)
    need(not p.exists(), "Refusing overwrite: " + str(p))
    p.write_text(json.dumps(data, indent=2, sort_keys=True, allow_nan=False) + "\n")


def table_write(df, path, columns=None):
    p = Path(path)
    need(not p.exists(), "Refusing overwrite: " + str(p))
    if columns is not None:
        df = df.reindex(columns=columns)
    options = {"method": "gzip", "mtime": 0} if p.suffix == ".gz" else None
    df.to_csv(p, sep="\t", index=False, na_rep="NA", float_format="%.17g", compression=options)


def pickle_write(path, value):
    need(not Path(path).exists(), "Refusing overwrite: " + str(path))
    with gzip.open(path, "wb", compresslevel=6) as stream:
        pickle.dump(value, stream, protocol=4)


def pickle_read(path):
    # Callers verify their own native stage checkpoint before reading this file.
    with gzip.open(path, "rb") as stream:
        return pickle.load(stream)


def fresh_dir(path):
    p = Path(path)
    need(not p.is_symlink(), "Symlink output directory rejected: " + str(p))
    if p.exists():
        need(p.is_dir() and not any(p.iterdir()), "Output stage is not strictly empty: " + str(p))
    else:
        p.mkdir(parents=True)
    return p


def settings(cfg):
    need(cfg.get("schema", "satb2-pyscenic/v1") == "satb2-pyscenic/v1", "Unexpected scientific configuration schema")
    need(cfg.get("allow_network", False) is False and cfg.get("overwrite_prior", False) is False,
         "Protocol forbids network acquisition and prior overwrites inside scientific tasks")
    for key, value in FIXED.items():
        need(cfg.get(key, value) == value, "This frozen protocol fixes " + key + "=" + str(value))
    need(set(cfg["lineages"]) == {"CA1", "IN"}, "Exactly CA1 and IN lineages required")
    need(len(cfg["references"]["rankings"]) == 2, "Exactly two matching mm10 ranking databases required")
    need("analysis_id" in cfg and cfg["analysis_id"], "Missing analysis_id")
    if "expected_group_counts" in cfg:
        need(cfg["expected_group_counts"] == {k: v[2] for k, v in ENDPOINTS.items()}, "Frozen configured endpoint counts differ")
    return dict(FIXED)


def runtime_versions():
    versions = {p: importlib.metadata.version(p) for p in
                ("pyscenic", "arboreto", "ctxcore", "numpy", "pandas", "scipy", "scikit-learn", "dask", "distributed", "pyarrow", "numba")}
    need(versions["pyscenic"] == "0.12.1", "Wrong pySCENIC version")
    need(versions["arboreto"] == "0.1.6" and versions["ctxcore"] == "0.2.0", "Wrong native inference versions")
    return {"python": sys.version, "executable": sys.executable, "packages": versions}


def configured_inputs(cfg):
    out = {"prepared_complete": cfg["prepared_complete"], "samples": cfg["samples"]}
    for lineage, records in cfg["lineages"].items():
        for kind in ("counts", "genes", "cells"):
            out[lineage + "/" + kind] = records[kind]
    for kind in ("tf_names", "motif_annotations"):
        out["reference/" + kind] = cfg["references"][kind]
    for i, r in enumerate(cfg["references"]["rankings"]):
        out["reference/rankings/" + str(i)] = r
    for key, r in cfg.get("additional_hash_inputs", {}).items():
        out["additional/" + key] = r
    return out


def verify_inputs(cfg):
    checked = {k: verify(r) for k, r in configured_inputs(cfg).items()}
    complete = json.loads(Path(cfg["prepared_complete"]["path"]).read_text())
    need(complete["status"] == "COMPLETE" and complete["source_files_unchanged"] is True,
         "Native preparation not complete/unchanged")
    need(complete["source_all_gene_pseudobulk_matches"] is True, "Original pseudobulk reconciliation missing")
    need(complete["group_counts"] == {k: v[2] for k, v in ENDPOINTS.items()}, "Frozen endpoint counts differ")
    roster = pd.read_csv(cfg["samples"]["path"], sep="\t", keep_default_na=False)
    need(len(roster) == 12 and roster["sample"].is_unique and dict(zip(roster["sample"], roster["condition"])) == ROSTER,
         "Hash-bound independent mouse roster differs")
    for lin, info in cfg["lineages"].items():
        for kind, filename in (("counts", "counts.mtx.gz"), ("genes", "genes.tsv"), ("cells", "cells.tsv")):
            native = complete["files"]["lineages/" + lin + "/" + filename]
            need(native["sha256"] == info[kind]["sha256"] and native["bytes"] == info[kind]["bytes"],
                 "Input disagrees with completed native preparation: " + lin + "/" + kind)
    return checked


def checkpoint(out, cfg_path, cfg, stage, lineage, summary):
    source = receipt(__file__)
    outputs = {str(p.relative_to(out)): receipt(p) for p in sorted(out.rglob("*")) if p.is_file()}
    json_write(out / "PENDING_COMPLETE.json", {"schema": "satb2-pyscenic-stage/v1", "status": "COMPLETE",
               "analysis_id": cfg["analysis_id"], "stage": stage, "lineage": lineage,
               "finished_UTC": utc(), "config": receipt(cfg_path), "source": source,
               "parameters": settings(cfg), "summary": summary, "output_files": outputs})


def verified_stage(root, lineage, stage, cfg_path):
    out = root / lineage / stage
    marker = out / "COMPLETE.json"
    d = json.loads(marker.read_text())
    need(d["schema"] == "satb2-pyscenic-stage/v1" and d["status"] == "COMPLETE", "Incomplete stage: " + str(out))
    need(d["lineage"] == lineage and d["stage"] == stage, "Checkpoint identity differs")
    need(d["config"]["sha256"] == receipt(cfg_path)["sha256"], "Configuration changed between stages")
    verify(d["source"])
    need(d["source"]["sha256"] == receipt(__file__)["sha256"], "Scientific source changed between stages")
    for r in d["output_files"].values():
        verify(r)
    return out, d


def bools(series):
    need(series.astype(str).isin(["True", "False", "TRUE", "FALSE"]).all(), "Invalid boolean metadata")
    return series.astype(str).str.lower().eq("true")


def prepare_matrix(counts, genes, cells, fraction=0.05):
    """Pure genotype-blind normalization/filter helper, also used by synthetic tests."""
    x = sparse.csc_matrix(counts, dtype=np.float64)
    need(x.shape == (len(genes), len(cells)), "MTX orientation/metadata mismatch")
    need(np.isfinite(x.data).all() and (x.data >= 0).all() and np.equal(x.data, np.floor(x.data)).all(), "Invalid raw integer RNA")
    library = np.asarray(x.sum(axis=0)).ravel()
    need((library > 0).all(), "Zero original full-gene library")
    need(np.array_equal(library, cells["library_umi_full_genes"].to_numpy(dtype=float)), "Full-library denominator mismatch")
    symbols = genes["gene_symbol"].astype(str)
    unique = ~symbols.duplicated(keep=False)
    valid = unique & ~bools(genes["duplicate_symbol"]) & ~bools(genes["symbol_missing"])
    valid &= symbols.ne("") & ~symbols.isin(["NA", "nan", "None"]) & symbols.eq(symbols.str.strip())
    audit = genes.copy()
    audit["valid_unambiguous_symbol"] = valid
    detected = []
    for ep in sorted(cells["endpoint"].unique()):
        cols = np.flatnonzero(cells["endpoint"].to_numpy() == ep)
        f = np.asarray((x[:, cols] > 0).sum(axis=1)).ravel() / len(cols)
        audit["detected_fraction__" + ep] = f
        detected.append(f)
    detection = np.max(np.vstack(detected), axis=0) >= fraction
    norm = x.multiply((10000.0 / library)[None, :]).tocsc()
    norm.data = np.log1p(norm.data)
    # Sparse moment check excludes constant profiles without densifying all RNA.
    means = np.asarray(norm.mean(axis=1)).ravel()
    variance = np.asarray(norm.power(2).mean(axis=1)).ravel() - means ** 2
    nonconstant = variance > np.finfo(float).eps * np.maximum(1, means ** 2)
    eligible = np.asarray(valid) & detection & nonconstant
    audit["genotype_blind_detection_eligible"] = detection
    audit["normalized_nonconstant"] = nonconstant
    audit["eligible_full_RNA_background"] = eligible
    audit["reason"] = np.select([~np.asarray(valid), ~detection, ~nonconstant],
                                ["MISSING_OR_AMBIGUOUS_SYMBOL", "DETECTED_IN_FEWER_THAN_5_PERCENT_IN_EVERY_ENDPOINT", "CONSTANT_NORMALIZED_RNA"], default="ELIGIBLE")
    return norm[eligible, :].T.tocsc(), symbols[eligible].tolist(), audit


def stage_prepare(root, lineage, cfg, cfg_path):
    out = fresh_dir(root / lineage / "prepared")
    records = cfg["lineages"][lineage]
    genes = pd.read_csv(records["genes"]["path"], sep="\t", keep_default_na=False)
    cells = pd.read_csv(records["cells"]["path"], sep="\t", keep_default_na=False)
    expected = {k: v[2] for k, v in ENDPOINTS.items() if v[3] == lineage}
    need(cells["endpoint"].value_counts().to_dict() == expected, "Exact endpoint membership counts differ")
    need(cells["nucleus_id"].is_unique and cells["barcode"].equals(cells["nucleus_id"]), "Nonunique/unaligned nucleus IDs")
    need(set(cells["sample"]) == set(ROSTER), "Independent-mouse roster changed")
    need(all(ROSTER[s] == c for s, c in zip(cells["sample"], cells["condition"])), "Genotype differs from frozen roster")
    for ep in expected:
        need(set(cells.loc[cells.endpoint == ep, "sample"]) == set(ROSTER), "Endpoint lacks a frozen mouse")
    with gzip.open(records["counts"]["path"], "rb") as stream:
        raw = io.mmread(stream)
    expression, names, audit = prepare_matrix(raw, genes, cells)
    del raw
    need(len(names) >= 20, "Insufficient eligible full RNA background")
    tf_all = [s.strip() for s in Path(cfg["references"]["tf_names"]["path"]).read_text().splitlines() if s.strip()]
    need(len(tf_all) == len(set(tf_all)) and len(tf_all) > 1000, "Invalid official mouse TF list")
    matched = sorted(set(names).intersection(tf_all))
    need(len(matched) >= 2, "Too few measured official TFs")
    sparse.save_npz(out / "expression_cells_by_genes.npz", expression, compressed=True)
    table_write(pd.DataFrame({"gene_symbol": names}), out / "expression_genes.tsv")
    table_write(cells, out / "nucleus_metadata.tsv.gz")
    table_write(audit, out / "gene_background_audit.tsv.gz")
    audit_symbols = {s: i for i, s in enumerate(genes.gene_symbol.astype(str))}
    candidate = []
    for tf in tf_all:
        i = audit_symbols.get(tf)
        candidate.append({"lineage": lineage, "TF": tf, "official_TF_candidate": True,
                          "assayed_symbol": i is not None, "eligible_RNA": tf in set(names),
                          "reason": audit.iloc[i]["reason"] if i is not None else "NOT_ASSAYED",
                          "forced_into_TF_list": False})
    table_write(pd.DataFrame(candidate), out / "TF_CANDIDATE_COVERAGE.tsv")
    (out / "eligible_official_TFs.txt").write_text("\n".join(matched) + "\n")
    detected_n = np.asarray((expression > 0).sum(axis=1)).ravel()
    low_n = int(np.ceil(len(names) * FIXED["aucell_threshold"]))
    auc_audit = pd.DataFrame({"nucleus_id": cells.nucleus_id, "detected_eligible_genes": detected_n,
                              "AUCell_top_fraction": FIXED["aucell_threshold"], "AUCell_rank_window_approx_genes": low_n,
                              "rank_window_includes_undetected_RNA": detected_n < low_n})
    table_write(auc_audit, out / "AUCell_detection_window_audit.tsv.gz")
    checkpoint(out, cfg_path, cfg, "prepared", lineage,
               {"n_nuclei": len(cells), "n_genes": len(names), "n_TF_candidates": len(matched),
                "Satb2_in_official_list": "Satb2" in tf_all, "Satb2_eligible_RNA": "Satb2" in matched,
                "AUCell_nuclei_window_includes_undetected_RNA": int((detected_n < low_n).sum()),
                "normalization": "log1p(10000*raw_count/full_original_all_gene_library)",
                "genotype_blind": True, "downsampled_nuclei": 0, "DEG_or_HVG_selection": False})


def load_expression(out):
    matrix = sparse.load_npz(out / "expression_cells_by_genes.npz").tocsc()
    genes = pd.read_csv(out / "expression_genes.tsv", sep="\t").gene_symbol.tolist()
    cells = pd.read_csv(out / "nucleus_metadata.tsv.gz", sep="\t", keep_default_na=False)
    need(matrix.shape == (len(cells), len(genes)), "Prepared expression axes changed")
    return matrix, genes, cells


def local_client(n):
    from distributed import Client, LocalCluster
    cluster = LocalCluster(n_workers=n, threads_per_worker=1, processes=False,
                           dashboard_address=None, memory_limit=0)
    return cluster, Client(cluster)


def native_grn(expression, genes, tf_names, seed, workers):
    from arboreto.algo import grnboost2
    cluster, client = local_client(workers)
    try:
        adj = grnboost2(expression_data=expression, gene_names=genes, tf_names=tf_names,
                       client_or_address=client, seed=seed, verbose=True)
    finally:
        client.close()
        cluster.close()
    need(list(adj.columns) == ["TF", "target", "importance"] and not adj.empty, "Unexpected native GRNBoost2 output")
    need(np.isfinite(adj.importance).all() and (adj.importance >= 0).all(), "Invalid native importance")
    need(not adj[["TF", "target"]].duplicated().any(), "Duplicate native adjacency")
    need(set(adj.TF).issubset(tf_names) and set(adj.target).issubset(genes), "Native GRN changed gene/TF universe")
    need(set(adj.target) == set(genes), "Native GRN returned no links for one or more nonconstant RNA targets; inspect warnings")
    return adj.sort_values(["importance", "TF", "target"], ascending=[False, True, True], kind="mergesort").reset_index(drop=True)


def stage_grn(root, lineage, cfg, cfg_path):
    prep, _ = verified_stage(root, lineage, "prepared", cfg_path)
    out = fresh_dir(root / lineage / "grn")
    x, genes, cells = load_expression(prep)
    tfs = (prep / "eligible_official_TFs.txt").read_text().splitlines()
    adj = native_grn(x, genes, tfs, FIXED["seed"], FIXED["num_workers"])
    table_write(adj, out / "adjacencies_native_GRNBoost2.tsv.gz")
    checkpoint(out, cfg_path, cfg, "grn", lineage,
               {"n_nuclei": len(cells), "n_genes": len(genes), "n_TFs": len(tfs), "n_edges": len(adj),
                "algorithm": "native_arboreto_grnboost2", "seed": FIXED["seed"],
                "execution": "LocalCluster_processesFalse_four_single_thread_workers", "positive_importance_is_not_regulatory_sign": True})


def motif_table_plain(df):
    if df.empty:
        return pd.DataFrame(columns=["TF", "MotifID", "NES", "AUC", "Context", "TargetGenes", "Annotation"])
    flat = df.copy()
    if flat.columns.nlevels == 2:
        flat.columns = flat.columns.droplevel(0)
    flat = flat.reset_index()
    for col in ("Context", "TargetGenes"):
        if col in flat:
            flat[col] = flat[col].map(lambda v: json.dumps(sorted(v), ensure_ascii=False))
    return flat


def stage_ctx(root, lineage, cfg, cfg_path):
    from ctxcore.rnkdb import FeatherRankingDatabase
    from pyscenic.utils import add_correlation, modules_from_adjacencies
    from pyscenic.prune import prune2df
    from pyscenic.transform import df2regulons
    prep, _ = verified_stage(root, lineage, "prepared", cfg_path)
    grn, _ = verified_stage(root, lineage, "grn", cfg_path)
    out = fresh_dir(root / lineage / "ctx")
    x, genes, cells = load_expression(prep)
    df = pd.DataFrame(x.toarray(), index=cells.nucleus_id, columns=genes)
    del x
    adj = pd.read_csv(grn / "adjacencies_native_GRNBoost2.tsv.gz", sep="\t")
    correlated = add_correlation(adj, df, rho_threshold=FIXED["rho_threshold"], mask_dropouts=False)
    need(np.isfinite(correlated.rho).all(), "Nonfinite TF/target correlation")
    modules = list(modules_from_adjacencies(
        correlated, df, thresholds=(0.75, 0.90), top_n_targets=(50,), top_n_regulators=(5, 10, 50),
        min_genes=20, absolute_thresholds=False, rho_dichotomize=True,
        keep_only_activating=True, rho_threshold=0.03, rho_mask_dropouts=False))
    need(all("activating" in m.context for m in modules), "Nonpositive candidate module")
    del df, adj
    gc.collect()
    databases = [FeatherRankingDatabase(r["path"], r.get("name", Path(r["path"]).stem)) for r in cfg["references"]["rankings"]]
    table_write(correlated, out / "adjacencies_native_with_Pearson_rho.tsv.gz")
    del correlated
    pickle_write(out / "candidate_modules_native.pkl.gz", modules)
    candidate_rows = []
    module_rows = []
    for idx, m in enumerate(modules):
        module_id = lineage + "::candidate_" + str(idx)
        module_rows.append({"module_id": module_id, "native_name": m.name, "TF": m.transcription_factor,
                            "n_genes": len(m), "context": json.dumps(sorted(m.context))})
        for db in databases:
            n = len(set(m.genes).intersection(db.geneset))
            candidate_rows.append({"module_id": module_id, "native_name": m.name, "TF": m.transcription_factor,
                                   "database": db.name, "n_module_genes": len(m), "n_mapped_genes": n,
                                   "fraction_mapped": n / len(m), "native_mapping_gate_pass": n / len(m) > 0.80})
    table_write(pd.DataFrame(module_rows, columns=["module_id", "native_name", "TF", "n_genes", "context"]), out / "CANDIDATE_MODULES.tsv")
    table_write(pd.DataFrame(candidate_rows, columns=["module_id", "native_name", "TF", "database", "n_module_genes", "n_mapped_genes", "fraction_mapped", "native_mapping_gate_pass"]), out / "MODULE_RANKING_DATABASE_COVERAGE.tsv")
    if modules:
        cluster, client = local_client(4)
        try:
            motifs = prune2df(databases, modules, cfg["references"]["motif_annotations"]["path"],
                              rank_threshold=1500, auc_threshold=0.05, nes_threshold=3.0,
                              motif_similarity_fdr=0.001, orthologuous_identity_threshold=0.0,
                              weighted_recovery=False, client_or_address=client, num_workers=4,
                              module_chunksize=1, filter_for_annotation=True)
        finally:
            client.close()
            cluster.close()
    else:
        motifs = pd.DataFrame()
    pickle_write(out / "motif_enrichment_native_dataframe.pkl.gz", motifs)
    # Keep native multi-index CSV alongside an explicit JSON-cell TSV adaptation.
    need(not (out / "motif_enrichment_native.csv.gz").exists(), "Refusing native motif overwrite")
    motifs.to_csv(out / "motif_enrichment_native.csv.gz", compression={"method": "gzip", "mtime": 0}, float_format="%.17g")
    table_write(motif_table_plain(motifs), out / "motif_enrichment_flat.tsv.gz")
    regulons = sorted(df2regulons(motifs, save_columns=["NES", "OrthologousIdentity", "MotifSimilarityQvalue", "Annotation"]), key=lambda r: r.name) if not motifs.empty else []
    need(all(r.name.endswith("(+)") and "activating" in r.context for r in regulons), "Nonpositive final regulon")
    need(len({r.name for r in regulons}) == len(regulons), "Duplicate native regulon names")
    need(all(set(r.genes).issubset(genes) for r in regulons), "Regulon genes outside measured RNA background")
    pickle_write(out / "regulons_native.pkl.gz", regulons)
    target_rows, regulon_rows = [], []
    for r in regulons:
        need(math.isfinite(float(r.score)) and all(math.isfinite(float(w)) and float(w) > 0 for w in r.gene2weight.values()),
             "Nonfinite/nonpositive native regulon score or target weight")
        rid = lineage + "::" + r.name
        regulon_rows.append({"lineage": lineage, "regulon_id": rid, "native_name": r.name, "TF": r.transcription_factor,
                             "n_native_targets": len(r), "n_targets": len(r), "n_background_genes": len(genes),
                             "n_self_targets": int(r.transcription_factor in r.genes), "eligible": len(r) >= 10,
                             "reason": "ELIGIBLE" if len(r) >= 10 else "FEWER_THAN_10_MEASURED_MOTIF_PRUNED_TARGETS",
                             "native_regulon_score": float(r.score), "native_context": json.dumps(sorted(r.context))})
        for gene, weight in sorted(r.gene2weight.items()):
            self_gene = gene == r.transcription_factor
            target_rows.append({"lineage": lineage, "regulon_id": rid, "native_name": r.name, "TF": r.transcription_factor,
                                "target": gene, "weight": float(weight), "self_target": self_gene,
                                "edge_origin": "native_module_added_TF_self_gene" if self_gene else "native_GRNBoost2_then_cisTarget_leading_edge_union",
                                "edge_basis": "positive_coexpression_plus_motif_resource_support_not_measured_binding_or_genotype_effect"})
    table_write(pd.DataFrame(regulon_rows, columns=["lineage", "regulon_id", "native_name", "TF", "n_native_targets", "n_targets", "n_background_genes", "n_self_targets", "eligible", "reason", "native_regulon_score", "native_context"]), out / "NATIVE_REGULON_COVERAGE.tsv")
    table_write(pd.DataFrame(target_rows, columns=TARGET_COLS), out / "REGULON_TARGETS.tsv")
    checkpoint(out, cfg_path, cfg, "ctx", lineage,
               {"n_nuclei": len(cells), "n_genes": len(genes), "n_candidate_modules": len(modules),
                "n_motif_enrichment_rows": len(motifs), "n_regulons": len(regulons),
                "n_score_eligible_regulons": sum(len(r) >= 10 for r in regulons),
                "status": "COMPLETE_WITH_REGULONS" if regulons else "EMPTY_NO_SUPPORTED_REGULONS",
                "direct_binding_proven": False, "positive_correlation_means_proven_activation": False})


def native_auc_batched(matrix, genes, cell_ids, regulons, batch_size=1000, seed=1337):
    """Native row-independent AUCell with identical seeded gene permutation per batch."""
    from pyscenic.aucell import aucell
    need(len(regulons) > 0, "No eligible native regulons to score")
    parts = []
    for start in range(0, matrix.shape[0], batch_size):
        stop = min(start + batch_size, matrix.shape[0])
        exp = pd.DataFrame(matrix[start:stop, :].toarray(), index=cell_ids[start:stop], columns=genes)
        auc = aucell(exp, regulons, auc_threshold=0.05, noweights=False, normalize=False, seed=seed, num_workers=1)
        need(auc.index.is_unique and set(auc.index) == set(cell_ids[start:stop]), "Native AUCell changed nucleus membership")
        need(set(auc.columns) == {r.name for r in regulons}, "Native AUCell changed regulon axes")
        # Native serial AUCell unstack sorts cell labels. Restore only their
        # verified labels/order; do not substitute coordinates or score values.
        auc = auc.reindex(index=list(cell_ids[start:stop]), columns=[r.name for r in regulons])
        need(np.isfinite(auc.to_numpy()).all() and (auc.to_numpy() >= 0).all() and (auc.to_numpy() <= 1 + 1e-12).all(),
             "Invalid native nucleus AUC")
        parts.append(auc)
        print("AUCell nuclei", stop, "/", matrix.shape[0], flush=True)
    joined = pd.concat(parts, axis=0)
    need(joined.index.is_unique and joined.index.tolist() == list(cell_ids), "Batched nucleus axes differ")
    return joined


def mouse_means(auc, cells, lineage, regulons, n_background, variant, suffix=""):
    rows = []
    for r in regulons:
        rid = lineage + "::" + r.name + suffix
        for ep in sorted(cells.endpoint.unique()):
            for sample in ROSTER:
                idx = np.flatnonzero((cells.endpoint == ep).to_numpy() & (cells["sample"] == sample).to_numpy())
                need(len(idx) > 0, "Missing independent-mouse endpoint group")
                v = auc.iloc[idx][r.name].to_numpy()
                need(np.isfinite(v).all(), "Nonfinite actual-nucleus AUC")
                info = ENDPOINTS[ep]
                rows.append({"lineage": lineage, "endpoint": ep, "endpoint_id": info[0], "source_endpoint_id": info[1],
                             "regulon_id": rid, "native_name": r.name, "TF": r.transcription_factor, "n_targets": len(r),
                             "n_background_genes": n_background, "sample": sample, "condition": ROSTER[sample],
                             "n_nuclei": len(idx), "aucell_mean": float(v.mean()), "score_unit": "native_AUCell_AUC",
                             "eligible": True, "network_variant": variant})
    return pd.DataFrame(rows, columns=MEAN_COLS)


def stage_aucell(root, lineage, cfg, cfg_path):
    prep, _ = verified_stage(root, lineage, "prepared", cfg_path)
    ctx, info = verified_stage(root, lineage, "ctx", cfg_path)
    out = fresh_dir(root / lineage / "aucell")
    x, genes, cells = load_expression(prep)
    regs = pickle_read(ctx / "regulons_native.pkl.gz")
    eligible = [r for r in regs if len(r) >= 10]
    # This explicit axis file is authoritative even for an empty-regulon lineage.
    axis = cells[["nucleus_id", "sample", "condition", "endpoint", "lineage", "sequencing_partition"]].copy()
    table_write(axis, out / "nucleus_axis_metadata.tsv.gz")
    if eligible:
        auc = native_auc_batched(x, genes, cells.nucleus_id.tolist(), eligible)
        scoring = pd.DataFrame({"regulon_id": [lineage + "::" + r.name for r in eligible], "native_name": [r.name for r in eligible],
                                "n_finite_nuclei": len(cells), "n_zero_AUC": (auc.to_numpy() == 0).sum(axis=0),
                                "min_AUC": auc.min().to_numpy(), "max_AUC": auc.max().to_numpy(),
                                "nonconstant_AUC": (auc.max() > auc.min()).to_numpy(),
                                "all_zero_AUC": (auc.max() == 0).to_numpy()})
        table_write(scoring, out / "AUCell_score_integrity.tsv")
        means = mouse_means(auc, cells, lineage, eligible, len(genes), "native_motif_pruned_positive_regulon")
    else:
        auc = pd.DataFrame(index=cells.nucleus_id)
        means = pd.DataFrame(columns=MEAN_COLS)
        table_write(pd.DataFrame(columns=["regulon_id", "native_name", "n_finite_nuclei", "n_zero_AUC", "min_AUC", "max_AUC", "nonconstant_AUC", "all_zero_AUC"]), out / "AUCell_score_integrity.tsv")
    export_auc = auc.rename(columns={r.name: lineage + "::" + r.name for r in eligible})
    table_write(export_auc.reset_index().rename(columns={"index": "nucleus_id", "Cell": "nucleus_id"}), out / "nucleus_AUC.tsv.gz")
    table_write(means, out / "REGULON_ACTIVITY_PER_MOUSE.tsv", MEAN_COLS)
    main_cov = pd.read_csv(ctx / "NATIVE_REGULON_COVERAGE.tsv", sep="\t", keep_default_na=False)
    cov_rows = []
    for r in main_cov.to_dict("records"):
        for ep in sorted(cells.endpoint.unique()):
            cov_rows.append(dict(r, endpoint=ep, endpoint_id=ENDPOINTS[ep][0], source_endpoint_id=ENDPOINTS[ep][1],
                                 network_variant="native_motif_pruned_positive_regulon"))
    table_write(pd.DataFrame(cov_rows, columns=COV_COLS), out / "REGULON_COVERAGE.tsv", COV_COLS)
    sat = next((r for r in regs if r.transcription_factor == "Satb2" and r.name == "Satb2(+)"), None)
    no_self = None
    reason = "NO_NATIVE_MOTIF_SUPPORTED_SATB2_REGULON"
    if sat is not None:
        if "Satb2" not in sat.genes:
            reason = "SELF_GENE_ABSENT"
        elif len(sat) - 1 < 10:
            reason = "FEWER_THAN_10_REMAINING_TARGETS_AFTER_SELF_EXCLUSION"
        else:
            no_self = sat.copy(gene2weight={g: w for g, w in sat.gene2weight.items() if g != "Satb2"})
            reason = "ELIGIBLE_SEPARATE_SATB2_SELF_EXCLUDED_SENSITIVITY"
    no_cov = []
    for ep in sorted(cells.endpoint.unique()):
        no_cov.append({"lineage": lineage, "endpoint": ep, "endpoint_id": ENDPOINTS[ep][0], "source_endpoint_id": ENDPOINTS[ep][1],
                       "regulon_id": lineage + "::Satb2(+)_NO_SELF_TARGET", "native_name": "Satb2(+)", "TF": "Satb2",
                       "n_targets": len(no_self) if no_self else max(0, len(sat) - int("Satb2" in sat.genes)) if sat else 0,
                       "n_native_targets": len(sat) if sat else 0, "n_background_genes": len(genes), "n_self_targets": 0,
                       "eligible": no_self is not None, "reason": reason, "network_variant": "SATB2_native_regulon_self_gene_excluded"})
    table_write(pd.DataFrame(no_cov, columns=COV_COLS), out / "SATB2_NO_SELF_TARGET_COVERAGE.tsv", COV_COLS)
    if no_self:
        sens_auc = native_auc_batched(x, genes, cells.nucleus_id.tolist(), [no_self])
        export_sens_auc = sens_auc.rename(columns={no_self.name: lineage + "::" + no_self.name + "_NO_SELF_TARGET"})
        table_write(export_sens_auc.reset_index().rename(columns={"index": "nucleus_id", "Cell": "nucleus_id"}), out / "SATB2_NO_SELF_TARGET_nucleus_AUC.tsv.gz")
        sens_mean = mouse_means(sens_auc, cells, lineage, [no_self], len(genes), "SATB2_native_regulon_self_gene_excluded", "_NO_SELF_TARGET")
    else:
        sens_mean = pd.DataFrame(columns=MEAN_COLS)
    table_write(sens_mean, out / "SATB2_NO_SELF_TARGET_ACTIVITY_PER_MOUSE.tsv", MEAN_COLS)
    checkpoint(out, cfg_path, cfg, "aucell", lineage,
               {"n_nuclei": len(cells), "n_genes": len(genes), "n_regulons": len(regs), "n_score_eligible_regulons": len(eligible),
                "n_mouse_endpoint_regulon_rows": len(means), "native_AUCell": True, "nucleus_axes_exact": True,
                "batch_size": 1000, "seeded_tie_breaking": True, "normalize_to_observed_max": False,
                "Satb2_self_excluded_status": reason,
                "status": "COMPLETE_WITH_REGULONS" if regs else "EMPTY_NO_SUPPORTED_REGULONS"})


def stage_assemble(root, cfg, cfg_path):
    need(not (root / "COMPLETE.json").exists(), "Native assembled outputs already complete")
    frames = {k: [] for k in ("REGULON_ACTIVITY_PER_MOUSE", "REGULON_COVERAGE", "REGULON_TARGETS",
                              "SATB2_NO_SELF_TARGET_ACTIVITY_PER_MOUSE", "SATB2_NO_SELF_TARGET_COVERAGE")}
    lineages = []
    paths = {}
    for lin in ("CA1", "IN"):
        stages = {s: verified_stage(root, lin, s, cfg_path) for s in ("prepared", "grn", "ctx", "aucell")}
        for key in frames:
            stage = "ctx" if key == "REGULON_TARGETS" else "aucell"
            frames[key].append(pd.read_csv(stages[stage][0] / (key + ".tsv"), sep="\t", keep_default_na=False))
        summary = stages["aucell"][1]["summary"]
        lineages.append({"lineage": lin, "status": summary["status"], "n_regulons": summary["n_regulons"],
                         "n_score_eligible_regulons": summary["n_score_eligible_regulons"], "n_nuclei": summary["n_nuclei"], "n_genes": summary["n_genes"]})
        paths[lin] = {"nucleus_AUC": receipt(stages["aucell"][0] / "nucleus_AUC.tsv.gz"),
                      "nucleus_axis_metadata": receipt(stages["aucell"][0] / "nucleus_axis_metadata.tsv.gz"),
                      "prepared_metadata": receipt(stages["prepared"][0] / "nucleus_metadata.tsv.gz"),
                      "native_regulons": receipt(stages["ctx"][0] / "regulons_native.pkl.gz"),
                      "native_motif_enrichment": receipt(stages["ctx"][0] / "motif_enrichment_native_dataframe.pkl.gz"),
                      "status": summary["status"]}
    for key, parts in frames.items():
        table_write(pd.concat(parts, ignore_index=True), root / (key + ".tsv"))
    table_write(pd.DataFrame(lineages), root / "LINEAGE_STATUS.tsv")
    means = pd.concat(frames["REGULON_ACTIVITY_PER_MOUSE"], ignore_index=True)
    if len(means):
        need(not means[["lineage", "endpoint", "regulon_id", "sample"]].duplicated().any(), "Duplicate mouse activity")
        for _, part in means.groupby(["endpoint", "regulon_id"]):
            need(len(part) == 12 and set(part["sample"]) == set(ROSTER), "Incomplete all12 mouse means")
    verified = verify_inputs(cfg)
    checks = {k: True for k in ["positive_motif_regulons_only", "full_genotype_neutral_RNA_background", "no_DEG_selection",
                               "native_nucleus_AUCell", "exact_nucleus_axes_verified", "all12_mouse_roster", "min_measured_targets_10",
                               "Satb2_self_target_sensitivity_separate", "inputs_config_source_SHA256_unchanged"]}
    provenance = {"schema": "satb2-pyscenic-native-provenance/v1", "analysis_id": cfg["analysis_id"], "finished_UTC": utc(),
                  "config": receipt(cfg_path), "source": receipt(__file__), "parameters": settings(cfg), "runtime": runtime_versions(),
                  "inputs": verified, "lineages": paths, "checks": checks,
                  "method": "native_arboreto_GRNBoost2_positive_modules_native_pySCENIC_cisTarget_pruning_AUCell",
                  "cohort": "pooled_genotypes_per_lineage_for_discovery_then_all_exact_nucleus_scoring_and_independent_mouse_means",
                  "no_new_gene_DE": True, "nuclei_are_not_genotype_replicates": True,
                  "regulon_gene_weight": "native_df2regulons_union_maximum_importance_not_binding_strength_or_genotype_logFC",
                  "self_gene": "Native modules add TF gene; retained and explicitly flagged; not inferred direct self_binding",
                  "caveats": ["Positive coexpression plus motif resource support does not prove activation, binding, or causality",
                              "CA1 labels are molecular layer programs, not physical anatomical assignments",
                              "IN_C8 has mixed identity and one mouse contributes only13 nuclei",
                              "Sequencing partition is partly confounded with genotype; 933 partition is unverified",
                              "Pooled data infer regulons; mouse comparisons on the same cohort are exploratory and not independent validation",
                              "AUCell background and regulon gene sets vary between lineages; absolute cross-lineage scores are not calibrated"]}
    json_write(root / "PROVENANCE.json", provenance)
    output_files = {str(p.relative_to(root)): receipt(p) for p in sorted(root.rglob("*")) if p.is_file()}
    json_write(root / "PENDING_COMPLETE.json", {"schema": "satb2-pyscenic-native-complete/v1", "status": "COMPLETE",
               "analysis_id": cfg["analysis_id"], "finished_UTC": utc(), "config": receipt(cfg_path),
               "source": receipt(__file__), "checks": checks, "lineages": lineages, "output_files": output_files})


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--stage", choices=["prepare", "grn", "ctx", "aucell", "assemble"], required=True)
    parser.add_argument("--lineage", choices=["CA1", "IN"])
    args = parser.parse_args(argv)
    cfg_path = Path(args.config).resolve()
    cfg = json.loads(cfg_path.read_text())
    settings(cfg)
    runtime_versions()
    need(args.stage == "assemble" or args.lineage is not None, "Lineage required for this stage")
    need(not Path(args.outdir).is_symlink(), "Symlink native output root rejected")
    root = Path(args.outdir).resolve()
    need(root == Path(cfg["run_root"]).resolve() / "results" / "pyscenic_native",
         "Native outdir must be the new configured run_root/results/pyscenic_native")
    root.mkdir(parents=True, exist_ok=True)
    before_config, before_source = receipt(cfg_path), receipt(__file__)
    before_inputs = verify_inputs(cfg)
    random.seed(1337)
    np.random.seed(1337)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    if args.stage == "assemble":
        stage_assemble(root, cfg, cfg_path)
    else:
        globals()["stage_" + args.stage](root, args.lineage, cfg, cfg_path)
    need(receipt(cfg_path) == before_config and receipt(__file__) == before_source, "Source/config changed during stage")
    need(verify_inputs(cfg) == before_inputs, "Configured inputs changed during stage")
    stage_dir = root if args.stage == "assemble" else root / args.lineage / ("prepared" if args.stage == "prepare" else args.stage)
    need(not (stage_dir / "COMPLETE.json").exists(), "Refusing completion marker overwrite")
    (stage_dir / "PENDING_COMPLETE.json").rename(stage_dir / "COMPLETE.json")


if __name__ == "__main__":
    main()
