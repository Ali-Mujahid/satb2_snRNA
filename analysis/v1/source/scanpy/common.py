"""Shared IO and integrity helpers; import performs no file IO or analysis."""
from __future__ import annotations

import hashlib
import importlib.metadata
import json
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse


def read_table(path):
    return pd.read_csv(path, sep="\t", dtype={"nucleus_id": str, "sample": str,
                                           "barcode": str, "condition": str,
                                           "gene_id": str}, low_memory=False)


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def sparse_digest(matrix):
    x = matrix.tocsr(copy=False)
    h = hashlib.sha256()
    h.update(str(x.shape).encode())
    for part in (x.data, x.indices, x.indptr):
        h.update(str(part.dtype).encode())
        h.update(memoryview(np.ascontiguousarray(part)).cast("B"))
    return h.hexdigest()


def boolean(values, name):
    s = values.astype(str).str.strip().str.lower()
    if not s.isin(["true", "false", "1", "0"]).all():
        raise ValueError(f"{name} requires explicit TRUE/FALSE or 1/0")
    return s.isin(["true", "1"]).to_numpy()


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + ".partial")
    temp.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    temp.replace(path)


def write_table(path, frame):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + ".partial")
    frame.to_csv(temp, sep="\t", index=False, na_rep="NA")
    temp.replace(path)


def write_h5ad(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.stem + ".partial.h5ad")
    data.write_h5ad(temp, compression="gzip", compression_opts=1)
    temp.replace(path)


def versions():
    names = ("scanpy", "anndata", "numpy", "scipy", "h5py", "pandas",
             "harmonypy", "igraph", "leidenalg", "scikit-misc", "scikit-learn")
    result = {}
    for name in names:
        try:
            result[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            result[name] = "not-installed"
    return result


def verify_mapmycells_original_query(query_root, sample, metadata, gene_ids, original_counts):
    """Verify old annotation was calculated from exactly these integer counts.

    One capture at a time; this reads a saved query and does not invoke the mapper.
    Newly computed hashes and exact-count evidence enter the new run's receipt.
    """
    folder = Path(query_root) / sample
    query_path = folder / "query_raw_counts.h5ad"
    config_path = folder / "config.json"
    completion_path = folder / "COMPLETE.json"
    config_hash = sha256(config_path)
    config = json.loads(config_path.read_text())
    if config.get("type_assignment", {}).get("normalization") != "raw":
        raise ValueError(f"{sample}: reused MapMyCells mapping was not configured for raw counts")
    if Path(config.get("query_path", "")).resolve() != query_path.resolve():
        raise ValueError(f"{sample}: mapper config does not name the verified original query")
    completion = json.loads(completion_path.read_text())
    if completion.get("sample") != sample or int(completion.get("n_nuclei", -1)) != len(metadata):
        raise ValueError(f"{sample}: old MapMyCells completion identity differs from the authoritative roster")
    if len(gene_ids) != 32285:
        raise ValueError(f"{sample}: original annotation reuse requires the exact 32,285-gene feature axis")
    query_hash = sha256(query_path)
    query = ad.read_h5ad(query_path, backed="r")
    try:
        if query.shape != original_counts.shape or query.var_names.tolist() != gene_ids:
            raise ValueError(f"{sample}: original MapMyCells query shape/gene IDs differ from current original counts")
        if query.obs_names.has_duplicates or set(query.obs_names) != set(metadata.nucleus_id):
            raise ValueError(f"{sample}: original MapMyCells query nucleus IDs differ from the authoritative roster")
        positions = query.obs_names.get_indexer(metadata.nucleus_id)
        query_counts = query.X[positions, :].tocsr()
        if np.any(query_counts.data < 0) or np.any(query_counts.data != np.floor(query_counts.data)):
            raise ValueError(f"{sample}: original MapMyCells query is not nonnegative integer counts")
        query_counts.sum_duplicates()
        query_counts.sort_indices()
        unequal = query_counts != original_counts
        if unequal.nnz != 0:
            raise ValueError(f"{sample}: reused MapMyCells query counts differ at {unequal.nnz} entries")
        verified_shape = list(query_counts.shape)
        del query_counts, unequal
    finally:
        query.file.close()
    if sha256(query_path) != query_hash or sha256(config_path) != config_hash:
        raise ValueError(f"{sample}: MapMyCells query/config changed during verification")
    return {"status": "VERIFIED_EXACT_ORIGINAL_INTEGER_COUNTS", "sample": sample,
            "query_path": str(query_path), "query_sha256": query_hash,
            "mapper_config_sha256": config_hash, "old_completion_sha256": sha256(completion_path),
            "normalization": "raw", "shape": verified_shape, "gene_IDs_exact": True,
            "nucleus_IDs_exact": True, "unequal_sparse_count_entries": 0,
            "mapper_fit_or_rerun_performed": False}


def load_original_counts(cfg):
    """Read the declared 10x GEX count branch, applying its explicit frozen roster."""
    branch = cfg.get("count_branch", "original_full")
    is_original = branch in ("original_full", "paired_original")
    sheet = read_table(cfg["sample_sheet"])
    roster = read_table(cfg["qc_roster"])
    required = {"sample", "condition", "source_h5", "source_h5_sha256"}
    if not required.issubset(sheet):
        raise ValueError(f"Sample sheet missing {required - set(sheet)}")
    if sheet["sample"].duplicated().any() or set(sheet.condition) != {"WT", "cKO"}:
        raise ValueError("Sample sheet needs unique sample IDs and WT/cKO conditions")
    if len(sheet) != cfg.get("expected_samples", 12):
        raise ValueError("Sample sheet does not contain the expected mice")
    condition_counts = sheet.groupby("condition", observed=True).size().to_dict()
    if condition_counts != {"WT": 6, "cKO": 6}:
        raise ValueError("Expected six WT and six cKO biological replicates")
    if not {"nucleus_id", "sample", "barcode", "condition"}.issubset(roster):
        raise ValueError("Roster lacks identity fields")
    keep_column = cfg.get("keep_column", "passes_primary_QC")
    if keep_column in roster:
        roster = roster.loc[boolean(roster[keep_column], keep_column)].copy()
    if roster.nucleus_id.duplicated().any() or roster.nucleus_id.isna().any():
        raise ValueError("Duplicate or missing roster nucleus IDs")
    if not (roster.nucleus_id == roster["sample"] + "_" + roster.barcode).all():
        raise ValueError("Preserve nucleus_id = sample + '_' + original barcode")
    if set(roster["sample"]) != set(sheet["sample"]):
        raise ValueError("Roster and sample sheet do not cover the same mice")
    if not (roster.condition == roster["sample"].map(sheet.set_index("sample").condition)).all():
        raise ValueError("Roster conditions conflict with sample sheet")
    if len(roster) != cfg.get("expected_nuclei", 108062):
        raise ValueError("Roster differs from the explicit expected nucleus count for this branch")
    # All metadata are retained, but no earlier clustering is read as a new result.
    forbidden = [c for c in roster if c.lower().startswith(("seurat_clusters", "leiden", "umap"))]
    if forbidden:
        roster = roster.drop(columns=forbidden)
    matrices, metadata, source_records = [], [], []
    reference_ids = reference_symbols = None
    for row in sheet.itertuples(index=False):
        p = Path(row.source_h5)
        observed_hash = sha256(p)
        if observed_hash != row.source_h5_sha256:
            raise ValueError(f"Original source checksum mismatch for {row.sample}")
        md = roster.loc[roster["sample"] == row.sample].copy()
        with h5py.File(p, "r") as f:
            g = f["matrix"]
            ids = g["features/id"].asstr()[()].tolist()
            symbols = g["features/name"].asstr()[()].tolist()
            kinds = g["features/feature_type"].asstr()[()].tolist()
            barcodes = g["barcodes"].asstr()[()].tolist()
            x = sparse.csc_matrix((g["data"][()], g["indices"][()], g["indptr"][()]),
                                  shape=tuple(g["shape"][()]))
        if len(set(ids)) != len(ids) or len(set(barcodes)) != len(barcodes):
            raise ValueError(f"Duplicate original gene ID/barcode in {row.sample}")
        if set(kinds) != {"Gene Expression"}:
            raise ValueError("Split non-GEX modalities explicitly before this RNA analysis")
        if reference_ids is None:
            reference_ids, reference_symbols = ids, symbols
        elif ids != reference_ids or symbols != reference_symbols:
            raise ValueError("Feature IDs and symbols must match exactly across all 12 inputs")
        lookup = pd.Index(barcodes).get_indexer(md.barcode)
        if (lookup < 0).any():
            raise ValueError(f"Frozen QC nuclei absent from source for {row.sample}")
        x = x[:, lookup].T.tocsr()
        x.sort_indices()
        if np.any(x.data < 0) or np.any(x.data != np.floor(x.data)):
            raise ValueError("Original-count branch requires nonnegative integer UMIs")
        totals = np.asarray(x.sum(axis=1)).ravel().astype(np.int64)
        if np.any(totals <= 0):
            raise ValueError("Frozen QC cohort contains a zero-UMI nucleus")
        if is_original and "nCount_RNA" in md and not np.array_equal(totals, md.nCount_RNA.to_numpy()):
            raise ValueError(f"Original per-nucleus UMI totals disagree for {row.sample}")
        if cfg.get("verify_sample_original_umis", branch == "original_full") and hasattr(row, "original_UMIs") and int(totals.sum()) != int(row.original_UMIs):
            raise ValueError(f"Sample original UMI totals disagree for {row.sample}")
        if "nCount_RNA" in md:
            md["frozen_original_nCount_RNA"] = md["nCount_RNA"].to_numpy()
        md["branch_UMIs_verified"] = totals
        if branch == "original_full" and cfg.get("mapmycells_query_root"):
            query_verification = verify_mapmycells_original_query(
                cfg["mapmycells_query_root"], row.sample, md, ids, x)
        else:
            query_verification = {"status": "NOT_REQUESTED" if branch == "original_full" else "PAIRED_BRANCH_SHARED_ORIGINAL_LABEL_PROVENANCE",
                                  "source_label_branch": "original_full", "mapper_fit_or_rerun_performed": False}
        matrices.append(x)
        metadata.append(md)
        source_records.append({"sample": row.sample, "source": str(p), "sha256": observed_hash,
                               "n_nuclei": len(md), "UMIs": int(totals.sum()),
                               "mapmycells_query_verification": query_verification})
        print(f"Loaded {branch} counts: {row.sample}, {len(md):,} nuclei", flush=True)
    x = sparse.vstack(matrices, format="csr")
    del matrices
    obs = pd.concat(metadata, ignore_index=True).set_index("nucleus_id", drop=False)
    var = pd.DataFrame({"gene_symbol": reference_symbols}, index=pd.Index(reference_ids, name="gene_id"))
    for c in ("sample", "condition", "sequencing_partition"):
        if c in obs:
            obs[c] = pd.Categorical(obs[c])
    # Ensure AnnData does not receive pandas StringArray values from new pandas versions.
    for c in obs:
        if obs[c].dtype == object:
            obs[c] = obs[c].fillna("NA").astype(str)
    data = ad.AnnData(X=x, obs=obs, var=var)
    data.uns["count_branch"] = branch
    return data, source_records


def join_mapmycells(data, prediction_path, allow_superset=False):
    predictions = read_table(prediction_path)
    if predictions.nucleus_id.duplicated().any():
        raise ValueError("MapMyCells prediction IDs are not unique")
    predictions = predictions.set_index("nucleus_id")
    if (not allow_superset and set(predictions.index) != set(data.obs_names)) or not set(data.obs_names).issubset(predictions.index):
        raise ValueError("MapMyCells labels must cover exactly the frozen QC cohort")
    predictions = predictions.loc[data.obs_names]
    for c in ("sample", "condition"):
        if c in predictions and not (predictions[c].astype(str).to_numpy() == data.obs[c].astype(str).to_numpy()).all():
            raise ValueError(f"MapMyCells {c} conflicts with the authoritative roster")
    for c in predictions:
        if c not in ("sample", "condition", "barcode"):
            data.obs["mmc_" + c] = predictions[c].to_numpy()
    data.uns["MapMyCells_annotation_evidence"] = (
        "Verified original-count query annotation reused; this is common annotation evidence, "
        "not a new independent mapping or new clustering outcome.")


def population_mask(data, population):
    if population == "allcells":
        return np.ones(data.n_obs, dtype=bool)
    column = {"Excitatory": "REF_EX_HS", "Inhibitory": "REF_IN_HS", "CA1_ProS": "REF_CA1_PROS_HS"}[population]
    if column not in data.obs:
        raise ValueError(f"Authoritative frozen population flag {column} absent")
    return boolean(data.obs[column], column)
