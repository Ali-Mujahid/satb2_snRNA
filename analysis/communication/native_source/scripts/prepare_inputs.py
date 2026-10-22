#!/usr/bin/env python3
"""Prepare frozen Seurat-defined, full-gene-background LIANA input.

No fitting, imputation, cell filtering, or differential expression is performed.
The CLI refuses nonempty output directories and verifies configured input hashes.
All collapsed RNA genes are retained: native LIANA SingleCellSignalR computes its
background mean before its own ligand/receptor-resource filtering.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem
from scipy import sparse


PROGRAMS = {
    "CA1_Deep": "Deep-associated expression program enriched",
    "CA1_Superficial": "Superficial-associated expression program enriched",
}
ANNOTATION_FIELDS = [
    "nucleus_id", "sample", "condition", "tool", "branch", "population",
    "cluster_scope", "cluster_id", "source_cluster_id", "current_cluster_id",
    "keep_for_analysis", "keep_for_plot", "exclusion_reason", "source_count_branch",
    "broad_label", "reviewed_label", "identity_label", "label_confidence",
    "layer_program_label", "expression_program_label", "program_evidence_status",
    "physical_layer_assignment", "identity_review_status", "uncertainty_reason",
    "annotation_evidence_status", "annotation_method", "full_cohort_marker_status",
    "reference_identity_review_status", "reference_broad_high_support",
    "reference_subclass_high_support", "reference_subclass_name",
]
REQUIRED_FIELDS = {
    "nucleus_id", "sample", "condition", "tool", "branch", "population",
    "cluster_scope", "cluster_id", "keep_for_analysis", "layer_program_label",
    "reviewed_label", "identity_review_status",
}


def require(ok, message):
    if not bool(ok):
        raise ValueError(message)


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path, data):
    Path(path).write_text(json.dumps(data, indent=2, sort_keys=True, allow_nan=False) + "\n")


def checked_path(value):
    p = Path(value)
    require(p.is_absolute(), f"Expected an absolute path: {p}")
    require(p.is_file(), f"Missing input file: {p}")
    return p


def configured_files(config):
    records = {}
    for name, item in config["inputs"].items():
        path = checked_path(item["path"])
        require(bool(item.get("sha256")), f"Missing input hash: {name}")
        records[f"input:{name}"] = {
            "path": str(path), "expected_sha256": item["sha256"],
            "expected_bytes": item.get("bytes"),
        }
    for ep in config["endpoints"]:
        tested = ep["status"] != "NOT_TESTED"
        for field in ("all_results_path", "deg_path", "stable_audit_path"):
            value = ep.get(field)
            expected = ep.get(field + "_sha256")
            if not value or not Path(value).is_file():
                require(not tested, f"Tested endpoint missing {field}: {ep['population_key']}")
                continue
            require(bool(expected), f"Existing endpoint file is not hash bound: {value}")
            records[f"endpoint:{ep['population_key']}:{field}"] = {
                "path": str(checked_path(value)), "expected_sha256": expected,
                "expected_bytes": None,
            }
    resource = checked_path(config["resource_file"])
    # The immutable reference can be bound directly or through the inputs map.
    configured_resource_hash = config.get("resource_sha256") or config.get("resource_file_sha256")
    records["reference:resource"] = {
        "path": str(resource), "expected_sha256": configured_resource_hash,
        "expected_bytes": None,
    }
    return records


def verify_hashes(records, phase):
    for key, record in records.items():
        path = Path(record["path"])
        if record["expected_bytes"] is not None:
            require(path.stat().st_size == int(record["expected_bytes"]), f"Input size changed: {key}")
        actual = sha256(path)
        if record["expected_sha256"]:
            require(actual == record["expected_sha256"], f"Input hash mismatch: {key}")
        if phase == "after":
            require(actual == record["before_sha256"], f"Input changed during preparation: {key}")
        record[phase + "_sha256"] = actual
        record["bytes"] = path.stat().st_size


def read_tsv_chunks(path, columns=None):
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False,
                       usecols=columns, chunksize=100_000)


def select_annotations(path):
    header = set(pd.read_csv(path, sep="\t", nrows=0).columns)
    require(REQUIRED_FIELDS <= header, f"Missing annotation fields: {REQUIRED_FIELDS - header}")
    fields = [c for c in ANNOTATION_FIELDS if c in header]
    selected = []
    for chunk in read_tsv_chunks(path, fields):
        common = ((chunk.tool == "seurat") & (chunk.branch == "original_full") &
                  (chunk.cluster_scope == "original_frozen") &
                  chunk.keep_for_analysis.str.lower().isin(["true", "1"]))
        in_scope = (chunk.population == "Inhibitory") | (
            (chunk.population == "CA1_ProS") & chunk.layer_program_label.isin(PROGRAMS.values()))
        subset = chunk.loc[common & in_scope].copy()
        if not subset.empty:
            selected.append(subset)
    require(bool(selected), "No requested Seurat nuclei were selected")
    result = pd.concat(selected, ignore_index=True)
    require(not result.nucleus_id.duplicated().any(), "Overlapping/duplicate target nucleus IDs")
    result["source_cluster"] = pd.to_numeric(result.cluster_id, errors="raise").astype(int).astype(str)
    result["ccc_group"] = ""
    for group, label in PROGRAMS.items():
        result.loc[(result.population == "CA1_ProS") & (result.layer_program_label == label), "ccc_group"] = group
    in_mask = result.population == "Inhibitory"
    result.loc[in_mask, "ccc_group"] = "IN_C" + result.loc[in_mask, "source_cluster"]
    require((result.ccc_group != "").all(), "Unassigned communication groups")
    require(set(result.loc[in_mask, "source_cluster"]) == {str(i) for i in range(1, 21)},
            "Expected exactly Seurat inhibitory clusters C1–C20")
    return result


def validate_categories(path, annotations):
    fields = ["nucleus_id", "sample", "condition", "tool", "branch", "population", "category_family", "label"]
    selected = []
    for chunk in read_tsv_chunks(path, fields):
        common = (chunk.tool == "seurat") & (chunk.branch == "original_full")
        ca = ((chunk.population == "CA1_ProS") &
              (chunk.category_family == "CA1_LAYER_EXPRESSION_PROGRAM") & chunk.label.isin(PROGRAMS.values()))
        inhib = ((chunk.population == "Inhibitory") &
                 (chunk.category_family == "REVIEWED_COARSE") & (chunk.label == "Inhibitory"))
        subset = chunk.loc[common & (ca | inhib)].copy()
        if not subset.empty:
            selected.append(subset)
    require(bool(selected), "No matching category memberships")
    categories = pd.concat(selected, ignore_index=True)
    require(not categories.nucleus_id.duplicated().any(), "Duplicated requested category memberships")
    require(set(categories.nucleus_id) == set(annotations.nucleus_id), "Category and nucleus annotation memberships differ")
    aligned = categories.set_index("nucleus_id").loc[annotations.nucleus_id]
    for field in ("sample", "condition", "population"):
        require(np.array_equal(aligned[field].to_numpy(), annotations[field].to_numpy()),
                f"Category metadata differs from annotation: {field}")
    for group, label in PROGRAMS.items():
        require(set(categories.loc[categories.label == label, "nucleus_id"]) ==
                set(annotations.loc[annotations.ccc_group == group, "nucleus_id"]),
                f"Frozen category membership mismatch: {group}")


def endpoint_group(ep):
    if ep["analysis_kind"] == "total" and "Inhibitory" in ep["population_key"]:
        return "IN_Total"
    if ep["analysis_kind"] == "cluster" and "Inhibitory" in ep["population_key"]:
        return "IN_C" + str(int(ep["analysis_label"]))
    for group, label in PROGRAMS.items():
        if ep["analysis_kind"] == "layer_program" and ep["analysis_label"] == label:
            return group
    raise ValueError(f"Unexpected endpoint in scoped configuration: {ep['population_key']}")


def group_metadata(config, annotations):
    by_group = {}
    for ep in config["endpoints"]:
        group = endpoint_group(ep)
        require(group not in by_group, f"Multiple endpoints for {group}")
        by_group[group] = ep
    required = set(PROGRAMS) | {f"IN_C{i}" for i in range(1, 21)} | {"IN_Total"}
    require(set(by_group) == required, "Endpoint map must cover both CA1 programs, all 20 IN clusters and IN total")
    rows = []
    for group in list(PROGRAMS) + [f"IN_C{i}" for i in range(1, 21)] + ["IN_Total"]:
        ep = by_group[group]
        subset = annotations.loc[annotations.population == "Inhibitory"] if group == "IN_Total" else annotations.loc[annotations.ccc_group == group]
        require(len(subset) == int(ep["n_nuclei"]), f"Frozen endpoint count mismatch: {group}")
        unique = lambda field: "; ".join(sorted(set(subset[field]) - {""})) if field in subset else ""
        rows.append({
            "group": group,
            "label": PROGRAMS.get(group, "Inhibitory total" if group == "IN_Total" else unique("reviewed_label")),
            "identity_review": unique("identity_review_status"),
            "endpoint_key": ep["population_key"],
            "n_nuclei": len(subset),
            "source_clusters": unique("source_cluster"),
            "reviewed_label": unique("reviewed_label"),
            "expression_program_label": unique("expression_program_label"),
            "label_confidence": unique("label_confidence"),
            "uncertainty_reason": unique("uncertainty_reason"),
            "physical_layer_assignment": unique("physical_layer_assignment"),
            "endpoint_status": ep["status"], "endpoint_reason": ep.get("reason", ""),
            "supplemental_resolution": group == "IN_Total",
        })
    return pd.DataFrame(rows)


def strings(values):
    return np.asarray(["" if pd.isna(v) else (v.decode() if isinstance(v, bytes) else str(v)) for v in values], dtype=object)


def resource_features(resource, raw_symbols, stable_ids):
    require({"ligand", "receptor"} <= set(resource.columns), "Resource must contain ligand and receptor columns")
    require(len(resource) > 0, "Empty ligand/receptor resource")
    require(resource[["ligand", "receptor"]].ne("").all().all(), "Blank ligand/receptor names")
    # LIANA's frozen protein complexes use '_' between required subunits.
    components = lambda value: list(dict.fromkeys(str(value).split("_")))
    resource_genes = set(gene for col in ("ligand", "receptor") for value in resource[col] for gene in components(value))
    symbols = [symbol.strip() if symbol.strip() else stable for symbol, stable in zip(raw_symbols, stable_ids)]
    require(len(set(stable_ids)) == len(stable_ids), "Duplicate stable feature IDs")
    all_symbols = list(dict.fromkeys(symbols))
    represented = set(all_symbols)
    overlapping_resource_symbols = [s for s in all_symbols if s in resource_genes]
    require(bool(overlapping_resource_symbols), "No resource genes overlap native features")
    # Keep the full collapsed universe for native SingleCellSignalR's mat_mean.
    # LIANA performs its resource restriction internally, after that background.
    retained_symbols = all_symbols
    symbol_column = {s: i for i, s in enumerate(retained_symbols)}
    feature = pd.DataFrame({"source_feature_index": np.arange(len(symbols)), "gene_id": stable_ids,
                            "source_gene_symbol": raw_symbols, "analysis_symbol": symbols})
    feature["missing_symbol_fallback"] = [not str(s).strip() for s in raw_symbols]
    sizes = feature.groupby("analysis_symbol", sort=False).size()
    feature["n_source_features_for_symbol"] = feature.analysis_symbol.map(sizes).astype(int)
    feature["duplicate_symbol"] = feature.n_source_features_for_symbol > 1
    feature["in_resource"] = feature.analysis_symbol.isin(resource_genes)
    feature["retained_in_h5ad"] = feature.analysis_symbol.isin(retained_symbols)
    feature["output_feature_index"] = feature.analysis_symbol.map(symbol_column).fillna(-1).astype(int)
    coverage = resource.copy()
    coverage.insert(0, "resource_row", np.arange(1, len(coverage) + 1))
    for role in ("ligand", "receptor"):
        subs = coverage[role].map(components)
        coverage[role + "_subunits"] = subs.map(lambda x: ";".join(x))
        coverage[role + "_missing_subunits"] = subs.map(lambda x: ";".join(g for g in x if g not in represented))
        coverage[role + "_all_subunits_present"] = coverage[role + "_missing_subunits"] == ""
    coverage["all_subunits_present"] = coverage.ligand_all_subunits_present & coverage.receptor_all_subunits_present
    return feature, coverage, retained_symbols


def read_sparse_counts(config, annotations, feature, retained_symbols):
    counts_path = config["inputs"]["counts"]["path"]
    chunks, libraries, source_indices = [], [], []
    stats = {"source_n_rows": 0, "source_n_features": len(feature), "source_count_dtype": "",
             "selected_umi_total_full_library": 0, "selected_umi_total_collapsed_all_genes": 0,
             "selected_umi_total_resource_features": 0,
             "source_raw_values_checked": 0, "library_sizes_verified_against_nCount_RNA": False}
    with h5py.File(counts_path, "r") as handle:
        source_ids = strings(read_elem(handle["obs/nucleus_id"]))
        require(len(set(source_ids)) == len(source_ids), "Duplicate native nucleus IDs")
        annotation_ids = set(annotations.nucleus_id)
        selected = np.asarray([v in annotation_ids for v in source_ids])
        require(int(selected.sum()) == len(annotations), "Not all exact target nuclei exist in source count matrix")
        ordered_ids = source_ids[selected]
        obs = annotations.set_index("nucleus_id", drop=False).loc[ordered_ids].copy()
        source_sample = strings(read_elem(handle["obs/sample"]))
        source_condition = strings(read_elem(handle["obs/condition"]))
        for field, native in [("sample", source_sample), ("condition", source_condition)]:
            require(np.array_equal(native[selected], obs[field].to_numpy()), f"Native {field} differs after exact-ID join")
        matrix = handle["X"]
        encoding = matrix.attrs.get("encoding-type", "")
        if isinstance(encoding, bytes):
            encoding = encoding.decode()
        require(encoding == "csr_matrix", "Expected native CSR count matrix")
        shape = tuple(int(v) for v in matrix.attrs["shape"])
        require(shape == (len(source_ids), len(feature)), "Count and feature metadata shape mismatch")
        require(np.issubdtype(matrix["data"].dtype, np.integer), "Native count matrix must store integers")
        stats.update(source_n_rows=shape[0], source_count_dtype=str(matrix["data"].dtype))
        keep_features = feature.retained_in_h5ad.to_numpy()
        aggregation = sparse.csr_matrix(
            (np.ones(int(keep_features.sum()), dtype=np.int64),
             (feature.loc[keep_features, "source_feature_index"].to_numpy(),
              feature.loc[keep_features, "output_feature_index"].to_numpy())),
            shape=(shape[1], len(retained_symbols)))
        resource_columns = feature.loc[feature.in_resource, "output_feature_index"].unique()
        chunk_size = 4096
        for start in range(0, shape[0], chunk_size):
            stop = min(start + chunk_size, shape[0])
            ptr = matrix["indptr"][start:stop + 1].astype(np.int64)
            low, high = int(ptr[0]), int(ptr[-1])
            data = matrix["data"][low:high]
            require((data >= 0).all(), f"Negative native count in block starting {start}")
            stats["source_raw_values_checked"] += int(len(data))
            if not selected[start:stop].any():
                continue
            block = sparse.csr_matrix((data.astype(np.int64), matrix["indices"][low:high], ptr - low),
                                      shape=(stop - start, shape[1]))
            block = block[selected[start:stop]]
            full_library = np.asarray(block.sum(axis=1)).ravel().astype(np.int64)
            require((full_library > 0).all(), "Zero-library target nucleus")
            collapsed = (block @ aggregation).tocsr()
            require(np.issubdtype(collapsed.dtype, np.integer), "Symbol collapsing lost integer counts")
            require(np.array_equal(np.asarray(collapsed.sum(axis=1)).ravel(), full_library),
                    "Full-gene symbol aggregation changed a nucleus's integer library")
            stats["selected_umi_total_full_library"] += int(full_library.sum())
            stats["selected_umi_total_collapsed_all_genes"] += int(collapsed.sum())
            stats["selected_umi_total_resource_features"] += int(collapsed[:, resource_columns].sum())
            normalized = collapsed.astype(np.float64).multiply((10_000.0 / full_library)[:, None]).tocsr()
            np.log1p(normalized.data, out=normalized.data)
            require(np.isfinite(normalized.data).all() and (normalized.data >= 0).all(), "Invalid normalized expression")
            chunks.append(normalized.astype(np.float32))
            libraries.append(full_library)
            source_indices.extend(np.flatnonzero(selected[start:stop]) + start)
        require(bool(chunks), "No target count chunks read")
        X = sparse.vstack(chunks, format="csr")
        library_sizes = np.concatenate(libraries)
        require(X.shape[0] == len(obs), "Row count changed during sparse preparation")
        if "nCount_RNA" in handle["obs"]:
            reported = np.asarray(read_elem(handle["obs/nCount_RNA"]), dtype=np.float64)[selected]
            require(np.array_equal(reported, library_sizes.astype(np.float64)), "Native nCount_RNA does not match complete integer library sums")
            stats["library_sizes_verified_against_nCount_RNA"] = True
        obs["source_row_index"] = np.asarray(source_indices, dtype=np.int64)
        obs["library_umi_full_genes"] = library_sizes
        obs["log1p_cp10k_scale_factor"] = 10_000.0 / library_sizes
    return X, obs, stats


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()
    config_path = checked_path(args.config)
    outdir = Path(args.outdir)
    require(outdir.is_absolute(), "--outdir must be absolute")
    config = json.loads(config_path.read_text())
    run_root = Path(config["run_root"])
    require(run_root.is_absolute() and outdir.is_relative_to(run_root), "Output must be inside the new configured run root")
    require(not outdir.exists() or (outdir.is_dir() and not any(outdir.iterdir())), "Refusing nonempty output directory")
    outdir.mkdir(parents=True, exist_ok=True)
    started = datetime.now(timezone.utc).isoformat()
    records = configured_files(config)
    records["configuration"] = {"path": str(config_path), "expected_sha256": sha256(config_path), "expected_bytes": config_path.stat().st_size}
    verify_hashes(records, "before")
    print("Verified configured sources and endpoint hashes", flush=True)
    annotations = select_annotations(config["inputs"]["nucleus_annotations"]["path"])
    require(len(annotations) == int(config["expected_nuclei"]), "Target nucleus count differs from plan")
    for group, n in config["expected_group_counts"].items():
        observed = (annotations.population == "Inhibitory").sum() if group == "IN_Total" else (annotations.ccc_group == group).sum()
        require(int(observed) == int(n), f"Target group count differs from plan: {group}")
    validate_categories(config["inputs"]["category_memberships"]["path"], annotations)
    samples = pd.read_csv(config["inputs"]["samples"]["path"], sep="\t", dtype=str, keep_default_na=False)
    require(not samples["sample"].duplicated().any(), "Sample roster has duplicates")
    require(set(samples["sample"]) == set(annotations["sample"]), "Target nuclei do not represent the exact 12-mouse roster")
    require(samples.groupby("condition")["sample"].nunique().to_dict() == {"WT": 6, "cKO": 6}, "Expected six independent WT and six cKO mice")
    sample_to_condition = samples.set_index("sample").condition
    require((annotations["sample"].map(sample_to_condition) == annotations.condition).all(), "Genotype differs from sample roster")
    groups = group_metadata(config, annotations)
    resource = pd.read_csv(config["resource_file"], sep="\t", dtype=str, keep_default_na=False)
    with h5py.File(config["inputs"]["counts"]["path"], "r") as handle:
        stable_ids = strings(read_elem(handle["var/gene_id"]))
        raw_symbols = strings(read_elem(handle["var/gene_symbol"]))
    feature, coverage, symbols = resource_features(resource, raw_symbols, stable_ids)
    print(f"Exact memberships verified: {len(annotations):,} nuclei; {len(symbols):,} full-background symbols", flush=True)
    X, obs, stats = read_sparse_counts(config, annotations, feature, symbols)
    group_order = list(PROGRAMS) + [f"IN_C{i}" for i in range(1, 21)]
    for col in obs.columns:
        if pd.api.types.is_object_dtype(obs[col].dtype):
            obs[col] = obs[col].fillna("").astype(str)
    obs["ccc_group"] = pd.Categorical(obs.ccc_group, categories=group_order, ordered=True)
    obs["condition"] = pd.Categorical(obs.condition, categories=["WT", "cKO"])
    obs["sample"] = pd.Categorical(obs["sample"], categories=samples["sample"].tolist())
    mapped = feature.groupby("analysis_symbol", sort=False)
    var = pd.DataFrame(index=pd.Index(symbols, name="gene_symbol"))
    var["gene_ids"] = [";".join(mapped.get_group(s).gene_id) for s in symbols]
    var["n_source_features"] = [len(mapped.get_group(s)) for s in symbols]
    var["duplicate_symbol"] = var.n_source_features > 1
    var["resource_gene_flag"] = [bool(mapped.get_group(s).in_resource.any()) for s in symbols]
    output = ad.AnnData(X=X, obs=obs, var=var)
    output.uns["normalization"] = {"method": "log1p(CP10000)", "target_sum": 10000,
                                   "denominator": "complete native integer RNA library",
                                   "duplicate_symbols": "sum integer features before normalization and log1p",
                                   "feature_universe": "all collapsed native RNA genes retained",
                                   "resource_filtering": "performed internally by native LIANA after background calculation",
                                   "imputation": "none", "expression_is_integrated": False}
    output.uns["source_counts_sha256"] = records["input:counts"]["before_sha256"]
    output.write_h5ad(outdir / "lr_expression.h5ad", compression="gzip")
    feature.to_csv(outdir / "feature_mapping.tsv", sep="\t", index=False)
    ambiguities = feature.loc[feature.duplicate_symbol].copy()
    ambiguities["DEG_support_policy"] = "AMBIGUOUS_SYMBOL_EXCLUDED_FROM_PRIMARY_DEG_SUPPORT"
    ambiguities.to_csv(outdir / "symbol_ambiguities.tsv", sep="\t", index=False)
    coverage.to_csv(outdir / "resource_coverage.tsv", sep="\t", index=False)
    groups.to_csv(outdir / "group_metadata.tsv", sep="\t", index=False)
    all_rows = []
    for sample in samples["sample"]:
        condition = sample_to_condition.loc[sample]
        for group in group_order + ["IN_Total"]:
            selected = (obs.population == "Inhibitory") if group == "IN_Total" else (obs.ccc_group == group)
            count = int(((obs["sample"] == sample) & selected).sum())
            all_rows.append({"sample": sample, "condition": condition, "group": group, "n_nuclei": count,
                             "eligible_mouse_support": count >= int(config["min_cells"]),
                             "supplemental_resolution": group == "IN_Total"})
    pd.DataFrame(all_rows).to_csv(outdir / "nuclei_per_mouse.tsv", sep="\t", index=False)
    pd.DataFrame({"nucleus_id": obs.nucleus_id.astype(str), "sample": obs["sample"].astype(str),
                  "condition": obs.condition.astype(str), "ccc_group": obs.ccc_group.astype(str),
                  "source_row_index": obs.source_row_index, "library_umi_full_genes": obs.library_umi_full_genes,
                  "log1p_cp10k_scale_factor": obs.log1p_cp10k_scale_factor}).to_csv(
                      outdir / "nucleus_library_sizes.tsv.gz", sep="\t", index=False, compression="gzip")
    stats.update({"method": "log1p(CP10000)", "target_sum": 10000,
                  "library_denominator": f"all {len(feature)} native features",
                  "n_selected_nuclei": int(output.n_obs), "n_retained_symbols": int(output.n_vars),
                  "n_resource_symbols": int(var.resource_gene_flag.sum()),
                  "output_sparse_format": "csr", "output_dtype": str(X.dtype),
                  "source_features_with_duplicate_symbol": int(feature.duplicate_symbol.sum()),
                  "ambiguous_symbols": int(ambiguities.analysis_symbol.nunique()),
                  "resource_interactions": int(len(coverage)),
                  "resource_interactions_all_subunits_present": int(coverage.all_subunits_present.sum()),
                  "full_collapsed_gene_background_retained": True,
                  "resource_subsetting_in_preparation": False,
                  "native_SingleCellSignalR_mat_mean_background": "all collapsed normalized RNA genes",
                  "no_imputation": True, "native_endpoint_refits": 0})
    write_json(outdir / "expression_normalization.json", stats)
    verify_hashes(records, "after")
    write_json(outdir / "source_hashes.json", records)
    files = {p.name: {"bytes": p.stat().st_size, "sha256": sha256(p)}
             for p in sorted(outdir.iterdir()) if p.is_file()}
    complete = {"status": "COMPLETE", "workflow_version": config["workflow_version"],
                "started_utc": started, "finished_utc": datetime.now(timezone.utc).isoformat(),
                "config_sha256": records["configuration"]["before_sha256"],
                "n_nuclei": int(output.n_obs), "n_genes": int(output.n_vars),
                "group_counts": {g: int((obs.ccc_group == g).sum()) for g in group_order},
                "supplemental_IN_Total_nuclei": int((obs.population == "Inhibitory").sum()),
                "prior_sources_unchanged": True, "files": files}
    write_json(outdir / "COMPLETE.json", complete)
    print(f"Preparation complete: {outdir}", flush=True)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"Preparation failed: {type(exc).__name__}: {exc}", file=sys.stderr, flush=True)
        raise
