#!/usr/bin/env python3
"""Prepare exact frozen Seurat whole-transcriptome inputs for a registered workflow.

Implementation draft. Creating this source does not execute preparation. The CLI
is intended to be called only by the approved Snakemake workflow. It never fits a
model, infers a regulon, selects DEGs/HVGs/LR genes, integrates expression, changes
an existing annotation, or overwrites a nonempty output directory.

Config contract and output schemas: docs/PREPARE_INPUTS_CONTRACT.md.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path
from datetime import datetime, timezone

import anndata as ad
from anndata.io import read_elem
import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmwrite


SCHEMA = "satb2-regulatory-preparation/v1"
ENDPOINTS = {
    "CA1_Deep": {"key": "P0244_seurat_CA1_ProS_LAYER_PROGRAM", "population": "CA1_ProS", "label": "Deep-associated expression program enriched", "clusters": {"1", "5"}, "n": 8750, "lineage": "CA1"},
    "CA1_Superficial": {"key": "P0247_seurat_CA1_ProS_LAYER_PROGRAM", "population": "CA1_ProS", "label": "Superficial-associated expression program enriched", "clusters": {"2", "3", "8"}, "n": 11706, "lineage": "CA1"},
    "IN_C1": {"key": "P0140_seurat_Inhibitory_CLUSTER", "population": "Inhibitory", "label": "1", "clusters": {"1"}, "n": 1718, "lineage": "IN"},
    "IN_C6": {"key": "P0145_seurat_Inhibitory_CLUSTER", "population": "Inhibitory", "label": "6", "clusters": {"6"}, "n": 501, "lineage": "IN"},
    "IN_C8": {"key": "P0147_seurat_Inhibitory_CLUSTER", "population": "Inhibitory", "label": "8", "clusters": {"8"}, "n": 468, "lineage": "IN"},
}
EXPECTED_LINEAGES = {"CA1": 20456, "IN": 2687}
REQUIRED_INPUTS = {"counts", "nucleus_annotations", "category_memberships", "population_index", "population_memberships", "samples"} | {"endpoint_counts_" + group for group in ENDPOINTS}
ANNOTATION_REQUIRED = {"nucleus_id", "sample", "condition", "tool", "branch", "population", "cluster_scope", "cluster_id", "keep_for_analysis", "layer_program_label", "source_count_branch"}
ANNOTATION_EXTRA = ["source_cluster_id", "current_cluster_id", "keep_for_plot", "exclusion_reason", "broad_label", "reviewed_label", "identity_label", "label_confidence", "expression_program_label", "program_evidence_status", "physical_layer_assignment", "identity_review_status", "uncertainty_reason", "annotation_evidence_status", "annotation_method", "full_cohort_marker_status", "reference_identity_review_status", "reference_broad_high_support", "reference_subclass_high_support", "reference_subclass_name"]


def require(ok, message):
    if not bool(ok):
        raise ValueError(message)


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def membership_sha(ids):
    values = [str(value) for value in ids]
    require(len(values) == len(set(values)), "Cannot hash duplicate nucleus IDs")
    return hashlib.sha256("".join(value + "\n" for value in sorted(values)).encode("utf-8")).hexdigest()


def strings(values):
    return np.asarray(["" if pd.isna(value) else value.decode() if isinstance(value, bytes) else str(value) for value in values], dtype=object)


def read_table(path, **kwargs):
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, **kwargs)


def table_chunks(path, columns):
    return read_table(path, usecols=columns, chunksize=100_000)


def json_write(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")


def table_write(table, path):
    path = Path(path)
    compression = {"method": "gzip", "mtime": 0} if path.suffix == ".gz" else None
    table.to_csv(path, sep="\t", index=False, compression=compression)


def checked_path(value):
    path = Path(value)
    require(path.is_absolute() and path.is_file(), f"Expected existing absolute input file: {path}")
    return path.resolve()


def input_records(config, config_path):
    require(config.get("schema") == SCHEMA, f"Expected configuration schema {SCHEMA}")
    require(REQUIRED_INPUTS <= set(config.get("inputs", {})), f"Missing configured inputs: {REQUIRED_INPUTS - set(config.get('inputs', {}))}")
    expected = config.get("expected_group_counts", {group: value["n"] for group, value in ENDPOINTS.items()})
    require(expected == {group: value["n"] for group, value in ENDPOINTS.items()}, "This frozen analysis requires DL8750/SL11706/C1 1718/C6 501/C8 468; changed membership requires a new plan")
    records = {}
    for name, item in config["inputs"].items():
        path = checked_path(item["path"])
        require(len(item.get("sha256", "")) == 64, f"Input must be hash-bound: {name}")
        records[name] = {"path": str(path), "expected_sha256": item["sha256"], "expected_bytes": item.get("bytes")}
    records["configuration"] = {"path": str(config_path), "expected_sha256": sha256(config_path), "expected_bytes": config_path.stat().st_size}
    return records


def verify_records(records, phase):
    for name, item in records.items():
        path = Path(item["path"])
        if item.get("expected_bytes") is not None:
            require(path.stat().st_size == int(item["expected_bytes"]), f"Input size changed: {name}")
        digest = sha256(path)
        require(digest == item["expected_sha256"], f"Input hash mismatch: {name}")
        if phase == "after":
            require(digest == item["before_sha256"], f"Input changed during preparation: {name}")
        item[phase + "_sha256"] = digest
        item["bytes"] = path.stat().st_size


def select_annotations(path):
    """Keep native flags unchanged; count requested members before and after them."""
    header = set(read_table(path, nrows=0).columns)
    require(ANNOTATION_REQUIRED <= header, f"Missing annotation columns: {ANNOTATION_REQUIRED - header}")
    columns = sorted(ANNOTATION_REQUIRED) + [field for field in ANNOTATION_EXTRA if field in header]
    candidate_parts = []
    for chunk in table_chunks(path, columns):
        frozen = chunk.loc[(chunk.tool == "seurat") & (chunk.branch == "original_full") & (chunk.cluster_scope == "original_frozen")].copy()
        if frozen.empty:
            continue
        frozen["endpoint"] = ""
        for group, spec in ENDPOINTS.items():
            mask = frozen.population.eq(spec["population"])
            mask &= frozen.layer_program_label.eq(spec["label"]) if spec["lineage"] == "CA1" else frozen.cluster_id.eq(spec["label"])
            frozen.loc[mask, "endpoint"] = group
        candidate = frozen.loc[frozen.endpoint.ne("")].copy()
        if not candidate.empty:
            candidate_parts.append(candidate)
    require(candidate_parts, "No requested frozen annotation rows")
    candidates = pd.concat(candidate_parts, ignore_index=True)
    require(not candidates.nucleus_id.duplicated().any(), "Requested frozen populations overlap or contain duplicate nucleus IDs")
    candidates["selected_native_keep_flag"] = candidates.keep_for_analysis.str.lower().isin(["true", "1", "yes"])
    selected = candidates.loc[candidates.selected_native_keep_flag].copy()
    require(len(selected) == 23143, "Expected exact five-endpoint total 23,143 nuclei")
    require(selected.source_count_branch.eq("original_full").all(), "Unexpected corrected/integrated count branch")
    rows = []
    for group, spec in ENDPOINTS.items():
        before = candidates.loc[candidates.endpoint.eq(group)]
        after = selected.loc[selected.endpoint.eq(group)]
        require(len(after) == spec["n"], f"Frozen retained count differs: {group}")
        require(set(after.cluster_id) == spec["clusters"], f"Frozen cluster composition differs: {group}")
        rows.append({"endpoint": group, "population_key": spec["key"], "lineage": spec["lineage"], "candidate_nuclei_before_native_keep_flag": len(before), "retained_nuclei": len(after), "excluded_by_existing_native_keep_flag": len(before) - len(after), "membership_sha256": membership_sha(after.nucleus_id), "source_cluster_ids": ";".join(sorted(spec["clusters"])), "scope_rule": "exact original_frozen program or cluster; no new QC or annotation"})
    selected["lineage"] = selected.endpoint.map({group: value["lineage"] for group, value in ENDPOINTS.items()})
    selected["population_key"] = selected.endpoint.map({group: value["key"] for group, value in ENDPOINTS.items()})
    return selected, pd.DataFrame(rows), candidates


def validate_frozen_memberships(config, records, annotations, summary):
    index = read_table(records["population_index"]["path"])
    require({"population_key", "tool", "branch", "analysis_kind", "analysis_label", "n_nuclei", "membership_sha256", "counts_file", "counts_sha256"} <= set(index), "Population index schema differs")
    desired = {spec["key"] for spec in ENDPOINTS.values()}
    parts = []
    for chunk in table_chunks(records["population_memberships"]["path"], ["population_key", "nucleus_id"]):
        subset = chunk.loc[chunk.population_key.isin(desired)]
        if not subset.empty:
            parts.append(subset)
    require(parts, "No exact native population memberships")
    members = pd.concat(parts, ignore_index=True)
    native_keys = {}
    for group, spec in ENDPOINTS.items():
        population = index.loc[index.population_key.eq(spec["key"])]
        require(len(population) == 1, f"Population index must have one endpoint row: {group}")
        row = population.iloc[0]
        require(row.tool == "seurat" and row.branch == "original_full" and row.analysis_label == spec["label"], f"Wrong native endpoint mapping: {group}")
        ids = members.loc[members.population_key.eq(spec["key"]), "nucleus_id"]
        chosen = annotations.loc[annotations.endpoint.eq(group), "nucleus_id"]
        require(set(ids) == set(chosen), f"Raw source endpoint roster and annotation-final roster differ: {group}")
        require(len(ids) == spec["n"] == int(row.n_nuclei), f"Native endpoint count differs: {group}")
        require(membership_sha(ids) == row.membership_sha256, f"Native population membership hash differs: {group}")
        counts = checked_path(row.counts_file)
        bound = records["endpoint_counts_" + group]
        require(counts == Path(bound["path"]), f"Native aggregate path differs from explicit workflow input: {group}")
        require(row.counts_sha256 == bound["expected_sha256"], f"Native aggregate hash differs from explicit workflow input: {group}")
        native_keys[group] = counts
        records["native_pseudobulk:" + group] = {"path": str(counts), "expected_sha256": row.counts_sha256, "expected_bytes": None}
        summary.loc[summary.endpoint.eq(group), "native_population_membership_sha256"] = row.membership_sha256
        summary.loc[summary.endpoint.eq(group), "native_integer_pseudobulk_path"] = str(counts)
        summary.loc[summary.endpoint.eq(group), "native_integer_pseudobulk_sha256"] = row.counts_sha256
    # Independent category-family memberships validate exact CA1 program labels,
    # while the inhibitory broad family verifies the selected nuclei's parent.
    columns = ["nucleus_id", "sample", "condition", "tool", "branch", "population", "category_family", "label"]
    wanted = set(annotations.nucleus_id)
    category_parts = []
    for chunk in table_chunks(records["category_memberships"]["path"], columns):
        common = chunk.tool.eq("seurat") & chunk.branch.eq("original_full") & chunk.nucleus_id.isin(wanted)
        family = ((chunk.population == "CA1_ProS") & (chunk.category_family == "CA1_LAYER_EXPRESSION_PROGRAM")) | ((chunk.population == "Inhibitory") & (chunk.category_family == "REVIEWED_COARSE") & (chunk.label == "Inhibitory"))
        subset = chunk.loc[common & family]
        if not subset.empty:
            category_parts.append(subset)
    require(category_parts, "No matching frozen category memberships")
    categories = pd.concat(category_parts, ignore_index=True)
    require(not categories.nucleus_id.duplicated().any(), "Duplicate scoped frozen category memberships")
    require(set(categories.nucleus_id) == wanted, "Category-family coverage differs from the requested roster")
    aligned = categories.set_index("nucleus_id").loc[annotations.nucleus_id]
    for field in ("sample", "condition", "population"):
        require(np.array_equal(aligned[field].to_numpy(), annotations[field].to_numpy()), f"Frozen category metadata differs: {field}")
    for group in ("CA1_Deep", "CA1_Superficial"):
        ids = set(categories.loc[categories.label.eq(ENDPOINTS[group]["label"]), "nucleus_id"])
        require(ids == set(annotations.loc[annotations.endpoint.eq(group), "nucleus_id"]), f"CA1 category label differs: {group}")
    # Optional additional atomic declarations are checked only when both supplied.
    optional = {"population_atoms", "atomic_population_memberships"} & set(records)
    require(not optional or len(optional) == 2, "Supply population_atoms and atomic_population_memberships together")
    if optional:
        atoms = read_table(records["population_atoms"]["path"])
        needed_atoms = set(atoms.loc[atoms.population_key.isin(desired), "atomic_population_key"])
        atomic_parts = []
        for chunk in table_chunks(records["atomic_population_memberships"]["path"], ["population_key", "nucleus_id", "tool", "branch"]):
            subset = chunk.loc[chunk.tool.eq("seurat") & chunk.branch.eq("original_full") & chunk.population_key.isin(needed_atoms)]
            if not subset.empty:
                atomic_parts.append(subset)
        require(atomic_parts, "No optional atomic memberships")
        atomic = pd.concat(atomic_parts, ignore_index=True)
        for group, spec in ENDPOINTS.items():
            declared = set(atoms.loc[atoms.population_key.eq(spec["key"]), "atomic_population_key"])
            ids = atomic.loc[atomic.population_key.isin(declared), "nucleus_id"]
            require(not ids.duplicated().any(), f"Declared atomic masks overlap: {group}")
            require(set(ids) == set(annotations.loc[annotations.endpoint.eq(group), "nucleus_id"]), f"Atomic-union roster differs: {group}")
    return native_keys


def feature_mapping(gene_ids, raw_symbols):
    """Sum stable-ID features sharing a real symbol; blank symbols never merge."""
    ids, symbols = strings(gene_ids), strings(raw_symbols)
    require(len(ids) == len(symbols) and len(set(ids)) == len(ids), "Stable gene IDs must be unique and match the symbol axis")
    symbols = np.asarray([value.strip() for value in symbols], dtype=object)
    missing = np.asarray([value in ("", "NA", "nan", "None") for value in symbols])
    keys = ["__gene_id__" + gene_id if absent else symbol for gene_id, symbol, absent in zip(ids, symbols, missing)]
    require(not any(symbol.startswith("__gene_id__") for symbol, absent in zip(symbols, missing) if not absent), "Real symbol collides with reserved missing-symbol namespace")
    unique = list(dict.fromkeys(keys))
    positions = {key: index for index, key in enumerate(unique)}
    mapping = pd.DataFrame({"source_feature_index": np.arange(len(ids)), "gene_id": ids, "source_gene_symbol": symbols, "analysis_feature_id": keys, "symbol_missing": missing})
    mapping["output_feature_index"] = mapping.analysis_feature_id.map(positions).astype(int)
    sizes = mapping.groupby("analysis_feature_id", sort=False).size()
    mapping["n_source_features"] = mapping.analysis_feature_id.map(sizes).astype(int)
    mapping["duplicate_symbol"] = mapping.n_source_features > 1
    groups = mapping.groupby("analysis_feature_id", sort=False)
    features = pd.DataFrame([{"feature_id": key, "gene_symbol": key, "gene_ids": ";".join(groups.get_group(key).gene_id), "n_source_features": len(groups.get_group(key)), "duplicate_symbol": len(groups.get_group(key)) > 1, "symbol_missing": bool(groups.get_group(key).symbol_missing.any())} for key in unique])
    projection = sparse.csr_matrix((np.ones(len(ids), dtype=np.int64), (np.arange(len(ids)), mapping.output_feature_index.to_numpy())), shape=(len(ids), len(unique)))
    return mapping, features, projection


def read_selected_raw_counts(config, records, annotations):
    with h5py.File(records["counts"]["path"], "r") as handle:
        ids = strings(read_elem(handle["obs/nucleus_id"]))
        require(len(ids) == len(set(ids)), "Duplicate native raw nucleus IDs")
        wanted = set(annotations.nucleus_id)
        keep = np.asarray([nucleus in wanted for nucleus in ids])
        require(int(keep.sum()) == len(annotations), "At least one frozen nucleus is absent from raw RNA counts")
        obs = annotations.set_index("nucleus_id", drop=False).loc[ids[keep]].copy()
        for field in ("sample", "condition"):
            values = strings(read_elem(handle["obs/" + field]))[keep]
            require(np.array_equal(values, obs[field].to_numpy()), f"Exact raw/annotation {field} differs")
        gene_ids = strings(read_elem(handle["var/gene_id"]))
        symbols = strings(read_elem(handle["var/gene_symbol"]))
        matrix = handle["X"]
        encoding = matrix.attrs.get("encoding-type", "")
        if isinstance(encoding, bytes):
            encoding = encoding.decode()
        require(encoding == "csr_matrix", "Expected native sparse CSR RNA X")
        shape = tuple(int(value) for value in matrix.attrs["shape"])
        expected_shape = tuple(config.get("expected_source_shape", [108062, 32285]))
        require(shape == expected_shape == (len(ids), len(gene_ids)), "Raw count shape differs from the approved source")
        require(np.issubdtype(matrix["data"].dtype, np.integer), "Use raw integer RNA counts, not SCT/Harmony/log-expression")
        source_dtype = str(matrix["data"].dtype)
        chunks, libraries, source_indices = [], [], []
        checked_values = 0
        block_size = int(config.get("row_block_size", 2048))
        require(64 <= block_size <= 8192, "row_block_size must be 64–8192")
        for start in range(0, shape[0], block_size):
            stop = min(start + block_size, shape[0])
            if not keep[start:stop].any():
                continue
            pointers = matrix["indptr"][start:stop + 1].astype(np.int64)
            low, high = int(pointers[0]), int(pointers[-1])
            data = matrix["data"][low:high]
            require((data >= 0).all(), f"Negative source integer counts in raw block {start}")
            checked_values += len(data)
            block = sparse.csr_matrix((data.astype(np.int64), matrix["indices"][low:high], pointers - low), shape=(stop - start, shape[1]))[keep[start:stop]]
            library = np.asarray(block.sum(axis=1)).ravel().astype(np.int64)
            require((library > 0).all(), "Frozen kept nucleus has zero complete RNA library")
            chunks.append(block)
            libraries.append(library)
            source_indices.extend(np.flatnonzero(keep[start:stop]) + start)
        require(chunks, "No exact frozen raw count blocks")
        counts = sparse.vstack(chunks, format="csr")
        full_libraries = np.concatenate(libraries)
        require(counts.shape == (len(obs), len(gene_ids)), "Selected count axes changed")
        if "nCount_RNA" in handle["obs"]:
            native_library = np.asarray(read_elem(handle["obs/nCount_RNA"]), dtype=np.float64)[keep]
            require(np.array_equal(native_library, full_libraries.astype(np.float64)), "Full original gene denominator differs from native nCount_RNA")
        obs["source_row_index"] = np.asarray(source_indices, dtype=np.int64)
        obs["library_umi_full_genes"] = full_libraries
        obs["log1p_cp10k_scale_factor"] = 10000.0 / full_libraries
    return counts, obs, gene_ids, symbols, {"source_shape": list(shape), "source_count_dtype": source_dtype, "selected_integer_dtype": str(counts.dtype), "source_sparse_values_checked": int(checked_values), "selected_UMIs": int(full_libraries.sum()), "selected_nuclei": len(obs)}


def verify_pseudobulks(counts, obs, gene_ids, native_paths, samples):
    """Every raw stable-ID gene × mouse entry must reproduce the native aggregate."""
    rows = []
    for group, path in native_paths.items():
        source = read_table(path)
        require(source.columns.tolist() == ["gene_id", "gene_symbol"] + samples, f"Raw aggregate sample axis differs: {group}")
        require(np.array_equal(source.gene_id.to_numpy(), gene_ids), f"Raw aggregate gene-ID axis differs: {group}")
        original = source[samples].to_numpy(dtype=np.int64)
        observed = np.column_stack([np.asarray(counts[(obs.endpoint == group).to_numpy() & (obs['sample'] == sample).to_numpy()].sum(axis=0)).ravel() for sample in samples]).astype(np.int64)
        require(np.array_equal(original, observed), f"Selected raw full-gene counts do not reproduce original all-mouse aggregate: {group}")
        rows.append({"endpoint": group, "population_key": ENDPOINTS[group]["key"], "gene_id_rows": len(source), "mouse_count_columns": len(samples), "compared_entries": int(original.size), "differing_entries": 0, "native_count_path": str(path), "native_count_sha256": sha256(path), "status": "PASS"})
    return pd.DataFrame(rows)


def attach_sample_metadata(obs, sample_table):
    require({"sample", "condition", "sequencing_partition"} <= set(sample_table), "Sample roster schema differs")
    require(not sample_table['sample'].duplicated().any(), "Duplicated mice in sample roster")
    require(sample_table.groupby("condition")['sample'].nunique().to_dict() == {"WT": 6, "cKO": 6}, "Expected six independent mice per genotype")
    require(set(obs['sample']) == set(sample_table['sample']), "Selected nuclei and 12-mouse roster differ")
    lookup = sample_table.set_index("sample")
    require((obs['sample'].map(lookup.condition) == obs.condition).all(), "Genotype differs from the mouse roster")
    obs["sample_id"] = obs['sample']
    obs["barcode"] = obs.nucleus_id
    obs["sequencing_partition"] = obs['sample'].map(lookup.sequencing_partition)
    obs["lane"] = obs.sequencing_partition
    return obs


def normalized_copy(counts, full_libraries):
    normalized = counts.astype(np.float64).multiply((10000.0 / np.asarray(full_libraries))[:, None]).tocsr()
    np.log1p(normalized.data, out=normalized.data)
    require(np.isfinite(normalized.data).all() and (normalized.data >= 0).all(), "Invalid full-library log1p(CP10k) expression")
    return normalized.astype(np.float32)


def write_lineage(outdir, lineage, counts, log_expression, obs, features, normalization):
    folder = outdir / "lineages" / lineage
    folder.mkdir(parents=True, exist_ok=False)
    take = obs.lineage.eq(lineage).to_numpy()
    require(int(take.sum()) == EXPECTED_LINEAGES[lineage], f"Exact {lineage} lineage count differs")
    cells = obs.loc[take].copy()
    front = ["barcode", "nucleus_id", "sample_id", "sample", "condition", "endpoint", "lineage", "sequencing_partition", "lane"]
    table_write(cells[front + [col for col in cells if col not in front]], folder / "cells.tsv")
    table_write(features, folder / "genes.tsv")
    raw = counts[take].tocsr()
    with gzip.open(folder / "counts.mtx.gz", "wb") as handle:
        mmwrite(handle, raw.T.tocoo(), field="integer", comment="Raw integer RNA UMIs; full duplicate-symbol-collapsed gene universe; rows=genes, columns=nuclei; exact order genes.tsv/cells.tsv")
    if log_expression is not None:
        var = features.set_index("feature_id", drop=False).copy()
        value = ad.AnnData(X=raw, obs=cells.copy(), var=var)
        value.layers["log1p_cp10k"] = log_expression[take].tocsr()
        value.uns["normalization"] = normalization
        value.uns["lineage"] = lineage
        value.uns["genotype_replication"] = "Six independent WT and six cKO mice; nuclei are not independent genotype replicates"
        value.write_h5ad(folder / "expression.h5ad", compression="gzip")
    return {"lineage": lineage, "n_nuclei": int(raw.shape[0]), "n_genes": int(raw.shape[1]), "raw_X": "integer RNA counts after sum by real gene symbol", "normalized_layer": "log1p_cp10k" if log_expression is not None else None, "normalized_h5ad_written": log_expression is not None, "MTX_axes": "gene rows × nucleus columns", "membership_sha256": membership_sha(cells.nucleus_id)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    config_path = checked_path(args.config)
    config = json.loads(config_path.read_text())
    outdir, run_root = args.outdir.resolve(), Path(config["run_root"]).resolve()
    require(args.outdir.is_absolute() and Path(config["run_root"]).is_absolute(), "Use absolute run/output directories")
    require(outdir != run_root and outdir.is_relative_to(run_root), "Preparation output must be below the new run root")
    require(not outdir.exists() or (outdir.is_dir() and not any(outdir.iterdir())), "Refusing to overwrite nonempty output directory")
    records = input_records(config, config_path)
    require(not any(outdir == Path(item["path"]) or outdir in Path(item["path"]).parents for item in records.values()), "Output directory contains a source input")
    verify_records(records, "before")
    started = datetime.now(timezone.utc).isoformat()
    annotations, membership_audit, candidate_audit = select_annotations(records["nucleus_annotations"]["path"])
    native_paths = validate_frozen_memberships(config, records, annotations, membership_audit)
    verify_records(records, "before")  # Also bind the population-index-authorized count files.
    counts_gene_id, obs, gene_ids, raw_symbols, count_stats = read_selected_raw_counts(config, records, annotations)
    sample_table = read_table(records["samples"]["path"])
    obs = attach_sample_metadata(obs, sample_table)
    source_count_audit = verify_pseudobulks(counts_gene_id, obs, gene_ids, native_paths, sample_table['sample'].tolist())
    mapping, features, projection = feature_mapping(gene_ids, raw_symbols)
    counts = (counts_gene_id @ projection).tocsr()
    require(np.issubdtype(counts.dtype, np.integer), "Gene-symbol sum lost integer RNA counts")
    libraries = obs.library_umi_full_genes.to_numpy(dtype=np.int64)
    require(np.array_equal(np.asarray(counts.sum(axis=1)).ravel(), libraries), "Duplicate-symbol collapse changed original per-nucleus library")
    write_normalized_h5ad = config.get("write_normalized_h5ad", False)
    write_raw_gene_id_h5ad = config.get("write_raw_gene_id_h5ad", False)
    require(isinstance(write_normalized_h5ad, bool) and isinstance(write_raw_gene_id_h5ad, bool), "Optional H5AD output flags must be JSON booleans")
    normalization = {"method_if_materialized": "log1p(CP10000)", "target_sum": 10000, "normalized_layer_materialized": write_normalized_h5ad, "library_denominator": "all original 32,285 stable-ID RNA features before any downstream gene filtering", "duplicate_symbol_policy": "sum raw integer counts before normalization/log1p; source IDs retained in feature_mapping.tsv", "missing_symbol_policy": "separate reserved __gene_id__<stableID> feature keys; no merging of missing symbols", "full_gene_background_retained": True, "DEG_HVG_LR_gene_filter": "none", "imputation": "none", "integration": "none", "new_QC_or_nucleus_filter": "none; exact existing native keep flag only", "normalized_dtype_if_materialized": "float32", "source_counts_sha256": records["counts"]["before_sha256"]}
    log_expression = normalized_copy(counts, libraries) if write_normalized_h5ad else None
    for field in obs:
        if pd.api.types.is_object_dtype(obs[field].dtype):
            obs[field] = obs[field].fillna("").astype(str)
    require(obs.index.is_unique, "Output nucleus axis contains duplicates")
    outdir.mkdir(parents=True, exist_ok=True)
    table_write(mapping, outdir / "feature_mapping.tsv")
    table_write(mapping.loc[mapping.duplicate_symbol].copy(), outdir / "duplicate_symbol_audit.tsv")
    table_write(membership_audit, outdir / "frozen_membership_audit.tsv")
    table_write(candidate_audit, outdir / "frozen_candidates_before_after_keep.tsv.gz")
    table_write(source_count_audit, outdir / "raw_count_source_verification.tsv")
    table_write(obs.reset_index(drop=True), outdir / "nucleus_metadata.tsv.gz")
    coverage = obs.groupby(["endpoint", "lineage", "sample", "condition", "sequencing_partition"], sort=False).agg(n_nuclei=("nucleus_id", "size"), total_UMIs=("library_umi_full_genes", "sum")).reset_index()
    require(len(coverage) == 60, "Each of five endpoints must have all twelve mouse entries")
    table_write(coverage, outdir / "mouse_coverage.tsv")
    # Preserve a sparse, exact stable-ID subset as an audit/input alternative.
    if write_raw_gene_id_h5ad:
        raw_var = pd.DataFrame({"gene_id": gene_ids, "gene_symbol": raw_symbols}, index=pd.Index(gene_ids, name="source_gene_id"))
        raw_value = ad.AnnData(X=counts_gene_id, obs=obs.copy(), var=raw_var)
        raw_value.uns["counts"] = "Untransformed original integer RNA counts; all 32,285 stable-ID features"
        raw_value.write_h5ad(outdir / "raw_gene_id_counts.h5ad", compression="gzip")
    lineage_receipts = [write_lineage(outdir, lineage, counts, log_expression, obs, features, normalization) for lineage in EXPECTED_LINEAGES]
    json_write(outdir / "expression_normalization.json", normalization)
    verify_records(records, "after")
    json_write(outdir / "source_hashes.json", records)
    files = {str(path.relative_to(outdir)): {"bytes": path.stat().st_size, "sha256": sha256(path)} for path in sorted(outdir.rglob("*")) if path.is_file()}
    receipt = {"schema": SCHEMA, "status": "COMPLETE", "started_utc": started, "finished_utc": datetime.now(timezone.utc).isoformat(), "n_nuclei": len(obs), "n_raw_stable_ID_features": len(gene_ids), "n_collapsed_features": len(features), "group_counts": {group: int(obs.endpoint.eq(group).sum()) for group in ENDPOINTS}, "lineages": lineage_receipts, "source_count_stats": count_stats, "source_all_gene_pseudobulk_matches": True, "source_files_unchanged": True, "model_fits": 0, "new_degrees_of_freedom_from_nuclei": 0, "count_correction": {"IN_C8": {"incorrect_draft_count": 568, "authoritative_count": 468, "source": "hash-verified native population index and exact frozen memberships", "artificial_nuclei_added": 0}}, "caveats": ["CA1 Deep/Superficial are molecular expression programs, not measured physical layers", "IN_C8 retains its native mixed-identity review flag", "TF inference must use full genes; LR selection remains separately DEG-only", "Sequencing partition is partially associated with genotype and cKO_933 partition is unverified"], "files": files}
    json_write(outdir / "COMPLETE.json", receipt)


if __name__ == "__main__":
    main()
