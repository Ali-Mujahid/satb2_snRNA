#!/usr/bin/env python3
"""Independent Python integer count aggregation; no normalization or model fitting.

Reads the declared branch's counts once in backed CSR row chunks. The shared R
DESeq2/edgeR stage owns biological inference; this export checks transformation
consistency for identical gene IDs, nuclei memberships and biological mice.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path
import re
import time

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


MICE = ["wt_1122", "wt_1123", "wt_1178", "wt_1220", "wt_1226", "wt_1228",
        "cKO_1138", "cKO_1139", "cKO_1160", "cKO_1165", "cKO_933", "cKO_991"]
ATLAS = [("ATLAS_TOTAL_ALL", "ALL", None), ("ATLAS_TOTAL_EX", "EX", "REF_EX_HS"),
         ("ATLAS_TOTAL_IN", "IN", "REF_IN_HS"),
         ("ATLAS_TOTAL_CA1_PROS", "CA1_PROS", "REF_CA1_PROS_HS")]
SCANPY = [("allcells", "SCANPY_ALLCELLS_harmony", None),
          ("Excitatory", "SCANPY_EXCITATORY_harmony", "REF_EX_HS"),
          ("Inhibitory", "SCANPY_INHIBITORY_harmony", "REF_IN_HS"),
          ("CA1_ProS", "SCANPY_CA1_PROS_harmony", "REF_CA1_PROS_HS")]


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def id_digest(ids):
    """Order-independent population identity, sorted IDs separated by LF."""
    return hashlib.sha256("".join(value + "\n" for value in sorted(map(str, ids))).encode()).hexdigest()


def write_json(path, value):
    temp = path.with_name(path.name + ".partial")
    temp.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    temp.replace(path)


def write_tsv(path, frame):
    temp = path.with_name(path.name + ".partial")
    if path.suffix == ".gz":
        # Deterministic gzip metadata avoids time-dependent output identities.
        with temp.open("wb") as raw:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as gz:
                frame.to_csv(gz, sep="\t", index=False, na_rep="NA", encoding="utf-8")
    else:
        frame.to_csv(temp, sep="\t", index=False, na_rep="NA")
    temp.replace(path)


def truth(values, name):
    values = values.astype(str).str.strip().str.lower()
    if not values.isin(["true", "false", "1", "0"]).all():
        raise ValueError(f"{name}: expected explicit TRUE/FALSE or 1/0")
    return values.isin(["true", "1"]).to_numpy()


def mask_for(obs, flag):
    if flag is None:
        return np.ones(len(obs), dtype=bool)
    if flag not in obs:
        raise ValueError(f"Authoritative atlas membership flag missing: {flag}")
    return truth(obs[flag], flag)


def make_populations(obs, scanpy_root, count_branch, label_branch):
    ids = pd.Index(obs.nucleus_id)
    source_ids = set(ids)
    populations = []
    cluster_sources = []
    frozen_all_ids = None
    for family, label, flag in ATLAS:
        positions = np.flatnonzero(mask_for(obs, flag))
        populations.append({"family": family, "label": label, "positions": positions})
    for population, family, flag in SCANPY:
        metadata_path = scanpy_root / population / "nucleus_metadata.tsv"
        receipt_path = scanpy_root / population / "COMPLETE.json"
        if not receipt_path.is_file():
            raise ValueError(f"Missing completed Scanpy population receipt: {receipt_path}")
        receipt = json.loads(receipt_path.read_text())
        if receipt.get("status") != "COMPLETE" or receipt.get("population") != population or receipt.get("count_branch") != label_branch:
            raise ValueError(f"Wrong/incomplete Scanpy result identity: {receipt_path}")
        md = pd.read_csv(metadata_path, sep="\t", dtype={"nucleus_id": str, "sample": str,
                         "condition": str, "cluster_harmony": str}, low_memory=False)
        if not {"nucleus_id", "sample", "barcode", "condition", "cluster_harmony"}.issubset(md):
            raise ValueError(f"Missing Scanpy pseudobulk identity columns: {metadata_path}")
        if md.nucleus_id.isna().any() or md.nucleus_id.duplicated().any():
            raise ValueError(f"Invalid Scanpy nucleus IDs: {metadata_path}")
        if not set(md["sample"]).issubset(MICE):
            raise ValueError("Frozen Scanpy labels contain unknown source capture IDs")
        if not (md.nucleus_id == md["sample"] + "_" + md.barcode).all():
            raise ValueError("Frozen Scanpy labels violate exact sample/barcode nucleus identity")
        expected_condition = md["sample"].map({s: "WT" if s.startswith("wt_") else "cKO" for s in MICE})
        if not md.condition.eq(expected_condition).all():
            raise ValueError("Frozen Scanpy labels conflict with the twelve-mouse genotype roster")
        if len(md) != int(receipt["n_nuclei"]):
            raise ValueError("Frozen Scanpy metadata size differs from its completed population receipt")
        label_ids = set(md.nucleus_id)
        if population == "allcells":
            frozen_all_ids = label_ids
            if not source_ids.issubset(frozen_all_ids):
                raise ValueError("Count branch contains nuclei absent from the frozen original-full labels")
        elif not label_ids.issubset(frozen_all_ids):
            raise ValueError("Frozen subset labels include nuclei absent from the frozen allcells capture cohort")
        expected_ids = set(ids[mask_for(obs, flag)])
        if count_branch == label_branch:
            if label_ids != expected_ids:
                raise ValueError(f"Primary Scanpy {population} membership differs from the exact count cohort")
        elif (label_ids & source_ids) != expected_ids:
            raise ValueError(f"Frozen Scanpy {population} does not cover exactly the paired atlas cohort")
        if md.cluster_harmony.isna().any() or md.cluster_harmony.eq("").any():
            raise ValueError("Missing frozen cluster labels must not become a population")
        frozen_cluster_labels = sorted(md.cluster_harmony.unique())
        # Freeze the original-full clustering, then intersect with paired source
        # nuclei. No correction-branch reclustering enters count sensitivity DE.
        original_label_nuclei = len(md)
        md = md.loc[md.nucleus_id.isin(source_ids)].copy()
        positions = ids.get_indexer(md.nucleus_id)
        if (positions < 0).any():
            raise ValueError("Scanpy nucleus absent from source-count cohort")
        aligned = obs.iloc[positions]
        for field in ("sample", "condition"):
            if not np.array_equal(md[field].to_numpy(), aligned[field].astype(str).to_numpy()):
                raise ValueError(f"Scanpy {field} conflicts with source-count metadata")
        for label in frozen_cluster_labels:
            positions = ids.get_indexer(md.loc[md.cluster_harmony.eq(label), "nucleus_id"])
            populations.append({"family": family, "label": str(label), "positions": positions})
        cluster_sources.append({"population": population, "metadata_path": str(metadata_path),
                                "metadata_sha256": sha256(metadata_path),
                                "receipt_sha256": sha256(receipt_path),
                                "source_label_branch": label_branch,
                                "original_label_nuclei": original_label_nuclei,
                                "count_branch_intersection_nuclei": len(md)})
    for i, population in enumerate(populations, 1):
        population["population_key"] = f"P{i:04d}_" + re.sub(r"[^A-Za-z0-9_-]", "_", population["family"])
        positions = population["positions"]
        if len(np.unique(positions)) != len(positions):
            raise ValueError("Duplicate nucleus within a pseudobulk membership")
        population["membership_sha256"] = id_digest(ids[positions])
    return populations, cluster_sources


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--shared-dir", required=True, type=Path)
    parser.add_argument("--scanpy-root", required=True, type=Path,
                        help="Branch directory containing allcells/Excitatory/Inhibitory/CA1_ProS")
    parser.add_argument("--branch", required=True,
                        choices=["original_full", "paired_original", "paired_cellbender"])
    parser.add_argument("--label-branch", default="original_full", choices=["original_full"],
                        help="Freeze original-full Scanpy memberships, then intersect paired count-source IDs")
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--chunk-size", type=int, default=4096)
    args = parser.parse_args()
    if args.chunk_size < 1:
        parser.error("chunk-size must be positive")
    started = time.monotonic()
    out = args.outdir.resolve()
    if out.exists() and (not out.is_dir() or any(out.iterdir())):
        raise ValueError("Output destination contains data; use workflow caching or a new output directory")
    # Snakemake creates output parent directories before launching this stage.
    out.mkdir(parents=True, exist_ok=True)
    shared_path = args.shared_dir / "branch_counts.h5ad"
    shared_receipt_path = args.shared_dir / "COMPLETE.json"
    receipt = json.loads(shared_receipt_path.read_text())
    if receipt.get("status") != "COMPLETE" or receipt.get("count_branch") != args.branch:
        raise ValueError("Shared source completion/branch identity mismatch")
    source_hash = sha256(shared_path)
    if source_hash != receipt["branch_counts_h5ad_sha256"]:
        raise ValueError("Shared count artifact no longer matches the prepared count-source receipt")
    data = ad.read_h5ad(shared_path, backed="r")
    try:
        if not hasattr(data.X, "format") or data.X.format != "csr":
            raise ValueError("Independent backed row-pass requires the shared CSR count artifact")
        obs = data.obs.copy()
        if "nucleus_id" not in obs:
            obs["nucleus_id"] = data.obs_names.astype(str)
        if not np.array_equal(obs.nucleus_id.astype(str), data.obs_names.astype(str)):
            raise ValueError("AnnData nucleus_id column/index differ")
        if obs.nucleus_id.duplicated().any() or set(obs["sample"].astype(str)) != set(MICE):
            raise ValueError("Counts must contain unique nuclei from exactly the twelve mice")
        if not (obs.nucleus_id.astype(str) == obs["sample"].astype(str) + "_" + obs.barcode.astype(str)).all():
            raise ValueError("Preserve exact sample + '_' + barcode identity")
        mouse = obs[[c for c in ("sample", "condition", "sequencing_partition") if c in obs]].drop_duplicates()
        if mouse["sample"].duplicated().any():
            raise ValueError("Mouse-level metadata conflict within the source cohort")
        mouse = mouse.set_index("sample").loc[MICE].reset_index()
        expected_conditions = ["WT"] * 6 + ["cKO"] * 6
        if mouse.condition.astype(str).tolist() != expected_conditions:
            raise ValueError("Six WT and six cKO mice required, with unchanged genotype identity")
        genes = pd.DataFrame({"gene_id": data.var_names.astype(str),
                             "gene_symbol": data.var.gene_symbol.astype(str).to_numpy()})
        if genes.gene_id.duplicated().any() or genes.gene_id.isna().any():
            raise ValueError("Gene IDs must be unique stable rows")
        genes["gene_symbol"] = genes.gene_symbol.mask(genes.gene_symbol.isin(["", "nan", "None"]), genes.gene_id)
        populations, cluster_sources = make_populations(obs, args.scanpy_root, args.branch, args.label_branch)
        samples = pd.Index(MICE).get_indexer(obs["sample"].astype(str))
        n_groups = len(populations) * len(MICE)
        estimated_aggregate_bytes = data.n_vars * n_groups * np.dtype(np.int64).itemsize
        if estimated_aggregate_bytes > 8 * 1024**3:
            raise ValueError("Pseudobulk aggregate exceeds the explicit 8GiB dense accumulator bound")
        group_rows, group_cols = [], []
        for p, population in enumerate(populations):
            positions = population["positions"]
            group_rows.append(positions)
            group_cols.append(p * len(MICE) + samples[positions])
        rows = np.concatenate(group_rows)
        cols = np.concatenate(group_cols)
        assignment = sparse.csr_matrix((np.ones(len(rows), dtype=np.int64), (rows, cols)),
                                       shape=(data.n_obs, n_groups))
        assignment.sum_duplicates()
        if np.any(assignment.data != 1):
            raise ValueError("Repeated membership within a population/sample column")
        aggregate = np.zeros((data.n_vars, n_groups), dtype=np.int64)
        expected_group_umis = np.zeros(n_groups, dtype=np.int64)
        nucleus_counts = np.asarray(assignment.sum(axis=0)).ravel().astype(np.int64)
        source_gene_totals = np.zeros(data.n_vars, dtype=np.int64)
        source_mouse_umis = np.zeros(len(MICE), dtype=np.int64)
        data_digest, indices_digest, indptr_digest = (hashlib.sha256() for _ in range(3))
        offset = 0
        for start in range(0, data.n_obs, args.chunk_size):
            end = min(start + args.chunk_size, data.n_obs)
            x = data.X[start:end].tocsr()
            if np.any(~np.isfinite(x.data)) or np.any(x.data < 0) or np.any(x.data != np.floor(x.data)):
                raise ValueError("Pseudobulk input must contain finite nonnegative integer counts")
            if np.any(x.data > np.iinfo(np.int64).max // max(1, data.n_obs)):
                raise OverflowError("Count magnitude threatens exact int64 aggregation")
            x = x.astype(np.int64, copy=False)
            x.sum_duplicates()
            x.sort_indices()
            row_totals = np.asarray(x.sum(axis=1)).ravel().astype(np.int64)
            if "branch_UMIs_verified" in obs and not np.array_equal(row_totals, obs.branch_UMIs_verified.iloc[start:end].to_numpy()):
                raise ValueError("Per-nucleus source UMI totals differ from shared preparation")
            np.add.at(source_mouse_umis, samples[start:end], row_totals)
            source_gene_totals += np.asarray(x.sum(axis=0)).ravel().astype(np.int64)
            local_assignment = assignment[start:end]
            product = (x.T @ local_assignment).tocoo()
            product.sum_duplicates()
            aggregate[product.row, product.col] += product.data
            local = local_assignment.tocoo()
            np.add.at(expected_group_umis, local.col, row_totals[local.row])
            # Canonical CSR component digests are streamed during this sole count
            # pass; canonical int64 types make them independent of storage dtype.
            data_digest.update(memoryview(np.ascontiguousarray(x.data, dtype=np.int64)).cast("B"))
            indices_digest.update(memoryview(np.ascontiguousarray(x.indices, dtype=np.int64)).cast("B"))
            global_indptr = x.indptr.astype(np.int64) + offset
            indptr_digest.update(memoryview(np.ascontiguousarray(global_indptr if start == 0 else global_indptr[1:])).cast("B"))
            offset += x.nnz
            print(f"Aggregated {args.branch}: {end:,}/{data.n_obs:,} nuclei", flush=True)
        if not np.array_equal(aggregate.sum(axis=0), expected_group_umis):
            raise ValueError("Gene-by-mouse sums do not conserve membership library totals")
        all_matrix = aggregate[:, :len(MICE)]
        if not np.array_equal(all_matrix.sum(axis=1), source_gene_totals) or not np.array_equal(all_matrix.sum(axis=0), source_mouse_umis):
            raise ValueError("ALL population fails gene-level/mouse-level source conservation")
        if np.any(aggregate > np.iinfo(np.int32).max):
            raise OverflowError("Aggregate exceeds the R integer-count contract")
        canonical = {"shape": [data.n_obs, data.n_vars], "nnz": offset,
                     "data_int64_sha256": data_digest.hexdigest(),
                     "indices_int64_sha256": indices_digest.hexdigest(),
                     "indptr_int64_sha256": indptr_digest.hexdigest()}
        canonical_digest = hashlib.sha256(json.dumps(canonical, sort_keys=True).encode()).hexdigest()
        index_rows, membership_frames, coverage_frames, eligibility_rows, count_hashes = [], [], [], [], {}
        for p, population in enumerate(populations):
            key = population["population_key"]
            folder = out / key
            folder.mkdir(exist_ok=True)
            y = aggregate[:, p*len(MICE):(p+1)*len(MICE)]
            coverage = mouse.copy()
            coverage["n_nuclei"] = nucleus_counts[p*len(MICE):(p+1)*len(MICE)]
            coverage["total_UMIs"] = y.sum(axis=0)
            if np.any((coverage.n_nuclei.to_numpy() == 0) != (coverage.total_UMIs.to_numpy() == 0)):
                raise ValueError("Absent population and zero library disagree")
            coverage["population_key"] = key
            coverage["eligible_primary"] = (coverage.n_nuclei >= 20) & (coverage.total_UMIs > 0)
            coverage["inclusion_or_exclusion_reason"] = np.where(coverage.n_nuclei == 0, "ABSENT_POPULATION",
                np.where(coverage.n_nuclei < 20, "FEWER_THAN_20_NUCLEI", "INCLUDED"))
            nwt = int(((coverage.condition.astype(str) == "WT") & coverage.eligible_primary).sum())
            nko = int(((coverage.condition.astype(str) == "cKO") & coverage.eligible_primary).sum())
            eligible_model = min(nwt, nko) >= 3
            eligibility_rows.append({"population_key": key, "family": population["family"],
                "label": population["label"], "n_eligible_WT": nwt, "n_eligible_cKO": nko,
                "eligible_primary_unadjusted_model": eligible_model,
                "unadjusted_design_full_rank": nwt > 0 and nko > 0,
                "unadjusted_residual_df": nwt+nko-2 if nwt > 0 and nko > 0 else 0,
                "status": "ELIGIBLE_COUNTS_NO_MODEL_FIT" if eligible_model else "FEWER_THAN_THREE_ELIGIBLE_MICE_PER_GENOTYPE"})
            count_frame = pd.concat([genes.reset_index(drop=True), pd.DataFrame(y, columns=MICE)], axis=1)
            count_path = folder / "SOURCE_INTEGER_COUNTS_ALL_GENES_ALL12.tsv.gz"
            write_tsv(count_path, count_frame)
            write_tsv(folder / "ALL12_POPULATION_COVERAGE.tsv", coverage)
            count_hashes[key] = sha256(count_path)
            membership_frames.append(pd.DataFrame({"population_key": key,
                "source_label_branch": args.label_branch,
                "nucleus_id": obs.nucleus_id.iloc[population["positions"]].astype(str).to_numpy()}))
            index_rows.append({"population_key": key, "family": population["family"], "label": population["label"],
                "n_nuclei": len(population["positions"]), "membership_sha256": population["membership_sha256"],
                "label_interpretation": "CA1/ProS_candidate_not_anatomical_CA1" if "CA1_PROS" in population["family"] else "declared_analysis_identity",
                "count_source": args.branch, "source_label_branch": args.label_branch,
                "membership_hash_scope": "sorted intersected source nucleus IDs; source label branch recorded separately",
                "counts_file": str(count_path)})
            coverage_frames.append(coverage)
        index_frame = pd.DataFrame(index_rows)
        write_tsv(out / "population_index.tsv", index_frame)
        write_tsv(out / "population_memberships.tsv.gz", pd.concat(membership_frames, ignore_index=True))
        write_tsv(out / "all12_population_coverage.tsv", pd.concat(coverage_frames, ignore_index=True))
        write_tsv(out / "population_eligibility.tsv", pd.DataFrame(eligibility_rows))
        write_tsv(out / "confirmed_all12_mouse_roster.tsv", mouse)
        write_tsv(out / "source_count_conservation.tsv", pd.DataFrame({"sample": MICE,
            "source_n_nuclei": [int((samples == i).sum()) for i in range(len(MICE))],
            "source_total_UMIs": source_mouse_umis, "atlas_ALL_total_UMIs": all_matrix.sum(axis=0),
            "exact_conservation": True}))
        write_json(out / "COMPLETE.json", {"status": "COMPLETE", "branch": args.branch,
            "source_label_branch": args.label_branch,
            "inference_scope": "Independent integer aggregation only; no normalization, gene test, or model fit",
            "shared_count_artifact": str(shared_path), "shared_count_artifact_sha256": source_hash,
            "shared_receipt_sha256": sha256(shared_receipt_path),
            "shared_preparation_count_digest": receipt.get("branch_counts_digest"),
            "canonical_CSR_int64_components": canonical, "canonical_count_digest": canonical_digest,
            "all_source_membership_sha256": id_digest(obs.nucleus_id),
            "gene_ID_order_sha256": hashlib.sha256(("\n".join(genes.gene_id) + "\n").encode()).hexdigest(),
            "n_nuclei": data.n_obs, "n_genes": data.n_vars, "n_populations": len(populations),
            "source_total_UMIs": int(source_mouse_umis.sum()), "counts_conservation_verified": True,
            "all12_sample_order": MICE, "cluster_sources": cluster_sources,
            "population_counts_file_sha256": count_hashes,
            "aggregate_accumulator_bytes": estimated_aggregate_bytes,
            "chunk_size": args.chunk_size, "seconds": time.monotonic()-started,
            "R_comparison_join": "family + label + sorted nucleus-ID membership_sha256; ordinal keys may differ"})
    finally:
        data.file.close()


if __name__ == "__main__":
    main()
