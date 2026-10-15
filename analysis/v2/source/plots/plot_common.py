"""Shared rendering/summary utilities. This module never opens nucleus counts."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

KEYS = ["tool", "branch", "population", "cluster_scope", "cluster_id"]
POPULATION_ALIASES = {
    "All": "allcells", "ALL": "allcells", "global": "allcells",
    "EX": "Excitatory", "IN": "Inhibitory", "CA1": "CA1_ProS",
    "CA1_PROS": "CA1_ProS", "CA1_PROS_CANDIDATE": "CA1_ProS",
}
BRANCHES = ["original_full", "paired_original", "paired_cellbender"]
POPULATIONS = ["allcells", "Excitatory", "Inhibitory", "CA1_ProS"]


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_table(path):
    return pd.read_csv(path, sep="\t", dtype={"cluster_id": str,
                       "source_cluster_id": str, "gene_id": str, "sample": str,
                       "nucleus_id": str}, low_memory=False)


def canonical_keys(frame):
    frame = frame.copy()
    if "branch" not in frame and "count_branch" in frame:
        frame = frame.rename(columns={"count_branch": "branch"})
    if "cluster_scope" not in frame and "label_scope" in frame:
        frame = frame.rename(columns={"label_scope": "cluster_scope"})
    if "population" in frame:
        frame["population"] = frame.population.replace(POPULATION_ALIASES)
    if "cluster_id" in frame:
        frame["cluster_id"] = frame.cluster_id.astype(str)
    if "tool" in frame:
        frame["tool"] = frame.tool.str.lower()
    if "cluster_scope" in frame:
        frame["cluster_scope"] = frame.cluster_scope.replace({
            "frozen_original_harmony_r0p4": "original_frozen",
            "fresh_branch_harmony_r0p4": "fresh_branch"})
    return frame


def bool_values(values):
    if pd.api.types.is_bool_dtype(values):
        return values.fillna(False)
    mapping = {"true": True, "false": False, "1": True, "0": False}
    z = values.astype(str).str.lower().map(mapping)
    if z.isna().any():
        raise ValueError("Boolean column contains missing/unrecognized values")
    return z.astype(bool)


def empty_outdir(path):
    path = Path(path).resolve()
    if path.exists() and any(path.iterdir()):
        raise FileExistsError(f"Destination already contains artifacts: {path}")
    path.mkdir(parents=True, exist_ok=True)
    return path


def write_table(frame, path):
    frame.to_csv(path, sep="\t", index=False, na_rep="NA")


def write_receipt(out, config_path, inputs, extra):
    receipt = dict(status="COMPLETE", scientific_fitting_performed=False,
                   raw_nucleus_count_aggregation_performed=False,
                   config_path=str(Path(config_path).resolve()),
                   config_sha256=sha256(config_path),
                   input_sources=[{"path": str(Path(p).resolve()), "sha256": sha256(p)}
                                  for p in dict.fromkeys(map(str, inputs))],
                   output_sha256={str(p.relative_to(out)): sha256(p) for p in sorted(out.rglob("*"))
                                  if p.is_file() and p.name != "COMPLETE.json"},
                   artifacts=[{"path": str(p.resolve()), "sha256": sha256(p)}
                              for p in sorted(out.rglob("*"))
                              if p.is_file() and p.name != "COMPLETE.json"])
    receipt.update(extra)
    (out / "COMPLETE.json").write_text(json.dumps(receipt, indent=2) + "\n")


def numeric_checks(frame, columns, integer=False):
    a = frame[columns].to_numpy(dtype=float)
    if not np.isfinite(a).all() or (a < 0).any():
        raise ValueError(f"Invalid nonnegative native aggregate: {columns}")
    if integer and not np.equal(a, np.floor(a)).all():
        raise ValueError(f"Noninteger native aggregate: {columns}")


def cluster_sort(value):
    try:
        return (0, float(value))
    except (ValueError, TypeError):
        return (1, str(value))
