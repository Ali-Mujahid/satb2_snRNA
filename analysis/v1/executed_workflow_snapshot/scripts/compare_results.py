#!/usr/bin/env python3
"""Compare fresh results, without fitting, reclustering or modifying inputs."""
from __future__ import annotations

import argparse
import csv
import functools
import gzip
import hashlib
import itertools
import json
import math
import re
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from decimal import Decimal, InvalidOperation
from pathlib import Path

import numpy as np

def initialize_plotting():
    global plt, PdfPages
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

BRANCHES = ("original_full", "paired_original", "paired_cellbender")
POPULATIONS = ("allcells", "Excitatory", "Inhibitory", "CA1_ProS")
SUBSETS = {"Excitatory": "EX", "Inhibitory": "IN", "CA1_ProS": "CA1_PROS_CANDIDATE"}
MISSING = {"", "NA", "NaN", "nan", "None", "<NA>"}
MICE = ("wt_1122", "wt_1123", "wt_1178", "wt_1220", "wt_1226", "wt_1228",
        "cKO_1138", "cKO_1139", "cKO_1160", "cKO_1165", "cKO_933", "cKO_991")
CONDITIONS = {sample: "WT" if sample.startswith("wt_") else "cKO" for sample in MICE}


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def member_digest(ids):
    return hashlib.sha256(("\n".join(sorted(ids)) + "\n").encode()).hexdigest()


def read_tsv(path):
    with Path(path).open(newline="", encoding="utf-8-sig") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def write_tsv(path, records, fields):
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(records)


def numeric(value):
    try:
        result = float(value)
        return result if math.isfinite(result) else None
    except (TypeError, ValueError):
        return None


def label_sort(label):
    value = numeric(label)
    return (0, value, str(label)) if value is not None else (1, 0, str(label))


def partition_metrics(counts):
    """ARI and arithmetic-mean normalized mutual information from a contingency."""
    counts = np.asarray(counts, dtype=np.int64)
    n = int(counts.sum())
    if n == 0:
        raise ValueError("No common labeled nuclei")
    if n < 2:
        return {"ARI": 1.0, "NMI_arithmetic": 1.0}
    rows, columns = counts.sum(axis=1), counts.sum(axis=0)
    choose2 = lambda values: float(np.sum(values.astype(np.float64) * (values - 1) / 2))
    within, row_pairs, column_pairs = choose2(counts), choose2(rows), choose2(columns)
    expected = row_pairs * column_pairs / (n * (n - 1) / 2)
    denominator = (row_pairs + column_pairs) / 2 - expected
    ari = (within - expected) / denominator if abs(denominator) > 1e-12 else 1.0
    p, q = rows / n, columns / n
    h_left, h_right = -np.sum(p[p > 0] * np.log(p[p > 0])), -np.sum(q[q > 0] * np.log(q[q > 0]))
    ii, jj = np.nonzero(counts)
    joint = counts[ii, jj] / n
    mutual = np.sum(joint * np.log(joint / (p[ii] * q[jj])))
    entropy = h_left + h_right
    nmi = 2 * mutual / entropy if entropy > 1e-12 else 1.0
    return {"ARI": float(ari), "NMI_arithmetic": float(np.clip(nmi, 0, 1))}


def hungarian_maximize(weights):
    """Exact one-to-one maximum-overlap assignment; rectangular arrays are padded."""
    weights = np.asarray(weights, dtype=np.int64)
    rows, columns = weights.shape
    size = max(rows, columns)
    if size == 0:
        return []
    padded = np.zeros((size, size), dtype=np.int64)
    padded[:rows, :columns] = weights
    cost = padded.max() - padded
    u, v = np.zeros(size + 1), np.zeros(size + 1)
    p, way = np.zeros(size + 1, dtype=int), np.zeros(size + 1, dtype=int)
    for row in range(1, size + 1):
        p[0] = row
        column = 0
        minimum = np.full(size + 1, np.inf)
        used = np.zeros(size + 1, dtype=bool)
        while True:
            used[column] = True
            current_row = p[column]
            delta, next_column = np.inf, 0
            for candidate in range(1, size + 1):
                if used[candidate]:
                    continue
                reduced = cost[current_row - 1, candidate - 1] - u[current_row] - v[candidate]
                if reduced < minimum[candidate]:
                    minimum[candidate], way[candidate] = reduced, column
                if minimum[candidate] < delta:
                    delta, next_column = minimum[candidate], candidate
            for candidate in range(size + 1):
                if used[candidate]:
                    u[p[candidate]] += delta
                    v[candidate] -= delta
                else:
                    minimum[candidate] -= delta
            column = next_column
            if p[column] == 0:
                break
        while True:
            previous = way[column]
            p[column] = p[previous]
            column = previous
            if column == 0:
                break
    return [(p[column] - 1, column - 1) for column in range(1, size + 1)
            if 0 < p[column] <= rows and column <= columns and weights[p[column] - 1, column - 1] > 0]


@dataclass
class Dataset:
    tool: str
    branch: str
    population: str
    path: Path
    ids: list
    index: dict
    samples: np.ndarray
    conditions: np.ndarray
    columns: dict
    partitions: dict
    partition_columns: dict
    umaps: dict
    sources: list = field(default_factory=list)
    by_animal: dict = field(default_factory=dict)
    source_digests: dict = field(default_factory=dict)
    features: set = field(default_factory=set)
    pca_features: set = field(default_factory=set)

    @property
    def key(self):
        return (self.tool, self.branch, self.population)


def metadata_path(root, tool, branch, population):
    base = root / "results" / ("seurat" if tool == "seurat_sct" else tool) / branch
    if tool == "seurat_sct":
        return base / "sct_global" / "nucleus_metadata.tsv"
    if tool == "scanpy":
        return base / population / "nucleus_metadata.tsv"
    if population == "allcells":
        return base / "log_global" / "nucleus_metadata.tsv"
    return base / "annotate_subsets" / SUBSETS[population] / "nucleus_metadata.tsv"


def load_dataset(root, tool, branch, population, main_resolution):
    path = metadata_path(root, tool, branch, population)
    if not path.is_file():
        raise FileNotFoundError(path)
    if not path.resolve().is_relative_to((root / "results").resolve()):
        raise ValueError(f"Historical/external metadata input refused: {path}")
    digests = {str(path): sha256(path)}
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        for key in ("nucleus_id", "sample", "condition"):
            if key not in fields:
                raise ValueError(f"{path}: missing {key}")
        wanted = [x for x in fields if re.search(r"(?:leiden|log|sct|EX|IN|CA1_PROS_CANDIDATE)_(?:unintegrated|harmony)_r", x)
                  or x.startswith(("cluster_", "UMAP_", "reference_", "mmc_"))]
        columns = {x: [] for x in wanted}
        ids, samples, conditions = [], [], []
        for row in reader:
            identity = row["nucleus_id"]
            if not identity or row["sample"] in MISSING or row["condition"] in MISSING:
                raise ValueError(f"{path}: missing nucleus/sample/condition")
            if row.get("count_source") not in (None, "", branch):
                raise ValueError(f"{path}: count_source does not match declared branch {branch}")
            ids.append(identity)
            samples.append(sys.intern(row["sample"]))
            conditions.append(sys.intern(row["condition"]))
            for column in wanted:
                columns[column].append(sys.intern(row.get(column, "")))
    if not ids or len(set(ids)) != len(ids):
        raise ValueError(f"{path}: empty or duplicated nucleus IDs")
    if set(samples) != set(MICE) or any(CONDITIONS.get(sample) != condition for sample, condition in zip(samples, conditions)):
        raise ValueError(f"{path}: not the declared six WT / six cKO mouse roster")
    columns = {k: np.asarray(v, dtype=object) for k, v in columns.items()}
    partitions, partition_columns = {}, {}
    for column, values in columns.items():
        match = re.fullmatch(r"(?:leiden|log|sct|EX|IN|CA1_PROS_CANDIDATE)_(unintegrated|harmony)_r([0-9.p]+)", column)
        if match:
            key = (match.group(1), float(match.group(2).replace("p", ".")))
        elif column in ("cluster_unintegrated", "cluster_harmony"):
            key = (column.removeprefix("cluster_"), main_resolution)
        else:
            continue
        if any(x in MISSING for x in values):
            raise ValueError(f"{path}: missing cluster label in {column}")
        if key in partitions and not np.array_equal(values, partitions[key]):
            raise ValueError(f"{path}: inconsistent primary/grid aliases for {key}")
        partitions[key], partition_columns[key] = values, column
    for representation in ("unintegrated", "harmony"):
        if (representation, main_resolution) not in partitions:
            raise ValueError(f"{path}: missing primary {representation} resolution {main_resolution}")
    sources = [path]
    umaps = {}
    for representation, native_prefix in (("unintegrated", "UMAPU"), ("harmony", "UMAPH")):
        pair = [f"UMAP_{representation}_{i}" for i in (1, 2)]
        if all(x in columns for x in pair):
            coordinates = np.column_stack([columns[x].astype(float) for x in pair])
        else:
            embedding = path.parent / f"umap_{representation}.tsv"
            if not embedding.is_file():
                raise FileNotFoundError(f"UMAP coordinates missing: {embedding}")
            digests[str(embedding)] = sha256(embedding)
            points = {}
            for row in read_tsv(embedding):
                identity = row["nucleus_id"]
                if identity in points:
                    raise ValueError(f"{embedding}: duplicate nucleus ID")
                points[identity] = [float(row[f"{native_prefix}_{i}"]) for i in (1, 2)]
            if set(points) != set(ids):
                raise ValueError(f"{embedding}: UMAP/metadata nuclei differ")
            coordinates = np.asarray([points[x] for x in ids])
            sources.append(embedding)
        if coordinates.shape != (len(ids), 2) or not np.isfinite(coordinates).all():
            raise ValueError(f"{path}: invalid {representation} UMAP coordinates")
        umaps[representation] = coordinates
    sentinel = next((p / "COMPLETE.json" for p in (path.parent, path.parent.parent, path.parent.parent.parent)
                     if (p / "COMPLETE.json").is_file()), None)
    if sentinel is None:
        raise FileNotFoundError(f"Missing completed producer receipt: {path}")
    digests[str(sentinel)] = sha256(sentinel)
    receipt = json.loads(sentinel.read_text())
    if tool == "scanpy":
        if receipt.get("status") != "COMPLETE" or receipt.get("count_branch") != branch or receipt.get("population") != population:
            raise ValueError(f"{sentinel}: wrong/incomplete Scanpy producer")
    else:
        expected_stage = "sct_global" if tool == "seurat_sct" else ("log_global" if population == "allcells" else "annotate_subsets")
        if receipt.get("status") != "EXECUTED_REVIEW_REQUIRED" or receipt.get("count_source") != branch or receipt.get("stage") != expected_stage:
            raise ValueError(f"{sentinel}: wrong/incomplete Seurat producer")
    sources.append(sentinel)
    shared = root / "results" / "shared" / branch
    membership_path, shared_receipt = shared / f"membership_{population}.tsv", shared / "COMPLETE.json"
    if not membership_path.is_file() or not shared_receipt.is_file():
        raise FileNotFoundError(f"Missing authoritative completed shared membership: {membership_path}")
    digests[str(membership_path)], digests[str(shared_receipt)] = sha256(membership_path), sha256(shared_receipt)
    shared_status = json.loads(shared_receipt.read_text())
    if shared_status.get("status") != "COMPLETE" or shared_status.get("count_branch") != branch:
        raise ValueError(f"{shared_receipt}: incomplete/wrong branch source")
    membership_rows = list(read_tsv(membership_path))
    authoritative = {r["nucleus_id"]: (r["sample"], r["condition"]) for r in membership_rows}
    if len(authoritative) != len(membership_rows):
        raise ValueError(f"{membership_path}: duplicate authoritative nucleus ID")
    if set(authoritative) != set(ids) or any(authoritative[x] != (sample, condition) for x, sample, condition in zip(ids, samples, conditions)):
        raise ValueError(f"{path}: fresh metadata differs from authoritative shared membership")
    sources.extend([membership_path, shared_receipt])
    feature_path = shared / f"common_hvg_{population}.tsv" if tool == "scanpy" else path.parent / "fixed_features.tsv"
    if not feature_path.is_file():
        raise FileNotFoundError(f"Missing actual representation feature panel: {feature_path}")
    digests[str(feature_path)] = sha256(feature_path)
    feature_ids = [r["gene_id"] for r in read_tsv(feature_path)]
    if not feature_ids or len(set(feature_ids)) != len(feature_ids):
        raise ValueError(f"{feature_path}: empty or duplicate modeled representation features")
    sources.append(feature_path)
    pca_features = set(feature_ids)
    if tool != "scanpy":
        pca_path = path.parent / "PCA_ACTUAL_FEATURES.tsv"
        if not pca_path.is_file():
            raise FileNotFoundError(f"Missing actual PCA feature roster: {pca_path}")
        digests[str(pca_path)] = sha256(pca_path)
        pca_ids = [r["gene_id"] for r in read_tsv(pca_path)]
        if not pca_ids or len(set(pca_ids)) != len(pca_ids) or not set(pca_ids).issubset(feature_ids):
            raise ValueError(f"{pca_path}: invalid actual PCA feature subset")
        pca_features = set(pca_ids)
        sources.append(pca_path)
    if any(sha256(p) != digest for p, digest in digests.items()):
        raise RuntimeError(f"Producer/shared input changed while loading {path}")
    by_animal = {sample: {x for x, animal in zip(ids, samples) if animal == sample} for sample in MICE}
    return Dataset(tool, branch, population, path, ids, {x: i for i, x in enumerate(ids)},
                   np.asarray(samples, dtype=object), np.asarray(conditions, dtype=object), columns,
                   partitions, partition_columns, umaps, sources, by_animal, digests, set(feature_ids), pca_features)


def comparison(left, right, representation, resolution, category, identifier):
    common = sorted(left.index.keys() & right.index.keys())
    li = np.asarray([left.index[x] for x in common], dtype=int)
    ri = np.asarray([right.index[x] for x in common], dtype=int)
    if not len(common):
        raise ValueError(f"{identifier}: no common nuclei")
    if not np.array_equal(left.samples[li], right.samples[ri]) or not np.array_equal(left.conditions[li], right.conditions[ri]):
        raise ValueError(f"{identifier}: common nuclei disagree on animal or condition")
    a, b = left.partitions[(representation, resolution)][li], right.partitions[(representation, resolution)][ri]
    labels_a, labels_b = sorted(set(a), key=label_sort), sorted(set(b), key=label_sort)
    ai, bi = {x: i for i, x in enumerate(labels_a)}, {x: i for i, x in enumerate(labels_b)}
    counts = np.zeros((len(ai), len(bi)), dtype=np.int64)
    for (x, y), count in Counter(zip(a, b)).items():
        counts[ai[x], bi[y]] = count
    pairs = hungarian_maximize(counts)
    record = {"comparison_id": identifier, "comparison_category": category, "population": left.population,
              "representation": representation, "resolution": resolution, "left_tool": left.tool,
              "left_branch": left.branch, "left_cluster_column": left.partition_columns[(representation, resolution)],
              "right_tool": right.tool, "right_branch": right.branch,
              "right_cluster_column": right.partition_columns[(representation, resolution)],
              "left_nuclei": len(left.ids), "right_nuclei": len(right.ids), "common_nuclei": len(common),
              "left_only_nuclei": len(left.ids) - len(common), "right_only_nuclei": len(right.ids) - len(common),
              "same_nucleus_set": set(left.ids) == set(right.ids), "left_clusters_on_common": len(ai),
              "right_clusters_on_common": len(bi), **partition_metrics(counts),
              "hungarian_matched_fraction": sum(int(counts[i, j]) for i, j in pairs) / len(common),
              "left_normalization": "SCTransform_v2" if left.tool == "seurat_sct" else "library_size_log1p",
              "right_normalization": "SCTransform_v2" if right.tool == "seurat_sct" else "library_size_log1p",
              "left_representation_features": len(left.features), "right_representation_features": len(right.features),
              "common_representation_features": len(left.features & right.features),
              "same_representation_feature_panel": left.features == right.features,
              "left_PCA_actual_features": len(left.pca_features), "right_PCA_actual_features": len(right.pca_features),
              "common_PCA_actual_features": len(left.pca_features & right.pca_features),
              "same_PCA_actual_feature_panel": left.pca_features == right.pca_features,
              "left_metadata": str(left.path), "right_metadata": str(right.path)}
    return {"record": record, "left": left, "right": right, "ids": common, "li": li, "ri": ri,
            "a": a, "b": b, "labels_a": labels_a, "labels_b": labels_b, "counts": counts,
            "pairs": pairs, "mapping": {labels_a[i]: labels_b[j] for i, j in pairs}}


def plot_colors(item):
    left = item["labels_a"]
    mapping = {item["labels_b"][j]: item["labels_a"][i] for i, j in item["pairs"]}
    unmatched = [x for x in item["labels_b"] if x not in mapping]
    canonical = left + [f"right:{x}" for x in unmatched]
    colors = {x: plt.get_cmap("turbo")((i + .5) / max(1, len(canonical))) for i, x in enumerate(canonical)}
    return colors, mapping


def plot_comparison(item, pdf):
    r, counts = item["record"], item["counts"]
    fig, ax = plt.subplots(figsize=(11.7, 8.3), constrained_layout=True)
    fractions = counts / counts.sum(axis=1, keepdims=True)
    image = ax.imshow(fractions, aspect="auto", cmap="Blues", vmin=0, vmax=1)
    ax.set_xticks(range(len(item["labels_b"])), item["labels_b"], fontsize=8, rotation=90)
    ax.set_yticks(range(len(item["labels_a"])), item["labels_a"], fontsize=8)
    ax.set_xlabel(f'{r["right_tool"]} {r["right_branch"]}: cluster')
    ax.set_ylabel(f'{r["left_tool"]} {r["left_branch"]}: cluster')
    ax.set_title(f'{r["population"]} / {r["representation"]} / resolution {r["resolution"]}\n'
                 f'{r["comparison_category"]}; common nuclei {r["common_nuclei"]:,}; '
                 f'ARI {r["ARI"]:.3f}; NMI {r["NMI_arithmetic"]:.3f}', fontsize=12)
    fig.colorbar(image, ax=ax, label="Fraction within left cluster")
    for i, j in item["pairs"]:
        ax.plot(j, i, marker="s", markersize=11, markerfacecolor="none", markeredgecolor="darkorange")
    fig.text(.01, .005, "Orange squares: one-to-one maximum-overlap matches. Splits and merges remain visible. IDs have no biological rank.", fontsize=8)
    pdf.savefig(fig); plt.close(fig)


def plot_umaps(item, pdf, seed, max_points):
    r = item["record"]
    colors, mapping = plot_colors(item)
    n = len(item["ids"])
    selected = np.sort(np.random.default_rng(seed).choice(n, min(n, max_points), replace=False))
    fig, axes = plt.subplots(1, 2, figsize=(11.7, 5.8), constrained_layout=True)
    for ax, side, indices, labels in ((axes[0], "left", item["li"], item["a"]),
                                     (axes[1], "right", item["ri"], item["b"])):
        dataset = item[side]
        coordinates = dataset.umaps[r["representation"]][indices[selected]]
        category = labels[selected]
        if side == "right":
            category = [mapping.get(x, f"right:{x}") for x in category]
        ax.scatter(coordinates[:, 0], coordinates[:, 1], c=[colors[x] for x in category], s=1.5,
                   alpha=.55, rasterized=True, linewidths=0)
        ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")
        ax.set_title(f'{dataset.tool}: {dataset.branch}\n{r["population"]} / {r["representation"]}')
    fig.suptitle(f'Maximum-overlap colors; {min(n, max_points):,}/{n:,} common nuclei plotted', fontsize=12)
    fig.text(.01, .005, "UMAP rotations, reflections and axis ranges are cosmetic. Coordinates are not used as an agreement test.", fontsize=8)
    pdf.savefig(fig); plt.close(fig)


def plot_animal_proportions(item, pdf):
    r = item["record"]
    colors, mapping = plot_colors(item)
    samples = sorted(set(item["left"].samples[item["li"]]))
    fig, axes = plt.subplots(2, 1, figsize=(11.7, 8.3), constrained_layout=True)
    for ax, side, indices, labels in ((axes[0], "left", item["li"], item["a"]),
                                     (axes[1], "right", item["ri"], item["b"])):
        dataset = item[side]
        categories = labels if side == "left" else np.asarray([mapping.get(x, f"right:{x}") for x in labels])
        animal = dataset.samples[indices]
        bottom = np.zeros(len(samples))
        for category in colors:
            values = np.asarray([np.sum((animal == sample) & (categories == category)) / np.sum(animal == sample) for sample in samples])
            ax.bar(samples, values, bottom=bottom, color=colors[category], width=.85, label=category)
            bottom += values
        ax.set_ylim(0, 1); ax.set_ylabel("Fraction within animal")
        ax.tick_params(axis="x", rotation=45, labelsize=8)
        ax.set_title(f'{dataset.tool} {dataset.branch}: same common-nucleus cohort')
    fig.suptitle(f'{r["population"]} / {r["representation"]}: descriptive animal proportions', fontsize=12)
    fig.text(.01, .005, "Colors use the contingency match; denominators are common nuclei per animal. This is not a genotype proportion test.", fontsize=8)
    pdf.savefig(fig); plt.close(fig)


def rank_average(values):
    values = np.asarray(values)
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(len(values), dtype=float)
    start = 0
    while start < len(order):
        end = start + 1
        while end < len(order) and values[order[end]] == values[order[start]]:
            end += 1
        ranks[order[start:end]] = (start + end - 1) / 2 + 1
        start = end
    return ranks


def correlation(a, b):
    if len(a) < 2 or np.std(a) == 0 or np.std(b) == 0:
        return None
    return float(np.corrcoef(a, b)[0, 1])


def open_count_table(path):
    return gzip.open(path, "rt", newline="") if Path(path).suffix == ".gz" else Path(path).open(newline="")


def integer_count(value):
    try:
        if re.fullmatch(r"[0-9]+", value):
            return int(value)
        decimal = Decimal(value)
        if not decimal.is_finite() or decimal < 0 or decimal != decimal.to_integral_value():
            raise ValueError("not a finite nonnegative integer")
        return int(decimal)
    except (TypeError, InvalidOperation) as error:
        raise ValueError(f"Invalid exact integer count: {value}") from error


def load_count_table(path):
    rows, totals = {}, [0] * len(MICE)
    with open_count_table(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        if "gene_id" not in fields:
            raise ValueError(f"{path}: stable gene_id column missing")
        animal_fields = set(fields) - {"gene_id", "gene_symbol", "display_symbol"}
        if animal_fields != set(MICE):
            raise ValueError(f"{path}: exact all-twelve mouse columns required; observed {sorted(animal_fields)}")
        for row in reader:
            gene = row["gene_id"]
            if gene in MISSING or gene in rows:
                raise ValueError(f"{path}: missing/duplicate stable gene ID")
            counts = tuple(integer_count(row[sample]) for sample in MICE)
            rows[gene] = counts
            totals = [a + b for a, b in zip(totals, counts)]
    if not rows:
        raise ValueError(f"{path}: empty aggregate count table")
    return rows, totals


def exact_pseudobulk_comparison(root, datasets, source_paths, source_digests, allow_partial, main_resolution=.4):
    """Independently compare actual R/Python gene-by-animal integer aggregates."""
    summaries, unmatched, disagreements, conservation, identities, branch_receipts, missing = [], [], [], [], [], [], []
    expected_families = {"ATLAS_TOTAL_ALL", "ATLAS_TOTAL_EX", "ATLAS_TOTAL_IN", "ATLAS_TOTAL_CA1_PROS",
                         "SCANPY_ALLCELLS_harmony", "SCANPY_EXCITATORY_harmony", "SCANPY_INHIBITORY_harmony", "SCANPY_CA1_PROS_harmony"}
    def bind(path):
        path = Path(path)
        if not path.is_file():
            raise FileNotFoundError(path)
        source_paths.add(path)
        source_digests.setdefault(str(path), sha256(path))
        return path
    for branch in BRANCHES:
        pybase, rbase = root / "results/python_pseudobulk" / branch, root / "results/seurat" / branch / "pseudobulk"
        try:
            py_index, r_index = bind(pybase / "population_index.tsv"), bind(rbase / "POPULATION_INDEX.tsv")
            py_sentinel, r_sentinel = bind(pybase / "COMPLETE.json"), bind(rbase / "COMPLETE.json")
            py_receipt, r_receipt = json.loads(py_sentinel.read_text()), json.loads(r_sentinel.read_text())
            if py_receipt.get("status") != "COMPLETE" or py_receipt.get("branch") != branch:
                raise ValueError(f"{py_sentinel}: wrong/incomplete independent exporter")
            if r_receipt.get("status") != "EXECUTED_REVIEW_REQUIRED" or r_receipt.get("stage") != "pseudobulk" or r_receipt.get("count_source") != branch:
                raise ValueError(f"{r_sentinel}: wrong/incomplete R aggregate/model stage")
            py_records, r_records = list(read_tsv(py_index)), list(read_tsv(r_index))
            if set(r["family"] for r in py_records) != expected_families:
                raise ValueError(f"{py_index}: expected four common atlas totals and all four fresh Scanpy Harmony families")
            r_members = defaultdict(list)
            membership_file = bind(rbase / "POPULATION_MEMBERSHIPS.tsv")
            for row in read_tsv(membership_file):
                r_members[row["population_key"]].append(row["nucleus_id"])
            expected_r = {}
            all_r = []
            for record in r_records:
                ids = r_members.get(record["population_key"])
                if ids is None and int(record["n_nuclei"]) == 0:
                    ids = []
                if ids is None or len(set(ids)) != len(ids):
                    raise ValueError(f"{membership_file}: missing/duplicate population membership")
                digest = member_digest(ids)
                if record.get("membership_sha256") and record["membership_sha256"] != digest:
                    raise ValueError(f"{r_index}: R membership SHA256 disagrees with actual membership")
                key = (record["family"], record["label"], digest)
                if key in expected_r:
                    raise ValueError(f"{r_index}: ambiguous duplicate population join key")
                expected_r[key] = record
                all_r.append((key, record))
            matched_keys, branch_ok, atlas_totals, expected_py_genes = set(), True, {}, None
            gene_order, atlas_arrays, cluster_arrays = None, {}, {}
            cluster_to_atlas = {"SCANPY_ALLCELLS_harmony": "ATLAS_TOTAL_ALL", "SCANPY_EXCITATORY_harmony": "ATLAS_TOTAL_EX",
                                "SCANPY_INHIBITORY_harmony": "ATLAS_TOTAL_IN", "SCANPY_CA1_PROS_harmony": "ATLAS_TOTAL_CA1_PROS"}
            seen_py = set()
            for py in py_records:
                digest = py.get("membership_sha256", "")
                if not re.fullmatch(r"[0-9a-f]{64}", digest):
                    raise ValueError(f"{py_index}: missing valid sorted-membership SHA256")
                atlas_population = {"ATLAS_TOTAL_ALL": "allcells", "ATLAS_TOTAL_EX": "Excitatory", "ATLAS_TOTAL_IN": "Inhibitory", "ATLAS_TOTAL_CA1_PROS": "CA1_ProS"}.get(py["family"])
                if atlas_population:
                    expected_dataset = datasets.get(("scanpy", branch, atlas_population))
                    expected_ids = expected_dataset.ids if expected_dataset else None
                else:
                    context = de_context(py["family"], main_resolution=main_resolution)
                    label_branch = py.get("source_label_branch") or py_receipt.get("source_label_branch")
                    if label_branch != "original_full":
                        raise ValueError(f"{py_index}: explicit frozen original-full label branch required")
                    expected_dataset = datasets.get(("scanpy", label_branch, context[0])) if context else None
                    count_dataset = datasets.get(("scanpy", branch, context[0])) if context else None
                    count_ids = set(count_dataset.ids) if count_dataset else None
                    expected_ids = [x for x, label in zip(expected_dataset.ids, expected_dataset.partitions[context[1:3]])
                                    if label == py["label"] and x in count_ids] if expected_dataset and count_ids is not None and context else None
                if expected_ids is not None and (member_digest(expected_ids) != digest or len(expected_ids) != int(py["n_nuclei"])):
                    raise ValueError(f"{py_index}: Python population differs from authoritative current atlas/fresh cluster membership")
                key = (py["family"], py["label"], digest)
                if key in seen_py:
                    raise ValueError(f"{py_index}: ambiguous duplicate population join key")
                seen_py.add(key)
                r = expected_r.get(key)
                if r is None:
                    unmatched.append({"branch": branch, "side": "Python", "population_key": py["population_key"], "family": py["family"], "label": py["label"], "membership_sha256": digest, "status": "NO_EXACT_R_FAMILY_LABEL_MEMBERSHIP_MATCH"})
                    branch_ok = False
                    continue
                matched_keys.add(key)
                pypath = Path(py.get("counts_file") or pybase / py["population_key"] / "SOURCE_INTEGER_COUNTS_ALL_GENES_ALL12.tsv.gz")
                rpath = rbase / r["population_key"] / "SOURCE_INTEGER_COUNTS_ALL_GENES_ALL12.tsv"
                for path in (pypath, rpath):
                    if not path.resolve().is_relative_to((root / "results").resolve()):
                        raise ValueError(f"Historical/external aggregate input refused: {path}")
                    bind(path)
                expected_hash = py_receipt.get("population_counts_file_sha256", {}).get(py["population_key"])
                if not expected_hash or sha256(pypath) != expected_hash:
                    raise ValueError(f"{pypath}: independent exporter count-file receipt hash mismatch")
                a, totals_a = load_count_table(pypath)
                b, totals_b = load_count_table(rpath)
                if len(a) != py_receipt.get("n_genes"):
                    raise ValueError(f"{pypath}: aggregate does not contain every declared source gene")
                if expected_py_genes is None:
                    expected_py_genes = set(a)
                    gene_order = sorted(a)
                elif set(a) != expected_py_genes:
                    raise ValueError(f"{pypath}: gene universe varies by population")
                common = a.keys() & b.keys()
                different = 0
                maximum_difference = 0
                for gene in common:
                    for sample, av, bv in zip(MICE, a[gene], b[gene]):
                        if av != bv:
                            different += 1
                            maximum_difference = max(maximum_difference, abs(av - bv))
                            if len(disagreements) < 1000:
                                disagreements.append({"branch": branch, "family": py["family"], "label": py["label"], "gene_id": gene,
                                                      "sample": sample, "Python_integer_count": av, "R_integer_count": bv})
                exact = set(a) == set(b) and different == 0
                branch_ok = branch_ok and exact
                pa = np.asarray([a[g] for g in gene_order], dtype=np.int64)
                ra = np.asarray([b.get(g, (0,) * len(MICE)) for g in gene_order], dtype=np.int64)
                if py["family"].startswith("ATLAS_TOTAL_"):
                    atlas_arrays[py["family"]] = (pa, ra)
                elif py["family"] in cluster_to_atlas:
                    if py["family"] not in cluster_arrays:
                        cluster_arrays[py["family"]] = (np.zeros_like(pa), np.zeros_like(ra))
                    cluster_arrays[py["family"]][0][:] += pa
                    cluster_arrays[py["family"]][1][:] += ra
                summaries.append({"branch": branch, "family": py["family"], "label": py["label"], "membership_sha256": digest,
                                  "source_label_branch": py.get("source_label_branch", py_receipt.get("source_label_branch", "")),
                                  "Python_population_key": py["population_key"], "R_population_key": r["population_key"],
                                  "Python_gene_rows": len(a), "R_gene_rows": len(b), "common_gene_rows": len(common),
                                  "Python_only_gene_rows": len(a.keys() - b.keys()), "R_only_gene_rows": len(b.keys() - a.keys()),
                                  "animal_columns": len(MICE), "exact_common_integer_entries": len(common) * len(MICE) - different,
                                  "different_common_integer_entries": different, "maximum_absolute_integer_difference": maximum_difference,
                                  "all_geneID_by_animal_counts_equal": exact, "Python_counts": str(pypath), "R_counts": str(rpath)})
                if py["family"] == "ATLAS_TOTAL_ALL":
                    atlas_totals = {sample: (a_total, b_total) for sample, a_total, b_total in zip(MICE, totals_a, totals_b)}
                del a, b
            for key, r in all_r:
                if key not in matched_keys:
                    in_scope = r["family"] in expected_families
                    unmatched.append({"branch": branch, "side": "R", "population_key": r["population_key"], "family": r["family"], "label": r["label"], "membership_sha256": key[2],
                                      "status": "UNMATCHED_R_WITHIN_PYTHON_SCOPE" if in_scope else "R_EXTRA_OUTSIDE_PYTHON_ATLAS_AND_SCANPY_HARMONY_SCOPE"})
                    if in_scope:
                        branch_ok = False
            for cluster_family, atlas_family in cluster_to_atlas.items():
                available = cluster_family in cluster_arrays and atlas_family in atlas_arrays
                py_mismatch = int(np.sum(cluster_arrays[cluster_family][0] != atlas_arrays[atlas_family][0])) if available else None
                r_mismatch = int(np.sum(cluster_arrays[cluster_family][1] != atlas_arrays[atlas_family][1])) if available else None
                passed = available and py_mismatch == r_mismatch == 0
                branch_ok = branch_ok and passed
                identities.append({"branch": branch, "cluster_family": cluster_family, "atlas_total_family": atlas_family,
                                   "gene_rows": len(gene_order) if gene_order else 0, "animal_columns": len(MICE),
                                   "Python_different_entries": py_mismatch, "R_different_entries": r_mismatch,
                                   "all_required_population_matrices_available": available, "cluster_sum_equals_atlas_total_gene_by_animal": passed})
            py_conservation = bind(pybase / "source_count_conservation.tsv")
            r_conservation = bind(rbase / "STREAMED_COUNT_CONSERVATION.tsv")
            pc_rows, rc_rows = list(read_tsv(py_conservation)), list(read_tsv(r_conservation))
            pc = {r["sample"]: r for r in pc_rows}
            rc = {r["sample"]: r for r in rc_rows}
            if len(pc_rows) != 12 or len(rc_rows) != 12:
                raise ValueError(f"{branch}: duplicate/missing source conservation rows")
            if set(pc) != set(MICE) or set(rc) != set(MICE):
                raise ValueError(f"{branch}: missing all-twelve source/ALL conservation coverage")
            for sample in MICE:
                pa, ra = atlas_totals.get(sample, (None, None))
                p_source, r_source = integer_count(pc[sample]["source_total_UMIs"]), integer_count(rc[sample]["source_total_UMIs"])
                p_nuclei, r_nuclei = integer_count(pc[sample]["source_n_nuclei"]), integer_count(rc[sample]["source_nuclei"])
                expected_dataset = datasets.get(("scanpy", branch, "allcells"))
                expected_nuclei = len(expected_dataset.by_animal[sample]) if expected_dataset else None
                passed = pa == ra == p_source == r_source and p_nuclei == r_nuclei and (expected_nuclei is None or p_nuclei == expected_nuclei)
                branch_ok = branch_ok and passed
                conservation.append({"branch": branch, "sample": sample, "Python_atlas_ALL_UMIs": pa,
                                     "R_atlas_ALL_UMIs": ra, "Python_source_UMIs": p_source, "R_source_UMIs": r_source,
                                     "Python_source_nuclei": p_nuclei, "R_source_nuclei": r_nuclei, "authoritative_nuclei": expected_nuclei,
                                     "exact_ALL_source_conservation": passed})
            branch_receipts.append({"branch": branch, "expected_Python_populations": len(py_records), "exactly_joined_populations": len(matched_keys),
                                    "all_twelve_source_conservation_checked": set(atlas_totals) == set(MICE), "exact_integer_aggregation_PASS": bool(branch_ok and len(matched_keys) == len(py_records))})
        except FileNotFoundError as error:
            if not allow_partial:
                raise
            missing.append({"branch": branch, "reason": str(error)})
            branch_receipts.append({"branch": branch, "expected_Python_populations": None, "exactly_joined_populations": 0,
                                    "all_twelve_source_conservation_checked": False, "exact_integer_aggregation_PASS": False})
    for path, digest in source_digests.items():
        if sha256(path) != digest:
            raise RuntimeError(f"Input changed during exact aggregate comparison: {path}")
    return {"summaries": summaries, "unmatched": unmatched, "disagreements": disagreements,
            "conservation": conservation, "identities": identities, "branches": branch_receipts, "missing": missing,
            "PASS": len(branch_receipts) == 3 and all(r["exact_integer_aggregation_PASS"] for r in branch_receipts)}


def de_context(family, declared_population="", declared_column="", main_resolution=.4):
    column = declared_column or family
    if column.startswith("SCANPY_"):
        match = re.fullmatch(r"SCANPY_(ALLCELLS|EXCITATORY|EX|INHIBITORY|IN|CA1_PROS)_(unintegrated|harmony)", column)
        if not match:
            return None
        populations = {"ALLCELLS": "allcells", "EXCITATORY": "Excitatory", "EX": "Excitatory",
                       "INHIBITORY": "Inhibitory", "IN": "Inhibitory", "CA1_PROS": "CA1_ProS"}
        return (populations[match.group(1)], match.group(2), main_resolution, "cluster_" + match.group(2))
    column = column.removeprefix("SEURAT_")
    match = re.search(r"(unintegrated|harmony)_r([0-9.p]+)$", column)
    if not match:
        return None
    if declared_population:
        population = declared_population
    elif column.startswith("EX_"):
        population = "Excitatory"
    elif column.startswith("IN_"):
        population = "Inhibitory"
    elif column.startswith("CA1_PROS_CANDIDATE_"):
        population = "CA1_ProS"
    elif column.startswith("log_"):
        population = "allcells"
    else:
        return None
    return (population, match.group(1), float(match.group(2).replace("p", ".")), column)


def discover_de(root, datasets, manifest, source_paths, source_digests, main_resolution=.4):
    records, notes = [], []
    for tool, branch in itertools.product(("scanpy", "seurat"), BRANCHES):
        base = root / "results" / tool / branch / "pseudobulk"
        index = base / "POPULATION_INDEX.tsv"
        if not index.is_file():
            notes.append({"tool": tool, "branch": branch, "status": "optional_pseudobulk_not_available", "path": str(index)})
            continue
        source_paths.add(index)
        source_digests.setdefault(str(index), sha256(index))
        memberships = defaultdict(list)
        membership_file = base / "POPULATION_MEMBERSHIPS.tsv"
        if membership_file.is_file():
            source_paths.add(membership_file)
            source_digests.setdefault(str(membership_file), sha256(membership_file))
            for row in read_tsv(membership_file):
                memberships[row["population_key"]].append(row["nucleus_id"])
        for population in read_tsv(index):
            context = de_context(population["family"], main_resolution=main_resolution)
            logical_tool = "scanpy" if population["family"].startswith("SCANPY_") else ("seurat" if population["family"].startswith("SEURAT_") else tool)
            members = memberships.get(population["population_key"])
            if members is None and int(population.get("n_nuclei", -1)) == 0:
                members = []
            digest = member_digest(members) if members is not None else ""
            verified = False
            fixed_verified = False
            if context and members is not None:
                pop, representation, resolution, column = context
                dataset = datasets.get((logical_tool, branch, pop))
                if dataset and column in dataset.columns:
                    expected = [x for x, label in zip(dataset.ids, dataset.columns[column]) if label == population["label"]]
                    verified = digest == member_digest(expected)
                if population.get("source_label_branch") == "original_full":
                    original = datasets.get((logical_tool, "original_full", pop))
                    count_dataset = datasets.get((logical_tool, branch, pop))
                    if original and count_dataset and column in original.columns:
                        source_ids = set(count_dataset.ids)
                        fixed_expected = [x for x, label in zip(original.ids, original.columns[column]) if label == population["label"] and x in source_ids]
                        fixed_verified = digest == member_digest(fixed_expected)
            for method, effect, fdr in (("DESeq2", "log2FoldChange", "padj"), ("edgeR", "logFC", "FDR")):
                path = base / population["population_key"] / "AUDIT" / method / "PRIMARY" / "ALL_RESULTS_STABLE_ID_AUDIT.tsv"
                if not path.is_file():
                    notes.append({"tool": tool, "branch": branch, "status": "optional_population_model_not_available", "path": str(path)})
                    continue
                records.append({"tool": logical_tool, "model_execution_tool": tool, "branch": branch, **population, "method": method, "endpoint": "PRIMARY",
                                "results_path": path, "effect_column": effect, "fdr_column": fdr, "gene_id_column": "gene_id",
                                "contrast": "cKO_vs_WT", "contrast_provenance": "declared_native_pseudobulk_module_contract",
                                "membership_digest": digest, "fresh_cluster_membership_verified": verified,
                                "fixed_original_cluster_membership_verified": fixed_verified, "context": context})
    if manifest:
        source_paths.add(manifest)
        source_digests.setdefault(str(manifest), sha256(manifest))
        for row in read_tsv(manifest):
            for required in ("tool", "branch", "family", "label", "method", "results_path"):
                if not row.get(required):
                    raise ValueError(f"{manifest}: missing required {required}")
            path = Path(row["results_path"])
            if not path.is_absolute():
                path = root / path
            if not path.resolve().is_relative_to((root / "results").resolve()):
                raise ValueError(f"Historical/external DE input refused: {path}")
            if not path.is_file():
                raise FileNotFoundError(path)
            context = de_context(row["family"], row.get("population", ""), row.get("cluster_column", ""), main_resolution)
            verified = False
            digest = row.get("membership_digest", "")
            if context and digest:
                pop, representation, resolution, column = context
                dataset = datasets.get((row["tool"], row["branch"], pop))
                if dataset and column in dataset.columns:
                    expected = [x for x, label in zip(dataset.ids, dataset.columns[column]) if label == row["label"]]
                    verified = digest == member_digest(expected)
            row.update(results_path=path, gene_id_column=row.get("gene_id_column") or "gene_id",
                       effect_column=row.get("effect_column") or ("logFC" if row["method"] == "edgeR" else "log2FoldChange"),
                       fdr_column=row.get("fdr_column") or ("FDR" if row["method"] == "edgeR" else "padj"),
                       endpoint=row.get("endpoint") or "PRIMARY", context=context, membership_digest=digest,
                       fresh_cluster_membership_verified=verified, contrast_provenance="explicit_manifest")
            records.append(row)
    unique = {}
    for record in records:
        key = (record["tool"], record["branch"], record["family"], record["label"], record["method"], record["endpoint"])
        if key in unique and unique[key]["results_path"].resolve() != record["results_path"].resolve():
            raise ValueError(f"Ambiguous duplicate DE input: {key}")
        unique[key] = record
        source_paths.add(record["results_path"])
        source_digests.setdefault(str(record["results_path"]), sha256(record["results_path"]))
    return list(unique.values()), notes


@functools.lru_cache(maxsize=4)
def load_de(path, gene_column, effect_column, fdr_column):
    output = {}
    for row in read_tsv(path):
        for key in (gene_column, effect_column, fdr_column):
            if key not in row:
                raise ValueError(f"{path}: missing stable-ID/effect/FDR column {key}")
        gene = row[gene_column]
        if gene in MISSING or gene in output:
            raise ValueError(f"{path}: missing or duplicated stable gene ID")
        effect, fdr = numeric(row[effect_column]), numeric(row[fdr_column])
        if fdr is not None and not 0 <= fdr <= 1:
            raise ValueError(f"{path}: FDR outside [0,1]")
        output[gene] = (effect, fdr)
    return output


def de_pairs(records, comparisons, minimum_overlap=.5, minimum_nuclei=20):
    maps = {}
    for item in comparisons:
        r = item["record"]
        counts = item["counts"]
        adequate = {item["labels_a"][i]: item["labels_b"][j] for i, j in item["pairs"]
                    if counts[i, j] >= minimum_nuclei and counts[i, j] / counts[i, :].sum() >= minimum_overlap
                    and counts[i, j] / counts[:, j].sum() >= minimum_overlap}
        maps[(r["left_tool"], r["left_branch"], r["right_tool"], r["right_branch"],
              r["population"], r["representation"], r["resolution"])] = adequate
    for left, right in itertools.combinations(sorted(records, key=lambda r: (r["tool"], r["branch"], r["family"], r["label"], r["method"])), 2):
        if left["endpoint"] != "PRIMARY" or right["endpoint"] != "PRIMARY":
            continue
        if not left.get("contrast") or left.get("contrast") != right.get("contrast"):
            continue
        same_tool_branch = (left["tool"], left["branch"]) == (right["tool"], right["branch"])
        same_population = (left["family"], left["label"]) == (right["family"], right["label"])
        if same_tool_branch and same_population and left["method"] != right["method"]:
            yield left, right, "DESeq2_vs_edgeR_same_population"
            continue
        if same_tool_branch or left["method"] != right["method"]:
            continue
        if same_population and (left["family"].startswith(("reference_", "mmc_", "TOTAL_", "ATLAS_"))
                                or (left.get("membership_digest") and left["membership_digest"] == right.get("membership_digest"))
                                or (left.get("fixed_original_cluster_membership_verified") and right.get("fixed_original_cluster_membership_verified"))):
            yield left, right, "shared_annotation_or_fixed_population_reproducibility"
            continue
        if not left.get("fresh_cluster_membership_verified") or not right.get("fresh_cluster_membership_verified"):
            continue
        lc, rc = left.get("context"), right.get("context")
        if not lc or not rc or lc[:3] != rc[:3]:
            continue
        population, representation, resolution = lc[:3]
        key = (left["tool"], left["branch"], right["tool"], right["branch"], population, representation, resolution)
        matched = maps.get(key, {})
        if matched.get(left["label"]) == right["label"]:
            yield left, right, "fresh_cluster_populations_matched_by_common_nucleus_overlap"
        else:
            reverse = maps.get((right["tool"], right["branch"], left["tool"], left["branch"], population, representation, resolution), {})
            if reverse.get(right["label"]) == left["label"]:
                yield left, right, "fresh_cluster_populations_matched_by_common_nucleus_overlap"


def compare_de(left, right, scope, identifier, gene_writer):
    a = load_de(str(left["results_path"]), left["gene_id_column"], left["effect_column"], left["fdr_column"])
    b = load_de(str(right["results_path"]), right["gene_id_column"], right["effect_column"], right["fdr_column"])
    genes = sorted(a.keys() & b.keys())
    finite = [g for g in genes if a[g][0] is not None and b[g][0] is not None]
    x, y = np.asarray([a[g][0] for g in finite]), np.asarray([b[g][0] for g in finite])
    tested = [g for g in finite if a[g][1] is not None and b[g][1] is not None]
    sig_a, sig_b = {g for g in tested if a[g][1] < .05}, {g for g in tested if b[g][1] < .05}
    nonzero = (x != 0) & (y != 0)
    union = sig_a | sig_b
    record = {"de_comparison_id": identifier, "comparison_scope": scope, "contrast": left["contrast"],
              "left_tool": left["tool"], "left_branch": left["branch"], "left_family": left["family"], "left_label": left["label"], "left_model": left["method"],
              "right_tool": right["tool"], "right_branch": right["branch"], "right_family": right["family"], "right_label": right["label"], "right_model": right["method"],
              "left_all_gene_rows": len(a), "right_all_gene_rows": len(b), "common_gene_ids": len(genes),
              "common_finite_effects": len(finite), "common_finite_FDR": len(tested),
              "effect_Pearson": correlation(x, y), "effect_Spearman": correlation(rank_average(x), rank_average(y)),
              "effect_median_absolute_difference": float(np.median(np.abs(x - y))) if len(x) else None,
              "same_nonzero_effect_sign_fraction": float(np.mean(np.sign(x[nonzero]) == np.sign(y[nonzero]))) if np.any(nonzero) else None,
              "nonzero_effect_sign_denominator": int(nonzero.sum()), "FDR05_left_on_common_tested": len(sig_a),
              "FDR05_right_on_common_tested": len(sig_b), "FDR05_overlap": len(sig_a & sig_b),
              "FDR05_jaccard_on_common_tested": len(sig_a & sig_b) / len(union) if union else None,
              "same_exact_population_membership": bool(left.get("membership_digest") and left["membership_digest"] == right.get("membership_digest")),
              "left_source_label_branch": left.get("source_label_branch", left["branch"]), "right_source_label_branch": right.get("source_label_branch", right["branch"]),
              "left_population_nuclei": left.get("n_nuclei", ""), "right_population_nuclei": right.get("n_nuclei", ""),
              "left_results": str(left["results_path"]), "right_results": str(right["results_path"])}
    for gene in genes:
        gene_writer.writerow({"de_comparison_id": identifier, "gene_id": gene, "left_effect": a[gene][0],
                              "right_effect": b[gene][0], "left_FDR": a[gene][1], "right_FDR": b[gene][1]})
    return record, x, y


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--main-resolution", type=float, default=.4)
    parser.add_argument("--seed", type=int, default=20261013)
    parser.add_argument("--max-umap-points", type=int, default=20000)
    parser.add_argument("--max-de-plot-pages", type=int, default=16)
    parser.add_argument("--min-de-cluster-overlap", type=float, default=.5,
                        help="Require this reciprocal common-cohort overlap fraction for matched-cluster DE")
    parser.add_argument("--min-de-cluster-nuclei", type=int, default=20)
    parser.add_argument("--de-manifest", type=Path)
    parser.add_argument("--allow-partial", action="store_true", help="Explicit exploratory mode; COMPLETE status becomes partial")
    args = parser.parse_args()
    root, out = args.run_root.resolve(), args.outdir.resolve()
    if out.exists() and any(out.iterdir()):
        raise FileExistsError(f"Refusing to overwrite nonempty comparison output: {out}")
    if args.max_umap_points < 1 or args.max_de_plot_pages < 0:
        raise ValueError("Invalid plot limit")
    if not 0 < args.min_de_cluster_overlap <= 1 or args.min_de_cluster_nuclei < 1:
        raise ValueError("Invalid matched-cluster DE overlap gate")
    initialize_plotting()
    datasets, missing, source_paths = {}, [], set()
    for tool, branch, population in itertools.product(("scanpy", "seurat"), BRANCHES, POPULATIONS):
        try:
            dataset = load_dataset(root, tool, branch, population, args.main_resolution)
        except FileNotFoundError as error:
            if not args.allow_partial:
                raise
            missing.append({"tool": tool, "branch": branch, "population": population, "reason": str(error)})
            continue
        datasets[dataset.key] = dataset
        source_paths.update(dataset.sources)
    try:
        sct = load_dataset(root, "seurat_sct", "original_full", "allcells", args.main_resolution)
        datasets[sct.key] = sct
        source_paths.update(sct.sources)
    except FileNotFoundError as error:
        if not args.allow_partial:
            raise
        missing.append({"tool": "seurat_sct", "branch": "original_full", "population": "allcells", "reason": str(error)})
    if not datasets:
        raise ValueError("No fresh result datasets")
    source_hashes_before = {}
    for dataset in datasets.values():
        for path, digest in dataset.source_digests.items():
            if path in source_hashes_before and source_hashes_before[path] != digest:
                raise RuntimeError(f"Shared input changed between producer loads: {path}")
            source_hashes_before[path] = digest
    exact = exact_pseudobulk_comparison(root, datasets, source_paths, source_hashes_before, args.allow_partial, args.main_resolution)
    sct_features = []
    log = datasets.get(("seurat", "original_full", "allcells"))
    sct = datasets.get(("seurat_sct", "original_full", "allcells"))
    if sct and log:
        coverage_path = sct.path.parent / "SCT_MODELED_FEATURE_COVERAGE.tsv"
        if not coverage_path.is_file():
            raise FileNotFoundError(coverage_path)
        source_paths.add(coverage_path)
        source_hashes_before[str(coverage_path)] = sha256(coverage_path)
        coverage_rows = list(read_tsv(coverage_path))
        coverage = {r["gene_id"]: r for r in coverage_rows}
        if len(coverage_rows) != len(coverage) or not sct.features.issubset(log.features) or set(coverage) != log.features:
            raise ValueError("SCT modeled-feature intersection differs from the declared original common panel")
        for gene in sorted(log.features):
            row = coverage[gene]
            included = row["included_in_SCT_representation"] == "TRUE"
            n_captures = integer_count(row["n_captures_with_modeled_residual"])
            if included != (gene in sct.features) or included != (n_captures == 12) or n_captures > 12:
                raise ValueError("SCT modeled coverage/actual feature panel mismatch")
            sct_features.append({"gene_id": gene, "in_log_panel": True, "in_SCT_panel": included,
                                 "n_captures_with_modeled_residual": n_captures, "reason": row["reason"]})
    comparisons = []
    def add(left, right, category):
        if category in ("cross_tool_same_count_branch", "paired_count_correction_sensitivity", "normalization_and_modeled_feature_panel_sensitivity") and set(left.ids) != set(right.ids) and not args.allow_partial:
            raise ValueError(f"Same-cohort comparison has different nuclei: {left.key} vs {right.key}")
        keys = sorted(left.partitions.keys() & right.partitions.keys())
        for representation, resolution in keys:
            identifier = f"C{len(comparisons) + 1:04d}"
            comparisons.append(comparison(left, right, representation, resolution, category, identifier))
    for branch, population in itertools.product(BRANCHES, POPULATIONS):
        left, right = datasets.get(("scanpy", branch, population)), datasets.get(("seurat", branch, population))
        if left and right:
            add(left, right, "cross_tool_same_count_branch")
    for tool, population in itertools.product(("scanpy", "seurat"), POPULATIONS):
        for a, b in itertools.combinations(BRANCHES, 2):
            left, right = datasets.get((tool, a, population)), datasets.get((tool, b, population))
            if left and right:
                scope = "paired_count_correction_sensitivity" if a == "paired_original" else (
                    "full_vs_intersection_cohort_and_reclustering" if b == "paired_original" else "cohort_and_count_correction_combined")
                add(left, right, scope)
    if log and sct:
        add(log, sct, "normalization_and_modeled_feature_panel_sensitivity")
    out.mkdir(parents=True, exist_ok=True)
    write_tsv(out / "normalization_sensitivity_features.tsv", sct_features, ["gene_id", "in_log_panel", "in_SCT_panel", "n_captures_with_modeled_residual", "reason"])
    write_tsv(out / "exact_pseudobulk_counts.tsv", exact["summaries"], ["branch", "family", "label", "membership_sha256", "source_label_branch", "Python_population_key", "R_population_key", "Python_gene_rows", "R_gene_rows", "common_gene_rows", "Python_only_gene_rows", "R_only_gene_rows", "animal_columns", "exact_common_integer_entries", "different_common_integer_entries", "maximum_absolute_integer_difference", "all_geneID_by_animal_counts_equal", "Python_counts", "R_counts"])
    write_tsv(out / "exact_pseudobulk_unmatched.tsv", exact["unmatched"], ["branch", "side", "population_key", "family", "label", "membership_sha256", "status"])
    write_tsv(out / "exact_pseudobulk_first_1000_disagreements.tsv", exact["disagreements"], ["branch", "family", "label", "gene_id", "sample", "Python_integer_count", "R_integer_count"])
    write_tsv(out / "exact_pseudobulk_conservation.tsv", exact["conservation"], ["branch", "sample", "Python_atlas_ALL_UMIs", "R_atlas_ALL_UMIs", "Python_source_UMIs", "R_source_UMIs", "Python_source_nuclei", "R_source_nuclei", "authoritative_nuclei", "exact_ALL_source_conservation"])
    write_tsv(out / "exact_pseudobulk_partition_identities.tsv", exact["identities"], ["branch", "cluster_family", "atlas_total_family", "gene_rows", "animal_columns", "Python_different_entries", "R_different_entries", "all_required_population_matrices_available", "cluster_sum_equals_atlas_total_gene_by_animal"])
    write_tsv(out / "exact_pseudobulk_branch_coverage.tsv", exact["branches"], ["branch", "expected_Python_populations", "exactly_joined_populations", "all_twelve_source_conservation_checked", "exact_integer_aggregation_PASS"])
    metrics = [item["record"] for item in comparisons]
    if metrics:
        write_tsv(out / "partition_metrics.tsv", metrics, list(metrics[0]))
    contingency, matches, attrition, proportions = [], [], [], []
    for item in comparisons:
        r, counts = item["record"], item["counts"]
        row_totals, column_totals = counts.sum(axis=1), counts.sum(axis=0)
        for i, j in zip(*np.nonzero(counts)):
            contingency.append({"comparison_id": r["comparison_id"], "left_cluster": item["labels_a"][i],
                                "right_cluster": item["labels_b"][j], "n_common_nuclei": int(counts[i, j])})
        matched_a, matched_b = set(), set()
        for i, j in item["pairs"]:
            matched_a.add(i); matched_b.add(j)
            n, a, b = int(counts[i, j]), int(row_totals[i]), int(column_totals[j])
            matches.append({"comparison_id": r["comparison_id"], "match_status": "one_to_one_matched",
                            "left_cluster": item["labels_a"][i], "right_cluster": item["labels_b"][j],
                            "overlap_nuclei": n, "left_cluster_nuclei": a, "right_cluster_nuclei": b,
                            "left_overlap_fraction": n / a, "right_overlap_fraction": n / b, "jaccard": n / (a + b - n)})
        for side, labels, totals, matched in (("left", item["labels_a"], row_totals, matched_a), ("right", item["labels_b"], column_totals, matched_b)):
            for index, label in enumerate(labels):
                if index not in matched:
                    matches.append({"comparison_id": r["comparison_id"], "match_status": f"unmatched_{side}",
                                    "left_cluster": label if side == "left" else "", "right_cluster": label if side == "right" else "",
                                    "overlap_nuclei": 0, "left_cluster_nuclei": int(totals[index]) if side == "left" else 0,
                                    "right_cluster_nuclei": int(totals[index]) if side == "right" else 0})
        animals = sorted(set(item["left"].samples) | set(item["right"].samples))
        for sample in animals:
            own = [item[side].by_animal.get(sample, set()) for side in ("left", "right")]
            common_n = len(own[0] & own[1])
            attrition.append({"comparison_id": r["comparison_id"], "sample": sample, "left_nuclei": len(own[0]),
                              "right_nuclei": len(own[1]), "common_nuclei": common_n,
                              "left_only_nuclei": len(own[0] - own[1]), "right_only_nuclei": len(own[1] - own[0])})
    for dataset in datasets.values():
        for (representation, resolution), labels in sorted(dataset.partitions.items()):
            cluster_labels = sorted(set(labels), key=label_sort)
            for sample in sorted(set(dataset.samples)):
                mask = dataset.samples == sample
                denominator = int(mask.sum())
                condition = set(dataset.conditions[mask])
                if len(condition) != 1:
                    raise ValueError(f"{dataset.path}: animal has conflicting conditions")
                counts = Counter(labels[mask])
                for cluster in cluster_labels:
                    proportions.append({"tool": dataset.tool, "branch": dataset.branch, "population": dataset.population,
                                        "representation": representation, "resolution": resolution, "sample": sample,
                                        "condition": next(iter(condition)), "cluster": cluster, "cluster_nuclei": counts[cluster],
                                        "animal_population_nuclei": denominator, "within_animal_fraction": counts[cluster] / denominator})
    write_tsv(out / "contingency_nonzero.tsv", contingency, ["comparison_id", "left_cluster", "right_cluster", "n_common_nuclei"])
    write_tsv(out / "cluster_matches.tsv", matches, ["comparison_id", "match_status", "left_cluster", "right_cluster", "overlap_nuclei", "left_cluster_nuclei", "right_cluster_nuclei", "left_overlap_fraction", "right_overlap_fraction", "jaccard"])
    write_tsv(out / "roster_attrition_per_animal.tsv", attrition, ["comparison_id", "sample", "left_nuclei", "right_nuclei", "common_nuclei", "left_only_nuclei", "right_only_nuclei"])
    write_tsv(out / "cluster_proportions_per_animal.tsv", proportions, ["tool", "branch", "population", "representation", "resolution", "sample", "condition", "cluster", "cluster_nuclei", "animal_population_nuclei", "within_animal_fraction"])
    selected = [x for x in comparisons if x["record"]["resolution"] == args.main_resolution
                and x["record"]["representation"] == "harmony"
                and x["record"]["comparison_category"] in ("cross_tool_same_count_branch", "paired_count_correction_sensitivity", "full_vs_intersection_cohort_and_reclustering", "normalization_and_modeled_feature_panel_sensitivity")]
    if selected:
        with PdfPages(out / "partition_comparison.pdf") as pdf:
            for item in selected:
                plot_comparison(item, pdf)
        with PdfPages(out / "umap_comparison.pdf") as pdf:
            for item in comparisons:
                if item["record"]["resolution"] == args.main_resolution and item["record"]["comparison_category"] in ("cross_tool_same_count_branch", "normalization_and_modeled_feature_panel_sensitivity"):
                    plot_umaps(item, pdf, args.seed, args.max_umap_points)
        with PdfPages(out / "animal_proportions.pdf") as pdf:
            for item in selected:
                if item["record"]["comparison_category"] in ("cross_tool_same_count_branch", "normalization_and_modeled_feature_panel_sensitivity"):
                    plot_animal_proportions(item, pdf)
    de_records, de_notes = discover_de(root, datasets, args.de_manifest, source_paths, source_hashes_before, args.main_resolution)
    de_metrics = []
    de_fields = ["de_comparison_id", "gene_id", "left_effect", "right_effect", "left_FDR", "right_FDR"]
    pairs = list(de_pairs(de_records, comparisons, args.min_de_cluster_overlap, args.min_de_cluster_nuclei))
    with gzip.open(out / "de_gene_pairs.tsv.gz", "wt", newline="") as compressed:
        writer = csv.DictWriter(compressed, fieldnames=de_fields, delimiter="\t")
        writer.writeheader()
        pdf = PdfPages(out / "de_effect_comparison.pdf") if pairs and args.max_de_plot_pages else None
        try:
            for index, (left, right, scope) in enumerate(pairs):
                record, x, y = compare_de(left, right, scope, f"D{index + 1:04d}", writer)
                de_metrics.append(record)
                if pdf and index < args.max_de_plot_pages:
                    fig, ax = plt.subplots(figsize=(8.3, 8.3), constrained_layout=True)
                    ax.scatter(x, y, s=2, alpha=.25, rasterized=True)
                    if len(x):
                        low, high = min(float(x.min()), float(y.min())), max(float(x.max()), float(y.max()))
                        ax.plot([low, high], [low, high], color="grey", linestyle="--", linewidth=1)
                    ax.set_xlabel(f'{left["tool"]} {left["branch"]} {left["method"]}: log2 effect')
                    ax.set_ylabel(f'{right["tool"]} {right["branch"]} {right["method"]}: log2 effect')
                    ax.set_title(f'{record["de_comparison_id"]}: {left["family"]} {left["label"]}\n{scope}\nFinite stable gene IDs: {len(x):,}', fontsize=10)
                    fig.text(.01, .005, "Raw model effect estimates; FDR families and cluster memberships may differ. Shared atlas labels are common evidence.", fontsize=8)
                    pdf.savefig(fig); plt.close(fig)
        finally:
            if pdf:
                pdf.close()
    if de_metrics:
        write_tsv(out / "de_effect_metrics.tsv", de_metrics, list(de_metrics[0]))
    else:
        write_tsv(out / "de_effect_metrics.tsv", [], ["de_comparison_id", "comparison_scope", "common_finite_effects"])
    write_tsv(out / "optional_de_status.tsv", de_notes, ["tool", "branch", "status", "path"])
    changed = [path for path, digest in source_hashes_before.items() if sha256(path) != digest]
    if changed:
        raise RuntimeError(f"Input files changed during comparison; COMPLETE withheld: {changed}")
    source_receipts = [{"path": path, "sha256": digest} for path, digest in sorted(source_hashes_before.items())]
    write_tsv(out / "INPUT_SOURCE_SHA256.tsv", source_receipts, ["path", "sha256"])
    flags = [f'{r["comparison_id"]}: compared nucleus sets differ' for r in metrics if not r["same_nucleus_set"]
             and r["comparison_category"] in ("cross_tool_same_count_branch", "paired_count_correction_sensitivity")]
    if not exact["PASS"]:
        flags.append("Independent R/Python exact integer pseudobulk aggregation did not satisfy complete-scope PASS; inspect mismatches/coverage.")
    notes = ["ARI/NMI use all common labeled nuclei; cluster matching does not alter assignments.",
             "Resolution-grid results are compared only at resolutions present in both sources.",
             "Full/intersection comparisons combine cohort restriction with reclustering; paired original/CellBender isolates count source only when rosters match.",
             "UMAP rotations/reflections and coordinate ranges are cosmetic and are not tested as numerical agreement.",
             "Per-animal proportions are descriptive and do not replace animal-level abundance inference.",
             "Common original MapMyCells/atlas annotations are shared evidence, not independent biological truth.",
             "Optional shared-atlas DE comparisons with identical inputs assess reproducibility, not independent annotation accuracy.",
             "Cross-cluster DE effect comparisons use overlap-matched populations whose memberships can differ; they do not isolate a count correction effect.",
             "Seurat log/SCT comparisons share original nuclei/count source; SCT uses the all-capture modeled-feature intersection, so differing panels accompany normalization differences.",
             "Finite effect/FDR denominators, original BH-family sizes and missing values are reported; a nonsignificant result is not proof of absence.",
             "No automatic consistency verdict or biological significance claim is assigned."]
    report = ["# Fresh Scanpy / Seurat comparison", "", f"Loaded {len(datasets)}/25 datasets; computed {len(metrics)} partition and {len(de_metrics)} optional DE comparisons.",
              f"Independent R/Python exact integer aggregation: {'PASS' if exact['PASS'] else 'NOT_PASS'}; scoped coverage is documented per branch.", "", *[f"- {note}" for note in notes]]
    if flags:
        report += ["", "Roster discrepancies requiring review:", *[f"- {flag}" for flag in flags]]
    if missing:
        report += ["", "Explicit partial mode: missing inputs are recorded in COMPLETE.json."]
    report += ["", "Metrics: partition_metrics.tsv, cluster_matches.tsv, contingency_nonzero.tsv, roster_attrition_per_animal.tsv and cluster_proportions_per_animal.tsv.",
               "Independent aggregates: exact_pseudobulk_counts.tsv, exact_pseudobulk_unmatched.tsv, exact_pseudobulk_conservation.tsv, exact_pseudobulk_branch_coverage.tsv and first 1000 disagreements; complete mismatch counts remain in the population metrics.",
               "Optional DE: de_effect_metrics.tsv, de_gene_pairs.tsv.gz and optional_de_status.tsv. All input source hashes are recorded in INPUT_SOURCE_SHA256.tsv.",
               "PDFs show the primary-resolution comparisons; all available matched resolution-grid metrics remain in the tables. Zero contingency entries are implicit."]
    (out / "COMPARISON_REVIEW.md").write_text("\n".join(report) + "\n")
    outputs = [{"path": str(path), "sha256": sha256(path)} for path in sorted(out.iterdir()) if path.is_file()]
    receipt = {"status": "PARTIAL_EXECUTED_REVIEW_REQUIRED" if missing or exact["missing"] else "EXECUTED_REVIEW_REQUIRED", "expected_datasets": 25,
               "observed_datasets": len(datasets), "missing_inputs": missing, "partition_comparisons": len(metrics),
               "optional_de_comparisons": len(de_metrics), "review_flags": flags, "notes": notes,
               "exact_integer_pseudobulk_aggregation_PASS": exact["PASS"], "exact_aggregation_branch_coverage": exact["branches"], "missing_exact_aggregation_inputs": exact["missing"],
               "parameters": {"main_resolution": args.main_resolution, "seed": args.seed, "max_umap_points": args.max_umap_points,
                              "min_de_cluster_overlap": args.min_de_cluster_overlap, "min_de_cluster_nuclei": args.min_de_cluster_nuclei},
               "source_script_sha256": sha256(Path(__file__)), "input_sources": source_receipts, "outputs": outputs,
               "no_scientific_fitting_performed": True, "automatic_consistency_verdict": None}
    (out / "COMPLETE.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps({k: receipt[k] for k in ("status", "observed_datasets", "partition_comparisons", "optional_de_comparisons", "review_flags")}))


if __name__ == "__main__":
    main()
