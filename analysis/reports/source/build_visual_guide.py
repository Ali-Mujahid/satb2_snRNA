#!/usr/bin/env python3
"""Offline guide and figures from completed produced data; no fits/count loading."""
from __future__ import annotations

import argparse
import base64
import colorsys
import hashlib
import html
import json
import os
from pathlib import Path
import subprocess
import tempfile
import textwrap
from urllib.parse import quote

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "satb2_visual_guide_matplotlib_cache"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.patheffects as pe
import numpy as np
import pandas as pd

plt.rcParams.update({"pdf.fonttype": 42, "svg.fonttype": "none", "font.size": 10})
KEYS = ["tool", "branch", "population", "cluster_scope", "cluster_id"]
BRANCHES = ["original_full", "paired_original", "paired_cellbender"]
POP_LABELS = {"allcells": "All nuclei", "Excitatory": "Excitatory neurons", "Inhibitory": "Inhibitory neurons", "CA1_ProS": "CA1-ProS candidates"}
CONFIDENCE_LABELS = {
    "REFERENCE_AND_FULL_COHORT_COARSE_MARKERS_CORROBORATE": "Atlas + full-cohort coarse marker support",
    "REFERENCE_AND_MARKERS_CORROBORATE_COARSE": "Atlas + initial coarse marker support",
    "REFERENCE_SUPPORTED_MARKER_REVIEW_REQUIRED": "Atlas supported; marker review needed",
    "RARE_REFERENCE_MARKER_EVIDENCE_REVIEW_REQUIRED": "Rare group; further marker review needed",
    "CONFLICTING_EVIDENCE": "Conflicting signatures",
    "UNCERTAIN": "Uncertain",
}


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for part in iter(lambda: f.read(8 * 1024 * 1024), b""): h.update(part)
    return h.hexdigest()


class Evidence:
    def __init__(self):
        self.inputs = {}
        self.expected_outputs = {}
    def bind(self, path):
        path = Path(path).resolve(strict=True)
        actual = sha(path)
        if str(path) in self.expected_outputs and self.expected_outputs[str(path)] != actual:
            raise ValueError(f"Producer output hash differs: {path}")
        if str(path) in self.inputs and self.inputs[str(path)] != actual: raise ValueError("Guide input changed")
        self.inputs[str(path)] = actual
        return path
    def table(self, path, **kwargs):
        path = self.bind(path)
        return pd.read_csv(path, sep="\t", dtype={"cluster_id": str, "source_cluster_id": str, "nucleus_id": str}, **kwargs)
    def stage(self, directory, status, receipt_name="COMPLETE.json"):
        directory = Path(directory).resolve(strict=True)
        receipt = json.loads(self.bind(directory / receipt_name).read_text())
        if receipt.get("status") != status: raise ValueError(f"Required producer is incomplete: {directory}")
        for relative, digest in receipt.get("output_sha256", {}).items():
            path = (directory / relative).resolve()
            if not path.is_relative_to(directory): raise ValueError("Producer receipt path escapes its directory")
            if str(path) in self.expected_outputs and self.expected_outputs[str(path)] != digest:
                raise ValueError("Producer receipts disagree")
            self.expected_outputs[str(path)] = digest
        return receipt
    def verify(self):
        for path, digest in self.inputs.items():
            if sha(path) != digest: raise ValueError("Input changed before guide publication")


def boolean(series):
    v = series.astype(str).str.lower().map({"true": True, "false": False, "1": True, "0": False})
    if v.isna().any(): raise ValueError("Invalid annotation mask boolean")
    return v.astype(bool)


def color(label):
    low = str(label).lower()
    if "conflicting" in low: return "#856a9c"
    if any(x in low for x in ["mixed", "weak", "unresolved", "not assigned"]): return "#a7afb5"
    if "superficial-associated" in low: return "#458ab3"
    if "deep-associated" in low: return "#b66942"
    if "pvalb" in low: return "#4379b6"
    if "sst" in low: return "#49956b"
    if "vip" in low: return "#b768a0"
    if "cck" in low or "sncg" in low: return "#8762a6"
    h = int(hashlib.sha256(str(label).encode()).hexdigest()[:8], 16) / 2**32
    return matplotlib.colors.to_hex(colorsys.hsv_to_rgb(h, .5, .68))


def local_link(path, out):
    return quote(os.path.relpath(Path(path).resolve(), out), safe="/")


def embed(path, evidence):
    path = evidence.bind(path)
    return "data:image/png;base64," + base64.b64encode(path.read_bytes()).decode()


def thumbnail(pdf_path, page, out, evidence):
    """Standard PDF-page raster preview; native figures and their data unchanged."""
    pdf_path = evidence.bind(pdf_path)
    key = hashlib.sha256((evidence.inputs[str(pdf_path)] + "__" + str(page)).encode()).hexdigest()[:20]
    folder = out / "thumbnails"; folder.mkdir(exist_ok=True)
    stem = folder / (key + "__page_" + str(page))
    path = stem.with_suffix(".png")
    if not path.exists():
        renderer = "/Users/mujahid357/.cache/codex-runtimes/codex-primary-runtime/dependencies/bin/override/pdftoppm"
        subprocess.run([renderer, "-f", str(page), "-singlefile", "-r", "90", "-png", str(pdf_path), str(stem)],
                       check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return embed(path, evidence)


def save_figure(fig, stem, out):
    records = {}
    for ext in ["pdf", "svg", "png"]:
        path = out / (stem + "." + ext)
        fig.savefig(path, dpi=170)
        records[ext] = str(path.resolve())
    plt.close(fig)
    return records


def diagnostic_figure(evidence, cfg, out):
    root = Path(cfg["ca1_diagnostic_dir"])
    receipt = evidence.stage(root, "COMPLETE")
    if receipt.get("primary_memberships_unchanged") is not True: raise ValueError("Diagnostic changed primary memberships")
    latent = evidence.table(root / "LATENT_SPACE_COMPARISON.tsv")
    graph = evidence.table(root / "NEIGHBOR_GRAPH_COMPARISON.tsv")
    controls = evidence.table(root / "LEIDEN_FIXED_GRAPH_CONTROLS.tsv")
    convergence = evidence.table(root / "PYTHON_HARMONY_CONVERGENCE.tsv")
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.subplots_adjust(left=.07, right=.96, bottom=.16, top=.88, hspace=.55, wspace=.30)
    labels = ["Original full", "Paired original", "Paired corrected"]
    x = np.arange(3)
    ax = axes[0, 0]
    for offset, rep, name, shade in [(-.18, "unintegrated", "Unintegrated PCs", "#4c83b5"), (.18, "harmony", "Harmony PCs", "#d18d44")]:
        values = [float(latent.loc[latent.branch.eq(b) & latent.representation.eq(rep), "euclidean_knn20_mean_overlap"].iloc[0]) for b in BRANCHES]
        ax.bar(x + offset, values, .34, label=name, color=shade)
        for p, value in zip(x + offset, values): ax.text(p, value + .013, f"{value:.3f}", ha="center", fontsize=9)
    ax.set(ylim=(0, 1.25), ylabel="Mean shared 20-neighbor fraction", title="1 · PCA neighborhoods are closer before integration")
    ax.set_xticks(x, labels, rotation=12); ax.legend(frameon=False, fontsize=9, loc="upper center", ncol=2)
    ax = axes[0, 1]
    z = graph[graph.representation.eq("harmony")].set_index("branch").loc[BRANCHES]
    ax.bar(x - .18, z.scanpy_mean_degree, .34, color="#4c83b5", label="Stored Scanpy connectivity")
    ax.bar(x + .18, z.seurat_mean_degree, .34, color="#d18d44", label="Stored Seurat SNN")
    for i, b in enumerate(BRANCHES): ax.text(i, float(z.loc[b, "seurat_mean_degree"]) + 3, f"edge Jaccard {z.loc[b, 'directed_stored_edge_jaccard']:.3f}", ha="center", fontsize=8)
    ax.set(ylim=(0, 87), ylabel="Mean stored neighbors / nucleus", title="2 · Neighborhood graph definitions differ")
    ax.set_xticks(x, labels, rotation=12); ax.legend(frameon=False, fontsize=9, loc="upper center", ncol=2)
    ax = axes[1, 0]
    original = controls[controls.branch.eq("original_full") & controls.representation.eq("harmony")]
    for graph_name, runner, label, shade in [("scanpy_connectivity", "scanpy_igraph", "Scanpy graph / runner", "#4c83b5"),
                                            ("seurat_SNN", "native_seurat_igraph_GroupSingletons", "Native Seurat SNN / runner", "#d18d44")]:
        q = original[original.graph.eq(graph_name) & original.runner.eq(runner)].sort_values("n_iterations")
        ax.plot(q.n_iterations, q.ARI_vs_original_Seurat, "o-", color=shade, label=label)
        for _, row in q.iterrows(): ax.text(row.n_iterations, row.ARI_vs_original_Seurat - .04, f"{row.ARI_vs_original_Seurat:.3f}", ha="center", fontsize=8)
    ax.set(ylim=(.5, 1.06), xlabel="Leiden iterations on each fixed graph", ylabel="ARI versus original Seurat partition",
           title="3 · More iterations make a small original-CA1 change")
    ax.set_xticks([2, 10, 30]); ax.legend(frameon=False, fontsize=8, loc="lower right")
    ax = axes[1, 1]
    for offset, cap, shade in [(-.18, 10, "#4c83b5"), (.18, 50, "#d18d44")]:
        values = [int(convergence.loc[convergence.branch.eq(b) & convergence.max_iter_harmony.eq(cap), "completed_harmony_iterations"].iloc[0]) for b in BRANCHES]
        ax.bar(x + offset, values, .34, color=shade, label=f"Maximum {cap}")
        for p, v in zip(x + offset, values): ax.text(p, v + .18, str(v), ha="center", fontsize=9)
    ax.set(ylim=(0, 9), ylabel="Completed Python Harmony iterations", title="4 · CA1 controls converge before either cap")
    ax.set_xticks(x, labels, rotation=12); ax.legend(frameon=False, fontsize=9)
    fig.suptitle("CA1–ProS: Harmony spaces and neighbor graphs differ between tools", fontsize=15, fontweight="bold")
    fig.text(.05, .045, "Produced diagnostic controls only; main memberships and coordinates remain unchanged. Python cap-10/cap-50 CA1 coordinates agree exactly after saved float32 matching.\nDegree/edge overlap exclude self-loops; native clustering controls use each runner's recorded conversion conventions. Upstream Harmony, graph construction, backend and random streams differ.\nThese controls do not provide an exact causal decomposition or establish global Harmony convergence, physical layers or biological correctness.", fontsize=9)
    records = save_figure(fig, "CA1_technical_diagnostics", out)
    return records, {"convergence_iterations": convergence[["branch", "max_iter_harmony", "completed_harmony_iterations"]].to_dict("records"),
                     "source_receipt": str((root / "COMPLETE.json").resolve())}


def coordinates(evidence, cfg, tool, branch, population):
    root = Path(cfg["previous_run_root"]) / "results"
    if tool == "scanpy":
        md = evidence.table(root / "scanpy" / branch / population / "nucleus_metadata.tsv",
                            usecols=["nucleus_id", "sample", "condition", "UMAP_harmony_1", "UMAP_harmony_2"])
    else:
        aliases = {"Excitatory": "EX", "Inhibitory": "IN", "CA1_ProS": "CA1_PROS_CANDIDATE"}
        stage = root / "seurat" / branch / ("log_global" if population == "allcells" else "annotate_subsets")
        if population != "allcells": stage = stage / aliases[population]
        md = evidence.table(stage / "nucleus_metadata.tsv", usecols=["nucleus_id", "sample", "condition"])
        xy = evidence.table(stage / "umap_harmony.tsv")
        cols = [c for c in xy if c != "nucleus_id"]
        if len(cols) != 2 or xy.nucleus_id.duplicated().any() or set(xy.nucleus_id) != set(md.nucleus_id): raise ValueError("Supplement coordinate axis differs")
        md = md.merge(xy.rename(columns={cols[0]: "UMAP_harmony_1", cols[1]: "UMAP_harmony_2"}), on="nucleus_id", validate="one_to_one")
    if md.nucleus_id.duplicated().any() or not np.isfinite(md[["UMAP_harmony_1", "UMAP_harmony_2"]]).all().all(): raise ValueError("Invalid supplement coordinates")
    return md[["nucleus_id", "sample", "condition", "UMAP_harmony_1", "UMAP_harmony_2"]]


def program_figure(frame, tool, branch, population, field):
    labels = frame[field].astype(object).fillna("Unresolved expression program").astype(str)
    fig, ax = plt.subplots(figsize=(15, 8.5))
    fig.subplots_adjust(left=.04, right=.65, top=.84, bottom=.12)
    order = np.random.default_rng(20261014).permutation(len(frame))
    xy = frame[["UMAP_harmony_1", "UMAP_harmony_2"]].to_numpy()
    colors = np.array([color(x) for x in labels])
    ax.scatter(xy[order, 0], xy[order, 1], c=colors[order], s=2.2, alpha=.9, linewidths=0, rasterized=True)
    for category in sorted(labels.unique()):
        ax.scatter([], [], s=38, color=color(category), label=textwrap.fill(f"{category} (n={labels.eq(category).sum():,})", 56))
    ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1), frameon=False, fontsize=8, labelspacing=.8)
    for cluster, group in frame.groupby("cluster_id"):
        center = group[["UMAP_harmony_1", "UMAP_harmony_2"]].median().to_numpy()
        text = ax.text(*center, f"C{cluster}", fontsize=9, fontweight="bold", ha="center")
        text.set_path_effects([pe.withStroke(linewidth=2.5, foreground="white")])
    ax.set(xlabel="Existing Harmony UMAP 1", ylabel="Existing Harmony UMAP 2"); ax.set_xticks([]); ax.set_yticks([])
    fig.suptitle(f"{tool.capitalize()} · {POP_LABELS[population]} · {branch}\nReviewed expression programs on {len(frame):,} nuclei", fontsize=14, fontweight="bold")
    fig.text(.03, .035, "Color shows the final produced expression-program label. Original-full WT evidence defines program calls; paired labels are carried unchanged.\nThese are expression-associated programs, not physical superficial/deep positions or interneuron morphology. Cluster IDs refer to frozen original Harmony 0.4 memberships.\nAll displayed nuclei use existing coordinates; no embedding, nucleus score or RNA-count aggregate was fitted or recomputed for this figure.", fontsize=9)
    return fig


def make_programs(evidence, cfg, nuclei, out):
    records = []
    pages = out / "programs"
    pages.mkdir()
    for population, field in [("CA1_ProS", "layer_program_label"), ("Inhibitory", "expression_program_label")]:
        with PdfPages(pages / (population + "_expression_programs_all_tools_branches.pdf")) as book:
            for branch in BRANCHES:
                for tool in ["scanpy", "seurat"]:
                    ann = nuclei[nuclei.tool.eq(tool) & nuclei.branch.eq(branch) & nuclei.population.eq(population) & nuclei.cluster_scope.eq("original_frozen")].copy()
                    md = coordinates(evidence, cfg, tool, branch, population)
                    if ann.nucleus_id.duplicated().any() or set(ann.nucleus_id) != set(md.nucleus_id): raise ValueError("Program annotation axis differs from source UMAP")
                    actual = ann.merge(md, on="nucleus_id", suffixes=("", "_coordinate"), validate="one_to_one")
                    for c in ["sample", "condition"]:
                        if not actual[c].eq(actual[c + "_coordinate"]).all(): raise ValueError("Program supplement sample/genotype changed")
                    actual = actual.loc[boolean(actual.keep_for_analysis)].copy()
                    fig = program_figure(actual, tool, branch, population, field)
                    book.savefig(fig, dpi=200)
                    stem = f"{tool}__{branch}__{population}__expression_program"
                    paths = save_figure(fig, stem, pages)
                    actual[["nucleus_id", "sample", "condition", "cluster_id", "reviewed_label", "label_confidence", field,
                            "UMAP_harmony_1", "UMAP_harmony_2"]].to_csv(pages / (stem + "_exact_metadata.tsv.gz"), sep="\t", index=False)
                    records.append(dict(tool=tool, branch=branch, population=population, program_field=field, **paths))
    return records


def make_crosswalks(nuclei, out):
    folder = out / "crosswalks"; folder.mkdir()
    records = []
    with PdfPages(folder / "CA1_tool_cluster_crosswalks.pdf") as book:
        for scope in ["original_frozen", "fresh_branch"]:
            for branch in BRANCHES:
                if scope == "fresh_branch" and branch == "original_full": continue
                z = nuclei[nuclei.population.eq("CA1_ProS") & nuclei.branch.eq(branch) & nuclei.cluster_scope.eq(scope)].copy()
                z = z.loc[boolean(z.keep_for_plot)]
                left = z[z.tool.eq("scanpy")][["nucleus_id", "cluster_id"]]
                right = z[z.tool.eq("seurat")][["nucleus_id", "cluster_id"]]
                if set(left.nucleus_id) != set(right.nucleus_id): raise ValueError("CA1 crosswalk lacks exact shared displayed nuclei")
                joined = left.merge(right, on="nucleus_id", validate="one_to_one", suffixes=("_scanpy", "_seurat"))
                counts = pd.crosstab(joined.cluster_id_scanpy, joined.cluster_id_seurat)
                order = lambda label: (0, int(label)) if str(label).isdigit() else (1, str(label))
                counts = counts.reindex(index=sorted(counts.index, key=order), columns=sorted(counts.columns, key=order))
                fraction = counts.div(counts.sum(axis=1), axis=0)
                fig, ax = plt.subplots(figsize=(12, 8))
                fig.subplots_adjust(left=.18, right=.90, bottom=.19, top=.84)
                im = ax.imshow(fraction.to_numpy(), vmin=0, vmax=1, cmap="Blues", aspect="auto")
                ax.set_xticks(range(len(counts.columns)), [f"C{x}" for x in counts.columns])
                ax.set_yticks(range(len(counts.index)), [f"C{x} · n={counts.loc[x].sum():,}" for x in counts.index])
                for i in range(len(counts.index)):
                    for j in range(len(counts.columns)):
                        q = fraction.iloc[i, j]
                        if q >= .04: ax.text(j, i, f"{q:.0%}", ha="center", va="center", color="white" if q > .6 else "#1d3e57", fontsize=9)
                ax.set(xlabel="Seurat cluster", ylabel="Scanpy cluster")
                fig.colorbar(im, ax=ax, fraction=.035, label="Fraction of this Scanpy group")
                fig.suptitle(f"CA1-ProS · {branch}\n{'Frozen original memberships' if scope == 'original_frozen' else 'Secondary fresh paired partitions'} · {len(joined):,} shared nuclei", fontsize=13, fontweight="bold")
                fig.text(.06, .035, "Rows describe where Scanpy group's nuclei occur in Seurat groups; numbers do not assert equivalent identities.\nFrozen paired membership maps are held fixed across raw/corrected counts. Fresh views are a separate graph sensitivity.\nThis is tabulation of produced metadata, not RNA-count aggregation, a statistical test or anatomical validation.", fontsize=9)
                stem = f"CA1__{branch}__{scope}__crosswalk"
                book.savefig(fig); paths = save_figure(fig, stem, folder)
                counts.to_csv(folder / (stem + "_nucleus_counts.tsv"), sep="\t")
                fraction.to_csv(folder / (stem + "_row_fractions.tsv"), sep="\t")
                records.append(dict(branch=branch, cluster_scope=scope, **paths))
    return records


def make_pipeline(retention, out):
    main = retention[retention.population.eq("allcells") & retention.tool.eq("scanpy") & retention.cluster_scope.eq("original_frozen")].set_index("branch")
    if main.index.duplicated().any() or set(main.index) != set(BRANCHES): raise ValueError("Pipeline retention scope is ambiguous")
    if main.loc["paired_original", "kept_nuclei"] != main.loc["paired_cellbender", "kept_nuclei"]:
        raise ValueError("Pipeline paired cohorts differ")
    fig, ax = plt.subplots(figsize=(14, 6)); ax.axis("off")
    labels = [f"Accepted original source\n{int(main.loc['original_full','n_nuclei']):,} nuclei", 
              f"Revised shared All mask\n{int(main.loc['original_full','kept_nuclei']):,} nuclei", 
              f"Paired original\n{int(main.loc['paired_original','kept_nuclei']):,} shared retained nuclei", 
              f"Paired CellBender\n{int(main.loc['paired_cellbender','kept_nuclei']):,} shared retained nuclei"]
    positions = [(.07, .65), (.34, .65), (.69, .75), (.69, .35)]
    for i, (label, pos) in enumerate(zip(labels, positions)):
        ax.text(*pos, label, transform=ax.transAxes, fontsize=13, ha="left", va="center",
                bbox=dict(boxstyle="round,pad=.8", facecolor=["#eaf0f5", "#eaf0f5", "#e5f0ec", "#fbefdf"][i], edgecolor="#c5d0d7"))
    for a, b in [((.26, .65), (.32, .65)), ((.56, .65), (.67, .75)), ((.56, .65), (.67, .35))]:
        ax.annotate("", xy=b, xytext=a, xycoords="axes fraction", arrowprops=dict(arrowstyle="->", color="#586f7f", lw=2))
    ax.text(.32, .40, "All/EX: iterative shared exclusion of\noriginal groups below 100 nuclei.\nCA1-ProS / IN: no 100-nucleus exclusion.", transform=ax.transAxes, fontsize=11)
    ax.text(.68, .10, "Compare these two branches to isolate\ncount source on identical frozen memberships.", transform=ax.transAxes, fontsize=11)
    fig.suptitle("How the three count branches relate", fontsize=16, fontweight="bold")
    fig.text(.05, .03, "All-population nucleus totals come from the produced retention audit. Original-full versus paired also changes the cohort; paired raw/corrected masks are identical.\nExisting QC and singlet calls are reused. Marker-confidence mouse coverage and genotype-model eligibility are separate rules.", fontsize=9)
    return save_figure(fig, "count_branches_and_filtering", out)


def programme_agreement(evidence, cfg, out):
    root = Path(cfg["study_root"]) / "reviews/CA1_programme_agreement_v2_0_0"
    receipt = evidence.stage(root, "COMPLETED_READ_ONLY_CA1_PROGRAMME_SENSITIVITY", "REVIEW_PROVENANCE.json")
    if receipt["primary_DE_memberships_modified"] or receipt["scientific_fits_performed"] or receipt["raw_RNA_aggregation_performed"]:
        raise ValueError("Programme supplement scope differs from produced metadata review")
    if receipt["actual_final_annotation_receipt_sha256"] != sha(Path(cfg["final_annotation_dir"]) / "COMPLETE.json"):
        raise ValueError("Programme review uses a different final annotation")
    semantic = evidence.table(root / "programme_semantic_agreement.tsv")
    overlap = evidence.table(root / "same_named_programme_overlap.tsv")
    original = semantic[semantic.branch.eq("original_full")].set_index("leading_WT_score_threshold")
    primary = overlap[overlap.branch.eq("original_full") & overlap.leading_WT_score_threshold.eq(.3)].set_index("programme")
    deep, superficial = primary.loc["Deep"], primary.loc["Superficial"]
    narrative = (f"At the primary WT-score threshold 0.3, same-name programme agreement is {original.loc[.3,'same_named_programme_fraction']:.1%}. "
                 f"Deep-associated groups have Jaccard {deep.same_named_Jaccard:.3f} and reciprocal overlaps {deep.Scanpy_overlap_fraction:.1%}/{deep.Seurat_overlap_fraction:.1%}, meeting 50% in both directions. "
                 f"Superficial-associated groups have Jaccard {superficial.same_named_Jaccard:.3f} and overlaps {superficial.Scanpy_overlap_fraction:.1%}/{superficial.Seurat_overlap_fraction:.1%}; they fail the reciprocal criterion. "
                 f"Thresholds 0.2/0.3/0.4 give {original.loc[.2,'same_named_programme_fraction']:.1%}/{original.loc[.3,'same_named_programme_fraction']:.1%}/{original.loc[.4,'same_named_programme_fraction']:.1%} same-name agreement. "
                 "These diagnostic recodings do not alter primary memberships. They do not establish equivalent superficial groups or physical layers.")
    pdf = root / "CA1_programme_agreement_and_threshold_sensitivity.pdf"
    cards = []
    for page, title in [(1, "Same-name programme overlap at the primary threshold"), (2, "Programme threshold sensitivity")]:
        png = evidence.bind(root / f"page-{page}.png")
        cards.append(dict(title=title, pdf=local_link(evidence.bind(pdf),out), png=local_link(png,out), thumbnail=thumbnail(pdf,page,out,evidence),caption="Produced metadata review; primary memberships unchanged and physical layers undetermined."))
    return dict(narrative=narrative, cards=cards, review=local_link(evidence.bind(root / "review.md"),out),
                semantic_table=local_link(root / "programme_semantic_agreement.tsv",out), overlap_table=local_link(root / "same_named_programme_overlap.tsv",out))


def write_html(out, data):
    payload = json.dumps(data, ensure_ascii=False).replace("<", "\\u003c").replace(">", "\\u003e").replace("&", "\\u0026")
    template = r'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>SATB2 · Visual results guide</title>
<style>:root{--ink:#203746;--muted:#607483;--blue:#3975ad;--bg:#f5f7f8;--line:#dbe3e8}*{box-sizing:border-box}body{margin:0;color:var(--ink);background:var(--bg);font:15px/1.5 system-ui,-apple-system,sans-serif}header{padding:42px max(24px,calc((100vw - 1250px)/2));background:#173446;color:white}h1{font-size:34px;margin:8px 0}h2{font-size:25px;margin-top:0}h3{margin:0 0 8px}p{margin:8px 0}header p{max-width:900px;color:#d3e0e8}.badge{display:inline-block;border-radius:16px;padding:4px 12px;background:#315667;color:#e4f1f7;font-size:12px}.wrap{max-width:1250px;margin:auto;padding:26px 22px}.toolbar{position:sticky;top:0;z-index:5;background:#fff;border:1px solid var(--line);border-radius:12px;padding:14px;display:flex;gap:14px;flex-wrap:wrap;box-shadow:0 6px 20px #1734460c}label{font-weight:600;font-size:12px;display:block}select{font:inherit;padding:8px 10px;border:1px solid #bdcbd5;border-radius:7px;min-width:160px}.note{padding:15px 18px;background:#eaf1f6;border-left:4px solid #628cab;border-radius:7px;margin:18px 0}.grid{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:22px}.card{background:#fff;border:1px solid var(--line);border-radius:13px;overflow:hidden;padding:20px;box-shadow:0 3px 14px #17344608}.card img{width:100%;height:270px;object-fit:contain;background:white;cursor:zoom-in}.links{display:flex;gap:8px;flex-wrap:wrap;margin:14px 0}a{color:#276c9b}a.button{display:inline-block;background:#e8f0f5;border:1px solid #d2e1eb;text-decoration:none;border-radius:7px;padding:7px 11px;font-weight:600;font-size:13px}.count{color:var(--muted);font-size:13px}.muted{color:var(--muted)}section{margin:28px 0}.hero-grid{display:grid;grid-template-columns:1fr 1fr;gap:20px}.hero-grid img{width:100%;object-fit:contain}.hero-grid .card{padding:16px}details{border-top:1px solid var(--line);padding:12px 0}summary{cursor:pointer;font-weight:650}.pages{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:9px}.pages a{display:block;font-size:13px}.tablebox{overflow:auto;border:1px solid var(--line);border-radius:8px;max-height:430px}table{border-collapse:collapse;width:100%;background:white;font-size:13px}th,td{padding:9px 12px;border-bottom:1px solid var(--line);text-align:left;vertical-align:top}th{background:#eaf0f4;position:sticky;top:0}.scope{font-size:11px;text-transform:uppercase;letter-spacing:.06em;color:#627887;margin:5px 0 12px}.program-grid{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:18px}.program-grid img{width:100%;object-fit:contain}.legend{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:14px}.legend div{background:white;border:1px solid var(--line);border-radius:9px;padding:14px}.anchor-links{display:flex;gap:18px;margin-top:18px}.anchor-links a{color:#d8e9f4;text-decoration:none}footer{font-size:12px;color:var(--muted);margin:28px 0}.empty{padding:25px;background:white;border-radius:10px}@media(max-width:850px){.grid,.hero-grid,.program-grid,.legend{grid-template-columns:1fr}.toolbar{position:static}.card img{height:auto}}@media print{.toolbar,.links{display:none}header{background:white;color:#203746}.card{break-inside:avoid}.grid{grid-template-columns:1fr}}
</style></head><body><header><span class="badge">12 independent mice · 6 WT + 6 cKO</span><h1>SATB2 · Visual results guide</h1><p>Start with original-full counts and compare both tools. Read the molecular labels together with marker expression and uncertainty, then open the paired count branches to check sensitivity.</p><nav class="anchor-links"><a href="#explore">Annotated UMAPs & markers</a><a href="#programs">CA1 / IN programs</a><a href="#technical">CA1 technical checks</a><a href="#effects">Genotype effect figures</a></nav></header><main class="wrap">
<section class="legend"><div><h3>Main population scope</h3><p>Frozen original Harmony 0.4 memberships with the reviewed shared analysis mask. Original-full is the default. A cluster number is local to its tool and declared scope.</p></div><div><h3>Confidence is evidence review</h3><p>Reference agreement and full-cohort markers support provisional molecular labels. These review tiers are not calibrated probabilities, anatomical layers or interneuron morphology.</p></div><div><h3>Read dotplots correctly</h3><p>Color: gene-scaled expression. Size: detected nuclei, with extra downweight below zero. Negative z dots are smaller and fainter; exact z and detection percentages remain in TSVs.</p></div></section>
<section id="pipeline" class="card"><h2>From accepted nuclei to matched count branches</h2><img id="pipeline-image" style="width:100%;max-height:370px;object-fit:contain" alt="Produced-data count branch and shared filtering diagram"><div class="links" id="pipeline-links"></div><details><summary>Show source cohort / revised mask counts</summary><div id="retention-table" class="tablebox"></div></details></section>
<section id="explore"><h2>Explore annotated UMAPs and marker evidence</h2><div class="toolbar"><div><label for="population">Population</label><select id="population"><option value="all">All four populations</option><option value="allcells">All nuclei</option><option value="Excitatory">Excitatory neurons</option><option value="Inhibitory">Inhibitory neurons</option><option value="CA1_ProS">CA1-ProS candidates</option></select></div><div><label for="branch">Count branch</label><select id="branch"><option value="original_full">Original full · start here</option><option value="paired_original">Paired original · same-profile control</option><option value="paired_cellbender">Paired CellBender · count sensitivity</option></select></div><div><label for="tool">Tool</label><select id="tool"><option value="all">Both tools</option><option value="scanpy">Scanpy</option><option value="seurat">Seurat v5</option></select></div><div><label for="scope">View scope</label><select id="scope"><option value="original_frozen">Main frozen population labels</option><option value="fresh_branch">Secondary fresh paired partitions</option></select></div></div><p class="note" id="scope-note"></p><div class="grid" id="books"></div></section>
<section id="programs"><h2>CA1 and inhibitory expression programs</h2><p class="note">These supplemental UMAPs color the final produced expression-program labels. Original-full WT evidence defines the calls; paired labels are carried unchanged. A superficial/deep-associated expression program is not a physical layer assignment. Mixed and unresolved labels stay visible.</p><div class="links" id="program-links"></div><div class="program-grid" id="program-cards"></div></section>
<section id="programme-agreement"><h2>CA1 programme agreement needs a separate check</h2><p id="agreement-note" class="note"></p><div id="agreement-cards" class="program-grid"></div><div id="agreement-links" class="links"></div></section>
<section id="crosswalks"><h2>CA1 cluster correspondence between tools</h2><p class="muted">A row shows where one Scanpy group's nuclei occur in Seurat groups. A high overlap does not establish the same anatomical identity; paired frozen maps are held fixed by design.</p><div class="links" id="crosswalk-links"></div><div class="program-grid" id="crosswalk-cards"></div></section>
<section id="technical" class="card"><h2>Why fine CA1 partitions differ</h2><img id="technical-image" style="width:100%;max-height:790px;object-fit:contain" alt="CA1 produced technical diagnostics"><p>Unintegrated PC neighborhoods are nearly identical. Harmony spaces and graph definitions differ more. Increasing Leiden iterations produces only a small original-CA1 agreement change. The reproduced CA1 Python Harmony controls converge before both tested caps.</p><p class="muted">These controls include upstream Harmony, graph construction, graph-to-runner conversion, backend and random-stream differences. They do not give an exact causal decomposition, validate anatomy or establish global convergence. Main memberships are unchanged.</p><div class="links" id="technical-links"></div></section>
<section id="effects"><h2>Mouse-level genotype effects</h2><p class="note">Cells are not independent genotype replicates. Both population definitions use the same R DESeq2 engine; Python independently validates count transformation. Unavailable FDR and untested populations stay distinct from nonsignificant results.</p><div class="program-grid" id="effect-cards"></div></section>
<section><h2>Current cluster labels and uncertainty</h2><p class="muted">The table follows the population, branch, tool and scope controls above. Downloadable exact source labels remain in the annotation table.</p><div id="annotation-table" class="tablebox"></div><div class="links" id="annotation-links"></div></section>
<footer id="provenance"></footer></main><script>
const DATA=__DATA__;
const el=id=>document.getElementById(id),make=(tag,text,cls)=>{const e=document.createElement(tag);if(text!==undefined)e.textContent=text;if(cls)e.className=cls;return e};
function link(text,path,download=false){const a=make('a',text,'button');a.href=path;if(download)a.setAttribute('download','');else{a.target='_blank';a.rel='noopener'}return a}
function table(rows,columns){const t=make('table'),head=make('thead'),hr=make('tr');columns.forEach(c=>hr.append(make('th',c[1])));head.append(hr);t.append(head);const b=make('tbody');rows.forEach(r=>{const tr=make('tr');columns.forEach(c=>tr.append(make('td',r[c[0]]===null||r[c[0]]===undefined?'—':String(r[c[0]]))));b.append(tr)});t.append(b);return t}
function imageCard(item,title){const d=make('div',undefined,'card');d.append(make('h3',title));const a=make('a'),im=make('img');a.href=item.png;a.target='_blank';im.src=item.thumbnail;im.alt=title;a.append(im);d.append(a);if(item.caption)d.append(make('p',item.caption,'muted'));const ls=make('div',undefined,'links');ls.append(link('Download PDF',item.pdf,true),link('Open full PNG',item.png));if(item.svg)ls.append(link('SVG',item.svg));d.append(ls);return d}
function render(){let branch=el('branch').value,scope=el('scope').value;if(scope==='fresh_branch'&&branch==='original_full'){branch='paired_original';el('branch').value=branch}const pop=el('population').value,tool=el('tool').value;const fits=x=>(pop==='all'||x.population===pop)&&(tool==='all'||x.tool===tool)&&x.branch===branch&&x.cluster_scope===scope;
 el('scope-note').textContent=scope==='original_frozen'?'Main analysis: reviewed original memberships, held fixed across paired raw/corrected counts. Existing UMAP coordinates are not refit after the shared mask. Open the complete book for markers, original-query MapMyCells labels, genotype views and all twelve mice.':'Secondary view: fresh paired graph labels can change membership. These pages are representation sensitivity only and do not borrow frozen-population dotplots or genotype-inference labels.';
 el('books').replaceChildren();DATA.books.filter(fits).forEach(x=>{const c=make('div',undefined,'card');c.append(make('h3',x.title));c.append(make('p',x.population_label+' · '+Number(x.displayed_nuclei).toLocaleString()+' displayed nuclei','count'));let selected=x.preview_png;const im=make('img');im.src=x.previews[0].thumbnail;im.alt=x.title+' reviewed UMAP';im.onclick=()=>window.open(selected,'_blank');c.append(im);const tabs=make('div',undefined,'links');x.previews.forEach(p=>{const b=make('button',p.label);b.type='button';b.style.cssText='border:1px solid #bdd0dd;background:#f1f6f9;border-radius:6px;padding:7px 10px;color:#265873;cursor:pointer';b.onclick=()=>{im.src=p.thumbnail;im.alt=x.title+' '+p.label;selected=p.png};tabs.append(b)});c.append(tabs);const ls=make('div',undefined,'links');ls.append(link('Download complete PDF',x.pdf,true),link('Full annotated UMAP',x.preview_png),link('Cluster IDs',x.cluster_png),link('All 12 mice',x.mouse_png));c.append(ls);const det=make('details');det.append(make('summary','Browse every figure in this book'));const pages=make('div',undefined,'pages');x.pages.forEach(p=>pages.append(link(p.label,p.png)));det.append(pages);c.append(det);el('books').append(c)});
 const rows=DATA.annotations.filter(fits);el('annotation-table').replaceChildren(table(rows,[['cluster_id','Cluster'],['n_nuclei','Retained nuclei'],['keep_for_analysis','In analysis'],['reviewed_label','Reviewed molecular label'],['confidence_display','Coarse identity review'],['expression_program_label','Expression program'],['program_evidence_display','Program evidence scope'],['uncertainty_reason','Uncertainty']]));
 el('program-cards').replaceChildren();if(scope==='fresh_branch')el('program-cards').append(make('div','Program maps use frozen original memberships. Choose Main frozen population labels to view them.','empty'));DATA.programs.filter(x=>scope==='original_frozen'&&x.branch===branch&&(pop==='all'||x.population===pop)&&(tool==='all'||x.tool===tool)).forEach(x=>el('program-cards').append(imageCard(x,x.title)));
 el('crosswalk-cards').replaceChildren();DATA.crosswalks.filter(x=>x.branch===branch&&x.cluster_scope===scope).forEach(x=>el('crosswalk-cards').append(imageCard(x,'CA1-ProS · '+branch+' · '+scope)));
}
['population','branch','tool','scope'].forEach(id=>el(id).addEventListener('change',render));
el('agreement-note').textContent=DATA.agreement.narrative;DATA.agreement.cards.forEach(x=>el('agreement-cards').append(imageCard(x,x.title)));el('agreement-links').append(link('Programme review',DATA.agreement.review),link('Same-name agreement table',DATA.agreement.semantic_table),link('Reciprocal overlap table',DATA.agreement.overlap_table));
el('pipeline-image').src=DATA.pipeline.thumbnail;el('technical-image').src=DATA.technical.thumbnail;for(const [id,item] of [['pipeline-links',DATA.pipeline],['technical-links',DATA.technical]])el(id).append(link('Download PDF',item.pdf,true),link('Open full PNG',item.png),link('SVG',item.svg));
el('retention-table').append(table(DATA.retention,[['tool','Tool'],['branch','Branch'],['population','Population'],['cluster_scope','Scope'],['n_nuclei','Source nuclei'],['kept_nuclei','Revised kept'],['excluded_nuclei','Explicitly excluded']]));
if(DATA.effects.length)DATA.effects.forEach(x=>el('effect-cards').append(imageCard(x,x.title)));else el('effect-cards').append(make('div','The final genotype-effect figure supplement was not yet published when this guide was built. Available annotated plots do not establish genotype effects.','empty'));
DATA.program_books.forEach(x=>el('program-links').append(link(x.title,x.pdf,true)));el('crosswalk-links').append(link('Download all CA1 crosswalks PDF',DATA.crosswalk_book,true));
el('annotation-links').append(link('Exact cluster annotations TSV',DATA.annotation_tsv),link('Exact nucleus annotations TSV.gz',DATA.nucleus_tsv));
el('provenance').textContent=DATA.provenance;render();
</script></body></html>'''
    (out / "visual_guide.html").write_text(template.replace("__DATA__", payload))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", required=True, type=Path)
    ap.add_argument("--outdir", required=True, type=Path)
    ap.add_argument("--effect-index", type=Path)
    ap.add_argument("--diagnostics-only", action="store_true")
    args = ap.parse_args()
    evidence = Evidence(); config_path = evidence.bind(args.config); cfg = json.loads(config_path.read_text())
    out = args.outdir.resolve()
    if out.exists() and any(out.iterdir()): raise ValueError("Refusing to overwrite an existing guide directory")
    if not args.diagnostics_only:
        evidence.stage(Path(cfg["final_annotation_dir"]), "ANNOTATION_FINAL_COMPLETE_REVIEW_REQUIRED")
        evidence.stage(Path(cfg["results_root"]) / "plots", "COMPLETE")
        manifest = evidence.table(Path(cfg["results_root"]) / "inputs/INPUT_SHA256.tsv")
        if manifest.path.duplicated().any(): raise ValueError("Frozen input manifest has duplicate paths")
        for row in manifest.itertuples():
            evidence.expected_outputs[str(Path(row.path).resolve())] = row.sha256
    out.mkdir(parents=True, exist_ok=True)
    technical, _ = diagnostic_figure(evidence, cfg, out)
    if not args.diagnostics_only:
        root = Path(cfg["results_root"]); ad = Path(cfg["final_annotation_dir"])
        ann = evidence.table(ad / "cluster_annotations.tsv")
        ann["confidence_display"] = ann.label_confidence.map(CONFIDENCE_LABELS).fillna(ann.label_confidence)
        ann["program_evidence_display"] = ann.program_evidence_status.str.replace("_", " ").str.lower()
        nuclei = evidence.table(ad / "nucleus_annotations.tsv.gz", usecols=KEYS + ["nucleus_id", "sample", "condition", "keep_for_analysis", "keep_for_plot",
                                                                                "reviewed_label", "label_confidence", "layer_program_label", "expression_program_label"])
        retention = evidence.table(ad / "annotation_retention.tsv")
        plot_books = evidence.table(root / "plots/PLOT_BOOK_INDEX.tsv")
        expected = {(t, b, p, "original_frozen") for t in cfg["tools"] for b in BRANCHES for p in cfg["populations"]}
        expected |= {(t, b, p, "fresh_branch") for t in cfg["tools"] for b in BRANCHES[1:] for p in cfg["populations"]}
        if plot_books.duplicated(KEYS[:-1]).any() or set(map(tuple, plot_books[KEYS[:-1]].to_numpy())) != expected:
            raise ValueError("Published plot books do not cover the declared main/secondary datasets")
        agreement = programme_agreement(evidence, cfg, out)
        programs = make_programs(evidence, cfg, nuclei, out)
        crosswalks = make_crosswalks(nuclei, out)
        pipeline = make_pipeline(retention, out)
        data = dict(books=[], programs=[], crosswalks=[], effects=[], agreement=agreement,
                    annotations=json.loads(ann.where(pd.notna(ann), None).to_json(orient="records")),
                    retention=json.loads(retention.to_json(orient="records")),
                    annotation_tsv=local_link(ad / "cluster_annotations.tsv", out),
                    nucleus_tsv=local_link(ad / "nucleus_annotations.tsv.gz", out))
        data["program_books"] = [dict(title="Download " + POP_LABELS[p] + " program book PDF", pdf=local_link(evidence.bind(out / "programs" / (p + "_expression_programs_all_tools_branches.pdf")), out))
                                 for p in ["CA1_ProS", "Inhibitory"]]
        data["crosswalk_book"] = local_link(evidence.bind(out / "crosswalks/CA1_tool_cluster_crosswalks.pdf"), out)
        for _, row in plot_books.iterrows():
            pdf = evidence.bind(row.pdf); stage = pdf.parent
            pages = evidence.table(stage / "PAGE_INDEX.tsv")
            preview = Path(pages.loc[pages.name.eq("reviewed_molecular_labels"), "png"].iloc[0])
            cluster = Path(pages.loc[pages.name.eq("cluster_ids_and_reviewed_labels"), "png"].iloc[0])
            mouse = Path(pages.loc[pages.name.eq("all_twelve_mouse_facets"), "png"].iloc[0])
            evidence.bind(cluster); evidence.bind(mouse)
            previews = [dict(label="Reviewed UMAP", png=local_link(preview, out), thumbnail=thumbnail(pdf, 1, out, evidence))]
            for prefix, label in [("representative_markers_", "Representative genes"), ("supplied_panel_", "Broad markers")]:
                q = pages[pages.name.str.startswith(prefix)]
                if len(q):
                    p = Path(q.iloc[0].png)
                    previews.append(dict(label=label, png=local_link(p, out), thumbnail=thumbnail(pdf, int(q.iloc[0].page), out, evidence)))
            wanted = []
            if row.population in ["Excitatory", "CA1_ProS"]:
                wanted.append(("CA1_population_candidates", "CA1 marker panel · first page"))
            if row.population == "CA1_ProS":
                wanted.extend([("Superficial_associated_VEH_program", "Superficial-associated program"),
                               ("Deep_associated_VEH_program", "Deep-associated program · first page")])
            if row.population == "Inhibitory":
                wanted.append(("Pvalb_Basket_cells", "Pvalb / Sst programs"))
                wanted.append(("Cck_Non-border", "Cck programs · first page"))
            for gene_set, label in wanted:
                for pr in pages[pages.name.str.startswith("supplied_panel_")].itertuples():
                    membership = evidence.table(stage / (pr.name + "_membership.tsv"))
                    if gene_set in set(membership["set"]):
                        p = Path(pr.png)
                        previews.append(dict(label=label, png=local_link(p, out), thumbnail=thumbnail(pdf, int(pr.page), out, evidence)))
                        break
            data["books"].append(dict(tool=row.tool, branch=row.branch, population=row.population, cluster_scope=row.cluster_scope,
                                      population_label=POP_LABELS[row.population], displayed_nuclei=int(row.displayed_nuclei),
                                      title=row.tool.capitalize() + " · " + row.branch.replace("_", " "),
                                      preview_png=local_link(preview, out),
                                      previews=previews,
                                      cluster_png=local_link(cluster, out), mouse_png=local_link(mouse, out), pdf=local_link(pdf, out),
                                      pages=[dict(label=str(p.page) + " · " + p.name.replace("_", " "), png=local_link(evidence.bind(p.png), out)) for p in pages.itertuples()]))
        for group, source in [("programs", programs), ("crosswalks", crosswalks)]:
            for item in source:
                data[group].append({**{k:v for k,v in item.items() if k not in ["png", "pdf", "svg"]},
                                    **{k:local_link(item[k], out) for k in ["png", "pdf", "svg"]},
                                    "thumbnail": thumbnail(item["pdf"], 1, out, evidence),
                                    "title": item.get("tool", "").capitalize() + " · " + POP_LABELS.get(item.get("population", "CA1_ProS"), "CA1-ProS"),
                                    "caption": ("Produced expression labels on existing coordinates; physical position and morphology remain undetermined."
                                                if group == "programs" else "Shared-nucleus membership tabulation; row percentages show overlap rather than equivalent identities.")})
        for name, item in [("pipeline", pipeline), ("technical", technical)]:
            data[name] = {**{k:local_link(item[k], out) for k in ["png", "pdf", "svg"]}, "thumbnail": thumbnail(item["pdf"], 1, out, evidence)}
        effect_index = args.effect_index or Path(cfg["study_root"]) / "reports/effect_figures/FIGURE_INDEX.tsv"
        if effect_index.exists():
            evidence.stage(effect_index.parent, "COMPLETE")
            effects = evidence.table(effect_index)
            for row in effects.to_dict("records"):
                row = {**row, **{ext: row.get(ext + "_path", row.get(ext)) for ext in ["pdf", "svg", "png"]}}
                if not all(Path(str(row.get(k, ""))).is_file() for k in ["pdf", "png"]): raise ValueError("Effect supplement index points to missing output")
                record = dict(title=str(row.get("title", row.get("name", "Mouse-level effect figure"))), caption=str(row.get("caption", "")),
                              pdf=local_link(evidence.bind(row["pdf"]), out), png=local_link(evidence.bind(row["png"]), out), thumbnail=thumbnail(row["pdf"], 1, out, evidence))
                if row.get("svg") and Path(str(row["svg"])).is_file(): record["svg"] = local_link(evidence.bind(row["svg"]), out)
                data["effects"].append(record)
        comparison = root / "comparison/COMPLETE.json"
        if comparison.exists(): evidence.bind(comparison)
        data["provenance"] = "Produced-data rendering only. No raw RNA matrix loading/aggregation or scientific fitting. Main scope: reviewed frozen original labels; secondary fresh partitions are separate. Exact source/output hashes are in GUIDE_RECEIPT.json."
        write_html(out, data)
    evidence.verify()
    receipt = dict(status="RENDERED_PRODUCED_DATA_REVIEW_REQUIRED", source_code_sha256=sha(__file__), config_sha256=sha(config_path),
                   input_sha256=evidence.inputs, scientific_fitting_performed=False, raw_nucleus_RNA_aggregation_performed=False,
                   coordinate_refit_performed=False, physical_layer_assignment_performed=False,
                   output_sha256={str(p.relative_to(out)):sha(p) for p in sorted(out.rglob("*")) if p.is_file()})
    (out / "GUIDE_RECEIPT.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps({"status":receipt["status"], "outdir":str(out), "diagnostics_only":args.diagnostics_only}))


if __name__ == "__main__": main()
