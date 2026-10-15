#!/usr/bin/env python3
"""Render reviewed annotations and native aggregate summaries; no count loading."""
from __future__ import annotations

import argparse
import colorsys
import hashlib
import json
import os
import textwrap
import tempfile
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "satb2_revision_matplotlib_cache"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import Normalize
import matplotlib.patheffects as pe
import numpy as np
import pandas as pd

from plot_common import (KEYS, bool_values, canonical_keys, cluster_sort, empty_outdir,
                         read_table, sha256, write_receipt, write_table)

plt.rcParams.update({"pdf.fonttype": 42, "svg.fonttype": "none", "font.size": 9})
BROAD_PALETTE = {"Excitatory": "#dc8b25", "Inhibitory": "#427cb6", "Astrocyte": "#8473b6",
                 "Oligodendrocyte": "#788c32", "OPC": "#bb8053", "Microglia": "#ad5077",
                 "Endothelial": "#35a080", "Pericyte": "#a46d68", "Unknown": "#9b9b9b",
                 "WT": "#3976ba", "cKO": "#db892b"}


def category_color(label):
    """Same named molecular category gets the same color in every book/tool."""
    label = str(label)
    if label in BROAD_PALETTE:
        return BROAD_PALETTE[label]
    h = int(hashlib.sha256(label.encode()).hexdigest()[:8], 16) / 2**32
    return matplotlib.colors.to_hex(colorsys.hsv_to_rgb(h, .55, .70))


def load_coordinates(cfg, tool, branch, population):
    previous = Path(cfg["previous_run_root"])
    base = previous if previous.name == "results" else previous / "results"
    matches = [x for x in cfg.get("plot_sources", [])
               if x["tool"] == tool and x["branch"] == branch and x["population"] == population]
    if matches:
        specification = matches[0]
        path = Path(specification["metadata_tsv"])
        md = read_table(path)
        hp = Path(specification["harmony_coordinates_tsv"]) if specification.get("harmony_coordinates_tsv") else None
        up = Path(specification["unintegrated_coordinates_tsv"]) if specification.get("unintegrated_coordinates_tsv") else None
    elif tool == "scanpy":
        path = base / "scanpy" / branch / population / "nucleus_metadata.tsv"
        md = read_table(path)
        hp = up = None
    elif tool == "seurat":
        aliases = {"Excitatory": "EX", "Inhibitory": "IN", "CA1_ProS": "CA1_PROS_CANDIDATE"}
        d = base / "seurat" / branch / ("log_global" if population == "allcells" else "annotate_subsets")
        if population != "allcells":
            d = d / aliases[population]
        path = d / "nucleus_metadata.tsv"
        md = read_table(path)
        hp = d / "umap_harmony.tsv"
        up = d / "umap_unintegrated.tsv"
    else:
        raise ValueError(f"Unsupported tool: {tool}")
    sources = [path]
    if md.nucleus_id.duplicated().any():
        raise ValueError("Coordinate-source metadata contains duplicate nucleus IDs")
    for basis, coordinates in [("harmony", hp), ("unintegrated", up)]:
        if coordinates:
            c = read_table(coordinates)
            columns = [x for x in c if x != "nucleus_id"]
            if len(columns) != 2 or c.nucleus_id.duplicated().any() or set(c.nucleus_id) != set(md.nucleus_id):
                raise ValueError("Coordinate TSV does not preserve the exact source nucleus axis")
            c = c.rename(columns={columns[0]: "UMAP_" + basis + "_1", columns[1]: "UMAP_" + basis + "_2"})
            md = md.merge(c, on="nucleus_id", validate="one_to_one")
            sources.append(coordinates)
        required = ["UMAP_" + basis + "_1", "UMAP_" + basis + "_2"]
        if set(required).difference(md.columns) or not np.isfinite(md[required].to_numpy()).all():
            raise ValueError("Missing/nonfinite source UMAP coordinates")
    return md, sources


def join_annotations(metadata, annotations):
    if annotations.nucleus_id.duplicated().any() or set(annotations.nucleus_id) != set(metadata.nucleus_id):
        raise ValueError("Per-scope annotations must exhaust the exact source metadata roster")
    annotations = annotations.set_index("nucleus_id").loc[metadata.nucleus_id].reset_index()
    for column in ["sample", "condition"]:
        if column in annotations:
            if not annotations[column].astype(str).eq(metadata[column].astype(str).reset_index(drop=True)).all():
                raise ValueError(f"Annotation/coordinate {column} identities conflict")
            annotations = annotations.drop(columns=column)
    overlapping = [c for c in annotations if c in metadata and c != "nucleus_id"]
    metadata = metadata.drop(columns=overlapping)
    merged = metadata.merge(annotations, on="nucleus_id", validate="one_to_one")
    for col in ["keep_for_plot", "keep_for_analysis"]:
        merged[col] = bool_values(merged[col])
    return merged


def title_scope(tool, branch, population, scope, n):
    text = "Frozen original Harmony 0.4 labels" if scope == "original_frozen" else "Secondary fresh paired Harmony 0.4 labels"
    return f"{tool.capitalize()} · {population} · {branch}\n{text}; {n:,} displayed nuclei"


def footer(fig, caption):
    fig.text(.02, .012, caption, fontsize=7, va="bottom", ha="left", wrap=True)


def scatter(ax, frame, category, basis="harmony", label_clusters=False, palette=None, seed=20261014,
            legend=True, point_size=1.8):
    labels = frame[category].astype(object).fillna("Unknown").astype(str)
    coords = frame[[f"UMAP_{basis}_1", f"UMAP_{basis}_2"]].to_numpy()
    order = np.random.default_rng(seed).permutation(len(frame))
    colors = [palette.get(x, category_color(x)) if palette else category_color(x) for x in labels]
    ax.scatter(coords[order, 0], coords[order, 1], c=np.array(colors)[order], s=point_size,
               linewidths=0, alpha=.85, rasterized=True)
    if legend:
        for label in sorted(labels.unique(), key=cluster_sort):
            color = palette.get(label, category_color(label)) if palette else category_color(label)
            ax.scatter([], [], s=24, c=color, linewidths=0, label=label)
        if labels.nunique() <= 40:
            ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1), frameon=False, fontsize=7,
                      labelspacing=.45, handletextpad=.3)
    if label_clusters:
        for label, cluster in frame.groupby("cluster_id", observed=True):
            center = cluster[[f"UMAP_{basis}_1", f"UMAP_{basis}_2"]].median().to_numpy()
            text = ax.text(*center, f"C{label}", fontsize=8, fontweight="bold", ha="center", va="center")
            text.set_path_effects([pe.withStroke(linewidth=2.5, foreground="white")])
    ax.set(xlabel="UMAP 1", ylabel="UMAP 2")
    ax.set_xticks([]); ax.set_yticks([])


def molecular_figure(frame, stage_title, category="reviewed_label", cluster_ids=True):
    fig, ax = plt.subplots(figsize=(14.5, 8))
    fig.subplots_adjust(left=.05, right=.65, bottom=.09, top=.88)
    scatter(ax, frame, category, label_clusters=cluster_ids)
    ax.set_title(stage_title + "\n" + category.replace("_", " "), fontsize=12)
    footer(fig, "All displayed nuclei on existing v1 coordinates; the declared analysis/display mask is applied after the coordinate fit. Molecular labels remain provisional; no physical layers.")
    return fig


def paired_view(frame, stage_title):
    fig, axes = plt.subplots(1, 2, figsize=(14.5, 7))
    fig.subplots_adjust(left=.04, right=.95, wspace=.15, bottom=.11, top=.85)
    palette = {x: category_color(x) for x in frame.reviewed_label.unique()}
    for ax, basis in zip(axes, ["unintegrated", "harmony"]):
        scatter(ax, frame, "reviewed_label", basis=basis, label_clusters=True, palette=palette, legend=False)
        ax.set_title(basis.capitalize())
    fig.suptitle(stage_title, fontsize=12)
    footer(fig, "The same frozen/reviewed Harmony labels are drawn on both existing UMAPs. Coordinates were not refit; rotations or shapes are not numerical agreement tests.")
    return fig


def facet_figure(frame, stage_title):
    samples = sorted(frame["sample"].unique())
    fig, axes = plt.subplots(3, 4, figsize=(15, 11))
    fig.subplots_adjust(left=.04, right=.98, bottom=.07, top=.88, hspace=.22, wspace=.12)
    coords = frame[["UMAP_harmony_1", "UMAP_harmony_2"]].to_numpy()
    for ax, sample in zip(axes.ravel(), samples):
        ax.scatter(coords[:, 0], coords[:, 1], s=.3, c="#b9b9b9", alpha=.20, linewidths=0, rasterized=True)
        this = frame.loc[frame["sample"].eq(sample)]
        condition = this.condition.iloc[0]
        ax.scatter(this.UMAP_harmony_1, this.UMAP_harmony_2, s=2.2, c=category_color(condition),
                   alpha=.9, linewidths=0, rasterized=True)
        ax.set_title(f"{sample} · {condition} · n={len(this):,}", fontsize=9)
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_xlim(coords[:, 0].min() - .3, coords[:, 0].max() + .3)
        ax.set_ylim(coords[:, 1].min() - .3, coords[:, 1].max() + .3)
    fig.suptitle(stage_title + "\nAll twelve mice; gray context is the complete displayed cohort", fontsize=12)
    footer(fig, "Full-cohort coordinates with common axes within this page. WT blue; cKO orange. No nuclei are sampled for the facets.")
    return fig


def dot_visual_values(values, max_area=100, negative_multiplier=.25, negative_alpha=.18):
    """Glyph downweight is explicit; percent and z fields remain unchanged."""
    values = values.copy()
    negative = values.scaled_expression_z.lt(0)
    values["negative_z_visual_downweight_applied"] = negative
    values["glyph_area_points2"] = values.fraction_detected * max_area
    values.loc[negative, "glyph_area_points2"] *= negative_multiplier
    values["glyph_alpha"] = np.where(negative, negative_alpha, .9)
    values["glyph_color_z_clipped"] = values.scaled_expression_z.clip(-2.5, 2.5)
    return values


def dot_figure(summary, annotations, genes, panel_title, stage_title, cfg):
    p = cfg.get("plots", {})
    chosen = summary.loc[summary.gene_id.isin(genes)].copy()
    clusters = sorted(annotations.loc[annotations.keep_for_plot, "cluster_id"].astype(str).unique(), key=cluster_sort)
    chosen = chosen.loc[chosen.cluster_id.isin(clusters)]
    chosen = dot_visual_values(chosen, float(p.get("dot_max_area", 100)),
                               float(p.get("negative_z_area_multiplier", .25)),
                               float(p.get("negative_z_alpha", .18)))
    if chosen.duplicated(["cluster_id", "gene_id"]).any():
        raise ValueError("Duplicate gene/cluster dot summary")
    x_index, y_index = {g: i for i, g in enumerate(genes)}, {c: i for i, c in enumerate(clusters)}
    chosen["dot_x"] = chosen.gene_id.map(x_index)
    chosen["dot_y"] = chosen.cluster_id.map(y_index)
    height = max(6.5, len(clusters) * .34 + 3.4)
    fig, ax = plt.subplots(figsize=(max(12.5, len(genes) * .34 + 5), height))
    fig.subplots_adjust(left=.28, right=.84, bottom=1.7 / height, top=1 - 1.5 / height)
    norm = Normalize(-2.5, 2.5)
    cmap = plt.colormaps["RdBu_r"]
    colors = cmap(norm(chosen.glyph_color_z_clipped.fillna(0).to_numpy()))
    colors[:, 3] = chosen.glyph_alpha.to_numpy()
    valid = chosen.scaled_expression_z.notna() & chosen.fraction_detected.notna()
    ax.scatter(chosen.loc[valid, "dot_x"], chosen.loc[valid, "dot_y"],
               s=chosen.loc[valid, "glyph_area_points2"], c=colors[valid.to_numpy()], linewidths=0)
    symbols = chosen.drop_duplicates("gene_id").set_index("gene_id").gene_symbol.to_dict()
    ax.set_xticks(range(len(genes)), [symbols.get(g, g) for g in genes], rotation=70, ha="right", fontsize=8)
    labels = annotations.set_index("cluster_id").reviewed_label.astype(str).to_dict()
    ax.set_yticks(range(len(clusters)), [textwrap.fill(f"C{c} · {labels.get(c, 'Unresolved')}", 38) for c in clusters], fontsize=8)
    ax.set(xlim=(-.7, len(genes) - .3), ylim=(len(clusters) - .5, -.5))
    ax.grid(axis="y", alpha=.12)
    ax.set_axisbelow(True)
    panel_caption = textwrap.shorten(panel_title, width=130, placeholder=" …; exact panel groups in membership TSV")
    fig.suptitle(stage_title + "\n" + textwrap.fill(panel_caption, 80), fontsize=11, y=.975)
    color_ax = fig.add_axes([.87, .47, .015, .22])
    fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap), cax=color_ax, label="Gene z (display clipped ±2.5)")
    handles = [ax.scatter([], [], s=float(p.get("dot_max_area", 100)) * q,
                          c="#666666", linewidths=0, label=f"{100*q:g}%") for q in [.1, .25, .5, 1.0]]
    ax.legend(handles=handles, title="Detected nuclei\n(z ≥ 0 size)", loc="upper left", bbox_to_anchor=(1.04, .30),
              frameon=False, fontsize=8, title_fontsize=8)
    ref = ",".join(sorted(chosen.scale_reference_branch.dropna().astype(str).unique()))
    footer(fig, f"All native cohort nuclei; detection uses count >0. z<0: area ×{p.get('negative_z_area_multiplier', .25)}, alpha={p.get('negative_z_alpha', .18)}; exact percent/z unchanged in TSV.\nGene z reference: {ref}; equal cluster means, population SD. Programs describe expression, not anatomical layers.")
    return fig, chosen


def save_page(fig, pdf, stage_out, number, stem, manifests):
    name = f"{number:02d}_{stem}"
    pdf.savefig(fig, dpi=200)
    for extension in ["svg", "png"]:
        path = stage_out / "pages" / (name + "." + extension)
        fig.savefig(path, dpi=180, metadata={"Creator": "SATB2 reviewed rendering"} if extension == "svg" else None)
    manifests.append({"page": number, "name": stem, "svg": str((stage_out / "pages" / (name + ".svg")).resolve()),
                      "png": str((stage_out / "pages" / (name + ".png")).resolve())})
    plt.close(fig)


def panel_pages(panel, maximum_genes=28):
    """Keep every supplied membership row auditable, with concise grouped pages."""
    if "group" not in panel:
        panel["group"] = panel.get("marker_group", panel.get("set", "Provided_marker_panel"))
    groups = []
    for group, x in panel.loc[panel.gene_id.notna()].groupby("group", sort=False):
        ids = x.gene_id.astype(str).drop_duplicates().tolist()
        groups.append((str(group), ids))
    pages, labels, ids = [], [], []
    for label, genes in groups:
        additional = [g for g in genes if g not in ids]
        if ids and len(ids) + len(additional) > maximum_genes:
            pages.append(("; ".join(labels), ids)); labels, ids = [], []
        if len(genes) > maximum_genes:
            for start in range(0, len(genes), maximum_genes):
                pages.append((label, genes[start:start + maximum_genes]))
        else:
            labels.append(label); ids.extend(g for g in genes if g not in ids)
    if ids:
        pages.append(("; ".join(labels), ids))
    return pages


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", type=Path, required=True)
    ap.add_argument("--summaries", type=Path, required=True)
    ap.add_argument("--annotations", type=Path, required=True)
    ap.add_argument("--panel-memberships", type=Path)
    ap.add_argument("--outdir", type=Path, required=True)
    args = ap.parse_args()
    cfg = json.loads(args.config.read_text())
    receipt_path = args.summaries / "COMPLETE.json"
    summary_receipt = json.loads(receipt_path.read_text())
    if summary_receipt.get("status") != "COMPLETE":
        raise ValueError("Rendering requires completed native aggregate summaries")
    annotation_receipt_path = args.annotations / "COMPLETE.json"
    annotation_receipt = json.loads(annotation_receipt_path.read_text())
    if annotation_receipt.get("status") != "ANNOTATION_FINAL_COMPLETE_REVIEW_REQUIRED" or annotation_receipt.get("stage") != "final":
        raise ValueError("Main rendering requires the completed reviewed final annotation producer")
    paths = [args.annotations / "cluster_annotations.tsv", args.annotations / "nucleus_annotations.tsv.gz",
             args.summaries / "gene_cluster_summary.tsv", args.summaries / "representative_markers.tsv",
             args.summaries / "cluster_mouse_counts.tsv"]
    membership_path = args.panel_memberships or Path(cfg.get("marker_panel_memberships_tsv", args.annotations.parent / "annotation_initial/panel_membership.tsv"))
    for path in paths[:2]:
        if annotation_receipt.get("output_sha256", {}).get(path.name) != sha256(path):
            raise ValueError("Final annotation artifact checksum mismatch")
    for path in paths[2:]:
        if summary_receipt.get("output_sha256", {}).get(path.name) != sha256(path):
            raise ValueError("Native marker summary artifact checksum mismatch")
    paths += [membership_path, receipt_path, annotation_receipt_path]
    annotations = canonical_keys(read_table(paths[0]))
    nuclei = canonical_keys(read_table(paths[1]))
    summary = canonical_keys(read_table(paths[2]))
    representatives = canonical_keys(read_table(paths[3]))
    counts = canonical_keys(read_table(paths[4]))
    memberships = read_table(membership_path)
    for frame in [annotations, nuclei]:
        if "keep_for_plot" not in frame:
            frame["keep_for_plot"] = frame.keep_for_analysis
        frame["keep_for_plot"] = bool_values(frame.keep_for_plot)
        frame["keep_for_analysis"] = bool_values(frame.keep_for_analysis)
    if annotations.duplicated(KEYS).any():
        raise ValueError("Cluster annotation definitions are not unique")
    # Paired inference masks, including iterative two-tool exclusions, must be identical.
    for (tool, population, scope), x in nuclei.groupby(["tool", "population", "cluster_scope"], observed=True):
        if scope != "original_frozen":
            continue
        paired = [set(x.loc[x.branch.eq(b) & x.keep_for_analysis, "nucleus_id"]) for b in ["paired_original", "paired_cellbender"]]
        if paired[0] != paired[1]:
            raise ValueError("Paired-original/corrected reviewed analysis masks differ")
    for (branch, population), x in nuclei.loc[nuclei.cluster_scope.eq("original_frozen")].groupby(["branch", "population"], observed=True):
        if population in ["allcells", "Excitatory"]:
            sets = [set(z.loc[z.keep_for_analysis, "nucleus_id"]) for _, z in x.groupby("tool", observed=True)]
            if len(sets) != 2 or sets[0] != sets[1]:
                raise ValueError("Main All/EX plotting masks must be shared between both tools")
    out = empty_outdir(args.outdir)
    books, palettes = [], []
    for keys, cluster_annotations in annotations.groupby(KEYS[:-1], observed=True, sort=False):
        tool, branch, population, scope = keys
        if scope == "fresh_branch" and branch == "original_full":
            continue
        selected_annotations = nuclei
        for key, value in zip(KEYS[:-1], keys):
            selected_annotations = selected_annotations.loc[selected_annotations[key].eq(value)]
        metadata, source_paths = load_coordinates(cfg, tool, branch, population)
        paths += source_paths
        merged = join_annotations(metadata, selected_annotations)
        displayed = merged.loc[merged.keep_for_plot].copy()
        if not len(displayed):
            raise ValueError("A declared stage has no displayed nuclei")
        if scope == "original_frozen":
            if not merged.keep_for_plot.eq(merged.keep_for_analysis).all():
                raise ValueError("Main original_frozen plots must use the actual analysis mask")
            count_rows = counts
            for key, value in zip(KEYS[:-1], keys):
                count_rows = count_rows.loc[count_rows[key].eq(value)]
            actual = merged.loc[merged.keep_for_analysis].groupby(["cluster_id", "sample"]).size()
            expected = count_rows.set_index(["cluster_id", "sample"]).n_nuclei
            if not actual.reindex(expected.index, fill_value=0).equals(expected.astype(actual.dtype)):
                raise ValueError("Displayed inference cohort conflicts with native cluster/mouse denominators")
        stage_out = out / tool / branch / population / scope
        (stage_out / "pages").mkdir(parents=True)
        write_table(merged, stage_out / "rendered_nucleus_metadata.tsv.gz")
        write_table(cluster_annotations, stage_out / "cluster_annotation_scope.tsv")
        write_table(displayed.groupby(["cluster_id", "sample", "condition"]).size().rename("displayed_nuclei").reset_index(), stage_out / "displayed_cluster_counts_per_mouse.tsv")
        title = title_scope(*keys, len(displayed))
        pages = []
        with PdfPages(stage_out / "annotated_umaps_and_marker_dotplots.pdf") as pdf:
            def add(fig, stem):
                save_page(fig, pdf, stage_out, len(pages) + 1, stem, pages)
            displayed["cluster_display_label"] = displayed.apply(lambda x: f"C{x.cluster_id} · {x.reviewed_label}", axis=1)
            add(molecular_figure(displayed, title), "reviewed_molecular_labels")
            add(molecular_figure(displayed, title, "cluster_display_label"), "cluster_ids_and_reviewed_labels")
            add(paired_view(displayed, title), "unintegrated_and_harmony_same_reviewed_labels")
            broad = "reference_broad_prediction" if "reference_broad_prediction" in displayed else "mmc_reference_broad_prediction"
            add(molecular_figure(displayed, title + "\nMapMyCells predictions from verified original query", broad), "original_MapMyCells_broad")
            subclass = "reference_subclass_name" if "reference_subclass_name" in displayed else "mmc_subclass_name"
            if subclass in displayed:
                common = displayed[subclass].value_counts().head(30).index
                displayed["MapMyCells_subclass_display"] = displayed[subclass].where(displayed[subclass].isin(common), "Other original subclasses; full values in TSV")
                add(molecular_figure(displayed, title + "\nOriginal MapMyCells subclass; up to thirty common categories", "MapMyCells_subclass_display"), "original_MapMyCells_subclasses")
            add(molecular_figure(displayed, title, "condition"), "WT_cKO_all_nuclei")
            add(facet_figure(displayed, title), "all_twelve_mouse_facets")
            stage_summary = summary
            for key, value in zip(KEYS[:-1], keys):
                stage_summary = stage_summary.loc[stage_summary[key].eq(value)]
            if len(stage_summary):
                reps = representatives.loc[representatives.tool.eq(tool) & representatives.population.eq(population)]
                genes = reps.sort_values(["cluster_id", "representative_rank"]).gene_id.drop_duplicates().tolist()
                for start in range(0, len(genes), 28):
                    block = genes[start:start + 28]
                    fig, values = dot_figure(stage_summary, cluster_annotations, block,
                                            "Original-full specificity/mouse-supported representatives", title, cfg)
                    stem = f"representative_markers_{start // 28 + 1:02d}"
                    add(fig, stem); write_table(values, stage_out / (stem + "_exact_values.tsv"))
                for i, (label, genes) in enumerate(panel_pages(memberships), 1):
                    fig, values = dot_figure(stage_summary, cluster_annotations, genes,
                                            f"Supplied expression panels · page {i}: " + label.replace("_", " "), title, cfg)
                    stem = f"supplied_panel_{i:02d}"
                    add(fig, stem); write_table(values, stage_out / (stem + "_exact_values.tsv"))
                    write_table(memberships.loc[memberships.gene_id.isin(genes)], stage_out / (stem + "_membership.tsv"))
            else:
                # Secondary fresh labels have no matching native frozen-population expression denominators.
                (stage_out / "EXPRESSION_SCOPE.txt").write_text("Fresh graph views only. Native aggregates use frozen original populations; no dotplot is relabeled as fresh-cluster expression.\n")
        for field in ["reviewed_label", "broad_label", "condition"]:
            if field in displayed:
                for category in displayed[field].dropna().astype(str).unique():
                    palettes.append(dict(field=field, category=category, hex_color=category_color(category)))
        write_table(pd.DataFrame(pages), stage_out / "PAGE_INDEX.tsv")
        write_receipt(stage_out, args.config, source_paths + paths[:6],
                      dict(stage="reviewed_annotation_and_native_expression_rendering",
                           **dict(zip(KEYS[:-1], keys)), displayed_nuclei=len(displayed),
                           n_pages=len(pages), all_displayed_nuclei_rendered=True,
                           native_expression_aggregates_available=bool(len(stage_summary))))
        books.append(dict(zip(KEYS[:-1], keys)) | {"displayed_nuclei": len(displayed), "pages": len(pages),
                     "pdf": str((stage_out / "annotated_umaps_and_marker_dotplots.pdf").resolve()),
                     "source_metadata": str(source_paths[0].resolve()), "main_inference_scope": scope == "original_frozen",
                     "coordinate_fit_nuclei": len(metadata), "coordinate_refit_performed": False})
    write_table(pd.DataFrame(books), out / "PLOT_BOOK_INDEX.tsv")
    write_table(pd.DataFrame(palettes).drop_duplicates(), out / "CATEGORY_PALETTE.tsv")
    write_table(memberships.loc[memberships.gene_id.isna()], out / "UNMAPPED_SUPPLIED_PANEL_ROWS.tsv")
    write_receipt(out, args.config, paths, dict(stage="reviewed_annotation_and_native_expression_rendering",
                  n_plot_books=len(books), n_pages=sum(x["pages"] for x in books),
                  inference_scope="original_frozen; same reviewed analysis masks in main plots",
                  secondary_scope="fresh_branch coordinates/labels only; no frozen expression relabeling",
                  all_displayed_nuclei_rendered=True,
                  negative_z_policy="fraction unchanged; area multiplied by 0.25 and alpha 0.18 unless config overrides"))


if __name__ == "__main__":
    main()
