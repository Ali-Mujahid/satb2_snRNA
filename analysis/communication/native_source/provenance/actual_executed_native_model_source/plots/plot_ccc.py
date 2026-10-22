#!/usr/bin/env python3
"""Render source-backed, portable LIANA figures without fitting any models.

The input is the completed interaction_evidence.tsv[.gz] table. The script
never changes that table, the count matrices, or previous analyses.
"""
from __future__ import annotations

import argparse
import hashlib
import html
import json
import math
from pathlib import Path
import re
import shutil
import textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
from matplotlib.cm import ScalarMappable
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, PathPatch, Wedge
from matplotlib.path import Path as MplPath
import numpy as np
import pandas as pd

CONDITIONS = ["WT", "cKO", "pooled"]
DIRECTIONS = ["CA1_to_IN", "IN_to_CA1", "bidirectional"]
EVIDENCES = ["DEG_supported", "all_candidates"]
MANDATORY = ["source", "target", "ligand_complex", "receptor_complex",
             "lr_logfc", "lr_gmeans", "magnitude_rank", "scope", "direction",
             "deg_supported", "ligand_props", "receptor_props"]
FOOTNOTE = ("Predicted transcript-based candidates; ranks are not genotype P values. "
            "CA1 Deep/Superficial are expression programs, not measured physical layers.")
BLUE, ORANGE, GREY = "#2166ac", "#d46a23", "#8a96a2"


def args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True,
                        help="Completed evidence TSV/TSV.GZ, or a results directory containing tables/interaction_evidence.tsv.gz")
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--pair-eligibility", type=Path,
                        help="Explicit evaluated-pair TSV; defaults to condition_pair_eligibility.tsv beside the evidence table")
    parser.add_argument("--rank-cutoff", type=float, default=.05)
    parser.add_argument("--top-lr", type=int, default=24)
    parser.add_argument("--pairs-per-page", type=int, default=20)
    parser.add_argument("--dpi", type=int, default=180)
    return parser.parse_args()


def bools(values):
    out = values.astype(str).str.lower().str.strip()
    valid = out.isin(["true", "false", "1", "0", "yes", "no", "t", "f"])
    if not valid.all():
        raise ValueError(f"Unrecognized deg_supported values: {out[~valid].unique().tolist()}")
    return out.isin(["true", "1", "yes", "t"])


def natural(value):
    return [int(x) if x.isdigit() else x.lower() for x in re.split(r"(\d+)", str(value))]


def label(group):
    if group == "CA1_Deep":
        return "CA1 Deep"
    if group == "CA1_Superficial":
        return "CA1 Superficial"
    if group == "IN_Total":
        return "IN Total"
    return str(group).replace("IN_", "IN ")


def short_review(row, endpoint):
    return str(row.get(endpoint + "_review", "")).strip().lower() not in [
        "", "nan", "false", "none", "ok", "resolved", "pass", "0", "reference_supported"]


def load_data(path):
    if path.is_dir():
        path = path / "tables" / "interaction_evidence.tsv.gz"
    if not path.is_file():
        raise FileNotFoundError(path)
    data = pd.read_csv(path, sep="\t", low_memory=False)
    missing = set(MANDATORY) - set(data)
    if missing:
        raise ValueError("Missing evidence columns: " + ", ".join(sorted(missing)))
    if "condition" not in data:
        if "dataset" in data:
            data["condition"] = data["dataset"]
        else:
            raise ValueError("An explicit condition column (WT, cKO, pooled) is required.")
    condition_map = {"wt": "WT", "cko": "cKO", "pooled": "pooled", "all": "pooled", "allpooled": "pooled"}
    data["condition"] = data["condition"].astype(str).str.lower().map(condition_map)
    if data["condition"].isna().any():
        raise ValueError("Unexpected or absent condition label; expected WT, cKO or pooled.")
    data["deg_supported"] = bools(data["deg_supported"])
    data["source"] = data["source"].astype(str)
    data["target"] = data["target"].astype(str)
    for col in ["lr_logfc", "lr_gmeans", "magnitude_rank", "specificity_rank", "cellphone_pvals", "lr_means", "ligand_props", "receptor_props", "lr_prop", "ligand_all_subunit_min_prop", "receptor_all_subunit_min_prop"]:
        if col in data:
            data[col] = pd.to_numeric(data[col], errors="coerce")
    if not data["magnitude_rank"].dropna().between(0, 1).all():
        raise ValueError("magnitude_rank must be between 0 and 1, or NA.")
    for col in ["ligand_props", "receptor_props", "ligand_all_subunit_min_prop", "receptor_all_subunit_min_prop"]:
        if col not in data:
            continue
        if not data[col].dropna().between(0, 1).all():
            raise ValueError(f"{col} must be a fraction from 0 to 1, or NA.")
    subunit_columns = ["ligand_all_subunit_min_prop", "receptor_all_subunit_min_prop"]
    if all(col in data for col in subunit_columns):
        # The weakest detected subunit across both complete complexes, rather
        # than the subunit selected by native mean-expression bottleneck rules.
        data["plot_min_detection_fraction"] = data[subunit_columns].min(axis=1, skipna=False)
        data["plot_detection_definition"] = "minimum across all ligand and receptor complex subunits"
    else:
        data["plot_min_detection_fraction"] = data[["ligand_props", "receptor_props"]].min(axis=1, skipna=False)
        data["plot_detection_definition"] = "native selected endpoint subunits only; full-complex minima unavailable"
    data["plot_lr"] = data["ligand_complex"].astype(str) + " → " + data["receptor_complex"].astype(str)
    data["plot_pair"] = data["source"] + " → " + data["target"]
    keys = ["scope", "condition", "source", "target", "ligand_complex", "receptor_complex"]
    if data.duplicated(keys).any():
        raise ValueError("Duplicate interaction identity within scope/condition; aggregation must be resolved upstream.")
    return data, path


def setup():
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
        "axes.titlesize": 14, "axes.labelsize": 11, "xtick.labelsize": 9,
        "ytick.labelsize": 9, "pdf.fonttype": 42, "ps.fonttype": 42,
        "svg.fonttype": "none", "savefig.facecolor": "white", "axes.spines.top": False,
        "axes.spines.right": False})


def subset(data, scope, direction, evidence, cutoff):
    out = data.loc[data["scope"].eq(scope)].copy()
    if direction != "bidirectional":
        out = out.loc[out["direction"].eq(direction)]
    out = out.loc[out["magnitude_rank"].notna() & out["magnitude_rank"].le(cutoff)]
    if evidence == "DEG_supported":
        out = out.loc[out["deg_supported"]]
    return out.copy()


def title_text(scope, direction, condition):
    direct = {"CA1_to_IN": "CA1 → inhibitory", "IN_to_CA1": "Inhibitory → CA1", "bidirectional": "CA1 ↔ inhibitory"}[direction]
    granularity = "inhibitory clusters" if scope == "IN_clusters" else "inhibitory total"
    return f"{direct} | {condition} | {granularity}"


def caption(cutoff, evidence):
    support = "Endpoint-specific DESeq2 support" if evidence == "DEG_supported" else "All expression-supported candidates"
    return f"{support}; LIANA magnitude rank ≤ {cutoff:g} (lower is better)."


def save_figure(fig, prefix, dpi):
    for extension in ["pdf", "png", "svg"]:
        fig.savefig(prefix.with_suffix("." + extension), dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def write_rows(data, path):
    data.to_csv(path, sep="\t", index=False, float_format="%.10g", na_rep="NA")


def new_record(kind, prefix, outdir, condition, scope, direction, evidence, page=1, pages=1):
    relative = prefix.relative_to(outdir).as_posix()
    return {"kind": kind, "prefix": relative, "condition": condition, "scope": scope,
            "direction": direction, "evidence": evidence, "page": page, "pages": pages}


def node_info(data):
    mapping = {}
    for endpoint in ["source", "target"]:
        for _, row in data.iterrows():
            group = row[endpoint]
            if group not in mapping:
                mapping[group] = {"group": group, "label": label(group), "annotation": "", "review": False}
            annotation = row.get(endpoint + "_label", "")
            if pd.notna(annotation):
                mapping[group]["annotation"] = str(annotation)
            mapping[group]["review"] |= short_review(row, endpoint)
    return mapping


def plot_dotplots(allrows, current, prefix, nodes, order, pair_order, a, meta):
    """Consistent global LR order and size/rank scales across all conditions."""
    page_records = []
    pages = max(1, math.ceil(len(pair_order) / a.pairs_per_page))
    for page in range(pages):
        pairs = pair_order[page * a.pairs_per_page:(page + 1) * a.pairs_per_page]
        rows = current.loc[current["plot_lr"].isin(order) & current["plot_pair"].isin(pairs)].copy()
        p = prefix.parent / (prefix.name + (f"_page{page+1:02d}" if pages > 1 else ""))
        height = max(6, len(order) * .29 + 2.8)
        width = max(9, len(pairs) * .63 + 5)
        fig, ax = plt.subplots(figsize=(width, height))
        if order and pairs:
            xs = {x: i for i, x in enumerate(pairs)}
            ys = {x: i for i, x in enumerate(order)}
            rows["plot_x_index"] = rows["plot_pair"].map(xs)
            rows["plot_y_index"] = rows["plot_lr"].map(ys)
            rows["plot_dot_area_points2"] = 6 + 170 * rows["plot_min_detection_fraction"].fillna(0)
            # A fixed rank color scale permits direct visual comparison across conditions.
            scatter = ax.scatter(rows["plot_x_index"], rows["plot_y_index"],
                s=rows["plot_dot_area_points2"], c=rows["magnitude_rank"],
                cmap="viridis_r", vmin=0, vmax=a.rank_cutoff, edgecolors="none")
            ax.set_xticks(range(len(pairs)), [f"{label(x.split(' → ')[0])}\n→ {label(x.split(' → ')[1])}" for x in pairs], rotation=50, ha="right")
            ax.set_yticks(range(len(order)), order)
            ax.set_xlim(-.6, len(pairs) - .4)
            ax.set_ylim(len(order) - .4, -.6)
            ax.grid(color="#e8edf0", linewidth=.55)
            ax.set_axisbelow(True)
            cb = fig.colorbar(scatter, ax=ax, fraction=.025, pad=.025)
            cb.set_label("LIANA magnitude rank\n(lower = higher priority)")
            handles = [plt.scatter([], [], s=6 + 170*x, color=GREY, label=f"{x:.0%}") for x in [.1, .5, 1.]]
            full_complex = allrows.get("plot_detection_definition",pd.Series(dtype=str)).eq("minimum across all ligand and receptor complex subunits").all()
            detection_label = "Minimum of all complex\nsubunit detection fractions" if full_complex else "Native endpoint-subunit\ndetection fraction"
            ax.legend(handles=handles, title=detection_label, loc="upper left", bbox_to_anchor=(1.04, -.01), frameon=False)
        else:
            ax.axis("off")
        if rows.empty:
            ax.text(.5, .5, "No qualifying interactions in this panel.\nBlank positions are absent / below the selected rank criterion, not measured zeros.",
                transform=ax.transAxes, ha="center", va="center", color="#536372", fontsize=12)
        ax.set_title(meta["title"] + (f" | page {page+1}/{pages}" if pages > 1 else ""), loc="left", pad=28)
        fig.text(.01, .025, caption(a.rank_cutoff, meta["evidence"]) + "\n" + FOOTNOTE, fontsize=9, color="#536372")
        fig.subplots_adjust(left=.27, bottom=.27, right=.85, top=.89)
        save_figure(fig, p, a.dpi)
        write_rows(rows, Path(str(p) + "_SourceData.tsv"))
        record = new_record("dotplot", p, a.outdir, meta["condition"], meta["scope"], meta["direction"], meta["evidence"], page+1, pages)
        record["rows_plotted"] = len(rows)
        record["selection"] = f"Top {len(order)} ligand–receptor identities by best rank across all three conditions; common order and cell-pair pages."
        page_records.append(record)
    return page_records


def edge_table(rows):
    columns = ["source", "target", "n_predicted_interactions", "best_magnitude_rank", "median_magnitude_rank", "n_deg_supported"]
    if rows.empty:
        return pd.DataFrame(columns=columns)
    return rows.groupby(["source", "target"], sort=False).agg(
        n_predicted_interactions=("ligand_complex", "size"),
        best_magnitude_rank=("magnitude_rank", "min"),
        median_magnitude_rank=("magnitude_rank", "median"),
        n_deg_supported=("deg_supported", "sum")).reset_index()


def plot_heatmap(rows, edges, prefix, nodes, vmax, a, meta):
    names = list(nodes)
    mat = pd.DataFrame(np.nan, index=names, columns=names)
    # A pair must occur in the unfiltered native evidence table before an empty
    # qualifying set can be shown as zero. Missing/unevaluated pairs remain NA.
    tested_pairs = meta.get("tested_pairs", set(zip(rows.source, rows.target)))
    for src in names:
        for tgt in names:
            src_ca1, tgt_ca1 = src.startswith("CA1_"), tgt.startswith("CA1_")
            permitted = src_ca1 != tgt_ca1
            permitted &= meta["direction"] == "bidirectional" or (src_ca1 if meta["direction"] == "CA1_to_IN" else tgt_ca1)
            if permitted and (src, tgt) in tested_pairs:
                mat.loc[src, tgt] = 0
    for _, row in edges.iterrows():
        mat.loc[row.source, row.target] = row.n_predicted_interactions
    size = max(7.5, len(names) * .40 + 2)
    fig, ax = plt.subplots(figsize=(size + 1.4, size))
    cmap = plt.colormaps["YlGnBu"].copy()
    cmap.set_bad("#eff2f5")
    image = ax.imshow(mat.values, cmap=cmap, vmin=0, vmax=max(1, vmax), aspect="equal")
    shown = [label(x) + (" *" if nodes[x]["review"] else "") for x in names]
    ax.set_xticks(range(len(names)), shown, rotation=55, ha="right")
    ax.set_yticks(range(len(names)), shown)
    ax.set_xticks(np.arange(-.5, len(names), 1), minor=True)
    ax.set_yticks(np.arange(-.5, len(names), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.2)
    ax.tick_params(which="minor", bottom=False, left=False)
    for i, j in zip(*np.where(np.isfinite(mat.values))):
        value = mat.iat[i, j]
        ax.text(j, i, str(int(value)), ha="center", va="center", fontsize=8,
            color="white" if value > max(1, vmax) * .55 else "#243746")
    ax.set_xlabel("Target (receptor-expressing group)")
    ax.set_ylabel("Source (ligand-expressing group)")
    ax.set_title(meta["title"], loc="left", pad=16)
    fig.colorbar(image, ax=ax, fraction=.035, pad=.035, label="Predicted ligand–receptor candidate count")
    fig.text(.015, .02, caption(a.rank_cutoff, meta["evidence"]) + "\nGray tiles: untested direction / no native pair evidence. * Identity review flag.\n" + FOOTNOTE, fontsize=9, color="#536372")
    fig.subplots_adjust(left=.20, bottom=.25, right=.86, top=.9)
    save_figure(fig, prefix, a.dpi)
    grid = pd.DataFrame({"source":np.repeat(names, len(names)),
                         "target":np.tile(names, len(names)),
                         "n_predicted_interactions":mat.values.ravel()})
    grid["tested_direction"] = grid.n_predicted_interactions.notna()
    write_rows(grid, Path(str(prefix) + "_SourceData.tsv"))
    return {"rows_plotted": len(edges), "selection": "All tested CA1–inhibitory cell-group directions; candidate counts use the stated rank criterion."}


def node_color(name):
    if name == "CA1_Deep":
        return "#875636"
    if name == "CA1_Superficial":
        return "#397fb0"
    if name == "IN_Total":
        return "#bb6b25"
    match = re.search(r"C(\d+)", name)
    index = int(match.group(1)) if match else 0
    return matplotlib.colors.to_hex(plt.colormaps["tab20"]((index - 1) % 20))


def circle_positions(names):
    # Fixed across conditions: family blocks and natural cluster order.
    return {name: (math.cos(math.pi/2 - 2*math.pi*i/len(names)), math.sin(math.pi/2 - 2*math.pi*i/len(names))) for i, name in enumerate(names)} if names else {}


def draw_node_labels(ax, names, positions, nodes, radius=1.12):
    for name in names:
        x, y = positions[name]
        text = label(name) + (" *" if nodes[name]["review"] else "")
        ax.text(radius*x, radius*y, text, ha="center" if abs(x) < .2 else ("left" if x > 0 else "right"),
                va="center", fontsize=10, color="#243746")


def plot_circle(rows, edges, prefix, nodes, positions, edge_max, node_max, a, meta):
    fig, ax = plt.subplots(figsize=(10.5, 10))
    ax.set_aspect("equal")
    names = list(nodes)
    outgoing = edges.groupby("source")["n_predicted_interactions"].sum().to_dict() if len(edges) else {}
    incoming = edges.groupby("target")["n_predicted_interactions"].sum().to_dict() if len(edges) else {}
    for _, row in edges.sort_values("n_predicted_interactions").iterrows():
        src, tgt = row.source, row.target
        x1, y1 = positions[src]
        x2, y2 = positions[tgt]
        width = .6 + 7.4 * row.n_predicted_interactions / max(1, edge_max)
        color = BLUE if src.startswith("CA1_") else ORANGE
        patch = FancyArrowPatch((.94*x1, .94*y1), (.94*x2, .94*y2),
            connectionstyle="arc3,rad=0.13", arrowstyle="-|>", mutation_scale=12 + width,
            linewidth=width, color=color, alpha=.60, shrinkA=8, shrinkB=10, zorder=1)
        ax.add_patch(patch)
    node_rows = []
    for name in names:
        total = outgoing.get(name, 0) + incoming.get(name, 0)
        area = 100 + 950 * total / max(1, node_max)
        x, y = positions[name]
        ax.scatter([x], [y], s=area, color=node_color(name), edgecolor="white", linewidth=1.5, zorder=3)
        node_rows.append({**nodes[name], "x": x, "y": y, "n_outgoing": outgoing.get(name, 0), "n_incoming": incoming.get(name, 0), "plot_node_area_points2": area})
    draw_node_labels(ax, names, positions, nodes)
    if edges.empty:
        ax.text(0, 0, "No qualifying interactions", ha="center", va="center", color="#536372", fontsize=14)
    ax.legend(handles=[Line2D([0], [0], color=BLUE, lw=3, label="CA1 → inhibitory"),
        Line2D([0], [0], color=ORANGE, lw=3, label="Inhibitory → CA1")], loc="lower center", bbox_to_anchor=(.5, -.08), ncol=2, frameon=False)
    ax.set_xlim(-1.6, 1.6)
    ax.set_ylim(-1.35, 1.45)
    ax.axis("off")
    ax.set_title(meta["title"], loc="left", pad=15)
    fig.text(.02, .025, caption(a.rank_cutoff, meta["evidence"]) + "\nArrow width = candidate count; node area = incoming + outgoing count. Common scales across conditions.\n* Identity review flag. " + FOOTNOTE, fontsize=9, color="#536372")
    fig.subplots_adjust(top=.9, bottom=.15, left=.1, right=.9)
    save_figure(fig, prefix, a.dpi)
    write_rows(edges, Path(str(prefix) + "_SourceData.tsv"))
    write_rows(pd.DataFrame(node_rows), Path(str(prefix) + "_Nodes.tsv"))
    return {"rows_plotted": len(edges), "selection": "All qualifying source–target edges; widths encode interaction identities, not nuclei or validated signaling."}


def arc_points(a0, a1, radius, n=12):
    t = np.linspace(a0, a1, n)
    return np.column_stack([radius*np.cos(t), radius*np.sin(t)])


def ribbon(ax, source_interval, target_interval, color):
    """A closed, finite-width chord ribbon; receiver is marked by a tip."""
    s0, s1 = source_interval
    t0, t1 = target_interval
    source = arc_points(s0, s1, .92)
    target = arc_points(t0, t1, .92)
    verts = list(source)
    codes = [MplPath.MOVETO] + [MplPath.LINETO] * (len(source)-1)
    verts += [source[-1]*.22, target[0]*.22, target[0]]
    codes += [MplPath.CURVE4]*3
    verts += list(target[1:])
    codes += [MplPath.LINETO]*(len(target)-1)
    verts += [target[-1]*.22, source[0]*.22, source[0], source[0]]
    codes += [MplPath.CURVE4]*3 + [MplPath.CLOSEPOLY]
    ax.add_patch(PathPatch(MplPath(np.asarray(verts), codes), facecolor=color, edgecolor="white", linewidth=.35, alpha=.50, zorder=1))
    midpoint = (t0+t1)/2
    ax.annotate("", xy=(.972*math.cos(midpoint), .972*math.sin(midpoint)),
        xytext=(.86*math.cos(midpoint), .86*math.sin(midpoint)),
        arrowprops={"arrowstyle": "-|>", "color": color, "lw": .7, "mutation_scale": 9}, zorder=4)


def plot_circos(rows, edges, prefix, nodes, a, meta):
    names = list(nodes)
    weights = {name: 0. for name in names}
    for _, row in edges.iterrows():
        weights[row.source] += row.n_predicted_interactions
        weights[row.target] += row.n_predicted_interactions
    # Equal sectors preserve node order across conditions; sub-arc lengths are
    # proportional within each sector. Ribbon width therefore is NOT globally
    # comparable across sectors; the exact candidate counts are printed/exported.
    fig, ax = plt.subplots(figsize=(10.5, 10))
    ax.set_aspect("equal")
    gap, theta = .025, math.pi/2
    span = 2*math.pi/max(1, len(names)) - gap
    intervals, cursors, positions = {}, {}, {}
    for name in names:
        intervals[name] = (theta, theta+span)
        cursors[name] = theta
        midpoint = theta + span/2
        positions[name] = (math.cos(midpoint), math.sin(midpoint))
        ax.add_patch(Wedge((0, 0), 1.02, math.degrees(theta), math.degrees(theta+span),
            width=.07, facecolor=node_color(name), edgecolor="white", linewidth=.7, zorder=3))
        theta += span + gap
    rendered_edges = []
    for _, row in edges.sort_values(["source", "target"], key=lambda c:c.map(str)).iterrows():
        s, t, count = row.source, row.target, row.n_predicted_interactions
        s0, t0 = cursors[s], cursors[t]
        s1 = s0 + span * count/max(1., weights[s])
        t1 = t0 + span * count/max(1., weights[t])
        color = BLUE if s.startswith("CA1_") else ORANGE
        ribbon(ax, (s0, s1), (t0, t1), color)
        cursors[s], cursors[t] = s1, t1
        rendered_edges.append({**row.to_dict(), "source_arc_start":s0, "source_arc_end":s1,
            "target_arc_start":t0, "target_arc_end":t1, "color":color})
    draw_node_labels(ax, names, positions, nodes, 1.16)
    if edges.empty:
        ax.text(0, 0, "No qualifying interactions", ha="center", va="center", color="#536372", fontsize=14)
    ax.legend(handles=[Line2D([0], [0], color=BLUE, lw=5, alpha=.6, label="CA1 → inhibitory"),
        Line2D([0], [0], color=ORANGE, lw=5, alpha=.6, label="Inhibitory → CA1")],
        loc="lower center", bbox_to_anchor=(.5, -.08), ncol=2, frameon=False)
    ax.set_xlim(-1.6, 1.6)
    ax.set_ylim(-1.4, 1.4)
    ax.axis("off")
    ax.set_title(meta["title"], loc="left", pad=15)
    fig.text(.02, .025, caption(a.rank_cutoff, meta["evidence"]) + "\nEqual sectors; ribbon width = within-sector fraction of endpoint candidate counts. Tips identify receivers.\nWidth across different sectors is not a common count scale. * Identity review flag.\n" + FOOTNOTE, fontsize=9, color="#536372")
    fig.subplots_adjust(top=.9, bottom=.16, left=.1, right=.9)
    save_figure(fig, prefix, a.dpi)
    write_rows(pd.DataFrame(rendered_edges, columns=list(edges.columns) + ["source_arc_start", "source_arc_end", "target_arc_start", "target_arc_end", "color"]), Path(str(prefix) + "_SourceData.tsv"))
    return {"rows_plotted": len(edges), "selection": "True chord ribbons aggregate candidate identities by source–target group; equal node sectors maintain legible common order."}


def html_report(records, summary, outdir, provenance):
    encoded = json.dumps({"figures": records, "summary": summary, "provenance": provenance}, ensure_ascii=False).replace("</", "<\\/")
    template = r'''<!DOCTYPE html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>SATB2 CA1–inhibitory interaction atlas</title><style>
    :root{font-family:system-ui,-apple-system,sans-serif;color:#213746;background:#edf3f6}*{box-sizing:border-box}body{margin:0}header{background:#173d58;color:white;padding:30px max(24px,5vw)}h1{font-size:clamp(24px,3vw,36px);margin:0 0 10px}header p{max-width:1000px;margin:8px 0;line-height:1.5}main{max-width:1550px;margin:24px auto;padding:0 22px}.warning,.card{background:white;border:1px solid #d1dde5;border-radius:12px;padding:20px;margin:16px 0}.warning{border-left:5px solid #c28433;background:#fffaf1}.controls{display:grid;grid-template-columns:repeat(auto-fit,minmax(190px,1fr));gap:15px}label{font-weight:650;font-size:14px;display:block}select{display:block;margin-top:7px;padding:10px;width:100%;border:1px solid #aebfcd;border-radius:6px;background:white;color:#213746;font:inherit}.downloads{display:flex;flex-wrap:wrap;gap:10px;margin:16px 0}a.button{padding:9px 14px;background:#eaf3fa;border:1px solid #b9cedd;border-radius:6px;text-decoration:none;color:#17496e}a{color:#155a8b}#figure{width:100%;height:auto;background:white}#status{font-weight:600}#description{color:#526877;line-height:1.5}.stats{display:flex;gap:24px;flex-wrap:wrap}.stats strong{font-size:23px;display:block}details{margin-top:18px}summary{cursor:pointer;font-weight:650}table{border-collapse:collapse;width:100%;font-size:13px;margin-top:15px}td,th{border-bottom:1px solid #d9e3e9;padding:9px;text-align:left}th{background:#f0f5f8}.scroll{overflow:auto;max-height:450px}.small{font-size:13px;color:#526877}.key{display:grid;grid-template-columns:repeat(auto-fit,minmax(240px,1fr));gap:15px}.key p{line-height:1.5}.empty{padding:50px;text-align:center}footer{padding:24px 0;color:#526877;font-size:13px}button:focus,select:focus,a:focus{outline:3px solid #cc8030;outline-offset:3px}@media print{header,.controls,.downloads,.warning,details,footer{display:none}main{padding:0}.card{border:none}#figure{max-height:90vh;object-fit:contain}}
    </style></head><body><header><h1>SATB2 CA1–inhibitory interaction atlas</h1><p>Seurat populations · CA1 Deep and Superficial expression programs · ligand–receptor predictions with LIANA and CellPhoneDB evidence</p><p>Portable report: open this HTML file locally. Keep the figures and data folders beside it. No internet or analysis server is required.</p></header><main>
    <div class="warning"><strong>Read as candidate signaling, not confirmed communication.</strong> RNA abundance does not establish protein availability, anatomical contact, synapses, or signaling activity. CA1 Deep/Superficial refer to molecular programs; physical layers were not measured. Condition-specific LIANA ranks and CellPhoneDB permutation P values do not test the WT–cKO difference.</div>
    <section class="card controls" aria-label="Figure controls"><label>Population resolution<select id="scope"><option value="IN_clusters">Inhibitory clusters (primary)</option><option value="IN_total">Inhibitory total (supplementary)</option></select></label><label>Interaction direction<select id="direction"><option value="CA1_to_IN">CA1 → inhibitory</option><option value="IN_to_CA1">Inhibitory → CA1</option><option value="bidirectional">Both directions</option></select></label><label>Condition<select id="condition"><option>WT</option><option>cKO</option><option value="pooled">Pooled (context)</option></select></label><label>Evidence view<select id="evidence"><option value="DEG_supported">Endpoint-specific DEG support</option><option value="all_candidates">All candidates (context)</option></select></label><label>Plot<select id="kind"><option value="dotplot">Ligand–receptor dotplot</option><option value="heatmap">Source–target tileplot</option><option value="circos">Circos / chord ribbons</option><option value="circle">Directed circle network</option></select></label><label id="pageLabel">Dotplot page<select id="page"></select></label></section>
    <section class="card"><div class="stats" id="stats"></div><p id="status" aria-live="polite"></p><p id="description"></p><div class="downloads" id="downloads"></div><img id="figure" alt="Selected ligand–receptor figure"><p id="missing" class="empty" hidden>No figure is available for these selections. Empty results are not evidence of absent biological communication.</p></section>
    <section class="card"><h2>How to read the plots</h2><div class="key"><p><strong>Dotplot</strong><br>Color is LIANA magnitude rank; lower ranks receive higher priority. Dot area preferably shows the minimum detection fraction across every subunit of the source ligand complex and target receptor complex. If only native selected-subunit fractions are available, the plot legend and manifest explicitly identify that fallback. Native audit columns remain unchanged. The same ligand–receptor order, rank scale and size scale are used across conditions. Missing dots are absent or below the displayed rank criterion.</p><p><strong>Tileplot</strong><br>Counts distinct predicted ligand–receptor identities for each tested source–target direction. Gray tiles were not tested. Zero is an empty qualifying set, not a measured absence of signaling. All eligible groups remain visible.</p><p><strong>Circos</strong><br>Closed ribbons connect source and receiver sectors. Tips mark receivers. Blue starts in CA1; orange starts in inhibitory neurons. Equal sectors keep the groups legible; ribbon width is the within-sector fraction and cannot be compared as an absolute count across sectors.</p><p><strong>Circle</strong><br>Arrows point from ligand source to receptor target. Width is the number of predicted interaction identities. Node area is incoming plus outgoing count. Scales and positions are fixed across conditions. These are neither cell counts nor significant genotype-change counts.</p></div>
    <p><strong>DEG-supported view:</strong> at least one ligand subunit has DESeq2 support in the matched source, or one receptor subunit has support in the matched target. This prioritizes candidates; it does not establish a statistically significant interaction change. The upstream evidence table documents the exact endpoint tests and complex rules.</p><p><strong>Table score:</strong> native LIANA <code>lr_logfc</code> compares group expression with the remaining cells; it is not the DESeq2 cKO/WT log2 fold change.</p>
    <p><strong>Identity review:</strong> an asterisk marks groups with provisional or mixed identity evidence. All-group tables retain these annotations. Annotation uncertainty should travel with any biological interpretation.</p>
    <details><summary>Interaction totals for all selections</summary><div class="scroll"><table><thead><tr><th>Scope</th><th>Direction</th><th>Condition</th><th>Evidence</th><th>Interactions</th><th>Source–target edges</th></tr></thead><tbody id="summaryBody"></tbody></table></div></details>
    <details><summary>Files and reproducibility</summary><p><a href="data/interaction_evidence.tsv.gz" download>Full unfiltered interaction evidence (TSV.GZ)</a> · <a href="data/plot_manifest.json" download>Plot manifest and selection rules</a> · <a href="data/node_annotations.tsv" download>Node annotations and identity-review flags</a> · <a href="README.md">Reading guide</a></p><p class="small" id="provenance"></p></details></section><footer>This report contains only completed saved result tables and precomputed figures. Changing controls does not rerun LIANA or DESeq2.</footer></main><script id="reportData" type="application/json">__DATA__</script><script>
    const report=JSON.parse(document.getElementById('reportData').textContent),ids=['scope','direction','condition','evidence','kind'],el=Object.fromEntries([...ids,'page'].map(id=>[id,document.getElementById(id)]));
    function selected(){return Object.fromEntries(ids.map(id=>[id,el[id].value]));}function same(r,s){return ids.every(id=>r[id]===s[id]);}
    function update(repage=true){const s=selected(),matches=report.figures.filter(r=>same(r,s));if(repage){el.page.replaceChildren();for(const r of matches){const o=document.createElement('option');o.value=r.page;o.textContent=`${r.page} / ${r.pages}`;el.page.append(o);}}document.getElementById('pageLabel').style.display=(s.kind==='dotplot'&&matches.length>1)?'block':'none';const r=matches.find(x=>x.page===+el.page.value)||matches[0],image=document.getElementById('figure'),downloads=document.getElementById('downloads');downloads.replaceChildren();image.hidden=!r;document.getElementById('missing').hidden=!!r;document.getElementById('status').textContent=r?`${s.condition} · ${s.scope==='IN_clusters'?'inhibitory clusters':'inhibitory total'} · ${s.direction.replaceAll('_',' ')}`:'';document.getElementById('description').textContent=r?r.selection:'';if(r){image.src=r.prefix+'.png';image.alt=`${s.kind}, ${s.direction}, ${s.condition}, ${s.evidence}; see source data and reading guide`;for(const ext of ['pdf','png','svg']){const a=document.createElement('a');a.className='button';a.href=r.prefix+'.'+ext;a.download='';a.textContent='Download '+ext.toUpperCase();downloads.append(a);}const a=document.createElement('a');a.className='button';a.href=r.prefix+'_SourceData.tsv';a.download='';a.textContent='Exact plot data (TSV)';downloads.append(a);}const stat=report.summary.find(r=>['scope','direction','condition','evidence'].every(k=>r[k]===s[k]));const stats=document.getElementById('stats');stats.replaceChildren();if(stat){for(const[k,t]of[['interactions','Predicted interaction identities'],['edges','Source–target edges']]){const d=document.createElement('div'),b=document.createElement('strong');b.textContent=stat[k].toLocaleString();d.append(b,document.createTextNode(t));stats.append(d);}}}
    for(const id of ids)el[id].addEventListener('change',()=>update(true));el.page.addEventListener('change',()=>update(false));for(const row of report.summary){const tr=document.createElement('tr');for(const k of ['scope','direction','condition','evidence','interactions','edges']){const td=document.createElement('td');td.textContent=row[k];tr.append(td);}document.getElementById('summaryBody').append(tr);}document.getElementById('provenance').textContent=`Input SHA-256: ${report.provenance.input_sha256}. Rank criterion: ≤ ${report.provenance.rank_cutoff}. Input rows: ${report.provenance.input_rows.toLocaleString()}.`;update();
    </script></body></html>'''
    (outdir / "index.html").write_text(template.replace("__DATA__", encoded), encoding="utf-8")


def main():
    a = args()
    if not 0 < a.rank_cutoff <= 1 or not 1 <= a.top_lr <= 60 or not 1 <= a.pairs_per_page <= 30:
        raise ValueError("Use rank cutoff in (0,1], 1–60 LR rows and 1–30 cell pairs per page.")
    a.outdir = a.outdir.resolve()
    data, input_path = load_data(a.results.resolve())
    eligibility_path = (a.pair_eligibility or input_path.parent / "condition_pair_eligibility.tsv").resolve()
    if not eligibility_path.is_file():
        raise FileNotFoundError(f"An explicit evaluated-pair table is required to distinguish zeros from untested pairs: {eligibility_path}")
    eligibility = pd.read_csv(eligibility_path, sep="\t")
    pair_columns = ["scope", "condition", "source", "target", "source_nuclei", "target_nuclei"]
    if set(pair_columns) - set(eligibility):
        raise ValueError("Missing columns in pair eligibility table: " + ", ".join(sorted(set(pair_columns) - set(eligibility))))
    if eligibility.duplicated(["scope", "condition", "source", "target"]).any():
        raise ValueError("Duplicate evaluated-pair identity in the eligibility table.")
    eligibility["condition"] = eligibility["condition"].astype(str).str.lower().map({"wt":"WT", "cko":"cKO", "pooled":"pooled"})
    if eligibility.condition.isna().any():
        raise ValueError("Unexpected condition in evaluated-pair table.")
    for col in ["source_nuclei", "target_nuclei"]:
        if not pd.to_numeric(eligibility[col], errors="coerce").ge(5).all():
            raise ValueError("The evaluated-pair table must contain eligible groups only (at least 5 nuclei per endpoint).")
    # A finished output directory is immutable by default; use a new path for a rerender.
    if (a.outdir / "data" / "plot_manifest.json").exists():
        raise FileExistsError(f"Refusing to overwrite completed report: {a.outdir}")
    (a.outdir / "figures").mkdir(parents=True, exist_ok=True)
    (a.outdir / "data").mkdir(parents=True, exist_ok=True)
    setup()
    available_scopes = set(data.scope) | set(eligibility.scope)
    scope_names = [x for x in ["IN_clusters", "IN_total"] if x in available_scopes]
    if available_scopes - set(scope_names):
        raise ValueError("Unexpected scope; expected IN_clusters or IN_total.")
    key_columns = ["scope", "condition", "source", "target"]
    pair_index = pd.MultiIndex.from_frame(eligibility[key_columns])
    if not pd.MultiIndex.from_frame(data[key_columns]).isin(pair_index).all():
        raise ValueError("Some interaction rows are absent from the explicit evaluated-pair table.")
    records, summaries = [], []
    allnodes = node_info(data)
    for name in set(eligibility.source) | set(eligibility.target):
        if name not in allnodes:
            allnodes[name] = {"group":name,"label":label(name),"annotation":"No native interaction rows; consult upstream group annotations.","review":True}
    write_rows(pd.DataFrame(list(allnodes.values())), a.outdir / "data" / "node_annotations.tsv")
    for scope in scope_names:
        base = data.loc[data.scope.eq(scope)]
        scope_eligible = eligibility.loc[eligibility.scope.eq(scope)]
        names = sorted(set(base.source) | set(base.target) | set(scope_eligible.source) | set(scope_eligible.target), key=lambda x:(not x.startswith("CA1_"), natural(x)))
        nodes = {name:allnodes[name] for name in names}
        positions = circle_positions(names)
        for direction in DIRECTIONS:
            for evidence in EVIDENCES:
                allrows = subset(data, scope, direction, evidence, a.rank_cutoff)
                ranked_lr = allrows.groupby("plot_lr").magnitude_rank.min().sort_values(kind="stable")
                order = ranked_lr.head(a.top_lr).index.tolist()
                pair_order = sorted(allrows.loc[allrows.plot_lr.isin(order), "plot_pair"].unique(), key=natural)
                edges_all = [edge_table(allrows.loc[allrows.condition.eq(c)]) for c in CONDITIONS]
                edge_max = max([1.] + [float(e.n_predicted_interactions.max()) for e in edges_all if len(e)])
                node_max = max([1.] + [float(pd.concat([e.groupby("source").n_predicted_interactions.sum(),e.groupby("target").n_predicted_interactions.sum()],axis=1).fillna(0).sum(axis=1).max()) for e in edges_all if len(e)])
                for condition, edges in zip(CONDITIONS, edges_all):
                    rows = allrows.loc[allrows.condition.eq(condition)].copy()
                    folder = a.outdir / "figures" / scope / direction / evidence / condition
                    folder.mkdir(parents=True, exist_ok=True)
                    meta = {"condition":condition,"scope":scope,"direction":direction,"evidence":evidence,"title":title_text(scope,direction,condition)}
                    tested = scope_eligible.loc[scope_eligible.condition.eq(condition)]
                    meta["tested_pairs"] = set(zip(tested.source, tested.target))
                    summaries.append({**{k:meta[k] for k in ["scope","direction","condition","evidence"]},"interactions":len(rows),"edges":len(edges)})
                    records.extend(plot_dotplots(allrows, rows, folder / "dotplot", nodes, order, pair_order, a, meta))
                    for kind, fun in [("heatmap",plot_heatmap),("circle",plot_circle),("circos",plot_circos)]:
                        prefix = folder / kind
                        if kind == "heatmap":
                            detail = fun(rows, edges, prefix, nodes, edge_max, a, meta)
                        elif kind == "circle":
                            detail = fun(rows, edges, prefix, nodes, positions, edge_max, node_max, a, meta)
                        else:
                            detail = fun(rows, edges, prefix, nodes, a, meta)
                        records.append({**new_record(kind,prefix,a.outdir,condition,scope,direction,evidence),**detail})
                    write_rows(rows, folder / "qualifying_interactions.tsv.gz")
                    write_rows(rows[["source", "target", "ligand_complex", "receptor_complex", "lr_logfc", "lr_gmeans", "magnitude_rank"]], folder / "liana_rank_aggregate_7columns.tsv")
                    write_rows(edges, folder / "source_target_counts.tsv")
                    print(f"Rendered {scope} / {direction} / {evidence} / {condition}: {len(rows)} interactions", flush=True)
    portable_input = a.outdir / "data" / "interaction_evidence.tsv.gz"
    if input_path.name.endswith(".gz"):
        shutil.copyfile(input_path, portable_input)
    else:
        pd.read_csv(input_path,sep="\t",low_memory=False).to_csv(portable_input,sep="\t",index=False,compression="gzip")
    shutil.copyfile(eligibility_path, a.outdir / "data" / "condition_pair_eligibility.tsv")
    provenance = {"input_sha256":hashlib.sha256(input_path.read_bytes()).hexdigest(),"input_rows":len(data),"input_file":str(input_path),"pair_eligibility_sha256":hashlib.sha256(eligibility_path.read_bytes()).hexdigest(),"rank_cutoff":a.rank_cutoff,"top_lr":a.top_lr,"pairs_per_page":a.pairs_per_page,"dot_detection_definition":data.plot_detection_definition.unique().tolist(),"methods":"Render completed saved LIANA/DESeq2 evidence only; no statistical fits or differential interaction tests."}
    manifest = {"provenance":provenance,"figures":records,"summary":summaries}
    (a.outdir / "data" / "plot_manifest.json").write_text(json.dumps(manifest,indent=2,ensure_ascii=False),encoding="utf-8")
    html_report(records,summaries,a.outdir,provenance)
    readme = Path(__file__).with_name("README.md")
    if readme.is_file():
        shutil.copyfile(readme,a.outdir / "README.md")
    print(f"Portable report: {a.outdir / 'index.html'}",flush=True)


if __name__ == "__main__":
    main()
