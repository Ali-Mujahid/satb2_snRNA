#!/usr/bin/env python3
"""Render completed DESeq2 results and existing mouse pseudobulk arrays only.

No fitting, raw-nucleus aggregation, membership changes or native-output edits.
The native final-comparison receipt is mandatory. --check performs readiness
inspection only and never creates the destination.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import textwrap

import numpy as np
import pandas as pd
from scipy.stats import rankdata

POPS = ["allcells", "Excitatory", "CA1_ProS", "Inhibitory"]
NAMES = {"allcells": "All selected nuclei", "Excitatory": "Excitatory neurons",
         "CA1_ProS": "CA1/ProS candidates", "Inhibitory": "Inhibitory neurons"}
BRANCHES = ["original_full", "paired_original", "paired_cellbender"]
BRANCH_LABELS = {"original_full": "Original full", "paired_original": "Paired raw",
                 "paired_cellbender": "Paired CellBender"}
COL = {"WT": "#0072B2", "cKO": "#E69F00", "robust": "#009E73", "hit": "#8064A2",
       "other": "#B8BEC6", "single": "#D55E00", "missing": "#3F4650", "unsupported": "#D7DBDF"}
STATS = ["baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda: f.read(8 * 1024 * 1024), b""): h.update(chunk)
    return h.hexdigest()


def truth(values):
    return values.astype(str).str.lower().isin(["true", "1"])


def readiness(cfg):
    root = Path(cfg["results_root"])
    required = [root / "comparison/COMPLETE.json", root / "count_validation/COMPLETE.json"]
    required += [root / "deseq2" / b / "COMPLETE.json" for b in BRANCHES]
    return {"ready": all(p.exists() for p in required),
            "missing_producer_receipts": [str(p) for p in required if not p.exists()]}


def setup_plotting():
    os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parents[1] / "reviews/.matplotlib_effect_figures"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
        "axes.titlesize": 12, "axes.titleweight": "bold", "axes.labelsize": 10,
        "xtick.labelsize": 8.5, "ytick.labelsize": 9, "legend.fontsize": 9,
        "axes.spines.top": False, "axes.spines.right": False,
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
        "svg.hashsalt": "SATB2_20261014_DESeq2_effect_figures", "savefig.dpi": 300})
    return plt, PdfPages


def caption(fig, title, text):
    fig.suptitle(title, x=.055, y=.98, ha="left", fontsize=16, fontweight="bold")
    fig.text(.055, .025, "\n".join(textwrap.wrap(text, 155)), ha="left", va="bottom", fontsize=8.6, color="#3F4650")


def style(ax):
    ax.grid(axis="y", color="#E8EBEE", linewidth=.6, zorder=0)
    ax.set_axisbelow(True)


class FigureWriter:
    def __init__(self, out, plt, PdfPages, synthetic=False):
        self.out, self.plt, self.rows = Path(out), plt, []
        self.book = PdfPages(self.out / "DESeq2_effect_overview.pdf", metadata={"Title": "SATB2 DESeq2 effect overview", "Author": "Source-bound analytical supplement"})
        self.synthetic = synthetic

    def save(self, fig, key, title, text, population="all_four", branch="mixed", scope="produced_result_tables_only"):
        if self.synthetic:
            fig.text(.5, .50, "FABRICATED LAYOUT FIXTURE", ha="center", va="center", fontsize=27, color="#B22222", alpha=.16, rotation=25)
        paths = {}
        for ext in ["pdf", "svg", "png"]:
            p = self.out / f"{key}.{ext}"
            fig.savefig(p, dpi=300, facecolor="white", metadata={"Creator": "build_effect_figures.py"} if ext == "svg" else None)
            paths[f"{ext}_path"] = str(p.resolve())
        self.book.savefig(fig, dpi=300)
        self.rows.append({"figure_id": key, "title": title, "population": population, "branch": branch,
                          "method": "DESeq2", "data_scope": scope, "caption": text, **paths})
        self.plt.close(fig)

    def close(self):
        self.book.close()
        pd.DataFrame(self.rows).to_csv(self.out / "FIGURE_INDEX.tsv", sep="\t", index=False)


def draw_status(plt, writer, audits):
    fig, axes = plt.subplots(1, 2, figsize=(11.6, 6.4))
    x = np.arange(len(POPS)); width = .23
    for ax, tool in zip(axes, ["scanpy", "seurat"]):
        for j, branch in enumerate(BRANCHES):
            d = audits[(audits.tool == tool) & (audits.branch == branch) & audits.endpoint.eq("PRIMARY")]
            tested = np.array([int((d.population.eq(p) & d.status.eq("TESTED_REVIEW_REQUIRED")).sum()) for p in POPS])
            unsupported = np.array([int((d.population.eq(p) & d.status.eq("NOT_TESTED")).sum()) for p in POPS])
            pos = x + (j - 1) * width
            ax.bar(pos, tested, width=width*.9, color=["#4D89BA", "#438E79", "#8766A6"][j], label=BRANCH_LABELS[branch])
            ax.bar(pos, unsupported, bottom=tested, width=width*.9, color=COL["unsupported"], hatch="///", edgecolor="white", linewidth=.4)
            for a, n, u in zip(pos, tested, unsupported):
                ax.text(a, n+u+.4, f"{n}/{u}", ha="center", va="bottom", fontsize=7.5)
        ax.set_xticks(x, ["All", "EX", "CA1/ProS", "IN"]); ax.set_title(tool.capitalize())
        ax.set_ylabel("Number of PRIMARY population models"); ax.set_ylim(0, max(1, ax.get_ylim()[1])*1.15); style(ax)
    axes[0].legend(frameon=False, loc="upper left", fontsize=8.5)
    axes[1].text(.01, .98, "Labels: tested / unsupported\nHatched gray: no fitted result", transform=axes[1].transAxes, va="top", fontsize=8.5)
    text = "Every requested PRIMARY endpoint is counted, including totals, clusters and overlapping reviewed categories/programs. Unsupported models include emptied or excluded groups and insufficient independent-mouse support. Bars count exploratory population definitions, not independent discoveries; BH correction is within each fitted endpoint."
    caption(fig, "01  Which population models have enough support?", text); fig.subplots_adjust(left=.075, right=.98, top=.86, bottom=.22, wspace=.28)
    writer.save(fig, "01_model_support", "PRIMARY models: tested and unsupported", text)


def draw_coverage(plt, writer, coverage):
    fig, axes = plt.subplots(2, 2, figsize=(11.6, 8.5))
    for ax, pop in zip(axes.flat, POPS):
        a, b = coverage[(pop, "original_full")], coverage[(pop, "paired_original")]
        b = b.set_index("sample").loc[a["sample"]].reset_index(); x = np.arange(len(a))
        colors = [COL[v] for v in a.condition]
        ax.bar(x-.18, a.n_nuclei, width=.34, color=colors, alpha=.9, label="Original full")
        ax.bar(x+.18, b.n_nuclei, width=.34, color="white", edgecolor=colors, hatch="///", linewidth=1, label="Paired raw / CellBender")
        ax.set_xticks(x, [s.split("_")[-1] for s in a["sample"]], rotation=45, ha="right")
        ax.set_title(NAMES[pop]); ax.set_ylabel("Nuclei contributing to the total"); ax.set_ylim(0, max(a.n_nuclei.max(), b.n_nuclei.max())*1.18)
        ax.axvline(5.5, color="#CCD1D6", linewidth=.8); ax.text(.02, .96, "WT (blue)                         cKO (orange)", transform=ax.transAxes, va="top", fontsize=8)
        style(ax)
    axes[0, 0].legend(frameon=False, loc="upper left", bbox_to_anchor=(0, 1.23), ncol=2, fontsize=8)
    text = "Bars use existing full-count population coverage tables for the same twelve independent mice. Paired raw and CellBender have identical nucleus memberships. All/EX use the shared iterative cluster-size filter; CA1/ProS and IN retain small original clusters. These are sampling/coverage counts, not a cell-abundance hypothesis test."
    caption(fig, "02  How many nuclei does each mouse contribute?", text); fig.subplots_adjust(left=.08, right=.98, top=.82, bottom=.18, hspace=.56, wspace=.26)
    writer.save(fig, "02_core_mouse_coverage", "Shared core totals: all twelve mice", text, scope="existing_all12_population_coverage_tables")


def label_genes(ax, frame, robust, maximum=4):
    selected = frame[frame.gene_id.isin(robust) & frame.display_symbol.notna() & np.isfinite(frame.padj)].sort_values(["padj", "gene_id"]).head(maximum)
    side_count = {True: 0, False: 0}
    for row in selected.itertuples(index=False):
        right = row.log2FoldChange >= 0
        rank = side_count[right]; side_count[right] += 1
        ax.annotate(str(row.display_symbol), xy=(row.log2FoldChange, -np.log10(max(row.padj, np.nextafter(0., 1.)))),
          xytext=(.96 if right else .04, .92-.11*rank), textcoords="axes fraction", ha="right" if right else "left", va="top",
          fontsize=9, color="#006B4E", bbox=dict(boxstyle="round,pad=.2", facecolor="white", edgecolor="none", alpha=.9),
          arrowprops=dict(arrowstyle="-", color="#648379", linewidth=.65))


def draw_volcano(plt, writer, effects, robust):
    fig, axes = plt.subplots(2, 2, figsize=(11.6, 8.5))
    for ax, pop in zip(axes.flat, POPS):
        z = effects[(pop, "original_full")]
        if z is None:
            ax.text(.5,.5,"Core model unsupported",ha="center",transform=ax.transAxes); ax.set_title(NAMES[pop]); continue
        ok = np.isfinite(z.log2FoldChange) & np.isfinite(z.padj)
        hit = z.padj.lt(.05) & ok; green = z.gene_id.isin(robust[pop]) & ok
        yy = -np.log10(z.padj.clip(lower=np.nextafter(0., 1.)))
        for mask, color, size in [(ok & ~hit, COL["other"], 5), (hit & ~green, COL["hit"], 8), (green, COL["robust"], 12)]:
            ax.scatter(z.loc[mask,"log2FoldChange"],yy[mask],c=color,s=size,alpha=.72,linewidths=0,rasterized=True)
        ax.axhline(-np.log10(.05),color="#8C939B",linestyle="--",linewidth=.8); ax.axvline(0,color="#CCD1D6",linewidth=.8)
        ax.set_title(NAMES[pop]); ax.set_xlabel("cKO / WT log₂ fold change"); ax.set_ylabel("−log₁₀ BH adjusted P value")
        ax.text(.02,.03,f"Modeled: {len(z):,} | FDR available: {int(np.isfinite(z.padj).sum()):,}\nBH<0.05: {int(z.padj.lt(.05).sum()):,} | supported: {int(green.sum()):,}",transform=ax.transAxes,fontsize=7.6,va="bottom",bbox=dict(facecolor="white",edgecolor="none",alpha=.85,pad=2))
        label_genes(ax,z,robust[pop]); style(ax)
    text = "Original-full PRIMARY DESeq2 totals, cKO versus WT. Green genes meet the native core sensitivity flag; purple genes meet this PRIMARY endpoint BH<0.05 only. Gray genes have available FDR≥0.05. Missing FDR/effects are counted but not plotted as nonsignificant. Labels rank supported genes by PRIMARY adjusted P value, without a fold-change cutoff. Gene-body RNA does not establish protein loss, exon retention or anatomy."
    caption(fig,"03  Which core RNA effects meet sensitivity checks?",text);fig.subplots_adjust(left=.08,right=.98,top=.87,bottom=.19,hspace=.38,wspace=.25)
    writer.save(fig,"03_core_volcano","Original-full core DESeq2 volcano plots",text,branch="original_full")


def paired_effect_frame(left,right,robust):
    z=left.set_index("gene_id").join(right.set_index("gene_id"),how="outer",lsuffix="_raw",rsuffix="_cb")
    finite=np.isfinite(z.log2FoldChange_raw)&np.isfinite(z.log2FoldChange_cb)
    af=np.isfinite(z.padj_raw);bf=np.isfinite(z.padj_cb);ah=af&z.padj_raw.lt(.05);bh=bf&z.padj_cb.lt(.05)
    missing_hit=(ah&~bf)|(bh&~af)
    z["plot_class"]="Other finite effects"
    z.loc[af&bf&(ah^bh),"plot_class"]="One BH hit; both FDR available"
    z.loc[ah&bh,"plot_class"]="Both BH hits"
    z.loc[missing_hit,"plot_class"]="BH hit; counterpart FDR unavailable"
    z.loc[z.index.isin(robust),"plot_class"]="Core sensitivity-supported"
    z["both_effects_finite"]=finite
    return z


def draw_concordance(plt,writer,effects,robust):
    fig,axes=plt.subplots(2,2,figsize=(11.6,8.8)); metrics=[]
    styles=[("Other finite effects",COL["other"],5,"o"),("One BH hit; both FDR available",COL["single"],9,"o"),
            ("Both BH hits",COL["hit"],10,"o"),("Core sensitivity-supported",COL["robust"],12,"o"),
            ("BH hit; counterpart FDR unavailable",COL["missing"],14,"x")]
    for ax,pop in zip(axes.flat,POPS):
        left,right=effects[(pop,"paired_original")],effects[(pop,"paired_cellbender")]
        if left is None or right is None:ax.text(.5,.5,"Paired model unsupported",ha="center",transform=ax.transAxes);continue
        z=paired_effect_frame(left,right,robust[pop]);finite=z.both_effects_finite
        for name,color,size,marker in styles:
            keep=finite&z.plot_class.eq(name)
            ax.scatter(z.loc[keep,"log2FoldChange_raw"],z.loc[keep,"log2FoldChange_cb"],s=size,c=color,marker=marker,alpha=.7,linewidths=.7 if marker=="x" else 0,rasterized=True,label=name)
        values=np.r_[z.loc[finite,"log2FoldChange_raw"],z.loc[finite,"log2FoldChange_cb"]]
        limit=max(1,float(np.max(np.abs(values))))*1.06 if len(values) else 1
        ax.plot([-limit,limit],[-limit,limit],color="#535B65",linestyle="--",linewidth=.8);ax.set_xlim(-limit,limit);ax.set_ylim(-limit,limit);ax.set_aspect("equal")
        # Descriptive rank correlation only; no new correlation P value.
        aa=z.loc[finite,"log2FoldChange_raw"].to_numpy(float);bb=z.loc[finite,"log2FoldChange_cb"].to_numpy(float)
        rho=float(np.corrcoef(rankdata(aa),rankdata(bb))[0,1]) if len(aa)>2 and np.std(aa)>0 and np.std(bb)>0 else np.nan
        row={"population":pop,"paired_raw_gene_rows":len(left),"paired_cb_gene_rows":len(right),"both_finite_effects":int(finite.sum()),"effect_Spearman":rho,
             "raw_FDR_available":int(left.padj.notna().sum()),"cb_FDR_available":int(right.padj.notna().sum()),
             "raw_BH_hits":int(left.padj.lt(.05).sum()),"cb_BH_hits":int(right.padj.lt(.05).sum()),
             "BH_hit_counterpart_FDR_unavailable_with_finite_effects":int((finite&z.plot_class.eq("BH hit; counterpart FDR unavailable")).sum()),
             "unplottable_effect_rows":int((~finite).sum())};metrics.append(row)
        ax.set_title(NAMES[pop]);ax.set_xlabel("Paired raw: cKO / WT log₂FC");ax.set_ylabel("Paired CellBender: cKO / WT log₂FC")
        ax.text(.025,.97,f"ρ={rho:.3f} | finite effects: {int(finite.sum()):,}\nBH hits raw→CB: {row['raw_BH_hits']:,}→{row['cb_BH_hits']:,}",transform=ax.transAxes,va="top",fontsize=8)
    handles,labels=axes.flat[0].get_legend_handles_labels();fig.legend(handles,labels,loc="lower center",bbox_to_anchor=(.5,.095),ncol=3,frameon=False,fontsize=7.7)
    text="Same paired nuclei, original frozen groups and paired-original gene universe. Every finite effect pair is plotted; the diagonal denotes unchanged magnitude. FDR-unavailable counterparts are explicit crosses, not lost-significance or nonsignificant genes. Green denotes the stricter native core flag across count branches plus partition/leave-one-mouse-out direction checks. Count correction is a sensitivity analysis; source availability and normalization affect uncertainty."
    caption(fig,"04  How much does CellBender change the estimated effect?",text);fig.subplots_adjust(left=.08,right=.98,top=.87,bottom=.23,hspace=.35,wspace=.30)
    writer.save(fig,"04_paired_count_effect_concordance","Paired raw versus CellBender DESeq2 effects",text,branch="paired_original_vs_paired_cellbender")
    return pd.DataFrame(metrics)


def draw_fine(plt,writer,crosswalk,comparisons):
    fig,axes=plt.subplots(2,2,figsize=(11.6,8.4)); saved=[]
    matches=crosswalk[crosswalk.branch.eq("original_full")&crosswalk.analysis_kind.eq("cluster")].copy()
    comps=comparisons[comparisons.comparison_scope.eq("FROZEN_POPULATION_DEFINITION_COMPARISON_SAME_DESEQ2_ENGINE")&comparisons.left_branch.eq("original_full")].copy()
    keys=["population","left_population_key","right_population_key"]
    z=matches.merge(comps,on=keys,how="left",suffixes=("_membership",""),validate="one_to_one")
    for ax,pop in zip(axes.flat,POPS):
        x=z[z.population.eq(pop)].copy(); ok=np.isfinite(x.effect_Spearman) if "effect_Spearman" in x else np.zeros(len(x),bool)
        good=truth(x.fine_population_interpretation_guard_pass.fillna(False)) if "fine_population_interpretation_guard_pass" in x else np.zeros(len(x),bool)
        low=x.get("left_coverage_flag",pd.Series("",index=x.index)).fillna("").str.contains("LOW_NUCLEUS")|x.get("right_coverage_flag",pd.Series("",index=x.index)).fillna("").str.contains("LOW_NUCLEUS")
        for keep,color,marker,label in [(ok&~good&~low,COL["other"],"o","Provisional label/overlap"),(ok&good,COL["robust"],"o","Stronger label/coverage evidence"),(ok&low,COL["single"],"^","IN mouse contribution <5 nuclei")]:
            ax.scatter(x.loc[keep,"jaccard"],x.loc[keep,"effect_Spearman"],s=38,c=color,marker=marker,label=label,alpha=.85,edgecolors="white",linewidths=.4)
        for i,row in enumerate(x.loc[ok].sort_values("effect_Spearman").head(2).itertuples(index=False)):
            ax.annotate(f"Scanpy {row.left_label} / Seurat {row.right_label}",(row.jaccard,row.effect_Spearman),
              xytext=(.04,.45-.10*i),textcoords="axes fraction",fontsize=8,
              bbox=dict(facecolor="white",edgecolor="none",alpha=.9,pad=2),
              arrowprops=dict(arrowstyle="-",color="#8C939B",linewidth=.6))
        ax.axvline(.5,color="#8C939B",linestyle="--",linewidth=.7);ax.set_xlim(0,1.04);ax.set_ylim(-1.02,1.06)
        ax.set_title(NAMES[pop]);ax.set_xlabel("Cluster membership Jaccard");ax.set_ylabel("Gene-effect Spearman ρ")
        ax.text(.03,.03,f"Assigned pairs: {len(x)} | effects compared: {int(ok.sum())}",transform=ax.transAxes,fontsize=8);style(ax)
        saved.append(x)
    handles,labels=axes.flat[0].get_legend_handles_labels();fig.legend(handles,labels,loc="lower center",bbox_to_anchor=(.5,.092),ncol=3,frameon=False,fontsize=8)
    text="Original Harmony0.4 cluster assignments are matched by nucleus overlap. Effect comparisons require ≥20 shared nuclei, ≥50% reciprocal overlap and both DESeq2 models tested. Empty/unmatched/unsupported groups remain in native tables. Green also requires stronger reference/marker evidence and Jaccard≥0.5; it does not validate anatomy. Both definitions use one DESeq2 engine. Fine groups have PRIMARY-only inference, not the core partition/LOO flag. S/R labels denote Scanpy/Seurat cluster IDs."
    caption(fig,"05  Which fine cluster comparisons can be interpreted?",text);fig.subplots_adjust(left=.09,right=.98,top=.87,bottom=.23,hspace=.34,wspace=.26)
    writer.save(fig,"05_matched_fine_effects","Matched original clusters: membership and effect agreement",text,branch="original_full")
    return pd.concat(saved,ignore_index=True)


def draw_reliability(plt,writer,effects,core):
    fig,ax=plt.subplots(figsize=(11.6,5.8));x=np.arange(len(POPS));width=.22
    counts=[]
    for pop in POPS:
        z=effects[(pop,"original_full")];r=core[(core.tool=="scanpy")&core.population.eq(pop)]
        counts.append([int(z.padj.lt(.05).sum()) if z is not None else 0,int(truth(r.PRIMARY_FDR05_all_three_count_branches_same_sign).sum()),int(truth(r.stringent_reliability_flag).sum())])
    counts=np.array(counts)
    for i,(name,color) in enumerate([("Original PRIMARY BH<0.05",COL["hit"]),("All3 PRIMARY BH<0.05; same sign","#4D89BA"),("Plus partition / all12 LOO direction",COL["robust"]) ]):
        ax.bar(x+(i-1)*width,counts[:,i],width*.92,color=color,label=name)
        for px,n in zip(x+(i-1)*width,counts[:,i]):ax.text(px,n+max(1,counts.max())*.015,f"{n:,}",ha="center",va="bottom",fontsize=9)
    ax.set_xticks(x,["All","EX","CA1/ProS","IN"]);ax.set_ylabel("Genes passing the stated core criterion");ax.set_ylim(0,max(1,counts.max())*1.28);style(ax);ax.legend(frameon=False,loc="upper left",fontsize=9)
    text="Native core flags require BH<0.05 and the same nonzero effect direction in all3 PRIMARY count branches, same-direction finite effects in both original-full partition sensitivities, and all12 tested leave-one-mouse-out effects. Partition/LOO checks require estimability and direction, not significance in every refit. Shared totals have identical counts across tools and one DESeq2 implementation. BH remains within each endpoint; this is sensitivity support, not a pooled discovery FDR or causal proof."
    caption(fig,"06  How many core effects pass each sensitivity criterion?",text);fig.subplots_adjust(left=.08,right=.98,top=.85,bottom=.25)
    writer.save(fig,"06_core_sensitivity_support","Core DESeq2 sensitivity support",text)


def draw_mouse_genes(plt,writer,mouse_values,selections):
    fig,axes=plt.subplots(4,2,figsize=(11.6,12.5))
    for row,pop in enumerate(POPS):
        genes=selections[selections.population.eq(pop)].sort_values("selection_rank").head(2)
        for col,ax in enumerate(axes[row]):
            if col>=len(genes):
                ax.text(.5,.5,"No additional gene meets\nthe core sensitivity flag",transform=ax.transAxes,ha="center",va="center",color="#6A717A");ax.set_axis_off();continue
            gene=genes.iloc[col];z=mouse_values[(mouse_values.population==pop)&mouse_values.gene_id.eq(gene.gene_id)]
            for j,condition in enumerate(["WT","cKO"]):
                m=z[z.condition.eq(condition)&np.isfinite(z.normalized_aggregate_UMIs)].copy();jitter=np.linspace(-.32,.32,len(m)) if len(m)>1 else np.array([0.])
                values=np.log10(1+m.normalized_aggregate_UMIs.to_numpy(float));ax.scatter(j+jitter,values,s=38,c=COL[condition],edgecolors="white",linewidths=.5,zorder=3)
                for k,(xx,yy,sample) in enumerate(zip(j+jitter,values,m["sample"])):
                    ax.annotate(sample.split("_")[-1],(xx,yy),xytext=(0,6 if k%2==0 else -11),textcoords="offset points",ha="center",fontsize=6.8,color="#48515B",
                      bbox=dict(facecolor="white",edgecolor="none",alpha=.8,pad=.2),zorder=5)
                if len(m):ax.plot([j-.35,j+.35],[np.log10(1+m.normalized_aggregate_UMIs.mean())]*2,color="#26323E",linewidth=1.7,zorder=2)
            ax.set_xticks([0,1],["WT","cKO"]);ax.set_xlim(-.45,1.45);ax.set_ylabel("log₁₀(1 + DESeq2 normalized UMIs)",fontsize=8.4)
            ax.set_title(f"{NAMES[pop]} · {gene.gene}\nPRIMARY log₂FC={gene.log2FoldChange:.2f}; BH={gene.padj:.2g}",fontsize=10.3);style(ax)
    text="Top2 core sensitivity-supported genes per population, ranked by original-full PRIMARY BH adjusted P value; no fold-change cutoff. Points are independent mice, labeled by mouse ID suffix. Existing full-gene aggregate UMIs are divided by that endpoint's recorded fitted size factor, then displayed as log10(1+x). Horizontal lines are transformed arithmetic genotype means of normalized counts, not a new test or confidence interval. No nuclei are pooled as replicates. These RNA measurements do not establish Satb2 exon retention/protein loss or a physical layer shift."
    caption(fig,"07  What do the independent mice show?",text);fig.subplots_adjust(left=.095,right=.98,top=.93,bottom=.14,hspace=.58,wspace=.29)
    writer.save(fig,"07_core_supported_gene_mouse_points","Independent mice: top core-supported RNA effects",text,branch="original_full",scope="existing_integer_gene_by_mouse_counts_divided_by_recorded_size_factors")


def review_all_primary_exports(e, data, module):
    """Read the completed endpoints and policies; never aggregate or fit counts."""
    records=[]
    for branch in BRANCHES:
        for row in data.indices[branch].itertuples(index=False):
            # The native reader validates the seven fields, stable-ID row binding,
            # settings, mouse eligibility and every reported result denominator.
            z,ad=data.endpoint(branch,row)
            ep=data.model_roots[branch]/row.tool/row.population_key/"PRIMARY"
            mice=e.table(ep/"AUDIT/MICE_ALL12_AUDIT.tsv",dtype={"sample":str,"condition":str})
            if len(mice)!=12 or mice["sample"].duplicated().any():raise ValueError("Incomplete all12 coverage audit")
            if set(mice["sample"])!=set(data.rosters[branch]) or (mice.condition!=mice["sample"].map(data.rosters[branch])).any():raise ValueError("All12 sample/genotype identities differ")
            included=module.truth(mice.included_in_model);selected=mice[included]
            minimum=1 if row.population=="Inhibitory" else 20
            if not np.array_equal(included,mice.n_nuclei.ge(minimum)&mice.total_UMIs.gt(0)):raise ValueError("PRIMARY mouse inclusion differs from declared full-coverage policy")
            lo=int(selected.n_nuclei.min()) if len(selected) else np.nan
            hi=int(selected.n_nuclei.max()) if len(selected) else np.nan
            if len(selected) and ((selected.n_nuclei<minimum).any() or (selected.total_UMIs<=0).any()):raise ValueError("Recorded included mouse violates declared coverage")
            if ad["status"]=="TESTED_REVIEW_REQUIRED":
                if lo!=ad["min_included_mouse_nuclei"] or hi!=ad["max_included_mouse_nuclei"]:raise ValueError("Included-mouse range disagrees with endpoint audit")
                low=minimum==1 and lo<5
                if low!=(str(ad.get("coverage_flag",""))=="LOW_NUCLEUS_COVERAGE_EXPLORATORY"):raise ValueError("Low-nucleus IN coverage flag mismatch")
            pf=e.table(ep.parent/"AUDIT/FROZEN_PRIMARY_GENE_PREFILTER_AUDIT.tsv",dtype={"gene_id":str})
            if pf.gene_id.duplicated().any():raise ValueError("Duplicated frozen gene ID")
            cutoff=5 if minimum==1 else 10
            required=3 if minimum==1 else min(int(selected.condition.eq("WT").sum()),int(selected.condition.eq("cKO").sum()))
            if pf.count_cutoff_per_mouse.ne(cutoff).any() or pf.required_mice.ne(required).any():raise ValueError("Frozen gene prefilter policy differs")
            universe=set(pf.loc[module.truth(pf.passes_primary_prefilter),"gene_id"])
            if len(universe)!=int(ad["n_frozen_primary_genes"]):raise ValueError("Frozen primary universe count differs")
            if branch=="paired_cellbender":
                raw=data.families["paired_original"].get(module.identity(pd.Series(row._asdict())))
                if raw is None or universe!=data.universe("paired_original",raw):raise ValueError("Corrected branch changed paired-original frozen universe")
            if z is not None and not set(z.gene_id).issubset(universe):raise ValueError("Modeled gene outside frozen source universe")
            warnings_path=ep/"AUDIT/MODEL_WARNINGS.txt";warnings=""
            if warnings_path.exists():warnings=e.bind(warnings_path).read_text().strip()
            fit_type_path=ep/"AUDIT/DISPERSION_FIT_TYPE.txt"
            fit_type=e.bind(fit_type_path).read_text().strip() if fit_type_path.exists() else ""
            records.append({"branch":branch,"tool":row.tool,"population":row.population,"analysis_kind":row.analysis_kind,
              "analysis_label":row.analysis_label,"population_key":row.population_key,"status":ad["status"],"reason":ad.get("reason",""),
              "n_WT":ad["n_WT"],"n_cKO":ad["n_cKO"],"residual_df":ad["residual_df"],
              "min_nuclei_policy":minimum,"min_included_mouse_nuclei":lo,"max_included_mouse_nuclei":hi,
              "included_mice_below5_nuclei":int(selected.n_nuclei.lt(5).sum()),"zero_nucleus_mice":int(mice.n_nuclei.eq(0).sum()),
              "n_frozen_genes":len(universe),"n_modeled_genes":len(z) if z is not None else 0,
              "n_pvalue_available":int(z.pvalue.notna().sum()) if z is not None else 0,
              "n_FDR_available":int(z.padj.notna().sum()) if z is not None else 0,"n_BH05":int(z.padj.lt(.05).sum()) if z is not None else 0,
              "missing_symbol_rows":int(z.display_symbol.isna().sum()) if z is not None else 0,
              "rows_with_duplicated_symbol":int(z.display_symbol.notna().where(z.display_symbol.duplicated(keep=False),False).sum()) if z is not None else 0,
              "coverage_flag":ad.get("coverage_flag",""),"warning_count":ad.get("warning_count",0),"recorded_warning_text":warnings,
              "actual_dispersion_fit_type":fit_type,
              "execution_origin":ad.get("execution_origin","UNKNOWN"),"cache_reason":ad.get("cache_reason",""),
              "seven_symbol_export_verified":z is not None,"result_path":str(ep/"ALL_RESULTS.txt") if z is not None else "",
              "endpoint_audit_path":str(ep/"AUDIT/EXECUTION_AUDIT.tsv"),"model_enabled":row.model_enabled})
    return pd.DataFrame(records)


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument("--config",required=True,type=Path);ap.add_argument("--outdir",required=True,type=Path);ap.add_argument("--check",action="store_true");args=ap.parse_args()
    cfg=json.loads(args.config.read_text());state=readiness(cfg)
    if args.check:print(json.dumps(state,indent=2));return
    if not state["ready"]:raise RuntimeError("Required native producers are not complete: "+str(state["missing_producer_receipts"]))
    root=Path(cfg["results_root"]);cmp=root/"comparison";receipt=json.loads((cmp/"COMPLETE.json").read_text())
    if receipt.get("status")!="EXECUTED_REVIEW_REQUIRED" or receipt.get("stage")!="final_produced_artifact_comparison" or receipt.get("exact_integer_counts_PASS") is not True or receipt.get("same_shared_total_counts_PASS") is not True or receipt.get("config_sha256")!=sha(args.config):raise ValueError("Final native acceptance/config binding failed")
    module_path=Path(cfg["study_root"])/"source/diagnostics/compare_revision.py"
    if sha(module_path)!=receipt.get("source_code_sha256"):raise ValueError("Native audit reader source differs from the completed comparator")
    spec=importlib.util.spec_from_file_location("satb2_completed_result_reader",module_path);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    class RoundTripEvidence(module.Evidence):
        def table(self,path,**kwargs):
            kwargs.setdefault("float_precision","round_trip")
            return super().table(path,**kwargs)
    e=RoundTripEvidence();e.bind(args.config);e.bind(module_path);e.bind(__file__);e.stage(cmp,"EXECUTED_REVIEW_REQUIRED","final_produced_artifact_comparison")
    data=module.Results(e,cfg)
    core=e.table(cmp/"CORE_WT_cKO_RELIABILITY_ALL_GENES.tsv.gz",dtype={"gene_id":str,"tool":str,"population":str})
    comparisons=e.table(cmp/"EFFECT_COMPARISONS.tsv");crosswalk=e.table(cmp/"FINE_POPULATION_CROSSWALK.tsv",dtype={"left_label":str,"right_label":str})
    status=[];effects={};coverage={};robust={};selections=[];mouse_values=[]
    for branch in BRANCHES:
        audit=e.table(data.model_roots[branch]/"PSEUDOBULK_EXECUTION_AUDIT.tsv",dtype={"population_key":str})
        primary=audit[audit.endpoint.eq("PRIMARY")]
        if primary.duplicated(["tool","population_key"]).any() or len(primary)!=len(data.indices[branch]) or primary.status.eq("ERROR").any():raise ValueError("Incomplete/failed all-population primary status registry")
        status.append(primary)
    primary_review=review_all_primary_exports(e,data,module)
    for pop in POPS:
        for branch in BRANCHES:
            ix=data.indices[branch];rows=ix[ix.population.eq(pop)&ix.analysis_kind.eq("total")];a=rows[rows.tool.eq("scanpy")];b=rows[rows.tool.eq("seurat")]
            if len(a)!=1 or len(b)!=1:raise ValueError("Missing unique shared total")
            z,ad=data.endpoint(branch,a.iloc[0]);zz,aad=data.endpoint(branch,b.iloc[0])
            if (z is None)!=(zz is None):raise ValueError("Shared total model support differs")
            if z is not None:
                if not np.array_equal(z.gene_id,zz.gene_id) or not np.allclose(z[STATS],zz[STATS],rtol=0,atol=0,equal_nan=True):raise ValueError("Shared total effects unexpectedly differ")
            effects[(pop,branch)]=z
            coverage[(pop,branch)]=e.table(Path(a.iloc[0].counts_file).parent/"ALL12_POPULATION_COVERAGE.tsv",dtype={"sample":str,"condition":str})
        rr=core[(core.tool=="scanpy")&core.population.eq(pop)]
        robust[pop]=set(rr.loc[truth(rr.stringent_reliability_flag),"gene_id"])
        z=effects[(pop,"original_full")]
        if z is None:continue
        chosen=z[z.gene_id.isin(robust[pop])&z.display_symbol.notna()].sort_values(["padj","gene_id"]).head(2)
        ix=data.indices["original_full"];row=ix[ix.tool.eq("scanpy")&ix.population.eq(pop)&ix.analysis_kind.eq("total")].iloc[0]
        count_path=Path(row.counts_file)
        if sha(e.bind(count_path))!=row.counts_sha256:raise ValueError("Normalized display source count array changed")
        counts=e.table(count_path,dtype={"gene_id":str}).set_index("gene_id")
        ep=data.model_roots["original_full"]/row.tool/row.population_key/"PRIMARY"
        sf=e.table(ep/"AUDIT/FITTED_SIZE_FACTORS.tsv",dtype={"sample":str}).set_index("sample").fitted_size_factor
        mice=e.table(ep/"AUDIT/MICE_ALL12_AUDIT.tsv",dtype={"sample":str,"condition":str})
        if sf.index.duplicated().any() or not np.isfinite(sf).all() or (sf<=0).any():raise ValueError("Invalid existing fitted normalization factors")
        for rank,gene in enumerate(chosen.itertuples(index=False),1):
            selections.append({"population":pop,"gene_id":gene.gene_id,"gene":gene.display_symbol,"selection_rank":rank,"log2FoldChange":gene.log2FoldChange,"padj":gene.padj,"selection_scope":"native_core_stringent_flag;original_full_primary_BH_rank"})
            for mouse in mice.itertuples(index=False):
                value=float(counts.loc[gene.gene_id,mouse.sample]);factor=float(sf.loc[mouse.sample]) if mouse.sample in sf.index else np.nan
                if not np.isfinite(value) or value<0 or value!=np.floor(value):raise ValueError("Noninteger existing aggregate gene count")
                mouse_values.append({"population":pop,"gene_id":gene.gene_id,"gene":gene.display_symbol,"sample":mouse.sample,"condition":mouse.condition,"n_nuclei":mouse.n_nuclei,"included_in_model":mouse.included_in_model,"source_integer_aggregate_UMIs":value,"fitted_size_factor":factor,"normalized_aggregate_UMIs":value/factor if np.isfinite(factor) else np.nan})
    e.finish()
    out=args.outdir.resolve()
    if out.exists() and any(out.iterdir()):raise ValueError("Refusing nonempty supplement destination")
    out.mkdir(parents=True,exist_ok=True);audit_dir=out/"audit_data";audit_dir.mkdir()
    plt,PdfPages=setup_plotting();writer=FigureWriter(out,plt,PdfPages)
    audits=pd.concat(status,ignore_index=True);selected=pd.DataFrame(selections,columns=["population","gene_id","gene","selection_rank","log2FoldChange","padj","selection_scope"])
    mice=pd.DataFrame(mouse_values,columns=["population","gene_id","gene","sample","condition","n_nuclei","included_in_model","source_integer_aggregate_UMIs","fitted_size_factor","normalized_aggregate_UMIs"])
    draw_status(plt,writer,audits);draw_coverage(plt,writer,coverage);draw_volcano(plt,writer,effects,robust)
    metrics=draw_concordance(plt,writer,effects,robust);fine=draw_fine(plt,writer,crosswalk,comparisons);draw_reliability(plt,writer,effects,core);draw_mouse_genes(plt,writer,mice,selected);writer.close()
    for name,frame in [("ALL_PRIMARY_MODEL_STATUS.tsv",audits),("ALL_PRIMARY_EXPORT_AND_COVERAGE_REVIEW.tsv",primary_review),("PAIRED_CORE_EFFECT_METRICS.tsv",metrics),("MATCHED_FINE_CLUSTER_DISPLAY_ROWS.tsv",fine),("TOP_SUPPORTED_GENE_SELECTION.tsv",selected),("NORMALIZED_MOUSE_VALUES.tsv",mice)]:frame.to_csv(audit_dir/name,sep="\t",index=False,na_rep="NA")
    e.finish();outputs={str(p.relative_to(out)):sha(p) for p in out.rglob("*") if p.is_file()}
    provenance={"status":"COMPLETE","stage":"completed_DESeq2_result_figure_supplement","config_sha256":sha(args.config),"source_code_sha256":sha(__file__),"input_sha256":e.hashes,"output_sha256":outputs,"scientific_fitting_performed":False,"raw_nucleus_count_aggregation_performed":False,"native_outputs_modified":False,"core_representative_tool":"scanpy_after_exact_shared_total_count_and_effect_verification","core_reliability_definition":"All3 PRIMARY BH<0.05 same sign; both original partition effects finite same sign; all12 tested LOO effects finite same sign. No requirement of significance in every sensitivity endpoint.","normalization_display":"existing integer gene-by-mouse aggregate / recorded endpoint fitted size factor; no refit","n_figures":len(writer.rows),"pvalue_policy":"no new tests; missing FDR retained unavailable","physical_layer_validation":False}
    (out/"COMPLETE.json").write_text(json.dumps(provenance,indent=2,allow_nan=False)+"\n")
    print(json.dumps({"status":"COMPLETE","n_figures":len(writer.rows),"figure_index":str(out/"FIGURE_INDEX.tsv")},indent=2))


if __name__=="__main__":main()
