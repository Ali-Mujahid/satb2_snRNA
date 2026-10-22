#!/usr/bin/env python3
"""Supplemental circular ligand–receptor complex networks from saved evidence.

No models or scientific workflows are executed. Each node is a unique
(population, complex, ligand/receptor role), not a gene collapsed across groups.
The top interaction identities are selected jointly across WT, cKO and pooled
contexts, then displayed with common node order and edge/node scales.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, Wedge
import numpy as np
import pandas as pd

CONDITIONS = ["WT", "cKO", "pooled"]
DIRECTIONS = ["CA1_to_IN", "IN_to_CA1", "bidirectional"]
KEYS = ["source", "target", "ligand_complex", "receptor_complex"]
MIN_COLUMNS = ["ligand_all_subunit_min_prop", "receptor_all_subunit_min_prop"]
BLUE, ORANGE, GRAY = "#2166ac", "#d46a23", "#8a96a2"


def cli():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--evidence", type=Path, required=True,
                   help="Completed interaction_evidence.tsv[.gz]")
    p.add_argument("--outdir", type=Path, required=True,
                   help="A new supplemental figure directory")
    p.add_argument("--top", type=int, default=12,
                   help="Maximum interaction identities per scope/direction; 1–20 (default 12)")
    p.add_argument("--rank-cutoff", type=float, default=.05)
    p.add_argument("--dpi", type=int, default=200)
    return p.parse_args()


def natural(value):
    return [int(x) if x.isdigit() else x.lower() for x in re.split(r"(\d+)", str(value))]


def group_label(name):
    return {"CA1_Deep":"CA1 Deep", "CA1_Superficial":"CA1 Superficial",
            "IN_Total":"IN Total"}.get(name, str(name).replace("IN_", "IN "))


def group_color(name):
    if name == "CA1_Deep":
        return "#875636"
    if name == "CA1_Superficial":
        return "#397fb0"
    if name == "IN_Total":
        return "#bb6b25"
    match = re.search(r"C(\d+)", name)
    i = int(match.group(1)) if match else 0
    return matplotlib.colors.to_hex(plt.colormaps["tab20"]((i-1) % 20))


def load(path):
    df = pd.read_csv(path, sep="\t", low_memory=False)
    required = {"scope", "condition", "direction", "deg_supported", "magnitude_rank",
                "lr_gmeans", *KEYS, *MIN_COLUMNS}
    if required - set(df):
        raise ValueError("Missing columns: " + ", ".join(sorted(required - set(df))))
    boolean = df.deg_supported.astype(str).str.lower().str.strip()
    if not boolean.isin(["true", "false", "1", "0", "t", "f", "yes", "no"]).all():
        raise ValueError("deg_supported must be explicitly boolean.")
    df["deg_supported"] = boolean.isin(["true", "1", "t", "yes"])
    df["condition"] = df.condition.astype(str).str.lower().map({"wt":"WT", "cko":"cKO", "pooled":"pooled"})
    if df.condition.isna().any():
        raise ValueError("Expected WT, cKO and pooled condition labels.")
    if not df.scope.isin(["IN_clusters", "IN_total"]).all():
        raise ValueError("Expected IN_clusters or IN_total scope.")
    if not df.direction.isin(["CA1_to_IN", "IN_to_CA1"]).all():
        raise ValueError("Unexpected direction in evidence table.")
    for col in ["magnitude_rank", "lr_gmeans", *MIN_COLUMNS]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    for col in ["magnitude_rank", *MIN_COLUMNS]:
        if not df[col].dropna().between(-1e-9, 1+1e-9).all():
            raise ValueError(f"{col} must be 0–1, or NA.")
    if not df.lr_gmeans.dropna().ge(0).all():
        raise ValueError("Native geometric-mean scores must be nonnegative.")
    if df.duplicated(["scope", "condition", *KEYS]).any():
        raise ValueError("Duplicate interaction identity within scope and condition.")
    return df


def select_keys(base, n, cutoff):
    candidates = base.loc[base.deg_supported & base.magnitude_rank.le(cutoff)].copy()
    ranked = candidates.groupby(KEYS, as_index=False).agg(
        selection_best_magnitude_rank=("magnitude_rank","min"),
        selection_max_lr_gmeans=("lr_gmeans","max"))
    ranked = ranked.sort_values(["selection_best_magnitude_rank", "selection_max_lr_gmeans", *KEYS],
        ascending=[True, False, True, True, True, True], kind="stable").head(n).reset_index(drop=True)
    ranked.insert(0,"selection_order",range(1,len(ranked)+1))
    return ranked


def node_table(selected):
    unique = set()
    for row in selected.itertuples():
        unique.add((row.source, row.ligand_complex, "ligand"))
        unique.add((row.target, row.receptor_complex, "receptor"))
    ordered = sorted(unique, key=lambda x:(not x[0].startswith("CA1_"), natural(x[0]), x[2], natural(x[1])))
    if len(ordered) > 40:
        raise ValueError("More than 40 endpoint-role nodes; reduce --top.")
    df = pd.DataFrame(ordered,columns=["group","complex","role"])
    df.insert(0,"node_id",[f"N{i+1:03d}" for i in range(len(df))])
    return df


def positions(nodes):
    """Group-block sectors and identical node positions across all conditions."""
    sectors, coordinates = [], {}
    groups = list(dict.fromkeys(nodes.group))
    gap = math.radians(4)
    unit = (2*math.pi-gap*len(groups))/max(1,len(nodes))
    cursor = math.pi/2
    for group in groups:
        sub = nodes.loc[nodes.group.eq(group)]
        start, end = cursor, cursor-unit*len(sub)
        sectors.append(dict(group=group,start_radians=start,end_radians=end))
        for i, row in enumerate(sub.itertuples()):
            angle = cursor-unit*(i+.5)
            coordinates[row.node_id] = dict(angle=angle,x=math.cos(angle),y=math.sin(angle))
        cursor=end-gap
    return coordinates,sectors


def review_needed(value):
    return str(value).strip().lower() not in ["", "nan", "none", "false", "0", "reference_supported", "resolved", "ok", "pass"]


def node_values(nodes, context, coordinates):
    values=[]
    for row in nodes.itertuples():
        endpoint="source" if row.role=="ligand" else "target"
        molecule="ligand_complex" if row.role=="ligand" else "receptor_complex"
        fraction=MIN_COLUMNS[0] if row.role=="ligand" else MIN_COLUMNS[1]
        native=context.loc[context[endpoint].eq(row.group)&context[molecule].eq(row.complex)]
        props=native[fraction].dropna()
        # Identical complex/group contexts should have identical fractions;
        # a conservative minimum plus range retains any unexpected discrepancy.
        prop=float(props.min()) if len(props) else np.nan
        review=any(review_needed(x) for x in native.get(endpoint+"_review",pd.Series(dtype=str)))
        annotations=native.get(endpoint+"_label",pd.Series(dtype=str)).dropna().astype(str).unique().tolist()
        values.append(dict(node_id=row.node_id,group=row.group,complex=row.complex,role=row.role,
            detection_fraction=prop,detection_fraction_max=float(props.max()) if len(props) else np.nan,
            detection_definition="minimum across all subunits of this complex in this population",
            plot_area_points2=25+340*prop if np.isfinite(prop) else 25,
            annotation="; ".join(annotations),review=review,color=group_color(row.group),
            **coordinates[row.node_id]))
    columns=["node_id","group","complex","role","detection_fraction","detection_fraction_max",
        "detection_definition","plot_area_points2","annotation","review","color","angle","x","y"]
    return pd.DataFrame(values,columns=columns)


def complex_label(name):
    # A newline between complex subunits preserves every symbol rather than
    # abbreviating a different biological complex into a similar label.
    parts=str(name).split("_")
    return "\n".join(parts) if len(parts)>1 else str(name)


def render(base, selected, nodes, coordinates, sectors, scope, direction, condition, folder, a, score_max):
    context=base.loc[base.condition.eq(condition)].copy()
    plotted=context.merge(selected,on=KEYS,how="inner",validate="one_to_one")
    plotted=plotted.loc[plotted.deg_supported & plotted.magnitude_rank.le(a.rank_cutoff)].copy()
    lookup={(r.group,r.complex,r.role):r.node_id for r in nodes.itertuples()}
    if len(plotted):
        plotted["source_node_id"]=[lookup[(s,l,"ligand")] for s,l in zip(plotted.source,plotted.ligand_complex)]
        plotted["target_node_id"]=[lookup[(t,r,"receptor")] for t,r in zip(plotted.target,plotted.receptor_complex)]
    else:
        plotted["source_node_id"]=pd.Series(dtype=str)
        plotted["target_node_id"]=pd.Series(dtype=str)
    plotted["plot_edge_width_points"]=.5+5*plotted.lr_gmeans.fillna(0)/max(score_max,1e-12)
    plotted["plot_edge_color"]=[BLUE if str(x).startswith("CA1_") else ORANGE for x in plotted.source]
    nodedata=node_values(nodes,context,coordinates)
    fig,ax=plt.subplots(figsize=(14,14))
    ax.set_aspect("equal")
    for sector in sectors:
        ax.add_patch(Wedge((0,0),1.035,math.degrees(sector["end_radians"]),math.degrees(sector["start_radians"]),
            width=.035,facecolor=group_color(sector["group"]),edgecolor="white",linewidth=.8,zorder=2))
        mid=(sector["start_radians"]+sector["end_radians"])/2
        flagged=bool(len(nodedata) and nodedata.loc[nodedata.group.eq(sector["group"]),"review"].any())
        ax.text(1.55*math.cos(mid),1.55*math.sin(mid),group_label(sector["group"])+(" *" if flagged else ""),
            ha="center",va="center",fontsize=10,fontweight="bold",color=group_color(sector["group"]))
    for row in plotted.sort_values("lr_gmeans",na_position="first").itertuples():
        first,last=coordinates[row.source_node_id],coordinates[row.target_node_id]
        ax.add_patch(FancyArrowPatch((.96*first["x"],.96*first["y"]),(.96*last["x"],.96*last["y"]),
            arrowstyle="-|>",connectionstyle="arc3,rad=0.12",linewidth=row.plot_edge_width_points,
            color=row.plot_edge_color,alpha=.56,mutation_scale=10+row.plot_edge_width_points,
            shrinkA=5,shrinkB=8,zorder=1))
    for row in nodedata.itertuples():
        available=np.isfinite(row.detection_fraction)
        ax.scatter(row.x,row.y,s=row.plot_area_points2,marker="o" if row.role=="ligand" else "s",
            color=row.color if available else "white",edgecolor="white" if available else GRAY,
            linewidth=1,zorder=4)
        deg=math.degrees(row.angle)%360
        rotate=deg if not 90<deg<270 else deg+180
        ha="left" if not 90<deg<270 else "right"
        ax.text(1.07*row.x,1.07*row.y,complex_label(row.complex),fontsize=8.5,ha=ha,va="center",
            rotation=rotate,rotation_mode="anchor",color="#243746")
    if plotted.empty:
        ax.text(0,0,"No selected interactions pass the display criteria\nin this expression context.",
            ha="center",va="center",fontsize=13,color="#536372")
    direct={"CA1_to_IN":"CA1 → inhibitory","IN_to_CA1":"Inhibitory → CA1","bidirectional":"CA1 ↔ inhibitory"}[direction]
    resolution="inhibitory clusters" if scope=="IN_clusters" else "inhibitory total"
    ax.set_title(f"{direct} | {condition} | {resolution}\nLigand–receptor complex circle network",loc="left",fontsize=18,pad=20)
    legends=[Line2D([0],[0],color=BLUE,lw=3,label="CA1 ligand source"),Line2D([0],[0],color=ORANGE,lw=3,label="Inhibitory ligand source"),
        Line2D([0],[0],color=GRAY,marker="o",ls="",markersize=9,label="Ligand complex"),Line2D([0],[0],color=GRAY,marker="s",ls="",markersize=9,label="Receptor complex")]
    ax.legend(handles=legends,loc="lower center",bbox_to_anchor=(.5,-.06),ncol=2,frameon=False,fontsize=10)
    ax.set_xlim(-1.85,1.85)
    ax.set_ylim(-1.8,1.8)
    ax.axis("off")
    footer=(f"Fixed top {len(selected)} interaction identities across WT, cKO and pooled; endpoint-specific DESeq2 support; magnitude rank ≤ {a.rank_cutoff:g}.\n"
        "Node area: minimum detection fraction across all complex subunits (shared scale). Hollow: fraction unavailable in this context.\n"
        "Arrow width: native LIANA geometric-mean expression score (shared scale); arrows identify source → receptor recipient, not activation/inhibition.\n"
        "Predicted transcript-based candidates; ranks are not genotype P values. CA1 Deep/Superficial are expression programs, not physical layers. * Identity review flag.")
    fig.text(.035,.018,footer,fontsize=9,color="#536372",linespacing=1.5)
    fig.subplots_adjust(left=.08,right=.92,top=.90,bottom=.17)
    folder.mkdir(parents=True,exist_ok=True)
    prefix=folder/"complex_circle"
    for ext in ["pdf","png","svg"]:
        fig.savefig(prefix.with_suffix("."+ext),dpi=a.dpi,bbox_inches="tight",facecolor="white")
    plt.close(fig)
    plotted.to_csv(folder/"complex_circle_SourceData.tsv",sep="\t",index=False,na_rep="NA",float_format="%.10g")
    nodedata.to_csv(folder/"complex_circle_Nodes.tsv",sep="\t",index=False,na_rep="NA",float_format="%.10g")
    pd.DataFrame(sectors).to_csv(folder/"complex_circle_Sectors.tsv",sep="\t",index=False,float_format="%.10g")
    return dict(scope=scope,direction=direction,condition=condition,prefix=prefix.relative_to(a.outdir).as_posix(),
        edges_plotted=len(plotted),fixed_interaction_identities=len(selected),nodes=len(nodes),edge_score_max=score_max)


def write_index(manifest,outdir):
    data=json.dumps(manifest,ensure_ascii=False).replace("</","<\\/")
    template=r'''<!DOCTYPE html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Supplemental complex circle networks</title><style>body{font:16px system-ui;color:#213746;background:#edf3f6;max-width:1400px;margin:30px auto;padding:0 20px}section{background:white;padding:20px;border:1px solid #cad9e3;border-radius:12px;margin:18px 0}select{font:inherit;padding:9px;margin:5px 15px 5px 5px}img{width:100%;height:auto}a{color:#155a8b}#downloads a{display:inline-block;margin:10px;padding:8px;border:1px solid #b6cbd9;border-radius:6px}p{line-height:1.5}.notice{border-left:5px solid #c28433;background:#fffaf1}</style><h1>Supplemental ligand–receptor complex circle networks</h1><section class="notice"><p>Nodes are distinct <strong>population × complex × ligand/receptor role</strong> endpoints. These figures prioritize saved transcript-based candidates; they do not establish anatomical contact, signaling activity, or a significant genotype difference. CA1 Deep/Superficial are expression programs, not measured physical layers.</p></section><section><label>Scope <select id="scope"><option>IN_clusters</option><option>IN_total</option></select></label><label>Direction <select id="direction"><option value="CA1_to_IN">CA1 → inhibitory</option><option value="IN_to_CA1">Inhibitory → CA1</option><option value="bidirectional">Both directions</option></select></label><label>Condition <select id="condition"><option>WT</option><option>cKO</option><option value="pooled">Pooled</option></select></label><p id="status"></p><div id="downloads"></div><img id="plot" alt="Complex-level circular ligand–receptor network"></section><section><h2>Reading guide</h2><p>The same top interaction identities, positions and scales are used across conditions. The default maximum is 12 identities, chosen by their best magnitude rank across WT, cKO and pooled contexts; the condition-specific panel shows those passing the rank criterion and endpoint DEG-support filter. This union-based selection is descriptive, not a differential communication test.</p><p>Outer sectors identify populations. Circles are ligand complexes; squares are receptor complexes. Node area is the minimum detection fraction over all subunits of that complex in that population; a hollow node indicates unavailable fraction data. Blue arrows originate in CA1 and orange arrows originate in inhibitory neurons. Arrow widths show native geometric-mean expression scores on a common scale, not P values or numbers of cells. An asterisk retains an identity-review flag. Full source/target annotations are in the node data.</p><p>All figures have PDF, PNG and SVG downloads; exact interaction rows, node values and positions are downloadable separately. Complex subunits remain joined in the tabular identifiers and appear on separate lines in labels.</p><p><a href="manifest.json" download>Manifest and input SHA-256</a> · <a href="README.md">Reproducibility guide</a></p></section><script id="data" type="application/json">__DATA__</script><script>const m=JSON.parse(document.getElementById('data').textContent),ids=['scope','direction','condition'];function draw(){let row=m.figures.find(r=>ids.every(k=>r[k]===document.getElementById(k).value)),image=document.getElementById('plot'),links=document.getElementById('downloads');links.replaceChildren();image.hidden=!row;if(!row){document.getElementById('status').textContent='No panel for this selection.';return;}image.src=row.prefix+'.png';document.getElementById('status').textContent=`${row.edges_plotted} displayed interactions · ${row.nodes} endpoint-role nodes · fixed selection across conditions`;for(const [suffix,title] of [['.pdf','PDF'],['.png','PNG'],['.svg','SVG'],['_SourceData.tsv','Exact interaction data'],['_Nodes.tsv','Node detection and positions'],['_Sectors.tsv','Group sectors']]){let a=document.createElement('a');a.href=row.prefix+suffix;a.download='';a.textContent=title;links.append(a);}}ids.forEach(k=>document.getElementById(k).addEventListener('change',draw));draw();</script></html>'''
    (outdir/"index.html").write_text(template.replace("__DATA__",data),encoding="utf-8")


def main():
    a=cli()
    a.evidence=a.evidence.resolve()
    a.outdir=a.outdir.resolve()
    if not 1<=a.top<=20 or not 0<a.rank_cutoff<=1:
        raise ValueError("Use --top from 1–20 and rank cutoff in (0,1].")
    if (a.outdir/"manifest.json").exists():
        raise FileExistsError("Refusing to overwrite a completed supplemental report.")
    data=load(a.evidence)
    a.outdir.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({"font.family":"DejaVu Sans","pdf.fonttype":42,"ps.fonttype":42,"svg.fonttype":"none"})
    records=[]
    for scope in [x for x in ["IN_clusters","IN_total"] if x in data.scope.unique()]:
        for direction in DIRECTIONS:
            base=data.loc[data.scope.eq(scope)].copy()
            if direction!="bidirectional":
                base=base.loc[base.direction.eq(direction)]
            selected=select_keys(base,a.top,a.rank_cutoff)
            nodes=node_table(selected)
            coordinates,sectors=positions(nodes)
            allselected=base.merge(selected,on=KEYS,how="inner",validate="many_to_one")
            finite_scores=allselected.loc[allselected.deg_supported & allselected.magnitude_rank.le(a.rank_cutoff),"lr_gmeans"].dropna()
            score_max=float(finite_scores.max()) if len(finite_scores) else 1.
            folder=a.outdir/scope/direction
            folder.mkdir(parents=True,exist_ok=True)
            selected.to_csv(folder/"fixed_selected_interaction_identities.tsv",sep="\t",index=False,float_format="%.10g")
            for condition in CONDITIONS:
                records.append(render(base,selected,nodes,coordinates,sectors,scope,direction,condition,folder/condition,a,score_max))
                print(f"Rendered supplemental {scope} / {direction} / {condition}",flush=True)
    manifest=dict(input_file=str(a.evidence),input_sha256=hashlib.sha256(a.evidence.read_bytes()).hexdigest(),
        input_rows=len(data),maximum_interaction_identities=a.top,maximum_endpoint_nodes=40,rank_cutoff=a.rank_cutoff,
        selection="Minimum magnitude rank across all conditions within scope/direction; ties by maximum geometric mean then identity; fixed selection in every condition.",
        methods="Pure rendering of saved tables. No fits, permutation tests, DESeq2 reruns or differential interaction statistics.",figures=records)
    (a.outdir/"manifest.json").write_text(json.dumps(manifest,indent=2,ensure_ascii=False),encoding="utf-8")
    write_index(manifest,a.outdir)
    (a.outdir/"README.md").write_text("# Supplemental complex circle networks\n\n"
        "Open index.html locally; keep all subfolders together when copying the report. No server or internet is required.\n\n"
        "This is a pure visualization of completed LIANA evidence. The registered primary workflow and previous analyses are unchanged.\n\n"
        "Nodes are population × ligand/receptor complex × role. At most 20 fixed interaction identities produce at most 40 nodes. "
        "Selection, coordinates, node areas, native scores, review flags and every plotted row are exported. "
        "Default selection is 12 interaction identities by best rank across all conditions; every condition uses that same selection and scales.\n\n"
        "Detection fractions use the minimum across all subunits of the complex in that population; native endpoint-subunit audit columns remain unmodified. "
        "If repeated rows unexpectedly supply different fractions for one node/context, the minimum is plotted and the observed maximum is also exported. "
        "Absent fractions are hollow nodes, not zero expression. Edge width uses native lr_gmeans, not a genotype log fold change.\n\n"
        "Ranks prioritize transcript-based candidates and are not genotype P values. Endpoint DESeq2 support is not a differential interaction test. "
        "CA1 Deep/Superficial are molecular programs, not measured physical layers. Ligand-to-receptor direction does not assign activation or inhibition.\n\n"
        "Reproduce with: python plot_gene_circos.py --evidence /absolute/interaction_evidence.tsv.gz --outdir /absolute/new_supplement_dir\n",encoding="utf-8")
    print(f"Supplemental report: {a.outdir/'index.html'}",flush=True)


if __name__=="__main__":
    main()
