#!/usr/bin/env python3
"""Synthetic-only arithmetic/render check: no study counts, clustering, or fit."""
from pathlib import Path
import argparse
import json
import subprocess
import sys

import numpy as np
import pandas as pd

from plot_common import BRANCHES, KEYS, sha256
from render_plots import dot_visual_values


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--outdir", type=Path, required=True)
    args = p.parse_args()
    out = args.outdir.resolve()
    if out.exists():
        raise FileExistsError("Synthetic fixture directory already exists")
    out.mkdir(parents=True)
    initial = out / "annotation_initial"
    initial.mkdir()
    previous = out / "synthetic_previous_run"
    source = Path(__file__).resolve().parent
    genes = [("gA", "SyntheticA"), ("gB", "SyntheticB"), ("gC", "SyntheticC")]
    mice = [(f"wt_{i}", "WT") for i in range(6)] + [(f"cKO_{i}", "cKO") for i in range(6)]
    ca, na = [], []
    for branch in BRANCHES:
        aggregate, coverage, metadata = [], [], []
        for cluster in ["0", "1", "2"]:
            total = 0
            for i, (sample, condition) in enumerate(mice):
                n = (20 if i == 0 else 10) if cluster == "0" else (6 if cluster == "1" else (0 if condition == "WT" else 5))
                key = dict(tool="scanpy", branch=branch, population="Inhibitory", cluster_scope="original_frozen", cluster_id=cluster)
                coverage.append(key | dict(sample=sample, condition=condition, n_nuclei=n))
                total += n
                for gene, symbol in genes:
                    high = (gene == "gA" and cluster == "0") or (gene == "gB" and cluster == "1") or (gene == "gC" and cluster == "2")
                    detected = n if high else (n // 10)
                    counts = detected * (3 if high else 1)
                    mean = (2.0 if high else .1) if n else 0.0
                    if gene == "gA" and cluster == "0":
                        mean = 2.0 if i == 0 else 1.0
                    if branch == "paired_cellbender" and gene == "gB":
                        counts = detected
                        mean *= .5
                    aggregate.append(key | dict(sample=sample, condition=condition, gene_id=gene, gene_symbol=symbol,
                                                n_nuclei=n, sum_counts=counts, n_detected=detected,
                                                sum_log1p_normalized=mean * n))
                for j in range(n):
                    nucleus = f"{sample}_C{cluster}_{j}"
                    annotation = key | dict(nucleus_id=nucleus, sample=sample, condition=condition,
                                            reviewed_label=f"Synthetic expression group {cluster}", broad_label="Inhibitory",
                                            label_confidence="SYNTHETIC_ONLY", layer_program_label="UNDETERMINED",
                                            keep_for_analysis=True, keep_for_plot=True)
                    na.append(annotation)
                    metadata.append(dict(nucleus_id=nucleus, sample=sample, condition=condition,
                                         reference_broad_prediction="Inhibitory", reference_subclass_name=f"Synthetic subclass {cluster}",
                                         UMAP_harmony_1=float(cluster) * 4 + np.sin(j),
                                         UMAP_harmony_2=i * .2 + np.cos(j),
                                         UMAP_unintegrated_1=float(cluster) * 3 + np.sin(j),
                                         UMAP_unintegrated_2=i * .3 + np.cos(j)))
            ca.append(key | dict(n_nuclei=total, reviewed_label=f"Synthetic expression group {cluster}", broad_label="Inhibitory",
                                 label_confidence="SYNTHETIC_ONLY", keep_for_analysis=True, keep_for_plot=True))
        prepared = out / "prepared" / branch
        prepared.mkdir(parents=True)
        apath, cpath = prepared / "gene_cluster_mouse_aggregates.tsv.gz", prepared / "cluster_mouse_counts.tsv"
        pd.DataFrame(aggregate).to_csv(apath, sep="\t", index=False)
        pd.DataFrame(coverage).to_csv(cpath, sep="\t", index=False)
        (prepared / "COMPLETE.json").write_text(json.dumps({"status": "EXECUTED_REVIEW_REQUIRED", "stage": "prepare", "branch": branch,
                                                            "scientific_fitting_performed": False, "raw_nucleus_count_aggregation_performed": True,
                                                            "n_panel_genes": len(genes), "synthetic_only": True,
                                                            "output_sha256": {apath.name: sha256(apath), cpath.name: sha256(cpath)}}))
        previous_stage = previous / "results/scanpy" / branch / "Inhibitory"
        previous_stage.mkdir(parents=True)
        pd.DataFrame(metadata).to_csv(previous_stage / "nucleus_metadata.tsv", sep="\t", index=False)
    pd.DataFrame(ca).to_csv(initial / "cluster_annotations.tsv", sep="\t", index=False)
    pd.DataFrame(na).to_csv(initial / "nucleus_annotations.tsv.gz", sep="\t", index=False)
    pd.DataFrame(genes, columns=["gene_id", "gene_symbol"]).to_csv(initial / "gene_panel.tsv", sep="\t", index=False)
    members = pd.DataFrame([dict(group="Synthetic superficial-associated panel", set="Synthetic superficial-associated panel", gene_id="gA", gene_symbol="SyntheticA"),
                            dict(group="Synthetic deep-associated panel", set="Synthetic deep-associated panel", gene_id="gB", gene_symbol="SyntheticB"),
                            dict(group="Synthetic missing-WT panel", set="Synthetic missing-WT panel", gene_id="gC", gene_symbol="SyntheticC"),
                            dict(group="Missing exact symbol", set="Missing exact symbol", gene_id=None, gene_symbol="UNMAPPED_SYNTHETIC")])
    members.to_csv(initial / "panel_membership.tsv", sep="\t", index=False)
    (initial / "COMPLETE.json").write_text(json.dumps({"status": "ANNOTATION_FINAL_COMPLETE_REVIEW_REQUIRED", "stage": "final",
                                                        "synthetic_only": True,
                                                        "output_sha256": {name: sha256(initial / name) for name in ["cluster_annotations.tsv", "nucleus_annotations.tsv.gz"]}}))
    cfg = dict(previous_run_root=str(previous), initial_annotation_dir=str(initial), prepared_root=str(out / "prepared"),
               plot_gene_panel_tsv=str(initial / "gene_panel.tsv"), marker_panel_memberships_tsv=str(initial / "panel_membership.tsv"))
    config = out / "synthetic_config.json"
    config.write_text(json.dumps(cfg, indent=2))
    summaries = out / "marker_summaries"
    subprocess.run([sys.executable, str(source / "summarize_aggregates.py"), "--config", str(config), "--outdir", str(summaries)], check=True)
    pooled = pd.read_csv(summaries / "gene_cluster_summary.tsv", sep="\t", dtype={"cluster_id": str})
    wt = pd.read_csv(summaries / "gene_cluster_summary_WT.tsv", sep="\t", dtype={"cluster_id": str})
    row = pooled[(pooled.branch == "original_full") & (pooled.cluster_id == "0") & (pooled.gene_id == "gA")].iloc[0]
    assert np.isclose(row.mean_log1p_normalized, 150 / 130)
    assert wt.loc[wt.cluster_id.eq("2"), "n_WT_nuclei"].eq(0).all()
    assert wt.loc[wt.cluster_id.eq("2"), "n_WT_mice"].eq(0).all()
    assert wt.loc[wt.cluster_id.eq("2"), ["mean_log1p_normalized", "fraction_detected"]].isna().all().all()
    ref = pooled[pooled.branch.eq("paired_cellbender")]
    assert ref.scale_reference_branch.eq("paired_original").all()
    x = pd.DataFrame({"scaled_expression_z": [-1.0, 1.0], "fraction_detected": [.8, .8]})
    glyph = dot_visual_values(x)
    assert glyph.fraction_detected.tolist() == [.8, .8]
    assert glyph.scaled_expression_z.tolist() == [-1.0, 1.0]
    assert glyph.glyph_area_points2.tolist() == [20, 80]
    assert glyph.glyph_alpha.tolist() == [.18, .9]
    plots = out / "plots"
    subprocess.run([sys.executable, str(source / "render_plots.py"), "--config", str(config),
                    "--summaries", str(summaries), "--annotations", str(initial), "--outdir", str(plots)], check=True)
    result = json.loads((plots / "COMPLETE.json").read_text())
    assert result["n_plot_books"] == 3 and result["scientific_fitting_performed"] is False
    checks = dict(synthetic_only=True, study_nucleus_counts_read=False, scientific_fitting_performed=False,
                  weighted_denominator_mean_verified=True, zero_WT_preserved_as_NA=True,
                  corrected_reference_anchored_to_paired_original=True,
                  negative_z_exact_percent_retained_and_glyph_downweight_verified=True,
                  render_books=result["n_plot_books"], render_pages=result["n_pages"])
    (out / "SYNTHETIC_TEST_RECEIPT.json").write_text(json.dumps(checks, indent=2) + "\n")
    print(json.dumps(checks))


if __name__ == "__main__":
    main()
