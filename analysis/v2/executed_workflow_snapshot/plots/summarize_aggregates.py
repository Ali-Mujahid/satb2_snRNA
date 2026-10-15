#!/usr/bin/env python3
"""Summarize native full-cohort aggregates. No nucleus counts or model fitting."""
from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd

from plot_common import (BRANCHES, KEYS, bool_values, canonical_keys, empty_outdir,
                         numeric_checks, read_table, sha256, write_receipt, write_table)

SAMPLE_KEY = KEYS + ["sample"]
GENE_KEY = KEYS + ["gene_id", "gene_symbol"]
SUMS = ["n_nuclei", "sum_counts", "n_detected", "sum_log1p_normalized"]


def summarize(frame, counts, mouse_n_min, mouse_detection_min):
    """Weighted native sums, with zero-sample rows distinct from expression NA."""
    frame = canonical_keys(frame)
    counts = canonical_keys(counts)
    required = SAMPLE_KEY + ["condition", "gene_id", "gene_symbol"] + SUMS
    if set(required).difference(frame.columns):
        raise ValueError(f"Missing native aggregate columns: {set(required).difference(frame.columns)}")
    if frame.duplicated(SAMPLE_KEY + ["gene_id"]).any() or counts.duplicated(SAMPLE_KEY).any():
        raise ValueError("Duplicate native gene/sample/cluster or count-denominator row")
    numeric_checks(frame, ["n_nuclei", "sum_counts", "n_detected"], integer=True)
    numeric_checks(frame, ["sum_log1p_normalized"])
    numeric_checks(counts, ["n_nuclei"], integer=True)
    if (frame.n_detected > frame.n_nuclei).any():
        raise ValueError("Detected-nucleus count exceeds its denominator")
    if ((frame.n_nuclei == 0) & frame[SUMS[1:]].ne(0).any(axis=1)).any():
        raise ValueError("Zero-nucleus mouse/cluster has nonzero expression sums")
    joined = frame.merge(counts[SAMPLE_KEY + ["condition", "n_nuclei"]],
                         on=SAMPLE_KEY, suffixes=("", "_denominator"),
                         how="left", validate="many_to_one")
    if joined.n_nuclei_denominator.isna().any() or not joined.n_nuclei.eq(joined.n_nuclei_denominator).all():
        raise ValueError("Native gene denominator conflicts with cluster_mouse_counts")
    if not joined.condition.eq(joined.condition_denominator).all():
        raise ValueError("Native gene/denominator genotype metadata conflicts")
    all_mice = sorted(counts["sample"].unique())
    if len(all_mice) != 12 or set(counts.condition) != {"WT", "cKO"}:
        raise ValueError("Native summaries require the exact twelve source mice and WT/cKO")
    assignments = counts[["sample", "condition"]].drop_duplicates()
    if assignments["sample"].duplicated().any() or not assignments.condition.value_counts().eq(6).all():
        raise ValueError("Native summary mouse/genotype assignments conflict or are not six per genotype")
    if not counts.groupby(KEYS, observed=True)["sample"].nunique().eq(12).all():
        raise ValueError("An atomic cluster lacks explicit all-twelve-mouse coverage")
    n_genes = frame.gene_id.nunique()
    if len(frame) != len(counts) * n_genes:
        raise ValueError("Missing explicit zero gene/sample/cluster aggregate rows")
    if not frame.groupby(SAMPLE_KEY, observed=True).gene_id.nunique().eq(n_genes).all():
        raise ValueError("Unequal gene universes among native mouse/cluster rows")
    positive = frame.n_nuclei.gt(0)
    frame["mean_log1p_normalized"] = frame.sum_log1p_normalized.div(frame.n_nuclei.where(positive))
    frame["fraction_detected"] = frame.n_detected.div(frame.n_nuclei.where(positive))
    frame["mean_source_counts"] = frame.sum_counts.div(frame.n_nuclei.where(positive))
    frame["mouse_eligible_for_expression_review"] = frame.n_nuclei.ge(mouse_n_min)
    frame["mouse_expression_supported"] = (frame.mouse_eligible_for_expression_review
                                              & frame.fraction_detected.ge(mouse_detection_min))
    frame["measurement_status"] = np.where(positive, "OBSERVED_NATIVE_FULL_COHORT", "NO_NUCLEI")
    aggregates = []
    for condition in [None, "WT", "cKO"]:
        selected = frame if condition is None else frame.loc[frame.condition.eq(condition)]
        pooled = selected.groupby(GENE_KEY, observed=True, sort=False)[SUMS].sum().reset_index()
        support = selected.groupby(GENE_KEY, observed=True, sort=False).agg(
            n_mice_with_nuclei=("n_nuclei", lambda x: int(x.gt(0).sum())),
            n_mice_eligible_for_expression_review=("mouse_eligible_for_expression_review", "sum"),
            n_mice_expression_supported=("mouse_expression_supported", "sum")).reset_index()
        pooled = pooled.merge(support, on=GENE_KEY, validate="one_to_one")
        positive = pooled.n_nuclei.gt(0)
        pooled["mean_log1p_normalized"] = pooled.sum_log1p_normalized.div(pooled.n_nuclei.where(positive))
        pooled["fraction_detected"] = pooled.n_detected.div(pooled.n_nuclei.where(positive))
        pooled["mean_source_counts"] = pooled.sum_counts.div(pooled.n_nuclei.where(positive))
        pooled["condition_scope"] = condition or "ALL_GENOTYPES"
        pooled["measurement_status"] = np.where(positive, "OBSERVED_NATIVE_FULL_COHORT", "NO_NUCLEI_IN_CONDITION")
        aggregates.append(pooled)
    return frame, aggregates


def add_mask_and_contrasts(frame, annotations):
    key_mask = annotations[KEYS + ["keep_for_plot"]].copy()
    if key_mask.duplicated(KEYS).any():
        raise ValueError("Initial annotation has duplicate cluster definitions")
    key_mask["keep_for_plot"] = bool_values(key_mask.keep_for_plot)
    z = frame.merge(key_mask, on=KEYS, validate="many_to_one", how="left")
    if z.keep_for_plot.isna().any():
        raise ValueError("A native cluster lacks an initial annotation/display-mask definition")
    totals = z.loc[z.keep_for_plot].groupby(KEYS[:-1] + ["gene_id"], observed=True)[SUMS].sum().reset_index()
    totals = totals.rename(columns={c: c + "_retained_population" for c in SUMS})
    z = z.merge(totals, on=KEYS[:-1] + ["gene_id"], how="left", validate="many_to_one")
    # Own contribution is subtracted only when it belongs to the displayed population.
    own = z.keep_for_plot.astype(int)
    rest_n = z.n_nuclei_retained_population - z.n_nuclei * own
    rest_detected = z.n_detected_retained_population - z.n_detected * own
    rest_log = z.sum_log1p_normalized_retained_population - z.sum_log1p_normalized * own
    z["fraction_detected_rest"] = rest_detected.div(rest_n.where(rest_n.gt(0)))
    z["mean_log1p_normalized_rest"] = rest_log.div(rest_n.where(rest_n.gt(0)))
    z["pct_difference_vs_rest"] = z.fraction_detected - z.fraction_detected_rest
    z["mean_log1p_difference_vs_rest"] = z.mean_log1p_normalized - z.mean_log1p_normalized_rest
    return z


def add_reference_scaling(frame):
    """Per-gene z; paired corrected uses the identical paired-original reference."""
    frame = frame.copy()
    reference = frame.loc[frame.keep_for_plot & frame.mean_log1p_normalized.notna()]
    scales = reference.groupby(KEYS[:-1] + ["gene_id"], observed=True).mean_log1p_normalized.agg(
        scale_reference_mean="mean", scale_reference_sd=lambda x: float(x.std(ddof=0))).reset_index()
    records = []
    for branch in frame.branch.unique():
        target = frame.loc[frame.branch.eq(branch)].copy()
        target["scale_reference_branch"] = np.where(
            target.cluster_scope.eq("original_frozen") & target.branch.eq("paired_cellbender"),
            "paired_original", target.branch)
        sources = scales.rename(columns={"branch": "scale_reference_branch"})
        target = target.merge(sources, on=["tool", "scale_reference_branch", "population", "cluster_scope", "gene_id"],
                              validate="many_to_one", how="left")
        if target.loc[target.mean_log1p_normalized.notna(), "scale_reference_mean"].isna().any():
            raise ValueError("Scaled expression lacks its declared native reference branch")
        sd = target.scale_reference_sd
        target["scaled_expression_z"] = (target.mean_log1p_normalized - target.scale_reference_mean).div(sd.where(sd.gt(0)))
        target.loc[sd.eq(0) & target.mean_log1p_normalized.notna(), "scaled_expression_z"] = 0.0
        target["scale_definition"] = "gene_z_of_full_cohort_mean_log1p; equal_cluster_weight; population_sd_ddof0"
        records.append(target)
    return pd.concat(records, ignore_index=True)


def choose_representatives(pooled, cfg):
    p = cfg.get("plots", {})
    selected = pooled.loc[pooled.branch.eq("original_full") & pooled.cluster_scope.eq("original_frozen") & pooled.keep_for_plot].copy()
    # Perturbation target and housekeeping/nuisance genes are not identity representatives.
    nuisance = selected.gene_symbol.astype(str).str.match(r"^(mt-|Rpl|Rps)") | selected.gene_symbol.eq("Satb2")
    supported = (selected.fraction_detected.ge(float(p.get("representative_min_detection", .1)))
                 & selected.pct_difference_vs_rest.ge(float(p.get("representative_min_pct_difference", .15)))
                 & selected.mean_log1p_difference_vs_rest.ge(float(p.get("representative_min_mean_difference", .25)))
                 & selected.n_mice_expression_supported.ge(int(p.get("representative_min_supported_mice", 3)))
                 & selected.n_mice_expression_supported.ge(selected.n_mice_eligible_for_expression_review * .5)
                 & ~nuisance)
    selected["representative_criteria_pass"] = supported
    selected = selected.sort_values(["tool", "population", "cluster_id", "representative_criteria_pass",
                                     "pct_difference_vs_rest", "mean_log1p_difference_vs_rest", "n_mice_expression_supported", "gene_id"],
                                    ascending=[True, True, True, False, False, False, False, True])
    candidates = selected.copy()
    selected = selected.loc[selected.representative_criteria_pass].groupby(KEYS, observed=True, sort=False).head(
        int(p.get("representative_genes_per_cluster", 7))).copy()
    selected["representative_rank"] = selected.groupby(KEYS, observed=True).cumcount() + 1
    selected["selection_reference_branch"] = "original_full"
    selected["selection_scope"] = "descriptive_native_specificity_and_mouse_support; fixed_genes_for_paired_comparisons"
    registry = pooled.loc[pooled.branch.eq("original_full") & pooled.cluster_scope.eq("original_frozen")].groupby(KEYS, observed=True).first().reset_index()
    n = selected.groupby(KEYS, observed=True).size().rename("n_supported_representative_genes").reset_index()
    registry = registry[KEYS + ["n_nuclei", "keep_for_plot"]].merge(n, on=KEYS, how="left", validate="one_to_one")
    registry["n_supported_representative_genes"] = registry.n_supported_representative_genes.fillna(0).astype(int)
    registry["marker_review_status"] = np.where(registry.n_supported_representative_genes.gt(0), "SUPPORTED_DESCRIPTIVE_CANDIDATES", "NO_GENE_MEETS_DECLARED_SUPPORT_CRITERIA")
    registry.loc[~registry.keep_for_plot, "marker_review_status"] = "EXCLUDED_FROM_PLOT_BY_EXPLICIT_ATOMIC_CLUSTER_MASK"
    return selected, candidates, registry


def program_evidence(wt, panel):
    if "group" not in panel:
        panel["group"] = panel.get("marker_group", panel.get("set", "Provided_marker_panel"))
    mapped = panel.loc[panel.gene_id.notna(), ["group", "gene_id", "gene_symbol"]].drop_duplicates()
    joined = wt.merge(mapped.rename(columns={"gene_symbol": "panel_gene_symbol"}), on="gene_id", how="inner", validate="many_to_many")
    outputs = []
    for keys, x in joined.groupby(KEYS + ["group"], observed=True, sort=False):
        row = dict(zip(KEYS + ["group"], keys))
        row.update(n_WT_nuclei=int(x.n_nuclei.iloc[0]), n_mapped_program_genes=int(x.gene_id.nunique()),
                   mean_WT_gene_z=float(x.scaled_expression_z.mean()),
                   mean_WT_fraction_detected=float(x.fraction_detected.mean()),
                   n_positive_WT_gene_z=int(x.scaled_expression_z.gt(0).sum()),
                   n_genes_supported_in_at_least_three_WT_mice=int(x.n_mice_expression_supported.ge(3).sum()),
                   member_gene_ids=";".join(x.gene_id.astype(str)), member_gene_symbols=";".join(x.panel_gene_symbol.astype(str)),
                   evidence_scope="WT_only_full_cohort_native_expression_program_arithmetic; not_anatomical_layer_or_cell_level_scoring",
                   physical_layer_assignment="UNDETERMINED")
        if row["n_WT_nuclei"] == 0:
            row["mean_WT_gene_z"] = row["mean_WT_fraction_detected"] = np.nan
        outputs.append(row)
    return pd.DataFrame(outputs)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", type=Path, required=True)
    ap.add_argument("--aggregates", type=Path, action="append")
    ap.add_argument("--cluster-counts", type=Path, action="append")
    ap.add_argument("--initial-annotations", type=Path)
    ap.add_argument("--panel", type=Path)
    ap.add_argument("--panel-memberships", type=Path)
    ap.add_argument("--outdir", type=Path, required=True)
    args = ap.parse_args()
    cfg = json.loads(args.config.read_text())
    prepared = Path(cfg.get("prepared_root", args.outdir.parent / "prepared"))
    branches = cfg.get("branches", BRANCHES)
    aggregates = args.aggregates or [prepared / b / "gene_cluster_mouse_aggregates.tsv.gz" for b in branches]
    counts_paths = args.cluster_counts or [prepared / b / "cluster_mouse_counts.tsv" for b in branches]
    if len(aggregates) != len(counts_paths):
        raise ValueError("Provide one cluster-counts table for every aggregate file")
    annotations_path = args.initial_annotations or Path(cfg.get("initial_annotation_dir", cfg.get("annotation_initial_dir", args.outdir.parent / "annotation_initial"))) / "cluster_annotations.tsv"
    panel_path = args.panel or Path(cfg.get("plot_gene_panel_tsv", annotations_path.parent / "gene_panel.tsv"))
    memberships_path = args.panel_memberships or Path(cfg.get("marker_panel_memberships_tsv", annotations_path.parent / "panel_membership.tsv"))
    annotations = canonical_keys(read_table(annotations_path))
    if "keep_for_plot" not in annotations:
        annotations["keep_for_plot"] = annotations.keep_for_analysis
    panel = read_table(panel_path)
    memberships = read_table(memberships_path)
    out = empty_outdir(args.outdir)
    p = cfg.get("plots", {})
    pooled, wt_rows, condition_rows, count_tables = [], [], [], []
    native_receipts = []
    with gzip.open(out / "gene_cluster_mouse_summary.tsv.gz", "wt") as handle:
        header = True
        for aggregate_path, counts_path in zip(aggregates, counts_paths):
            if aggregate_path.parent.resolve() != counts_path.parent.resolve():
                raise ValueError("Native aggregate and denominator tables must share their prepare stage")
            native_receipt_path = aggregate_path.parent / "COMPLETE.json"
            native = json.loads(native_receipt_path.read_text())
            if native.get("status") != "EXECUTED_REVIEW_REQUIRED" or native.get("stage") != "prepare":
                raise ValueError("Native full-cohort aggregate prepare producer is not complete")
            if native.get("scientific_fitting_performed") is not False or native.get("raw_nucleus_count_aggregation_performed") is not True:
                raise ValueError("Native prepare receipt lacks its authorized count-aggregation scope")
            for path in [aggregate_path, counts_path]:
                if native.get("output_sha256", {}).get(path.name) != sha256(path):
                    raise ValueError(f"Native prepared aggregate hash mismatch: {path.name}")
            frame = canonical_keys(read_table(aggregate_path))
            if set(frame.branch) != {native.get("branch")}:
                raise ValueError("Native aggregate branch disagrees with its prepare receipt")
            if "n_panel_genes" in native and frame.gene_id.nunique() != int(native["n_panel_genes"]):
                raise ValueError("Native aggregate gene-panel size disagrees with its receipt")
            native_receipts.append(native_receipt_path)
            counts = canonical_keys(read_table(counts_path))
            mouse, summaries = summarize(frame, counts,
                                         int(p.get("marker_review_mouse_n_min", 5)),
                                         float(p.get("marker_review_mouse_detection_min", .1)))
            mouse.to_csv(handle, sep="\t", index=False, na_rep="NA", header=header)
            header = False
            pooled.append(summaries[0]); wt_rows.append(summaries[1]); condition_rows += summaries[1:]
            count_tables.append(counts)
    pooled = add_reference_scaling(add_mask_and_contrasts(pd.concat(pooled, ignore_index=True), annotations))
    wt = add_reference_scaling(add_mask_and_contrasts(pd.concat(wt_rows, ignore_index=True), annotations))
    wt["n_WT_nuclei"] = wt.n_nuclei
    wt["n_WT_mice"] = wt.n_mice_with_nuclei
    wt["WT_program_annotation_reference_branch"] = "original_full"
    conditions = pd.concat(condition_rows, ignore_index=True)
    if pooled.duplicated(GENE_KEY).any() or wt.duplicated(GENE_KEY).any():
        raise ValueError("Input aggregate files repeat branch/atomic-cluster definitions")
    representatives, candidates, marker_registry = choose_representatives(pooled, cfg)
    write_table(pooled, out / "gene_cluster_summary.tsv")
    write_table(wt, out / "gene_cluster_summary_WT.tsv")
    write_table(conditions, out / "gene_cluster_condition_summary.tsv")
    write_table(pd.concat(count_tables, ignore_index=True), out / "cluster_mouse_counts.tsv")
    write_table(representatives, out / "representative_markers.tsv")
    write_table(candidates, out / "representative_candidate_audit.tsv.gz")
    write_table(marker_registry, out / "representative_cluster_registry.tsv")
    write_table(program_evidence(wt, memberships), out / "program_evidence_WT.tsv")
    write_table(memberships.loc[memberships.gene_id.isna()], out / "unmapped_panel_entries.tsv")
    write_receipt(out, args.config, aggregates + counts_paths + native_receipts + [annotations_path, panel_path, memberships_path],
                  dict(stage="native_aggregate_summary", n_gene_cluster_rows=len(pooled),
                       mean_expression_definition="sum_native_log1p(count/full_source_library*10000)/n_nuclei",
                       percent_definition="native_n_detected/n_nuclei; no_downsampling",
                       representative_settings=p, WT_program_evidence_only=True,
                       negative_z_visual_policy="handled_by_renderer; exact_z_and_percent_preserved"))


if __name__ == "__main__":
    main()
