#!/usr/bin/env python3
"""Freeze reviewed original/paired inputs and workbook-derived review panels."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024**2), b""):
            h.update(chunk)
    return h.hexdigest()


def truth(values):
    values = values.astype(str).str.lower()
    if not values.isin(["true", "false", "1", "0"]).all():
        raise ValueError("Boolean fields must be explicit")
    return values.isin(["true", "1"])


def save_table(frame, path):
    if path.exists():
        raise RuntimeError(f"Refusing to overwrite {path}")
    frame.to_csv(path, sep="\t", index=False, na_rep="NA")


def workbook_table(cells, workbook, sheet):
    frame = cells.loc[(cells.workbook == workbook) & (cells.sheet == sheet)].copy()
    frame["column"] = frame.cell.str.extract(r"^([A-Z]+)")[0]
    frame["row"] = frame.cell.str.extract(r"([0-9]+)$")[0].astype(int)
    table = frame.pivot(index="row", columns="column", values="value")
    return table


def review_panels(cells_path):
    cells = pd.read_csv(cells_path, sep="\t", dtype=str, keep_default_na=False)
    entries = []
    canonical = {
        "Neuron": ["Snap25", "Rbfox3", "Syt1"],
        "Excitatory": ["Slc17a7", "Slc17a6", "Camk2a"],
        "Inhibitory": ["Gad1", "Gad2", "Slc32a1"],
        "Oligodendrocyte": ["Mbp", "Plp1", "Mog"],
        "Astrocyte": ["Aqp4", "Slc1a3", "Aldoc"],
        "OPC": ["Pdgfra", "Cspg4"],
        "Microglia": ["C1qa", "C1qb", "P2ry12"],
        "Endothelial": ["Cldn5", "Pecam1"],
    }
    for group, genes in canonical.items():
        entries.extend({"group": group, "gene_symbol": gene,
                        "source": "canonical coarse-lineage expression review",
                        "selection": "diagnostic expression only; does not define cohort"}
                       for gene in genes)
    ca1 = workbook_table(cells, "CA1.markers.xlsx", "PyrCA1_markers").drop(index=1)
    for col in ["B", "C", "D", "F"]:
        ca1[col] = pd.to_numeric(ca1[col], errors="coerce")
    ca1["specificity"] = ca1.C - ca1.D
    ca1 = ca1.loc[(ca1.F <= .05) & (ca1.B > 0) & (ca1.specificity > .1)]
    ca1 = ca1.sort_values(["specificity", "B"], ascending=False).head(20)
    for row, item in ca1.iterrows():
        entries.append({"group": "CA1_population_candidates", "gene_symbol": item.A,
                        "source": f"CA1.markers.xlsx:PyrCA1_markers!A{row}:F{row}",
                        "selection": "positive FC,FDR<=.05,pct difference>.1; top20 specificity"})
    layer = workbook_table(cells, "CA1.markers.xlsx", "DEGs_PyrCA1Sup_vs_PyrCA1Deep").drop(index=1)
    for col in ["B", "C", "D", "F"]:
        layer[col] = pd.to_numeric(layer[col], errors="coerce")
    layer = layer.loc[(layer.F <= .05) & (abs(layer.C-layer.D) > .1) & (layer.A != "Satb2")]
    layer = layer.loc[~layer.A.str.match(r"^(mt-|Rpl|Rps)")]
    for positive, group in [(True, "Superficial_associated_VEH_program"),
                            (False, "Deep_associated_VEH_program")]:
        rows = layer.loc[layer.B > 0] if positive else layer.loc[layer.B < 0]
        rows = rows.assign(abs_fc=abs(rows.B)).sort_values("abs_fc", ascending=False).head(15)
        for row, item in rows.iterrows():
            entries.append({"group": group, "gene_symbol": item.A,
                            "source": f"CA1.markers.xlsx:DEGs_PyrCA1Sup_vs_PyrCA1Deep!A{row}:F{row}",
                            "selection": "VEH-only program review; top15 absoluteFC; Satb2/mt/ribo excluded; no anatomical calls"})
    interneuron = workbook_table(cells, "Interneuron classification marker.xlsx", "Sheet1")
    for col in ["C", "D", "E", "F"]:
        interneuron[col] = pd.to_numeric(interneuron[col], errors="coerce")
    interneuron = interneuron.loc[(interneuron.F <= .05) & (interneuron.C > 0)]
    interneuron["specificity"] = interneuron.D-interneuron.E
    for name, rows in interneuron.groupby("G", dropna=False):
        rows = rows.loc[rows.specificity > .1].sort_values(["specificity", "C"], ascending=False).head(5)
        for row, item in rows.iterrows():
            entries.append({"group": str(name), "gene_symbol": item.A,
                            "source": f"Interneuron classification marker.xlsx:Sheet1!A{row}:G{row}",
                            "selection": "top5 specificity,positiveFC,FDR<=.05; Module labels are programs, not cell types"})
    result = pd.DataFrame(entries).drop_duplicates(["group", "gene_symbol"])
    result["set"] = result.group
    result["eligible"] = True
    result["annotation_scope"] = "expression review; no forced cell/layer identity"
    return result


def qc_report(full, primary, pairs, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    samples = sorted(full["sample"].unique())
    summary = []
    for sample in samples:
        source = full.loc[full["sample"] == sample]
        retained = primary.loc[primary["sample"] == sample]
        paired = pairs.loc[pairs["sample"] == sample]
        summary.append({"sample": sample, "condition": source.condition.iloc[0],
                        "called_nuclei": len(source), "primary_QC_nuclei": len(retained),
                        "recorded_doublet_calls": int((source.call_new == "doublet").sum()),
                        "paired_count_nuclei": len(paired),
                        "primary_original_UMIs": int(retained.nCount_RNA.sum()),
                        "paired_original_UMIs": int(paired.nCount_RNA.sum()),
                        "paired_corrected_UMIs": int(paired.cellbender_nCount_RNA.sum())})
    summary = pd.DataFrame(summary)
    save_table(summary, out / "QC_RETENTION_AND_COUNT_SUMMARY.tsv")
    with PdfPages(out / "QC_REVIEW.pdf") as pdf:
        fig, ax = plt.subplots(figsize=(11, 5), constrained_layout=True)
        pos = np.arange(len(samples))
        ax.bar(pos-.25, summary.called_nuclei, .25, label="Cell Ranger called")
        ax.bar(pos, summary.primary_QC_nuclei, .25, label="Original audited QC/singlets")
        ax.bar(pos+.25, summary.paired_count_nuclei, .25, label="Matched corrected-count support")
        ax.set_xticks(pos, samples, rotation=45, ha="right")
        ax.set(ylabel="Nuclei", title="All twelve mice: frozen QC policy and paired-count coverage")
        ax.legend(frameon=False)
        pdf.savefig(fig); plt.close(fig)
        for column, title, logarithm in [("nCount_RNA", "Original UMIs per nucleus", True),
                                         ("nFeature_RNA", "Original detected genes per nucleus", False),
                                         ("pct_mt", "Original mitochondrial UMI fraction (%)", False)]:
            fig, ax = plt.subplots(figsize=(12, 5), constrained_layout=True)
            before = [full.loc[full["sample"] == sample, column].to_numpy() for sample in samples]
            after = [primary.loc[primary["sample"] == sample, column].to_numpy() for sample in samples]
            if logarithm:
                before = [np.log10(v+1) for v in before]
                after = [np.log10(v+1) for v in after]
            for values, shift, color in [(before, -.18, "#8ba6ba"), (after, .18, "#305f82")]:
                box = ax.boxplot(values, positions=pos+shift, widths=.3, showfliers=False,
                                 patch_artist=True, manage_ticks=False)
                for item in box["boxes"]:
                    item.set_facecolor(color)
            ax.set_xticks(pos, samples, rotation=45, ha="right")
            ax.set(ylabel="log10(UMIs+1)" if logarithm else title,
                   title=title + "; light=called, dark=retained")
            pdf.savefig(fig); plt.close(fig)
        fig, axes = plt.subplots(3, 4, figsize=(14, 10), constrained_layout=True)
        retained_ids = set(primary.nucleus_id)
        for ax, sample in zip(axes.ravel(), samples):
            source = full.loc[full["sample"] == sample]
            keep = source.nucleus_id.isin(retained_ids)
            ax.scatter(np.log10(source.nCount_RNA+1), source.nFeature_RNA,
                       c=np.where(keep, "#305f82", "#c98982"), s=1, rasterized=True, alpha=.5)
            ax.set(title=sample, xlabel="log10 original UMIs+1", ylabel="Detected genes")
        fig.suptitle("Original nuclear complexity; red=outside frozen primary QC")
        pdf.savefig(fig); plt.close(fig)
        fig, ax = plt.subplots(figsize=(11, 5), constrained_layout=True)
        ax.boxplot([pairs.loc[pairs["sample"] == sample, "umi_fraction_removed"].to_numpy()
                    for sample in samples], showfliers=False, tick_labels=samples)
        ax.tick_params(axis="x", rotation=45)
        ax.set(ylabel="Fraction of original UMIs removed", title="CellBender: identical paired nuclei per mouse")
        pdf.savefig(fig); plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    args = parser.parse_args()
    cfg = json.loads(args.config.read_text())
    accepted_cb = json.loads(Path(cfg["cellbender_receipt"]).read_text())
    for key, expected in [("original_source_hashes_matched", 24), ("primary_QC_nuclei", 108062),
                          ("paired_supported_nuclei", 104298),
                          ("all_filtered_corrected_nonnegative_integers", True),
                          ("all_features_match_original_order", True)]:
        if accepted_cb.get(key) != expected:
            raise ValueError(f"Reviewed CellBender audit gate changed: {key}")
    out = Path(cfg["run_root"]) / "results" / "inputs"
    if out.exists() and (not out.is_dir() or any(out.iterdir())):
        raise RuntimeError("Input snapshot directory contains data; use a fresh run")
    # Snakemake creates parent folders for declared output files before execution.
    out.mkdir(parents=True, exist_ok=True)
    original_sheet = pd.read_csv(cfg["original_manifest"], sep="\t")
    for row in original_sheet.itertuples(index=False):
        if digest(row.source_h5) != row.source_h5_sha256:
            raise ValueError(f"Previously accepted original source hash changed: {row.sample}")
    original = pd.read_csv(cfg["original_roster"], sep="\t", low_memory=False)
    full_qc = pd.read_csv(cfg["full_called_qc"], sep="\t", low_memory=False)
    predictions = pd.read_csv(cfg["mapmycells_predictions"], sep="\t", low_memory=False)
    genes = pd.read_csv(cfg["gene_manifest"], sep="\t")
    cb = pd.read_csv(cfg["cellbender_manifest"], sep="\t")
    if len(original) != 108062 or original.nucleus_id.duplicated().any():
        raise ValueError("Frozen original roster contract failed")
    if len(full_qc) != 121486 or full_qc.nucleus_id.duplicated().any():
        raise ValueError("Audited called-nucleus QC roster contract failed")
    if set(full_qc.loc[truth(full_qc.passes_primary_QC), "nucleus_id"]) != set(original.nucleus_id):
        raise ValueError("Full audited QC and retained roster do not reconcile")
    if set(predictions.nucleus_id) != set(original.nucleus_id):
        raise ValueError("MapMyCells source query roster mismatch")
    if not truth(original.passes_primary_QC).all():
        raise ValueError("Original roster has unretained profiles")
    availability = pd.concat([pd.read_csv(path, sep="\t", low_memory=False)
                              for path in cb.primary_paired_roster_path], ignore_index=True)
    if availability.nucleus_id.duplicated().any() or set(availability.nucleus_id) != set(original.nucleus_id):
        raise ValueError("CellBender availability roster mismatch")
    availability = availability.set_index("nucleus_id").loc[original.nucleus_id]
    paired_ids = availability.loc[truth(availability.paired_sensitivity_usable)].index
    if len(paired_ids) != 104298:
        raise ValueError("Audited paired count cohort changed")
    pairs = original.loc[original.nucleus_id.isin(paired_ids)].copy()
    for column in ["cellbender_analyzed", "cellbender_nonzero_full", "cellbender_filtered_present",
                   "cellbender_nCount_RNA", "cellbender_nFeature_RNA", "paired_sensitivity_usable",
                   "umi_fraction_removed"]:
        original[column] = availability[column].to_numpy()
        pairs[column] = availability.loc[pairs.nucleus_id, column].to_numpy()
    save_table(original, out / "roster_original_full.tsv")
    save_table(pairs, out / "roster_paired_original.tsv")
    save_table(pairs, out / "roster_paired_cellbender.tsv")
    # Reference fields are provenance-preserved external annotation, not inherited clusters.
    fields = ["nucleus_id"] + [c for c in original if c.startswith(("reference_", "REF_"))]
    save_table(original[fields], out / "labels_original_full.tsv")
    save_table(pairs[fields], out / "labels_paired_original.tsv")
    save_table(pairs[fields], out / "labels_paired_cellbender.tsv")
    save_table(predictions, out / "MapMyCells_original_predictions.tsv")
    save_table(genes, out / "gene_features.tsv")
    coverage = []
    for branch, roster in [("original_full", original), ("paired_original", pairs), ("paired_cellbender", pairs)]:
        sheet = original_sheet.copy()
        sheet["count_branch"] = branch
        totals = roster.groupby("sample").nCount_RNA.sum()
        sheet["original_UMIs"] = sheet["sample"].map(totals).astype(np.int64)
        sheet["n_nuclei"] = sheet["sample"].map(roster.groupby("sample").size()).astype(int)
        if branch == "paired_cellbender":
            paths = cb.set_index("sample").corrected_full_path
            sheet["source_h5"] = sheet["sample"].map(paths)
            sheet["branch_UMIs"] = sheet["sample"].map(roster.groupby("sample").cellbender_nCount_RNA.sum()).astype(np.int64)
        else:
            sheet["branch_UMIs"] = sheet.original_UMIs
        sheet["filtered_h5"] = sheet.source_h5
        sheet["source_h5_sha256"] = [digest(path) for path in sheet.source_h5]
        save_table(sheet, out / f"samples_{branch}.tsv")
        for row in sheet.itertuples(index=False):
            coverage.append({"branch": branch, "sample": row.sample, "condition": row.condition,
                             "n_nuclei": row.n_nuclei, "original_UMIs": row.original_UMIs,
                             "branch_UMIs": row.branch_UMIs, "source_h5_sha256": row.source_h5_sha256})
    save_table(pd.DataFrame(coverage), out / "COUNT_BRANCH_COVERAGE.tsv")
    save_table(review_panels(cfg["marker_cells"]), out / "marker_review_panels.tsv")
    qc_report(full_qc, original, pairs, out)
    run = Path(cfg["run_root"])
    source_counts = {"original_full": 108062, "paired_original": 104298, "paired_cellbender": 104298}
    scanpy_paths = {population: str(run / "results" / "scanpy" / "original_full" / population / "nucleus_metadata.tsv")
                    for population in ["allcells", "Excitatory", "Inhibitory", "CA1_ProS"]}
    for branch, expected in source_counts.items():
        shared = run / "results" / "shared" / branch
        scanpy = dict(cfg["scanpy_parameters"])
        scanpy.update(sample_sheet=str(out / f"samples_{branch}.tsv"),
                      qc_roster=str(out / f"roster_{branch}.tsv"),
                      mapmycells_predictions=str(out / "MapMyCells_original_predictions.tsv"),
                      markers_tsv=str(out / "marker_review_panels.tsv"),
                      count_branch=branch, expected_nuclei=expected,
                      mapmycells_allow_superset=branch != "original_full",
                      annotation_provenance="existing original-query MapMyCells labels reused; paired denoising uses fixed original labels")
        accepted_shared = cfg.get("reuse_prepared_shared", {}).get(branch)
        if accepted_shared:
            scanpy["reuse_prepared_shared_dir"] = accepted_shared
        accepted_populations = cfg.get("reuse_completed_scanpy", {}).get(branch)
        if accepted_populations:
            scanpy["reuse_completed_populations"] = accepted_populations
        if branch != "original_full":
            scanpy["hvg_source_dir"] = str(run / "results" / "shared" / "original_full")
        (out / f"scanpy_{branch}.json").write_text(json.dumps(scanpy, indent=2) + "\n")
        seurat = dict(cfg["seurat_parameters"])
        seurat.update(sample_manifest_tsv=str(out / f"samples_{branch}.tsv"),
                      roster_tsv=str(out / f"roster_{branch}.tsv"),
                      gene_manifest_tsv=str(out / "gene_features.tsv"),
                      fixed_features_tsv=str(shared / "common_hvg_allcells.tsv"),
                      labels_tsv=str(out / f"labels_{branch}.tsv"),
                      marker_tsv=str(out / "marker_review_panels.tsv"),
                      output_root=str(run / "results" / "seurat" / branch),
                      count_source=branch, expected_nuclei=expected,
                      subset_features_tsv={"EX": str(shared / "common_hvg_Excitatory.tsv"),
                                           "IN": str(shared / "common_hvg_Inhibitory.tsv"),
                                           "CA1_PROS_CANDIDATE": str(shared / "common_hvg_CA1_ProS.tsv")},
                      scanpy_metadata_paths=scanpy_paths)
        if branch != "original_full":
            seurat.update(cluster_metadata_root=str(run / "results" / "seurat" / "original_full"),
                          prefilter_universe_root=str(run / "results" / "seurat" / "original_full" / "pseudobulk"),
                          run_loo=False, run_partition_sensitivities=False)
        (out / f"seurat_{branch}.json").write_text(json.dumps(seurat, indent=2) + "\n")
    source_keys = ["original_manifest", "original_roster", "mapmycells_predictions", "gene_manifest",
                   "cellbender_manifest", "cellbender_receipt", "marker_cells", "full_called_qc"]
    provenance = [{"role": key, "path": cfg[key], "sha256": digest(cfg[key])} for key in source_keys]
    save_table(pd.DataFrame(provenance), out / "SOURCE_PROVENANCE.tsv")
    receipt = {"status": "COMPLETE", "primary_nuclei": 108062, "paired_nuclei": 104298,
               "genotype_replication": "six WT and six cKO independent mice in every branch",
               "MapMyCells": "existing same-original-count query labels, verified exact nucleus IDs; no identical mapping refit",
               "QC": "audited original primary QC and singlet calls reused; no new exclusions",
               "files": {p.name: digest(p) for p in out.iterdir() if p.is_file()}}
    (out / "COMPLETE.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps({"status": "COMPLETE", "primary": 108062, "paired": 104298}))


if __name__ == "__main__":
    main()
