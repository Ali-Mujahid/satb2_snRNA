#!/usr/bin/env python3
"""Check completed produced plot artifacts/metadata; no fits or RNA aggregation."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd


def sha(p):
    h = hashlib.sha256()
    with Path(p).open("rb") as f:
        for x in iter(lambda: f.read(8 * 1024 * 1024), b""): h.update(x)
    return h.hexdigest()


def truth(s):
    z = s.astype(str).str.lower().map({"true": True, "false": False, "1": True, "0": False})
    if z.isna().any(): raise ValueError("Invalid mask")
    return z.astype(bool)


def read(p, **kwargs): return pd.read_csv(p, sep="\t", dtype={"cluster_id": str, "nucleus_id": str}, **kwargs)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--config", type=Path, required=True); ap.add_argument("--outdir", type=Path, required=True)
    args = ap.parse_args(); cfg = json.loads(args.config.read_text()); root = Path(cfg["results_root"]); out = args.outdir
    if out.exists(): raise ValueError("Review destination already exists")
    receipt_path = root / "plots/COMPLETE.json"
    receipt = json.loads(receipt_path.read_text())
    if receipt.get("status") != "COMPLETE": raise ValueError("Global plot producer incomplete")
    checked = 0
    for relative, digest in receipt["output_sha256"].items():
        if sha(root / "plots" / relative) != digest: raise ValueError(f"Plot artifact hash differs: {relative}")
        checked += 1
    ann = read(Path(cfg["final_annotation_dir"]) / "nucleus_annotations.tsv.gz", usecols=["tool", "branch", "population", "cluster_scope", "cluster_id", "nucleus_id", "sample", "condition", "keep_for_plot", "keep_for_analysis"])
    books = read(root / "plots/PLOT_BOOK_INDEX.tsv")
    corekeys = ["tool", "branch", "population", "cluster_scope"]
    expected = {(t,b,p,"original_frozen") for t in cfg["tools"] for b in cfg["branches"] for p in cfg["populations"]}
    expected |= {(t,b,p,"fresh_branch") for t in cfg["tools"] for b in cfg["branches"][1:] for p in cfg["populations"]}
    assert set(map(tuple, books[corekeys].to_numpy())) == expected and not books.duplicated(corekeys).any()
    panels = read(Path(cfg["initial_annotation_dir"]) / "panel_membership.tsv")
    panel_rows, results = [], []
    values_count = negative_count = 0
    for row in books.itertuples():
        stage = Path(row.pdf).parent
        r = json.loads((stage / "COMPLETE.json").read_text()); idx = read(stage / "PAGE_INDEX.tsv")
        assert r["status"] == "COMPLETE" and r["n_pages"] == len(idx) == row.pages
        assert r["scientific_fitting_performed"] is False and r["raw_nucleus_count_aggregation_performed"] is False
        assert r["all_displayed_nuclei_rendered"] is True and not bool(row.coordinate_refit_performed)
        a = ann
        for k in corekeys: a = a[a[k].eq(getattr(row,k))]
        actual = read(stage / "rendered_nucleus_metadata.tsv.gz", usecols=["nucleus_id", "sample", "condition", "cluster_id", "keep_for_plot", "keep_for_analysis"])
        assert not actual.nucleus_id.duplicated().any() and set(actual.nucleus_id) == set(a.nucleus_id)
        joined = actual.merge(a[["nucleus_id", "sample", "condition", "cluster_id", "keep_for_plot", "keep_for_analysis"]], on="nucleus_id", suffixes=("_rendered", "_annotation"), validate="one_to_one")
        for k in ["sample", "condition", "cluster_id"]: assert joined[k+"_rendered"].eq(joined[k+"_annotation"]).all()
        for k in ["keep_for_plot", "keep_for_analysis"]: assert truth(joined[k+"_rendered"]).equals(truth(joined[k+"_annotation"]))
        displayed = actual[truth(actual.keep_for_plot)]
        assert len(displayed) == row.displayed_nuclei == r["displayed_nuclei"]
        assert displayed["sample"].nunique() == 12 and displayed["condition"].nunique() == 2
        native = read(stage / "displayed_cluster_counts_per_mouse.tsv")
        got = displayed.groupby(["cluster_id", "sample", "condition"]).size().rename("n").reset_index()
        test = native.merge(got, on=["cluster_id", "sample", "condition"], validate="one_to_one", how="outer")
        assert test.displayed_nuclei.eq(test.n).all()
        for path in sorted(stage.glob("*_exact_values.tsv")):
            v = read(path); f = v.fraction_detected.to_numpy(); z = v.scaled_expression_z.to_numpy(); neg = np.isfinite(z) & (z < 0)
            policy = cfg["plots"]
            expected_area = policy["dot_max_area"] * f * np.where(neg, policy["negative_z_area_multiplier"], 1)
            assert np.all((f >= 0) & (f <= 1)) and np.isfinite(f).all()
            assert np.allclose(v.glyph_area_points2, expected_area, rtol=1e-10, atol=1e-10)
            assert np.allclose(v.glyph_alpha, np.where(neg,policy["negative_z_alpha"],.9), rtol=1e-10, atol=1e-10)
            assert np.allclose(v.glyph_color_z_clipped,np.clip(z,-2.5,2.5),equal_nan=True)
            assert np.array_equal(truth(v.negative_z_visual_downweight_applied).to_numpy(),neg)
            assert np.allclose(f,v.n_detected / v.n_nuclei,rtol=1e-10,atol=1e-10)
            assert np.allclose(v.mean_log1p_normalized,v.sum_log1p_normalized/v.n_nuclei,rtol=1e-10,atol=1e-10)
            assert v.n_mice_expression_supported.ge(0).all()
            assert set(v.scale_reference_branch) == {"original_full" if row.branch == "original_full" else "paired_original"}
            values_count += len(v); negative_count += int(neg.sum())
        local_members = pd.concat([read(p) for p in sorted(stage.glob("supplied_panel_*_membership.tsv"))], ignore_index=True) if row.cluster_scope == "original_frozen" else pd.DataFrame()
        if row.cluster_scope == "original_frozen":
            keys = ["set", "group", "gene_symbol", "gene_id", "source", "selection"]
            actual_rows = set(map(tuple,local_members[keys].fillna("").to_numpy()))
            expected_rows = set(map(tuple,panels[panels.gene_id.notna()][keys].fillna("").to_numpy()))
            assert actual_rows == expected_rows
            panel_rows.append(dict(**{k:getattr(row,k) for k in corekeys}, mapped_panel_rows=len(expected_rows),unique_panel_genes=local_members.gene_id.nunique(),full_supplied_panel_membership="PASS"))
        else: assert not list(stage.glob("*_exact_values.tsv"))
        results.append(dict(**{k:getattr(row,k) for k in corekeys},displayed_nuclei=len(displayed),n_clusters=displayed.cluster_id.nunique(),n_pages=len(idx),metadata_and_masks="PASS",no_coordinate_refit="PASS",native_all_mouse_views="PASS"))
    out.mkdir(parents=True)
    pd.DataFrame(results).to_csv(out / "BOOK_CHECKS.tsv",sep="\t",index=False)
    pd.DataFrame(panel_rows).to_csv(out / "SUPPLIED_PANEL_COVERAGE.tsv",sep="\t",index=False)
    result=dict(status="PASS",stage="completed_native_plot_artifact_review",n_books=len(results),n_primary_books=sum(x["cluster_scope"]=="original_frozen" for x in results),
                n_pages=int(books.pages.sum()),n_hash_checked_artifacts=checked,n_dotplot_values_checked=values_count,n_negative_z_values_checked=negative_count,
                native_input_receipt_sha256=sha(receipt_path),annotation_receipt_sha256=sha(Path(cfg["final_annotation_dir"])/"COMPLETE.json"),source_sha256=sha(__file__),
                visual_review_scope="SEPARATE_MANUAL_PNG_REVIEW",browser_interaction_QA="UNVERIFIED_URL_SECURITY_POLICY",scientific_fitting_performed=False,raw_nucleus_RNA_count_aggregation_performed=False)
    (out / "COMPLETE.json").write_text(json.dumps(result,indent=2)+"\n");print(json.dumps(result))


if __name__ == "__main__": main()
