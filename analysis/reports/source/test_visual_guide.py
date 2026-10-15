#!/usr/bin/env python3
"""Synthetic renderer checks only; no scientific input or model fitting."""
from pathlib import Path
import argparse
import base64
import json
import tempfile
import pandas as pd
from build_visual_guide import Evidence, sha, boolean, write_html, make_crosswalks
import matplotlib.pyplot as plt


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--outdir", type=Path, required=True)
    out = ap.parse_args().outdir
    if out.exists(): raise ValueError("Synthetic fixture directory already exists")
    out.mkdir(parents=True)
    fig, ax = plt.subplots(figsize=(8, 4)); ax.text(.5, .5, "SYNTHETIC LAYOUT FIXTURE\nNot scientific results", ha="center", va="center", fontsize=22); ax.axis("off")
    fig.savefig(out / "fixture.png", dpi=100); fig.savefig(out / "fixture.pdf"); fig.savefig(out / "fixture.svg"); plt.close(fig)
    image = "data:image/png;base64," + base64.b64encode((out / "fixture.png").read_bytes()).decode()
    artwork = {"thumbnail": image, "png": "fixture.png", "pdf": "fixture.pdf", "svg": "fixture.svg"}
    books = []
    for b in ["original_full", "paired_original", "paired_cellbender"]:
        for t in ["scanpy", "seurat"]:
            for p in ["allcells", "Excitatory", "Inhibitory", "CA1_ProS"]:
                for scope in ["original_frozen", "fresh_branch"]:
                    if b == "original_full" and scope == "fresh_branch": continue
                    books.append(dict(tool=t, branch=b, population=p, population_label=p, cluster_scope=scope,
                                      title="SYNTHETIC " + t + " " + b, displayed_nuclei=10, thumbnail=image,
                                      preview_png="fixture.png", cluster_png="fixture.png", mouse_png="fixture.png", pdf="fixture.pdf",
                                      previews=[dict(label="Reviewed UMAP", **artwork), dict(label="Representative genes", **artwork)],
                                      pages=[dict(label="1 · synthetic page", png="fixture.png")]))
    data = dict(books=books, programs=[], crosswalks=[], effects=[], annotations=[], retention=[],
                annotation_tsv="synthetic_annotations.tsv", nucleus_tsv="synthetic_nuclei.tsv.gz", pipeline=artwork,
                technical=artwork, programme_books=[], program_books=[], crosswalk_book="fixture.pdf",
                agreement=dict(narrative="SYNTHETIC agreement fixture only",cards=[],review="synthetic_annotations.tsv",semantic_table="synthetic_annotations.tsv",overlap_table="synthetic_annotations.tsv"),
                provenance="SYNTHETIC FIXTURE ONLY: no scientific results or fits.")
    write_html(out, data)
    text = (out / "visual_guide.html").read_text().replace("SATB2 · Visual results guide", "SYNTHETIC visual guide test")
    (out / "visual_guide.html").write_text(text)
    pd.DataFrame([{"synthetic": True}]).to_csv(out / "synthetic_annotations.tsv", sep="\t", index=False)
    pd.DataFrame([{"synthetic": True}]).to_csv(out / "synthetic_nuclei.tsv.gz", sep="\t", index=False)
    assert boolean(pd.Series(["True", "false", 1, "0"])).tolist() == [True, False, True, False]
    try: boolean(pd.Series(["unknown"]))
    except ValueError: pass
    else: raise AssertionError("Invalid mask accepted")
    with tempfile.TemporaryDirectory() as tmp:
        r = Path(tmp); (r / "file.tsv").write_text("value\n1\n")
        (r / "COMPLETE.json").write_text(json.dumps({"status": "COMPLETE", "output_sha256": {"file.tsv": sha(r / "file.tsv")}}))
        e = Evidence(); e.stage(r, "COMPLETE"); e.bind(r / "file.tsv")
        (r / "file.tsv").write_text("value\n2\n")
        try: e.bind(r / "file.tsv")
        except ValueError: pass
        else: raise AssertionError("Changed producer artifact accepted")
    receipt = dict(status="PASS", scope="SYNTHETIC_ONLY", n_synthetic_books=len(books),
                   tests=["invalid_mask_rejected", "producer_hash_drift_rejected", "offline_HTML_payload_generated"],
                   source_sha256=sha(Path(__file__).with_name("build_visual_guide.py")), scientific_fitting_performed=False)
    (out / "SYNTHETIC_TEST_RECEIPT.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt))


if __name__ == "__main__": main()
