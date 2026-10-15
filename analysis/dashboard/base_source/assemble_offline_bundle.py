#!/usr/bin/env python3
"""Package the offline viewer, without fitting or changing scientific results."""
from pathlib import Path
import argparse
import datetime as dt
import hashlib
import json
import re
import shutil
import zipfile

BASE = Path(__file__).resolve().parents[1]
SOURCE = BASE / "offline_dashboard_source"
DEFAULT_OUT = BASE / "reports" / "offline_dashboard_v1"

README = """# SATB2: offline Scanpy–Seurat comparison

Unzip the complete folder, then double-click **index.html** in the unzipped folder.
Keep index.html, app.js, styles.css, data/, vendor/ and assets/ together. The page
uses local files and needs no web server, installation, account, or internet.

Use the four population tabs to explore all nuclei, excitatory nuclei,
inhibitory nuclei, and CA1–ProS candidates. Scanpy and Seurat views are shown
side by side. Select a saved embedding, annotation, marker set and palette;
include or exclude clusters using each tool's controls. The source analyses
and their original labels/counts remain unchanged. Edited display labels,
colors and selections can be saved as a preset and exported for review.

The main population masks are the completed analysis defaults: clusters below
100 nuclei were excluded by a shared iterative mask in the all-nucleus and
excitatory views. Inhibitory and CA1–ProS groups have no such size cutoff.
Hiding clusters does not recompute PCA, Harmony, UMAP, clustering or markers.

Marker dot plots use already-produced full-cohort statistics for 878 available
gene IDs, with gene symbols for display. Provided marker-sheet panels and
representative-marker presets are included. A gene absent from that panel has
no dot-plot statistics; it is not represented as a zero. Dot area represents
detection fraction and color represents the chosen expression metric. Negative
scaled values are smaller and fainter by default. Gene overlays on UMAP show
cluster-mean expression, not individual-nucleus gene expression.

The differential-expression viewer contains native DESeq2 outputs only.
Positive log2 fold change means cKO greater than WT. The exported result fields
are gene, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj; gene uses the source
symbol. Stable IDs, coverage and model settings are available separately for
auditing. NA adjusted p-values remain NA. Significance is within each fitted
endpoint, not across every population explored in the viewer.

An unchanged cluster retains its own precomputed fit when other clusters are
hidden. A selected total or edited annotation group can show a precomputed fit
only if its exact members match an existing fitted population. Otherwise the
page states **NO MATCHING FIT** and lets you download a request describing the
selection. An offline page cannot run DESeq2. It never combines p-values,
averages fold changes, or presents the old whole-population fit as a new fit.

Inhibitory pseudobulk includes every available nucleus in each cluster. All
nonempty positive-library mice are eligible, with at least three independent
WT and three cKO mice required. Contributions below five nuclei are flagged as
exploratory. A retained cluster is not automatically a reliable cell identity.
CA1 superficial/deep labels describe marker-expression programs; these data
do not determine physical anatomical layers. Rare/mixed/conflicting groups
remain explicit in the annotation and coverage views.

Use a plot's PDF button for the current custom view. Native precomputed PDF
books are separate downloads and retain their source settings. Save presets
as JSON before moving the folder to another computer; browser local storage
may differ between computers or be unavailable for file:// pages.

Validation information is in DASHBOARD_REVIEW.md and BUNDLE_MANIFEST.json.
The viewer was checked through source, local-file integrity, simulated-DOM
interactions and exported artifacts. Interactive browser inspection was
unavailable because the automation browser rejected local file:// navigation.
That restriction was respected; no server or alternate browser was used.

This separate delivery does not overwrite any earlier analysis. Original
source paths and SHA-256 hashes are retained in the provenance records. The
vendored JavaScript libraries and their licenses are in vendor/.
"""

def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--zip", action="store_true")
    args = parser.parse_args()
    out = args.out.resolve()
    if "runs" in out.relative_to(BASE).parts or out == BASE:
        raise SystemExit("Viewer packaging must use a separate reports subfolder")
    out.mkdir(parents=True, exist_ok=True)
    for name in ("index.html", "app.js", "styles.css"):
        src = SOURCE / name
        if not src.is_file():
            raise SystemExit(f"Missing frontend: {name}")
        shutil.copy2(src, out / name)
    shutil.copytree(SOURCE / "vendor", out / "vendor", dirs_exist_ok=True)
    (out / "README.md").write_text(README)
    html = (out / "index.html").read_text()
    refs = re.findall(r'(?:src|href)=["\']([^"\']+)', html)
    errors = []
    for ref in refs:
        if ref.startswith(("http:", "https:", "//")) and re.search(
            r'<(?:script|link)\b[^>]*(?:src|href)=["\']' + re.escape(ref), html
        ):
            errors.append("External runtime dependency: " + ref)
        elif not ref.startswith(("http:", "https:", "#", "data:", "mailto:")):
            local = out / ref.split("#")[0].split("?")[0]
            if not local.exists():
                errors.append("Missing HTML reference: " + ref)
    if errors:
        raise SystemExit("\n".join(errors))
    files = [p for p in sorted(out.rglob("*")) if p.is_file() and p.name != "BUNDLE_MANIFEST.json"]
    manifest = {
        "schema": "satb2-offline-bundle-1",
        "generated_at": dt.datetime.now(dt.timezone.utc).isoformat(),
        "scientific_fitting_performed": False,
        "raw_nucleus_count_aggregation_performed": False,
        "previous_analysis_files_modified": False,
        "entrypoint": "index.html",
        "file_count": len(files),
        "bytes": sum(p.stat().st_size for p in files),
        "files": [{"path": str(p.relative_to(out)), "bytes": p.stat().st_size,
                   "sha256": sha(p)} for p in files],
    }
    (out / "BUNDLE_MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\n")
    if args.zip:
        dest = out.with_suffix(".zip")
        with zipfile.ZipFile(dest, "w", compression=zipfile.ZIP_DEFLATED,
                             compresslevel=6, allowZip64=True) as archive:
            for path in sorted(out.rglob("*")):
                if path.is_file():
                    archive.write(path, Path(out.name) / path.relative_to(out))
        print(json.dumps({"zip": str(dest), "bytes": dest.stat().st_size,
                          "sha256": sha(dest)}))
    print(json.dumps({"folder": str(out), "files": len(files), "bytes": manifest["bytes"]}))

if __name__ == "__main__":
    main()
