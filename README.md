# SATB2 single-nucleus RNA-seq atlas

An interactive Scanpy–Seurat v5 atlas of six WT and six Satb2 cKO mice, with
MapMyCells reference annotation and saved mouse-level DESeq2 comparisons.

**[Open the web atlas](https://ali-mujahid.github.io/satb2_snRNA/)** ·
**[Data, figures and source downloads](https://ali-mujahid.github.io/satb2_snRNA/downloads.html)** ·
**[Release assets](https://github.com/Ali-Mujahid/satb2_snRNA/releases/tag/atlas-v2.1.0)**

The website runs in the browser on other computers without R, Python or a
local server. It supports annotated/gene-colored UMAPs, 60–180-gene dotplot
pages, adjustable expression/z color limits, label editing, cluster exclusions,
mouse expression boxplots, exact saved DESeq2 results and PDF downloads.
Labels and selections stay in that browser; export settings to transfer them.

## Study and interpretation

Views cover all nuclei, excitatory nuclei, inhibitory nuclei, CA1/ProS
candidates and CA1 layer-associated RNA programs. Marker dotplots and gene
UMAPs use 878 already-produced cluster summaries. Mouse boxplots can select
all 32,285 source genes and use the selected exact endpoint's recorded size
factors. Each boxplot point is an independent mouse.

Gene UMAP colors are cluster summaries, not measured per-nucleus expression.
The optional inhibitory spacing control changes display geometry only.
Seventy original cluster records across overlapping views retain uncertain
fine identities, with nonbinding parent/“-like” label options. Deep/superficial
CA1-associated programs do not establish physical anatomical layers.

The saved model index has 948 PRIMARY endpoints (852 tested, 96 unsupported)
and 112 partition/leave-one-mouse-out endpoints. DESeq2 exports retain the
standard seven columns and gene symbols, with stable IDs in separate audits.
An edited total/union only reuses a fit when the exact nucleus membership,
tool, count branch and endpoint match. Both tools' fits use one DESeq2 engine;
identical shared-count results are not independent statistical replication.

## Repository layout

* `site/`: GitHub Pages frontend and losslessly transported interactive data.
* `data_downloads/`: exact source-download payloads, outside the Pages size
  budget, fetched from this release's pinned Git tag with SHA-256 checks.
* `analysis/`: original and revision workflow/configuration snapshots, marker
  workbooks, annotation/reference evidence, package records and provenance.
* `docs/analysis/`: reproduction and processed-object guidance.
* `scripts/`: transport, integrity and release packaging helpers.

Large portable archives, figure PDFs and compatible processed objects are
release assets rather than ordinary Git blobs. The 28-object processed set
contains 12 Scanpy AnnData files, 13 Seurat HVG-only representation RDS files,
and three shared full-gene integer-count AnnData files. Seurat representations
need the matching count object and final masks/annotation for count-based
analyses; these are not falsely described as full-gene standalone R objects.

## Offline use

Download the **recommended corrected runtime** `SATB2_atlas_v2_1_runtime.zip`
and the unchanged companion `SATB2_atlas_v2_saved_results.zip`. Extract both
into the same parent folder, allowing their common `offline_dashboard_v2/`
folder to merge, then open its `index.html`. They are independent ordinary
ZIPs, not spanned archive volumes.

The v2.1 runtime corrects tiny P-value scientific-notation display and retains
the selected overlay gene when using the gene-search controls. Scientific
payloads, model results and the classic offline data loader are unchanged.
The original `SATB2_atlas_v2_runtime.zip` remains available as a **legacy
original**. Use the recommended v2.1 runtime for corrected offline exploration.

Together, the corrected runtime and unchanged saved-results companion retain
the original source tables, 67 PDF books and 450 individually indexed figure
pages. The combined bundle manifest updates the two UI-file hashes while
preserving all other file hashes. See [UI update provenance](docs/offline_ui_update.json)
and [the derivative archive builder](scripts/build_corrected_offline_runtime.py).

`SATB2_analysis_source_v2.zip` supplies an additional source/reference snapshot.
The website does not need processed AnnData/RDS files to operate. Those files
are separate downloads for R/Scanpy work. Original CellRanger/raw-droplet inputs
and large external reference models are not recreated by this repository;
follow the source/input provenance for a complete scientific rerun.

## Validation and deployment

Original analysis outputs and the offline v1/v2 deliveries were preserved.
The hosted transport keeps native JSON bytes and exact source text, verifies
compressed and decoded SHA-256, and performs no new RNA extraction, embedding
or differential model fitting. Validation evidence is under `docs/` and
`analysis/provenance/`; the hosted browser checks are recorded separately.

GitHub Pages deploys the prebuilt `site/` directory. Source-download URLs are
pinned to `atlas-v2.1.0`; preserve that tag. Publication helpers verify data
integrity and keep the published site below its size limit.

Vendored software retains its existing licenses. No additional blanket reuse
license is granted to the study data or other source materials by this README.

