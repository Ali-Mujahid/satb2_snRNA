Current hosted update: **[atlas v3.7.0](https://ali-mujahid.github.io/satb2_snRNA/)** — adds the **[hdWGCNA and curated TF target-signature explorer](https://ali-mujahid.github.io/satb2_snRNA/regulatory/)**, mouse-level scores, native figure gallery, adjustable networks and offline/PDF exports. The **[CA1–inhibitory interaction explorer](https://ali-mujahid.github.io/satb2_snRNA/communication/)** retains strict DEGs at both exact endpoints with signed ligand/receptor effects. Native scores, original gene fits and prior analyses are preserved. [New findings and limitations](analysis/reviews/20261026_regulatory_findings.md) · [Corrected communication results](docs/DEG_ONLY_RESULTS.md).

The regulatory release includes **10 primary CA1 modules, 5 inhibitory modules and 23 primary TF target-signature associations**. None of the TF associations passes the smaller matched-partition sensitivity FDR; SATB2 activity is not FDR-significant. The pinned v3.7 release records the earlier unavailable motif references. Matching official mouse references became reachable on 28 October; a separate [native pySCENIC run and source snapshot](analysis/regulatory_20261028/README.md) are now in progress, with no completed motif inference claimed yet. [Regulatory figures and source downloads](https://ali-mujahid.github.io/satb2_snRNA/regulatory/downloads.html) · [Standalone regulatory HTML](https://github.com/Ali-Mujahid/satb2_snRNA/releases/download/atlas-v3.7.0/SATB2_regulatory_interactive_offline.html).

The [gene-level communication/module/TF interpretation addendum](analysis/reviews/20261028_completed_LR_HD_TF_integration_guide.md) traces five selected LR examples to original DESeq2 signs, native module membership and actual saved TF-prior links. These are joins of completed results, without new fits or changes to the full strict communication datasets.

# SATB2 single-nucleus RNA-seq atlas

An interactive Scanpy–Seurat v5 atlas of six WT and six Satb2 cKO mice, with
MapMyCells reference annotation and saved mouse-level DESeq2 comparisons.

**[Open the web atlas](https://ali-mujahid.github.io/satb2_snRNA/)** ·
**[Data, figures and source downloads](https://ali-mujahid.github.io/satb2_snRNA/downloads.html)** ·
**[Portable atlas files](https://github.com/Ali-Mujahid/satb2_snRNA/releases/tag/atlas-v3.0.0)**

The website runs in the browser on other computers without R, Python or a
local server. It supports annotated/gene-colored UMAPs, 60–180-gene dotplot
pages, adjustable expression/z color limits, label editing, cluster exclusions,
mouse expression boxplots, exact saved DESeq2 results and PDF downloads.
Labels and selections stay in that browser; export settings to transfer them.

The v3 marker tools add a **Screenshot · all 75 genes** button, preserving the
image's left-to-right order and fifteen blocks of five in both tools. The preset
places all 75 genes on one 90-gene page. Genes can also be grouped by cluster
expression pattern, ranked for a focus cluster, pasted in a custom order or
moved individually. Scanpy and Seurat share the chosen column order using one
selected reference tool; their cluster row orders can be adjusted separately.
See the [marker preset and ordering guide](docs/marker_preset_v3.md).

The v3.1 web controls add individual-nucleus gene expression, separate category
and expression color palettes, compact/custom display labels, reverse/A–Z
marker ordering, and volcano axes with adjustable DEG labels and colors.
See the [v3.1 control guide](docs/ATLAS_V3_1_CONTROLS.md). These new web controls
do not change the preserved v3.0.0 portable release archives.

The v3.2 volcano controls make non-DEGs visible with stronger light-gray dots
and add a full-range overview plus a separate low-y detail panel. The shared
view selector, low-y maximum, gray-dot size and opacity apply to both tools.
PDF downloads include the selected full/low/both view; combined plots use a
portrait page. These display changes preserve every saved DESeq2 value.

## Study and interpretation

Views cover all nuclei, excitatory nuclei, inhibitory nuclei, CA1/ProS
candidates and CA1 layer-associated RNA programs. Marker dotplots and gene
UMAPs use **889** marker genes: the original 878 summaries plus eleven added
descriptive summaries from the same frozen native integer-count populations.
The original marker values, source memberships and saved fits are unchanged.
Automatic gene ordering describes expression contrasts; it is not a new marker
test or evidence that screenshot cluster numbers match this study's clusters.
Mouse boxplots can select
all 32,285 source genes and use the selected exact endpoint's recorded size
factors. Each boxplot point is an independent mouse.

Gene UMAP resolution is explicit. **Individual nuclei** uses each nucleus's
exact existing same-branch source count or
`log1p(10000 × gene count / complete source library)` for all 889 panel genes.
**Saved cluster summaries** retains the prior cluster-mean/detection/z views;
uniform colors in that mode do not represent individual-nucleus RNA.
Both tools join the same branch counts by exact stable nucleus ID while
retaining their existing coordinates, populations and exclusions. Source
count zero is distinct from an unavailable source nucleus. The
[gene-export provenance](site/data/GENE_EXPRESSION_PROVENANCE.json) records
the three native count-object hashes, full-library checks and ID validation.
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

* `site/`: GitHub Pages frontend, verified original interactive data and the
  separately documented marker and individual-nucleus expression extensions.
* `gene_expression/`: sparse source-gene vectors and stable nucleus/library
  axes, fetched from the preserved `atlas-v3.1.0` tag with SHA-256 checks.
* `data_downloads/`: exact source-download payloads, outside the Pages size
  budget, fetched from the preserved `atlas-v2.1.0` tag with SHA-256 checks.
* `analysis/`: original and revision workflow/configuration snapshots, marker
  workbooks, annotation/reference evidence, package records and provenance.
  `analysis/dashboard/v3_marker_extension/` includes the new descriptive sums,
  source hashes, ordered marker panel and validation evidence.
* `docs/analysis/`: reproduction and processed-object guidance.
* `scripts/`: transport, integrity and release packaging helpers.

Large portable archives, figure PDFs and compatible processed objects are
release assets rather than ordinary Git blobs. The 28-object processed set
contains 12 Scanpy AnnData files, 13 Seurat HVG-only representation RDS files,
and three shared full-gene integer-count AnnData files. Seurat representations
need the matching count object and final masks/annotation for count-based
analyses; these are not falsely described as full-gene standalone R objects.

## Offline use

Download the **recommended v3 runtime** `SATB2_atlas_v3_runtime.zip` and
`SATB2_atlas_v3_saved_results.zip` from
[atlas-v3.0.0](https://github.com/Ali-Mujahid/satb2_snRNA/releases/tag/atlas-v3.0.0).
Extract both into the same parent folder, allowing their common `offline_dashboard_v3/`
folder to merge, then open its `index.html`. They are independent ordinary
ZIPs, not spanned archive volumes.

The v3 runtime includes the 75-gene preset, eleven additional marker summaries,
shared gene-column ordering and manual row/gene ordering. It retains the v2.1
tiny-P-value and selected-overlay-gene corrections and the classic offline
data loader. The saved results companion places the existing source results
under the common v3 folder; it does not refit the analysis. The earlier
[v2/v2.1 archives and original assets remain available](https://github.com/Ali-Mujahid/satb2_snRNA/releases/tag/atlas-v2.1.0).

Together, the v3 runtime and saved-results companion retain
the original source tables, 67 PDF books and 450 individually indexed figure
pages. The new runtime's marker-extension provenance distinguishes newly
aggregated display summaries from preserved source results. The prior v2.1
[UI update provenance](docs/offline_ui_update.json) remains available.

`SATB2_analysis_source_v2.zip` supplies an additional source/reference snapshot.
The website does not need processed AnnData/RDS files to operate. Those files
are separate downloads for R/Scanpy work. Original CellRanger/raw-droplet inputs
and large external reference models are not recreated by this repository;
follow the source/input provenance for a complete scientific rerun.

## Validation and deployment

Original analysis outputs and the offline v1/v2 deliveries were preserved.
The hosted transport verifies compressed and decoded SHA-256 and preserves
original numeric values, nulls and exact source text. The authorized v3 marker
extension computes eleven descriptive gene summaries from existing counts,
without changing the original 878 genes or fitting embeddings, clusters,
annotations or differential models. All 64,548 added gene × atom × mouse
integer sums match the 489 existing full-gene native count arrays exactly.
Validation evidence is under `docs/` and
`analysis/provenance/`; the hosted browser checks are recorded separately.

The v3.1 individual-nucleus exporter reads the three existing full-gene count
objects without modifying them. It writes 2,667 gene vectors (889 × three
branches) and three nucleus/library axes; it does not rerun QC, CellBender,
normalization-model fitting, clustering, embeddings, annotation or DESeq2.
The 3,140 prior base-atlas shard descriptors and source fingerprint are unchanged.

GitHub Pages deploys the prebuilt `site/` directory. Source-download URLs are
pinned to `atlas-v2.1.0`; preserve that tag. Publication helpers verify data
integrity and keep the published site below its size limit.

Vendored software retains its existing licenses. No additional blanket reuse
license is granted to the study data or other source materials by this README.


Hosted compatibility patch `atlas-v3.0.1` gives each compressed data request a content-specific cache URL, so browsers with older marker files load the correct extended payload. The v3.0.0 portable archives and all scientific data are unchanged.

The v3.3 presentation controls add optional compact layouts, per-cluster placement, rotation/reflection, a clearly marked hippocampus-like schematic, reversible individual dot styling, and clean volcano PDFs without label leader lines. See [the presentation guide](docs/ATLAS_V3_3_PRESENTATION.md). Original coordinates, measurements, memberships and DESeq2 results remain unchanged.
