# Seurat CA1 top-two marker dotplot

Published reader: `site/ca1_sublayers/top_markers/`. The same component is embedded near the top of the existing CA1 sublayer reader and linked from the main atlas.

## Selection and interpretation

The first two entries by `representative_rank` were copied from the existing `representative_markers.tsv`, restricting to Seurat, original_full counts, CA1_ProS and original_frozen clusters. These are descriptive, specificity- and mouse-support-screened representatives within the saved curated marker panel, not genome-wide DEGs or newly fitted marker tests. The chosen genes stay fixed when switching expression display between pooled, WT and cKO nuclei.

Ten clusters have two supported candidates. C7 has zero candidates passing the original criteria. Its two exploratory entries, Ndrg2 and Gm19951, come from the existing candidate audit sorted by descending detected-fraction difference versus the rest, then mean log1p-expression difference, then stable gene ID. Their weak specificity and unsupported status are disclosed in the chart, evidence table and gene download. The supported-only control removes them without hiding the C7 row. C8 is CA1-prosubiculum; C11 remains identity-conflicting, without changing the original assignments.

All 11 original cluster rows and 22 candidate genes are present by default. Dot area is proportional to the saved nonzero fraction. Colour is the saved gene-wise expression z-score with an adjustable clipped display range, initially -0.7 to +2. No negative-z area suppression, row-dependent rescaling, UMAP change, annotation change or model refit is applied. Low or zero detection is not an absent-data estimate; missing values are shown with crosses.

## Controls and downloads

The component supports pooled/WT/cKO expression, colour limits and three colour stops, dot radius and opacity, gene and cluster label size, gene removal, row selection, selection-block visibility, and supported-only filtering. Downloads include vector SVG/PDF, 3x PNG, selected values TSV, ranked gene-list TSV, and settings JSON with validated reload. Default SVG/PDF files are also stored under the reader's `figures/` directory.

## Sources and validation

`selection_and_source_audit.json` binds the source files by SHA256. The source summaries are under the original 20261014 v2 run's `results/marker_summaries/`; the existing marker shard supplies every displayed value. For each of three conditions, 242 cluster-by-gene cells are copied, with 10,164 source fields checked for equality. The original 3,140 atlas shard descriptors and source files remain unchanged.

Run the compact data checks with `node tests/ca1_top_markers.cjs`. Local end-to-end controls/export checks use `node tests/ca1_top_markers.cjs --browser`; its Playwright location can be provided via `PLAYWRIGHT_MODULE`. The browser acceptance record includes genotype nucleus counts, settings roundtrip, all downloads, and the parent-page embed. SVG/PNG were inspected at default and large labels; PDF was independently rendered and inspected after fixing rotated right-aligned labels.

## Pages size budget

The existing site was at its 950 MB guard. The byte-identical older `ca1_sublayers_offline_20261031.zip` has moved to this folder's `archives/` directory, outside the Pages bundle. Its link in the CA1 reader now points to the repository's raw download. `archive_relocation.json` records its unchanged SHA256 and public path. The size guard is unchanged, and the archive, original reader results, workflows and scientific payloads are retained.
