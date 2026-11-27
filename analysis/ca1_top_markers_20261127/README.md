# Seurat CA1 top-two marker dotplot

Published reader: `site/ca1_sublayers/top_markers/`. The same component is embedded near the top of the existing CA1 sublayer reader and linked from the main atlas.

## Selection and interpretation

The first two entries by `representative_rank` were copied from the existing `representative_markers.tsv`, restricting to Seurat, original_full counts, CA1_ProS and original_frozen clusters. These are descriptive, specificity- and mouse-support-screened representatives within the saved curated marker panel, not genome-wide DEGs or newly fitted marker tests. The chosen genes stay fixed when switching expression display between pooled, WT and cKO nuclei.

C7 and its exploratory genes are excluded at the user's request. The ten retained original clusters each supply two supported representatives (20 unique genes). C8 is CA1-prosubiculum; C11 remains identity-conflicting. The original native assignments and underlying expression values are unchanged.

Rows run C1, C2, C3, C4, C5, C6, C8, C9, C10, C11 from top to bottom. Gene columns are grouped by their pooled peak mean-expression cluster in that same ascending order. The pooled column order is fixed across WT/cKO views and when rows are hidden, so comparisons retain aligned columns. Selection origin and expression peak are separate fields: a marker selected for C4 or C5 can peak in a smaller other cluster. The evidence table and gene-list export retain both. The optional original-selection order restores two-gene source-cluster blocks.

Dot area is proportional to the saved nonzero fraction. Colour uses each saved condition mean standardized against the same saved pooled gene-wise mean and SD, with an adjustable clipped display range initially -0.7 to +2. This shared display scale is stored in display_scaled_expression_z; original scaled_expression_z values remain unchanged, including the unavailable native cKO z-scores. The reference includes all original 11 clusters and is not recalculated after excluding C7. No negative-z area suppression, row-dependent rescaling, UMAP change, annotation change or model refit is applied. Low or zero detection is distinct from unavailable data. The displayed nuclei total 26,733 (14,360 WT and 12,373 cKO); the original 27,703-nucleus CA1/ProS population remains unchanged.

## Controls and downloads

The component supports pooled/WT/cKO expression, colour limits and three colour stops, dot radius and opacity, gene and cluster label size, gene removal, row selection, expression-block visibility and pooled-peak versus original-selection order. Downloads include vector SVG/PDF, 3x PNG, selected values TSV, ranked gene-list TSV, and settings JSON with validated reload. Default SVG/PDF files are also stored under the reader's `figures/` directory.

## Sources and validation

`selection_and_source_audit.json` binds the source files by SHA256. The source summaries are under the original 20261014 v2 run's `results/marker_summaries/`; the existing marker shard supplies every displayed value. For each of three conditions, 200 cluster-by-gene cells are copied, with 8,400 source fields checked for equality. The original 3,140 atlas shard descriptors and source files remain unchanged. The added display z values are arithmetic standardizations of existing summaries, not newly fitted normalization or differential-expression models.

Run the compact data checks with `node tests/ca1_top_markers.cjs`. Local end-to-end controls/export checks use `node tests/ca1_top_markers.cjs --browser`; its Playwright location can be provided via `PLAYWRIGHT_MODULE`. The browser acceptance record includes genotype nucleus counts, settings roundtrip, all downloads, and the parent-page embed. SVG/PNG were inspected at default and large labels; PDF was independently rendered and inspected after fixing rotated right-aligned labels.

## Pages size budget

The existing site was at its 950 MB guard. The byte-identical older `ca1_sublayers_offline_20261031.zip` has moved to this folder's `archives/` directory, outside the Pages bundle. Its link in the CA1 reader now points to the repository's raw download. `archive_relocation.json` records its unchanged SHA256 and public path. The size guard is unchanged, and the archive, original reader results, workflows and scientific payloads are retained.
