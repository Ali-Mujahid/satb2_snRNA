# Seurat CA1 and inhibitory marker editors

Readers:
- `site/ca1_sublayers/top_markers/` (retains the existing CA1 URL and embedded CA1 section).
- `site/inhibitory_markers/` (new inhibitory reader).
- Both load `site/marker_explorer/app.js` and the same namespaced stylesheet.

## Evidence and scope

The exact 88 CA1 and 73 inhibitory screenshot gene names map uniquely to existing source stable IDs. The ordered lists are in `scripts/build_marker_explorers.py`, mapping JSONs in this folder, and each reader's `marker_selection.tsv`. No genes were inferred from image dot colours. The original Seurat `original_full / original_frozen` summaries supply expression and detected fractions. The base marker shard and already-published `atlas-v3.4.0` extension provide 1,055 cached genes per population. The reader verifies both compressed and decoded SHA256 values before plotting. No new counts, DESeq2 fits, cluster assignments or embeddings are produced.

All 11 CA1 clusters are available, with C7 hidden initially to preserve the preceding user preference; checking C7 or Add all clusters restores it. The C8 label remains CA1-prosubiculum. All 20 original inhibitory clusters are available and visible by default. Inhibitory subtype labels are shortened working annotations; the parent C6 is not silently renamed C6.0 or split.

The screenshot gene set is selected by default. The earlier CA1 top-two panel is retained as a preset, and a second preset uses saved supported representatives for inhibitory clusters where available. The complete cached catalog can be searched; typed genes must resolve to a source symbol or stable ID, and unavailable names are reported rather than plotted as zeros.

## Encoding

Dot area is proportional to the measured nonzero-nucleus fraction. Colour can show shared-reference expression z, mean log1p CP10k, or detection. Shared-reference z standardizes the saved condition mean against the existing pooled gene reference; source scaled-expression fields, including unavailable native cKO z values, remain unchanged in value exports. The original reference includes all 11 CA1 or 20 inhibitory clusters. Display filtering never recomputes it. These pooled displays are not mouse-level genotype tests.

Dotplots and pheatmap-style heatmaps share genes, row order and colour scales. Heatmaps use coloured tiles without a detection-size encoding. Peak-expression ordering is fixed from pooled data; screenshot/manual and alphabetical orders are available. Existing classifications remain unchanged. The heatmap is an ordered matrix display, not a new clustering analysis.

## Requested controls

Both readers support gene and cluster removal/restoration, all-screenshot and representative presets, extra cached genes and bulk paste, three-stop colour palettes, measurement and colour limits, label size, row label style, column spacing, dot radius, mark opacity, plot borders, dot outlines, heatmap cell borders, and horizontal/vertical guides. Guide lines pass through dot centres or form cell grids in heatmap view; they are presentation guides, not biological interactions. Their colours, widths and opacity are editable.

SVG and PNG export the displayed gene page; vector PDF exports all selected gene pages. Pagination keeps large selections readable. PNG resolution is limited to safe browser canvas dimensions; vector formats preserve full detail. TSV exports contain all selected genes/rows regardless of pagination, including original values and explicitly named display fields. Settings JSON validates dataset identity, settings, gene IDs and cluster IDs. The preceding CA1 v2 settings are also accepted.

## Validation

`node tests/marker_explorers.cjs` checks both configs, exact screenshot mapping, all original source hashes, restored CA1 C7 availability, and 7,284 screenshot-gene by cluster by condition source cells. The existing Pages workflow calls this through `tests/ca1_top_markers.cjs`.

`node tests/marker_explorers.cjs --browser` tests cluster removal/restoration, gene removal/restoration and paste, border/outline/guide operations, heatmap switching, cKO views, all downloads, settings roundtrip, multi-page PDF export, and the parent CA1 embed. The browser test substitutes byte-identical local extension shards for the pinned remote URL; live checks verify the actual remote transport after deployment. `PLAYWRIGHT_MODULE` can override the local runtime location.

Browser exports and QA evidence are retained in `browser_checks/`. Independent PDF rendering checks complete gene labels and page counts. Source-bound atlas validation and the unchanged 950 MB Pages guard are run before publication. Earlier standalone top-two files and their source audit remain in `analysis/ca1_top_markers_20261127/` and Git history.
