# Screenshot marker extension (15 October 2026)

This adds **11 descriptive expression markers** to the existing 878-gene atlas, for 889 available dotplot/cluster-mean UMAP genes. The ordered screenshot preset contains 75 genes, arranged in the supplied 15 blocks of five. The display-order preset is a user-supplied marker panel; it does not assign corresponding cluster numbers to this study.

Added symbols: `Unc5c Cpne5 Nrxn1 Tcf4 Fgf12 Ptpre Gria4 Gria1 Ctnnd2 Mgat4c Ptprd`. Each maps uniquely to one native stable gene ID across all three source count branches. Their exact mappings, the complete ordered 75-gene list, and ordered five-gene blocks are in `manifest_patch.json`.

The extension reads existing **full-gene integer-count AnnData objects**, rather than the Seurat HVG-only representation. It uses the same exact frozen native atomic cluster memberships and existing shared population masks as the completed analysis. All 24 combinations of two tools, three count branches and four source populations are covered, with ALL, WT and cKO display scopes. CA1 programme/sublayer display groupings reuse these same CA1-ProS source atoms.

For each nucleus, expression is `log1p(10000 × gene count / complete selected-source library)`. Cluster means and detection fractions use all nuclei in that frozen atom; per-mouse sums and coverage are retained explicitly. The library denominator is recomputed from all 32,285 source genes and checked exactly against each nucleus's frozen verified source library. CellBender uses the corrected full-library denominator.

The original native Python summary functions calculate the extension's ALL and WT contrasts, mouse support and z-scores. Z-scores give equal weight to retained cluster means, use population SD (`ddof=0`), and corrected paired CellBender uses the paired-original reference. cKO is derived from the ALL-minus-WT numerators and denominators; **cKO z-scores stay unavailable**, following the original atlas contract. No unavailable value is replaced with zero.

The existing 878 marker rows in every dataset and genotype scope are preserved exactly, including numeric bits, text and missingness. The new eleven rows per atomic cluster are appended. QC, exclusions, cluster memberships, embeddings, annotations, DESeq2 results and full-gene mouse boxplot data remain unchanged. This is new descriptive marker aggregation, with **no new model fitting**.

`native/` contains the new eleven-gene summaries only:

- `extension_gene_cluster_mouse_summary.tsv.gz`: exact atom × mouse coverage, counts, detection and log-expression sums, including zero-nucleus mouse/atom rows.
- `extension_gene_cluster_summary.tsv`: ALL summaries, support, rest-of-population contrasts and declared reference scaling.
- `extension_gene_cluster_summary_WT.tsv`: WT-only summaries with independently calculated WT reference scaling.

`validation/native_aggregate_comparison.json` verifies the new producer against three already-produced marker genes (`Ccnd2`, `Pvalb`, `Sst`) in every branch × atomic cluster × mouse. Counts and detection are exact; the largest floating log-sum difference is below `6 × 10⁻¹³`.

`validation/exact_native_count_comparison.json` is an additional independent read-only check of **all 64,548 new gene × atom × mouse integer sums** against every one of the 489 pre-existing full-gene native pseudobulk arrays. Every value agrees exactly. This supplemental receipt was generated after the extension's immutable completion receipt and binds to that receipt's hash.

`validation/existing_marker_preservation.json` records preservation hashes for all original marker rows. `COMPLETE.json` records source object/member hashes, output hashes, executed code hashes and scope. Its output paths are relative to the original extension build directory; hosted payloads live in `site/data/shards/` after integration.

Reproduction requires the existing processed full-gene count objects, frozen native prepared memberships and native results from the recorded completed run. The two source scripts are [`extend_screenshot_markers.py`](../../../scripts/extend_screenshot_markers.py) and [`validate_marker_extension_counts.py`](../../../scripts/validate_marker_extension_counts.py). Both require explicit input/output paths and refuse to overwrite their outputs. No analysis is automatically launched by opening this repository.
