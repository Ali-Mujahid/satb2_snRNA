# Exact per-nucleus expression display vectors

These files provide measured source counts for all **889 atlas marker genes** in
the three existing count branches. They support single-gene UMAP overlays for
both Scanpy and Seurat without estimating expression from a cluster mean.

The source objects are the unchanged full-gene integer-count AnnData files from
the completed native analysis, already distributed in the `atlas-v2.1.0` release.
Their original SHA-256 hashes and release asset names are recorded in
`PROVENANCE.json`. The original source cohorts contain 108,062 nuclei for
`original_full` and 104,298 each for `paired_original` and `paired_cellbender`.
Rendering uses the atlas's existing coordinate subset and display exclusions;
the broader source axis does not restore excluded nuclei to a plot.

Each branch has one `axis__BRANCH.json.gz` containing ordered unique nucleus IDs
and their complete source-library counts. Every source library was recomputed
from all 32,285 genes and agrees exactly with the native verified library field.

Each `gene__BRANCH__GENE_ID.json.gz` contains sorted source-axis indices with
strictly positive integer gene counts. A nucleus present in the source axis but
absent from that sparse vector has a measured count of zero. A nucleus absent
from the source axis has **unavailable expression**, not zero. The payload's
`axis_payload_sha256` binds it to the correct branch and axis order.

Normalized expression is evaluated per nucleus as:

```
log1p(10000 * exact_source_gene_count / complete_source_library_count)
```

This uses full-library denominators, including in the corrected CellBender
branch. It does not use a marker-only library, a cluster mean, imputation, an
integrated assay, an SCT residual or a new statistical model. Exact source counts
are also retained for a raw-count display.

The small catalog in `site/data/GENE_EXPRESSION_CATALOG.json` maps branch and
stable gene ID to these files. It records both compressed-byte and decompressed
payload SHA-256 hashes and byte sizes. Its raw GitHub URLs are pinned to the
immutable `atlas-v3.1.0` tag. Keeping the vectors outside `site/` avoids inflating
the GitHub Pages deployment; only requested genes and source axes are fetched.

Validation includes every coordinate ID in all 24 unchanged datasets; all 889
per-gene source-count sums and detected-nucleus totals per branch; independently
read, reproducible 257 × 61 source cell/gene spot checks per branch; and exact
counts, detection and full-library normalized means for `Ccnd2`, `Pvalb` and `Sst`
in every displayed atomic cluster against the previously packaged native marker
summaries. Existing QC, coordinates, annotations, masks, cluster memberships,
DESeq2 fits and marker-summary payloads are not modified by this export.

Producer: `scripts/build_nucleus_gene_vectors.py`. Focused representation and
integrity tests: `tests/test_nucleus_gene_vectors.py`. The producer refuses to
overwrite an existing output directory or catalog. It requires the existing
scientific Python environment, explicit native object manifest and output path.
