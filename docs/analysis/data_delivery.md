# Data deliveries and processed native objects

The current atlas has 24 main datasets: 2 tools × 3 count branches × 4 populations. It has 878 marker genes; the produced mouse catalog covers 32,285 stable source genes and 1,060 endpoint descriptors. Statistical exports preserve exactly `gene,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj`, with separate stable-ID audits for duplicate/missing symbols and NA values.

The viewer uses existing produced results. Gene UMAP colors are cluster means, not per-nucleus RNA measurements; mouse values are produced aggregates with declared normalization/size factors. Exact results, model eligibility and all 12-mouse coverage remain in the saved-result delivery.

## Processed objects

The [JSON](../../analysis/provenance/PROCESSED_OBJECTS_RELEASE_MANIFEST.json) and [TSV](../../analysis/provenance/PROCESSED_OBJECTS_RELEASE_MANIFEST.tsv) manifests name 28 collision-free assets for `atlas-v2.1.0`, with exact streaming SHA256, byte size, source-relative path, branch/population and native nucleus count. URLs target declared assets; the publication task records completed transfers. Initial `PREPARED_FOR_UPLOAD` status alone does not establish network completion.

| Family | Files | Logical GiB |
|---|---:|---:|
| Scanpy `scanpy.h5ad` |12|8.879|
| Seurat `representation_HVG_only.rds` (12 log +1 original SCT)|13|10.327|
| Shared full-gene `branch_counts.h5ad` |3|1.817|
| Total |28|21.023|

Every individual object is below 2 GiB. Objects are byte-identical, with no alias insertion, fitting, format conversion or RNA recomputation. Publication hashed existing bytes without loading the objects.

Original native objects have 108,062 nuclei; paired branches 104,298. They precede revised All/EX shared filtering. Current original/paired retained cohorts are 107,747/103,983 All nuclei and 60,908/60,772 EX nuclei. IN/CA1 have no 100-nucleus cutoff. Join revision labels/masks by stable nucleus ID; do not assume native full objects equal displayed rosters.

Seurat HVG objects are representations, not full-gene count stores. Use same-branch shared full-gene integer AnnData and stable-ID gene metadata separately. Corrected/original paired masks and original-cluster identities remain fixed, separately from fresh paired cluster labels. Reviewed category/program maps are under `analysis/data/annotation_final/`.

## Included and excluded files

Included: original/executed workflows, actual configs/runtime records, source-bound sample/roster/reference metadata, exact workbooks, panel mappings, atomic/category memberships, completed producer receipts, endpoint indices and curated reviews. Large text metadata use documented lossless gzip. Panel coverage is 147 source/canonical rows, 145 mapped rows / 139 unique IDs/two unmapped rows; workbook-only coverage is 123 mapped rows / 117 IDs. This is separate from the full 878-gene panel union.

Explicitly excluded from Git and the28-object set:

- Raw Cell Ranger/CellBender and original MapMyCells query matrices. Source metadata/hashes are preserved; full upstream reproduction requires original study inputs.
- 13 redundant pre-graph checkpoints (8.424 GiB) and 26 separate neighbor-graph R objects (1.275 GiB). Completed technical control metrics remain included.
- Large native nucleus annotation/aggregate exports and historical recovery/build duplicates. Compatible labels, arrays, DE tables and figures are delivered through the viewer/result archives; exact stage reproduction can require additional source-bound files.
- Virtual environments, `.DS_Store`, `__pycache__`, controller caches, authentication files and unrelated data.

Original source locators are provenance. Relative public paths, explicit release URLs, hashes and semantic membership identify delivered files; a path matching the author's computer is not required for browsing.
