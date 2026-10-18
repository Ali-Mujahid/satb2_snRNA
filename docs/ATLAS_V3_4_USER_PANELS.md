# Your marker panels and visible cluster arranger

Open the [live atlas](https://ali-mujahid.github.io/satb2_snRNA/). Version 3.4 adds your three panels to **Marker expression**. Each is available in every population and selects the same ordered gene columns for Scanpy and Seurat. Each complete verified selection fits on one 120-gene plot page and can be downloaded using each plot's PDF button. Click its preset again to recover the supplied sequence after manual or expression-based reordering.

| Preset | Requested entries | Unique names | Available unique genes | Unavailable unique names |
|---|---:|---:|---:|---:|
| All-cell panel | 100 | 100 | 99 | 1 |
| CA1 / ProS panel | 158 | 152 | 118 | 34 |
| CA1 sublayer panel | 100 | 100 | 97 | 3 |

The original literal list, every position, repeated names and stable-ID mappings are available beside each preset and through **Download name mappings**. Unavailable names are shown explicitly and are never plotted as zero. Repeated names are selected once at their first verified occurrence. No spelling guesses were applied. `AC129186.1` is present in your original workbook but absent from the native gene annotation.

Four documented older names are disclosed as requested → stored source symbol: Nov → Ccn3, Fam19a1 → Tafa1, Fam19a2 → Tafa2, and Neph1 → Kirrel (current NCBI symbol Kirrel1). Their authoritative evidence links appear in the preset notes. Native source symbols remain on plotted axes. The audit retains your original requested names.

The registry now offers 1,055 genes, including 166 additive exports. Their marker summaries and individual-nucleus expression use the existing raw-count branches and full source library sizes. cKO-only z values remain unavailable as in the prior atlas; use mean expression for that scope. Descriptive marker grouping does not establish identity or create a new DE result.

## Arrange UMAPs

The controls are now directly **above the two UMAP plots**, with four large buttons:

- **Restore original UMAP** recovers the exact original coordinates and clears this population/branch's manual placement, rotation and reflection.
- **Bring clusters closer** packs cluster islands using their 95% core outlines. Every nucleus remains present, including peripheral nuclei that can overlap other islands. Internal geometry stays unchanged unless you manually resize a cluster; inter-cluster distances and neighbors are artificial.
- **Rounded cluster schematic** places the same nuclei in deterministic rounded discs. Disc area approximates nucleus count with a minimum readable radius. Both internal and external geometry are artificial. This is explicitly marked SCHEMATIC on the plot and PDF.
- **Hippocampus-like schematic** arranges clusters along a C-shaped guide in cluster-ID order. This is an illustration, without anatomical or spatial placement evidence.

Choose **Tool** and **Cluster**, or click a nucleus to select its cluster. The arrow buttons immediately move that cluster; **Move step** sets the distance as a percentage of the original source span. Arrows follow the screen direction after rotation and mirroring. Exact X/Y shift and scale controls are also visible. Advanced rotation, reflection and island gap remain in the left sidebar. One-click presets clear manual adjustments in this population/branch; switching the sidebar display mode retains them.

PDFs show the current display layout and its geometry warning. Display-coordinate CSVs include both original and adjusted coordinates, IDs and a nonspatial flag. Settings export/import preserves these display edits. The original UMAP, memberships, exclusions, biological labels and DESeq2 results remain intact.

## Validation and reproducibility

- Existing 889-gene catalogs and all original scientific shards remain byte-identical.
- 974,088 independent new gene × atomic cluster × mouse count sums match frozen native arrays exactly.
- 24 additive marker shards and 498 sparse nucleus-expression vectors are hash-verified and use the existing full-library source axes.
- `scripts/verify_v34_extensions.py` checks all new hashes, source bindings, ALL/WT/cKO sum conservation and produced scaling, plus literal panel order and unresolved symbols.
- Tests exercise all panels in all four underlying population scopes, both tools, all three count branches, actual layout controls, PDF paths and exact reset.

New data are pinned to immutable `atlas-v3.4.0` paths outside the Pages bundle. The previous portable v3.0 archives are preserved; this update targets the hosted GitHub atlas. No model fit, clustering, annotation, QC or differential analysis was rerun. The stopped inhibitory downsampling analysis remains stopped.
