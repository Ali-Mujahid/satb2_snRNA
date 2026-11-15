# C6 and SST cluster comparison: data provenance

This descriptive display uses frozen Seurat original-full inhibitory cluster memberships (all available 6 WT + 6 cKO mice), not the 6 WT / 3 cKO sensitivity analysis. `data.json` contains per-cluster and per-mouse nuclei and saved MapMyCells reference assignment counts, and pooled WT+cKO marker summaries. The supertype predictions are not physical anatomy, tracing, or physiological classifications.

Marker dots use detected fraction (`n_detected/n_nuclei`) and mean log1p-normalized original UMI counts per 10,000 full-library UMIs. `scaled_expression_z` uses the pre-existing z definition: equal-weight all original inhibitory cluster means, population standard deviation (ddof 0). It is not recomputed over only the four selected clusters. Detection is not expression intensity and a zero does not prove biological absence.

For Elfn1, Chodl, Chrna2, Ctnnd2 and Htr3a, which were missing in the published saved summary, the original-full 8,911-nucleus h5ad `layers/counts` matrix was matched by `nucleus_id` to frozen Seurat labels. The shared original count layer is used solely for quantification; Scanpy clusters are not used. Counts per nucleus were normalized using its original full library, not the 32,285-feature subtotal. The reconstructed quantities for Cdh7, Sst, Reln, Trpc6, Pcdh11x, and Grin3a agreed with saved Seurat summaries across all four groups (means and z within 1e-6; fractions within 1e-9).

Rebuild locally, in this order, from repository root:

```
/Users/mujahid357/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3 analysis/inhibitory_c6_review_20261116/build_data.py
/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/.venv_scanpy/bin/python analysis/inhibitory_c6_review_20261116/extract_unavailable_markers.py
```

The inputs are absolute paths recorded in `data.json.provenance`; the 20261013 source files and 20261014 frozen Seurat outputs are preserved. No DE, cluster memberships, normalization model, label or UMAP coordinates are changed.
