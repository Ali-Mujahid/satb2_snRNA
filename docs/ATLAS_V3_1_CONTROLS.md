# Atlas v3.1 plot controls

Open the [web atlas](https://ali-mujahid.github.io/satb2_snRNA/), select a population and count branch, then use the controls below. Export settings to transfer the same view to another browser. Earlier portable archives remain preserved; this guide describes the v3.1 web viewer.

## UMAP colors and single-gene expression

**Category palette** changes categorical labels, clusters, mice or genotypes. Select **Category to recolor**, choose a color, then **Apply color** to override that category in the current coloring field. **Reset color** removes that override. These colors are independent of the **Expression / z-score continuous palette**, whose low/mid/high colors are used for numeric expression and dotplots.

Choose a gene under **Find an overlay gene**. With **Show gene when selected** enabled, selection opens its expression UMAP automatically; the **Show selected gene on UMAP** button does the same explicitly. Both Scanpy and Seurat use the same selected gene and count branch.

- **Individual nuclei · source counts** shows each nucleus's exact existing branch count or `log1p(10000 × gene count / complete source library)`. It covers all **889** panel genes on all three branches. Zero-expression nuclei use the adjustable zero color. Select a genotype under **Nuclei shown on UMAP** to display WT or cKO; the cluster-summary scope selector does not apply in individual mode.
- **Saved cluster summaries** shows existing cluster means, detected fractions or produced gene-scaled z values. Every nucleus in a cluster receives that summary color; this mode does not show individual RNA. WT/cKO summary selection applies in this mode. cKO produced z remains unavailable.

**Gene color limits** supports automatic or manual minimum/maximum values. **Same color scale in both tools** aligns numeric colors. Limits clip display colors only. Points still use existing UMAP coordinates, existing tool-specific memberships and current exclusions.

The exporter reads the three unchanged native full-gene count objects and joins exact stable nucleus IDs. Full-library denominators include all source genes, not only the 889-gene panel. A present source nucleus with no sparse entry is zero; an absent source ID is unavailable. The [export receipt](../site/data/GENE_EXPRESSION_PROVENANCE.json) records source hashes, full-library and ID checks. No new clustering or statistical model is fitted.

## Short labels

**Full original labels** retains source wording. **Compact** uses names such as `Sup.assoc.exp.` and `Deep.assoc.exp.`. **Short** uses `SL-like` and `DL-like`. These mean superficial- and deep-associated **expression**; physical anatomical layers remain undetermined. Mixed, weak, conflict, review and provisional qualifiers remain visible.

For your own wording, select an **Original label to shorten**, enter **Custom display text**, and apply it. Display aliases affect UMAP centroids/legends, dotplot rows, selectors and PDF display text. They preserve original labels, grouping keys, nucleus membership and saved DE targets. They do not create custom DE unions. The separate **Labels & exclusions** panel retains its explicit annotation/union workflow.

For dotplots, **Dotplot row labels** chooses cluster IDs, display labels, or both. Cluster IDs remain independent between tools.

## Marker sequence

**Screenshot · all 75 genes** restores the screenshot's exact left-to-right sequence and fifteen blocks of five, on one 90-gene page. The screenshot cluster numbers are not identities assigned to this dataset.

**Shared gene-column order** supports preset/selected order, grouping by cluster expression pattern, or a focused cluster's expression contrast. Select Scanpy or Seurat as the reference to give both plots the same gene-column sequence. Automatic ordering describes mean-expression contrasts against nucleus-weighted other retained clusters; it is not a new differential test or validated representative-marker selection.

The open manual section can move a selected gene, apply a pasted full/partial order, reverse the order, or sort symbols A–Z. Freeze an automatic order before editing if desired. Cluster rows have separate ascending, descending or custom order for each tool. PDF and selected-value exports retain the chosen gene sequence.

## Volcano plots

Axes are **log2 fold change (cKO / WT)** and **−log10 adjusted P**, with numeric ticks. Positive effects mean higher expression in cKO. **Maximum adjusted P** and **Minimum absolute log2 fold change** define which saved rows qualify as DEGs for highlighting and labels.

**Label top DEGs** defaults to up to **20** eligible genes per plot. Eligible genes are ranked by smallest adjusted P, then greatest absolute log2 fold change, then stable ID. Increase/decrease the maximum or turn labels off. If fewer genes pass, only those genes are eligible; missing/invalid adjusted P or fold change is not replaced. Labels may be omitted if the plot lacks space, with the number reported below the plot. Extremely small/zero adjusted P values have a documented display cap of 320 on the y axis; source values are unchanged.

Choose separate colors for **Higher in cKO**, **Lower in cKO** and **Below selected DEG thresholds**. Threshold guides can be toggled. **Volcano PDF** retains the displayed axes, gene labels, colors and measurement notes. These controls use the exact saved DESeq2 result and never refit it. Source seven-column DEG files remain unchanged.
