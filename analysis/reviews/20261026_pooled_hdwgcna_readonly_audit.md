# Supplemental read-only audit of the completed pooled hdWGCNA checkpoints

Reviewed 2026-10-26; final checkpoint rehash at **10:40:12 UTC**. This is a supplemental inspection of saved **pooled primary** outputs. It is **not** a terminal analysis/registry report, a consensus-network audit, a new fit or a final genotype-significance report. The parent workflow will assess global association BH FDR and enrichment after those outputs complete.

No registered scientific source, count matrix, annotation, prior result or hdWGCNA output was modified. Computation in this audit consisted of reading tables, checking identities/checksums, counting existing rows, summarizing saved statistics and viewing existing figures. No network, linear model, DESeq2 analysis or preservation test was launched.

## Exact checkpoints and integrity

Run root: `/Users/mujahid357/Documents/codex_new/SATB2_REGULATORY_20261026/runs/20261026_hdwgcna_v1_0_1`

| Checkpoint | Native completed time | Declared files rehashed | MD5 mismatches | completed.json SHA256 |
| --- | --- | ---: | ---: | --- |
| CA1 pooled |2026-10-26 10:32:08 UTC|56|0|`c45524ad4837ddb42239a82ac77b67ac24983bc4b3fd39ba3e8166780546acc2`|
| IN pooled |2026-10-26 10:36:57 UTC|56|0|`42f050d624a98a4ec95d901b36fabdedda42dceb0bacf430b8196858b8867dd4`|

Native manifests:

- [CA1 completed.json](/Users/mujahid357/Documents/codex_new/SATB2_REGULATORY_20261026/runs/20261026_hdwgcna_v1_0_1/results/hdwgcna/CA1/pooled/completed.json), manifest MD5 `e6c13aff2e54d062fa3498a323c3358c`.
- [IN completed.json](/Users/mujahid357/Documents/codex_new/SATB2_REGULATORY_20261026/runs/20261026_hdwgcna_v1_0_1/results/hdwgcna/IN/pooled/completed.json), manifest MD5 `9770fd5b50648f839c5e2180a30016bf`.

All files listed by each native manifest existed and matched its stored MD5. This demonstrates file consistency with the completed checkpoint; it does not by itself validate a biological interpretation.

## Coverage and score checks

| Quantity | CA1 pooled | IN pooled |
| --- | ---: | ---: |
| Unique tested gene symbols / assignment rows |10,618|11,065|
| Assigned genes |4,898|2,304|
| Non-grey modules |10|5|
| Grey/unassigned genes |5,720 (53.87%)|8,761 (79.18%)|
| Original nuclei with finite scores |20,456|2,687|
| Distinct independent mice |12|12|
| Original endpoint counts |DL8,750; SL11,706|C1 1,718; C6 501; C8 468|
| Mouse × endpoint × module summary rows |240|180|
| Duplicated summary keys |0|0|
| Metacells |1,200|359|
| Metacell member rows |30,000|3,590|
| Nuclei per metacell |25 throughout|10 throughout|
| Mixed mouse/endpoint/genotype metacells |0|0|
| Nonfinite nuclear eigengene values |0|0|
| PC1 correlation with mean module RNA |0.97448–0.99281|0.97024–0.99792|

CA1 has50 metacells in every one of24 mouse × program strata. All strata are included in discovery. IN has35 included discovery strata and one explicitly excluded stratum: **cKO_933 × IN_C8,13 nuclei, below min_cells20**. All13 of these nuclei are present in the saved nuclear eigengene score table. Discovery metacells are C1=259, C6=50 and C8=50; therefore C1 contributes about72% of IN discovery metacells. Sparse C6/C8 strata constrain how comprehensively their covariance is learned even though all mice remain in downstream score summaries.

The high grey fractions are explicit unassigned **genes**, not missing nuclei or unresolved cell clusters. They must not be hidden or forced into modules. Every assigned gene has a finite nonnegative saved nuclear kME; grey kMEs are intentionally NA. However, median assigned nuclear kME is only0.095–0.185 across CA1 modules and0.133–0.213 across IN modules. Thus strong hub behavior should not be generalized to equally strong membership evidence for every assigned gene.

Source tables:

- [CA1 assignments](/Users/mujahid357/Documents/codex_new/SATB2_REGULATORY_20261026/runs/20261026_hdwgcna_v1_0_1/results/hdwgcna/CA1/pooled/module_membership.tsv), [CA1 discovery coverage](/Users/mujahid357/Documents/codex_new/SATB2_REGULATORY_20261026/runs/20261026_hdwgcna_v1_0_1/results/hdwgcna/CA1/pooled/metacell_coverage.tsv), [CA1 membership](/Users/mujahid357/Documents/codex_new/SATB2_REGULATORY_20261026/runs/20261026_hdwgcna_v1_0_1/results/hdwgcna/CA1/pooled/metacell_membership.tsv.gz), [CA1 nuclear scores](/Users/mujahid357/Documents/codex_new/SATB2_REGULATORY_20261026/runs/20261026_hdwgcna_v1_0_1/results/hdwgcna/CA1/pooled/nucleus_eigengenes.tsv.gz), [CA1 orientation](/Users/mujahid357/Documents/codex_new/SATB2_REGULATORY_20261026/runs/20261026_hdwgcna_v1_0_1/results/hdwgcna/CA1/pooled/eigengene_orientation.tsv).
- [IN assignments](/Users/mujahid357/Documents/codex_new/SATB2_REGULATORY_20261026/runs/20261026_hdwgcna_v1_0_1/results/hdwgcna/IN/pooled/module_membership.tsv), [IN discovery coverage](/Users/mujahid357/Documents/codex_new/SATB2_REGULATORY_20261026/runs/20261026_hdwgcna_v1_0_1/results/hdwgcna/IN/pooled/metacell_coverage.tsv), [IN membership](/Users/mujahid357/Documents/codex_new/SATB2_REGULATORY_20261026/runs/20261026_hdwgcna_v1_0_1/results/hdwgcna/IN/pooled/metacell_membership.tsv.gz), [IN nuclear scores](/Users/mujahid357/Documents/codex_new/SATB2_REGULATORY_20261026/runs/20261026_hdwgcna_v1_0_1/results/hdwgcna/IN/pooled/nucleus_eigengenes.tsv.gz), [IN orientation](/Users/mujahid357/Documents/codex_new/SATB2_REGULATORY_20261026/runs/20261026_hdwgcna_v1_0_1/results/hdwgcna/IN/pooled/eigengene_orientation.tsv).

## Module sizes and hub evidence

Color names are arbitrary native labels. **CA1 blue and IN blue are separate modules and are not homologous simply because their colors match.** Report a lineage prefix with each module name.

| CA1 module | Genes | Median nuclear kME | Selected saved top hubs |
| --- | ---: | ---: | --- |
|turquoise|1,667|0.115|Lingo2, Ak5, Pcdh7, Cadm1, Fgf14, Ntng1|
|blue|1,117|0.095|Opcml, Ano3, Cntnap5a, Ryr2, Rbfox1, Epha3|
|brown|575|0.102|Pde10a, R3hdm1, Ntrk2, Ndst4, Sorcs3, Homer1|
|yellow|512|0.130|Nrxn3, Kcnh7, Chrm3, Epha6, Cadm2, Kcnd2|
|green|333|0.118|mt-Co3, mt-Co1, Gm42418, mt-Nd1, mt-Atp6, mt-Co2|
|red|170|0.165|Epha7, Ppp3ca, Htr1f, Frmpd4, Epha5, Kcnq3|
|black|152|0.185|Fgf12, Calb1, Slc35f4, Adamts3, Ncam2, Ntm, Dcc|
|pink|132|0.179|Lrrtm4, Neto1, Cnih3, Sorbs1, Nrg1, Htr2c|
|magenta|131|0.147|Snhg11, Atp6v0b, Slc17a7, Ahi1, Hnrnpa2b1, Hspa5|
|purple|109|0.163|Astn2, Cntnap5c, Gabrg3, Khdrbs3, Kcnip4, Ptn|

| IN module | Genes | Median nuclear kME | Selected saved top hubs |
| --- | ---: | ---: | --- |
|turquoise|979|0.198|Ntrk3, Dock4, Lrrc4c, Adgrb3, Gabrg3, Tmem108|
|blue|900|0.133|Gm20754, Grin3a, Dcc, Sema6d, Grm7, Nrxn1, Reln|
|brown|210|0.213|Nrg3, Dpp10, Rbfox1, Gria3, Cntnap5a, Lsamp|
|yellow|122|0.186|mt-Co3, mt-Co1, mt-Atp6, mt-Nd1, mt-Co2, mt-Nd4|
|green|93|0.147|Necab1, Gria4, Cntnap5c, Kirrel3, Il1rapl2, Stxbp6|

These are RNA co-expression hub rankings from [CA1 hub genes](/Users/mujahid357/Documents/codex_new/SATB2_REGULATORY_20261026/runs/20261026_hdwgcna_v1_0_1/results/hdwgcna/CA1/pooled/hub_genes.tsv) and [IN hub genes](/Users/mujahid357/Documents/codex_new/SATB2_REGULATORY_20261026/runs/20261026_hdwgcna_v1_0_1/results/hdwgcna/IN/pooled/hub_genes.tsv). They do not identify verified TF-target edges or direct ligand–receptor binding. Exported network edges are a bounded top-hub bicor subset, not a complete saved TOM.

## Stability and descriptive direction checks

All120 CA1 and60 IN fixed-module leave-one-mouse-out stability rows are finite and have unique module × omitted-mouse keys. Minimum full-versus-omitted connectivity correlation is0.95091 in CA1 and0.98276 in IN. Minimum top10 hub retention is80% in each network. These summaries hold modules and nuclear PC1 scores fixed; they are **not de novo leave-one-mouse-out discovery, independent validation or formal WGCNA modulePreservation**. Shared original nuclei within metacells also mean metacells are not genotype replicates.

The saved genotype bootstrap tables contain20 CA1 contrasts and15 IN contrasts, each with200 independent-mouse resamples and correctly ordered interval endpoints. Primary score summaries use6WT/6cKO. The same-partition summaries contain6WT/3cKO throughout.

Saved effect directions were compared descriptively, without selecting on unadjusted P-values or calculating a substitute FDR:

- **CA1:**18/20 contrasts retain direction in all12 omitted-mouse fits. Blue/SL retains it11/12 and reverses in the same-partition subset. Green/DL retains it7/12 and has a near-zero full-data effect (−0.018 raw PC1 units). These two contrasts are fragile.
- **IN:** C8 turquoise retains direction6/12 and C8 green8/12, with near-zero full-data effects (+0.0032 and+0.0259 raw PC1 units). C8 brown retains it11/12 and reverses in the same-partition subset. C8 yellow retains it12/12 but also reverses in that subset. C6 turquoise and yellow retain direction11/12. These require extra caution; stable full-sample sign alone does not resolve partition confounding.

No branch `effects.tsv` or `*_unadjusted.tsv` was interpreted as final significance. Raw PC1 units differ with the module; absolute coefficients across different modules are not a direct ranking of biological effect strength. Final inference must use the parent workflow's global BH family across tested module × endpoint associations, along with design limitations and enrichment evidence.

Stability source examples: [CA1 fixed-module stability](/Users/mujahid357/Documents/codex_new/SATB2_REGULATORY_20261026/runs/20261026_hdwgcna_v1_0_1/results/hdwgcna/CA1/pooled/fixed_module_connectivity_stability.tsv), [IN fixed-module stability](/Users/mujahid357/Documents/codex_new/SATB2_REGULATORY_20261026/runs/20261026_hdwgcna_v1_0_1/results/hdwgcna/IN/pooled/fixed_module_connectivity_stability.tsv).

## Candidate interpretation and cautions

CA1 black is a **Calb1/adhesion-associated candidate program**, while brown is an **Ndst4-associated candidate program**. Saved equal-mouse mean scores place black higher in SL than DL within each genotype and brown higher in DL. This is internal consistency with pre-existing marker-defined programs; it is not independent validation or grounds to relabel nuclei anatomically. Spatial work shows that CA1 layer markers change their relative anatomical position along the longitudinal axis, including Ndst4-associated layers. [Pachicano et al., Nature Communications2025](https://www.nature.com/articles/s41467-025-66613-y)

CA1 yellow's Nrxn3/Cadm2/Chrm3/Kcnd2 and pink's Lrrtm4/Neto1/Cnih3/Nrg1 nominate broad synapse/excitability hypotheses. IN blue's Grin3a/Grm7/Nrxn1/Reln/Dcc/Sema6d pattern is mainly C6-associated; IN turquoise's Ntrk3/Lrrc4c/Gabrg3/Adgrb3 pattern is mainly C1-associated in the saved heatmap. Thus subclass/endpoint identity can drive covariance. These provisional descriptions come from saved hub content; **no formal pathway-enrichment claim is made here**. Wait for the completed native enrichment tables before attaching specific significant GO terms.

Mitochondrial content deserves particular caution. CA1 green contains13 mitochondrial genes and9 of its top10 hubs are mt genes. IN yellow contains12 mitochondrial genes and9 of its top10 hubs are mt genes. In nuclei data, such programs can reflect residual cytoplasm, contamination or sample quality; they must not be presented as a SATB2 regulatory mechanism solely because they form a stable module. Official snRNA-seq guidance discusses mitochondrial signal as a stripping/QC indicator, and primary nuclei-QC work documents technical noise requiring explicit assessment. [Bioconductor OSCA snRNA-seq guidance](https://bioconductor.org/books/3.16/OSCA.advanced/single-nuclei-rna-seq-processing.html), [primary nuclei QC study](https://pmc.ncbi.nlm.nih.gov/articles/PMC5338304/)

Neither all-mouse linear models nor6WT/3cKO sensitivity fully separate genotype from the partly confounded technical design. C8's small strata and pre-existing mixed-identity caveats persist. These modules are co-expression programs, not SCENIC motif-supported regulons or independent corroboration of communication. Reusing the same genes in DEG, LR and module views does not create independent mechanistic evidence.

## Actual visual inspection

The following native PNGs were opened and inspected. Titles, axes, labels, legends and data marks were visible without clipping or broken characters:

- [CA1 mouse eigengene heatmap PNG](/Users/mujahid357/Documents/codex_new/SATB2_REGULATORY_20261026/runs/20261026_hdwgcna_v1_0_1/results/hdwgcna/CA1/pooled/figures/mouse_eigengene_heatmap.png)
- [CA1 module RNA dotplot PNG](/Users/mujahid357/Documents/codex_new/SATB2_REGULATORY_20261026/runs/20261026_hdwgcna_v1_0_1/results/hdwgcna/CA1/pooled/figures/module_RNA_dotplot.png)
- [CA1 top hubs PNG](/Users/mujahid357/Documents/codex_new/SATB2_REGULATORY_20261026/runs/20261026_hdwgcna_v1_0_1/results/hdwgcna/CA1/pooled/figures/top_hubs.png)
- [CA1 soft-power diagnostics PNG](/Users/mujahid357/Documents/codex_new/SATB2_REGULATORY_20261026/runs/20261026_hdwgcna_v1_0_1/results/hdwgcna/CA1/pooled/figures/soft_power_diagnostics.png)
- [CA1 block1 dendrogram PNG](/Users/mujahid357/Documents/codex_new/SATB2_REGULATORY_20261026/runs/20261026_hdwgcna_v1_0_1/results/hdwgcna/CA1/pooled/figures/gene_dendrogram_block_1.png)
- [IN top hubs PNG](/Users/mujahid357/Documents/codex_new/SATB2_REGULATORY_20261026/runs/20261026_hdwgcna_v1_0_1/results/hdwgcna/IN/pooled/figures/top_hubs.png)
- [IN mouse eigengene heatmap PNG](/Users/mujahid357/Documents/codex_new/SATB2_REGULATORY_20261026/runs/20261026_hdwgcna_v1_0_1/results/hdwgcna/IN/pooled/figures/mouse_eigengene_heatmap.png)
- [IN metacell coverage PNG](/Users/mujahid357/Documents/codex_new/SATB2_REGULATORY_20261026/runs/20261026_hdwgcna_v1_0_1/results/hdwgcna/IN/pooled/figures/metacell_coverage.png)

The actual [CA1 top-hubs PDF](/Users/mujahid357/Documents/codex_new/SATB2_REGULATORY_20261026/runs/20261026_hdwgcna_v1_0_1/results/hdwgcna/CA1/pooled/figures/top_hubs.pdf) was also checked with pdfinfo (one page,720×864pt, Cairo1.17.6, PDF1.7) and rasterized at90dpi using pdftoppm. The rendered page was opened and the gene labels, numeric axes, module panels and subtitle were intact. The inspection-only temporary raster is `/tmp/satb2_ca1_20261026_top_hubs_pdf_review.png`; it is not a new scientific result or a modification of the native PDF.

The heatmap clearly labels its values as display z-scores. The RNA dotplot legends correctly identify mean expression for both size and color; they are not percent-expressing-cell legends. Native module colors include bright yellow, green and magenta, which are visually harsh against white; this is a presentation limitation, not grounds to alter the frozen analysis. Interactive display palettes can be changed later while retaining the original scientific values and lineage-qualified labels.

## Method-source context

The basic and consensus hdWGCNA workflows distinguish gene/module inference from downstream eigengene summaries. The pinned connectivity implementation's sparse branch uses corSparse; this run's explicit sparse=FALSE setting was chosen for the stated bicor kME. Native raw module PC1 is oriented toward positive mean RNA; no genotype/sample Harmony was applied to these outcome scores. [Official basic workflow](https://smorabit.github.io/hdWGCNA/articles/basic_tutorial.html), [pinned ModuleConnectivity implementation](https://github.com/smorabit/hdWGCNA/blob/e3344d1f7bbac4264adf94f7aa31e0802fa8282d/R/ModuleConnectivity.R), [pinned ModuleEigengenes implementation](https://github.com/smorabit/hdWGCNA/blob/e3344d1f7bbac4264adf94f7aa31e0802fa8282d/R/ModuleEigengenes.R)

This audit leaves consensus evaluation, final global association FDR, enrichment, terminal registry reporting and the subsequent native TF-activity analysis to their respective workflow stages.
