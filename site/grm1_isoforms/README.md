# Grm1 alpha/beta: twelve-mouse RNA evidence

Open `index.html` in Chrome or Safari. All plot data and plotting code are embedded and work offline; static schematics and downloads require the extracted `figures/` folder alongside the HTML.

The analysis uses the original Seurat populations: superficial CA1 C2+C3, deep CA1 C1+C5, all inhibitory clusters, and especially inhibitory C2 and the original parent C6. C8 is excluded from superficial CA1. All six WT and six cKO mice are included. No previous clustering, labels, counts or DESeq2 models were changed.

## Start with the right names

- Grm1-201 encodes the long mGluR1 alpha protein (1199 amino acids).
- Grm1-202 and Grm1-203 encode the same shorter mGluR1 beta protein (906 amino acids) with different untranslated-region models.
- Grm1-204/-205/-206 are noncoding RNA models retained as competitors for shared reads.
- Grm7 encodes the separate mGlu7A/B receptors. It does not encode mGluR1 alpha/beta.
- Figure 4B in Ferraguti et al.1998 is a hippocampal synaptosome protein blot; “4B” is a panel label, not an isoform. It does not identify the contributing subregions or cell types.

## Choose the appropriate evidence

`alpha_beta_evidence_per_mouse.tsv`: diagnostic RNA features, conservative alpha/beta assignments and compatible upper bounds. Conservative assignments require family-exclusive model compatibility, diagnostic features and compatible splicing. The splice may be constitutive and need not be the alpha/beta-discriminating terminal junction. These are evidence lower bounds, not full-length transcript counts.

`alpha_beta_specific_junction_counts_per_mouse.tsv`: direct counts of the three alpha/beta-discriminating terminal junctions. These counts are much sparser and are the appropriate evidence for local terminal-splicing events. Zero support does not prove absence. Adjacent junction counts can share molecules and must not be summed into transcript abundance.

`isoform_evidence_per_mouse.tsv`: all six individual transcript models. Beta is often shared between 202/203; absent unique model assignments do not imply absent beta-family RNA. Exon-only signal can be pre-mRNA.

`gene_UMI_counts_per_mouse.tsv`, exon/junction/CDS/UTR tables and `population_genotype_master_summary.tsv` retain source values and denominators. Normalization uses original full-library UMIs, not the rendered HVG-only nCount_RNA. Mouse is the independent unit. Raw totals are not genotype fold changes or new DESeq2 statistics.

## Plots and controls

The explorer has 35 populations and 12 evidence views. Choose population and raw, per-nucleus, full-library or gene-relative values. Open the appearance controls to change colors, label size, column width, heatmap/dotplot style, dot size, displayed feature filter and printed values. Downloads include SVG, PNG, PDF, TSV and plot-state/data JSON.

`figures/` contains 28 static figures, each in SVG, PNG and PDF. The transcript and protein schematics explain the verified alpha/beta mapping. A protein schematic is not a protein measurement.

## Validation and reproducibility

`VALIDATION.json` records count checks; `EXPLORER_QA.json` records browser and export checks. `ANALYSIS_REPORT.md` explains the results, sparse terminal junctions and limits. `COMPLETE.json`, `annotation.json`, the Grm1-only reference GTF, transcript/feature tables and source hashes record the native run. `source/` contains the corrected counting workflow and reviewed export code; raw BAMs and individual nucleus barcodes are not distributed in this package.

The corrected native run completed successfully after an initial reference-only preflight attempt identified that noncoding models also had to be retained. That failed preflight produced no biological results. No full-length transcript reconstruction, differential isoform-switching test, protein abundance assay or new functional inference was performed.
