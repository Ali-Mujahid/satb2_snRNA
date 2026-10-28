Completed hdWGCNA enrichment and atlas packaging while preserving all four previously finished network fits.

This analysis examined coexpression programs in the frozen Seurat CA1 deep-associated and superficial-associated populations and inhibitory clusters 1, 6 and 8. It used 23,143 selected nuclei from six WT and six cKO mice, with each mouse as the unit for genotype inference. The primary networks contain 10 CA1 and 5 inhibitory modules, with separately fitted consensus sensitivities containing 8 and 4 modules. The earlier controller stopped during SVG export after the network fits had finished, so this recovery verified and copied the completed preparation and network files into a fresh directory and ran only unchanged-specification enrichment and packaging. The recovery completed with exit code 0 and a successful FINAL_VALIDATION.json. These outputs support exploratory coexpression and gene-set associations, while technical partition confounding and several module-quality limitations restrict mechanistic interpretation.

## 1. Analysis context

Registry 54ace3f1-80e2-439c-b810-fae72950596d; native run snakemake-satb2-hdwgcna-20261026-20261026-110747-aeebc0a6. Existing user-selected Snakemake 9.27 controller, immutable recovery workflow 1.0.2, hdWGCNA 0.4.12 and WGCNA 1.74. Frozen endpoint counts: CA1_Deep 8,750; CA1_Superficial 11,706; IN_C1 1,718; IN_C6 501; IN_C8 468. Raw RNA, memberships and original DESeq2 gene fits were reused unchanged. Full genotype-blind eligible RNA is required for coexpression; the separate communication analysis retains its strict both-endpoint DEG-only rule.

## 2. Lifecycle and verification

Durable native state: completed, return code 0. All four copy/enrichment/package/final tasks completed. All 250 original prepared/network files, totaling 812,883,283 bytes, match the frozen source inventory and recovered SHA-256 receipts; each network checkpoint MD5 seal and the original fitting R source are verified. Scientific configuration is unchanged except new output root, workflow version and explicit copy provenance. All six enrichment TSVs are byte-identical to those written before the previous export failure. SVG and PDF now use explicit working graphics devices. Scientific result receipts exclude mutable controller metadata.

## 3. Findings

Primary CA1: 20,456 nuclei, 1,200 metacells, 10,618 tested genes, 4,898 assigned genes and 5,720 grey genes. Primary inhibitory: 2,687 nuclei, 359 metacells, 11,065 tested genes, 2,304 assigned genes and 8,761 grey genes. All module–endpoint scores contain twelve independent mice. The 13-nucleus IN_C8/cKO933 stratum contributes to scoring but falls below the 20-cell metacell-discovery minimum. Primary global BH yields 12/35 associations (8 CA1, 4 inhibitory); the separately corrected same-partition sensitivity yields 5/35, all CA1. The latter are CA1 black higher in cKO in both programs, purple lower in both, and brown higher in superficial-associated nuclei. Consensus yields 5/28 primary associations and 0/28 in the same-partition sensitivity. Enrichment contains 993/72,805 GO BP tests and 40/105 original-DEG overlap tests below FDR 0.05; no original-DEG module overlap passes in C6 or C8.

## 4. Interpretation

CA1 black and purple are the strongest within-dataset genotype-associated module leads because their direction and FDR support persist in the matched-partition subset and agree with directional original-DEG enrichment. This is convergence within reused data, not independent replication. Module eigengenes are continuous native PC1 scores oriented to average RNA expression; their differences are not gene log2 fold changes. GO enrichment describes a gene set and does not establish pathway activation. Consensus modules must be matched by genes, not color names, and are neither independent replication nor formal preservation tests. Network edges are saved undirected RNA coexpression, not TF binding.

## 5. Artifacts inspected

Run root: /Users/mujahid357/Documents/codex_new/SATB2_REGULATORY_20261026/runs/20261026_hdwgcna_v1_0_2_export_recovery. Inspected results/FINAL_VALIDATION.json, REUSED_NATIVE_RECEIPT.json, prepared membership/coverage audits, four completed network checkpoints and adjusted association tables, enrichment TSVs, atlas transport and paired native GO PNG/PDF. Native figure labels and axes were checked visually; detailed read-only audits retain exact quality diagnostics. Publication staging has no scientific fitting and remains a separate operation.

## 6. Limits

Genotype partly confounds sequencing partition; the prespecified same-partition subset is six WT and three cKO mice. IN_C8 remains mixed in molecular identity, C1 dominates inhibitory metacell discovery, and 79.2% of tested inhibitory genes remain grey in the pooled fit. Primary CA1 green and inhibitory yellow are mitochondrial-heavy and should not be promoted as SATB2 mechanisms without QC review. Sparse nuclear RNA causes material use of the configured individual Pearson fallback in signed-bicor connectivity; median nuclear kME is modest. The original R session reported 50 or more warnings without retaining their full list, so their exact causes and count are unknown. Original Satb2 gene RNA is higher in cKO in both CA1 programs; these data do not verify functional SATB2 protein deletion, targeted exon loss or direct regulation. No motif-supported regulons or nucleus-level TF activity UMAPs were produced.

## 7. Next action

Publish verified saved tables, figures, native state/source archive and the interactive explorer additively to the existing GitHub atlas. Keep primary and sensitivity FDR families distinct, expose exact scope and missing values, and explain that display changes never refit the saved results. Retain all previous analyses and the export-failure records. Matching mouse motif references remain unavailable; no substitute database or fabricated regulon is used.
