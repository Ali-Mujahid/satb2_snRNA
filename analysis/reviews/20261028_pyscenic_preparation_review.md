# pySCENIC preparation review — 28 October 2026

**Passed preparation checks; this does not establish a completed GRN or motif result.** Both frozen native preparation stages in `runs/20261028_pyscenic_v1/results/pyscenic_native` are complete. This independent review used existing audit tables, metadata, sparse-file headers, reference schemas and recorded file receipts. It performed no model fitting, normalization, raw-matrix reload, source edits or changes to prior analyses.

## Cohort, RNA background and axes

| Lineage | Exact nuclei retained | Original prepared features | Valid unambiguous symbols | Eligible nonconstant RNA | Official TF predictors retained |
|---|---:|---:|---:|---:|---:|
| CA1 | 20,456 | 32,245 | 32,207 | 10,606 | 938 |
| Inhibitory | 2,687 | 32,245 | 32,207 | 11,051 | 973 |

The stored sparse CSC shapes are exactly `20456×10606` and `2687×11051`, with nucleus rows and eligible gene columns. All nucleus IDs, ordering, sample/genotype labels, exact endpoint assignments, sequencing partitions and original full-library UMI denominators match the frozen input metadata. There are six WT and six cKO mice at every endpoint. All 23,143 original selected nuclei remain included; none was downsampled or added. The endpoint counts remain deep 8,750, superficial 11,706, C1 1,718, C6 501 and C8 468.

Every one of the 32,245 feature eligibility decisions per lineage reconciles to the same rule: valid unambiguous symbol, RNA detected in at least 5% of nuclei in any exact endpoint, and a nonconstant normalized profile. No condition, DEG status or HVG ranking is used in this rule. The executed normalization is `log1p(10000 × raw count / original full-gene library UMI sum)`; it does not change the denominator after selecting genes. This auditor checked the denominator fields against frozen metadata and the upstream completed count-reconciliation evidence, rather than independently repeating normalization.

All 14 sealed native preparation output files pass independent byte-count/SHA256 verification. Source SHA256 matches the tested frozen producer (`91a2d534c59fa843f067fc66038ae7eb2039bc9327515aac5d396072d33a6e90`). All included nuclei have more detected eligible RNA genes than the prescribed 5% AUCell rank window; the detection-window audit flags zero nuclei whose window would include undetected RNA.

## TF candidate attrition

The actual unchanged official TF resource contains **1,860 unique names**, including Satb2. For each lineage, 1,830 are assayed, 30 are absent from the assay, and three assayed symbols are ambiguous/invalid. The final RNA detection gate removes another 889 CA1 or 854 inhibitory candidates, leaving 938 and 973 predictors. None was manually added.

Reference coverage imposes a further distinct limit: 645 of the 938 CA1 predictors, and 674 of the 973 inhibitory predictors, have native-filtered annotations for at least one motif actually ranked in these databases. These are possible annotation-supported sources, not predicted final regulon counts. Positive-module support, native database mapping gates, motif enrichment and final target coverage must still be met. A missing final regulon therefore cannot be equated with absent TF activity.

## Why inhibitory SATB2 cannot be inferred in this run

Satb2 is an official candidate, has a single valid assayed symbol, and is nonconstant after normalization in both lineages. It is retained in CA1 but fails the prespecified 5% detection gate in **every** inhibitory endpoint:

| Exact population | Satb2-detected nuclei / total | Detection | Minimum detected nuclei required | Eligibility |
|---|---:|---:|---:|---|
| CA1 deep-associated program | 5,339 / 8,750 | 61.0171% | 438 | Pass |
| CA1 superficial-associated program | 9,096 / 11,706 | 77.7037% | 586 | Pass |
| Inhibitory C1 | 65 / 1,718 | 3.7835% | 86 | Fail |
| Inhibitory C6 | 10 / 501 | 1.9960% | 26 | Fail |
| Inhibitory C8 | 20 / 468 | 4.2735% | 24 | Fail |

The recorded inhibitory reason is `DETECTED_IN_FEWER_THAN_5_PERCENT_IN_EVERY_ENDPOINT`. The native inhibitory expression background and TF predictor list therefore both exclude Satb2. **This fixed run cannot infer an inhibitory SATB2 regulon**, regardless of subsequent motif resources. This is a coverage limitation, not evidence that SATB2 protein or function is absent. The independently preserved CollecTRI target-signature analysis asks a different question and does not require the TF source transcript to pass this nucleus-detection gate. No gate was relaxed after inspecting SATB2.

## Matching reference compatibility

The four selected reference files pass independent SHA256/size checks. The two exact mm10 v10-cluster gene-based ranking resources contain 24,068 and 24,130 mouse-symbol genes, respectively, and the same 5,032 ranked motifs. All ranked motifs have compatible annotations under the declared native thresholds. Approximately 91.84–91.97% of eligible CA1 RNA genes and 92.00–92.15% of eligible inhibitory genes map to the ranking gene axes. Per-module coverage remains the later native pruning decision; overall mapping fractions do not guarantee that any particular module passes.

Satb2 is present in both ranking gene axes. Its native-filtered annotation file has 207 Satb2 rows; 77 rows annotate 23 motifs actually ranked here. Of those 77 descriptions, three report direct annotation and 74 involve orthology and/or motif similarity. This establishes resource availability, not an inferred regulon, mouse-specific occupancy, proven transcriptional activation or SATB2 causality. No full ranking matrix was loaded during this audit; only Arrow schemas and the small motif-index columns were read.

Detailed counts and receipts are in `20261028_pyscenic_preparation_review.json`. GRN/ctx/AUCell outputs require separate review once their actual native stages complete. CA1 labels remain molecular layer programs, C8 retains its mixed-identity warning, and genotype–sequencing-partition confounding remains.
