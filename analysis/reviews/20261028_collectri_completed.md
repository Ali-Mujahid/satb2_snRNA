Completed curated CollecTRI TF target-signature analysis in all five populations; motif inference remains unavailable.

This analysis estimated curated TF target-signature scores for the exact two CA1 programs and inhibitory clusters 1, 6 and 8. It used original integer RNA pseudobulk from all twelve mice, unchanged original DESeq2 size factors and a genotype-blind full eligible RNA background. Native decoupler ULM scores were then compared with separate continuous mouse-level models, preserving six WT and six cKO independent replicates. There are 2,768 tested population–TF combinations and 23 global FDR-significant associations, with no FDR-significant associations in the smaller matched-partition sensitivity. SATB2 target-signature activity is not FDR-significant in any population, including the separate self-target-exclusion sensitivity. The curated analysis and its packaging completed successfully, while pySCENIC motif inference was not run because matching mouse references were unavailable and the user requested continuation of available analyses.

## 1. Analysis context

Registry aa3fe6be-fb68-40a0-8f30-8f4db2c5e934; run snakemake-satb2-collectri-tf-20261026-20261026-103502-75b1b8b7. Existing user-selected Snakemake controller and isolated decoupler 2.3.0 runtime. The frozen mouse-mapped CollecTRI prior contains 43,226 links, 1,165 sources and 6,582 targets. Exact selected memberships match hdWGCNA: 8,750 deep-associated CA1, 11,706 superficial-associated CA1, and 1,718/501/468 inhibitory C1/C6/C8 nuclei. The original gene-DE method remains DESeq2; no score is fitted with DESeq2 and no edgeR or new gene-DE fit is used.

## 2. Lifecycle and verification

Durable native state is completed, return code 0; all four activity/mouse-association/package/final tasks completed. FINAL_VALIDATION.json reports COMPLETE_WITH_BLOCKED_MOTIF, verified original source/input hashes, finite native mouse scores, exact memberships and unique mouse keys. The final native/package receipts cover the scientific outputs; controller logs are external mutable metadata. Native score, target coverage, RNA/Wald input audits and all output tables were inspected. An independent read-only audit verified 227 native/package SHA receipts and reproduced saved linear-model arithmetic and BH adjustment within floating-point precision without fitting new scientific models.

## 3. Findings

Tested signatures: deep CA1 622, superficial CA1 638, C1 563, C6 466, C8 479. Global BH across all 2,768 tests gives respectively 7, 13, 1, 2 and 0 associations at FDR below 0.05. The inhibitory candidates are Maf lower in C1 and Klf15/Rest lower in C6. All 23 retain the effect direction in twelve leave-one-mouse-out fits and in the same-partition subset, but none passes that subset's separately corrected FDR. Full-prior Satb2 FDR is 0.999/0.381/0.943/0.902/0.611 for DL/SL/C1/C6/C8. Removing its self-target link yields no significant association across the five-endpoint sensitivity family; DL FDR is 0.162 and SL 0.107 (SL nominal P 0.0213). No claim of proven SATB2 activity loss is supported.

## 4. Interpretation

Scores are native ULM t values measuring correspondence between a signed curated target prior and the measured RNA pattern. A score is distinct from source-TF RNA, source-TF protein, direct motif binding and a genotype log2 fold change. Native within-mouse gene-set FDR is separate from the genotype-association FDR. Thirteen of 23 significant signatures have source-TF RNA below the genotype-blind expression threshold, including five with zero source counts across all twelve mice in superficial CA1 (Dlx4, Foxh1, Hoxa5, Prop1, Tcf21); their target-pattern scores do not establish TF activity in those neurons. Several prior sources have identical signed target sets, including Dmtf1/Dmtf1l where modeled, and cannot be independently distinguished by these scores. None of the 23 significant signatures belongs to an exact duplicate-target-set group.

## 5. Artifacts inspected

Run root: /Users/mujahid357/Documents/codex_new/SATB2_REGULATORY_20261026/runs/20261026_collectri_v1. Key outputs: tf_activity/TF_ACTIVITY_PER_MOUSE.tsv; TF_TARGET_COVERAGE.tsv; TF_WALD_CONTRAST_ENRICHMENT.tsv; tf_mouse/TF_ASSOCIATIONS.tsv, TF_SAME_PARTITION_ASSOCIATIONS.tsv, TF_LOO_STABILITY.tsv and SATB2_NO_SELF_TARGET_ASSOCIATIONS.tsv; native used signed prior edges, eligibility/normalization audits, COMPLETE and FINAL_VALIDATION receipts. Native SATB2 self-target and mouse-score figures were inspected. The sensitivity figure retains native tick-label precision; interactive exports format display labels independently without changing values.

## 6. Limits

The human-derived curated prior is mapped to mouse using the documented native resource procedure; links are not all direct mouse-neuron binding evidence and some sources are regulatory cofactors. Target coverage varies substantially between endpoints. Partial genotype/partition confounding and only three cKO mice in the matched-partition check limit inference. Shared target genes and reused RNA mean agreement with gene DE or modules is not independent validation. Removing a self-edge is a sensitivity, not an independent assay. IN_C8 remains mixed. No per-nucleus activity was computed, so no activity UMAP is substituted by TF RNA or cluster averages. The prepared pySCENIC environment does not constitute motif-supported inference; both official reference hosts refused connections and no matching local cache exists.

## 7. Next action

Integrate curated scores, actual target networks, coverage, mouse plots, original saved associations and sensitivity figures into the additive GitHub/offline explorer. Preserve explicit unavailable statuses and NA values. Present the 23 associations as exploratory target-signature leads, with the matched-partition and source-RNA caveats immediately available. Keep the three communication directions strictly restricted to original DEGs at both exact endpoints and every complex subunit; the full eligible background here serves a different statistical method.
