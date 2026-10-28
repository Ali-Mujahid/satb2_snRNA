# CollecTRI TF activity: independent final audit

Completed: 23 primary target-signature associations were found; none passed the same-partition sensitivity FDR.

Audit date: 2026-10-28. Scope: read-only inspection of the completed native run `runs/20261026_collectri_v1/results`; no study fits, reference changes, native file edits, or exclusions were performed. This is a curated-prior target-pattern analysis, not completed motif-supported regulon inference. The statistical replicate is the independent mouse. Primary comparisons contain six WT and six cKO mice at each exact endpoint. The unavailable mouse cisTarget references still prevent pySCENIC motif inference, and the user authorized continuing the available analyses.

## Completion and numerical checks

- `FINAL_VALIDATION.json` reports `COMPLETE_WITH_BLOCKED_MOTIF`; TF activity ran, motif regulons did not, and no nucleus activity UMAP was packaged.
- All **227** recorded SHA-256/byte receipts checked independently: **222 native files and five packaged atlas files**, with no mismatch.
- **33,216** finite native mouse scores form **2,768** unique endpoint×TF groups. Every modeled group has exactly 12 unique mice, six per genotype; there are no duplicate endpoint×TF×mouse keys.
- The all-prior association table retains **5,825** rows (1,165 candidate sources × five endpoints), including **3,057** unmodeled endpoint×TF rows with unavailable statistics. These are missing values, not zero activity or evidence of no effect.
- Primary coefficient, pooled-variance SE, t statistic, two-sided p value, and 95% CI were independently checked against the saved mouse scores using arithmetic, without refitting. Maximum errors were respectively 6.44e-15, 1.64e-15, 1.72e-13, 8.78e-14, and 7.93e-15.
- BH adjustment independently recomputed from the native p values agreed within 3.00e-15 for both 2,768-test families. Both separate five-endpoint SATB2 sensitivity families agreed within 9.44e-16.

## Primary results and robustness

Effects are cKO minus WT in **native ULM t-score units**, not gene log2 fold changes. FDR is adjusted over all 2,768 tested endpoint×TF associations, rather than separately within each endpoint.

| Exact population | Nuclei | Background genes | Tested sources | Primary FDR < 0.05 | Higher/lower in cKO | Same-partition FDR < 0.05 |
|---|---:|---:|---:|---:|---:|---:|
| CA1 deep-associated program | 8,750 | 15,899 | 622 | 7 | 2 / 5 | 0 |
| CA1 superficial-associated program | 11,706 | 16,969 | 638 | 13 | 3 / 10 | 0 |
| Inhibitory C1 | 1,718 | 12,749 | 563 | 1 | 0 / 1 | 0 |
| Inhibitory C6 | 501 | 8,980 | 466 | 2 | 0 / 2 | 0 |
| Inhibitory C8 | 468 | 9,203 | 479 | 0 | 0 / 0 | 0 |

All **23** primary hits retain their effect direction in **12/12** leave-one-mouse-out comparisons, all have the same direction in the six-WT/three-cKO same-sequencing-partition analysis, and all 23 exploratory 200-resample mouse bootstrap intervals exclude zero. However, **none** passes the separate global same-partition BH threshold; its smallest adjusted p value is 0.101437. Directional consistency and bootstrap intervals are sensitivity summaries, not independent replication or multiplicity-adjusted confirmation. The smaller sensitivity sample also has reduced power, so failure there does not prove the primary effects absent.

The primary candidate signatures are:

- **CA1 deep:** lower Myrf, Pgr, Tfcp2, Mxd1, Barx2; higher Cdx1, Smad7. The strongest adjusted p values are Myrf 0.00737, Pgr 0.00773, and Tfcp2 0.00805.
- **CA1 superficial:** lower Barx2, Tcf21, Mxd1, Dlx4, Hoxa5, Foxh1, Nfatc4, Gli1, Prop1, Ewsr1; higher Irf7, Tbr1, Arid3a. Barx2 has adjusted p = 0.00473.
- **Inhibitory C1:** Maf lower, effect −0.38009, adjusted p = 0.03990.
- **Inhibitory C6:** Klf15 lower, effect −0.80478, adjusted p = 0.02841; Rest lower, effect −0.70399, adjusted p = 0.03684.
- **Inhibitory C8:** no primary global-FDR hit. This is not evidence that no regulation changes; C8 also has mixed identity and only 13 nuclei in cKO_933, although that mouse was correctly retained for scoring.

## SATB2 and self-target sensitivity

The curated prior includes a Satb2→Satb2 edge. Eligible full-prior target counts are 17/23/13/9/11 for deep/superficial/C1/C6/C8. Removing this edge leaves 16/22/12/9/11 targets; the source transcript remains in the otherwise unchanged RNA background. This is a controlled score sensitivity, not independent biological evidence.

| Endpoint | Full-prior effect | Full-prior global FDR | No-self effect | No-self five-endpoint FDR | No-self LOO direction |
|---|---:|---:|---:|---:|---:|
| CA1 deep | −0.00056 | 0.99921 | −0.12924 | 0.16157 | 12/12 |
| CA1 superficial | −0.22335 | 0.38085 | −0.31337 | 0.10659 | 12/12 |
| C1 | +0.03342 | 0.94312 | +0.04083 | 0.76893 | 10/12 |
| C6 | +0.06229 | 0.90203 | +0.06229 | 0.76893 | 11/12 |
| C8 | +0.24332 | 0.61121 | +0.24332 | 0.35997 | 12/12 |

**No SATB2 association is FDR-significant in either variant or its same-partition sensitivity.** Superficial no-self scoring has nominal p = 0.02132, but adjusted p = 0.10659; it is an exploratory downward target-signature tendency. Deep scoring is materially sensitive to the self-edge: the full-prior effect is essentially zero and its LOO direction is only 6/12, whereas the no-self effect is negative. The different FDR families mean a smaller no-self adjusted p value alone cannot be interpreted as stronger replication. Inhibitory scores do not support a significant SATB2 activity increase or decrease.

## Identifiability and source-expression guardrails

Exact endpoint-specific signed target sets were grouped independently. There are **7/8/5/4/4 duplicate pairs** in deep/superficial/C1/C6/C8, leaving 615/630/558/462/475 distinct signed target sets. Every duplicate pair has exactly identical saved mouse scores. Dmtf1/Dmtf1l share nine/seven/six targets in deep/superficial/C1; other repeated pairs include Eno1/Eno1b, Hmga1/Hmga1b, Kdm5c/Kdm5d, and Mycn/Mycs. None of the 23 primary FDR hits belongs to these duplicate groups. Preserve the native records, but describe such pairs as indistinguishable prior signatures, not two independently identified regulators.

**13 of 23 significant signatures have source RNA below the eligibility threshold** of ≥5 counts in ≥3 mice. Five sources have zero raw counts in all 12 superficial-program pseudobulks: **Dlx4, Foxh1, Hoxa5, Prop1, and Tcf21**. ULM can score their curated targets even when the source is undetected; that is valid target-pattern inference but cannot establish that these TF proteins are active in these neurons. Barx2 in both CA1 programs, Cdx1 in deep, Gli1/Irf7/Nfatc4 in superficial, and Klf15/Rest in C6 are also below the source RNA threshold. Myrf's RNA is eligible, but its CA1 target signature still requires cell-identity/contamination review before a neuron-intrinsic mechanistic claim. No source RNA filter was added after seeing the results.

## Methodological interpretation

The executed source correctly uses the matching original DESeq2 size factors, full eligible measured RNA backgrounds including non-DEGs/non-target genes, genotype-blind coverage filtering, native decoupler ULM with `tmin=5`, and mouse-level continuous-score linear models. It does not fit DESeq2 to continuous scores, re-estimate size factors after target selection, treat nuclei as replicates, or select scores by within-mouse gene-set FDR. Native gene-set enrichment FDR and genotype-association FDR are different quantities. The supplementary original-DESeq2-Wald enrichment uses the same original experiment and is not independent validation. Differences in available targets/backgrounds also limit absolute score comparisons between cell populations.

The defensible conclusions are the saved target-pattern estimates, their mouse coverage, the 23 primary associations under the stated model, and the explicit uncertainty revealed by sequencing-partition and SATB2 sensitivities. All mechanistic regulator claims remain exploratory: the prior partly uses cross-species orthology, overlapping signatures can identify the same expression pattern, motif-supported regulons are unavailable, and sequencing partition is partly confounded with genotype. The DEG-only requirement continues to apply to the separately preserved ligand–receptor analyses; it would be statistically inappropriate to apply it to the gene background of this regulatory analysis.

## Artifacts inspected

Native tables: `tf_mouse/TF_ASSOCIATIONS.tsv`, `TF_SAME_PARTITION_ASSOCIATIONS.tsv`, `TF_MOUSE_SCORES.tsv`, `TF_LOO_STABILITY.tsv`, `TF_MOUSE_BOOTSTRAP.tsv`, `SATB2_NO_SELF_TARGET_*`, and `MOUSE_PARTITION_AUDIT.tsv`; `tf_activity/TF_TARGET_COVERAGE.tsv`, endpoint-specific gene-background audits and used signed prior edges; final and step completion receipts. The actual `tf_mouse/figures/Satb2_self_target_sensitivity.png` was visually inspected and matches the native coefficients/CIs. Its long floating-point x tick labels are cosmetic; the underlying values and uncertainty are intact. The source and native outputs were left unchanged.
