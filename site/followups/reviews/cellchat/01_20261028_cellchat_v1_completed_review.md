Takeaway: SATB2 CellChat found three exploratory lower-cKO score associations; primary partition sensitivity was inconclusive.

This analysis asked whether the original DEG-supported CA1 and inhibitory-neuron partners show altered predicted communication in Satb2 cKO mice. CellChat completed both the primary triMean calculation and the planned truncated-mean sensitivity using the same twelve mice. Strict matching-endpoint DEG eligibility retained six molecular pairs across nine source–target contexts. Three contexts had lower cKO scores at the primary full-cohort adjusted threshold, with the same direction under the alternative averaging method. Their primary shared-partition comparisons did not meet that threshold, so they remain exploratory leads rather than independently confirmed differential signaling. The numerical outputs passed independent checks, while display defects require a separate figure-layout correction before featured presentation.

## Scientific context and question

The study compares six independent WT and six Satb2 cKO P60 mice. The requested communication directions are CA1→inhibitory, inhibitory→CA1 and a combined view. Existing Seurat endpoints are CA1 deep-associated and superficial-associated molecular programs, plus inhibitory clusters C1, C6 and C8. These labels do not establish anatomical position, temporal order or a new resolution of C8's mixed identity. Genotype partly overlaps sequencing partition; a separate six-WT/three-cKO common-partition comparison is retained.

## Question, lifecycle, and endpoint

The native run completed all four workflow steps with return code zero. Registry: `7e8d72cb-c013-43b2-a7fe-89e271430333`; native run: `snakemake-satb2-cellchat-deg-pathways-20261028-20261028-111445-3b934f90`. Its scientific output directory is `/Users/mujahid357/Documents/codex_new/SATB2_FOLLOWUPS_20261028/runs/20261028_cellchat_v1/results/cellchat`.

CellChat 2.1.2 and its pinned bundled mouse database were used with full-library log1p-CP10k RNA, no imputation, `population.size=FALSE`, minimum group size ten, seed 1337 and one hundred native specificity permutations. A common sixty-mouse×endpoint object keeps expression scaling comparable; only matching within-mouse CA1↔inhibitory edges are reported. Native cell-label specificity probabilities are distinct from genotype tests. Genotype comparisons use unthresholded eligible scores and exact independent-mouse label permutations, with separate BH families for each averaging method, level and cohort. No original DESeq2 model was refitted.

## Key findings

- Exact inputs: 23,143 nuclei, comprising 20,456 CA1-program nuclei and 2,687 inhibitory nuclei; all five endpoints occur in all twelve mice. The smallest mouse×endpoint group has thirteen nuclei. Every integer-matrix column sum matches its saved full original RNA library denominator.
- All 40,548 candidate-context eligibility decisions were checked independently. Nine contexts, representing six molecular pairs and five pathway labels, pass the original matching-source ligand and matching-target receptor DEG rule. Every required complex subunit must have original DESeq2 `padj<0.05`; none of the retained contexts has a cofactor modifier.
- The retained pairs and original gene directions are Ptn↓–Ptprz1↑, Nrxn1↓–Nlgn1↑, Slit1↓–Robo1↓, Sema6d↓–Plxna1↓, Tenm3↓–Adgrl2↓ and Slit2↓–Robo1↓. The complete endpoint mapping and component evidence are preserved in the eligibility table and independent component-evidence audit.
- Primary triMean: three of nine full-cohort context tests pass BH 0.05; zero of nine common-partition tests pass. The planned truncated-mean sensitivity yields seven of nine full-cohort and five of nine common-partition associations. These are separate exploratory families, not independent validations.

| Primary full-cohort context | Mean score difference, cKO−WT | Full-cohort BH | Common-partition BH |
|---|---:|---:|---:|
| CA1 deep-associated program → IN C6, Slit1–Robo1 | −0.018037387 | 0.0064935065 | 0.071428571 |
| IN C6 → CA1 deep-associated program, Tenm3–Adgrl2 | −0.033767803 | 0.0064935065 | 0.071428571 |
| IN C6 → CA1 superficial-associated program, Tenm3–Adgrl2 | −0.042087752 | 0.0064935065 | 0.071428571 |

All three retain their direction and full-cohort BH support under both averaging methods. Each retained pathway context contains exactly one qualifying pair in this strict subset, so identical LR and pathway outcomes are the same evidence. The combined direction view transports the two directional results and does not create another test family. C8 has no eligible strict-DEG edge; this does not establish absence of communication.

## Interpretation

Slit1–Robo1 from the deep-associated CA1 program to C6 and Tenm3–Adgrl2 from C6 to both CA1 programs are the clearest conditional leads. Their common-partition primary BH values remain above 0.05, despite concordant direction, and the sensitivity calculation is more permissive in observed outcomes. Thus neither robust partition-independent genotype causality nor confirmed differential protein signaling is established. Other sensitivity-only associations should remain explicitly method-sensitive. Component gene up/down directions and the composite interaction-score direction are separate quantities; one does not mechanically determine the other.

The analysis meets the requested literal both-end-DEG restriction. That restriction deliberately excludes most possible signaling and reuses the same study's differential-expression selection. Subsequent score tests are therefore exploratory post-selection comparisons. The pathway labels organize database-supported machinery; they are not evidence of pathway activation or the full pathway's behavior.

## Artifacts

All paths below are under the scientific output directory stated above unless noted otherwise.

**Results:** `interaction_endpoint_eligibility_and_signs.tsv`; `triMean/strict_DEG_within_mouse_LR_scores.tsv`; `triMean/strict_DEG_within_mouse_unthresholded_pathway_scores.tsv`; the two corresponding `exploratory_*_mouse_permutation.tsv` tables; and the same files under `truncatedMean10/`. All three direction views and native RDS objects are preserved. Each method has 108 per-mouse LR records with exact native score/probability transport.

**QC and visual reports:** `FINAL_VALIDATION.json` records complete native output validation. Thirty-eight one-page PDFs and thirty-eight PNGs cover dotplots, pathway tiles, circos, circle networks and mouse-level boxplots. These original figures are preserved with the display limitations below; they are not all presentation-approved. Independent audits and rendered PDF checks are in `/Users/mujahid357/Documents/codex_new/SATB2_FOLLOWUPS_20261028/reviews/CellChat_actual_readonly/`.

**Provenance:** `METHODS.json`, `executed_config.json`, `sessionInfo.txt`, database tables/native RDS, `mouse_roster.tsv`, `nucleus_metadata.tsv`, feature mapping, frozen source/config and the native plan/run receipts. The terminal Workbench record is saved in `/Users/mujahid357/Documents/codex_new/SATB2_FOLLOWUPS_20261028/receipts/cellchat_terminal_record.json`. Independent review verified all 115 final file receipts, original DEG gates, exact 924/84 permutation arithmetic, separate BH families and RDS/TSV/pathway transports.

**Model synthesis:** This review and `20261028_final_scientific_review.md` in the independent audit directory explain the findings, uncertainty and figure defects. Native completion alone did not establish acceptance.

## Limitations and blockers

Partial genotype/partition confounding, the smaller common-partition cKO cohort and same-data DEG selection limit inference. RNA-based scores do not prove protein expression, binding, synaptic transmission, functional E/I balance, anatomical layers or causal genotype effects. Existing original DESeq2 fits and label uncertainty remain unchanged.

The original figure outputs have presentation defects. All twelve circos headers extend above the PDF page bounds. Long heatmap context labels are clipped or too tight in PNG exports. Both mouse-boxplot groups use default vertical jitter, which can draw genuinely zero scores below zero; the actual TSV/RDS scores are nonnegative. The circle figures need an explicit readable context caption. These issues affect display, not the independently checked numerical outputs. Original files must remain intact, and corrected figures must be saved additively from the exact existing values without new inference.

## Recommended next step

Produce the minimal separate presentation-only correction: zero vertical jitter, adequate device margins/canvas and in-page headers/captions, preserving all numerical values, scales, eligibility, statistics and original files. Review the corrected PDFs before linking them as the featured atlas downloads. Present the three primary leads as exploratory candidates for independent follow-up, with full and common-partition outcomes side by side; do not expand the DEG rule or tune thresholds to obtain additional interactions.
