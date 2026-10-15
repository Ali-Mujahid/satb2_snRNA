# Actual input and QC review

The active v1_0_5 input freeze supports the matched twelve-mouse reanalysis. All 22 declared input artifacts and eight recorded provenance sources match their hashes; current primary, paired and mapper-prediction IDs reconcile exactly. The review covers frozen QC/roster consistency, supplied marker mappings and the QC report. Source-bound cached evidence confirms all twelve MapMyCells queries used these exact original integer counts and feature/nucleus axes. Existing QC/singlet calls and mapper labels were reused; no new doublet model or identical mapper fit was performed. This input acceptance does not establish anatomical annotation, cluster accuracy or differential-expression validity.

Reviewed active run: `runs/20261013_v1_0_5`, native registry `9702cbef-00e6-4f64-aa37-e9b62357bd71`. The [input completion receipt](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/runs/20261013_v1_0_5/results/inputs/COMPLETE.json) records `COMPLETE`. Current checks are recorded in the [active review provenance](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/reviews/input_qc_review_active_v1_0_5.provenance.json) and [22-artifact hash table](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/reviews/input_qc_active_artifact_hashes.tsv). The earlier v1_0_1 source-cell and figure inspections remain historical evidence; current numerical QC, roster, marker-panel and feature tables are byte-identical to those checked earlier. The current PDF has a different file hash, but all six re-rendered pages are pixel- and PNG-byte-identical to the earlier inspected figures.

## Counts and frozen QC

The actual [QC summary](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/runs/20261013_v1_0_5/results/inputs/QC_RETENTION_AND_COUNT_SUMMARY.tsv) totals 121,486 Cell Ranger called nuclei, 108,062 retained primary nuclei and 7,338 recorded doublet calls. The primary exclusion count is 13,424. Doublet calls and other QC criteria can overlap, so recorded doublets are not a separate additive explanation of all exclusions. The accepted earlier QC/singlet calls were reused; this stage did not fit a new doublet model or introduce fresh exclusions.

Both paired rosters contain exactly 104,298 unique nucleus IDs, cover the same twelve mice and are exact subsets of the 108,062-nucleus primary roster. Their roster files and label files are byte-identical between paired-original and paired-CellBender branches. Six independent WT mice and six independent cKO mice remain in every branch.

| Genotype | Called | Primary QC | Paired support | Corrected-profile attrition | Attrition / primary |
| --- | ---: | ---: | ---: | ---: | ---: |
| WT | 61,849 | 54,931 | 53,253 | 1,678 | 3.055% |
| cKO | 59,637 | 53,131 | 51,045 | 2,086 | 3.926% |
| Total | 121,486 | 108,062 | 104,298 | 3,764 | 3.483% |

| Mouse | Primary QC | Paired support | Attrition | Attrition / primary |
| --- | ---: | ---: | ---: | ---: |
| cKO_1138 | 7,943 | 7,687 | 256 | 3.223% |
| cKO_1139 | 9,750 | 9,465 | 285 | 2.923% |
| cKO_1160 | 8,221 | 7,874 | 347 | 4.221% |
| cKO_1165 | 9,851 | 9,444 | 407 | 4.132% |
| cKO_933 | 8,896 | 8,625 | 271 | 3.046% |
| cKO_991 | 8,470 | 7,950 | 520 | 6.139% |
| wt_1122 | 8,546 | 8,184 | 362 | 4.236% |
| wt_1123 | 9,201 | 8,968 | 233 | 2.532% |
| wt_1178 | 6,660 | 6,407 | 253 | 3.799% |
| wt_1220 | 10,430 | 10,176 | 254 | 2.435% |
| wt_1226 | 9,733 | 9,362 | 371 | 3.812% |
| wt_1228 | 10,361 | 10,156 | 205 | 1.979% |

Original primary UMIs total 849,122,137. The same paired profiles have 846,367,233 original UMIs and 793,300,571 corrected UMIs, a weighted reduction of 6.270%. Within paired WT the reduction is 6.454%; within paired cKO it is 6.116%. These are descriptive capture/count-source comparisons, not genotype tests or evidence of differential gene expression. The [36-row branch coverage table](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/runs/20261013_v1_0_5/results/inputs/COUNT_BRANCH_COVERAGE.tsv) preserves each sample's count-source hash, nucleus count and library total.

Corrected-profile attrition agrees with the separately completed [CellBender audit](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/audit_cellbender/CELLBENDER_AUDIT_REVIEW.md): 3,764 primary profiles are zero corrected, including 707 not analyzed. The paired branches preserve their supported overlap; missing corrected profiles remain documented in the full primary cohort. They must not be filled with original values inside a purported corrected matrix.

## QC figure inspection

All six pages of the active [QC_REVIEW.pdf](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/runs/20261013_v1_0_5/results/inputs/QC_REVIEW.pdf) were re-rendered from the current PDF and matched exactly to the historically inspected pages; [render checks](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/reviews/input_qc_png_active_v1_0_5/render_checks.json) retain both image hashes. Axes, sample names and numerical scales are readable, with no clipping or missing panels. Pages cover called/retained/paired counts, original UMI and gene-complexity distributions, original mitochondrial percentages, original complexity scatterplots with frozen-QC exclusions highlighted, and same-profile CellBender removal distributions.

The complexity scatterplots show high-complexity tails and excluded profiles in multiple mice, including a particularly high detected-gene tail for wt_1178. These plots support review of the existing policy; they do not establish that any remaining outlier is a doublet. The boxplots summarize distributions rather than displaying every outlier. The CellBender removal distribution is broad in each mouse and is widest for wt_1122, consistent with the earlier correction diagnostics. The PDF describes earlier audited QC applied to verified counts, not an independent new QC fit.

## Supplied workbook marker mappings

The actual [review panel](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/runs/20261013_v1_0_5/results/inputs/marker_review_panels.tsv) has 147 rows and 140 unique gene symbols: 22 canonical coarse-lineage genes, 20 CA1 population candidates, 15 superficial-associated VEH genes, 15 deep-associated VEH genes, and five genes for each of the 15 supplied interneuron/program labels. All 125 workbook-derived rows point to the exact matching gene cell in the original extracted workbook contents. All three original XLSX files still match their previously recorded binary hashes.

The two CA1 workbooks are duplicate source tables across all 13 worksheets, so only `CA1.markers.xlsx` is used for selection. CA1 population candidates come from `PyrCA1_markers` with positive fold change, adjusted P <=0.05 and a detection-rate difference >0.1, ranked by specificity. Layer programs come from control/VEH superficial-versus-deep comparisons in `DEGs_PyrCA1Sup_vs_PyrCA1Deep`, preserving the sign: positive is superficial-associated, negative deep-associated. They are explicitly VEH-associated expression programs rather than anatomical layer assignments. Satb2, mitochondrial and ribosomal genes are excluded from those layer programs. Epilepsy-state/trajectory DEG sheets are not treated as general identity marker panels.

Interneuron entries retain their exact `Sheet1` row locators and require positive fold change, adjusted P <=0.05 and detection-rate difference >0.1. The two `Module` labels remain expression programs. The supplied `Cholinergic neurons` label is a source-program review category; its presence must not force a GABAergic cell identity. All panels declare expression review only and do not determine the frozen cohorts or impose cell/layer labels.

**Mapping flag:** `AC129186.1` has no exact symbol in the actual [gene feature table](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/runs/20261013_v1_0_5/results/inputs/gene_features.tsv). It occurs twice in the panels: `CA1_population_candidates` and `Superficial_associated_VEH_program`. Thus 145 panel rows covering 139 unique symbols map to one stable feature ID each; two rows remain unavailable for exact-symbol expression plotting. No alias or stable gene ID was invented. The supplied workbook row must remain recorded even when the gene is unavailable in this reference. Panel availability/plotting outputs should expose this limitation.

Original workbook sources: [CA1.markers.xlsx](</Volumes/Neurowissenschaften/-Ali_on_server/snr_p60/CA1.markers.xlsx>), [CA1 superficial/deep source](</Volumes/Neurowissenschaften/-Ali_on_server/snr_p60/CA1.Superfacial vs deep layer markers.xlsx>), and [interneuron source](</Volumes/Neurowissenschaften/-Ali_on_server/snr_p60/Interneuron classification marker.xlsx>). Exact per-gene sheet/row references are preserved in the panel table.

## Exact MapMyCells query provenance

The active [original shared-stage receipt](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/runs/20261013_v1_0_5/results/shared/original_full/COMPLETE.json) records twelve `mapmycells_query_verification` results with status `VERIFIED_EXACT_ORIGINAL_INTEGER_COUNTS`, normalization `raw`, exact gene and nucleus IDs, each complete sample-by-32,285-gene shape and zero unequal sparse count entries. The sample shapes sum to 108,062 nuclei; this comparison covers the full gene-by-nucleus payload rather than only library totals or a selected feature panel. Each result binds the original source-count SHA, saved query SHA, mapper-config SHA and historical mapper completion SHA. The active prediction table has 108,062 unique IDs, equals the authoritative primary roster, and has the same verified prediction hash recorded by the shared-stage receipt. Per-capture bindings are retained in the [query evidence table](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/reviews/input_qc_MapMyCells_exact_query_evidence.tsv).

This is accepted source-bound cache evidence from the preserved v1_0_4 shared stage. The active run records `raw_count_sources_reloaded=false` and `scientific_preparation_refit_performed=false`; this review did not reread the large saved query matrices. The [original verifier](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/runs/20261013_v1_0_5/workflow/scanpy/common.py:83) checks raw mapper settings, the exact ordered gene axis, unique/exact nucleus IDs, nonnegative integer query values and all sparse entries after nucleus alignment. The [active reuse guard](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/runs/20261013_v1_0_5/workflow/scanpy/prepare_common_inputs.py:125) requires the accepted zero-mismatch evidence and full shape while matching current source hashes, mouse identities, library totals and prediction hash. The parent shared-stage receipt hash also matched. This supports reuse of the same original-query labels; it does not provide an independent annotation truth set or a new denoised-query mapping.

## Remaining interpretation limits

Existing same-original-query MapMyCells predictions were verified and reused as common annotation evidence; the input stage did not rerun an identical mapping. These labels are not independent biological truth. This review confirms the matched starting cohorts and report provenance. New normalization/integration, fresh clustering, tool agreement, subtype/layer support and mouse-level inference require their own actual outputs and reviews.
