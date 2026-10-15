# CellBender audit: source and correction diagnostics

The complete twelve-sample archive supports a **paired corrected-count sensitivity**, with original counts retained for the full primary comparison. Use the archive sources recorded in `usable_count_sources.tsv`. Older June outputs under `Desktop/data/snr/cellBender` are historical and differ from this archive: for example, their cKO_933 filtered matrix has 9,946 barcodes, versus 11,402 in the archive.

All twelve archive sample folders contain full and filtered corrected matrices, training logs, metrics, PDF summaries, cell-barcode CSVs and posterior files. The full/filtered/log/metrics/PDF/barcode source files are SHA256-bound in `source_artifact_inventory.tsv` (72 artifacts). All 24 archive original Cell Ranger raw/filtered hashes match the previously accepted original-source audit. That accepted count-integrity audit was reused rather than rerun. This task read local corrected headers, diagnostic metadata, sparse count chunks and fixed-roster marker totals; no scientific model was fitted or original file changed.

## Fixed-roster coverage

The frozen primary cohort contains 108,062 nuclei. Local full-output sparse pointers and filtered counts reproduce the earlier compatibility result exactly: **104,298 have usable analyzed, nonzero corrected profiles; 3,764 are zeroed, including 707 not analyzed**. The filtered outputs have 138,494 cell calls; their extra barcodes are outside the fixed primary comparison.

| Sample | Primary nuclei | Paired usable | Zero corrected | Not analyzed |
| --- | ---: | ---: | ---: | ---: |
| cKO_1138 | 7,943 | 7,687 | 256 | 11 |
| cKO_1139 | 9,750 | 9,465 | 285 | 16 |
| cKO_1160 | 8,221 | 7,874 | 347 | 200 |
| cKO_1165 | 9,851 | 9,444 | 407 | 17 |
| cKO_933 | 8,896 | 8,625 | 271 | 72 |
| cKO_991 | 8,470 | 7,950 | 520 | 88 |
| wt_1122 | 8,546 | 8,184 | 362 | 168 |
| wt_1123 | 9,201 | 8,968 | 233 | 0 |
| wt_1178 | 6,660 | 6,407 | 253 | 135 |
| wt_1220 | 10,430 | 10,176 | 254 | 0 |
| wt_1226 | 9,733 | 9,362 | 371 | 0 |
| wt_1228 | 10,361 | 10,156 | 205 | 0 |

The zeroed profiles had original sample medians of 623–927 UMIs and 454–662 detected genes. These descriptive observations do not justify retrospectively changing QC. All filtered corrected counts are finite, nonnegative integers, and ordered stable feature IDs match the original matrices. Each sample's complete primary roster, including unusable corrected profiles, is preserved in its `*_original_cellbender_primary_paired_roster.tsv`.

## Identity-specific diagnostic observations

Frozen high-support reference memberships retain 8,903/8,911 IN nuclei, 60,849/60,985 EX nuclei and 27,672/27,703 nested CA1-ProS nuclei. These are prior reference predictions, not new cell-type conclusions or anatomical validation. The full coverage table includes lower-support broad predictions and subclasses; their losses are larger than high-support group losses.

Across the identical 104,298 supported nuclei, original UMIs total 846,367,233 and corrected UMIs total 793,300,571 (6.27% removed). Whole-cohort marker removal is strongly affected by expression outside the expected identity, so within-identity comparisons are more informative:

- In the paired high-support IN cohort, Gad1/Gad2 lose approximately 0.1% of UMIs, Pvalb 1.6%, Sst 2.0%, Vip 1.1% and Reln 0.18%. The off-identity excitatory marker Slc17a7 loses 99.1%.
- In the paired CA1-ProS cohort, Slc17a7 loses 2.27%, Satb2 0.35%, Ndst4 0.76%, Col11a1 0.06%, Asic2 2.62% and Calb1 2.22%. Gad1/Gad2/Pvalb/Sst lose 96–99% of their small off-identity UMI totals.
- **Cck remains a specific review flag**: IN Cck UMIs fall from 4,301 to 1,905 (55.7%), with detected nuclei decreasing from 2,880 to 1,135. CA1-ProS Cck loses 22.6%. This low-copy signal needs subtype-aware comparison using other Cck-class markers and original counts; its removal alone cannot establish either successful correction or loss of true identity.

This pattern supports plausible ambient suppression while retaining most canonical neuronal signal. It does not establish gene-specific truth. Marker counts and detection totals are saved both per mouse and aggregated in the `marker_removal*` tables.

## Training, calling and provenance flags

All logs describe CellBender 0.3.0, 150 epochs, GPU use, FPR 0.01 and MCKP output estimation. All twelve training/test HDF5 curves are finite through epoch 150. The rendered PDF panels show late plateaus with no terminal collapse. Several probability tails are broad; wt_1122 and wt_1178 have the largest found/expected ratios (4.56 and 4.53) and recorded within-called-cell removal fractions (14.0% and 12.7%). Treat these metrics as review flags and retain the frozen roster.

The recorded `convergence_indicator` spans 0.154–2.041. Version 0.3.0 defines it from recent training-ELBO changes relative to recent variability; it is not a universal numerical acceptance threshold. Its denominator can also make comparisons misleading. Count-removal metrics refer to the CellBender-defined called roster, which differs from the paired fixed roster. [Version 0.3.0 metric definitions and output ordering](https://raw.githubusercontent.com/broadinstitute/CellBender/v0.3.0/cellbender/remove_background/run.py)

All twelve logs record HTML-report generation errors followed by completion. Matrices were written before report generation; these errors alone do not establish count corruption. Current archive logs show one successful checkpoint load (cKO_1138), eight rejected checkpoint-workflow-hash mismatches and three absent checkpoints (wt_1123/wt_1220/wt_1226). The historical compact note's statement of four loaded checkpoints is not reproduced by the actual local logs. Preserve cKO_1138's resumed-training provenance caution; no checkpoint was unpickled or executed here.

## Use in the requested reanalysis

Use the original 108,062-nucleus roster as the primary starting cohort for both Scanpy and Seurat. Run a separately labeled correction sensitivity on the **same 104,298 supported barcode overlap in both original and corrected branches**. Keep all twelve mouse replicates and their sample IDs. Do not substitute originals for zero corrected profiles within a purported corrected matrix, import extra CellBender cell calls automatically or silently drop unavailable profiles. Compare Cck subtype annotation and markers explicitly before accepting correction-dependent conclusions.

The audit receipt, count-source manifest, source hashes, complete rosters, reference-coverage tables, recorded log/metric flags, full learning curves and PDF source locators are in this directory. The full/filtered sparse-count equivalence and historical remote source-compatibility result remain the earlier accepted compatibility audit; this task did not repeat that entire accepted audit. Current local source digests now bind the exact files used for new analysis.
