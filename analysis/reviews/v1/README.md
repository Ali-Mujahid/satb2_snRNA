# SATB2 all-12 Scanpy / Seurat comparison

Status: COMPLETED and reviewed, 13 October 2026. All 34 native tasks finished; broad computational agreement is high and fine CA1/ProS identities remain provisional.

This project preserves the previous analyses and works from the audited common primary QC/singlet roster: 108,062 nuclei from six WT and six Satb2 cKO mice. The original-count branch is primary. CellBender sensitivity uses the same 104,298 usable nuclei in its original-count and corrected-count branches, without replacing unavailable corrected profiles with original counts.

MapMyCells annotation reuses the existing local Allen WMB 20240831 mapping after exact query-count, stable-gene-ID and nucleus checks across all 12 samples. These labels define common cohorts and accompany fresh Scanpy and Seurat representations. The corrected branches use the frozen original-query labels. Annotation is therefore shared evidence, rather than an independent classifier agreement test.

- [Main analysis report](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/run_analysis_summary.md)
- [Results navigation](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/RESULTS_INDEX.md)
- [Current status](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/RUN_STATUS.md)
- [CellBender review](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/audit_cellbender/CELLBENDER_AUDIT_REVIEW.md)
- [Input/QC review](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/reviews/input_qc_review.md)
- [Scanpy output review](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/reviews/scanpy_result_review.md)
- [Seurat output review](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/reviews/seurat_result_review.md)
- [Final cross-tool comparison](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/reviews/comparison_result_review.md)
- [Readable completed corrected-count Seurat coordinates](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/reviews/figures/seurat_v1_0_5_paired_cellbender_visible_coordinates.pdf)
- [Active results](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/runs/20261013_v1_0_5/results)
- [Active stage logs](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/runs/20261013_v1_0_5/logs)
- [Exact run configuration](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/config/analysis_v1_0_5.json)

The native registered run is 9702cbef-00e6-4f64-aa37-e9b62357bd71. Its immutable workflow source is version-88efd06f4ac5ec0eeb88ec3c8ee4d1e6, SHA-256 acac5a94719db38829bba09b15acad85fac6bcaadf282b0ae05ba0154fcec062. Earlier failed/canceled task runs are retained. Verified byte-identical duplicates within those task-created runs use APFS filesystem clones; the storage receipt records unchanged content hashes, sizes and modification times.

The Seurat RDS files contain the fixed variable-feature representation; full-gene counts remain available in the shared H5ADs and are independently imported for full-gene marker and mouse-level pseudobulk exports. Genotype differential expression uses integer counts with the mouse as the replicate, rather than integrated values or nuclei as independent replicates.

Transparent APFS compression is applied to task-owned closed text outputs, gated by completed producers or finalized endpoint audit publication. Their paths, logical contents, hashes, permissions and modification times are preserved; reviews/task_text_compression_receipt.jsonl records the storage operations. No previous user analysis or scientific source/configuration is changed.
