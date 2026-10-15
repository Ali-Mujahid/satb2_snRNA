# Results navigation

This index links task-owned outputs. Producer completion and independent scientific acceptance are separate; consult the main report for the current lifecycle and conclusions.

[Main analysis report](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/run_analysis_summary.md)

## Representations

The three count branches use the original full cohort (108,062 nuclei), paired original (104,298) and paired CellBender (104,298). Scanpy objects retain full-gene count layers. Seurat RDS files store the variable-feature representation; use full-source QC sidecars and full-gene inference exports for library/QC/count interpretation. SCT PCA uses 2,731 modeled genes.

| Count branch | Cohort | Scanpy | Seurat v5 |
| --- | --- | --- | --- |
| original_full | allcells | [H5AD + tables + plots](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/runs/20261013_v1_0_5/results/scanpy/original_full/allcells) | [HVG RDS + tables + plots](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/runs/20261013_v1_0_5/results/seurat/original_full/log_global) |
| original_full | Excitatory | [H5AD + tables + plots](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/runs/20261013_v1_0_5/results/scanpy/original_full/Excitatory) | [HVG RDS + tables + plots](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/runs/20261013_v1_0_5/results/seurat/original_full/annotate_subsets/EX) |
| original_full | Inhibitory | [H5AD + tables + plots](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/runs/20261013_v1_0_5/results/scanpy/original_full/Inhibitory) | [HVG RDS + tables + plots](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/runs/20261013_v1_0_5/results/seurat/original_full/annotate_subsets/IN) |
| original_full | CA1_ProS | [H5AD + tables + plots](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/runs/20261013_v1_0_5/results/scanpy/original_full/CA1_ProS) | [HVG RDS + tables + plots](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/runs/20261013_v1_0_5/results/seurat/original_full/annotate_subsets/CA1_PROS_CANDIDATE) |
| paired_original | allcells | [H5AD + tables + plots](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/runs/20261013_v1_0_5/results/scanpy/paired_original/allcells) | [HVG RDS + tables + plots](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/runs/20261013_v1_0_5/results/seurat/paired_original/log_global) |
| paired_original | Excitatory | [H5AD + tables + plots](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/runs/20261013_v1_0_5/results/scanpy/paired_original/Excitatory) | [HVG RDS + tables + plots](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/runs/20261013_v1_0_5/results/seurat/paired_original/annotate_subsets/EX) |
| paired_original | Inhibitory | [H5AD + tables + plots](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/runs/20261013_v1_0_5/results/scanpy/paired_original/Inhibitory) | [HVG RDS + tables + plots](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/runs/20261013_v1_0_5/results/seurat/paired_original/annotate_subsets/IN) |
| paired_original | CA1_ProS | [H5AD + tables + plots](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/runs/20261013_v1_0_5/results/scanpy/paired_original/CA1_ProS) | [HVG RDS + tables + plots](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/runs/20261013_v1_0_5/results/seurat/paired_original/annotate_subsets/CA1_PROS_CANDIDATE) |
| paired_cellbender | allcells | [H5AD + tables + plots](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/runs/20261013_v1_0_5/results/scanpy/paired_cellbender/allcells) | [HVG RDS + tables + plots](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/runs/20261013_v1_0_5/results/seurat/paired_cellbender/log_global) |
| paired_cellbender | Excitatory | [H5AD + tables + plots](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/runs/20261013_v1_0_5/results/scanpy/paired_cellbender/Excitatory) | [HVG RDS + tables + plots](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/runs/20261013_v1_0_5/results/seurat/paired_cellbender/annotate_subsets/EX) |
| paired_cellbender | Inhibitory | [H5AD + tables + plots](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/runs/20261013_v1_0_5/results/scanpy/paired_cellbender/Inhibitory) | [HVG RDS + tables + plots](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/runs/20261013_v1_0_5/results/seurat/paired_cellbender/annotate_subsets/IN) |
| paired_cellbender | CA1_ProS | [H5AD + tables + plots](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/runs/20261013_v1_0_5/results/scanpy/paired_cellbender/CA1_ProS) | [HVG RDS + tables + plots](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/runs/20261013_v1_0_5/results/seurat/paired_cellbender/annotate_subsets/CA1_PROS_CANDIDATE) |

[Original global SCT representation](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/runs/20261013_v1_0_5/results/seurat/original_full/sct_global)

## Full-gene markers and mouse-level inference

Marker tables describe cluster identity; genotype inference uses mouse-level integer pseudobulk with six WT and six cKO animals. Eligibility and model/filtering audits are retained for small or untested populations. Paired R marker/DE groups use frozen original memberships restricted to the paired roster, while paired representations above are reclustered. Numeric cluster IDs must not be exchanged across those scopes.

| Count branch | Full-gene Seurat markers | Population registry and model tree | Independent Python counts |
| --- | --- | --- | --- |
| original_full | [Markers](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/runs/20261013_v1_0_5/results/seurat/original_full/markers) | [Population index](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/runs/20261013_v1_0_5/results/seurat/original_full/pseudobulk/POPULATION_INDEX.tsv); [Models/audits](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/runs/20261013_v1_0_5/results/seurat/original_full/pseudobulk) | [Integer gene-by-mouse exports](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/runs/20261013_v1_0_5/results/python_pseudobulk/original_full) |
| paired_original | [Markers](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/runs/20261013_v1_0_5/results/seurat/paired_original/markers) | [Population index](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/runs/20261013_v1_0_5/results/seurat/paired_original/pseudobulk/POPULATION_INDEX.tsv); [Models/audits](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/runs/20261013_v1_0_5/results/seurat/paired_original/pseudobulk) | [Integer gene-by-mouse exports](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/runs/20261013_v1_0_5/results/python_pseudobulk/paired_original) |
| paired_cellbender | [Markers](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/runs/20261013_v1_0_5/results/seurat/paired_cellbender/markers) | [Population index](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/runs/20261013_v1_0_5/results/seurat/paired_cellbender/pseudobulk/POPULATION_INDEX.tsv); [Models/audits](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/runs/20261013_v1_0_5/results/seurat/paired_cellbender/pseudobulk) | [Integer gene-by-mouse exports](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/runs/20261013_v1_0_5/results/python_pseudobulk/paired_cellbender) |

Four broad population keys are P0001_ATLAS_TOTAL_ALL, P0002_ATLAS_TOTAL_EX, P0003_ATLAS_TOTAL_IN and P0004_ATLAS_TOTAL_CA1_PROS. In each model folder, ALL_RESULTS.txt and DEG_FDR05.txt are standard gene-symbol exports; the corresponding AUDIT folder carries stable IDs, full mouse coverage, model settings and preserved NA/eligibility information. Condition effects are cKO/WT. Python independently exports the four atlas totals and original Scanpy Harmony populations; Seurat-defined cluster aggregates are R-only.

## Reviews and readable figures

- [Input and QC review](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/reviews/input_qc_review.md)
- [CellBender audit](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/audit_cellbender/CELLBENDER_AUDIT_REVIEW.md)
- [Scanpy review](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/reviews/scanpy_result_review.md)
- [Seurat review](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/reviews/seurat_result_review.md)
- [Final comparison review](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/reviews/comparison_result_review.md)
- [All primary model navigation](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/reviews/pseudobulk_endpoint_navigation_20261013_v1_0_5_v1/PRIMARY_ENDPOINT_NAVIGATION.tsv)
- [Core same-roster count sensitivity](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/reviews/pseudobulk_first4_atlas_count_sensitivity_v1/effect_comparisons.tsv)
- [Supported frozen-population cross-tool effects](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/reviews/paired_frozen_cluster_DE_v1_0_5/effect_metrics.tsv)
- [Completed Cck strata](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/reviews/Cck_completed_strata_v1_0_5/Cck_stratum_summary.tsv)
- [Seurat marker label and sampling scope](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/reviews/seurat_20261013_v1_0_5_MARKER_SCOPE_AND_SOURCE_LABEL_MANIFEST.tsv)
- [Full-gene marker review](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/reviews/seurat_completed_marker_output_review_v1)
- [Cluster-agreement heatmaps](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/reviews/final_cluster_agreement_heatmaps_v1_0_5/cluster_agreement_heatmaps.pdf)
- [Same-nucleus CellBender effect figure](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/reviews/figures/paired_CellBender_DESeq2_effect_sensitivity_v1.pdf)
- [Original DESeq2 versus edgeR effect figure](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/reviews/figures/pseudobulk_first4_atlas_DESeq2_edgeR_effect_comparison_v1.pdf)
- [Full-source Seurat QC sidecar index](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/reviews/metadata/20261013_v1_0_5__QC_SIDECAR_INDEX.tsv)
- [Original overview and mouse facets](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/reviews/figures/seurat_v1_0_5_original_full_visible_coordinates.pdf)
- [Original CA1/ProS overview and mouse facets](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/reviews/figures/seurat_v1_0_5_original_full_CA1_PROS_visible_coordinates.pdf)
- [Matched CellBender overview and mouse facets](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/reviews/figures/seurat_v1_0_5_paired_cellbender_visible_coordinates.pdf)
- [Original SCT overview and mouse facets](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/reviews/figures/seurat_v1_0_5_original_full_SCT_visible_coordinates.pdf)

## Reproducibility

- [Frozen analysis configuration](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/config/analysis_v1_0_5.json)
- [Scientific source](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/source)
- [Native stage logs](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/runs/20261013_v1_0_5/logs)
- [Python package lock](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/scanpy_requirements_lock.txt)

Registry: 9702cbef-00e6-4f64-aa37-e9b62357bd71. Workflow snapshot: version-88efd06f4ac5ec0eeb88ec3c8ee4d1e6. Previous user analyses are preserved.
