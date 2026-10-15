The 14 October module defines a DESeq2-only revision. It preserves the previous run, uses original Harmony resolution0.4 clusters, aggregates every included nucleus, and exports exactly `gene`, `baseMean`, `log2FoldChange`, `lfcSE`, `stat`, `pvalue`, `padj`. `gene` contains the source gene symbol; duplicate symbols remain separate stable-ID modeled rows and missing symbols remain NA. Stable-ID identity and row correspondence are explicit in AUDIT.

Native CLI (strict argument order):

```
Rscript --vanilla run_deseq2_revision.R --config MASTER.json --stage prepare --branch original_full --outdir STAGE_DIRECTORY
```

Stages: `prepare`, `finalize_populations`, `deseq2`. Branches: `original_full`, `paired_original`, `paired_cellbender`. Output directories may exist only if genuinely empty. A failed directory is preserved and requires a new run root. No dependency installation occurs. Source/parse inspection does not execute the workflow.

Master JSON fields consumed:

| Field | Value |
| --- | --- |
| `previous_run_root` | Completed old native run root |
| `study_root` | Old study project root (legacy core fallback only) |
| `plot_gene_panel_tsv` | Initial annotation `gene_panel.tsv`, unique `gene_id`, `gene_symbol` |
| `deseq2.sample_manifest_tsv_by_branch` | Map branch→TSV with all12 `sample`, `condition`, `sequencing_partition`, `filtered_h5`, `source_h5_sha256` |
| `deseq2.roster_tsv_by_branch` | Map branch→frozen count-roster TSV with `nucleus_id`, `barcode`, `sample`, `condition`, `nCount_RNA`; corrected branch also `cellbender_nCount_RNA` |
| `deseq2.gene_manifest_tsv` | Source-order TSV `gene_id`, `gene_symbol` |
| `deseq2.population_index_tsv_by_branch` | Map branch→initial atomic index; one combined file may serve all branches |
| `deseq2.population_memberships_tsv_by_branch` | Map branch→initial atomic membership TSV(.gz); combined file is filtered by `branch` |
| `deseq2.prepared_root_by_branch` | Map branch→`results/prepared/BRANCH` |
| `deseq2.category_definition_tsv` | Final annotation `category_definition.tsv` |
| `deseq2.category_membership_tsv` | Final annotation `category_memberships.tsv.gz`, exact union check |
| `deseq2.finalized_root_by_branch` | Map branch→`results/populations/BRANCH` |
| `deseq2.count_acceptance_json_by_branch` | Map branch→independent `results/count_validation/COMPLETE.json`; may be one all-branch receipt |
| `deseq2.count_acceptance_bool_key` | `exact_integer_counts_PASS` (default); must be literal true before any model |
| `deseq2.prefilter_universe_root_by_branch` | Corrected branch obligatorily maps to paired-original **model output root**. Original-full/paired-original derive their own original-count primary universe; omit those keys. |
| `deseq2.cache_pseudobulk_root_by_branch` | Optional old pseudobulk root overrides; fallback `previous_run_root/results/seurat/BRANCH/pseudobulk` |
| `deseq2.legacy_statistical_core` | Audited old `source/seurat/pseudobulk_core.R` with recognized hash |
| `deseq2.memory_probe_python` | Installed psutil Python; top-level `memory_probe_python` fallback |
| `deseq2.minimum_available_gib` | Stage-start guard; default8GiB |
| `deseq2.input_hashes_tsv` | Optional frozen source-input `path`, `sha256` table; every listed input is verified |
| `deseq2.run_original_core_sensitivities` | true(default): known-partition/main-partition and12 leave-one-mouse-out endpoints only for original-full core totals |

Atomic index required fields: `population_key`, `tool`, `population`, `cluster_id`, `n_nuclei`, `membership_sha256`, `model_enabled`, and `branch` if combined. Membership required fields: `population_key`, `nucleus_id`, plus `branch` if combined. IDs use the authoritative underscore sample-barcode identity. The initial annotation owns the iterative common All/EX masks; preparation never changes them. All original cluster levels, including emptied/excluded groups, remain in the registry. Empty membership hashes use SHA256(empty bytes); nonempty hashes use sorted UTF-8 IDs, one LF after each ID.

Optional atomic reuse fields: `source_counts_file`, `source_counts_sha256`, `source_membership_sha256`, `source_count_branch`, `cache_population_key`. Exact source membership and branch match are required for count reuse; otherwise the native capture pass aggregates the current membership. Any supplied cached file hash must match before use. `cache_population_key` identifies a possible old DESeq2 endpoint; it never authorizes model reuse alone.

Preparation produces `POPULATION_INDEX.tsv`, `POPULATION_MEMBERSHIPS.tsv.gz`, `GENE_ID_SYMBOL_MANIFEST.tsv`; each atom has `AUDIT/SOURCE_INTEGER_COUNTS_ALL_GENES_ALL12.tsv.gz` and `AUDIT/ALL12_POPULATION_COVERAGE.tsv`. `cluster_mouse_counts.tsv` and `gene_cluster_mouse_aggregates.tsv.gz` enumerate every original atom×all12 mice (including zeros). Panel sums are `sum_counts`, `n_detected`, `sum_log1p_normalized`, where each nucleus is normalized using its actual source full-gene library to10000 before log1p. No marker sampling enters these summaries or pseudobulk.

Finalization sums disjoint atomic arrays, builds each tool/population total and cluster endpoints, and adds reviewed category/program unions. Annotation columns: `tool`, `branch`, `population`, `category_family`, `label`, JSON `member_cluster_ids`, `cluster_scope=original_frozen`, `source_branch=original_full`. Supplied category memberships must exactly match those unions. Category families cannot assign an atom to two categories. It writes `POPULATION_INDEX.tsv`, `POPULATION_ATOMS.tsv`, `POPULATION_MEMBERSHIPS.tsv.gz`, per-population AUDIT count/coverage arrays, and mask registry. Raw and corrected paired labels/masks remain fixed.

DESeq2 outputs `TOOL/POPULATION_KEY/PRIMARY/{ALL_RESULTS.txt,DEG_FDR05.txt,AUDIT/...}`; sensitivity endpoints use analogous names. No result files are fabricated for unsupported models. `PSEUDOBULK_EXECUTION_AUDIT.tsv` and `POPULATION_INDEX.tsv` provide navigation. An endpoint is published only when `AUDIT/EXECUTION_AUDIT.tsv` exists; a stage is published only by `COMPLETE.json`. No repeated normalized-all-gene-by12 tables are written.

IN includes every nonempty positive-library mouse aggregate (`min_nuclei=1`), gene counts≥5 in≥3 eligible mice; other populations require≥20 nuclei/mouse and gene counts≥10 in at least min(eligibleWT,eligiblecKO) mice. All models require≥3 independent mice/genotype, full rank and≥2 residual df. Included IN mice with<5 nuclei cause `LOW_NUCLEUS_COVERAGE_EXPLORATORY`, with min/max coverage and full all12 mouse audits. Original gene universes are frozen before sensitivities; corrected uses paired-original universe. Default Cook NAs and unavailable FDR values remain NA. BH is within each endpoint; overlapping categories are exploratory and do not share a global correction.

The contrast is cKO/WT; normalization poscounts, Wald test, parametric dispersion request, betaPrior false, count replacement disabled, sequential fits. Leave-one-mouse-out disables independent filtering for influence evidence. Known-partition comparisons exclude unverified cKO933; lane8 has only cKO, so the genotype coefficient is chiefly supported by the shared main partition. These sensitivities do not justify discarding mice.

Cache reuse requires all12 integer counts, primary and realized gene universes, selected mice/partitions/design, contrast/result policies, reviewed source code and matching runtime versions. Reuse is labeled explicitly with source hashes; it is never described as a fresh fit. Same-run identical endpoints may reuse their exact count/model evidence. A failed cache guard causes a native fit. Native count validation is mandatory first.

`test_revision_io.R` uses fabricated four-gene H5s and does not read study inputs or fit any model. It verifies explicit CSC import including latent root groups, all12 conservation, full-library panel arithmetic, IN1 metadata, exact symbol exports, disjoint unions and duplicate-atom rejection.
