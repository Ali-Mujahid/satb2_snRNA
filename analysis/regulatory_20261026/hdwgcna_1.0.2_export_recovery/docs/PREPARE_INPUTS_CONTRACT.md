# Frozen whole-transcriptome preparation draft

Status: source drafted and syntax compiled only. No study preparation, count transformation, normalization, analysis, installation, download or model execution has been performed by this draft. The registered, approved Snakemake workflow must own scientific execution and output paths.

## Entry point and inputs

`source/scripts/prepare_inputs.py --config ABS_CONFIG_JSON --outdir ABS_RUN_RESULTS_PREPARED`

The output must be a new/empty directory strictly below the configured absolute `run_root`. A nonempty output is rejected. Input files are opened read-only, verified against configured SHA256/size before use and again after preparation. All generated files and their hashes are recorded in `COMPLETE.json`.

Configuration schema: `satb2-regulatory-preparation/v1`. The actual native input paths, byte sizes and freshly verified hashes are in `config/preparation.inputs.example.json`. Its `run_root` contains a `REPLACE_WITH_REGISTERED_RUN_ROOT` placeholder and must be set by the registered workflow. The example is not a record of scientific execution.

Required `inputs` keys, each with `path`, `sha256`, and optional `bytes`:

- `counts`: shared `original_full/branch_counts.h5ad`, sparse raw integer X with 108,062 nuclei × 32,285 stable-ID features. `obs/nucleus_id`, `obs/sample`, `obs/condition`, `var/gene_id`, `var/gene_symbol` are read using the actual AnnData/HDF5 encoding.
- `nucleus_annotations`: final annotation table. Select `tool=seurat`, `branch=original_full`, `cluster_scope=original_frozen`, the exact endpoint label/cluster, and the existing native `keep_for_analysis` flag.
- `category_memberships`: independently checks exact CA1 program labels and the inhibitory broad parent for the selected nuclei.
- `population_index`: native finalized original-full endpoint definitions, membership hashes, and source raw-pseudobulk paths/hashes.
- `population_memberships`: exact finalized nucleus IDs for each original-full endpoint.
- `samples`: authoritative 12-mouse roster with sample, genotype and sequencing partition.
- `endpoint_counts_CA1_Deep`, `endpoint_counts_CA1_Superficial`, `endpoint_counts_IN_C1`, `endpoint_counts_IN_C6`, `endpoint_counts_IN_C8`: all five native all-gene/all-mouse integer pseudobulk files are explicit workflow readiness dependencies. Their paths and hashes must also agree exactly with the independently hash-bound finalized population index.

Optional `population_atoms` and `atomic_population_memberships` must be supplied together. The script checks declared atomic unions using only `tool=seurat`, `branch=original_full`, because atomic population keys repeat across branches. Other configured hash-bound inputs are preserved in provenance, but are not implicitly used to change memberships.

The finalized population index independently authorizes the same five native `SOURCE_INTEGER_COUNTS_ALL_GENES_ALL12.tsv.gz` files. Their paths and native hashes must match the explicit workflow input records and file contents. Every gene-ID × mouse integer count recomputed from the exact raw selected nuclei must equal its saved original pseudobulk count. No DESeq2 model is fitted by this check.

## Exact scope and corrected count

| Analysis endpoint | Native key | Exact frozen clusters | Nuclei | Lineage |
|---|---|---|---:|---|
| CA1_Deep | P0244_seurat_CA1_ProS_LAYER_PROGRAM | 1;5 | 8,750 | CA1 |
| CA1_Superficial | P0247_seurat_CA1_ProS_LAYER_PROGRAM | 2;3;8 | 11,706 | CA1 |
| IN_C1 | P0140_seurat_Inhibitory_CLUSTER | 1 | 1,718 | IN |
| IN_C6 | P0145_seurat_Inhibitory_CLUSTER | 6 | 501 | IN |
| IN_C8 | P0147_seurat_Inhibitory_CLUSTER | 8 | 468 | IN |

CA1 total = 20,456; IN total = 2,687; combined = 23,143 nuclei. The requested draft count of 568 for C8 was corrected to the authoritative 468, consistent with the original population index, all 12 mouse coverage rows and completed native communication provenance. The script cannot add artificial nuclei or silently accept another count. The configured expected-count dictionary must match these exact values.

`frozen_membership_audit.tsv` records scoped candidate counts before the native keep flag, retained counts, exclusions under the existing flag, selected/native membership hashes, raw-pseudobulk paths and hashes. `frozen_candidates_before_after_keep.tsv.gz` retains the candidate-level original flags and exclusion reasons. No new QC threshold or identity filter is introduced. C8's mixed-identity flag is retained; it is not upgraded or silently omitted.

## Feature and normalization policy

Raw stable-ID gene features sharing a real gene symbol are summed **before** normalization. Their original IDs and source positions remain in `feature_mapping.tsv`; duplicate symbols remain flagged in `duplicate_symbol_audit.tsv`. Gene IDs are never assigned a different biological symbol. Missing symbols receive separate reserved `__gene_id__<stableID>` analysis keys so unnamed genes cannot merge or collide with named genes. Downstream TF/network resources must exclude/handle these fallback keys explicitly rather than treat them as recognized symbols.

All genes are retained after duplicate-symbol collapse, including genes that are zero in a selected lineage. No DEG, highly variable gene, ligand–receptor or regulon target filter is applied during preparation. Original full-library integer UMIs must be conserved separately for every nucleus. This full RNA background is needed for regulatory inference; it does not loosen the separate DEG-only LR selection.

`write_normalized_h5ad` defaults to `false`. With this setting no normalized expression matrix is computed or written. hdWGCNA can read full raw MTX counts and apply its planned normalization using the complete collapsed library. `cells.tsv` preserves the original full-gene library and `10000 / library_umi_full_genes` scale factor.

If explicitly enabled in the registered configuration, `write_normalized_h5ad=true` additionally materializes log1p(CP10k) in each lineage H5AD layer using **all original 32,285 features** as the denominator before downstream filtering. X remains raw integer, duplicate-symbol-collapsed RNA counts; `layers['log1p_cp10k']` is float32 normalized expression. The script verifies that collapsing symbols did not change this denominator. No Harmony/SCT residual, imputation, integrated assay or UMAP coordinate enters expression preparation.

`write_raw_gene_id_h5ad` independently defaults to `false`; enabling it writes the sparse exact selected stable-ID raw matrix for auditing or methods requiring uncollapsed source features.

## Outputs and hdWGCNA interface

Always:

- `lineages/CA1/counts.mtx.gz` and `lineages/IN/counts.mtx.gz`: sparse Matrix Market **gene rows × nucleus columns**, raw integer UMIs after symbol collapse, no gene filtering.
- `lineages/{CA1,IN}/genes.tsv`: one row per MTX gene, exact order. Columns: `feature_id`, `gene_symbol`, `gene_ids`, `n_source_features`, `duplicate_symbol`, `symbol_missing`. `feature_id` and `gene_symbol` are unique analysis feature keys; `gene_ids` contains the original stable IDs joined by semicolons and is not guaranteed to represent one gene ID.
- `lineages/{CA1,IN}/cells.tsv`: one row per MTX column, exact order. Initial columns: `barcode`, `nucleus_id`, `sample_id`, `sample`, `condition`, `endpoint`, `lineage`, `sequencing_partition`, `lane`. Additional columns include original identity/review flags, `population_key`, `source_row_index`, `library_umi_full_genes`, `log1p_cp10k_scale_factor`. `barcode` is the complete sample-prefixed nucleus ID. `sample_id` aliases the original mouse ID. `lane` aliases the original sequencing-partition string, including `unverified` for cKO_933.
- `nucleus_metadata.tsv.gz`, `mouse_coverage.tsv`, `feature_mapping.tsv`, `duplicate_symbol_audit.tsv`, `frozen_membership_audit.tsv`, `frozen_candidates_before_after_keep.tsv.gz`, `raw_count_source_verification.tsv`, `expression_normalization.json`, `source_hashes.json`, `COMPLETE.json`.

Only when enabled:

- `lineages/{CA1,IN}/expression.h5ad`: X raw counts and normalized `log1p_cp10k` layer.
- `raw_gene_id_counts.h5ad`: selected raw stable-ID count matrix, uncollapsed, 23,143 × 32,285.

The R consumer should use `--input-dir ABS_PREP/lineages --lineages CA1,IN`. It must read exact MTX axes from `genes.tsv` and `cells.tsv`; use `feature_id`/`gene_symbol` as unique row keys and `barcode` as unique column keys. Preserve all original mouse and endpoint metadata. Network gene filtering, metacell policy, soft-power selection and TF testing belong to separately specified workflow stages.

## Scientific limits retained

Independent genotype replication is six WT and six cKO mice, not the number of nuclei or metacells. CA1 DL/SL denote molecular programs; no physical layer assignment is asserted. C6/C8 have sparse per-mouse coverage, and C8 retains mixed identity. Sequencing partitions are partially associated with genotype; cKO_933 partition is unverified. An inferred module or TF activity is not direct SATB2 binding or a causal regulatory mechanism.

## Validation status

Source syntax compilation passed. Following separate authorization for synthetic fixture validation, all six standard-library `unittest` fixtures in `tests/test_prepare_inputs.py` passed in the existing project Python environment. They cover only tiny synthetic feature-count arithmetic, safe blank/duplicate-symbol handling, exact fixed counts, membership hashing and mouse metadata aliases. No study expression or membership extraction was performed by those tests. Run with `python -B -m unittest discover -s tests -p test_prepare_inputs.py -v`; pytest is not required or installed in that environment.

Real input extraction, feature collapsing, normalizing, MTX writing, pseudobulk verification and frozen membership assertions remain **not executed** pending the registered approved workflow.
