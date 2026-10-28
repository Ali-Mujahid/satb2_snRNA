# Configuration contract

Supply the existing absolute `config/hdwgcna.v1.json` as `--configfile` through the Workbench plan. It contains `run_root`, `python_executable`, `Rscript_executable`, `threads=4`, `allow_network=false`, `overwrite_prior=false` and the `satb2-regulatory-preparation/v1` input schema. Every input is an absolute path with exact byte count and SHA-256. The workflow binds frozen counts, Seurat annotations/memberships, all five native raw pseudobulks, all five original ALL_RESULTS files, the corrected strict LR table and the isolated R runtime manifest.

The `hdwgcna` object uses the exact native implementation keys in `run_hdwgcna.R::defaults`. Primary settings: seed1337, threads4, expression_fraction0.05, signed bicor, max_p_outliers0.05, powers1–10+12–30 by2, topology_fit_min0.8, negative slope, mean_connectivity_min5, min_module_size50, merge_cut_height0.2, deep_split4, max_block_size3000. Neighborhood PCA uses 2000 HVGs; the network gene background uses full eligible RNA expression. CA1 and IN metacell parameters are explicit under `CA1` and `IN`. `execute=true` belongs to the immutable approval-bound config; it must never be flipped after approval. Consensus is separately labeled and conditional on coverage/topology gates.

The same-partition sensitivity must use the observed exact value `23HJT3LT4_lanes6_7`. It has six WT and three cKO and is exploratory. Keep the unverified cKO933 partition explicit.

No installation or reference download can be embedded here. Optional whole-transcriptome H5AD writes are false for this hdWGCNA-only run to reduce duplicated storage. Output extraction/normalization and all fitting happen only after native approval.
