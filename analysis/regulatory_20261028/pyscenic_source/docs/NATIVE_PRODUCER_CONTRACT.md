# Native pySCENIC producer contract, v1.0.0

This is an additive scientific implementation for the user's original selection of full pySCENIC, curated TF activity, and hdWGCNA. The official mouse reference service became reachable again on 28 October 2026. Download receipts must establish the actual matching files before execution. Earlier blocked-reference records remain historical evidence; they are not silently replaced. Only the root agent launches an immutable NGS Workbench Snakemake plan. Importing or inspecting this source does not execute study analysis.

## Inputs and immutable configuration

`scripts/run_pyscenic_native.py` consumes `--config /absolute/config.json --outdir /absolute/results/pyscenic_native --stage STAGE`, with `--lineage CA1|IN` except for `assemble`.

Required JSON keys:

```json
{
  "analysis_id": "a fresh immutable analysis identity",
  "run_root": "/absolute/new/immutable/run",
  "samples": {"path": "/absolute/frozen/twelve_mouse_samples.tsv", "bytes": 1, "sha256": "actual digest"},
  "prepared_complete": {"path": "/absolute/native/prepared/COMPLETE.json", "bytes": 1, "sha256": "actual digest"},
  "lineages": {
    "CA1": {"counts": {"path": "/absolute/counts.mtx.gz", "bytes": 1, "sha256": "actual"}, "genes": {"path": "/absolute/genes.tsv", "bytes": 1, "sha256": "actual"}, "cells": {"path": "/absolute/cells.tsv", "bytes": 1, "sha256": "actual"}},
    "IN": {"counts": {"path": "/absolute/counts.mtx.gz", "bytes": 1, "sha256": "actual"}, "genes": {"path": "/absolute/genes.tsv", "bytes": 1, "sha256": "actual"}, "cells": {"path": "/absolute/cells.tsv", "bytes": 1, "sha256": "actual"}}
  },
  "references": {
    "tf_names": {"path": "/absolute/allTFs_mm.txt", "bytes": 1, "sha256": "actual"},
    "motif_annotations": {"path": "/absolute/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl", "bytes": 1, "sha256": "actual"},
    "rankings": [
      {"path": "/absolute/mm10_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather", "name": "mm10_v10_500bp", "bytes": 1, "sha256": "actual"},
      {"path": "/absolute/mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather", "name": "mm10_v10_10kbp", "bytes": 1, "sha256": "actual"}
    ]
  },
  "additional_hash_inputs": {"frozen_membership_audit": {"path": "/absolute/audit.tsv", "bytes": 1, "sha256": "actual"}, "mouse_coverage": {"path": "/absolute/coverage.tsv", "bytes": 1, "sha256": "actual"}}
}
```

Example digest/size fields above are illustrative schema placeholders, never valid execution inputs. Root supplies actual verified receipts. The source verifies every input before and after every stage. Prepared inputs must match the SHA256/byte records in the completed original hdWGCNA preparation and its all-gene pseudobulk reconciliation. Exact fixed endpoint membership is 8,750 deep-program, 11,706 superficial-program, 1,718 C1, 501 C6, and 468 C8 nuclei. Every endpoint must contain all twelve named independent mice, with exact frozen genotype assignments. Counts are gene rows × nucleus columns, nonnegative integers, with full-library sums equal to saved nucleus metadata. Missing or ambiguous gene symbols are excluded, not guessed.

Fixed scalar parameters may appear in config but cannot differ from the protocol: `seed=1337`, `num_workers=4`, `min_detection_fraction=0.05`, `normalization_cp10k=10000`, `module_min_genes=20`, `rho_threshold=0.03`, `rank_threshold=1500`, `ctx_auc_threshold=0.05`, `nes_threshold=3.0`, `motif_similarity_fdr=0.001`, `orthologous_identity_threshold=0.0`, `min_regulon_targets=10`, `aucell_threshold=0.05`, `aucell_batch_size=1000`, `aucell_num_workers=1`.

## Stages and scientific semantics

1. **prepare** creates `<lineage>/prepared/`. Use full RNA with positive detection in at least 5% of nuclei of any exact endpoint in that lineage, without using genotype or DEGs/HVGs for selection. Normalize `log1p(10000 × count / original full-library UMI sum)` and exclude constant normalized profiles. Preserve a row for every original prepared feature in the audit. Retain only actually eligible genes from the unchanged official TF candidate list; record SATB2's actual list membership and measured RNA eligibility, with no forced addition. Store sparse CSC expression only. All nuclei remain included; no artificial nuclei or discovery downsampling.
2. **grn** creates `<lineage>/grn/` through native `arboreto.algo.grnboost2` with the pinned GRNBoost2 profile and seed. Sparse CSC and a four-worker threaded `LocalCluster` avoid large per-process dense copies. Save native adjacency weights without scaling or arbitrary edge truncation. Native importance is nonnegative and is not a regulatory sign or genotype effect. Every nonconstant target must return native inference links; absent target output is an explicit workflow failure for review, not a silently imputed zero.
3. **ctx** creates `<lineage>/ctx/`. Native `add_correlation` uses Pearson correlation on log-normalized RNA with dropout masking disabled. Native `modules_from_adjacencies` uses 75th/90th importance percentiles, top50 targets, top5/10/50 regulators, minimum20 genes, and positive correlation >0.03. These are positive coexpression modules, not proven activating mechanisms. Candidate modules, full correlations, and per-ranking-database mapping coverage are preserved. Native `prune2df` uses both exact mm10 gene-based ranking resources, NES≥3, recovery rank1500, AUC0.05, matching MGI annotations, motif similarity FDR≤0.001, orthology identity≥0, unweighted motif recovery, and annotation filtering. The native module mapping gate requires more than80% of module genes mapped. Native `df2regulons` unions the supported leading-edge targets per TF and positive context, retaining maximum native GRNBoost2 importance per target across surviving modules. Save the complete native multi-index dataframe, native CSV, native modules, and native regulons, plus labeled TSV adaptations. Preserve every returned native positive regulon; score eligibility requires at least10 measured retained targets.
4. **aucell** creates `<lineage>/aucell/`. Run native pySCENIC AUCell on each exact nucleus, with the same full eligible RNA background and gene order used for discovery. Use native target weights, top5% of ranked genes, seed1337 to break expression ties, and `normalize=False`; no observed-max normalization or AUC threshold/binarization. Batches of1,000 nuclei call the serial native kernel to propagate exceptions and avoid silently lost multiprocessing children. Native serial output sorts nucleus labels; the adapter verifies exact unique label membership and restores original label order without changing values. Synthetic native tests establish exact whole-versus-batched equality with zeros and ties. Save exact nucleus AUC, metadata axes, coverage, and score-integrity summaries; flag constant/all-zero native score columns rather than manufacture variability. Means are one arithmetic mean per exact endpoint × mouse × regulon using all included nuclei; no nucleus is a genotype replicate.
5. **assemble** runs after all eight upstream stage markers. Verify all source/config/output bytes, combine tables, and write `LINEAGE_STATUS.tsv`, `PROVENANCE.json`, and final native `COMPLETE.json`. A lineage with no supported regulons has `EMPTY_NO_SUPPORTED_REGULONS`; an empty finding is a completed computation, never fabricated zero AUC. A lineage with returned regulons has `COMPLETE_WITH_REGULONS` even if low target coverage prevents some/all score tests.

All stages refuse nonempty output directories or existing file overwrites; strictly empty directories created by Snakemake are allowed. Completion markers are first pending and are only promoted after unchanged-input/source checks pass. A failed partial stage must be recovered through a new immutable plan, not edited in place. Every fitted lineage rule reserves four Snakemake cores so lineages/stages do not overrun the machine simultaneously. Set OMP/BLAS/MKL/NUMBA thread counts to1 externally to avoid nested thread oversubscription. The maximum dense correlation step is approximately nuclei × eligible genes plus genes², held in memory rather than persisted as dense expression text. References and prepared inputs are reused read-only.

## Table and axis outputs for downstream R/packager

Final native marker schema: `satb2-pyscenic-native-complete/v1`; status `COMPLETE`; `config` and single producer `source` receipts; `output_files` maps relative file path to absolute `{path,bytes,sha256}`. Required checks explicitly state positive motif regulons, full genotype-neutral RNA, no DEG selection, actual nucleus AUCell, exact axes, all12 roster, min10 targets, separate SATB2 sensitivity, and unchanged hashes.

- `REGULON_ACTIVITY_PER_MOUSE.tsv`: lineage, endpoint, endpoint_id, source_endpoint_id, regulon_id, native_name, TF, n_targets, n_background_genes, sample, condition, n_nuclei, aucell_mean, score_unit=`native_AUCell_AUC`, eligible, network_variant.
- `REGULON_COVERAGE.tsv`: same endpoint/regulon identifiers, n_targets, n_background_genes, n_native_targets, n_self_targets, eligible, reason, network_variant. All returned regulons appear even when fewer than10 targets prevents scores.
- `REGULON_TARGETS.tsv`: lineage, regulon_id, native_name, TF, target, weight, self_target, edge_origin, edge_basis.
- `LINEAGE_STATUS.tsv`: lineage, status, n_regulons, n_score_eligible_regulons, n_nuclei, n_genes; CA1 and IN exactly once.
- `<lineage>/aucell/nucleus_AUC.tsv.gz`: first column `nucleus_id`, followed by lineage-distinct regulon IDs such as `CA1::Satb2(+)`. All exact nucleus rows occur in original prepared order. Empty-regulon lineages retain nucleus IDs with no fake activity columns.
- `<lineage>/aucell/nucleus_axis_metadata.tsv.gz`: nucleus_id, sample, condition, endpoint, lineage, sequencing_partition; full original labeled metadata is also retained at `<lineage>/prepared/nucleus_metadata.tsv.gz`.
- `PROVENANCE.json` supplies exact paths/receipts under `lineages.CA1|IN.nucleus_AUC`, `nucleus_axis_metadata`, `prepared_metadata`, `native_regulons`, and `native_motif_enrichment`.

The unique display ID is `<lineage>::<native regulon name>`, with `native_name` preserved verbatim. The two lineage regulons sharing a TF name are not silently merged. Gene coverage/backgrounds may differ; absolute AUC levels across lineages are not calibrated measurements of TF activity.

## Native self gene and conditional SATB2 sensitivity

Native pySCENIC adds the TF's own RNA gene to candidate modules. If retained after motif pruning it remains in the untouched native regulon, and its exported `self_target=True`/`edge_origin=native_module_added_TF_self_gene` explicitly prevents treating it as an inferred binding edge. Native importance weights are not binding strengths.

Only if an actual positive motif-supported Satb2 regulon exists, contains Satb2 itself, and retains at least10 genes after removing this self gene, run a **separate** native AUCell sensitivity using precisely the same rank background and seed. Do not modify the original regulon. `SATB2_NO_SELF_TARGET_ACTIVITY_PER_MOUSE.tsv` and `SATB2_NO_SELF_TARGET_COVERAGE.tsv` use the main schemas, with ID `<lineage>::Satb2(+)_NO_SELF_TARGET`, original `native_name=Satb2(+)`, and separate network_variant. Coverage always includes all five exact endpoints; missing regulon, absent self gene, or insufficient remaining targets is explicit with no fabricated scores. The sensitivity remains separate from the main global BH family and is not independent validation. If produced, actual nucleus sensitivity AUC is retained separately.

## Interpretation boundary

Both genotypes jointly discover lineage-specific regulons and all frozen nuclei are scored; downstream mouse comparisons are exploratory inference on this same cohort. Subtype composition, technical effects, and indirect coexpression can drive edges. Motif resource support is more specific than coexpression alone but does not measure chromatin accessibility, occupancy, protein activity, physical contacts, direction of regulation, or SATB2 causality. CA1 layer labels remain expression programs; C8 remains mixed identity. Partial genotype–sequencing-partition confounding remains. Original DESeq2 gene analysis, curated-prior TF activity, hdWGCNA, strict DEG-only ligand–receptor results and earlier files are preserved.

Official method references: [pySCENIC tutorial](https://pyscenic.readthedocs.io/en/latest/tutorial.html), [official pySCENIC source](https://github.com/aertslab/pySCENIC), and [Arboreto source](https://github.com/aertslab/arboreto). The executed version is the already pinned pySCENIC0.12.1/Arboreto0.1.6/ctxcore0.2.0 environment, not the moving documentation version.
