# SATB2 analysis source and reference data

This directory preserves compatible v1/v2 scientific workflows, original configurations, executed workflow snapshots, actual runtime records, supplied marker workbooks, reviewed annotations, membership masks and provenance. It does not install software or automatically execute analysis.

Start with the [reproducibility guide](../docs/analysis/reproducibility.md) and [data delivery guide](../docs/analysis/data_delivery.md).

| Directory | Contents |
|---|---|
| `v1/source/` | Scanpy/Seurat representations, preparation, independent integer aggregation and original comparison |
| `v1/executed_workflow_snapshot/` | Exact source attached to completed run 20261013_v1_0_5 |
| `v1/config/`, `v1/environment/` | Original configuration; Python lock/freeze/wheel resolution; actual R sessions/package records |
| `v1/cellbender_audit/` | Source bindings and existing-correction audit source/results; no new correction fit |
| `v2/source/` | Revised annotation/filtering, native aggregation, DESeq2, plots and CA1 diagnostics |
| `v2/executed_workflow_snapshot/` | Exact source attached to completed run 20261014_v2_0_0 |
| `v2/config/`, `v2/environment/` | Original revision configuration and actual runtime/session records |
| `dashboard/`, `reports/source/` | Viewer/catalog/build/QA source, vendor licenses and report/figure builders |
| `data/` | Source-bound input/reference metadata, exact workbooks, panels, shared masks and category memberships |
| `reviews/`, `provenance/` | Completed reviews, nonbinding suggestions, copied-file hashes, receipts, indices and object release manifests |
| `portable/`, `tools/` | Explicit templates and a JSON-only path-materialization helper |

[The copied-file ledger](provenance/PUBLISHED_SOURCE_FILES.json) records original and published SHA256. Most files are byte-identical; large text metadata use documented lossless deterministic gzip. Decompression recovers exact original bytes. Existing gzip files remain unchanged. No scientific source/result is rewritten by this packaging.

Native v1 objects precede the revision's shared All/EX size mask. Apply v2 masks/annotations by stable nucleus ID. A Seurat `representation_HVG_only.rds` is an HVG representation, not a full-gene raw count store; use the corresponding shared integer-count AnnData separately. Large processed objects are release downloads, not Git payloads.

Raw source matrices, redundant checkpoints/graphs, virtual environments and controller caches are not copied here. A clone alone is not a complete raw-to-final reproduction environment.
