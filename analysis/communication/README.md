# Completed SATB2 CA1–inhibitory communication analysis

This folder preserves the executed source and runtime provenance of completed native LIANA scoring and figure production. The atlas uses finished saved results; publication does not rerun this workflow.

- Final completed registered run: `d397e9d9-03bc-43b0-9bd8-901409248fd6`.
- LIANA: 2.0.0, official immutable commit `c59472ccc9de8360dbbf5016db75f8abde08dd3e`.
- `native_source/`: final executed rendering/validation source and configuration. Its `provenance/` subfolder preserves the actual executed native scientific model source.
- `environment/`: frozen dependencies, original artifact hashes, API/runtime audits and environment notes.
- `references/mouseconsensus.tsv`: frozen native LR resource.
- Exact result tables and original matching DESeq2 exports are under repository-root `communication_data/`.
- Compact offline interactive payload and verification manifests are under `site/communication/data/`.

Some original executed configuration and runtime notes contain absolute paths to the source Mac. They are preserved as execution evidence, not presented as a portable one-command rerun. Input hashes identify original count matrices and frozen Seurat memberships; those large matrices are not duplicated here. Reproducing native scoring requires those original identified inputs and a compatible runtime.

For the static publication data contract, score semantics, checksums and regeneration instructions, see [CCC_DATA_NOTES.md](../../docs/CCC_DATA_NOTES.md). The standard-library publication packager is [package_communication_atlas.py](../../scripts/package_communication_atlas.py).
