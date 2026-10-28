# Official mouse cisTarget references, retrieved 28 October 2026

`READY_MANIFEST.json` records the exact official primary HTTPS URLs, response lengths, local SHA256 hashes, and download provenance. The historical `../blocked_resource_manifest.json` remains unchanged and describes the earlier connection failure only.

All four resources were downloaded with normal certificate and hostname validation into fresh files, then renamed atomically after their response lengths were verified. Both Feather databases also match freshly retrieved official SHA1 sidecars. The optional official SHA256 registry and TF-list/annotation checksum sidecars returned HTTP 404; no upstream checksum claim is made for those files. Their exact HTTPS response lengths and locally computed SHA256 hashes are recorded.

The pinned pySCENIC Python environment accepted both Feather V2 databases through `ctxcore.rnkdb.FeatherRankingDatabase`, and `pyscenic.utils.load_motif_annotations` accepted the MGI v10 annotation table. `SCHEMA_INSPECTION.json` contains their dimensions, ID coverage, parser results, and sample ranking-range checks. All 5,032 ranked motif IDs have native-filtered MGI annotations. The 500bp database contains 24,068 genes and the 10kb database contains 24,130 genes; each retains its native gene axis.

The current official mouse TF list contains 1,860 unique TFs, including `Satb2`. It is preserved exactly, without manual additions. This differs from the older GitHub protocol example audited during the outage.

This folder establishes reference availability and compatibility only. It contains no study-data inference, fitted GRN, regulon, enrichment result, or activity score. Workflow execution readiness is a separate check.
