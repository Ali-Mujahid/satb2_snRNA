# Reviewed follow-up archive transport

This source packages existing files for the four separately requested SATB2 follow-up tasks. It performs **no fitting, scoring, normalization, eligibility selection, permutation test or new scientific calculation**. No actual study archive has been created by source preparation or its unit tests.

The disabled `config.template.json` names six components. Copy it outside this source, supply the remaining actual strict-sampler/reader terminal seals and review receipts, review the component definitions, then explicitly enable the copied configuration. The completed CellChat and Monocle3 components and the shared source/input component can be packaged separately while the strict sampler is active. Every invocation uses a different fresh destination directory on the server share. It cannot append to or overwrite an existing delivery.

```bash
python3 package_archives.py --configuration /absolute/approved-archives.json --component cellchat --component monocle3 --validate-only
python3 package_archives.py --configuration /absolute/approved-archives.json --component cellchat --component monocle3
```

Selections apply terminal/review gates only to the selected components. An unselected active strict run is neither read nor packaged. Omitting `--component` requires all six components. The program validates all selected gates, reviews, canonical paths, source hashes and conservative size budgets **before creating the destination**. The global approval field acknowledges the root agent's actual reviews; this transport code does not independently certify biological interpretation or visual correctness.

## Components and exact-byte preservation

1. **Initial scCODA v1.0.2:** every sealed count/proportion table, native chain/combined H5AD, authentic posterior/transition array, diagnostic, plot and status is retained. It completed computationally, but all sixteen models were inferentially inconclusive; native selected coefficients are not validated effects. The earlier failed attempts remain preserved locally, with their review/history records in the shared provenance selection.
2. **Strict same-model sampler v1.1.0:** separate archive after its actual native terminal seal and completed review exist. It preserves the same inputs/references/model and the single predeclared stricter integration check. A legitimate inconclusive outcome is accepted and remains explicit; this code does not change its gates or infer effects.
3. **CellChat and readable layouts:** original 115 sealed scientific files, all native objects and scores, plus the separate display-only correction and its native-to-corrected mapping. Original figures remain archival outputs with the known clipping/vertical-jitter caveats. The corrected layout must have an actual display-acceptance receipt bound to its exact seal. Native cell-label specificity probabilities are distinct from mouse-level genotype score tests.
4. **Monocle3:** all sealed four-model native CDS/RDS files, required nearest-neighbor indexes, graph/root tables, graph-associated genes, modules, raw/display aggregates, plots and support/QC tables. The output describes unrooted expression-state graphs; arbitrary-reference distance is not established developmental time. Original DESeq2 fits remain gene-genotype evidence.
5. **Shared frozen inputs and source:** exact frozen composition metadata/annotations, the prepared CA1+IN integer Matrix Market counts and cell/feature files, the five original matching-endpoint DESeq2 tables, executed/custom source/config snapshots, reference/runtime/install provenance and known native plan/review receipts. Native CellChat and Monocle3 source distributions are retained; external dependency URLs/digests/locks are provided rather than bundling mutable installed environments. Byte-identical shared files are stored once with explicit aliases.
6. **Portable offline reader:** the reviewed complete folder, with plots/PDF/table downloads and all relative links intact. This is separate from the small online index that references release-hosted assets. Do not substitute a partially built reader for a reviewed portable package.

Every component retains original native file bytes. The ZIP metadata uses a fixed date, safe portable permissions and sorted member order, so repeated packaging of the same bound inputs is deterministic. Executable script wrappers retain a portable executable bit; archived code is never run. The ZIP's `ARCHIVE_MANIFEST.json` is a **delivery receipt**, not a new scientific result schema. It records each requested relative member, original canonical source, size, SHA256 and stored member. Native scientific schemas, labels, tables, statuses and model objects are preserved. Original workflow snapshots/configurations and immutable run-relative paths should remain together.

## Integrity and safety

- The source accepts only explicitly allowed canonical roots and rejects symlinks, escapes, traversal, duplicate member collisions and unsafe/nonportable member names.
- Mutable software environments/caches are excluded. Workflow Python bytecode, `.git` and `.snakemake` state are not scientific source and are skipped from source trees.
- Each compressed archive must remain below **1,900,000,000 bytes**. A conservative uncompressed preflight leaves margin and avoids relying on an optimistic compression estimate. Every source is checked while streaming; every stored ZIP member is then read back and verified by SHA256. Large files stream in 1 MiB blocks.
- Staging goes to a fresh folder on the server share, not into a native run or the Mac's limited working volume. A failed fresh partial delivery is retained for inspection; no scientific file or prior delivery is deleted.
- `ARCHIVE_DELIVERY_RECEIPT.json` records actual ZIP sizes/digests and `BUILT_AND_READBACK_VERIFIED_NOT_UPLOADED`. It makes no GitHub or publication claim. Upload/server-digest/live checks remain with the parent agent.

Final release validation must require every configured component, even if they were built in separate invocations:

```bash
python3 package_archives.py --configuration /absolute/final-approved-archives.json \
  --verify-deliveries /server/delivery-A/ARCHIVE_DELIVERY_RECEIPT.json \
                       /server/delivery-B/ARCHIVE_DELIVERY_RECEIPT.json
```

This verifies exact expected component configuration, all current native/source receipts, each full ZIP and each decompressed member, rejecting duplicates or missing components. It is local transport verification, not publication approval. The final approved configuration must retain the already-packaged component definitions unchanged; only destination selection is operational and can differ between invocations.

## Opening and restoring the files

Open the portable reader's `index.html` after keeping its full folder together; no Python or R is needed for offline viewing. Native `.tsv/.txt/.csv` tables are text, `.mtx.gz` contains sparse integer counts, `.npz` contains authentic arrays, `.h5ad` uses AnnData, and `.rds`/native Monocle directories require compatible R/package versions. Gene symbols and original mapping tables are preserved. Biological mice, not nuclei, are the replicate unit. CA1 layer names are molecular programs and C8 uncertainty remains visible. There are no fabricated scVelo vectors or latent times: velocity remains input-blocked pending genuine splicing counts and LiSC access.

For complete path restoration, including deduplicated shared-input aliases, use the standalone standard-library extraction helper with the full ZIP SHA256 from the delivery receipt:

```bash
python3 extract_reviewed_archive.py /absolute/component.zip /absolute/new-extraction-folder \
  --expected-sha256 EXPECTED_FULL_ZIP_SHA256
```

The helper verifies safe members, exact sizes and SHA256 values, writes only into a fresh destination and executes no archived code. Standard ZIP tools also work for native result components; the shared-input archive needs alias restoration to recreate every requested file path. Source configurations retain the original paths for provenance. For a future reanalysis on another system, edit copies of those configurations to the new extracted paths and follow the recorded environment/reference requirements; do not edit the archived originals or describe the stored Mac runtime as a universal installer. Reproducing inference is a separately approved scientific execution.

## Source verification

Run only the tiny invented transport fixtures:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s tests -v
```

They check unchanged-byte round trips, deterministic archives, fresh-only output, approval/partial-state/hash/review guards, path/symlink and size limits, shared aliases, selected components, final all-component completeness and safe extraction. They do not open or package study matrices/models.
