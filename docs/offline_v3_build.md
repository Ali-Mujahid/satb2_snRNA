# Offline atlas v3 packaging

`scripts/package_marker_atlas_v3.py` packages the completed classic-script offline atlas into two release assets. It does not build the UI, extract RNA, calculate marker summaries, fit DE models or recompute embeddings. The final marker-extension receipt records the separate descriptive aggregation of the added genes. Previously produced DE, counts, size factors, embeddings and saved-result files retain their bytes.

Both new archives use the same top-level folder:

- `SATB2_atlas_v3_runtime.zip` contains `offline_dashboard_v3/` except its `saved_results/` subtree.
- `SATB2_atlas_v3_saved_results.zip` contains `offline_dashboard_v3/saved_results/`.

Extract both archives into one parent directory, then open `offline_dashboard_v3/index.html`. They are independent v3 downloads; the original v2 runtime, corrected v2.1 runtime and v2 saved-results archive remain unchanged. Each release asset must be **strictly below 2 GiB**. A ZIP's file bytes and SHA are established only after actual packaging finishes.

Wait until the v3 folder contains the final classic UI, added marker payloads, fixed-order preset, browser examples and metadata. No writer should modify that folder between the snapshot and packaging commands. The output directory must be outside the bundle. The directory name must be exactly `offline_dashboard_v3`.

```bash
python3 scripts/package_marker_atlas_v3.py --self-test

python3 scripts/package_marker_atlas_v3.py --prepare-manifest \
  --bundle-root /path/to/offline_dashboard_v3 \
  --extension-receipt /path/to/marker_extension/COMPLETE.json

python3 scripts/package_marker_atlas_v3.py --package \
  --bundle-root /path/to/offline_dashboard_v3 \
  --extension-receipt /path/to/marker_extension/COMPLETE.json \
  --runtime-source-zip /path/to/SATB2_atlas_v2_1_runtime.zip \
  --saved-source-zip /path/to/SATB2_atlas_v2_saved_results.zip \
  --outdir /path/to/new_v3_archives
```

The explicit `--prepare-manifest` command replaces the copied folder's `BUNDLE_MANIFEST.json` with a faithful v3 file snapshot; its prior SHA is recorded. It enumerates every file, including hidden files, with relative path, bytes and SHA256, and excludes only the manifest itself. Symlinks and unsafe names are refused. The manifest binds the completed extension receipt and distinguishes added-marker expression aggregation from DE fitting and embedding recomputation. Packaging rescans all files against that snapshot before and after archive creation. A changed file, missing file, new file or changed receipt stops acceptance.

The helper clones the old archives on APFS when available, otherwise makes a byte copy. It then renames the common folder prefix from `offline_dashboard_v2/` to `offline_dashboard_v3/`, retaining existing compressed payloads. `zipnote` is attempted first. On this macOS runtime, `zipnote -w` failed a fabricated single-file archive with a temporary-file error. The tested fallback updates only equal-length local and central filename headers; offsets, CRCs and compressed payloads stay unchanged. It uses Python's ZIP64-aware offsets and refuses unsupported Unicode-path override extra fields. A failed `zipnote` attempt must leave the staged archive unchanged before fallback is permitted.

Only changed or added runtime files are recompressed by Info-ZIP. Replacement does not depend on mtime, so an older timestamp cannot cause a changed UI file to be skipped. Runtime entries present in the owned baseline but absent from the frozen v3 folder are removed. The saved-results companion is stricter: any added, removed or changed file causes refusal before its archive copy is made. It is a folder-prefix rename of the original saved-result bytes.

Every final ZIP file is read through EOF to check CRC and decompressed SHA256 against the folder snapshot. Every unchanged member additionally has compression method, compressed size and a streaming SHA256 of its compressed payload compared with the original archive. Thus unchanged data are checked for both logical identity and compressed-byte identity. Original source archive SHA256 values are rechecked at the end. Existing v3 output assets or receipts are never overwritten. Failed staging is retained under the output directory for review.

The package output directory receives:

- `OFFLINE_V3_PACKAGE_RECEIPT.json`: complete local audit, including original archive paths and every unchanged member's compressed-byte proof. This is a private uploader receipt.
- `OFFLINE_V3_PUBLIC_PROVENANCE.json`: public-safe summary with relative member names, archive hashes, sizes, changes, extension binding and acceptance results. It contains no private absolute source paths.
- The two verified v3 ZIP files.

The synthetic test uses temporary fabricated ZIPs only. It covers automatic and forced fallback prefix rename; unchanged compressed payload identity; changed files with older mtimes; added and deleted runtime files; all-member CRC/SHA; stale snapshots; and refusal of a changed saved-results file. It does not read the study data or the real v3 folder. Actual release readiness additionally requires the finished UI/data validation and the root reviewer's browser acceptance; passing this fixture alone does not establish that the release archives were built.
