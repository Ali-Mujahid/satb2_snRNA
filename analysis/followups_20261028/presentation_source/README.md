# Portable native-results reader builder

This source packages a **read-only** offline HTML reader for the four separately requested follow-up tasks. It does not infer communication, trajectories, composition, RNA velocity or new statistics. The embedded JSON is only a navigation index of native files/statuses and explicitly root-approved interpretation; it does not replace a native scientific schema.

No final study reader is built yet. `config.template.json` is disabled. Actual study building requires all three native terminal seals and completed independent/root reviews, then a separately saved approved configuration. The reader has four task sections: scCODA composition, CellChat, Monocle3, and an explicit **INPUT_BLOCKED** velocity assessment with no fabricated velocity/latent-time plots.

## Read-only packaging

Use `build_reader.py --configuration APPROVED_JSON --destination FRESH_DIRECTORY` after the root supplies final approved findings/limitations and exact terminal/review receipts. The destination must not exist. The builder only writes its fresh presentation directory; it never alters run trees, source native tables, models, GitHub or releases.

- Composition: `satb2-composition-final-validation/v1`, accepting both COMPLETE and COMPLETE_WITH_INCONCLUSIVE_OR_UNTESTABLE_SCOPES; preserves original per-model statuses and Bayesian versus frequentist distinctions in approved interpretation.
- CellChat: `satb2-cellchat-final-validation/v1`, COMPLETE. Actual scientific root is nested `results/cellchat`.
- Monocle3: `satb2-monocle3-final-validation/v1`, COMPLETE. Scope statuses, lack of biological temporal root and genotype-DE flags remain native and visible.
- Velocity: current `satb2-scvelo-starting-point-assessment/v1` must retain INPUT_BLOCKED_NOT_EXECUTED and scientific_jobs_started=false.

Every native receipt is verified, including archive-only objects. Symlink/path escapes, missing/changed/duplicate outputs, incomplete seals, synthetic fixture roots, missing root approval/reviews and unsafe archive URLs fail. The reader does not independently certify scientific correctness; it requires root-approved interpretation plus actual review files.

Native PNG/PDF bytes and small tables/supporting files are copied unchanged, with relative links. The HTML needs no server, libraries/CDNs or network. Keep the complete folder together. Scope/cohort/method/direction/plot selectors intersect available native files. PNG plots link to their original paired PDFs. Table previews display up to60 existing rows as raw strings; downloads retain complete unchanged native tables. No display filter changes any fit, exclusion, threshold or statistical call. No calculations or custom plotting happen in the browser.

All original native status strings remain in scope cards and the provenance panel. Missing/inconclusive outcomes are visible, not converted to negative findings. Velocity has assessment/review links and a clear blocked panel, never an empty fake plot gallery. Model/raw-axis files and large tables remain separately preserved native artifacts, listed with hashes/bytes in NATIVE_ARCHIVE_ONLY_FILES.json. Optional approved HTTPS full-archive links appear alongside the reader's own asset receipts; absent archive URLs are honestly marked as separately preserved, not represented as downloadable.

## Optional verified CellChat display repair

The CellChat task can supply `presentation_override` only after the separate plot-repair native run completes and an independent visual review is approved. It must contain `approved_for_display:true`, `presentation_root`, exact `terminal_seal` receipt and nonempty `visual_review_receipts`. The builder requires schema `satb2-cellchat-layout-only-final/v1`, a terminal COMPLETE or COMPLETE_PENDING_VISUAL_REVIEW status plus the explicit approval/review, all preservation flags false, and the repair's exact original scientific-seal hash. The mapping file must itself be sealed, every original/corrected PNG/PDF must match its own receipt, relative identities must match, and every original native figure must be covered exactly once. A missing review or altered mapping fails before creating a reader.

With a verified repair, “Readable layout” is the default CellChat figure rendering. An explicit rendering selector switches between readable and “Archival native (known display issues)” figures; corrected cards link both their readable PDF and the clearly labeled archival original PDF. Tables remain the unchanged original scientific sources under either figure rendering. Both original and corrected assets, terminal seals, mapping and visual reviews are preserved with relative portable links. An always-visible CellChat archival warning explains clipped labels/titles and that original mouse-boxplot points had vertical jitter, so their y positions do not equal raw saved scores. The absence of an override preserves the archival originals and warning; no unavailable repair is fabricated. This optional step never changes native results or statistics.

Default package budget45MB and per-table budget2MB keep the reader modest for a GitHub site. All generated PNG/PDF outputs are included; exceeding total budget fails rather than silently dropping scientific figures. Raw nucleus/feature-axis tables and per-chain model objects are archive-only. The root can explicitly adjust these presentation budgets, never scientific parameters. The final package includes READER_MANIFEST.json, unchanged terminal seals/reviews, approved presentation configuration and source/copy hashes.

## Verification

`tests/test_builder.py` uses tiny fake **test-only** sealed file fixtures, no study data, to check terminal/approval guards, receipts, path containment, unchanged native-byte transport, table preview, task preservation, velocity blocking, local PDF links, bounded size and no overwrite. No real reader is produced by these tests; their temporary fixture folders are removed.

Rendered browser checking must follow actual final build: all four task selectors, intersecting plot controls, empty-filter messages, original PDF/table links, table previews, blocked velocity panel, narrow/wide layouts and offline relative-link behavior. That UI review is separate from the native scientific review and from source unit tests. No GitHub deployment is attempted by this source.

## Complete large Monocle3 gene tables

Only oversized native `graph_test_all.tsv` tables receive deterministic gzip download transport (`mtime=0`, empty embedded filename). No rows, columns, values or formatting are changed. Compression occurs in a temporary presentation staging directory, outside all run trees. The builder verifies each decompressed SHA256 against the sealed original, records original and compressed receipts in the reader manifest, and supplies a relative `.tsv.gz` download link with the original byte count/hash. Native models, nucleus axes and all other oversized files remain separate archive artifacts. A compressed table contributes its actual compressed size to the explicit package budget.

## Small online transport beside the complete offline package

The existing atlas already uses approximately930MB, and the regulatory publication gate is950,000,000bytes. Do not copy a60MB portable reader into that site or raise the existing gate merely to fit it. `build_online_reader.py` instead stages a **separate small online transport of the same reviewed portable package**. Large PNG/PDF/gzip/table assets can reside as individually downloadable GitHub release assets; the online index retains the same native status, interpretation, controls, previews, hashes and original/readable figure distinction. The full portable ZIP keeps every local relative asset for offline use. No rendering, native-data transformation or scientific analysis occurs in this transport.

The online builder never uploads, mutates the existing site or deploys. It first verifies every portable file receipt. It requires an explicit approved online configuration, portable visual-review receipts and independent remote-upload verification receipts. Every external mapping must match its sealed local byte size/SHA256 and a separately recorded matching remote byte size/SHA256. Accepted verification methods are `github_server_digest` or `download_sha256`; a precomputed URL alone is not verification. Canonical download links are restricted to the authorized `Ali-Mujahid/satb2_snRNA` repository. An incomplete, changed or unverified mapping fails before the staging directory is created.

After the actual portable reader has passed review, the following **read-only proposal** lists uniquely named release assets for upload. It does not upload anything and marks every remote receipt unverified:

```text
python build_online_reader.py plan-assets --portable-root /absolute/reviewed/portable/folder --tag REVIEWED_TAG
```

The default proposal sends every PNG/PDF/gzip asset and every other displayed asset larger than20KB to the release; smaller source/status tables stay local. Asset names contain a stable path hash so files from different methods/scopes cannot collide. The root should upload the exact proposed bytes, verify server bytes/digests or downloaded SHA256 values, and save the actual tool/API receipts before completing `online.config.template.json` in a separately approved configuration. Each `verified_release_assets` row has `portable_path`, canonical `url`, `bytes`, `sha256`, `remote_verified_bytes`, `remote_verified_sha256` and `verification_method`. It must retain the native exact-byte identity. The separate `portable_zip` row binds a local ZIP receipt, its verified remote receipt and optional `archive_prefix`. The ZIP's complete file list and streamed decompressed hashes must match the reviewed portable package, including its portable manifest, with no missing or extra files.

```text
python build_online_reader.py stage --configuration /absolute/approved-online-config.json --destination /absolute/fresh-online-stage
```

The destination must be separate from both the complete portable package and the existing site. Both an explicit online-size cap (default10MB) and the existing site's current bytes plus the new stage (default950,000,000bytes) are checked, including final receipt overhead. Native assets mapped remotely are omitted only from the small online transport, never from the complete offline ZIP; other small files are copied unchanged. `ONLINE_TRANSPORT_MANIFEST.json` records original and remote receipts, package identity and size evidence. It is a delivery/provenance record, not a replacement scientific schema. The page clearly explains that the online preview needs the network for some files and provides the verified complete offline ZIP link.

Publication still requires actual release-asset checks and live browser verification: decoded images, all four task sections/filters, corrected versus archival links, full native table downloads, PDF downloads, blocked velocity text and the offline ZIP. An asset-delivery error is visibly distinguished from scientific unavailability. Source fixture tests do not establish that GitHub-hosted image delivery works; test a real uploaded PNG in the intended browser before deploying this transport.

The current source regression suite has21groups (13portable-builder groups, including lossless gzip, and8online-transport groups). Online tests use only tiny explicitly fake fixtures. They cover exact local/remote transport identity, independent-review and approval requirements, unsafe URLs/paths, incomplete ZIPs, symbolic-link ZIP metadata, site/online size limits, no existing-site writes, no overwrite and preserved inconclusive/blocked status. No actual study reader or online site was built by these tests.
