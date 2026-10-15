# Deployment readiness review

Reviewed on 15 October 2026. **Local deployment checks pass.** Remote publication and HTTPS browser acceptance remain separate checks; this document does not claim the site has already deployed.

## Corrections verified

Two concrete blockers were identified and corrected before deployment:

- The validator now reads `site/data/DE_CATALOG.json` with the exact filename case. The former lowercase spelling could pass on a case-insensitive macOS filesystem and fail on the Ubuntu Actions runner.
- The `github-pages` environment now permits the `atlas-v*` tag pattern, in addition to the `main` branch. The workflow's version-tag trigger therefore matches an allowed deployment policy. The public environment-policy API returned the tag policy after the correction.

The publisher separately confirmed Pages is enabled with the GitHub Actions build type. This is required because the pinned `configure-pages` action defaults to leaving enablement unchanged.

## Validator and site paths

Running `scripts/verify_site.py` with the final catalog spelling produced:

| Check | Result |
|---|---:|
| Scientific payload shards, compressed and decoded hashes | 3140 verified |
| Current complete `site/` tree | 892,020,921 bytes |
| Deployment budget | Less than 950,000,000 bytes |
| Datasets | 24 |
| Full mouse-display gene axis | 32285 genes |
| Mouse display endpoint records | 1060 |
| Scientific execution performed by the validator | No |

An independent review checked 3166 HTML, catalog and shard paths against the actual spelling of every filesystem component. All resolved with exact case. In particular, the page loads `data/de_catalog.js`, which exists with that lowercase JavaScript filename; the JSON companion has the distinct uppercase name `DE_CATALOG.json`.

The Pages artifact contains only `site/`. Large exact-table payloads remain in the repository's `data_downloads/` directory and are fetched through the fixed `atlas-v2.1.0` raw-content URL. Both locations are checked by the validator. Source archives, local upload staging, raw input matrices and virtual environments are excluded from the Pages tree.

## Workflow review

The single-job workflow has the required `contents: read`, `pages: write` and `id-token: write` permissions, a `github-pages` environment, a deployment URL output, and an explicit `site` artifact path. Its `atlas-v*` trigger matches the permitted tag policy; manual dispatch remains available. A separate build-job dependency is unnecessary because validation, artifact upload and deployment are sequential steps in one job. This layout is supported by the [official custom Pages workflow documentation](https://docs.github.com/en/pages/getting-started-with-github-pages/using-custom-workflows-with-github-pages).

All four pinned action revisions resolved successfully to official action definitions:

| Action | Pinned revision | Implementation |
|---|---|---|
| Checkout | `3d3c42e5aac5ba805825da76410c181273ba90b1` | Node 24 |
| Configure Pages | `45bfe0192ca1faeb007ade9deae92b16b8254a0d` | Node 24 |
| Upload Pages artifact | `fc324d3547104276b827a68afc52ff2a11cc49c9` | Composite action |
| Deploy Pages | `368f82528645a54fb793d4d04e342629a3f51346` | Node 24 |

The final site remains below the [documented Pages size limit](https://docs.github.com/en/pages/getting-started-with-github-pages/github-pages-limits). Uploads use the official artifact action, which prepares the artifact format expected by Pages.

## Hosted runtime contract

The loader accepts only project-relative `data/shards/` filenames or the exact pinned prefix:

`https://raw.githubusercontent.com/Ali-Mujahid/satb2_snRNA/atlas-v2.1.0/data_downloads/`

It verifies the received gzip SHA256, decompresses explicitly, verifies decoded JSON SHA256, and parses data without evaluating downloaded code. HTTPS provides the secure context required by the browser's hashing API. The original offline download package retains its separate classic-script loader for local-file use.

The canonical data globals, dataset IDs, fit IDs, source count branches, exact membership guards, unavailable values and source-table text are unchanged. No scientific calculation, fitting, count extraction or annotation reassignment is performed by this publication layer.

## Final publication acceptance

At review time the complete data directories and immutable version tag still needed to be pushed. Therefore acceptance should finish by checking:

1. The `atlas-v2.1.0` tag identifies the complete data commit, including all referenced `data_downloads/` files.
2. The Pages workflow succeeds from that version or an explicitly selected permitted reference.
3. Actual HTTPS dataset and exact-table downloads return the expected bytes and hashes.
4. Release PDFs/documents and all complete archives have the recorded filenames, hashes and size limits.
5. The hosted browser controls, exports and cross-origin source-table downloads pass their separate browser review.

The expected project site is [the SATB2 atlas](https://ali-mujahid.github.io/satb2_snRNA/). Those final checks remain the publisher's responsibility and are not inferred from local success.

## Reviewed file identities

- `site/app.js`: `d462d5fb7447286f0c59640d31f56ec15fff773055d994292d7a9329f927d809`
- `scripts/verify_site.py`: `522fa1325793b018834d5f1f3b061fe34b3dd26ae487a7b8ec23fab91ee3448c`
- `.github/workflows/pages.yml`: `24277a2aa69347a98544084920cc5b40b7918a9cbcb4a3dba070636a1abb5ff7`

Later changes to these files need corresponding acceptance evidence. This review made no browser-control calls and did not modify the application, workflow, validator, original offline atlas or native analysis results.
