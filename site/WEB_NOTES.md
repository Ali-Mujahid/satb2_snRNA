# Hosted atlas engineering

This is a separate web adaptation of the preserved offline v2 explorer. Scientific arrays, source numeric values, nulls, endpoint identities and fit membership rules are preserved. Controls change views and local display state; this page performs no embedding fit, raw-nucleus RNA extraction or differential-expression model.

Interactive gzip binary shards use same-origin project-relative paths. Exact source-download payloads use only the pinned `atlas-v2.1.0` path under `raw.githubusercontent.com/Ali-Mujahid/satb2_snRNA/data_downloads`. The loader verifies compressed bytes against their SHA256, decompresses explicit gzip magic with a browser-native stream or bundled fflate, verifies the decoded JSON payload hash and parses JSON. Remote JavaScript is not evaluated. HTTPS/WebCrypto is required for the hosted loader. The decoded cache is bounded to six entries and two simultaneous decodes.

The source-data manifest, DE catalog and mouse catalog remain classic local JavaScript globals. A new data version requires a new pinned tag and regenerated descriptors. `published_downloads.js` holds native/PDF/archive release links and integrity metadata; `ADDITIONAL_GALLERY_ASSET_INDEX.tsv` lists six matching example PDFs for the release publisher. The publisher must fill final archive SHA256/size fields after packaging. Native download and review-document links are rewritten by the data builder; none should point to the omitted multi-GB offline folder.

The original offline dashboard remains the version to open without a network. The hosted page downloads public static data; local labels/settings are not uploaded. UMAP-only genotype/programme filters do not alter the saved inference population. Novel exclusions/unions without an exact saved membership show NO MATCHING FIT. Missing or unsupported size factors remain NA.

Engineering receipts and the actual-browser checklist are maintained in the analysis project's `reviews/github_publication_20261015/` directory. A VM test is not a browser witness. Root owns deployed-browser testing, release publication, final archive integrity and deployment acceptance.
