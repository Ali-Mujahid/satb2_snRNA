# Reproducibility and portable paths

V1 creates fresh Scanpy/Seurat representations from previously audited integer counts and reuses source-verified Allen MapMyCells labels. V2 consumes accepted v1 outputs and adds exhaustive reviewed annotations, shared filtering, category/program memberships, revised DESeq2 inference, full-cohort marker figures and technical controls. V2 does not repeat upstream preprocessing or independently remap the atlas.

| Goal | Required material |
|---|---|
| Explore existing UMAPs, marker dots and cached DE | Complete viewer/runtime and saved-result delivery; no Python/R execution |
| Inspect exact gene-by-mouse values/DE | Produced full-gene arrays, size-factor/model settings, seven-column tables, stable-ID and eligibility audits |
| Repeat aggregate inference | Exact arrays, mouse/condition map, fixed membership/gene universe, model policy, matching code/runtime; this is new scientific execution |
| Load existing native representations | 28 processed-object assets, their manifests and v2 annotations/masks |
| Repeat original representations | V1 code/runtime plus all source-bound matrices, rosters/QC, query/reference evidence and feature manifests |
| Repeat the complete revision | Restored accepted v1 tree, v2 code/runtime, source-bound ledgers and a fresh output root |

This repository alone does not promise a one-command raw-to-final rerun. The original workflow begins with audited matrix/count artifacts, not FASTQ. Cell Ranger and existing reviewed CellBender corrections are upstream inputs; source queries/matrices are not silently regenerated.

## New configuration on another machine

The [v1 original](../../analysis/v1/config/analysis_v1_0_5.original.json) and [v2 original](../../analysis/v2/config/analysis_v2_0_0.original.json) are unchanged. Their absolute paths preserve actual provenance. Portable example copies under [analysis/portable](../../analysis/portable) have explicit placeholders; the clean-v1 example removes historical recovery-cache keys.

1. Restore needed data into a separate working tree, checking branch, stable IDs and source SHA rather than filename alone. Release/result manifests identify exact produced objects.
2. Decompress packaged reference metadata where expected. The publication ledger records original locators/SHA and any lossless gzip transform.
3. Copy `analysis/portable/path_map.template.json` outside the immutable snapshot and fill every absolute root/executable placeholder. Prepare compatible environments separately from actual recorded versions.
4. Create a new JSON config with the standard-library helper:

   ```bash
   python analysis/tools/materialize_config.py --config analysis/v2/config/analysis_v2_0_0.original.json --path-map /absolute/filled_path_map.json --output /absolute/new_run/config.json
   ```

   For clean v1 use its original config and add `--drop-v1-reuse`. Cache reuse instead requires exact completed bound files/receipts; similar names are insufficient.
5. Remap paths inside referenced manifests separately and deliberately. The helper edits JSON only; it does not edit TSVs, receipts, AnnData/R objects or lazy payloads. Preserve original hashes, and record transformed publication/working copies separately.
6. Check required inputs and actual runtime before execution. The helper writes a config/remapping receipt; it neither establishes scientific readiness nor launches Snakemake/R/Python fits.

The revision binds 80 frozen inputs, including upstream counts, metadata, graphs, PCs and UMAPs. Some are separate object downloads; original count/query data and certain native dependencies still require the original study inputs. Not every original path exists after cloning Git.

Use actual `PACKAGE_VERSIONS.tsv` / `SESSION_INFO.txt`, not planning estimates. Recorded R 4.6.1/DESeq2 1.52.0 and Seurat 5.5.1 are preserved, as are distinct Python/R Harmony versions and Python freeze/lock/wheel records. A newer environment does not imply identical stochastic graphs, clusters, UMAPs or fits.

V1 retains unintegrated/Harmony and SCT controls; v2 reported inference uses DESeq2. Identical shared arrays/settings reproducing the same R engine is not an independent statistical implementation test. IN includes all available nuclei with positive-library nonempty mouse aggregates and at least 3 independent mice/genotype; sparse coverage remains exploratory. Other groups retain the documented 20-nuclei minimum. Controlled raw/corrected comparisons fix membership and gene universes.

Display suggestions remain nonbinding. Programmes describe expression with physical layers undetermined. An alias or custom selection does not authorize reuse of an unrelated cached fit: exact produced membership and eligibility are required.
