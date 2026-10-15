# Seurat v5 implementation for the all-twelve comparison

Authored 13 October 2026. Only R syntax and installed APIs were inspected during
authoring; no data normalization, SCT, Harmony, clustering or model fit was run.
The parent controller owns execution, scientific review and completion.

Run one stage per R process, sequentially, through the approved controller:

```text
Rscript --vanilla source/seurat/run_seurat.R CONFIG.json log_global
Rscript --vanilla source/seurat/run_seurat.R CONFIG.json annotate_subsets
Rscript --vanilla source/seurat/run_seurat.R CONFIG.json markers
Rscript --vanilla source/seurat/run_seurat.R CONFIG.json sct_global
Rscript --vanilla source/seurat/run_seurat.R CONFIG.json pseudobulk
```

Each stage writes `output_root/STAGE/COMPLETE.json` only after successful execution.
`FAILED.json` retains its error. A controller-created empty stage directory is
allowed. Any existing content, including hidden files or partial/failure evidence,
is refused; use a new root for an authorized retry and preserve failed evidence. Sources, old
checkpoints and prior results are read-only. There is no installer or upload.

## Required contract

Config fields: `sample_manifest_tsv`, `roster_tsv`, `gene_manifest_tsv`,
`fixed_features_tsv`, `labels_tsv`, `output_root`, `resource_budget_gib`,
`memory_probe_python`, `minimum_available_gib`, `expected_nuclei`, `count_source`
and `seed`.

- Sample manifest: `sample`, `condition`, `sequencing_partition`, `filtered_h5`,
  `source_h5_sha256`. Accepted original-count RDSs can instead use `counts_rds`
  and `counts_rds_sha256`. H5 inputs take precedence. Each source is hash-checked
  before reading; source paths/hashes are retained in representation audits.
  H5 import reads only the declared `/matrix` CSC arrays, stable feature IDs and
  barcodes; extra CellBender latent root groups are ignored as nonmatrix metadata.
  Shape, pointers, indices, unique axes and nonnegative integer counts are checked.
- Frozen roster: `nucleus_id`, `barcode`, `sample`, `condition`, `nCount_RNA`,
  `passes_primary_QC`, `call_new`. IDs are exactly `sample_barcode`; all twelve
  samples must remain. Original full-library UMIs must equal `nCount_RNA` for
  original branches. CellBender libraries are freshly summed and must be positive;
  original `nCount_RNA` is preserved separately, never asserted equal to corrected
  totals. Noninteger counts are rejected, never rounded.
- Gene manifest: unique `gene_id`, `gene_symbol` for every source gene. Duplicate
  symbols remain separate stable gene IDs.
- Common features: TSV `gene_id`, exactly 3,000 count-based batch-aware HVGs from
  the original-full Scanpy preparation. All log/count branches use this panel.
- Labels: original source-verified MapMyCells predictions keyed `nucleus_id`,
  `reference_*` fields and `REF_EX_HS`, `REF_IN_HS`, `REF_CA1_PROS_HS`. Predictions
  are explicitly reused, not described as a new mapping. The original-full label
  roster must match exactly; paired branches may explicitly select from that
  original label superset. `CA1_supported` is optional and must come from a
  separately reviewed identity contract. Otherwise the CA1 subanalysis is named
  `CA1_PROS_CANDIDATE`; no anatomical layer assignment is imposed.
- Optional curated marker TSV: `set`, `gene_symbol`. It affects marker figures,
  never QC membership, clusters, labels or genotype tests.

`count_source` is `original_full` (108,062), `paired_original` (104,298), or
`paired_cellbender` (104,298). The controller supplies separate per-branch roots.
For paired count inference, set `cluster_metadata_root` to original-full Seurat
results and `prefilter_universe_root` to original-full `pseudobulk`; fixed original
memberships and frozen gene universes then isolate count-source effects. Original
full is primary; CellBender is a qualified count sensitivity.

## Representation and resources

Each capture is loaded individually. Full source libraries normalize the selected
HVG counts: `log1p(count * 10000/full_library_UMIs)`. Only HVG counts/data are
merged into the durable Seurat representation. This avoids incorrectly using
the HVG library total as the normalization denominator. Dense scaling is limited
to the selected panel; all genes are reopened from source for markers/pseudobulk.

Default log settings: centered/scaled RNA, clipping max 10, 30 PCs, exact RANN
neighbors `k=20`, Seurat SNN pruning `1/15`, igraph Leiden modularity with 10
iterations, resolutions 0.2/0.4/0.6, primary 0.4. UMAP uses uwot, cosine metric,
30 neighbors and min-dist 0.3. Harmony 2 explicitly uses sample (never genotype),
theta 2, sigma 0.1, lambda 1, nclust 100, max-iter 10 and four cores. Seurat SNN
graphs and Scanpy connectivity graphs remain native algorithm differences;
matching settings does not imply bit-identical graphs or UMAPs.

The graph adapter constructs a weighted undirected igraph object directly from
the sparse SNN lower triangle, excluding self edges and retaining every nucleus.
It then calls native Seurat `RunLeiden` and `GroupSingletons`, preserving the
declared modularity objective, random seed, iterations and singleton grouping.
This avoids the installed igraph 2.3.3 `mode='lower'` NA check, which allocates a
dense `lower.tri` index before sparse graph conversion. A small synthetic fixture
matches native FindClusters memberships and factor levels at all three
resolutions, including isolates; a 104,298-vertex sparse synthetic conversion
also succeeds without a nucleus-by-nucleus dense matrix. This is an API repair,
not a change to the biological cohort or clustering settings.

A 3,000 × 108,062 double scale matrix alone occupies approximately 2.42 GiB.
`PRE_FIT_RESOURCE_PLAN.tsv` is a conservative dense/capture planning bound before
fitting; `RESOURCE_ESTIMATE.tsv` refines it using actual HVG sparse sizes before
PCA/graphs. These are estimates, not measured peak RSS. A psutil stage-start guard
requires the configured available-memory minimum (normally 24 GiB). The parent
should also enforce sequential stages and record actual memory. Dense scale data
are removed immediately after PCA/Harmony, before graphs. `PRE_GRAPH_HVG_ONLY.rds`,
PCA/Harmony TSVs and count-conservation audits retain fitted evidence before
graph/clustering. Sparse SNN checkpoints and graph-size audits are exported per
branch. These artifacts preserve failure diagnostics; they are not COMPLETE
stage receipts or automatically accepted scientific results.

SCT runs only for original-full. It fits each capture on all genes with SCT v2,
no nuisance/genotype regression and min-cells 5. Residuals are requested for
eligible common features; the final SCT PCA panel is the intersection modeled in
all twelve captures. Missing residuals are never fabricated as zeros. Requested,
modeled, dropped and retained gene IDs are exported. SCT consequently represents
a model/normalization/feature-eligibility sensitivity, not a pure software change.
The primary matched log panel is unchanged.

## Outputs and biological inference

Global stage exports `nucleus_metadata.tsv`, `pca.tsv`, `harmony.tsv`,
`umap_unintegrated.tsv`, `umap_harmony.tsv`, `representation_HVG_only.rds`,
per-mouse/count-conservation audits and readable rasterized embedding PDFs.
Cluster fields are `log_unintegrated_r0p4` and `log_harmony_r0p4` at the primary
resolution. Subsets retain their own `EX_*`, `IN_*` or `CA1_PROS_CANDIDATE_*`
fields and local `log_*` aliases; the merged subset sidecar retains unique fields.

Marker discovery uses a recorded seeded sample-balanced subsample, normally at
most 600 nuclei per cluster and 30,000 overall. It cannot change the full cluster
or inferential roster. `MARKER_SUBSAMPLING_MANIFEST.tsv`, complete exploratory
marker TSV/CSV, top-seven marker tables, curated-marker PDFs and all dotplot
numerical values are exported. Nucleus-level marker P values are exploratory
rank information and are not animal-level WT/cKO evidence.

Pseudobulk streams each full-gene capture once into gene × mouse sums. Optional
`scanpy_metadata_paths` maps allcells/Excitatory/Inhibitory/CA1_ProS to Scanpy
metadata TSVs with `nucleus_id`, `cluster_harmony`, `cluster_unintegrated`.
Both tools' fresh cluster memberships are then fitted by the same trusted R
models, distinct from common atlas populations. Families use `ATLAS_*`,
`SEURAT_*` and `SCANPY_*` namespaces. `POPULATION_INDEX.tsv` binds each ordinal
directory key to its family/label; `POPULATION_MEMBERSHIPS.tsv` binds exact nuclei.
The index records `source_label_branch` to distinguish fixed original identities
from fresh count-branch cluster assignments.
`SOURCE_POPULATION_LEVEL_REGISTRY.tsv` preserves original full-cohort cluster and
external label names before paired selection. A source-defined endpoint emptied
by that intersection remains present with twelve zero counts, zero nuclei and
an ordinary ineligible status; names/counts are never fabricated. The external
high-support membership policy still applies to the retained raw label registry.

Each population retains all twelve count columns and explicit coverage. Eligibility
requires at least 20 nuclei/animal, three independent mice/genotype, full rank and
two residual degrees of freedom. Primary DESeq2 uses unpaired `~condition`, WT
baseline, cKO/WT effect, poscounts, Wald, no count replacement, default Cook
filtering, independent filtering and endpoint-specific BH 0.05. Gene prefilter:
10 UMIs in at least min(nWT,ncKO) included mice. edgeR robust QL uses the same
input gene universe. All unavailable/filtered/outlier values and ineligible
statuses remain visible. No fold-change significance gate is used.

Original-full partition sensitivities are conditional on recorded partition
provenance; cKO_933 stays unverified. `loo_families` limits leave-one-mouse-out to
common atlas total families, normally ALL/EX/IN/CA1_PROS. Paired branches should
set `run_loo=false`, `run_partition_sensitivities=false`; they run primary DESeq2
and edgeR on fixed original identities/universes. Primary MA/volcano PDFs and
LOO effect-range/direction audits accompany complete stable-ID/statistical tables.
Volcano plotting floors numerical zero adjusted P values for display only.

`pseudobulk_core.R` is copied from the preserved reviewed twelve-mouse module.
Its statistical and eligibility functions are reused; the new wrapper replaces
whole-dataset count assembly with streamed sums and supplies truthful count-kind
metadata for corrected-count sensitivity. Its original whole-object runner is
not called. Agreement on identical counts/memberships checks reproducibility;
it does not supply independent biological validation.
