Scanpy source preparation for the newly authorized matched reanalysis. These files
have been authored and syntax checked; scientific stages must run through the
parent's reviewed workflow execution. Existing projects and the mapper runtime
are read-only inputs.

IO contract

- `sample_sheet`: TSV with `sample`, `condition`, `source_h5`,
  `source_h5_sha256`; the accepted ingest table also supplies `original_UMIs`.
- `qc_roster`: accepted branch-specific TSV with `nucleus_id`, `barcode`,
  `sample`, `condition`, `passes_primary_QC`, and population flags `REF_EX_HS`,
  `REF_IN_HS`, `REF_CA1_PROS_HS`. IDs remain `sample + '_' + unchanged barcode`.
- `mapmycells_predictions`: original-count prediction TSV with exactly the same
  nucleus IDs. Imported columns receive `mmc_` prefix. Reuse is recorded as common
  annotation evidence, not as a newly executed independent mapping.
- Optional `markers_tsv`: columns `group`, `gene_symbol`, optional `eligible`;
  descriptive expression evidence only, no forced anatomical subtype calls.

Commands, to be launched by the reviewed Snakemake workflow

```
<new-venv>/bin/python <source>/prepare_common_inputs.py --config <config.json> --outdir <run>/results/shared
<new-venv>/bin/python <source>/run_scanpy.py --config <config.json> --shared-dir <run>/results/shared --population allcells --outdir <run>/results/scanpy/allcells
```

Run the second command separately and sequentially for `Excitatory`, `Inhibitory`,
and `CA1_ProS`. Frozen reference flags select those matched populations. CA1-ProS
remains a provisional reference population that includes ProS; no prior cluster
labels become new results.

Preparation outputs under the passed `--outdir`

- `branch_counts.h5ad`: CSR declared-branch integer UMIs in X, all genes, selected
  nuclei, authoritative metadata, verified common MapMyCells annotations.
- `nucleus_metadata.tsv`, `gene_features.tsv`.
- `membership_<population>.tsv` and `common_hvg_<population>.tsv` for each population.
- `COMPLETE.json`: exact input and output hashes, count digest, versions, coverage.

Representation outputs under each population's passed `--outdir`

- `scanpy.h5ad`: sparse log-normalized X; untouched input-branch UMIs in
  `layers['counts']`; PCA, Harmony2 PCs, unintegrated/Harmony2 graphs and UMAPs,
  fresh Leiden grid and main clusters. Scaling is restricted to common HVGs.
- `nucleus_metadata.tsv`, `cluster_counts_per_mouse.tsv`.
- `cluster_markers_unintegrated.tsv`, `cluster_markers_harmony.tsv`: descriptive
  Wilcoxon cluster-vs-rest annotation markers, with equal per-mouse sampling within
  each cluster capped at100; not mouse-level genotype tests.
- `marker_selected_roster_unintegrated.tsv`, `marker_selected_roster_harmony.tsv`.
- `representation_review.pdf`, optional marker-panel expression/audit TSVs.
- `COMPLETE.json`: count-preservation assertion, versions and output hash.

Configuration parameters

```
{
  "sample_sheet": "/absolute/accepted/PER_SAMPLE_ORIGINAL_COUNT_ARTIFACTS.tsv",
  "qc_roster": "/absolute/accepted/SELECTED_NUCLEUS_ACTIVE_METADATA.tsv",
  "mapmycells_predictions": "/absolute/original/reference_predictions_all_nuclei.tsv",
  "expected_samples": 12,
  "expected_nuclei": 108062,
  "count_branch": "original_full",
  "keep_column": "passes_primary_QC",
  "populations": ["allcells", "Excitatory", "Inhibitory", "CA1_ProS"],
  "n_hvg": 3000,
  "hvg_flavor": "seurat_v3_paper",
  "feature_discovery_exclude_patterns": [],
  "target_sum": 10000,
  "scale_max": 10,
  "n_pcs": 30,
  "n_neighbors": 20,
  "umap_n_neighbors": 30,
  "umap_metric": "cosine",
  "umap_min_dist": 0.3,
  "resolutions": [0.2, 0.4, 0.6],
  "main_resolution": 0.4,
  "harmony_key": "sample",
  "harmony_theta": 2,
  "harmony_sigma": 0.1,
  "harmony_lambda": 1.0,
  "harmony_nclust": 100,
  "harmony_max_iter": 10,
  "harmony_max_iter_kmeans": 4,
  "harmony_epsilon": 0.01,
  "seed": 20261013,
  "threads": 4,
  "minimum_available_gib": 24,
  "minimum_marker_cluster_nuclei": 20,
  "marker_nuclei_per_mouse_cluster": 100,
  "marker_top_n": 100
}
```

`seurat_v3_paper` selects batch-aware integration features from original counts,
then the same gene IDs are fixed in both tool branches. Optional discovery
exclusions must be explicitly set by the agreed scientific design; all genes
and all original counts remain saved. Harmony2 parameters should match R's
reviewed configuration; this implementation rejects condition as its batch key.
With fixed lambda1 and 100 Harmony clusters, the explicit defaults above replace
the package's automatic lambda choice. The unintegrated control is retained.

Paired denoising sensitivity branches set `count_branch` to `paired_original` or
`paired_cellbender`, `expected_nuclei` to104298, and pass their explicitly audited
paired roster and source H5/hash sheet. Set `hvg_source_dir` to the original-full
shared-output directory; matching per-population panels are copied byte-for-byte.
`mapmycells_allow_superset=true` permits the original-full predictions to be
subset to the paired original cohort. A corrected branch must pass its declared
annotation source and `annotation_provenance`, stating fresh corrected-query
mapping or explicitly shared original-count annotation. Clustering populations
remain fixed by the common frozen roster for the paired comparison. Actual branch
library totals are recomputed; `frozen_original_nCount_RNA` preserves old metadata.
The original-full branch verifies per-mouse original UMI totals; paired-original
verifies per-nucleus original counts rather than comparing to full-cohort totals.

Runtime plan

Use bundled CPython3.12 arm64 to create a new project-only `.venv_scanpy`.
`requirements.in` pins the method-defining stable versions. Resolve ancillary
stable dependencies using `pip install --dry-run --ignore-installed
--only-binary=:all: --report <resolution.json> -r requirements.in`. Inspect
the report for CPython3.12/macOS arm64 or universal wheels and record every
resolved version/hash before actual installation. Then install into this venv
only and save `pip freeze`, `pip check`, and module import/version evidence.
Never modify bundled/system Python or the existing mapper environment.

Reserve approximately1.5GB for the environment plus0.5GB wheel/cache space,
at least10GB for sparse outputs (provisional estimates), and max20GiB RAM for
one stage. Require24GiB available RAM before each process; execute stages
sequentially with max4threads. Fine-grained peak RAM remains to be measured from
actual execution. Scientific fitting/importing the data has not run here.

Primary API sources reviewed

- https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.highly_variable_genes.html
- https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.leiden.html
- https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.rank_genes_groups.html
- https://pypi.org/project/scanpy/1.12.4/
- https://pypi.org/project/anndata/0.13.4/
- https://pypi.org/project/harmonypy/2.1.0/
- https://pypi.org/project/igraph/1.0.0/
- https://pypi.org/project/leidenalg/0.12.0/
- https://pypi.org/project/scikit-misc/0.5.3/

Independent Python pseudobulk count check

```
<new-venv>/bin/python <source>/export_pseudobulk.py --shared-dir <run>/results/shared/BRANCH --scanpy-root <run>/results/scanpy/BRANCH --branch BRANCH --outdir <run>/results/python_pseudobulk/BRANCH
```

This standalone export performs no normalization, statistical test or model fit.
It verifies the shared count artifact identity, reads CSR counts in4096-nucleus
backed chunks, and accumulates exact int64 gene-by-mouse sums. All12 mouse columns
are retained, including explicitly absent/ineligible population columns.
Common atlas families are `ATLAS_TOTAL_ALL`, `ATLAS_TOTAL_EX`, `ATLAS_TOTAL_IN`,
`ATLAS_TOTAL_CA1_PROS`, with labels `ALL`, `EX`, `IN`, `CA1_PROS`. Fresh Scanpy
Harmony cluster families are `SCANPY_ALLCELLS_harmony`,
`SCANPY_EXCITATORY_harmony`, `SCANPY_INHIBITORY_harmony`,
`SCANPY_CA1_PROS_harmony`; their labels are actual main-resolution cluster IDs.
These agree with the shared R pseudobulk population semantics.

Outputs: `population_index.tsv`, `population_memberships.tsv.gz`,
`all12_population_coverage.tsv`, `population_eligibility.tsv`,
`confirmed_all12_mouse_roster.tsv`, `source_count_conservation.tsv`,
`COMPLETE.json`; each population directory contains
`SOURCE_INTEGER_COUNTS_ALL_GENES_ALL12.tsv.gz` and `ALL12_POPULATION_COVERAGE.tsv`.
Count TSV columns are `gene_id`, `gene_symbol`, WT6 then cKO6 in the fixed R order.
Python and R aggregation comparisons should join by family, label and sorted
ID-membership hash, then assert every gene-by-mouse count exactly equal; ordinal
population keys can differ when R includes additional unintegrated populations.
Eligibility records the same20nuclei/mouse and3mice/genotype count requirements;
the trusted R stage owns model fitting, partition designs and gene prefilters.
Sorted ID hashes use UTF-8 newline-delimited IDs plus one final newline.

MapMyCells native annotation route (separate mapper environment)

Existing package: `cell_type_mapper1.7.4`; native entrypoint:
`python -m cell_type_mapper.cli.from_specified_markers --input_json <config.json>`.
Verified current local source schemas accept `query_path`, `precomputed_stats.path`,
`query_markers.serialized_lookup`, `type_assignment.normalization='raw'`,
`type_assignment.algorithm='hierarchical'`, bootstrap parameters, and output
paths. Pinned WMB reference20240831/CCN20230722 is already cached; no new data
upload/download is required for original-query annotation reuse.
