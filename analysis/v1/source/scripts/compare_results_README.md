# Fresh result comparison source

`compare_results.py` compares completed fresh outputs. It does not fit models,
change assignments, integrate data, alter inputs or reuse earlier result folders.
It requires Python 3.10+, NumPy and Matplotlib. All computations use metadata
and exported stable-ID model tables; expression matrices are not loaded.

Run through the parent workflow after every required branch/population stage:

```bash
<python> <source>/scripts/compare_results.py \
  --run-root <new-run-root> --outdir <new-run-root>/results/comparison
```

The output directory must be empty. `COMPLETE.json` is written last, after input
hashes are rechecked and output hashes are recorded. Missing core datasets cause
failure; optional pseudobulk absence is reported without pretending it was run.
`--allow-partial` is explicit exploratory mode and gives a distinct partial status.
The final normal status is `EXECUTED_REVIEW_REQUIRED`; it documents comparison
execution, with no automatic
consistency verdict. Inspect the exported metrics, source coverage and PDF plots.

Required core inputs

- Scanpy: `results/scanpy/BRANCH/POPULATION/nucleus_metadata.tsv`.
- Seurat all cells: `results/seurat/BRANCH/log_global/nucleus_metadata.tsv`.
- Seurat subsets: `results/seurat/BRANCH/annotate_subsets/SUBSET/nucleus_metadata.tsv`.
- Branches: `original_full`, `paired_original`, `paired_cellbender`.
- Populations: `allcells`, `Excitatory`, `Inhibitory`, `CA1_ProS`; corresponding
  Seurat subsets are `EX`, `IN`, `CA1_PROS_CANDIDATE`.
- A 25th mandatory result is `results/seurat/original_full/sct_global` with
  `sct_REP_r0p4` and its completed producer receipt. Its nucleus roster is the
  same original-full cohort. `fixed_features.tsv` and
  `SCT_MODELED_FEATURE_COVERAGE.tsv` must agree on the all-twelve modeled residual
  intersection. The comparator exports `normalization_sensitivity_features.tsv`
  and quantitative log-versus-SCT ARI/NMI, contingencies and animal proportions
  at every shared grid resolution. Feature-panel sizes/intersection are explicit:
  when panels differ, the sensitivity combines normalization and modeled features.
- Every metadata table has unique `nucleus_id`, consistent `sample`/`condition`,
  unintegrated and Harmony clusters at primary resolution 0.4.
- Producer receipts are mandatory: Scanpy `COMPLETE` with the correct branch and
  population; Seurat `EXECUTED_REVIEW_REQUIRED` with the correct branch and stage.
  Each result must match `results/shared/BRANCH/membership_POPULATION.tsv` and its
  completed shared receipt, including the exact six WT/six cKO mouse identity.
- Scanpy primary cluster aliases `cluster_unintegrated`/`cluster_harmony` and
  grid `leiden_REP_r0.2`, `r0.4`, `r0.6` are supported.
- Seurat `log_REP_r0p4` or unique subset prefixes are supported. Decimal and
  `p` resolution notation are normalized; conflicting aliases cause failure.
- Scanpy UMAP coordinates are `UMAP_REP_1`/`UMAP_REP_2` in metadata. Seurat
  coordinates are adjacent `umap_unintegrated.tsv` (`UMAPU_1`, `UMAPU_2`) and
  `umap_harmony.tsv` (`UMAPH_1`, `UMAPH_2`), each keyed by `nucleus_id`.

Core comparisons and outputs

- Cross-tool agreement within each count branch and population uses every
  common nucleus at every shared representation/resolution. ARI and arithmetic
  normalized mutual information are label-invariant. Roster differences remain
  visible rather than being hidden by the common-nucleus calculation. Normal mode
  refuses different rosters for same-count cross-tool or paired correction tests;
  full-versus-intersection comparisons explicitly permit the declared restriction.
- Within-tool branch comparisons distinguish paired count correction,
  full-versus-intersection cohort/reclustering changes and their combined effect.
- Exact Hungarian assignment maximizes one-to-one cluster overlap. Cluster
  splits/merges and unmatched clusters remain in the contingency tables.
- `partition_metrics.tsv`, `cluster_matches.tsv`, `contingency_nonzero.tsv` and
  `roster_attrition_per_animal.tsv` contain metrics and their denominators.
- `cluster_proportions_per_animal.tsv` includes zero counts for a represented
  animal/cluster combination and each animal's branch/population denominator.
  These are descriptive proportions, not abundance tests.
- `partition_comparison.pdf` shows primary-resolution Harmony contingencies;
  `umap_comparison.pdf` shows both representations across tools using matched
  colors; `animal_proportions.pdf` uses the same common nuclei per comparison.
- UMAP panels sample at most 20,000 common nuclei deterministically; metrics
  always use all common nuclei. Rotations/reflections and axis scales are
  cosmetic, so coordinates are not tested as an agreement score.
- `COMPARISON_REVIEW.md`, `INPUT_SOURCE_SHA256.tsv` and `COMPLETE.json` record
  source paths, hashes, explicit missing/optional statuses and review flags.

Required independent integer pseudobulk comparison

The stage requires all three completed Python exports at
`results/python_pseudobulk/BRANCH/population_index.tsv` and R stages at
`results/seurat/BRANCH/pseudobulk/POPULATION_INDEX.tsv`. It joins by exact
`family`, `label` and SHA256 of sorted nucleus membership; ordinal population
keys may differ. The intended scope is the four common atlas totals and the
four original-full fresh Scanpy Harmony cluster families, fixed and intersected
with each branch's source nuclei for count sensitivity. The explicit
`source_label_branch=original_full` is preserved; these inference populations
are separate from fresh reclustering of paired branches. Empty intersected
clusters retain full zero count matrices and an empty-membership digest rather
than disappearing. Extra R Seurat/unintegrated families
are recorded as outside that independent Python scope, not silently discarded.

For every matched population, compare actual
`SOURCE_INTEGER_COUNTS_ALL_GENES_ALL12.tsv.gz` (Python) with
`SOURCE_INTEGER_COUNTS_ALL_GENES_ALL12.tsv` (R) on stable `gene_id` and the twelve
sample columns. Symbol formatting is excluded from the test. Counts are parsed
as exact nonnegative integers, without floating point rounding; reordered gene
rows or animal columns are aligned by identity. Duplicate gene IDs, invalid
counts, incorrect mouse rosters or inconsistent membership hashes cause failure.

`exact_pseudobulk_counts.tsv` reports every matched population's complete entry
count, mismatch count and gene coverage. `exact_pseudobulk_unmatched.tsv` reports
all unmatched/inapplicable populations; `exact_pseudobulk_conservation.tsv` binds
actual ALL-matrix per-mouse totals and nucleus counts to both independent source
conservation exports. `exact_pseudobulk_branch_coverage.tsv` records each branch's
expected/matched population coverage. A diagnostic table retains the first 1,000
count disagreements; full mismatch counts remain in the main numerical table.
No exact-aggregation PASS is assigned unless all intended populations in all
three branches match and all actual integer arrays and animal-level conservation
checks agree. This is a computational reproducibility check, not annotation truth.

Optional pseudobulk comparison

Automatically discovered under `results/TOOL/BRANCH/pseudobulk`:

- `POPULATION_INDEX.tsv`: `population_key`, `family`, `label`.
- `POPULATION_MEMBERSHIPS.tsv`: `population_key`, `nucleus_id`, used to bind
  the fitted population to fresh cluster membership when available.
- `POPULATION_KEY/AUDIT/DESeq2/PRIMARY/ALL_RESULTS_STABLE_ID_AUDIT.tsv`:
  `gene_id`, `log2FoldChange`, `padj`.
- `POPULATION_KEY/AUDIT/edgeR/PRIMARY/ALL_RESULTS_STABLE_ID_AUDIT.tsv`:
  `gene_id`, `logFC`, `FDR`.

The R stage can export models for both `SCANPY_*` and `SEURAT_*` fresh cluster
families. The comparator follows the population family to the correct logical
clustering tool even when both models were executed by R in the Seurat folder.
`ATLAS_*` families retain their explicitly shared annotation interpretation.
Paired count sensitivities use verified fixed original-full cluster memberships;
fresh branch cluster results remain separate from those inference memberships.

The native module contract declares cKO-versus-WT effects. Optional custom
inputs can be supplied with `--de-manifest /absolute/new-run/manifest.tsv`.
All custom result paths must remain under the new run's `results` directory.
Required manifest columns are `tool`, `branch`, `family`, `label`, `method`,
`results_path`; also supply `contrast` to enable comparison. Optional columns:
`population`, `cluster_column`, `population_key`, `endpoint`, `gene_id_column`,
`effect_column`, `fdr_column`, `membership_digest`. A membership digest is SHA256
of lexicographically sorted nucleus IDs, one ID plus newline per record.

Only `PRIMARY` endpoints and matching explicitly declared contrasts are paired.
Gene IDs must be unique; gene symbols are never collapsed. Pearson/Spearman
effect correlations, median absolute effect difference, nonzero sign agreement
and FDR overlap are accompanied by finite-gene denominators. Missing/NA effects
and FDRs stay unavailable. No p-values are compared by numerical identity.

DESeq2-versus-edgeR comparisons use the same declared population. Cross-tool or
cross-branch common-atlas/fixed-population comparisons are labeled reproducibility
checks using shared annotation evidence. Fresh cluster populations can be paired
through their common-nucleus Hungarian match only after fitted population
membership is verified against that fresh metadata. By default, a matched pair
needs at least 20 overlapping nuclei and at least 50% reciprocal overlap in the
common comparison cohort (`--min-de-cluster-nuclei` and
`--min-de-cluster-overlap`). These explicit comparability gates do not establish
cell identity. Their different memberships
can affect effect sizes, so such comparisons do not isolate correction alone.

Outputs are `de_effect_metrics.tsv`, compressed complete stable-ID gene pairs
in `de_gene_pairs.tsv.gz`, `optional_de_status.tsv` and up to 16 pages in
`de_effect_comparison.pdf` (`--max-de-plot-pages` changes that display limit).
All available valid comparisons remain in numerical tables. Model-specific BH
families can differ even with a common starting gene universe; nonsignificance
does not establish absence of a biological effect.

Source preparation is not evidence of analysis completion. Actual PDFs require
visual inspection after the parent workflow executes this stage.
