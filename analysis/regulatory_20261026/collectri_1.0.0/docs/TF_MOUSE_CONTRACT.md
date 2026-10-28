# Continuous independent-mouse TF score associations

`scripts/analyze_tf_mouse.R` consumes a **completed, receipt-verified** native CollecTRI ULM activity analysis. It does not run gene differential expression, change original DESeq2 fits, infer motifs/regulons, or substitute TF RNA expression for TF target-set activity. Before launch, validation used synthetic scores only.

## CLI and inputs

```sh
/absolute/path/environment/run_R.sh /absolute/path/analyze_tf_mouse.R \
  --config /absolute/path/same_native_activity_config.json \
  --input-dir /absolute/path/results/tf_activity \
  --outdir /absolute/path/results/tf_mouse
```

The output directory may be nonexistent or strictly empty (to accommodate Snakemake parent preparation); nonempty directories, files and symlinks are rejected. The same native `satb2-tf-activity/v1` config is used. `samples` and optional bound mouse-coverage metadata supply each mouse's exact partition. Native output checksums, source, configuration and producer completion are verified before use and again after the analysis.

The primary input is `TF_ACTIVITY_PER_MOUSE.tsv`, filtered by its exact producer mode `mouse_log1p_original_sizefactor_expression`. Native `TF_TARGET_COVERAGE.tsv` defines the all-TF universe and target coverage. Every modeled endpoint/TF must contain exactly the12 frozen mice, six WT/six cKO, finite native ULM t-values and native gene-set FDR, with ≥5 eligible targets. Unknown/duplicated/missing mice, conflicting genotypes or coverage, wrong score units and unreceipted files fail validation. No native gene-set FDR filter is applied.

Optional policy overrides are in `mouse_associations`; unrecognized policy keys are rejected. The random seed is taken from `mouse_associations.seed`, then top-level `random_seed`/`seed`, otherwise1337, and is recorded. All model/plot policies are captured in output provenance. No study data are processed merely by sourcing the script.

## Models and eligibility

Each TF × endpoint gets an unweighted continuous linear model `escore ~ condition`, reference WT. The coefficient is the mean cKO−WT difference in native ULM t-value units; all inference uses independent mice. Exactly6WT/6cKO are required for a primary test. There is no nucleus/metacell pseudoreplication and no count-model/DESeq2/edgeR analysis.

Constant outcomes are not tested. Specifically, score range≤1e−12×max(1,max|score|) is numerically constant; within-genotype residual sum of squares≤n×[1e−12×max(1,max|score|)]² is numerically zero residual variance. These cases receive explicit model statuses and NA statistic/P/CI/FDR. Group means and a descriptive mean difference may remain available. There must be an estimable rank2 linear model with positive residual variance; all modeled outcomes have all12 mice. Unmodeled TFs with <5 targets retain explicit coverage/eligibility rows and receive no manufactured activity score.

The primary **association FDR** is BH across every tested TF × all five endpoints jointly. `association_fdr` and its exact alias `padj` are distinct from `gene_set_fdr`. The latter is the native within-mouse TF target-set enrichment FDR; it is preserved in raw scores and never used for eligibility or genotype significance. The raw producer field `genotype_association_tested=FALSE`, when present, describes the upstream producer only; downstream model tables explicitly record each performed association.

The same-partition sensitivity selects only exact partition `23HJT3LT4_lanes6_7`; it requires6WT/3cKO and has a separate global BH family. Missing or unexpected coverage is explicitly not tested. Genotype partly confounds technical partition; neither the full nor sensitivity model establishes a causally isolated SATB2 effect.

For each tested primary contrast, leave-one-mouse-out models summarize effect sign/range after omitting each of12 mice (at least3 mice/genotype still required). Two hundred independent-mouse bootstrap resamples are drawn **within each genotype**; their percentile95% effect interval is exploratory and not multiplicity-adjusted. Bootstrap sign fractions are not P-values or independent replications.

## Satb2 self-target sensitivity

`SATB2_NO_SELF_TARGET_ACTIVITY_PER_MOUSE.tsv` and `SATB2_NO_SELF_TARGET_COVERAGE.tsv` supply a separate native score variant: same eligible RNA background, only Satb2→Satb2 removed from the prior. Its exact expression mode is `mouse_log1p_original_sizefactor_expression_no_Satb2_self_target`, TF=`Satb2`, network_variant=`SATB2_only_self_target_excluded`. All five endpoints must remain in coverage, including explicitly unavailable endpoints. Header-complete empty score tables are supported. This variant is required by default in the approved protocol.

No-self-target continuous associations use their own five-endpoint BH family, with separate same-partition, leave-one-out and bootstrap outputs. Its native target-set FDR is based on the single Satb2 set and must not be confused with primary across-TF enrichment FDR or genotype-association FDR. Full-prior scores remain unchanged. Target-set scoring with/without a self-target is not independent proof of SATB2 binding or functional activity.

## Output tables

| File | Content |
| --- | --- |
| `TF_ASSOCIATIONS.tsv` | all prior endpoint/TF rows, exact eligibility/status and primary coefficients/P/CI/global association FDR |
| `TF_SAME_PARTITION_ASSOCIATIONS.tsv` | separately adjusted same-partition sensitivity |
| `TF_MODEL_ELIGIBILITY.tsv` | modeled/unmodeled/constant statuses, target/background counts, n_scored_mice, n_expected_mice, n_WT/n_cKO, tested |
| `TF_LOO_MODELS.tsv`, `TF_LOO_STABILITY.tsv` | individual omitted-mouse fits and effect sign/range summaries |
| `TF_MOUSE_BOOTSTRAP.tsv` |200 within-genotype mouse-resampling summaries |
| `TF_MOUSE_SCORES.tsv` | unchanged native long scores/FDR plus exact sample partition |
| `SATB2_NO_SELF_TARGET_*` | corresponding independently named sensitivity results and raw scores |
| `MOUSE_PARTITION_AUDIT.tsv` | the12 frozen samples, genotype and exact partition |
| `TF_PLOT_SELECTION.tsv`, `SATB2_STATUS.tsv` | transparent plotting selection and available/unavailable Satb2 activity |
| `MODEL_SUMMARY.tsv` | tested and association-FDR<0.05 counts by endpoint |
| `PROVENANCE.json`, `sessionInfo.txt`, `COMPLETE.json` | methods, policy, all input/source/software records and output SHA256; completion written last |

Association columns include `endpoint`, `endpoint_id`, `source_endpoint_id`, `TF`, `network_variant`, `n_used_targets`, `n_background_genes`, `n_scored_mice`, `n_expected_mice`, `n_WT`, `n_cKO`, `mean_WT`, `mean_cKO`, `effect_cKO_minus_WT`, `SE`, `statistic`, `pvalue`, `association_fdr`, `padj`, `CI_low`, `CI_high`, `residual_df`, `tested`, `model_status`, `fdr_family`, `score_unit` and `gene_set_fdr_used_for_filtering=FALSE`.

`COMPLETE.json` schema is `satb2-tf-mouse-complete/v1`, with status COMPLETE, input/config/source `{path,bytes,sha256}` records and `outputs` mapping each relative path to `{bytes,sha256}`. Partial/failed results have no completion receipt.

## Bounded figures

Nine plot views (each PDF/SVG/PNG): top30 TF-effect heatmap, full-prior and same-partition association volcanoes, five endpoint-specific box/dot panels (prespecified Satb2 plus five post hoc top TFs), and one full-prior/no-self Satb2 association comparison. Volcanoes use light-gray nonsignificant points, cKO higher/lower colors, native score differences and **association** FDR axes. Unavailable/untested heatmap cells are gray; no zero score is fabricated. Selection by association FDR/P/effect is explicitly post hoc and exploratory. Native score magnitudes depend on endpoint-specific target/background coverage; the plots do not establish directly comparable absolute TF activity between cell populations.

Synthetic tests are in `tests/test_tf_mouse.R`. They verify mean-difference sign/units, constant outcomes, exact12-mouse coverage, separate BH families,6WT/3cKO sensitivity, leave-one-out/bootstrap structure, unavailable no-self targets and independence from native gene-set FDR filtering. PDF skill execution bookkeeping/visual review is performed by the parent workflow when actual figures are authored.
