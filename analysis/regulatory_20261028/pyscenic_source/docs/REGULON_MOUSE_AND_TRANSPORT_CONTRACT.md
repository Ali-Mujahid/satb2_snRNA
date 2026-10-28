# Native positive motif regulon scores, mouse associations, and atlas transport

This is an additive pySCENIC branch. It does not alter earlier hdWGCNA, CollecTRI, DESeq2 or strict-DEG LIANA outputs. Execution of scientific scripts is reserved to the reviewed immutable native workflow. Before launch, only synthetic scores may exercise these scripts.

## Shared immutable configuration

Schema `satb2-pyscenic/v1`; absolute `run_root`, `analysis_id`, seed1337. Bound records are `{path,bytes,sha256}`. Required inputs are:

- `samples`: original twelve-mouse roster with sample/condition columns.
- `prepared_complete`: validated HD preparation COMPLETE receipt.
- `lineages.CA1` and `.IN`: raw integer `counts`, `genes`, `cells` receipts from that completed preparation; cells contain exact `nucleus_id,endpoint,sample,condition` and partition metadata.
- `references.tf_names`, `references.motif_annotations`, `references.rankings` (the two matched mouse mm10 gene-based v10-cluster ranking files). The actual current official TF list is audited, including Satb2 presence; no TF or regulon is forced into the result.
- `additional_hash_inputs.frozen_membership_audit` and `.mouse_coverage`: bound validated preparation cohort records. Every other supplied record is also revalidated at packaging.

Optional `mouse_associations` overrides allow only recognized policy keys. Protocol fixes six WT/six cKO, six WT/three cKO in exact partition `23HJT3LT4_lanes6_7`, 200 genotype-stratified mouse bootstraps, target minimum10, and no AUCell thresholding. Plot limits affect display only.

## Native producer output: `results/pyscenic_native/`

Producer `scripts/run_pyscenic_native.py` performs genuine GRNBoost2, positive-correlation module construction, matching cisTarget motif pruning, native df2regulons and nucleus-level AUCell. Full genotype-neutral eligible RNA ranks are used, with no DEG-only, genotype-enriched, or deliberately selected target universe. Native ctx/motif enrichment, regulator candidates, target weights, raw nucleus AUC, exact metadata, parameters and source/input receipts remain available.

Regulon IDs are lineage-specific, e.g. `CA1::Satb2(+)`; `native_name` retains unmodified native `Satb2(+)`. CA1 and IN definitions must not be joined because they share a TF symbol. Only positive motif-supported native regulons enter this branch. Native gene2weight usually inherits the maximum surviving GRNBoost2 importance; it is not a binding affinity, causal effect, or interaction strength. Native pySCENIC adds the TF itself to a regulon's gene2weight map; this selfgene is preserved and flagged, not advertised as self-binding evidence.

| Native table | Required columns / meaning |
|---|---|
| `REGULON_ACTIVITY_PER_MOUSE.tsv` | `lineage,endpoint,endpoint_id,source_endpoint_id,regulon_id,native_name,TF,n_targets,n_background_genes,eligible,sample,condition,n_nuclei,aucell_mean,score_unit` |
| `REGULON_COVERAGE.tsv` | same identity/target/background/eligible columns, plus `reason`; every discovered positive native regulon appears for every endpoint of its lineage |
| `REGULON_TARGETS.tsv` | `lineage,regulon_id,TF,target,weight`; native target members/weights without re-ranking; optional `self_target` must match TF==target |
| `LINEAGE_STATUS.tsv` | one CA1/IN row each: `lineage,status,n_regulons,n_nuclei,n_genes`; status `COMPLETE_WITH_REGULONS` or `EMPTY_NO_SUPPORTED_REGULONS` |
| `SATB2_NO_SELF_TARGET_ACTIVITY_PER_MOUSE.tsv` | main score columns plus `network_variant`; supported no-self scores only |
| `SATB2_NO_SELF_TARGET_COVERAGE.tsv` | coverage columns plus `network_variant`; all five exact endpoint states explicitly, including unsupported |

The score unit is exactly `native_AUCell_AUC`. Each eligible endpoint–regulon has twelve finite independent-mouse means, six of each genotype; each mean comes from all exact selected nuclei in that endpoint/mouse. No AUC threshold/binarization, nucleus downsampling, filtering by P value, or pseudo-replicated nucleus test is performed. All native positive sets and targets are retained, but a set with fewer than ten measured targets has `eligible=FALSE`, a reason, and no manufactured score rows. Completed native discovery can legitimately yield no supported regulons; empty tables remain header-complete and that lineage has an explicit unavailable activity state.

Raw matrix `<lineage>/aucell/nucleus_AUC.tsv.gz` has `nucleus_id` as its first column followed by lineage-distinct eligible regulon IDs. Matching `<lineage>/aucell/nucleus_axis_metadata.tsv.gz` retains exact ordered `nucleus_id,endpoint,sample,condition,lineage,sequencing_partition`. Even an empty-regulon lineage retains all nucleus IDs with a nucleus-ID-only matrix. PROVENANCE identifies both paths/receipts and the native mapping. A complete, exact full nucleus axis is required; cluster means are never substituted for nucleus measurements.

`COMPLETE.json` schema is `satb2-pyscenic-native-complete/v1`, status COMPLETE, analysis_id, same config receipt, single producer source receipt, inputs, and `output_files` mapping every relative native path to absolute `{path,bytes,sha256}`. Completion is written last; no unsealed files may silently enter the atlas/archive. Required true checks are `positive_motif_regulons_only`, `full_genotype_neutral_RNA_background`, `no_DEG_selection`, `native_nucleus_AUCell`, `exact_nucleus_axes_verified`, `all12_mouse_roster`, `min_measured_targets_10`, `Satb2_self_target_sensitivity_separate`, `inputs_config_source_SHA256_unchanged`.

## Mouse modeling: `scripts/analyze_regulon_mouse.R`

CLI uses absolute `--config`, `--input-dir` (completed native producer folder), and `--outdir` (new or strictly empty folder, normally `results/regulon_mouse`). A nonempty output, symlink or file is rejected. Configuration, producer completion, source/input receipts and every native output are verified; no partial native scores are accepted. Source import is inert.

For each eligible lineage-distinct regulon × matching exact endpoint, the unweighted continuous model is `mean_native_nucleus_AUCell_AUC ~ condition`, reference WT. Its effect is the cKO−WT mean difference in native AUCell AUC units. Exactly six independent mice per genotype are required. No count model, DESeq2, edgeR, new gene DE, or nucleus/metacell pseudoreplication is used. Constant outcomes and numerically zero within-group residual variance remain explicit untested statuses with NA P/FDR/CI; descriptive means/differences may remain.

Primary association FDR is BH jointly across **every tested regulon × matching endpoint, both lineages**. The same-partition six-WT/three-cKO model has a separate joint BH family. Leave-one-mouse-out models summarize twelve omitted-mouse effects; 200 mouse bootstrap resamples are drawn independently within each genotype and give exploratory percentile95% intervals. These bootstrap/sign summaries are not additional P values or independent replication. Genotype partially confounds sequencing partition, so primary and sensitivity effects remain observational.

Outputs are `REGULON_ASSOCIATIONS.tsv`, `REGULON_SAME_PARTITION_ASSOCIATIONS.tsv`, `REGULON_LOO_MODELS.tsv`, `REGULON_LOO_STABILITY.tsv`, `REGULON_MOUSE_BOOTSTRAP.tsv`, `REGULON_MODEL_ELIGIBILITY.tsv`, `REGULON_MOUSE_SCORES.tsv`, `MOUSE_PARTITION_AUDIT.tsv`, `MODEL_SUMMARY.tsv`, `SATB2_STATUS.tsv` (explicit native/no-self available or absent status for all five endpoints), plotting selection, provenance and `sessionInfo.txt`. Effect records retain native definition/TF/lineage/endpoint IDs, all target/background/coverage counts, group means, cKO−WT difference, SE/statistic/P/CI, association_fdr/padj, tested status and named FDR family. No native AUC significance score is invented; AUCell AUC is not an enrichment P value.

Nine bounded views are saved as PDF/SVG/PNG: top30 post hoc effect heatmap, primary and same-partition volcanoes (gray below-threshold points), five endpoint-specific top-six post hoc mouse box/dot plots, and the Satb2 full/no-self association comparison. Each point is an independent mouse. Missing activities receive no fabricated zero. Export uses explicit Cairo PDF/base SVG devices, avoiding the earlier optional svglite dependency.

`COMPLETE.json` schema `satb2-regulon-mouse-complete/v1` seals all output bytes after source/input rechecks. True checks include all12 roster, native means preserved, both-lineage global BH, separate sensitivity BH, no activity thresholding, no nucleus pseudoreplication, no new DESeq2, positive-only support, and source/input stability. Partial failures have no complete seal.

## Satb2 selfgene sensitivity

If a genuine native Satb2 positive motif regulon exists and includes Satb2 itself, the native producer separately scores that exact regulon minus its selfgene using unchanged full expression rankings/background, provided at least ten targets remain. This sensitivity does not change the primary full native regulon. ID is `<lineage>::Satb2(+)_NO_SELF_TARGET`, native_name stays `Satb2(+)`. Absent native regulon, absent selfgene, or fewer than ten remaining targets produce explicit unavailable coverage, not an RNA proxy or zero AUC. Coverage states exist for all five endpoints even if scores are absent.

Corresponding mouse outputs use `SATB2_NO_SELF_TARGET_` prefixes and a separate available-five-endpoint BH family, plus its own same-partition BH, LOO and bootstrap tables. This is a sensitivity to native membership bookkeeping, not independent binding/functional validation. All native regulons are discovered separately per lineage; a CA1 Satb2 set need not resemble an inhibitory Satb2 set.

## Transport: `scripts/package_pyscenic_results.py`

CLI `--config ABS --results ABS_RUN/results`. The packager performs no GRN inference, motif pruning, scoring, averaging, normalization, model, P value or FDR adjustment. It independently validates producer/mouse completion, hashes, matching source/config, frozen endpoint counts/membership hashes, twelve-mouse roster, exact native-to-R means, native target weights, per-nucleus score/membership axes, min10 eligibility, unique keys, and the separate no-self contract. It rejects modified/unsealed data and nonempty atlas/previous FINAL_VALIDATION. It revalidates input/native/source bytes before sealing.

Transport schema `satb2-regulatory-atlas/v1` has separate `pyscenic_motif_CA1` and `pyscenic_motif_IN` methods, family `motif_regulon`, activity_unit mouse. Features, raw values, native primary effects, full native weighted target networks, eligibility coverage and download receipts preserve lineage-specific definitions. Native target edges are candidate TF→target orientation with GRNBoost2-derived gene2weight; arrows are not proof of binding or causality. Selfgene loops retain a warning. Unmodelable native sets retain explicit null values, no zeros. No activity UMAP is packaged without independently bound actual original coordinates and exact nucleus alignment; this branch does not require one. Native per-nucleus scores remain downloadable.

Writes `atlas/data.json`, deterministic gzip, data.js, downloads.json and manifest.json. `FINAL_VALIDATION.json` schema `satb2-pyscenic-regulatory-final-validation/v1`, status COMPLETE, seals all native/transport receipts and exact axes. COMPLETE describes successful execution/validation, not a positive scientific result or proof of function. Earlier pySCENIC-reference-blocked reports and published v3.7 outputs remain historical; a new additive release integrates this component after actual result inspection.

Synthetic checks: `tests/test_regulon_mouse.R` verifies score identities/coverage, model units, constants, min10, global/separate BH, LOO/200 bootstrap, no-self unsupported states and truly empty lineages. `tests/test_package_pyscenic_synthetic.py` exercises a wholly synthetic producer→R→transport and rejects changed input references, native scores, raw nucleus axes and downstream means. It never reads study expression or runs GRN inference.
