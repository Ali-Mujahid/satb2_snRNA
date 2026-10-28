Takeaway: SATB2 scCODA completed all 16 composition models, but sampling diagnostics withhold every inferential call.

This analysis asked whether recovered cell-type and cluster proportions differ between the twelve WT and Satb2 cKO mice. Native scCODA completed all sixteen predeclared tool, population and cohort models with four genuine chains per model. Frozen nucleus membership, count denominators, posterior transport and native statistical arithmetic passed independent checks. Every model nevertheless retained divergent sampler transitions, and one also failed the prespecified mixing and effective-sample criteria. Therefore the raw mouse-level proportion plots are available, but neither a differential-proportion claim nor a conclusion of no difference is validated. One globally specified stricter sampling check of the same model is justified by these numerical diagnostics, while preserving every original result and unchanged biological and inferential settings.

## Scientific context and question

The P60 study contains six independent WT and six Satb2 cKO mice. The requested compositional method was scProportionTest or scCODA; scCODA was selected to retain one multivariate count vector per mouse and the compositional denominator. Frozen Scanpy and Seurat labels are compared in four separate populations. Full-cohort genotype associations are partly confounded with sequencing partition, so a separate common-partition six-WT/three-cKO sensitivity was predeclared. Relative recovered-nucleus composition cannot establish absolute tissue abundance or cell loss.

## Question, lifecycle, and endpoint

The third implementation attempt completed all nineteen native workflow steps with return code zero. Registry: `faf0eb93-87e6-415d-803a-18b4da6b885e`; native run: `snakemake-satb2-sccoda-composition-20261028-20261028-112519-832757e3`. Root: `/Users/mujahid357/Documents/codex_new/SATB2_FOLLOWUPS_20261028/runs/20261028_sccoda_v1_0_2`.

The first two preserved attempts failed at count-layout and native-summary coordinate compatibility checks, respectively; their failed-run reviews were already submitted. Version 1.0.2 corrected only those transport/metadata compatibility issues, with exact count values, sample identities and model preserved. This completed workflow is not scientifically accepted for inferential composition calls because its frozen MCMC gates fail.

Native pertpy 1.4.0 scCODA used four serial independent NUTS chains with seeds 1337, 2718, 31415 and 65537, 1,500 warmup and 3,000 retained draws per chain, target acceptance 0.90 and maximum depth twelve. The formula compares cKO with WT. Native expected Bayesian FDR 0.05 is separate within each scope/cohort, not a frequentist BH gene-DE threshold. All stochastic parameters require split Rhat≤1.01, effective samples≥400 and zero retained divergences before calls are validated.

## Key findings

| Frozen population, each tool | Retained nuclei | Categories |
|---|---:|---|
| All-nucleus broad composition | 107,747 | 8 |
| Excitatory subclusters | 60,908 | Scanpy 19; Seurat 23 |
| Inhibitory subclusters | 8,911 | Scanpy 17; Seurat 20 |
| CA1 molecular programs | 27,703 | 4, including mixed/weak and conflicting/unassigned |

The two tools use identical kept nucleus sets within each corresponding population, with their own frozen labels. Each population uses its own saved exclusions and denominator; overlapping populations must not be added together or silently intersected. The 315 all-nucleus exclusions and 77 excitatory exclusions were retained, while no new inhibitory minimum or downsampling was imposed.

All sixteen models completed, providing sixty-four genuine chains and 192,000 retained draws. Independent review checked 672 final output receipts, six bound input identities, 512 authentic chain-array slices and 206 effect-table rows. Count/sample/category axes, zero handling, native posterior/PIP/expected-FDR arithmetic and diagnostic calculations matched. No data or serialization mismatch explains the numerical warnings.

All sixteen models failed the unchanged zero-divergence gate: 1,770 retained divergent transitions in total, ranging from four to 504 per model. The Seurat CA1 full-cohort model additionally had maximum Rhat 1.01447 and minimum effective samples 296.81. There was no observed retained tree-depth saturation. Three native-selected coefficients are preserved in the tables, but all are withheld; there are **zero validated inferential calls**. This is an inconclusive inference result, not evidence of no compositional change.

References were chosen without genotype using the fixed pooled-dispersion/presence rule and retained for the common-partition sensitivity. Broad composition uses Inhibitory; inhibitory subclusters use C5, annotated 052 Pvalb Gaba; CA1 programs use the superficial-associated program. Excitatory references are Scanpy C14 and Seurat C21, whose saved annotations carry immune/glial-associated signature-review flags. These latter references are not established as biologically invariant and remain an interpretation limitation even if numerical sampling improves.

## Interpretation

The raw mouse-level counts, denominators and descriptive proportion plots are trustworthy transports of the frozen annotations. The current posterior cannot support accepted genotype proportion findings under its prespecified numerical standard. Native divergence flags originate in sampler transitions before serialization and reproduce exactly in the saved native objects/arrays. They are consistent with integration or posterior-exploration difficulty; the audit does not establish a unique causal parameter or prior.

Relaxing the diagnostic standard, selecting favorable clusters or treating a zero validated count as a negative biological conclusion would be inappropriate. Cross-tool cluster numbers are not automatically equivalent identities. CA1 class proportions reflect molecular-program classification, which may change with expression; they are not direct measurements of anatomical layer-cell abundance. Better sampling cannot remove genotype/partition confounding, annotation uncertainty or reference dependence.

## Artifacts

**Results:** `results/prepared/all_mouse_counts_and_proportions.tsv`, `composition.json`, sample roster and reference-candidate tables; sixteen `results/fits/<tool__population>/<cohort>/` directories containing native per-chain/combined H5ADs, posterior arrays, transition fields, native effects and `reviewable_effect_calls.tsv`; and `results/report/all_native_effect_calls.tsv` plus the full/common-partition descriptive comparison. Native selections and withheld status remain distinct.

**QC and visual reports:** Each fit has `RESULT_STATUS.json`, `MCMC_VALIDATION.json`, full stochastic-parameter diagnostics and transport audits, with raw stacked proportions, individual-mouse proportions and clearly marked inconclusive effect-interval PNG/PDF figures. `results/FINAL_VALIDATION.json` records sixteen inconclusive models and zero diagnostic passes. Read-only terminal audit: `/Users/mujahid357/Documents/codex_new/SATB2_FOLLOWUPS_20261028/reviews/composition_readonly_audit/20261028_terminal_v1_0_2/TERMINAL_AUDIT.json`, with model-level diagnostics and native arithmetic evidence beside it.

**Provenance:** Frozen source/config, original sample/annotation input receipts, isolated environment lock/status, native plan and earlier failed-run receipts remain intact. The terminal Workbench record is `/Users/mujahid357/Documents/codex_new/SATB2_FOLLOWUPS_20261028/receipts/composition_v1_0_2_terminal_record.json`. No original DESeq2 fit, annotation, previous analysis or concurrent regulatory run was altered.

**Model synthesis:** This review distinguishes a completed computational workflow and verified descriptive data from unsupported inferential claims. The external diagnostic-only repair proposal records one globally specified next check, rather than effect-driven tuning.

## Limitations and blockers

Current MCMC diagnostics block inferential use of every model. Small common-partition cKO sample size, partial technical confounding, relative rather than absolute abundance, uncertain labels and excitatory reference flags remain. The native reference coefficient is structurally zero, not proof that its population is unchanged. Separate scope/cohort Bayesian-FDR families are exploratory and do not provide a global confirmatory FDR guarantee. Completed preparation or tiny synthetic software tests do not prove adequate real-data sampling.

## Recommended next step

Use one separate immutable native plan for the same sixteen models with target acceptance 0.99, 4,000 warmup and 6,000 retained draws per chain. Preserve the same four seeds, depth twelve, code, counts, categories, references, priors, formula, expected-FDR rule and quality gates; assess success by diagnostics, not by obtaining more effects. The production scripts for the prepared candidate are byte-identical. Use two independent one-thread fit processes with chains serial and a fresh server-share output directory to preserve local disk reserve. Retain both runs transparently. If that single stricter check also fails, report the affected models as inconclusive without further automatic tuning or method substitution.
