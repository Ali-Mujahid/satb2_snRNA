Takeaway: Strict scCODA finished; four fits pass sampling gates, while twelve remain inconclusive and references limit interpretation.

The stricter composition sensitivity finished the same planned models for recovered nucleus proportions between wild-type and Satb2 knockout mice. Independent read-only checks verified the frozen counts, annotations, references, native posterior arrays and selection arithmetic. Four fits pass the unchanged sampling gates, while twelve retain divergent transitions and remain inconclusive. Two mixed excitatory categories receive native selections in the shared sequencing partition, relative to references that retain glial-associated identity warnings. This is a within-study, reference-sensitive result and does not establish a robust cell-type abundance change. No further sampling or tuning is proposed, and the original completed analysis is preserved separately.

## Scientific context and question

This additive mouse snRNA-seq analysis asks whether relative recovery of annotated nucleus populations differs between six independent WT and six independent Satb2 cKO mice. The four scopes are all-nucleus broad categories, excitatory subclusters, inhibitory subclusters and all four CA1 molecular expression programs, assessed separately using saved Scanpy and Seurat labels. Independent mice are the biological units; nuclei are not independent condition replicates. The full cohort has 6 WT / 6 cKO mice; the shared sequencing partition has 6 WT / 3 cKO mice.

## Question, lifecycle, and endpoint

The durable native run was observed completed with return code zero at registry `cf9eb1f5-7165-40be-ad0e-4ca04faa8df4`, run `snakemake-satb2-sccoda-composition-20261028-20261028-122900-5660fb24`, on the registered local target with server-share outputs. Its final seal records all 16 models and `COMPLETE_WITH_INCONCLUSIVE_OR_UNTESTABLE_SCOPES`: 4 diagnostic passes and 12 inconclusive models. This is terminal artifact completion with partial inferential validity, not validation of the whole family. Root owns the agent-authored terminal Workbench submission; this independent audit supplies supporting evidence.

## Key findings

- Every one of 672 final scientific output receipts and all 6 bound input hashes was rehashed successfully. Exact mouse rosters, counts, category order, membership, references, priors/contrast and inference gates agree with the initial frozen family. The five production source files are byte-identical to the initial workflow. Only the authorized global sampler settings differ: target acceptance 0.99, 4,000 warmup and 6,000 retained draws per chain; the four saved seeds are 1337, 2718, 31415 and 65537, with depth 12.
- All 64 genuine native chains and 384,000 retained draws are accounted for; 512 saved chain-array slices exactly match their corresponding pooled arrays. Native posterior inclusion probabilities, expected Bayesian FDR, effect values and MCMC diagnostics independently agree with the saved tables and gate decisions. The report contains 16 model rows and 206 effect rows. No model, sampler, native resummary or scoring API was invoked by this audit.
- Four models pass the unchanged gates of maximum split Rhat ≤ 1.01, minimum effective sample size ≥ 400 and zero retained divergences: Scanpy excitatory shared partition, Scanpy inhibitory shared partition, Seurat excitatory shared partition and Seurat inhibitory all12. Across all 16 models, maximum split Rhat is 1.002041 and minimum effective sample size is 2,118.87. The other 12 fail the zero-divergence gate, with 81 retained divergences in total. Scanpy CA1 all12 also reaches the saved depth cap on 2 of its 24,000 retained transitions (maximum leapfrog steps 4,095); it is already inconclusive. No passing model reaches that cap.
- The same-model sensitivity reduces retained divergences from 1,770 to 81 and increases numerical passes from 0 to 4. Native selected coefficients change from 3 to 2; qualified coefficients change from 0 to 2. These numbers assess diagnostic behavior, not a reason to prefer a run for yielding discoveries. The initial family remains separately available with its original inconclusive status.
- Two native coefficients survive the numerical gate, both within the shared partition: Scanpy excitatory C9, labeled `Excitatory; mixed atlas subclasses`, has compositional log2 fold change +0.646680 and PIP 0.993333; Seurat excitatory C4, also mixed atlas subclasses, has +0.625121 and PIP 0.995958. They are relative changes against Scanpy C14 and Seurat C21, respectively. The positive signs mean higher relative recovered proportions in cKO under those models; they do not establish absolute increases or resolved identities. No expected-FDR rounding-above-target flag was found.

|Tool / scope|Cohort|Divergences|Max split Rhat|Min ESS|Numerical gate|Native / qualified selections|
|---|---|---:|---:|---:|---|---:|
|scanpy__allcells_coarse|all12|11|1.001235|3731.42|INCONCLUSIVE|0 / 0|
|scanpy__allcells_coarse|shared_partition|1|1.000605|4419.28|INCONCLUSIVE|0 / 0|
|scanpy__ca1_layer_programs|all12|11|1.001196|4001.42|INCONCLUSIVE|0 / 0|
|scanpy__ca1_layer_programs|shared_partition|3|1.000618|5776.34|INCONCLUSIVE|0 / 0|
|scanpy__excitatory_clusters|all12|2|1.001509|3981.94|INCONCLUSIVE|0 / 0|
|scanpy__excitatory_clusters|shared_partition|0|1.000837|2555.48|PASS|1 / 1|
|scanpy__inhibitory_clusters|all12|5|1.000936|2668.04|INCONCLUSIVE|0 / 0|
|scanpy__inhibitory_clusters|shared_partition|0|1.001056|2541.37|PASS|0 / 0|
|seurat__allcells_coarse|all12|13|1.000667|4299.50|INCONCLUSIVE|0 / 0|
|seurat__allcells_coarse|shared_partition|18|1.000744|4330.52|INCONCLUSIVE|0 / 0|
|seurat__ca1_layer_programs|all12|3|1.001046|4053.99|INCONCLUSIVE|0 / 0|
|seurat__ca1_layer_programs|shared_partition|7|1.000783|4420.71|INCONCLUSIVE|0 / 0|
|seurat__excitatory_clusters|all12|1|1.001380|2871.87|INCONCLUSIVE|0 / 0|
|seurat__excitatory_clusters|shared_partition|0|1.001601|2118.87|PASS|1 / 1|
|seurat__inhibitory_clusters|all12|0|1.000893|2534.64|PASS|0 / 0|
|seurat__inhibitory_clusters|shared_partition|6|1.002041|2345.90|INCONCLUSIVE|0 / 0|

Full-cohort denominators are identical between tools: 107,747 broad nuclei, 60,908 excitatory, 8,911 inhibitory and 27,703 CA1/ProS. Shared-partition denominators are 82,504 / 45,966 / 6,817 / 22,065, respectively. All 16 exact family identities agree with the approved identity receipt and original preparation.

## Interpretation

Stricter global sampler settings improved diagnostics without resolving every model. A correct receipt, high average acceptance, satisfactory Rhat or satisfactory ESS cannot override retained divergences. All outputs from failed models remain visible for audit while inferential calls stay withheld. The two numerically qualified mixed excitatory outputs are reference-sensitive exploratory associations, not resolved cell-type findings. Their apparent agreement across tools does not prove matching population membership and is not independent replication. The passing inhibitory models have no selected coefficients; this is not proof that inhibitory proportions are unchanged. There is no qualified CA1 layer-program or broad cell-class genotype association in this family. Native PIP and expected Bayesian FDR are composition-model quantities, separate from BH-adjusted gene DESeq2 results; these scopes and tools remain separate exploratory FDR families.

## Artifacts

- Native scientific outputs: `/Volumes/Neurowissenschaften/-Ali_on_server/snr_p60/SATB2_FOLLOWUPS_20261028_sccoda_strict_v1_1_0/results`, including `FINAL_VALIDATION.json`, exact prepared counts, all native model H5AD/NPZ draws and transition fields, effect/call tables, diagnostics, per-mouse proportions and original PNG/PDF figures.
- Independent data/QC audits: `TERMINAL_AUDIT.json`, `NATIVE_ARITHMETIC_AUDIT.json`, `APPROVED_IDENTITY_COMPARISON.json` and `MODEL_DIAGNOSTICS.tsv` beside this review. The separate `INITIAL_VS_STRICT_DIAGNOSTICS.json` and `.tsv` compare saved diagnostics only. `PLOT_SPOTCHECK.json` records hashes, sizes, PDF pages and PNG dimensions for eight original plot files. The inspected excitatory plot retains the mixed target and flagged reference labels; the inspected failed CA1 plot explicitly withholds inferential calls. No figures or underlying data were transformed.
- Execution/source provenance: `/Users/mujahid357/Documents/codex_new/SATB2_FOLLOWUPS_20261028/composition_sampler_strict_source_v1_1_0`, the executed frozen config, runtime-contract/lock receipts and exact initial-versus-strict family identity proof. The independent wrapper is `/Users/mujahid357/Documents/codex_new/SATB2_FOLLOWUPS_20261028/reviews/strict_sampler_actual_readonly/audit_strict_terminal_family.py`; it delegates numeric read-only checks to the preserved pure NumPy/H5 auditor. Historical generic geometry prose inside that reused raw audit is not evidence that a passing strict fit failed.

## Limitations and blockers

Twelve sampling-gate failures block interpretation of those models. Full-cohort genotype association is partially confounded with sequencing partition; the shared-partition sensitivity has only three cKO mice. The reference choices were genotype blind but are not established biologically invariant. Broad models use Inhibitory, inhibitory models use C5 Pvalb, and CA1 models use the superficial-associated molecular program. Excitatory Scanpy C14 retains a Microglia/Other_immune-associated signature-review label; Seurat C21 retains Astrocyte/Microglia/Other_immune-associated signature review. Numerical convergence cannot resolve these reference identity issues or the mixed target labels. No anatomical layer localization, tissue cell count, cell death or causal genotype effect is demonstrated by molecular-program/recovered-nucleus proportions. No cross-tool cluster equivalence was inferred merely from a similar coefficient sign.

## Recommended next step

Preserve both original and strict terminal families with exact native artifacts and their distinct diagnostic status. Present the four numerical passes, the twelve inconclusive models and the two reference-sensitive mixed-category outputs explicitly in the reader; do not promote them to robust resolved cell-type changes. Root should complete the native terminal review before packaging this component, then verify every publication receipt and downloadable native file. Do not refit, tune thresholds, change references, relabel categories or modify prior analyses as part of this completed sensitivity. Any future independently designed biological validation or reference-identity study requires separate scope rather than an automatic search for additional selections.
