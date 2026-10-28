Takeaway: SATB2 stricter scCODA sampling qualified four fits; twelve remain inconclusive and mixed-cluster effects need caution.

This analysis checked whether stricter numerical sampling could support the requested WT–cKO proportion comparisons without changing the biological model. All sixteen models completed using the same frozen nuclei, counts, labels, references and inferential rules. Four models passed the unchanged sampling checks, while twelve retained divergences and remain inconclusive. Two mixed excitatory categories had numerically qualified relative increases in the shared sequencing partition, but their target identities and references limit biological interpretation. The passing inhibitory models selected no effects, which does not establish equivalence between genotypes. The appropriate next step is to deliver the complete qualified and inconclusive outputs together, without further automatic tuning or relabeling.

## Scientific context and question

The P60 study contains six independent WT and six Satb2 cKO mice. The requested scCODA analysis addresses relative recovered-nucleus proportions using frozen Scanpy and Seurat annotations in broad all-nucleus classes, excitatory clusters, inhibitory clusters and CA1 molecular programs. Separate full-cohort and common-sequencing-partition comparisons were predeclared. Mice are the biological units; nuclei, tools and repeated sampler runs are not independent biological replications. Molecular CA1 programs do not establish anatomical layers.

## Question, lifecycle, and endpoint

The separate same-model stricter-sampler run completed all nineteen native steps with return code zero. Registry: `cf9eb1f5-7165-40be-ad0e-4ca04faa8df4`; native run: `snakemake-satb2-sccoda-composition-20261028-20261028-122900-5660fb24`. Root: `/Volumes/Neurowissenschaften/-Ali_on_server/snr_p60/SATB2_FOLLOWUPS_20261028_sccoda_strict_v1_1_0`.

The initial completed family failed all sampling gates and was separately reviewed and preserved. This single global numerical check used native target acceptance 0.99, 4,000 warmup and 6,000 retained draws for each of four genuine independent chains. Seeds 1337, 2718, 31415 and 65537, depth twelve, priors, likelihood, formula, references, counts, expected-FDR rule and quality gates were unchanged. All five production files are byte-identical. No selection of favorable models or biological labels determined which fits were repeated. No further automatic tuning is planned.

## Key findings

All sixteen models retain their exact initial mouse/count/category/reference/member identities. Full-cohort denominators are 107,747 all-nucleus, 60,908 excitatory, 8,911 inhibitory and 27,703 CA1/ProS nuclei; corresponding common-partition totals are 82,504, 45,966, 6,817 and 22,065. Each scope partitions its own kept nuclei, and overlapping scopes must not be added together. Original exclusions and all inhibitory nuclei remain preserved.

Independent read-only audit verified 672 final output receipts, six bound input identities, sixty-four authentic native chains, 384,000 retained draws and 512 exact chain-array slices. Native PIP, expected Bayesian FDR, effect-table arithmetic and diagnostic decisions matched. Native expected FDR remains 0.05 within each scope/cohort, distinct from BH-adjusted gene-DE results.

| Numerically passing model | Cohort | Maximum split Rhat | Minimum ESS | Selected coefficients |
|---|---|---:|---:|---:|
| Scanpy excitatory clusters | 6 WT / 3 cKO common partition | 1.000837 | 2,555.48 | 1 |
| Seurat excitatory clusters | 6 WT / 3 cKO common partition | 1.001601 | 2,118.87 | 1 |
| Scanpy inhibitory clusters | 6 WT / 3 cKO common partition | 1.001056 | 2,541.37 | 0 |
| Seurat inhibitory clusters | Full 6 WT / 6 cKO | 1.000893 | 2,534.64 | 0 |

All four have zero retained divergences. The other twelve models fail the unchanged zero-divergence criterion, with 81 retained divergences in total. All models meet the declared Rhat and ESS criteria; those diagnostics alone do not override divergences. Scanpy CA1 full-cohort sampling also reached the depth cap on two of 24,000 retained transitions and is already inconclusive. No broad-class or CA1-program model qualifies for inferential use.

The two numerically qualified selections are:

| Saved target | Saved reference | Cohort | Native PIP | Native compositional log2 fold change |
|---|---|---|---:|---:|
| Scanpy excitatory C9, mixed atlas subclasses | C14 | Common partition | 0.993333 | +0.646680 |
| Seurat excitatory C4, mixed atlas subclasses | C21 | Common partition | 0.995958 | +0.625121 |

Their native coefficient values are +0.484993 and +0.467923 respectively; those coefficients and compositional log2 fold changes are different quantities. Positive signs indicate higher relative recovered proportions in cKO under the stated models. Scanpy C14 retains an immune-associated signature-review annotation, and Seurat C21 retains astrocyte/immune-associated signature review. Neither reference is established as invariant, and neither mixed target is a resolved cell identity. Similar coefficient signs do not establish equal cluster membership or independent replication.

Compared with the initial family, retained divergences decreased from 1,770 to 81 and numerical passes increased from zero to four. Native selected coefficients changed from three to two; numerically qualified selections changed from zero to two. This is a diagnostic sensitivity of the same study/model, not additional biological evidence or an attempt to obtain more findings.

## Interpretation

The defensible result is limited inferential support, with complete descriptive proportions and twelve explicitly inconclusive models. The two shared-partition mixed excitatory outputs are exploratory and reference-sensitive. They do not establish a reproducible increase in a named neuronal cell type or absolute tissue abundance. Passing inhibitory fits did not select effects at the native expected-FDR rule; this is not proof that inhibitory proportions are unchanged. No robust broad-cell-class or CA1-layer abundance conclusion follows from failed-gate fits.

Better numerical sampling cannot remove partial genotype/technical-partition confounding, limited common-partition cKO replication, annotation uncertainty or compositional reference dependence. A structural zero reference coefficient is not proof of unchanged reference biology. CA1 expression-program reassignment can reflect RNA regulation rather than physical cell-number changes. Original DESeq2 fits remain the separate evidence for genotype gene expression.

## Artifacts

**Results:** Under `results/`, `prepared/all_mouse_counts_and_proportions.tsv`, the exact prepared composition/roster/reference tables, sixteen `fits/<tool__population>/<cohort>/` folders with native H5AD/NPZ arrays and native/gated effect tables, and `report/all_native_effect_calls.tsv`, `model_statuses.tsv` and full/common-partition descriptive comparisons. Every failed model preserves its native values with inferential calls withheld.

**QC and visual reports:** Per-fit `RESULT_STATUS.json`, `MCMC_VALIDATION.json`, full stochastic-parameter and transition diagnostics, raw stacked/individual-mouse proportion plots and native effect-interval PNG/PDFs. The final seal records four diagnostic passes and twelve inconclusive models. Passing excitatory and inconclusive CA1 effect figures, plus raw mouse-proportion plots, were visually inspected with legible saved labels and explicit gate warnings.

**Provenance:** Executed source/config and runtime/input identities are preserved alongside both sampler families. The final Workbench record is `/Users/mujahid357/Documents/codex_new/SATB2_FOLLOWUPS_20261028/receipts/composition_strict_terminal_record.json`. Independent audit, sixteen-row diagnostics, initial-versus-strict comparisons and figure checks are under `/Users/mujahid357/Documents/codex_new/SATB2_FOLLOWUPS_20261028/reviews/strict_sampler_actual_readonly/terminal_v1_1_0/`.

**Model synthesis:** This review and the independent `FINAL_SCIENTIFIC_AUDIT.md` explain numerical versus biological qualification. The initial family and its failed gates remain separately available, rather than being overwritten or combined as additional replication.

## Limitations and blockers

Twelve models remain inferentially inconclusive under the prespecified standard. Full-cohort and common-partition outcomes, reference uncertainty, mixed identities and relative rather than absolute abundance must remain visible. Native high-acceptance warnings are generic sampling messages; high acceptance alone neither validates nor invalidates a posterior. No gate was relaxed to accommodate these warnings or the remaining divergences. This analysis does not resolve C8 identity, infer anatomical position, demonstrate cell loss, or replace the missing kinetic inputs for RNA velocity.

## Recommended next step

Deliver all native data, plots, uncertainty and both sampler families transparently in the offline/online results reader. Label the two mixed-category outputs as numerically qualified exploratory relative changes, with their flagged references and restricted cohort. Preserve the twelve inconclusive outcomes and stop automatic sampling adjustments. Any new reference-identity or independently replicated biological validation requires a separate scientific question; it should not be substituted for this completed prespecified check.
