# Independent read-only review of actual native Monocle3 results

Reviewed 28 October 2026, after both native lineage seals and `results/FINAL_VALIDATION.json` existed. Registry: `081614e2-5ca8-48a8-9b33-b431673f0eee`. Native run: `/Users/mujahid357/Documents/codex_new/SATB2_FOLLOWUPS_20261028/runs/monocle3_20261028_v1`. This review concerns the actual study outputs, not synthetic fixture results. Its scope is an independent source, data, native-object, numerical and figure audit; the parent agent retains responsibility for the Workbench terminal review and publication.

**Verdict: PASS for the implemented exploratory graph/module analysis, with explicit scientific interpretation limits.** All four native scope models completed. There is no independently established biological temporal root, no demonstrated developmental direction, and no genotype differential-expression test in this Monocle3 analysis.

## Independent evidence

The read-only R checker verified the actual runtime-lock closure, source/config/input/output receipts, retained sparse integer counts against the original frozen matrices, exact cell/feature axes, all mouse/endpoint axes, original/prepared/retained denominator reconciliation, native constructor size factors, saved full native objects and indexes, principal graph exports, rooted-vector exports, native graph-test BH values and module candidate/loading/assignment coverage. It called no fitting, graph-test, ordering, aggregation or scientific workflow. Arithmetic comparisons allowed ordinary TSV floating-point roundoff only. Neither the run tree, scientific source, environments, inputs nor previous analyses was modified.

Both lineage audits passed: 48 CA1 plus 126 inhibitory output hashes (174 total); both native lineage seals matched the final validation. All 32 final native PDFs parsed as nonempty one-page documents. Four representative actual downloadable PDFs were rasterized and visually checked (CA1 module heatmap and gene overlays; inhibitory C8 condition graph and C1 module heatmap). The corresponding native PNGs and other selected PNGs were also inspected. Axes, legends, values and titles were readable. Poppler emitted a fontconfig setup warning but successfully rendered all four inspected PDF pages; native upstream ggplot2 aesthetic deprecations were cosmetic warnings and did not invalidate scientific results.

## Actual observed scope results

| Native scope | Nuclei | Independent mice | Retained features tested | Graph components / partitions | Native q ≤ 0.05 graph-associated genes | Native modules |
|---|---:|---:|---:|---:|---:|---:|
| CA1 DL + SL | 20,456 | 12 (6 WT, 6 cKO) | 10,618 | 2 / 2 | 9,056 | 47 |
| Inhibitory C1 | 1,718 | 12 (6 WT, 6 cKO) | 10,539 | 1 / 1 | 3,669 | 19 |
| Inhibitory C6 | 501 | 12 (6 WT, 6 cKO) | 9,849 | 1 / 1 | 1,439 | 8 |
| Inhibitory C8 | 468 | 12 (6 WT, 6 cKO) | 10,258 | 1 / 1 | 2,002 | 11 |

Every tested native graph-test row had successful finite output. Every selected graph-associated gene was present on the PCA loading axis and assigned to a native module. All selected nuclei had finite distance from the scope/component-specific arbitrary reference; no infinite/unreachable values were substituted as zero. The primary saved CDSs remain unrooted, and biological-pseudotime exports remain unavailable. Inhibitory identities were modeled separately; no CA1–inhibitory or C1–C6–C8 trajectory was forced.

The prepared universe is 32,245 symbol-summed feature keys derived from the original 32,285 stable-ID gene universe. Eligibility was genotype blind (at least 5% detection in any included frozen endpoint and positive raw variance), with no DEG or HVG selection. The native model uses size factors on the retained eligible genes, not original full-library counts. These matched retained nucleus totals divided by their geometric mean within each scope to below 6 × 10⁻¹⁵ absolute arithmetic error. Retained/prepared count fractions ranged from 0.9701–0.9967 in CA1, 0.9664–0.9935 in C1, 0.9460–0.9882 in C6 and 0.9746–0.9888 in C8. All three library denominators are exported; they must not be described as identical normalizations.

## Observed CA1 structure and genes

The large component contains 19,609 nuclei (8,748 DL; 10,861 SL), with 10,186 WT and 9,423 cKO nuclei. The second component contains 847 nuclei (845 SL; 2 DL), with 577 WT and 270 cKO nuclei. Both components contain all 12 mice. This establishes observed molecular separation with mouse coverage, not a validated developmental branch, anatomical location or new cell type. The condition-colored graph shows substantial genotype-associated geometry within the main component; partly confounded sequencing partition and other technical/molecular differences remain alternative contributors.

The highest native graph-autocorrelation genes include Chrm3 (Moran I 0.7254), Zbtb20 (0.6991), C1qtnf7 (0.6971), Nrxn3 (0.6950), Mctp1 (0.6872), Spock1 (0.6768), Galntl6 (0.6720), Epha6 (0.6664), Cadm2 (0.6614) and Kcnh7 (0.6317). These are graph-associated genes, not newly established WT–cKO DEGs or time-dependent genes. In the native overlay, C1qtnf7 is prominent in the small component; that alone does not resolve its identity.

The 47 module sizes range from 90 to 316 genes. Examples ranked by graph Moran I include module 9 (Chrm3, Nrxn3, Galntl6, Cadm2, Cntnap2), module 20 (Kcnip4, Ndst4, Gabrg3, Col11a1, Cntn4), and module 21 (Dcc, Thsd7b, Grm1, Ntm, Calb1). Module 26 includes activity-associated candidates such as Nr4a2, Homer1, Nr4a3 and Egr1. Module 32 contains 13 mt-prefixed genes; module 42 includes Slc1a2, S100b, Aldoc and Sparcl1. Such modules require biological/quality review before a neuron-specific mechanism is claimed. The review does not diagnose contamination, alter labels or exclude cells from these patterns. “Top module genes” here means highest native graph Moran I, not network hub centrality. Native Monocle3 module IDs are arbitrary and do not correspond to prior hdWGCNA module IDs.

## Observed inhibitory genes and reference sensitivity

C1 graph-associated leads include Cntn5, Pcdh15, Il1rapl2, Chrm2, Necab1, Gpc6 and Kirrel3. C6 leads include Galntl6, Adamtsl1, Rgs6, Hs3st5, Ncam2, Mgat4c and Pcdh9. C8 leads include Grm1, Ptprt, Grin3a, Mgat4c, Galntl6, Cntnap2 and Pth2r. These observations characterize within-scope molecular geometry and require independent validation. C8's mixed/uncertain atlas identity is preserved; learning its graph does not resolve that label.

Changing the arbitrary reference from a component medoid to a farthest tip materially changes cell ordering. Spearman correlations between the two saved reference-distance rankings are −0.2951 and −0.2235 for the two CA1 components, −0.4111 for C1, −0.4064 for C6, and −0.1437 for C8. These are descriptive comparisons of actual saved native distances, not re-fits or biological root validation. They do not establish a clean reversal or biological transition. Distances from distinct components/models do not share a universal time axis.

## Required interpretation and publication wording

- Label primary outputs **exploratory unrooted molecular state graph**. Label optional ordering **arbitrary-reference graph distance**, with the reference and alternative visible. Biological temporal direction and latent developmental time remain unavailable in this P60 single-snapshot analysis.
- Native graph_test measures same-data graph autocorrelation, with a separate BH family across successful tested features in each scope. It is not genotype DE, an independent-mouse hypothesis test, or proof of trajectory dynamics. The 9,056 CA1 discoveries (85.3% of tested features) should be interpreted in that exploratory context.
- Native zero p/q values (CA1 1,124 rows; C1 80; C6 6; C8 7) reflect finite numerical precision/underflow. Preserve the native tables but do not claim exact zero biological probability.
- Native gene-overlay colors are per-gene percentage of maximum expression; they are not detection fraction and do not provide a common absolute expression scale across genes.
- The module heatmap contains mouse/endpoint means of native module aggregation with display z scores clipped to [−3, 3], not genotype log fold changes. Preserve the raw unscaled aggregation table beside the display table. Mice remain the independent biological units.
- Same-partition sensitivity outputs are a descriptive filter of the same fitted graph, not an independent topology re-fit. All-mouse support is coverage, not demonstrated leave-one-mouse-out or bootstrap stability.
- No new original DESeq2 fit, strict DEG-only communication record, cell annotation, hdWGCNA result or pySCENIC method was changed by this run or audit.

## Durable receipts and primary method sources

The machine-readable final receipt is `INDEPENDENT_FINAL_REVIEW_RECEIPT_20261028.json`. Detailed proofs are `ACTUAL_CA1_COMPLETED_AUDIT_20261028.json`, `ACTUAL_IN_COMPLETED_AUDIT_20261028.json`, both `*_PDF_PARSE_20261028.json`, and the actual observed-results JSON/TSV files beside this review. The earlier partial CA1 checkpoint is historical and superseded by these completed audits. Synthetic checker-calibration receipts remain separate from study results.

Official documentation: [Monocle3 trajectories](https://cole-trapnell-lab.github.io/monocle3/docs/trajectories/), [Monocle3 differential-expression and graph tests](https://cole-trapnell-lab.github.io/monocle3/docs/differential/), and [exact native source commit](https://github.com/cole-trapnell-lab/monocle3/tree/536f1033d6de7c957f26a1f403f81efbd825e0db). No independent temporal root is supplied by the observed study design.
