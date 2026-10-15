Takeaway: Twelve-mouse SATB2 snRNA-seq supports total-population RNA effects and provisional fine CA1 program labels.

Twelve independent WT and Satb2 cKO mice were compared to assess consistency between Scanpy and Seurat v5 single-nucleus results. The native revision completed its annotation, count validation, figures and mouse-level DESeq2 comparisons while preserving the earlier analyses. Shared population totals have identical integer counts, and the prespecified sensitivity checks support 604 all-nucleus, 785 excitatory, three inhibitory and 1,805 CA1–ProS population-gene effects. Finer CA1 differences involve divergent Harmony neighborhoods, graph partitions and threshold-sensitive marker programs, while superficial program equivalence remains unsupported. These results support relative RNA differences under the recorded models, with anatomical identity, population composition and partially confounded sequencing partitions limiting finer interpretation. The justified next step is to explore the qualified evidence in the separate offline viewer and reserve edited-population inference for an exact saved fit or a subsequent native analysis.

## Scientific context and question

The study compares six WT and six Satb2 cKO independent mice using 10X single-nucleus RNA-seq. Existing study metadata identifies P60 male mice, one GEM-X 3′ v4 capture per mouse, and CA1 tissue including ventral sampling. The previously validated genotype/deletion context is retained. The requested comparison covers all nuclei, excitatory nuclei, inhibitory nuclei and CA1–ProS candidates, including original clusters and reviewed marker-associated cluster unions.

Existing CellRanger 10/mm10 matrices, verified CellBender outputs and the exact prior MapMyCells query/results were reused. The MapMyCells reference is Allen WMB 20240831/CCN20230722, with mapper 1.7.4. The exact CellRanger GTF and intron-counting setting remain unresolved; this revision did not rerun alignment. A shared reference, source features and statistical engine provide controlled comparisons rather than independent methodological validation of every stage.

## Question, lifecycle, and endpoint

The Workbench registry run `30bef3d0-c984-4432-92c4-4059580c30b5` is **completed**, with 19/19 tasks completed and no reported failure. The run directory is [20261014_v2_0_0](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261014_revision/runs/20261014_v2_0_0). Its immutable scientific workflow/configuration and the earlier 13 October run remain preserved.

The native revision performs the authorized fixed-cohort/cluster annotation and count/model comparison. Original counts are primary. Paired original versus paired CellBender counts hold nuclei, labels and the paired-original gene universe fixed. Main memberships use the original Harmony 0.4 clusters; saved embeddings are not recomputed after the requested size mask. Fresh paired views are explicitly secondary. The subsequently requested webpage uses produced metadata, aggregate statistics and fitted results only; its engineering performs no raw-nucleus RNA aggregation or scientific fitting.

## Key findings

The shared iterative minimum-100 cluster policy for All/EX removes any nucleus belonging to a small cluster in either tool and converges after one iteration. Every original cluster remains audited, including empty excluded groups. IN and CA1–ProS have no such cutoff.

| Population | Original retained nuclei | Paired retained nuclei | Original PRIMARY BH<0.05 | Core sensitivity-supported genes |
|---|---:|---:|---:|---:|
| All nuclei | 107,747 | 103,983 | 645 | 604 |
| Excitatory | 60,908 | 60,772 | 875 | 785 |
| Inhibitory | 8,911 | 8,903 | 4 | 3 |
| CA1–ProS candidates | 27,703 | 27,672 | 2,025 | 1,805 |

These are separate, overlapping population-gene results; the rows do not sum to a unique gene-discovery count. Core support requires BH<0.05 and the same nonzero sign in all three PRIMARY count branches, both original partition estimates tested/finite/same-sign, and all twelve leave-one-mouse-out estimates tested/finite/same-sign. Sensitivity fits need not each be significant. Shared totals use one R DESeq2 implementation, so their equality demonstrates count/implementation reproducibility.

The size mask removes 315 All and 77 EX nuclei per branch. Original All removals comprise 174 WT/141 cKO; EX removals comprise 74 WT/three cKO. These are consequences of the requested shared size rule and affect the defined cohort, rather than proving that excluded cells are invalid or differentially abundant.

Independent Python/R validation checked **556,722,540 integer entries in 1,437 atomic/final population arrays**, with zero differing entries. All twelve shared-total tool pairs have identical nuclei and the full 32,285-gene axis. The final comparison records 642 effect pairs, 420 fine membership assignments and 163 explicit skips.

Three DESeq2 producers publish 1,060 endpoints: 964 tested and 96 unsupported, comprising 360 fresh fits and 604 exact unchanged/repeated endpoint reuses. The 948 PRIMARY endpoints comprise 852 tested and 96 unsupported; the remaining 112 are core sensitivity endpoints. Recorded warnings/errors are zero. PRIMARY dispersion fitting is parametric in 844 cases and local in eight; actual settings are retained.

IN pseudobulk sums every available nucleus, with all nonempty positive-library mouse aggregates eligible. It requires at least three independent mice per genotype, residual degrees of freedom of at least two, and gene counts of at least five in at least three eligible mice. Seventeen tested IN PRIMARY endpoints per branch include an animal with fewer than five nuclei and remain exploratory. Other populations require at least twenty nuclei per included mouse. Zero/absent contributions remain visible in all-twelve-mouse audits.

The three core-supported IN effects are **Apool↓, Arhgap31↓ and Pdia5↑**, with original log2 fold changes −1.318, −0.815 and +1.064. Illustrative supported All/EX effects include Adamts3↑ and Sema3a↑; CA1–ProS includes Col5a2↓ and Sema3a↑. These describe relative RNA abundance in the defined population, not independent evidence of cell loss or a pure within-type mechanism.

Satb2 gene-body RNA is relatively higher in cKO in All (+0.971), EX (+0.956) and CA1–ProS (+1.166), each meeting the core flag; IN Satb2 is not significant. The cKO/WT contrast and stable-ID binding were verified. The 3′ gene-body assay does not measure deleted exon 4 or protein, and the result does not reopen the validated genotype or establish a feedback mechanism.

CA1 original rotation-aligned PCA kNN20 overlap is 0.9965; Harmony overlap is 0.8481. Source graphs have mean off-diagonal degrees 28.90 versus 66.79 and edge Jaccard 0.3150. Increasing controlled Scanpy Leiden iterations changes cross-tool ARI only from 0.6818 to 0.6923. Tested Python Harmony caps 10/50 converge to identical representations, with actual completion after 5/5/6 iterations; the cap does not explain the tested discrepancy. These controls do not isolate every graph/implementation contribution or establish a biological optimum.

WT-defined same-name CA1 program agreement is 64.5% at the primary threshold and changes to 83.6%/85.5% at the two reviewed threshold alternatives. Same-name deep-associated overlap has reciprocal coverage 95.6%/68.8%, while superficial-associated coverage 89.7%/44.9% fails the declared reciprocal guard. Different-label matches are population-definition sensitivity. All 163 original clusters have review records; supported atlas subclass evidence, rare/mixed/conflicting groups and physical-layer uncertainty remain separate.

## Interpretation

The strongest current evidence concerns the shared total populations and the prespecified count/partition/mouse-omission sensitivity checks. Tool-specific fine groups differ in membership, and apparent fine CA1 effects cannot be interpreted through a matching name alone. Supported CA1–ProS subclass records mostly corroborate the shared atlas subclass rather than distinct validated layers or states. The marker sheets supply expression-associated programs, including perturbation-context candidates; they do not establish physical superficial/deep position or inhibitory morphology.

The main sequencing partition contains six WT and three cKO mice; partition-sensitive agreement has substantial shared support and is not independent replication. Remaining batch/sample-composition effects are not ruled out. BH correction is within each endpoint, not globally across overlapping category/cluster hypotheses. CellBender is a count sensitivity, with incomplete discovery-list identity and no direct ground truth for the origin of removed RNA.

## Artifacts

**Results:** [native comparison](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261014_revision/runs/20261014_v2_0_0/results/comparison/comparison_review.md), the three native DESeq2 result branches, exact seven-column symbol-based result files, stable-ID audits, all-mouse coverage and core/fine comparison tables.

**QC and visual review:** 40 native annotated UMAP/marker books containing 550 pages, exact dot values and source memberships; [seven-page effect overview](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261014_revision/reports/effect_figures/DESeq2_effect_overview.pdf); CA1 technical/program sensitivity PDFs; [annotation review](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261014_revision/reviews/annotation_result_review.md).

**Provenance:** immutable workflow/configuration, stage completion receipts, producer/model audits, exact-count acceptance, output hashes, package receipts and stable row/ID bindings. The offline data/catalog cover 24 main views, 878 precomputed marker IDs and 1,060 model endpoints. The 147-row reviewed supplied/canonical panel contains 145 mapped rows/139 unique IDs; its workbook-only subset contains 123 mapped rows/117 unique IDs. Both unmapped AC129186.1 rows remain explicit.

**Model-authored synthesis and viewer:** [DESeq2 interpretation](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261014_revision/reviews/final_DESeq2_effect_review.md) and [portable offline explorer](/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261014_revision/reports/offline_dashboard_v1/index.html). It supports side-by-side saved UMAPs, reference/reviewed/program annotations, cluster exclusions, marker/color choices, shared continuous scales, local labels/presets, exact-fit DESeq2 browsing, qualified reliability and PDF/source downloads.

## Limitations and blockers

There is no outstanding native execution blocker. Fine identity, anatomical layers, inhibitory morphology, cell abundance and mechanisms remain qualified as above. Unsupported model endpoints retain their reasons. Missing FDR/effects remain unavailable. Corrected counts, shared atlas assignments and sensitivity flags are not independent biological replication or causal proof.

Offline exclusions/label unions cannot run DESeq2. An edited total shows statistics only for an exact existing fit on that source branch; otherwise it exports a request and shows no matching fit. Marker browsing covers the existing 878-gene panel. Gene UMAP overlays are cluster means, not individual-nucleus expression. Real browser/CSS and screen-reader inspection remain unverified after the automation browser explicitly rejected local file navigation; source, actual-data state, simulated canvas/PDF and local integrity checks were used without bypassing that policy.

## Recommended next step

Use the portable viewer to inspect source-supported labels, mouse coverage, exact fits and their sensitivity limits. Preserve an exported settings/selection JSON for any proposed annotation union that lacks an exact saved fit. Such a population can then be reviewed and fitted through the native workflow if subsequently requested; no additional analysis is needed merely to use the page.
