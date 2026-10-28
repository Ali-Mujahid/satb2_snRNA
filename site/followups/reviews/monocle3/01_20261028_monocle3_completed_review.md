Takeaway: SATB2 Monocle3 recovered four expression-state graphs and gene modules; biological temporal direction remains unresolved.

This analysis examined molecular state structure within the selected CA1 programs and inhibitory clusters in the WT and Satb2 cKO study. Native Monocle3 completed one pooled CA1 analysis and separate analyses of inhibitory C1, C6 and C8. All four models retained twelve mice and produced graph-associated genes, gene modules, native objects and visual reports. The CA1 analysis recovered two disconnected components, including a small mostly superficial-associated component represented in every mouse. Ordering cells from different arbitrary reference nodes materially changed their rankings, so these results support expression-state maps rather than a biological clock or a demonstrated deep-to-superficial transition. The immediate use is to examine the exported gene/module patterns alongside the original independent-mouse DESeq2 evidence while preserving label and technical uncertainty.

## Scientific context and question

The existing P60 single-nucleus study has six independent WT and six Satb2 cKO mice, sampled at one age. The requested Monocle3 endpoints were trajectories/pseudotime, graph-associated genes and gene modules. No independent temporal root or validated developmental progression is supplied by this design. CA1 deep/superficial-associated labels are molecular expression programs; C8's mixed atlas identity remains unresolved. Genotype partly overlaps sequencing partition.

## Question, lifecycle, and endpoint

The registered native workflow completed all four steps with return code zero: registry `081614e2-5ca8-48a8-9b33-b431673f0eee`, native run `snakemake-satb2-monocle3-states-20261028-20261028-112530-8575a0c7`. Root: `/Users/mujahid357/Documents/codex_new/SATB2_FOLLOWUPS_20261028/runs/monocle3_20261028_v1`; scientific results are under `results/monocle3`.

Pinned native Monocle3 1.4.27 used the exact frozen prepared integer RNA counts. Features were selected without genotype or DEG/HVG selection using detection and variance criteria. Native normalization uses retained-feature size factors; original, prepared-all-feature and retained-feature UMI totals are exported separately. Native PCA/UMAP, partitions and principal graphs were fitted without changing original atlas coordinates. CA1 DL/SL were pooled, while C1/C6/C8 used separate CDSs; no cross-lineage or cross-interneuron bridge was forced. Native principal-graph tests and modules characterize same-data molecular structure; no genotype gene-DE model was fitted.

## Key findings

| Model | Nuclei | Features tested | Components / partitions | Graph-associated genes at native BH≤0.05 | Modules |
|---|---:|---:|---:|---:|---:|
| CA1 DL/SL programs | 20,456 | 10,618 | 2 / 2 | 9,056 | 47 |
| Inhibitory C1 | 1,718 | 10,539 | 1 / 1 | 3,669 | 19 |
| Inhibitory C6 | 501 | 9,849 | 1 / 1 | 1,439 | 8 |
| Inhibitory C8 | 468 | 10,258 | 1 / 1 | 2,002 | 11 |

Every model contains all twelve mice. All tested native rows succeeded, and every selected graph-associated gene had PCA-loading coverage and a native module assignment. Native counts, cell/feature axes and size-factor arithmetic matched the frozen source matrices independently. The two lineage seals and 174 actual output receipts passed review; all 32 one-page PDFs parsed, and representative downloadable PDFs and PNGs were inspected.

CA1 has 194 principal nodes and 192 edges. Its main component contains 19,609 nuclei: 8,748 DL and 10,861 SL, comprising 10,186 WT and 9,423 cKO nuclei. The smaller component has 847 nuclei: 845 SL and two DL, comprising 577 WT and 270 cKO nuclei. Both components include all twelve mice. Thus the small component is not attributable to a single missing or isolated mouse, but its biological identity or temporal role is not established.

CA1 graph-associated examples include Chrm3, Zbtb20, Nrxn3, Mctp1, Spock1, Galntl6, Epha6 and Cadm2. Module 9 contains Chrm3/Nrxn3/Galntl6/Cadm2/Cntnap2; module 20 contains Kcnip4/Ndst4/Gabrg3/Col11a1/Cntn4. Module 32 contains thirteen mt-prefixed genes, and module 42 includes Slc1a2/S100b/Aldoc/Sparcl1. These latter patterns warrant quality/background interpretation before proposing a neuron-specific mechanism; they did not trigger new exclusions or diagnoses. Native module IDs are arbitrary and do not correspond to hdWGCNA module numbers.

All nuclei have finite native distance from a component-specific arbitrary reference; no infinite distance was replaced with zero. Alternative-root rank correlations were −0.295 and −0.223 in the two CA1 components, and −0.411, −0.406 and −0.144 in C1, C6 and C8. This is substantial ordering sensitivity, not a verified temporal reversal. Primary saved CDSs remain unrooted, and biological pseudotime remains explicitly unavailable.

## Interpretation

The analyses recover within-scope molecular organization and gene sets that vary across it. They provide a useful way to examine CA1-program and inhibitory heterogeneity, including candidate modules containing supplied CA1 marker genes. They do not demonstrate developmental transitions, elapsed time, WT→cKO progression, C1→C6→C8 conversion or a newly resolved C8 identity. Strong genotype/program structure is visible in the CA1 geometry; partly confounded sequencing partition and other technical or molecular differences remain alternative contributors.

Native graph association is a cell-level same-data autocorrelation screen, with a separate native BH family within each model. It is not an independent-mouse genotype test, nor independent confirmation of the original DEGs. The large CA1 screen result—9,056 of 10,618 tested genes—must be interpreted accordingly. Numerical native p/q zeros reflect finite precision, not exact zero biological probability. Preserve original DESeq2 results for genotype gene inference.

The module heatmaps show mouse/endpoint aggregate expression with display z scores clipped to [−3,3], rather than genotype fold changes. Raw aggregates are provided beside them. Native gene overlays use each gene's percentage of its own maximum expression; this is neither detection fraction nor a shared absolute scale across genes. Reference distances from distinct components or models do not share a universal time axis.

## Artifacts

**Results:** Each model directory under `results/monocle3/{CA1,IN}/models/` contains `graph_test_all.tsv`, `gene_modules.tsv`, `module_candidate_reconciliation.tsv`, `module_mouse_endpoint_raw.tsv`, `module_mouse_endpoint_display_z.tsv`, `principal_nodes.tsv`, `principal_edges.tsv`, `embedding_cells.tsv`, `reference_nodes.tsv` and `reference_distances.tsv`. Native unrooted objects and their required nearest-neighbor indexes, plus separately saved reference-ordering objects, are retained.

**QC and visual reports:** `SCOPE_VALIDATION.json`, `gene_eligibility.tsv`, `library_denominators.tsv`, component/partition/reference-bin mouse support and `SAME_PARTITION_DISPLAY.json` retain coverage and interpretation. Each model has eight original PNG/PDF figure pairs covering endpoint/genotype/mouse/technical partition, native partitions, arbitrary-reference distance, gene overlays and module heatmaps. Representative original figures and downloadable PDF renders were readable. Upstream ggplot aesthetic deprecations were cosmetic and did not stop output generation.

**Provenance:** Native final/lineage seals, input/config/source/runtime receipts, session information and installed-file closure are preserved. The final Workbench record is `/Users/mujahid357/Documents/codex_new/SATB2_FOLLOWUPS_20261028/receipts/monocle3_terminal_record.json`. Independent audit code and actual-model receipts are under `/Users/mujahid357/Documents/codex_new/SATB2_FOLLOWUPS_20261028/reviews/monocle3_readonly_audit/`.

**Model synthesis:** This review and `INDEPENDENT_FINAL_REVIEW_20261028.md` in the independent audit directory distinguish observed molecular structure from unsupported temporal/genotype claims and identify quality-associated modules.

## Limitations and blockers

Biological temporal direction is unestablished. Same-partition exports are a descriptive filter of the same fitted graph, not a new graph fit. Twelve-mouse coverage is not demonstrated leave-one-mouse-out, seed or bootstrap topology stability. Genotype/partition confounding, expression-program identity, ambient/background patterns and C8 label uncertainty remain. Different reference choices materially alter the ordering, and distances cannot establish anatomical transitions. The analysis does not replace RNA velocity or supply the missing spliced/unspliced counts.

## Recommended next step

Publish the reviewed outputs as exploratory unrooted state graphs, with arbitrary-reference and alternative-reference distances clearly labeled. Use raw mouse/module tables and original DESeq2 fits to assess candidate gene patterns, keeping mitochondrial/glial-associated modules visible as review signals. A biological chronology requires independent temporal or kinetic evidence; do not force a root from genotype, anatomical assumptions or an attractive UMAP.
