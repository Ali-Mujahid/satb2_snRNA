# Native expression summaries and reviewed plots

These modules read completed native aggregate tables and existing metadata/UMAP
coordinates. They never load or aggregate raw nucleus count matrices and never
fit a normalization, embedding, clustering or genotype model.

`summarize_aggregates.py --config MASTER --outdir RESULTS/marker_summaries`
reads each `prepared/BRANCH/gene_cluster_mouse_aggregates.tsv.gz` and
`cluster_mouse_counts.tsv`. The native prepare receipt must have stage `prepare`,
status `EXECUTED_REVIEW_REQUIRED`, no fitting, declared native count aggregation,
the matching branch and exact hashes for both tables. Alternative inputs can be
provided with repeated `--aggregates` and `--cluster-counts` arguments.

Aggregate keys are `tool,branch,population,cluster_scope,cluster_id,sample` plus
`gene_id,gene_symbol`; values are `condition,n_nuclei,sum_counts,n_detected,
sum_log1p_normalized`. Every selected gene has explicit rows for every atomic
cluster and all twelve mice, including zero-nucleus rows. Log sums are produced
natively from `log1p(count/full_source_nucleus_library*10000)` before summation.
Means and detection fractions use every native cohort nucleus; marker-review
mouse-support thresholds never remove nuclei from these means or from DE.

Outputs include pooled `gene_cluster_summary.tsv`, explicit-denominator
`gene_cluster_summary_WT.tsv`, all-genotype per-mouse and condition summaries,
`program_evidence_WT.tsv`, representative marker selections and their exhaustive
criteria audit. WT absence remains NA; `n_WT_nuclei` and `n_WT_mice` are explicit.
Annotation uses original-full WT evidence to define expression programs and
carries those calls into the paired branches. Programs never establish physical
layers. This module's program table is descriptive arithmetic, not per-nucleus
scoring or a new model.

Representative genes require native detection/specificity, mean-expression
contrast and cross-mouse support. Original-full selected genes are fixed for
paired comparisons. Clusters without a supported candidate remain present in
the registry and plots; no marker is invented or selected by a tiny P value.
The default marker-review rule is at least five nuclei per mouse and ten-percent
detection, with support in at least three mice. This confidence rule is separate
from the explicitly configured pseudobulk rules.

Gene scaling uses equally weighted cluster mean-log1p values and population SD
(ddof 0), after the explicit atomic-cluster mask. Original-full plots use their
own scale; paired corrected plots use the paired-original scale on the same
frozen population roster. Exact means, z values and detection fractions are
saved. This is an explicit common full-cohort dotplot definition, rather than a
claim of numerical identity with Seurat's default sampled DotPlot calculation.

`render_plots.py --config MASTER --summaries RESULTS/marker_summaries
--annotations RESULTS/annotation_final --outdir RESULTS/plots` requires complete,
hash-bound final annotation and summaries. Every main book uses the reviewed
`original_frozen` inference mask. Shared All/EX iterative filters are checked
between tools; raw/corrected masks are checked pairwise. Secondary `fresh_branch`
books contain coordinate/annotation views only when no matching native fresh
aggregate exists. Frozen expression is never relabeled as a fresh-cluster result.

Books live under `plots/TOOL/BRANCH/POPULATION/CLUSTER_SCOPE/`, with PDF, individual
SVG/PNG pages, exact dot TSVs, panel memberships, rendered metadata, count tables,
page indices and per-book COMPLETE receipts. `PLOT_BOOK_INDEX.tsv` and a global
COMPLETE index the 24 main books and 16 secondary paired books. All supplied panel
memberships are represented; the two unmapped supplied rows remain an explicit
unmapped-row table. Actual page counts depend on representative-gene unions.

Negative scaled expression changes glyph area by the configured multiplier
(default 0.25) and alpha (default 0.18). Detection percentages and exact z values
remain unchanged. Every dotplot explains this additional visual downweight;
glyph area/alpha are also saved in its numerical table. Zero detection remains
zero glyph area. UMAPs retain all displayed nuclei, readable cluster identifiers,
reviewed molecular and original-query atlas labels, WT/cKO colors and twelve-mouse
facets. Existing coordinate fits are not rerun after a display/analysis mask.

`test_synthetic.py --outdir NEW_FIXTURE_DIRECTORY` tests weighted denominators,
zero-WT handling, paired scale anchoring, glyph semantics and end-to-end rendering
using synthetic sums/coordinates only. It never reads study count matrices.
