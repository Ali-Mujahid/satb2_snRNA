# Independent review of saved primary communication figures

Reviewed the saved PNGs under
`runs/20261022_liana_figures_final/results/report/figures` at full rendered size.
No code, scientific inputs, native results, or figure outputs were changed.
No models or tests were run during this review.

## Outcome

**Pass for the representative panels inspected. No material rendering defect
requiring a change was found.** The four requested chart types are present,
readable, and carry appropriate prediction/context caveats. This is a visual
and source-data spot check, not a claim that every exported panel or PDF was
individually reviewed.

## Panels inspected

- `IN_clusters/CA1_to_IN/replicated_specific/WT/circle.png`
- `IN_clusters/IN_to_CA1/replicated_specific/WT/dotplot_page01.png`
- `IN_clusters/IN_to_CA1/replicated_specific/cKO/dotplot_page01.png`
- `IN_clusters/IN_to_CA1/replicated_specific/pooled/circos.png`
- `IN_clusters/bidirectional/DEG_supported/WT/heatmap.png`
- `IN_clusters/bidirectional/DEG_supported/cKO/heatmap.png`
- `IN_clusters/bidirectional/DEG_supported/pooled/heatmap.png`
- `IN_clusters/bidirectional/DEG_supported/cKO/circle.png`
- `IN_total/CA1_to_IN/DEG_supported/WT/dotplot.png`
- `IN_total/IN_to_CA1/replicated_specific/cKO/heatmap.png`

## Checks and findings

- Titles identify direction, genotype/context, inhibitory resolution, and
  dotplot page number correctly. No inspected title, gene label, group label,
  legend or footnote was clipped. The dense panels remain readable at full
  size; the cluster circle is intentionally an overview with individual edge
  counts available in its source data.
- WT and cKO strict IN→CA1 dotplots retain the same LR row order and the same
  cell-pair columns on page 1. Both use magnitude rank 0–0.05 and matching
  10%, 50%, 100% dot-area legends. The legend explicitly says the minimum over
  all complex subunit detection fractions.
- WT, cKO and pooled bidirectional DEG-supported heatmaps share the same count
  color scale (maximum 33), group order and axes. The y-axis is ligand source;
  the x-axis is receptor target. Cross-category directions are populated,
  whereas CA1→CA1 and IN→IN tiles are gray/NA rather than false numeric zeros.
  Spot-checked WT source-data rows agree with the displayed first CA1 Deep
  counts (12, 12, 12, 10, 9) and with within-family NA values.
- IN review flags remain visible as asterisks for C8, C12, C16, C18, C19 and
  C20. Strict views retain these groups as isolated nodes or zero-count tested
  endpoints while excluding their uncertain-identity interaction rows, as the
  documented selection requires.
- Circle arrows and their blue CA1-source/orange inhibitory-source legends
  match the requested directions. The strict CA1→IN WT circle contains only
  blue arrows. The bidirectional cKO circle contains both colors; node labels
  remain outside the network.
- The strict IN→CA1 pooled circos has closed ribbons, orange inhibitory-source
  direction, and receiver tips at the CA1 sectors. Its caption correctly warns
  that equal-sector ribbon widths are within-sector fractions, not absolute
  cross-sector counts. Broad ribbons from a small source sector must therefore
  not be ranked visually as larger biological effects.
- The strict inhibitory-total heatmap correctly shows zero qualifying
  interactions for its tested directions. Its combined identity does not meet
  the exact `REFERENCE_SUPPORTED` filter; the report/README explains this
  exclusion. Those zeros do not imply absent biological communication.
- Strict captions explicitly list the rank, within-condition CellPhoneDB,
  mouse-support and reference-identity requirements. They state that the extra
  priority is exploratory context specificity, not genotype FDR. General
  captions distinguish prediction from validation and molecular CA1 programs
  from measured physical layers.

Supplemental gene/complex circle networks were not yet available at the time
of this primary review and are outside this review's coverage.
