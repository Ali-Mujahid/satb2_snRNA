Takeaway: CellChat's figure-only correction passed; 76 exports preserve native data and fix labels and boxplot heights.

This review checked whether the completed CellChat figures can be presented without changing the underlying communication results. The separate layout-only workflow completed its corrected PNG and PDF exports from the existing saved tables. All original scientific file receipts remain unchanged, and every corrected figure maps to its preserved native original. Both averaging methods now plot all individual-mouse boxplot points at their exact saved score, including genuine zeros. PDF-bound checks and visual inspection found readable headers, rotated labels, axes and legends in the corrected exports. The corrected figures are acceptable for display, while the original exploratory interpretation and its common-partition limitations remain unchanged.

## Scientific context and question

The underlying study contains six WT and six Satb2 cKO P60 mice. The completed CellChat analysis used the selected CA1 deep-associated and superficial-associated molecular programs and inhibitory clusters C1, C6 and C8. Both endpoints' required ligand/receptor subunits had to be original matching-population DESeq2 DEGs. The layout correction addresses display defects only: clipped circos headers, tight/clipped heatmap PNG context labels and default vertical boxplot jitter. It does not change the biological question, eligibility, labels, scores, comparison groups, statistical tests or limits of RNA-based communication inference.

## Question, lifecycle, and endpoint

The actual display-only workflow is complete: registry `b78aeedd-a3b4-429a-8d24-ce40d4e787a4`, native run `snakemake-satb2-cellchat-layout-only-20261028-20261028-115316-3dd73358`, scientific output `/Users/mujahid357/Documents/codex_new/SATB2_FOLLOWUPS_20261028/runs/20261028_cellchat_layout_v1/results/layout`.

Its native seal intentionally remains `COMPLETE_PENDING_VISUAL_REVIEW`. This external, read-only review supplies the actual display acceptance without modifying that seal. The original scientific run registry is `7e8d72cb-c013-43b2-a7fe-89e271430333`. No CellChat inference, expression normalization, interaction rescoring, permutation test or gene-DE fit was run for this audit. The in-memory `ggplot_build` check captured point coordinates without opening a figure device or saving a new scientific plot.

## Key findings

| Check | Actual result |
|---|---:|
| Original scientific output receipts reverified unchanged | 115 |
| Corrected layout output receipts verified | 79 |
| Original-to-corrected figure exports verified | 76 |
| Corrected one-page PDFs parsed and bounds checked | 38 |
| Corrected PNG files structurally verified and contact-sheet reviewed | 38 |
| PDF text outside page bounds | 0 |
| Raw mouse scores exactly matching plotted boxplot heights | 216 |
| Nonzero vertical point displacements / negative plotted points | 0 / 0 |
| Genuine zero triMean points retained at zero | 21 |
| Frozen installed-runtime file hashes independently verified | 12,643 |

The methods each contain 108 plotted mouse points. The prior saved-coordinate proof was independently compared with the native pathway tables by exact mouse/source/target/pathway key; all keys and values matched, with zero numeric error. A separate audit copy executed the same frozen display helper with `cc_save_pair` replaced by a no-render capture. It independently required exact equality between each built point-layer y coordinate and its raw saved score. The native summary, network-matrix, endpoint-palette and save-helper bodies were also byte-equivalent at the function-body level to their original counterparts.

All 76 corrected exports preserve their original relative figure identities. Each PDF has one page and no text character outside the page, with the required context headers/captions present. Both methods' long-label heatmaps, circos headers and mouse boxplots were visually reviewed through six actual PDF renders. Full-resolution corrected dotplot and circle PNGs were also inspected, and the two complete contact sheets cover all 38 PNGs. Headers, legends and labels are readable, and the original zero-score plotting artifact is resolved.

## Interpretation

The correction is a faithful presentation derivative. It adds margins/heatmap canvas space, puts circle/circos context text inside the page and sets vertical jitter to zero while retaining the same horizontal jitter. Original scores, summary formulas, colors, numeric scales, pathway sums and statistical tables are preserved. The fixed pictures should be the default reader downloads; the original pictures should remain available as explicitly archival layouts with known display issues.

This display acceptance does not strengthen the underlying biological evidence. The original primary triMean full-cohort analysis had three adjusted associations, but their primary shared-partition BH values were approximately 0.07143. The original sensitivity method gives different significance coverage. These remain exploratory, same-study DEG-selected communication leads, not independent validation or proof of physical signaling. No additional statistical finding was created by this correction.

## Artifacts

**Results:** The exact native tables remain under `runs/20261028_cellchat_v1/results/cellchat`. Corrected plot files remain under `runs/20261028_cellchat_layout_v1/results/layout/figures/{triMean,truncatedMean10}`. `NATIVE_TO_PRESENTATION_PLOT_MAPPING.json` binds every corrected file to its original path, byte length and SHA256.

**QC and visual reports:** `INDEPENDENT_INTEGRITY_AUDIT.json`, `20261028_actual_point_y_audit.json` and `20261028_actual_point_y_proof.tsv` are in this external review directory. Actual PDF-render review images and complete PNG contact sheets are preserved in the earlier `reviews/CellChat_layout_actual_readonly` directory. All native scientific and corrected-output hashes were checked again during this independent review.

**Provenance:** `DISPLAY_ACCEPTANCE.json` records the original and layout terminal-seal identities, mapping, runtime closure, independent numerical proofs and visual evidence receipts. The original scientific seal SHA256 is `c1dbdbe57c1606ab37d71cbf630905731bb630667038a25bd5ec0cd6bfd16683`; the layout seal SHA256 is `e15002fe6b056c9efe7d02dbe3d137cd3f041fb1f4e26a1e533ab2a655abda92`; the mapping SHA256 is `4ff34c93f451ca0baaca3259a3099760a713e3e24033cf3be20d4fb343ea6e5a`. Frozen source, configuration and runtime paths remain bound by the native workflow receipts.

**Model synthesis:** This review explains why the correction is acceptable for display while preserving the original biological and statistical caveats. It does not change Workbench lifecycle or replace the parent agent's run review/submission.

## Limitations and blockers

No remaining display blocker was identified in the reviewed exports. PDF text bounds alone cannot prove visual quality; representative full-size PDF/PNG checks and all-PNG contact-sheet inspection complement that test. The audit did not fit or independently regenerate any scientific result.

Circos sector/ribbon geometry is relative within a condition; absolute cross-condition score changes must not be read from ribbon width. Circle edge-width maxima are shared across views/conditions within an averaging method, rather than necessarily across methods. Dot/tile colors represent native per-figure magnitude and its legend, not gene up/down regulation or pathway activation. Empty/ineligible contexts and C8's absent eligible edge do not prove absent biological communication. Genotype/technical-partition confounding, same-data DEG selection, small common-partition cKO sample size and molecular-program identity limits remain.

## Recommended next step

Feature the corrected layouts in the portable/online reader using this display-acceptance receipt and the sealed original-to-corrected mapping. Preserve explicit original-layout downloads with their display-defect warning. Keep the unchanged scientific tables and original completed-run interpretation beside the plots. The parent agent can submit its terminal layout-run review through Workbench and then proceed with the already authorized additive publication gates.
