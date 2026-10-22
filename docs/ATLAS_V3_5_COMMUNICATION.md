# Interactive CA1–inhibitory ligand–receptor explorer

The [communication explorer](https://ali-mujahid.github.io/satb2_snRNA/communication/) is linked from the original nuclei atlas and its downloads page. It adds display controls around the completed Seurat original_full analysis. No biological model was rerun and no previous analysis was overwritten.

Choose individual inhibitory clusters or inhibitory total, pooled/WT/cKO data, either directed set or their union, and an evidence tier. Source/target filters, search, matching gene-level DEG support and rank cutoffs apply together. WT/cKO comparisons use common population positions, dotplot rows/columns and numeric/size scales. Counts and score differences remain descriptive.

## Figures and controls

- **LR dotplot:** magnitude rank, geometric-mean expression, group-versus-rest LR logFC or native CellPhoneDB specificity P; adjustable palette, low/mid/high colors, color limits, opacity, dot radius, text size, order, up to 120 LR identities and paginated population pairs. Dot area always means minimum detection over all complex subunits.
- **Heatmap:** distinct ordered LR counts or mean geometric-mean expression. Untested directions are gray NA. Eligible pairs with no qualifying records have count 0; this does not establish absent communication.
- **Directed circle:** population nodes with count-scaled areas; counts determine arrow width. Arrowhead size, width, curvature, opacity, node size and population/edge colors are adjustable.
- **Circos:** real directed ribbons, adjustable ribbon multiplier, curvature, arrowheads, arc gap/thickness and colors. Widths are within-sector fractions; absolute counts remain in the heatmap/circle/tooltips.
- **Gene/complex circle:** population × complex × role nodes; fixed selected interaction keys across conditions, geometric-mean edge widths and detection-fraction node areas. Missing condition records retain unavailable evidence rather than fabricated zeros.

Every plotted view exports its current appearance as PDF, SVG or PNG. PDF uses a high-resolution raster of the rendered SVG; SVG retains vector geometry. Plot-state JSON includes current settings, filters, exact native rows and provenance. Manual color limits only clip displayed colors, never the saved scores.

## Evidence and complete tables

All 64,936 native records and matching LR-subunit DESeq2 evidence are embedded in the interactive transport. The results table has all filtered records independent of plot limits. Inspect a row for the original matching sender-ligand and receiver-receptor DESeq2 effects, stable IDs, adjusted P values, ambiguity/NOT_TESTED flags and eligible-mouse support.

Seven-column LR downloads retain `source,target,ligand_complex,receptor_complex,lr_logfc,lr_gmeans,magnitude_rank`. Extended downloads add the native evidence fields. Original seven-column gene DESeq2 files are linked separately, without modifying their bytes. Large native downloads live outside the Pages artifact under the immutable `atlas-v3.5.0` repository tag.

The standalone `SATB2_LR_interactive_offline.html` release asset embeds all render dependencies, six fits and matching LR-component evidence. Filters, figures, table/component downloads and appearance controls operate without a server. Links to repository tables need internet; the complete 156 MB original-analysis ZIP contains those tables plus 276 precomputed figure panels for local use.

## Scientific interpretation

The default pooled DEG-supported tier has 1,589 population-specific records. The stricter replicated/context-specific tier has 324 (112 CA1→IN, 212 IN→CA1), with 322 CA1-only and 2 both-end matching gene-level DEG-supported records. C6 Cntn3→CA1 Ptprg is the both-end lead; opposing component RNA changes do not specify a net change in communication.

`lr_logfc` means native cell-group-versus-rest specificity, not cKO/WT fold change. LIANA ranks are priorities, and native CellPhoneDB P values concern cell-label specificity. No interaction-level genotype significance test was fitted. Deep/superficial are expression programs, not measured physical layers. RNA and curated priors do not establish protein localization, contact, activation/inhibition or synaptic function. Full interpretation and limitations are available from the explorer.

## Verification

`scripts/verify_communication_atlas.py` verifies the compressed payload, all six dataset sizes, unique native keys, 164 transport hashes, native-table preservation, links and Pages size. Four Node suites check numerical filters/native evidence, charts, network geometry and actual app code with VM/DOM mocks. These checks do not open native expression objects or fit models. Browser downloads and rendering are verified separately against the deployed HTTPS atlas.

`scripts/build_communication_portable.py` builds the standalone HTML and a source-SHA receipt. `scripts/package_communication_atlas.py` reproducibly transports the completed report; it is not a scoring workflow.

Version3.5.1 preserves the exact visible plot’s settings, native rows, title and filename during asynchronous PDF/PNG generation, including changes made before the pending redraw. Dotplot captions retain the full filtered LR and population-pair totals. The six scientific data fits and all native data URLs remain pinned to3.5.0.
