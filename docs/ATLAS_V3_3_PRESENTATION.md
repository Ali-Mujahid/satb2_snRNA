# Atlas v3.3 presentation controls

Open the [SATB2 atlas](https://ali-mujahid.github.io/satb2_snRNA/). These controls change how existing results are displayed. Counts, cluster membership, annotations, original embeddings and saved DESeq2 fits remain unchanged. The exploratory inhibitory downsampling analysis remains stopped.

## Arrange the nuclei

The UMAP panel offers three layouts:

- **Original** uses the saved embedding. With spacing 1× and no manual transforms, this recovers the original coordinates.
- **Compact** brings cluster centers closer while retaining the point pattern within each cluster. The spacing control works across all categories.
- **Hippocampus C schematic** arranges clusters into a C-shaped illustration. It is an artificial layout: these nuclei have no measured tissue positions, and the picture does not establish hippocampal location or layer.

Rotation and flips change orientation. A tool-specific cluster editor also lets you translate one cluster along X/Y and adjust its visual scale. Manual changes and schematic layouts are marked as presentation layouts, including exported figures; distances in these views must not be interpreted as the original UMAP distances. Reset layout controls to recover the saved embedding. **Download display coordinates** records the transformed coordinates separately from the original coordinates.

To bring the distant Seurat inhibitory **C4** closer:

1. Select **Inhibitory** and open the UMAP panel.
2. Try **Compact** and reduce cluster spacing until the overall arrangement is useful.
3. For precise placement, choose **Seurat**, then **C4**, in the cluster editor. Move its X/Y position toward the other clusters and apply. Keep scale at 100% to retain the cluster's internal size.
4. Inspect both tools and download the figure. Use reset if the new arrangement hides points or suggests relationships you did not intend to show.

## Make selected marker dots smaller or lighter

In the marker appearance panel, choose **Both tools**, **Scanpy** or **Seurat**, then select a gene from the selected markers. Choose **All clusters** for a whole gene column, or an original cluster ID for one dot. Set the dot **area %**, **opacity %**, or **Hide**, then apply.

Area changes dot area: 25% area produces half the original diameter. Opacity changes transparency. These adjustments multiply the existing dot appearance, including the controls for negative z-scores. They do not modify expression, z-scores or detected fractions. Hidden or faded dots must not be treated as evidence that a gene is absent.

Whole-gene settings apply to every cluster unless an individual dot has an override. **Reset selected** removes that override and restores the inherited appearance; **Reset panel** clears styles for the chosen tool(s), branch and population. With merged-label grouping, only whole-gene styles apply because one displayed dot can summarize several original clusters.

Edited appearance is disclosed in figure/export information. Marker CSVs retain the saved scientific values and include presentation settings. **Export marker styles** saves the exact stable gene IDs, original cluster IDs, scope and multipliers.

## Download a clean volcano figure

The default volcano PDF keeps the nucleus category title, axes, ticks, plotted points and selected gene names. Gene labels have no leader lines. The full view and low-y view use the same source values, significance thresholds and colors as the atlas. Select detailed export when you want the additional legend, status and provenance text in the PDF.

Clean figure formatting does not remove provenance from the underlying results. Keep the exact source DEG file and exported settings with a figure so its tool, branch, target, thresholds and color settings can be recovered.

## Save and recover a presentation

Export settings before making a large set of edits. Settings export/import preserves the presentation controls alongside existing selections. Marker-style export and display-coordinate download provide focused records of those edits. To recover an unedited view, reset marker styles in each edited scope, select **Original**, restore spacing to 1×, and reset manual position, scale, rotation and flips. No recomputation of QC, clustering, annotation or differential expression is needed.
