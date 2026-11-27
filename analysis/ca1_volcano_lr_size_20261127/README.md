# CA1 volcano LR highlights and fixed-mark resizing

The existing CA1 revised-layer reader is updated in place. Native fit payloads, their release-pinned URLs and SHA256 values, sample definitions, DESeq2 estimates and LIANA ranks remain unchanged. Superficial = C2+C3; deep = C1+C5; Seurat original-full 6WT/6cKO.

## Source lock

`scripts/build_ca1_volcano_lr.py` joins the 14 already-reviewed rows in `analysis/synaptic_evidence_20261125/results/CA1_IN_directed_pair_roles.tsv` to the uniquely mapped stable gene IDs in the exact current fits. It checks the native CA1 fold changes and adjusted P values against the rounded archived table. There are nine superficial rows/seven CA1 genes and five deep rows/five CA1 genes. Both endpoints pass the original FDR <0.05 criterion. Gene-level outlines do not merge or re-score LR rows. The full source mapping, roles, partners, ranks, binding status, corrected molecular/synaptic interpretations, compartments and reference URLs remain downloadable. IN C6 refers to the original parent cluster, not C6.0. No new fits, LR prediction or literature review is performed.

APP–GRM7 and PENK–GRM7 remain unverified. DCC–SLIT2 remains indirect; FSTL5–CNTN5 retains binding-screen support. NRXN1–NLGN1 reciprocal resource entries retain their separate resource roles and their shared reviewed CA1-to-IN predicted synaptic orientation. Local contacts remain untested. A dagger indicates that a gene has at least one currently included screen, indirect or unverified association; the table distinguishes these per pair. The candidate highlight is not a direct-binding or functional-contact claim. Only the CA1 endpoint is plotted on its own CA1 volcano; the inhibitory partner is table metadata.

## Figure behavior and visual specification

Claim: locate exact CA1-end genes from the saved LR candidates within the existing genotype-effect estimates, at a chosen final output size. X/log2FC and Y/−log10 P coordinates remain exact. Stroke outlines identify candidates while fill retains the genotype-effect colour. Filters choose resource direction, inhibitory partner, binding-evidence level and native LIANA magnitude rank. Checkbox selections choose LR text labels independently of top-N and custom labels. LR labels are placed first, with collision avoidance and leader lines. Omitted labels, unavailable coordinates and cropped points are counted. Zoom to LR genes changes axis limits only and switches off automatic top-N labels; Auto axis limits restores the full extent.

Width/height are expressed in millimetres or CSS pixels. Small (480×380), presentation (900×650), and poster (1200×850) presets reflow the plot; none changes dot radius or gene-name size. Titles, legends and long notes are optional, with a minimal whole-model/cropped-view note retained. Legend/header text wraps instead of shrinking. Invalid combinations that leave insufficient data space produce an explicit message and disable exports. Selecting fewer gene labels remains useful for compact plots.

SVG has explicit physical millimetre dimensions and unchanged internal CSS-pixel mark sizes. PDF uses jsPDF px_scaling (72/96 point per CSS pixel), so it has the same physical dimensions and dot/glyph size as SVG. PNG is 3× with a validated pHYs chunk at 288dpi, preserving the same intended physical size. Whole-object scaling in a downstream slide editor scales every element; the reader therefore instructs the user to choose final size before export. Gene label size remains a separate control.

The established red/blue/grey atlas palette is retained; dark teal outlines and redundant text/dagger flags add the candidate context. Scientific Visuals & Tables governs the scatter, labeling and source-integrity checks. The prior scatter references and exact numerical data remain applicable.

## Validation

`node tests/ca1_volcano_lr_size.cjs` checks source hashes, exact CA1-ID mapping, all row/gene totals, stable data and DEG values, small/large dot-radius and gene-font invariance, label subsets, filters, unverified flags and SVG physical sizing. The Pages workflow runs it alongside the existing fit validator. `--browser` checks both populations, input millimetres, all presets, filter and checkbox actions, focus, exports and settings roundtrip using byte-identical local copies of the immutable model payloads.

Final SVG/PNG and independently rendered vector PDF outputs were checked for clipping, colliding labels, notes, grayscale interpretation and correct model identity. PDF MediaBox and PNG dimensions/DPI were checked independently with pypdf/Pillow. The archived browser exports and acceptance record are in `browser_checks/`. The unchanged 950MB Pages guard and all frozen source shards are checked before publication. Live source transport and PDF downloads are checked after deployment.
