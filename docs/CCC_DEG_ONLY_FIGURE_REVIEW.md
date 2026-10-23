# Corrected DEG-only figures: verification

Status: PASS after final spacing correction, 2026-10-23.

- Exactly 15 views were generated: five chart types for each of the three requested directions. Each view has matching SVG, 192-dpi PNG, and one-page PDF.
- Native pooled records were filtered without refitting: 6 CA1-to-inhibitory, 10 inhibitory-to-CA1, and 16 bidirectional records. Gene circles show every retained population-specific record. The C8 zero-edge population remains visible in population-level plots.
- All 15 PDFs were reopened with pypdf, checked for one page and the intended dimensions, and rendered with Poppler. The complete PDF contact sheet was inspected.
- The final direct PNG and PDF contact sheets show no clipped labels, missing glyphs, or overlap between population labels and the signed endpoint summary. The previously close C6 header label was moved below the summary.
- The CA1-to-inhibitory gene-circle PDF was additionally inspected at full preview scale. Up/down arrows, colored endpoint outlines, population fills, native-expression edge widths, and legends are intact.
- Dot interior area still encodes minimum all-subunit detection. Metric fill is independent of signed left/right rims. Heatmap zero cells remain distinct from untested NA cells.
- Up/down denotes matching endpoint DESeq2 log2(cKO / WT), not inferred up/down interaction activity. Network arrows denote the listed ligand sender and receptor receiver.
- Network renderer regression checks (16 checks) and the separate signed-visual suite passed after the final correction. Population-label header clearance was checked at font sizes 8, 12, 20, and 30.

Figure manifest SHA-256: `e150a4e26a780d8186816b3287f175a88a2add50a368985482e035e0508011d9`

PDF roundtrip contact sheet SHA-256: `35259c72991aaae37df586ebdf93eca388a922625426609be3a7c05685b0f872`

The manifest contains the full 45 figure-file hashes, source code/data hashes, all plotting defaults, source TSV references, signed endpoint-combination counts, and the 15 PDF roundtrip checks. `tmp/pdfs/` contains intermediate roundtrip PNGs and can be omitted from a user download archive.
