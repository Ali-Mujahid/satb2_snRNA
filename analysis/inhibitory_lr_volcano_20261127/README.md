# Original inhibitory C1/C6 volcanoes with CA1 partner highlights

The atlas reader at `site/ca1_layer_volcano/` now offers four saved Seurat populations: superficial CA1 C2+C3 (C8 excluded), deep CA1 C1+C5, original inhibitory C1, and original inhibitory C6. This addition packages existing fits for display; no classification, DESeq2 model or LIANA analysis was rerun.

| Saved inhibitory population | Modeled genes | Nuclei | WT / cKO mice included | DEGs, FDR <0.05 | Higher / lower in cKO |
|---|---:|---:|---:|---:|---:|
| C1, PV reference | 12,769 | 1,718 | 6 / 6 | 27 | 11 / 16 |
| C6, original SST parent cluster | 8,989 | 501 | 6 / 6 | 22 | 4 / 18 |

`PROVENANCE.json` identifies the original stable-ID audit tables and the immutable public payloads. The native result tables, all-12-mouse inclusion audits and original model settings are saved alongside the compact JSON payloads. Missing estimates and adjusted P values remain null, and the browser verifies both compressed and decoded checksums. C6.0/C6.1/C6.2 and mouse-subset fits are separate analyses and are not substituted here.

## Ligand–receptor mapping

The existing 14 reviewed, directed population-pair records are rendered from both endpoint perspectives. `LR_VIEW_DATA.json` therefore contains 28 view records, representing the same 14 candidates, not 28 independent interaction predictions. `IN_CA1_LR_endpoint_mapping.tsv` records both stable IDs, both exact fitted effects, source/target roles, native LIANA scores, and the existing molecular and synaptic evidence annotations.

- **IN C1:** Nlgn1, Grm7, Ptprz1 and Fstl5; six superficial-partner rows and one deep-partner row. All four inhibitory genes are higher in cKO in the saved fit.
- **IN C6:** Cntn3, Robo1 and Slit2; three superficial-partner rows and four deep-partner rows. All three inhibitory genes are lower in cKO in the saved fit.

Inhibitory volcanoes plot the inhibitory gene's own fold change and P value. Their CA1 partner table uses the revised C2+C3 superficial or C1+C5 deep effects. Resource orientation is retained separately from synaptic direction. APP–GRM7 and PENK–GRM7 remain unverified; screen-supported, indirect and context-dependent relationships retain their original qualifications. A highlight does not establish a physical contact in these mice.

## Controls and exports

All four populations share physical-size and pixel-size controls, independent dot and label sizes, top-N/custom/selected-LR labels, CA1/IN partner filters, reviewed-evidence and native-rank filters, editable colors, opacity, axis limits and ticks, fold-change cutoffs and independent dotted guides, and plot boundaries. SVG, vector PDF, 3× PNG with physical-resolution metadata, DEG/all/label tables, filtered LR evidence and settings JSON use the currently selected population. Labels and settings use stable gene IDs.

## Reproduction and verification

Run `scripts/build_inhibitory_lr_volcano.py` after the existing CA1 payload and LR builders; it extends the frontend manifest while leaving the two frozen CA1 payloads unchanged. The current release tag is recorded in `PROVENANCE.json`.

- `node tests/inhibitory_lr_volcano.cjs`: native checksums and inclusion audits, exact inhibitory endpoint coordinates, counts, partner filters, binding flags and source immutability.
- Add `--browser` for Chromium checks of C1/C6, all eight downloads, label selection, partner filters, small-panel reflow, JSON roundtrip and both CA1 regression panels.
- Existing `tests/ca1_layer_volcano.cjs` and `tests/ca1_volcano_lr_size.cjs` pass with the expanded reader.
- `scripts/verify_site.py` passes with all 3,140 original shards preserved, no scientific model fitting, and the original 950 MB site guard intact.
- Independently inspected vector PDF renders and exported PNG/SVG dimensions: the 480×380 px small preset exports as 360×285 PDF points, 127×100.54 mm SVG, and 1440×1140 PNG pixels at approximately 288 dpi. Dot radius remains 2.4 px and gene labels 15 px through figure resizing. The zoomed LR views disclose the visible-gene denominator and full-fit DEG counts.

Browser and export acceptance records and representative compact LR figure exports are retained under `browser_checks/`.
