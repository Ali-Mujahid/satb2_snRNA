# Hosted atlas browser acceptance — 15 October 2026

The live HTTPS atlas was operated in the Codex in-app browser at its normal 1280-pixel desktop viewport. This is actual UI/download testing, separate from the VM, byte-identity and scientific-source checks. The original local analysis and offline v1/v2 directories were not edited.

## Passed UI workflows

- Loaded both original-coordinate tool views, including all nuclei (107,747 retained), excitatory (60,908), CA1/ProS (27,703) and inhibitory (8,911). Checked paired-original All (103,983) and paired-CellBender inhibitory (8,903) views as well.
- Switched Harmony/unintegrated coordinates, MapMyCells labels, genotype filters, continuous palettes, and saved gene-summary UMAPs. The final UI keeps Calb1 selected when its UMAP is requested. It explicitly uses cluster means, not individual-nucleus RNA.
- Selected supplied-workbook/canonical marker presets, searched the 878-gene panel, and used 60, 90, 120 and 180 genes per page. A 139-gene selection exported every gene range for both tools. Manual color limits [-1,2], WT/cKO summaries and smaller/fainter negative-z dots were exercised; source values are unchanged.
- Selected Malat1 outside the dotplot panel from the 32,285-source-gene mouse selector, added Calb1, and exported raw mouse boxplots and values. Also exported normalized Seurat inhibitory boxplots. One point represents one independent mouse; factors stay bound to the exact saved endpoint.
- Edited labels with actual keyboard input in both tools. Custom label unions appear in the DE selectors. Applied provisional suggestions, verified their mixed/“-like” wording, and restored source labels.
- Exported settings, reset the view, imported the exported JSON, and verified the custom union and its inference status returned. Annotation CSV export also downloaded.
- Excluded and restored a cluster in every population for both tools (matrix below), checked exclude-all gives NO NUCLEI, and checked an unmatched edited total/two-cluster union gives NO MATCHING FIT rather than another fit's statistics.
- Verified total, cluster, CA1 layer-associated program and leave-one-mouse-out DE selectors. The latter showed 6 WT/5 cKO mice when omitting cKO 1138. An 82-nucleus inhibitory cluster retained all 82 with 6 WT/6 cKO and an exploratory coverage flag; a 126-nucleus marker-defined category had 3 WT/3 cKO.
- Downloaded the exact CA1 total DESeq2 source through the browser's cross-origin, pinned data path. It has 16,255 rows and exactly the standard seven fields with gene symbols. Checked gene search and FDR/fold-change controls. The final displayed Adamts3 P value is 6.929e-110, with the exponent intact.
- Downloaded UMAP PDFs for both tools, both multipage dotplot PDFs, mouse boxplots, volcano/comparison plots, and an inhibitory spacing view. At 0.65× the spacing export explicitly says DISPLAY LAYOUT / distances altered; reset restores original geometry.
- Operated the download index: 96 total assets, 28 processed objects, 13 Seurat objects after searching, and four archives. All six gallery images loaded; a gallery PDF downloaded from the public GitHub release and matched its declared hash.

## Exclusion matrix

Original-full baseline cohorts. Each exclusion was undone using Select baseline.

| Population | Baseline, both tools | Scanpy after first-cluster exclusion | Seurat after first-cluster exclusion |
|---|---:|---:|---:|
| All nuclei | 107,747 | 105,296 | 92,577 |
| Excitatory | 60,908 | 58,467 | 52,534 |
| CA1/ProS | 27,703 | 25,248 | 21,378 |
| Inhibitory | 8,911 | 8,431 | 7,193 |

## Corrections and evidence

The audit found and corrected selected-gene replacement, scientific-notation exponent trimming, category-shortcut persistence, and stale frontend asset caching. On the final deployed app (content hash cefbc5906fee), clicking Inhibitory and reloading retained that category and 8,911 nuclei in both tools. Scientific payload bytes and exact-source text were not changed. Hosted JavaScript URLs now contain content hashes. The corrected offline runtime has the corresponding UI fixes; its original classic transport remains intact. Original archive assets remain available separately.

Screenshots: [gene UMAP comparison](browser/atlas_gene_umaps.jpg), [figure gallery](browser/figure_gallery.jpg). PDF/source file integrity is recorded in [browser_download_integrity.json](browser_download_integrity.json). Data/version checks are in [remote_data_integrity.json](remote_data_integrity.json), [release_integrity.json](release_integrity.json), and [offline_ui_update.json](offline_ui_update.json). The final unauthenticated HTTP and Actions checks are in [deployment_acceptance.json](deployment_acceptance.json).

## Practical scope

This tests control families and representative saved endpoints, not every possible gene/endpoint/settings combination or a separate Windows/Chrome/Firefox installation. The full-nuclei redraw can take several seconds; two automation click deadlines were followed by checking the actual panel state, which confirmed the actions completed. No browser JavaScript error was observed in the final tested session. Large legends continue onto additional PDF pages, so a three-marker-page selection may produce more than three PDF pages.

The browser policy prevented inspecting the local file URL directly; the requested GitHub-hosted HTTPS atlas was tested instead. The corrected offline derivative was checked through its classic-loader, UI regression, ZIP CRC and unchanged-file identity tests, not by controlling a local file tab. New exclusions or unions without saved memberships cannot run a new DESeq2 fit inside this static webpage. Fine identities and physical CA1 layers retain the source uncertainty.
