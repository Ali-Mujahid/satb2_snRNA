# Atlas v3.4.1 live acceptance

Verified the deployed [GitHub atlas](https://ali-mujahid.github.io/satb2_snRNA/) after [successful deployment](https://github.com/Ali-Mujahid/satb2_snRNA/actions/runs/39128737990).

- All three one-click presets select 99, 118 and 97 available genes in supplied verified order. Unavailable symbols, four documented legacy-name mappings and duplicates remain auditable.
- Downloaded marker CSVs contain 3,663 inhibitory all-cell-panel rows, 2,360 CA1/ProS-panel rows and 1,940 sublayer-panel rows. Every count, detection fraction, expression mean and produced z equals its native saved summary exactly.
- Nine PDFs downloaded successfully: marker panels, Wfs1 measured-expression UMAP, packed C4 and both rounded schematic UMAPs. Each complete plot occupies one wide figure page. Long label legends can continue onto a separate text-only page.
- Wfs1 measured individual-nucleus expression loads in both tools with the same full-library source-count normalization and shared range.
- Packing, Seurat C4 arrow movement and rounded schematics operate on the live page; adjusted geometry is labeled.
- Display-coordinate downloads retain exact original coordinates for all 17,822 inhibitory nuclei in both tools. Reset recovers them exactly. The new screen-position regression verifies C4 moves while the other 8,318 Seurat nuclei stay fixed.
- Browser console reports no warnings or errors. Prior settings are backed up and restored exactly. Previous scientific results remain unchanged; no models were fitted.

[Preset controls screenshot](screenshots/atlas_v341_one_click_panels.png) · [Packed layout screenshot](screenshots/atlas_v341_packed_inhibitory.png) · [Rounded schematic screenshot](screenshots/atlas_v341_rounded_schematic.png)

See the [usage guide](ATLAS_V3_4_USER_PANELS.md), [machine-readable browser receipt](atlas_v34_browser_acceptance.json), [export validation](../analysis/dashboard/v34_user_panels/export/VALIDATION_REPORT.json) and tests for reproducible details.
