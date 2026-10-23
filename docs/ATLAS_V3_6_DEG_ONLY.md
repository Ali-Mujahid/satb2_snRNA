# Atlas v3.6.0: corrected DEG-only communication

The active communication explorer now requires every ligand AND receptor complex subunit to be an unambiguous DEG (`padj < 0.05`) in its own exact source or target endpoint. Only P0244 CA1 deep-program, P0247 CA1 superficial-program, P0140 inhibitory C1, P0145 C6 and P0147 C8 are eligible. The old OR rule and inhibitory-total/other-cluster selectors cannot be restored using saved settings.

The corrected selection has 16 pooled, 15 WT and 16 cKO native records; rank≤0.05 retains 8/9/7. The additional replicated/context-specific tier retains two C6 Cntn3→CA1 Ptprg records. C8 remains explicit with no qualifying edge because its only DEG Calb1 is absent from mouseconsensus.

Gene up/down directions use the original matching DESeq2 log2(cKO/WT), separately for ligand/source and receptor/target. Controls filter each direction or the four ordered combinations, hide/show signed annotations, and recolor up/down/mixed. Dotplot numeric color and detection area are preserved; paired rims carry gene signs. Gene-circle outlines and labels show gene signs. Aggregate heatmap/circle/circos summaries count sign combinations without assigning net interaction activity.

Seven-column native LR tables are unchanged projections; signed extended tables add gene effects and adjusted P. Original seven-column DEG copies are byte-identical, with separate up/down subsets. The standalone HTML embeds only the corrected 47 records and exact component evidence. Corrected precomputed figures are packaged separately from historical broad figures.

This update is a deterministic projection of completed LIANA results. It does not refit DESeq2/LIANA or renormalize a truncated DEG matrix. Native scores/ranks retain their original full-gene and population-background meaning. The previous data payload, native score fields, original gene results, code and archives remain preserved.

See [detailed results](DEG_ONLY_RESULTS.md), [endpoint audit](DEG_ONLY_ENDPOINT_AUDIT.json), and the source/read-only checks in `tests/communication_deg_only.cjs`, `tests/communication_deg_ui.cjs`, `tests/communication_deg_visuals.cjs` and `tests/communication_app.cjs`.

Validation before publication: complete native data/transport and original atlas integrity checks passed; strict-data, policy, chart, network, sign-visual and application acceptance suites passed. Fifteen precomputed PDF panels passed parsing, roundtrip rendering and visual inspection. See the figure, portable and release receipts. Live deployed-browser checks are recorded separately after release.

Live HTTPS acceptance passed after deployment: strict default16, native rank≤0.05 eight, WT/cKO15/16, CA1→IN6 and IN→CA1 ten, four signed combinations3/2/4/7, C8zero with Calb1 explanation, exact Cntn3↓/Ptprg↑ component evidence, editable sign colors and annotations, all five PDF downloads and exact full filtered table downloads. Five live PDFs parsed/rendered successfully; dotplot/gene-circle visual inspection passed alongside all15 precomputed panels. No browser errors/warnings. See [live acceptance](CCC_DEG_ONLY_LIVE_ACCEPTANCE.json), [PDF checks](CCC_DEG_ONLY_LIVE_PDFS.json), and [screenshot](images/communication_v360_live.png).
