# Revised CA1 layer volcano reader

The public entry is `site/ca1_layer_volcano/index.html`, with `#superficial` and `#deep` direct links. The main atlas, CA1 program page and historical LR volcano page link to it.

## Frozen scientific inputs

Seurat only; original-full cohort; cKO versus WT; six independent WT and six cKO mice. Superficial is the existing **25 November C2+C3 refit**, excluding CA1 C8, with 14,206 modeled genes and 1,775 DEGs at adjusted P <0.05 (938 up, 837 down). Deep remains the original C1+C5 fit, with 13,439 modeled genes and 742 DEGs (360 up, 382 down). Complete original TSV text is archived in `*_source.tsv.gz`; coverage tables and all source/decoded/compressed SHA256 values are retained. Native unavailable values are preserved as null/NA, never zero. Stable gene IDs are the row key. No statistical model is fitted by the builder or reader.

`scripts/build_ca1_layer_volcano.py` reads these existing files and packages them into JSON gzip assets. The reader fetches assets from the immutable release tag in `PROVENANCE.json`, verifying compressed and decoded SHA256. Keeping model assets outside `site/` preserves the unchanged Pages size guard. The older LR explorer retains its historical C2+C3+C8 definition; its new link identifies this distinction.

## Visual specification

Question: which genes change in cKO versus WT within each exact CA1 layer program? Scatter: X is signed log2(cKO/WT), Y is −log10 adjusted P by default, optionally raw P. Colours use gene-level adjusted P and the selected absolute log2FC cutoff even with raw-P Y. Default |log2FC| cutoff 0 preserves the saved DEG counts; ±1 dotted guides are independent reference positions, explicitly labelled. Auto axes are fit-specific. For panel comparison, choose equal manual limits. Points outside axes are omitted and counted, rather than clamped to boundaries. Zero P values use a stated numerical display floor of 1e−300. This is a display of estimates and hypothesis tests, not a new test between layers.

Web/vector figure 900×650 by default, adjustable with independent 15px gene-name size. Blue/red/grey uses the established atlas/user mapping; signed X position redundantly encodes direction. White field, quiet axes and optional enclosing border. Top 10 visible significant genes are ranked by adjusted P then |log2FC|; deterministic displacement and leader lines avoid label overlap. All requested/shown/omitted label counts are exposed. Additional named genes are optional. CSV-style TSVs retain exact model values; no genome-wide re-ranking statistics are computed.

SVG, vector PDF, 3× PNG, current DEG lists, complete model tables, labelled gene lists and validated plot-state JSON are supported. Plot colour cutoff, FDR cutoff and independent guide position are recorded in settings. Model/source hashes prevent loading settings against a different fit.

## QA

`node tests/ca1_layer_volcano.cjs` checks source and payload hashes, both cohort definitions, all saved DEG totals, exact default label counts, colour and threshold settings, axis/border/guide geometry, missing/zero semantics and immutability. The Pages workflow invokes it. `--browser` uses headless Chrome with byte-identical local payload substitution for the immutable raw-GitHub URLs, and exercises both datasets, exports, invalid-axis rejection, reset and state roundtrip. Actual remote loading and PDF downloads are checked after deployment.

The default two figures were inspected at native delivery size, including extrema and label crowding, plus independent PDF rendering and grayscale inspection. Unchanged CA1 reassignment provenance validation passes. No pre-existing fit, embedding, LIANA result or annotation is altered. Browser acceptance artifacts remain under `browser_checks/`. The design follows Scientific Visuals & Tables scatter/label/threshold grammar and the user's established red/blue palette; style-reference retrieval selected OAI-VIS-008, OAI-VIS-085 and OAI-VIS-089.
