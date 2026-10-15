# Produced-data visual guide

`build_visual_guide.py` reads completed native annotation, plot books and CA1 diagnostic tables. It uses existing exported Harmony UMAP coordinates for supplementary program views. It does not load a raw RNA matrix, aggregate nucleus RNA counts, fit an embedding or run a statistical model.

The main HTML starts with original-full counts, both tools and frozen original memberships. It provides population, tool, branch and scope controls, UMAP/marker thumbnail tabs and direct links to every PDF/PNG page. Reviewed confidence tiers describe molecular evidence; they are not calibrated identity probabilities. CA1 superficial/deep-associated and inhibitory workbook programs remain expression programs with undetermined physical position and morphology.

Paired original and CellBender views use the same retained nuclei and frozen original labels. Fresh paired partitions are separate secondary views. Native marker plots use all retained nuclei, with full-source library normalization. The negative-z glyph policy reduces dot area to one quarter and alpha to 0.18 without changing saved z or detection percentages.

The supplementary program UMAPs and CA1 correspondence heatmaps tabulate produced annotation/coordinate metadata. The branching diagram uses the produced retention audit. The CA1 diagnostic graphic uses the completed native diagnostic controls, excludes self-loops in degree/edge overlap, and does not claim an exact causal decomposition or global Harmony convergence.

Publication refuses existing nonempty output directories and requires the final annotation and global plot receipts. Bound annotation/plot outputs are checked against producer hashes, and coordinate inputs are checked against the native frozen-input manifest. `GUIDE_RECEIPT.json` preserves source, configuration, input and output SHA256 values. A completed effect-figure supplement is linked when its index and receipt are available; otherwise the guide states its absence.

Build with the dedicated Scanpy Python runtime, passing `--config config/analysis_v2_0_0.json --outdir reports/visual_guide`. Use absolute paths from the project root. `--diagnostics-only` creates only the technical figure from completed diagnostic tables.

The HTML is self-contained for chart thumbnails, labels and filtering and uses relative local artifact links. No CDN, analytics or external data transmission is used. The PDF/PNG/SVG artifacts remain independently usable. Browser interaction QA is unverified because the available browser tool rejected local `file://` access under its URL security policy; the builder has Python/JavaScript syntax checks, synthetic mask/hash checks, link checks and independent PNG/PDF inspection.
