# SATB2 TF activity and motif-supported regulons — 26 October 2026

The user selected full pySCENIC, CollecTRI TF activity and hdWGCNA, then confirmed that no cached mouse cisTarget references are available and instructed us to continue the available analyses. This addendum does not change the already registered hdWGCNA scientific plan or existing strict DEG-only communication analysis.

## Cohort and distinct questions

Use exactly the frozen Seurat CA1 deep-associated and superficial-associated expression programs and inhibitory clusters 1, 6 and 8, across all six WT and six cKO mice. The frozen native population membership identities and original integer pseudobulk counts must match the completed hdWGCNA preparation audit. Layer expression programs do not establish physical laminar location; IN_C8 retains its mixed-identity warning.

Ligand–receptor views remain restricted to significant original DESeq2 genes at both exact endpoints, including every complex subunit. Regulatory activity and coexpression require the full eligible measured gene background; using only significant DEGs would bias these methods. No original gene DE fit, annotation, count or communication result is replaced.

## Available CollecTRI analysis

Use native decoupler 2.3.0 ULM with its frozen mouse-mapped CollecTRI signed network, obtained through native `op.collectri(organism='mouse', remove_complexes=False, license='academic')`. Preserve the downloaded human reference and HCOP orthology files, actual checksums and mapping policy. This reference includes 1,165 TFs and 43,226 signed links; SATB2 has 33 source links before measurement and eligibility filtering. Cross-species curation is not evidence of direct binding in these mice.

For each endpoint and mouse, normalize full original integer counts with the saved size factor from that exact original DESeq2 fit and apply log1p. Retain genotype-neutral genes with at least five counts in at least three of all twelve mice. Exclude missing or ambiguous symbols instead of guessing mappings. Keep non-target genes in the ULM background. Require at least five eligible prior targets per TF. Preserve native ULM t-values and native within-input gene-set FDR separately.

Compare continuous mouse activity scores using `score ~ condition`, with one observation per independent mouse and a cKO-minus-WT effect. Adjust across all tested endpoint–TF associations using one global BH family. Report a prespecified same-partition sensitivity using the six WT and three cKO mice from `23HJT3LT4_lanes6_7`, with its own BH family, plus leave-one-mouse-out effect signs and 200 stratified mouse bootstrap draws. Scores are not integer RNA counts: these continuous-score models are separate from gene differential expression and are never called DEGs or fitted with DESeq2. Original DESeq2 remains the sole gene DE method.

As a separate supplement, native ULM can use all finite, unambiguous original DESeq2 Wald statistics, including zeros and non-DEGs. Its gene-set enrichment FDR is not the independent-mouse activity-association FDR. Display and download these result types separately.

The frozen prior includes a Satb2-to-Satb2 self-target link. Export explicit self-target flags. Preserve the full-prior main analysis and run a separately named SATB2-only sensitivity that removes this single link, retains the same full eligible RNA background and requires at least five remaining eligible targets. Record unavailability if coverage fails. Compare its mouse scores in a separate five-endpoint BH family. Full-prior SATB2 activity is partly dependent on its own RNA and is not independent confirmation of SATB2 regulation.

## Requested pySCENIC branch

The isolated pinned pySCENIC 0.12.1 runtime is ready. Full motif-supported inference requires GRNBoost2, matching mm10 v10-cluster gene-based cisTarget rankings and matching MGI motif-to-TF annotations, followed by AUCell. Both official cisTarget hosts refused connections before TLS or HTTP; no reference bytes were retrieved and the user has no cache. This branch is blocked on reference acquisition. Do not substitute human, mm9, region-based or unverified files, or label curated TF activity as de novo motif-supported regulons. Do not force a SATB2 regulon where the resource does not support one. No ATAC data or direct binding evidence are available.

## Execution, deliverables and interpretation

Reuse the user-selected project Snakemake controller. Register a new immutable, hash-bound TF workflow and launch only through the NGS Workbench native execution path. Keep installed software and resource acquisition outside run preparation. Use fresh outputs, no network access inside scientific tasks, and exact before/after source input checks. Record Python, R and package versions and the limits of the existing controller/base library closure.

Export native scores, target coverage and signed prior edges, mouse-level association and sensitivity tables, original-Wald enrichment, source/configuration receipts, PDF/SVG/PNG figures and a separately validated offline atlas payload. Add results to the GitHub atlas only after actual output and figure review. Keep native result units unchanged under display scaling. Show pySCENIC's blocked status explicitly. Partial genotype/technical-partition confounding limits causal attribution; target-based inferred activity is not TF RNA abundance, TF occupancy or demonstrated regulation.

Method documentation: https://decoupler.scverse.org/en/stable/api/generated/decoupler.mt.ulm.html and https://decoupler.scverse.org/en/stable/api/generated/decoupler.op.collectri.html .
