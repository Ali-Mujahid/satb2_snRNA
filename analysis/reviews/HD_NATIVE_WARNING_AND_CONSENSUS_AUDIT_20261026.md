Takeaway: Failed: SATB2 export stopped after four complete coexpression fits; native warnings remain partly unclassified.

This follow-up examines coexpression programs in the frozen CA1 and inhibitory populations from six wild-type and six Satb2 conditional-knockout mice. All four pooled and consensus network fits produced complete native checkpoints, but the encompassing workflow failed during an SVG enrichment-figure export. The completed consensus outputs retain every selected nucleus and all twelve independent mice per endpoint and module, with finite scores oriented positively to mean module RNA. Numerous sparse nuclear gene vectors necessarily trigger the configured individual Pearson fallback, so nuclear hub correlations cannot all be described as purely robust biweight correlations. The saved global association tables contain exploratory mouse-level differences, while partial technical confounding, mixed inhibitory identity, changing module definitions, and unclassified native warnings limit biological interpretation. The justified recovery is an explicitly documented copy of verified completed outputs followed by corrected plotting and packaging, without refitting networks or changing their scores.

## Scientific context and question

The biological objective is an exploratory follow-up of WT versus Satb2 cKO snRNA-seq using exact frozen Seurat endpoints: CA1 deep-associated program (8,750 nuclei), superficial-associated program (11,706), inhibitory C1 (1,718), C6 (501), and C8 (468). Mice are the independent genotype replicates; nuclei and overlapping metacells are not independent genotype replicates. CA1 programs are molecular expression programs, not measured physical locations. C8 retains mixed atlas identity. This review does not establish direct transcription-factor binding, causal regulation, cell-cell signaling, or a new gene-level differential-expression result.

## Question, lifecycle, and endpoint

Reviewed source run: `/Users/mujahid357/Documents/codex_new/SATB2_REGULATORY_20261026/runs/20261026_hdwgcna_v1_0_1`.

The parent reports the registered run as FAILED. This agrees with `logs/enrich.log`: `ggsave` required unavailable `svglite` for SVG output and halted before enrichment completion or final workflow validation. The four native network checkpoints completed before that failure. No FINAL_VALIDATION receipt existed at this review. Any future recovery must be a separate run and preserve the failed run and its original configuration provenance.

## Key findings

Consensus checkpoint verification was read-only. Every recorded MD5 matched for CA1 (49 files, 213,406,806 bytes) and IN (49 files, 56,160,712 bytes); each checkpoint fingerprint matched its run identity. Every consensus nucleus ID matched the corresponding prepared lineage cells exactly, without duplication. CA1 has 192 unique mouse-by-endpoint-by-module rows (2 endpoints x 8 modules x 12 mice); IN has 144 (3 x 4 x 12). All native mouse scores are finite. Positive PC1 correlations with mean module RNA range from 0.9466 to 0.9982 in CA1 and 0.9380 to 0.9981 in IN; no eigengene harmonization was applied.

| Native network | Genes considered | Assigned modules | Unassigned grey genes | Metacells and genotype coverage |
|---|---:|---:|---:|---|
| CA1 pooled | 10,618 | 10 | 5,720 | Shared source metacells; native primary analysis |
| CA1 consensus | 10,618 | 8 | 6,442 (60.67%) | 600 WT / 600 cKO; six mice each |
| IN pooled | 11,065 | 5 | 8,761 | Shared source metacells; native primary analysis |
| IN consensus | 11,065 | 4 | 8,172 (73.85%) | 174 WT / 185 cKO; six mice each |

CA1 consensus sizes are turquoise 2,343; blue 477; brown 441; yellow 385; green 289; red 111; black 79; pink 51. IN consensus sizes are turquoise 1,228; blue 878; brown 604; yellow 183. Consensus soft-power topology gates passed with CA1 WT power 7 / cKO 6 and IN WT/cKO power 6.

The native pooled-versus-consensus contingency tables were independently checked against exact gene-set intersections. These are descriptive correspondences, not formal preservation statistics. For example, CA1 consensus red overlaps pooled magenta in 86/111 genes; consensus black overlaps pooled yellow in 69/79; consensus blue mixes pooled turquoise (245 genes) and pink (110). CA1 green has the closest descriptive one-to-one correspondence (250 shared genes; Jaccard 0.672). IN consensus turquoise overlaps pooled blue in 735 genes, whereas consensus blue overlaps pooled turquoise in 697; consensus brown contains 527/604 genes formerly grey. Module colors are arbitrary identifiers within a particular network and cannot be treated as matched biological labels across variants or lineages.

The four native global BH tables were checked arithmetically against their exported P values; maximum absolute discrepancy was below 4.5e-16. No statistical model or P value was refitted.

| Native table family | Tests | FDR <0.05 rows | Cohort |
|---|---:|---:|---|
| Pooled module x endpoint associations | 35 | 12 | 6 WT + 6 cKO |
| Pooled same-partition sensitivity | 35 | 5 | 6 WT + 3 cKO |
| Consensus module x endpoint associations | 28 | 5 | 6 WT + 6 cKO |
| Consensus same-partition sensitivity | 28 | 0 | 6 WT + 3 cKO |

The five pooled same-partition FDR-positive rows are CA1 black in both programs (positive cKO-minus-WT), CA1 purple in both programs (negative), and CA1 brown in the superficial program (positive). There are no FDR-positive IN same-partition rows. These continuous-score associations are exploratory and are not gene DEGs; absence of a sensitivity FDR hit does not prove the primary effect false because its cohort is smaller. Absolute PC1 magnitudes are not directly comparable between different modules or alternative module definitions.

## Interpretation of nuclear hub-correlation fallback

The existing `network_genes.tsv` exports provide endpoint detection fractions. These were converted back to positive-nucleus counts by multiplying each endpoint fraction by its authoritative nucleus count and rounding; the maximum rounding discrepancy was below 7.3e-12 nuclei. Counts were summed across endpoints, using 8,750/11,706 for CA1 and 1,718/501/468 for IN, to match the concatenated nucleus axis selected by native ModuleConnectivity. This is arithmetic QA of saved detection summaries, not an expression transformation or fit.

| Lineage/network | Assigned genes with detection <50% | Top-20 hubs per module with detection <50% |
|---|---:|---:|
| CA1 pooled | 3,330/4,898 (67.99%) | 34/200 (17.00%) |
| CA1 consensus | 2,855/4,176 (68.37%) | 33/160 (20.63%) |
| IN pooled | 1,424/2,304 (61.81%) | 23/100 (23.00%) |
| IN consensus | 1,928/2,893 (66.64%) | 26/80 (32.50%) |

There were no assigned genes exactly at 50% detection in these comparisons. Strictly below 50% detection means more than half of the nonnegative log-normalized nuclear values are zero, which guarantees zero gene-vector median absolute deviation. Other distributions can also have zero MAD, so these fractions are lower-bound diagnostics rather than exact event counts.

The pinned installed hdWGCNA ModuleConnectivity implementation with `sparse=FALSE` calls WGCNA signedKME on nuclear RNA and eigengenes. The frozen workflow supplies `corFnc='bicor'`, `maxPOutliers=0.05`, and `pearsonFallback='individual'`; native bicor defaults to robust weighting on both axes and issues a zero-MAD warning when that fallback is used. For a zero-MAD gene vector, the gene-side robust treatment falls back as configured; the eigengene side can remain robust. The proper description is **nuclear kME calculated with bicor and individual Pearson fallback**, not uniformly pure robust bicor and not necessarily pure Pearson on both axes. The nuclear detection audit does not quantify fallback use in aggregated metacell network discovery or metacell graph edges.

Some module hub rankings are particularly affected by this diagnostic: CA1 pooled pink has 11/20 low-detection hubs, IN pooled green 11/20, CA1 consensus blue and brown each 10/20, and IN consensus brown 14/20 and turquoise 11/20. This limits claims of uniform robustness; it does not by itself invalidate these modules or justify retrospectively changing their ranks.

## Artifacts and provenance

Results reviewed include each consensus `module_membership.tsv`, `module_eigengenes_mouse.tsv`, `nucleus_eigengenes.tsv.gz`, `eigengene_orientation.tsv`, `consensus_coverage.tsv`, `selected_power.json`, `pooled_consensus_module_overlap.tsv`, and corresponding pooled membership/detection/hub tables. Cohorts were checked against `results/prepared/lineages/{CA1,IN}/cells.tsv` and the frozen membership audit. Native seals are each branch's `completed.json` and `run_identity.json`; no RDS objects were loaded or modified. Global comparisons are the four `results/hdwgcna/{pooled,consensus}_module_genotype_{associations,same_partition}.tsv` tables. Existing figures and all native source/results remain unchanged by this review.

The installed implementations were inspected by dumping function definitions only; no expression calculation, network fit, genotype association fit, normalization, or new enrichment was invoked.

## Limitations and blockers

`logs/hdwgcna.log` ends with `There were 50 or more warnings` but does not preserve individual warning types or a precise count. The post-exit warning list cannot be reconstructed from the saved log. A separate read-only namespace-inspection session emitted ten import-replacement warnings; source also contains an expected TOM-not-found warning path when TOMs are intentionally not stored. These observations do not establish which warnings occurred in the original run and must not be used to classify its entire warning total. Original native warnings remain partly unclassified.

The enrichment log separately records an already-existing output-directory warning and a text-glyph substitution for the Unicode minus sign, followed by the missing-svglite SVG error. These export warnings are distinct from the earlier unknown network warning set. Enrichment tables/one PDF may exist, but the failed run has no successful enrichment completion seal and is not a publication-ready complete workflow. GO and DEG-enrichment conclusions require a successfully completed, validated recovery. GO terms from primary pooled module definitions cannot be transferred to consensus modules solely by matching color names.

The native module cross-reference also retains an important original gene-level observation: Satb2 RNA is higher in cKO in both CA1 programs (deep log2FC +1.34993, original padj 2.08e-48; superficial +1.09451, padj 8.64e-19), and Satb2 maps to pooled black. This is reused original DESeq2 evidence, not a newly fitted contrast, and does not verify SATB2 protein loss. Targeted exon coverage, residual transcripts, Cre/flox targeting or mosaicism, and protein depletion remain independent validation questions; no specific explanation follows from gene-level RNA alone. This also reinforces the separate caveat for curated activity priors that include a Satb2 self-target edge.

Consensus uses the same mice and source metacells as the primary network and is not independent validation. Genotype is partially confounded with technical partition; same-partition sensitivity retains only three cKO mice. C8 has mixed identity and one 13-nucleus mouse stratum excluded from metacell discovery, while all selected nuclei remain projected/scored. High grey fractions and changed assignments remain explicit. RNA coexpression and curated prior activity cannot prove anatomical layer, direct binding, causal SATB2 regulation, receptor activation, or physical cell-cell communication.

## Recommended next step

Preserve this failed run. Recover export/packaging in a new native workflow that first proves byte-identical reuse of prepared inputs, all four completed network trees, the original scientific configuration, and the original fitting source. Correct plotting through explicit installed graphics devices and rerun only downstream enrichment export/packaging, then validate the final receipt and source-aware archive before publishing. Keep the original unknown-warning disclosure and fallback QA in the result review; do not refit merely to recover missing warning messages.
