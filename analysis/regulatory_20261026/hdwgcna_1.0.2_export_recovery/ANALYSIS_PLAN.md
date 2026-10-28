# SATB2 regulatory and co-expression follow-up — 2026-10-26

## Objective and decision

Identify reproducible co-expression modules and candidate regulatory programs in the frozen Seurat CA1 layer-program and inhibitory populations, and relate them to the existing strict DEG-only CA1→IN, IN→CA1 and bidirectional communication results. This is an additive exploratory follow-up; prior counts, labels, DESeq2 fits and published communication results are immutable inputs.

## Observed starting point

- Twelve independent mice: six WT and six Satb2 cKO; nuclei and metacells are observations within mice, not independent genotype replicates.
- Original integer UMI counts and frozen tool-specific memberships are available. Exact populations: CA1 deep program P0244, superficial program P0247, inhibitory C1 P0140, C6 P0145 and C8 P0147.
- Expected exact memberships from the native audit: DL 8750, SL 11706, C1 1718, C6 501 and C8 468 (23143 total). A prior informal C8 count of 568 was inconsistent with the authoritative sample audit and is superseded. Preparation must verify rather than force these totals.
- Mouse reference is mm10. Exact original GTF release is not established. Stable-ID/symbol joins and duplicated-symbol collapsing must be audited.
- All five original DESeq2 fits include all 12 mice. Reuse their saved gene statistics/size factors; no new gene differential-expression model is proposed.
- Technical partition partially confounds genotype. Six WT and only three cKO share the known lane6/7 partition; two cKO are lane8 and one is unverified. Primary continuous-score associations cannot be interpreted as causally separated from all technical effects.

## Methods and comparison structure

### hdWGCNA — selected by the user

Pin hdWGCNA 0.4.12, immutable upstream commit e3344d1f7bbac4264adf94f7aa31e0802fa8282d. Build separate CA1 (DL+SL) and inhibitory (C1+C6+C8) networks using normalized RNA expression and eligible expressed genes, not a DEG-only or HVG-only network universe. Existing Harmony/UMAP coordinates are not expression inputs.

Normalize using the full raw gene-count denominator. Select network genes expressed in at least 5% of nuclei in any corresponding endpoint. Build metacells within mouse × endpoint; never mix mice or genotypes. Prespecified CA1 settings: k25, max_shared10, min_cells50, target50 per group. Inhibitory: k10, max_shared3, min_cells20, target25. Small-stratum discovery exclusion is recorded, while all original nuclei remain available for fixed-module projection and mouse summaries. In particular the 13-nucleus C8/cKO933 stratum is too small for this metacell discovery setting.

Use a signed bicor network, inspect soft-power diagnostics, and bound block size at 3000–4000 genes for the available memory/storage. The selected power and the full diagnostic table must be exported. A poor topology fit is a failed validity gate, not permission to hide diagnostics or force attractive modules. Keep grey/unassigned genes explicit. Export module membership, kME/hub ranks, dendrograms, eigengenes and network-edge summaries. Module edge weights describe co-expression, not direct TF binding or physical contact.

Primary genotype-blind pooled modules are fixed across WT/cKO. Consensus/preservation analyses are conditional on adequate metacell coverage and resources; do not label bootstrap/leave-one-mouse-out correlations as formal WGCNA preservation statistics. Report the actual performed method. Scores are projected to all nuclei, then summarized once per mouse × endpoint. Continuous module-eigengene genotype associations use a clearly separate mouse-level linear model with FDR over all tested module/endpoint associations; they are not DEGs and must never be supplied to DESeq2. Add same-partition sensitivity (6WT/3cKO) and leave-one-mouse-out direction/stability summaries where identifiable. These are exploratory and do not remove the design confounding.

### TF activity and regulon inference — method choice pending

The proposed full branch uses motif-supported pySCENIC regulons (GRNBoost2, matching mouse cisTarget pruning and AUCell) plus a separately labeled CollecTRI/decoupler activity analysis. The smaller staged alternative runs TF activity plus hdWGCNA first. Do not substitute curated activity for de novo inferred regulons. Full measured gene background is required for both; the DEG-only selection policy applies specifically to communication endpoints. Keep WT/cKO regulon definitions fixed for comparison. Audit actual Satb2 motif/prior coverage before claiming a Satb2 regulon. No ATAC data are present; enhancer-binding or SCENIC+ claims are unsupported.

### Communication — reuse exact completed results

Keep only ligand and receptor components that pass the matching original population's DESeq2 padj<0.05, requiring every complex subunit. Retain original native LIANA ranks/scores from the full scoring universe and original component-specific cKO/WT effect signs. Native lr_logfc is group-versus-rest, not a genotype interaction effect. Supply CA1→IN, IN→CA1 and both-direction views. C8 has no qualifying LR-resource overlap; this is not evidence of absent biological communication.

## Evidence contract and outputs

- Input hashes, exact membership rosters, all12 coverage, stable-ID/symbol mapping, package/source versions, selected network genes and metacell coverage.
- Soft-power diagnostics, modules including grey, hub/kME tables, continuous mouse-level eigengenes and effect/sensitivity tables; no nucleus-level pseudoreplication.
- Clearly distinguish TF RNA expression, curated TF activity, inferred motif-supported regulon activity and co-expression module membership.
- PDF/SVG/PNG figures, machine-readable complete tables, a compact portable offline regulatory explorer, and an additive link in the GitHub atlas only after results pass validation.
- Keep raw scores in downloads and mark display-only z-scores, colors, node/edge widths and layout edits. Artificial display coordinates cannot enter inference.

## Limits and implementation status

Exploratory modules and activities nominate hypotheses; they do not establish direct binding, receptor activation, synaptic direction, or a causal SATB2 mechanism. C8's mixed identity and small counts limit interpretation. No experimental validation is performed. New isolated dependencies are required, available disk space is limited/variable, and the remote registered compute host was unreachable during inspection. No scientific workflow has started at this plan's creation. Workflow readiness and execution approval are separate subsequent artifacts; the existing project Snakemake controller was already selected by the user.
