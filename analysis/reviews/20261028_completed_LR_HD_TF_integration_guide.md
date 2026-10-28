# Reading the completed communication, module and TF evidence together

This addendum supplies **exact gene-level joins** for the three requested example families, rather than adding a new fitted analysis or changing the existing findings guide. It uses completed strict DEG-only LIANA results, pooled hdWGCNA membership and mouse-score associations, original endpoint-specific DESeq2 results, and saved CollecTRI target links. No new communication effect, regulon, causal pathway or significance test is inferred here. Results from the later motif branch are deliberately excluded.

## Keep three different directions separate

“CA1→inhibitory” means that the ligand's RNA comes from the CA1 population and the receptor's RNA from the inhibitory population. The gene up/down signs refer to **cKO versus WT within those exact populations**. A module effect is a separate difference in a multi-gene PC1 score. None determines the others, and a LIANA native `lr_logfc` is group-versus-rest specificity, not genotype log2FC.

The five ordered pairs below are **selected interpretation examples**, not an exhaustive interaction list or a new significance filter. They occur in pooled, WT and cKO native scoring datasets, giving 15 example records. The complete strict communication datasets remain unchanged: **16 pooled, 15 WT and 16 cKO records**. Each selected pair satisfies the unchanged rule that both exact-endpoint component genes are unambiguous original DESeq2 DEGs. All six mice of each genotype meet their saved expression-support gate. This repeated expression qualification is useful coverage information; it is not a test that communication differs between genotypes.

| Ordered RNA candidate | Ligand log2(cKO/WT), original FDR | Receptor log2(cKO/WT), original FDR | Pooled native magnitude rank |
|---|---|---|---:|
| C6 Cntn3 → CA1 DL Ptprg | −0.53835; 0.0037006 | +0.39666; 0.016489 | 8.86983e-5 |
| C6 Cntn3 → CA1 SL Ptprg | −0.53835; 0.0037006 | +0.49044; 2.92705e-8 | 8.29423e-5 |
| CA1 DL Slit1 → C6 Robo1 | −0.30242; 0.019179 | −0.45823; 3.38377e-5 | 0.00744070 |
| CA1 DL Ptn → C1 Ptprz1 | −0.64643; 5.28170e-5 | +0.59000; 0.0038241 | 0.202279 |
| CA1 SL Ptn → C1 Ptprz1 | −1.01569; 4.67338e-6 | +0.59000; 0.0038241 | 0.382691 |

These ranks retain their original full-background meaning and were not recomputed after DEG selection. Smaller ranks indicate higher native priority, not smaller genotype p values. In the combined view, the two ordered directions remain separate records; it does not establish that each pair operates in both biological directions.

## What the module context actually adds

Module identity below includes the lineage and **pooled** network. A CA1 blue module and inhibitory blue module are different gene sets. Primary module tests use six WT/six cKO mice; same-partition sensitivity uses six WT/three cKO, with its own global FDR family.

| Exact component gene | Native pooled module | Gene's nuclear kME / within-module rank | Module cKO−WT effect | Primary module FDR | Same-partition module FDR | What is supported |
|---|---|---|---:|---:|---:|---|
| C6 Cntn3 | IN blue | 0.51369 / 32 | +1.20295 | 0.01214 | 0.15885 | Gene decreases while the module score increases; module context cannot be used as an interaction sign |
| DL Ptprg | CA1 turquoise | 0.43558 / 75 | +0.95713 | 0.29253 | 0.48109 | Receptor gene increase is not accompanied by a significant module effect |
| SL Ptprg | CA1 turquoise | 0.43558 / 75 | +1.27978 | 0.08527 | 0.27624 | Same-sign descriptive module tendency; it does not pass FDR |
| DL Slit1 | CA1 blue | 0.27019 / 40 | −0.44336 | 0.62047 | 0.76124 | Ligand and module directions agree descriptively, without a significant module association |
| C6 Robo1 | IN turquoise | 0.53748 / 44 | +0.14913 | 0.77317 | 0.98115 | Receptor gene decreases; near-zero positive module effect is nonsignificant |
| DL Ptn | CA1 purple | 0.49064 / 7 | −2.07299 | 7.10227e-6 | 1.15801e-4 | Coherent ligand-RNA and broader program decrease in this dataset |
| SL Ptn | CA1 purple | 0.49064 / 7 | −2.35614 | 3.54419e-8 | 2.86037e-5 | Same coherent gene/program decrease, retaining partition sensitivity support |
| C1 Ptprz1 | IN blue | 0.51140 / 34 | +0.73001 | 0.01236 | 0.16552 | Receptor RNA and primary module score increase, without partition-FDR confirmation |

All listed module effects retain direction in 12/12 native leave-one-mouse-out fits except C6 turquoise, which retains 11/12. This is a sensitivity check on the same experiment, not independent replication. A positive within-nucleus kME describes cross-nucleus correlation with the module; it does not require the gene's genotype contrast to match the module-score contrast. PC1 orientation and units also differ between modules, so their numerical effect magnitudes are not calibrated biological effect sizes across programs.

For presentation, the **Cntn3→Ptprg** examples are opposing component-RNA changes with strong saved native rank and full mouse expression support. Do not call the receptor increase compensatory or assign one net signaling sign. The **Slit1→Robo1** example has both RNAs lower, but its modules add no FDR-supported confirmation; lower signaling is still a hypothesis. **Ptn's decrease** has the clearest matching program context among these example ligands, but that does not elevate its lower-priority LR ranks to a significant genotype communication result.

## Why Ptn–Ptprz1 keeps an identity and functional caveat

PTN–PTPRZ has experimentally demonstrated oligodendrocyte-lineage biology, including PTN-mediated inhibition of PTPRZ and effects on precursor differentiation. This makes cell-identity context and ligand function important before calling a C1 RNA association neuron-to-interneuron signaling. It does **not** prove that the C1 signal is contamination, nor does it justify relabeling the cluster from these data alone. [Primary mouse/OPC study](https://pmc.ncbi.nlm.nih.gov/articles/PMC6605312/).

The exact receptor cell identity, protein/isoform localization and physical contact are not established by a population RNA mean. Increasing receptor RNA and decreasing ligand RNA cannot determine net receptor activity, especially when ligand binding need not activate the receptor's phosphatase function. Preserve C1's original annotation and the caveat rather than treating a single receptor as a new cell label.

## Genuine TF joins, including the negative findings

The traceable prior-overlap table contains 35 genuine saved CollecTRI source→target links to these eight exact LR component rows. It includes every source at these targets, not only statistically interesting ones. A further eight prior links concern the separately labeled CA1 Satb2 gene-context rows. No edge or target was invented, and no new rank or fitted effect was calculated.

**None of these eight LR components is a Satb2 target in its saved used CollecTRI network.** This resource-specific absence cannot exclude an uncurated or indirect mechanism, but it provides no supported direct SATB2→LR link to draw in a causal diagram.

Only one primary-FDR-positive TF signature joins to these specific components: a saved **Hoxa5→Ptn** positive prior in superficial CA1. Its target-signature effect is −0.23651 (global FDR 0.02841), but same-partition FDR is 0.10144 and **Hoxa5 RNA is zero in all 12 mice**. The overlap is real as a curated-prior record; it does not establish active HOXA5 regulation in these neurons and should not be promoted as the mechanistic explanation for Ptn loss.

Satb2 itself belongs to CA1 pooled black (nuclear kME 0.29471, rank 31). Its original gene RNA is higher in cKO in DL (+1.34993, FDR 2.08256e-48) and SL (+1.09451, FDR 8.63648e-19), alongside higher black-module scores with partition support. However, the full-prior SATB2 target-signature associations are nonsignificant, and the prior includes a self-target edge. These observations do not demonstrate functional SATB2 gain, successful protein deletion, direct binding, or a causal link to any LR pair. Gene-level counts alone do not identify retained exons/transcripts or protein state.

## Files and evidence limits

- `20261028_LR_HD_TF_components.tsv`: 10 exact endpoint–gene rows, original DE values, native module membership/kME/rank, primary and partition effects/FDR, leave-one-mouse-out signs, and source paths.
- `20261028_LR_HD_TF_directional_pairs.tsv`: the 15 unchanged condition-specific ordered LR records and component gene directions.
- `20261028_LR_HD_TF_saved_prior_overlap.tsv`: 43 actual used prior links, source-RNA eligibility/absence, and saved mouse-score/partition statistics.
- `20261028_LR_HD_TF_sources.json`: source and output byte-count/SHA256 receipts, exact scope and explicit noncausal interpretation boundaries.

The same RNA experiment contributes gene DE, modules, target signatures and expression-qualified LR records. Their agreement is **within-dataset convergence**, not independent validation. No gene-level partition-adjusted DE test or genotype interaction-activity test was fitted for this addendum. Original model confounding remains, C8's lack of strict resource-qualified DEG interactions is not biological silence, and future motif outputs require separate review before being added as another evidence layer.
