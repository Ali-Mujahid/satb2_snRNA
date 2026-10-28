# CA1–inhibitory communication, gene programs and TF target signatures

Analysis completed 26 October 2026; final audit and atlas publication prepared 28 October 2026. The original Scanpy/Seurat analyses, annotations, gene-level DESeq2 fits and strict DEG-only communication results are preserved.

The most consistent new findings are **higher pooled CA1 black-module scores and lower purple-module scores in cKO in both CA1 programs**, with superficial brown higher as a secondary lead. Curated TF analysis identifies 23 primary target-signature associations, but **none passes the matched sequencing-partition sensitivity FDR**. **SATB2 target-signature activity is not FDR-significant**, with or without its self-target link. Motif-supported pySCENIC regulons could not be inferred because matching mouse reference files were unavailable; the prepared environment is not a completed regulon analysis.

[Interactive regulatory explorer](https://ali-mujahid.github.io/satb2_snRNA/regulatory/) · [Native figure gallery](https://ali-mujahid.github.io/satb2_snRNA/regulatory/gallery.html) · [Tables, source and offline downloads](https://ali-mujahid.github.io/satb2_snRNA/regulatory/downloads.html) · [Strict DEG-only communication explorer](https://ali-mujahid.github.io/satb2_snRNA/communication/)

## Exact scope

All genotype inference uses the twelve independent mice, six WT and six cKO. Nuclei and overlapping metacells are not genotype replicates. These analyses cover five frozen Seurat endpoints:

| Display population | Original endpoint | Nuclei |
|---|---|---:|
| CA1 DL-associated expression program | P0244_seurat_CA1_ProS_LAYER_PROGRAM | 8,750 |
| CA1 SL-associated expression program | P0247_seurat_CA1_ProS_LAYER_PROGRAM | 11,706 |
| Inhibitory C1 | P0140_seurat_Inhibitory_CLUSTER | 1,718 |
| Inhibitory C6 | P0145_seurat_Inhibitory_CLUSTER | 501 |
| Inhibitory C8 | P0147_seurat_Inhibitory_CLUSTER | 468 |

The total is 23,143 selected nuclei. CA1 DL/SL are molecular expression programs; measured anatomical layer position is unavailable. C8 retains mixed molecular identity. Findings should not be generalized automatically to every CA1 or inhibitory nucleus in the broader atlas.

The communication restriction remains **original DESeq2 adjusted P < 0.05 for every ligand subunit in the exact source and every receptor subunit in the exact target**. In contrast, coexpression and TF target-signature inference require a full, genotype-blind eligible measured-RNA background, including non-DEGs. Their different backgrounds and statistics serve different questions; no communication eligibility was relaxed.

## How to read the three communication directions

The existing corrected analysis contains 16 pooled ordered population–interaction records: six CA1→inhibitory and ten inhibitory→CA1. WT contains 15 and cKO 16 qualifying saved records. Counts across these separate expression-based datasets are descriptive; they are not tests of a genotype change in signaling. The combined direction retains both ordered sets and does not make an interaction biologically bidirectional.

| Direction | Example signed lead | Defensible interpretation |
|---|---|---|
| CA1→inhibitory | Slit1 lower in DL → Robo1 lower in C6 | Both component RNAs decrease in their matching populations; a candidate for follow-up, not proven reduced signaling |
| Inhibitory→CA1 | Cntn3 lower in C6 → Ptprg higher in DL and SL | Opposing component changes; do not collapse them into one up/down communication label |
| Both directed sets | Preserve the two source→target orientations | Compare exact component DEGs, expression coverage and native ranks without treating reverse orientation as the same record |

C8's only original DEG, Calb1, has no qualifying mouseconsensus LR subunit overlap. No strict eligible C8 interaction does not mean C8 has no biological communication. Ptn–Ptprz1 and adhesion/synaptic pairs retain their context and identity caveats. The full signed component tables remain in the communication explorer and the native archive.

Native LIANA `lr_logfc` is a group-versus-rest statistic, **not** gene log2(cKO/WT). LIANA rank is a context-priority score, **not** a genotype FDR. A receptor increase plus ligand decrease is not evidence of compensation or of the net direction of signaling. The new module cross-reference organizes the same measured RNA; it does not supply independent communication validation.

## Detailed hdWGCNA findings

Signed coexpression networks were fitted separately within the selected CA1 and inhibitory lineages using hdWGCNA 0.4.12/WGCNA 1.74, full eligible RNA, genotype-blind coverage and prespecified topology gates. Native nucleus PC1 projections were oriented toward average module RNA, then summarized per mouse and exact population. Separate continuous-score linear models test cKO minus WT. These score effects are not DESeq2 gene fold changes or new DEG lists.

| Network | Tested genes | Modules | Assigned genes | Grey/unassigned genes |
|---|---:|---:|---:|---:|
| Primary CA1 pooled | 10,618 | 10 | 4,898 | 5,720 (53.87%) |
| Primary inhibitory pooled | 11,065 | 5 | 2,304 | 8,761 (79.18%) |
| CA1 consensus sensitivity | 10,618 | 8 | 4,176 | 6,442 (60.67%) |
| Inhibitory consensus sensitivity | 11,065 | 4 | 2,893 | 8,172 (73.85%) |

There are 1,200 CA1 and 359 inhibitory metacells. C1 contributes 259/359 inhibitory metacells, so discovery is weighted toward C1. The 13-cell C8/cKO933 stratum was below the metacell-discovery minimum but all its nuclei remain in twelve-mouse scoring. Grey means unassigned under these settings, not biological absence.

The primary global BH family contains 35 population–module tests: 12 pass FDR < 0.05, eight CA1 and four inhibitory. The matched sequencing-partition family has five significant associations, all CA1:

| Primary CA1 module / program | cKO − WT native eigengene difference | Primary FDR | Same-partition FDR |
|---|---:|---:|---:|
| black / DL | +2.806 | 3.70e-7 | 2.91e-5 |
| black / SL | +4.030 | 8.29e-11 | 9.18e-8 |
| purple / DL | −2.073 | 7.10e-6 | 1.16e-4 |
| purple / SL | −2.356 | 3.54e-8 | 2.86e-5 |
| brown / SL | +1.529 | 0.00110 | 0.0422 |

Each retains its direction in twelve leave-one-mouse-out fits, and its exploratory mouse-bootstrap interval excludes zero. Black contains hubs such as Fgf12, Calb1, Ncam2, Ntm and Dcc; purple includes Astn2, Cntnap5c, Gabrg3, Kcnip4 and Ptn; brown includes Pde10a, Ndst4, Ntrk2 and Sorcs3. Their hubs, original directional DEGs and mouse distributions provide coherent within-dataset hypotheses. Differences in unrelated module PC1 units should not be compared as biological effect magnitudes.

Inhibitory blue and green are higher in C1 and C6 in the primary model, but none passes same-partition FDR. C8 has no primary FDR-positive module. Consensus has five of 28 primary associations and none of 28 same-partition associations. Consensus uses the same mice and is neither independent replication nor formal preservation. Colors are arbitrary across fits: consensus black is not the same module as pooled black. Compare actual member genes.

GO BP enrichment has 993/72,805 tests below global FDR 0.05; original-DEG overlap has 40/105. These families include all declared tests across both lineages. Illustrative annotations include synaptic signaling in CA1 black, cell-junction assembly in purple, and synapse organization in inhibitory green. Black is strongly enriched for original CA1 up DEGs, and purple for down DEGs. Such agreement reuses the same RNA and is contextual consistency, not independent confirmation. No DEG-overlap test passes FDR in C6 or C8.

Primary CA1 green and inhibitory yellow have nine mitochondrial genes among their top ten hubs. Treat these as QC/metabolic/possible carryover programs until independently reviewed. Sparse nuclear RNA materially invokes the configured individual Pearson fallback for bicor kME; median nucleus connectivity is modest. The native R session retained only a “50 or more warnings” notice, so the exact warning types/count remain unknown. The source and audit retain these limitations.

## Curated TF target-signature results

Native decoupler 2.3.0 ULM uses the recorded mouse-mapped CollecTRI signed prior, original DESeq2 size factors, a full eligible RNA background and at least five eligible targets per source. Genotype association is a separate mouse-level continuous-score model. Per-mouse gene-set FDR and genotype-association FDR are distinct.

| Exact population | Eligible RNA genes | Tested sources | Primary FDR < 0.05 | Same-partition FDR < 0.05 |
|---|---:|---:|---:|---:|
| CA1 DL | 15,899 | 622 | 7 | 0 |
| CA1 SL | 16,969 | 638 | 13 | 0 |
| Inhibitory C1 | 12,749 | 563 | 1 | 0 |
| Inhibitory C6 | 8,980 | 466 | 2 | 0 |
| Inhibitory C8 | 9,203 | 479 | 0 | 0 |

The primary FDR family comprises all 2,768 tested population–source combinations. The three inhibitory candidates are Maf lower in C1, and Klf15/Rest lower in C6. All 23 primary associations retain direction in twelve leave-one-mouse-out fits and in the smaller partition subset, but no association passes the latter's separately corrected FDR. Reduced power in that subset matters; it does not prove the primary effects absent.

Thirteen of the 23 significant signatures have source-TF RNA below the genotype-blind expression threshold. In superficial CA1, Dlx4, Foxh1, Hoxa5, Prop1 and Tcf21 have zero source counts across all twelve mice. A nonzero target-pattern score can be computed without source RNA; it is not confirmation that the named TF protein is active in these neurons. Myrf's CA1 signature also warrants cell-identity/contamination review. No post-hoc source-RNA exclusion was made. Some prior sources, including modeled Dmtf1/Dmtf1l, have identical signed target sets and identical scores; those regulators cannot be distinguished by this assay. None of the 23 significant signatures belongs to an exact duplicate group. See the full [independent TF audit](20261028_TF_final_audit.md) for all candidates and source-expression cautions.

**SATB2:** full-prior target-signature FDR is 0.999 in DL and 0.381 in SL; all inhibitory populations are also nonsignificant. Removing the Satb2→Satb2 self-link gives DL effect −0.129/FDR 0.162 and SL −0.313/FDR 0.107 across the five-endpoint sensitivity family. The SL nominal P is 0.0213, which does not pass FDR. Do not describe these results as established loss or gain of SATB2 activity. Original gene-level Satb2 RNA is actually higher in cKO CA1 (DL log2FC +1.350; SL +1.095). Gene-level RNA does not establish which exons/transcripts remain or whether SATB2 protein is depleted; deletion design and functional loss require independent clarification.

## What to prioritize and how to use the atlas

Start with primary CA1 black, purple and superficial brown, inspect the twelve mouse points, then their actual hub genes, directional original DEGs and enrichment. Retain inhibitory blue/green and the 23 curated target signatures as exploratory leads with visible sensitivity limits. Keep individual LR gene directions: a module's direction does not apply to every member. For example, Cntn3 decreases in C6 even though its pooled blue-module score increases.

In the explorer choose a saved method, exact population and module/TF. The **Show SATB2 activity in all five populations** button opens the measured mouse-score panels. Switch native units/display z, per-mouse/genotype summaries, colors and numeric limits, point sizes, node layout, saved edge thresholds/widths, arrow sizes for directed priors, and arc appearance. Genotype display filters never refit the saved effects or alter the twelve-mouse z reference. NA remains missing; it is not zero. Network arrows show curated-prior orientation, not a newly inferred causal network. hdWGCNA edges remain undirected.

PDF/PNG reproduce the current plot, SVG retains vector geometry, and CSV/settings/plot-provenance JSON downloads preserve selection and source identity. The native gallery contains 58 figure groups with original PDF/PNG; SVGs and full native tables/state are in the complete archive. The standalone HTML embeds its rendering dependencies and precomputed payload; it runs without R, Python or a server. Native/source downloads remain separate release files. No per-nucleus TF activity was computed, so activity UMAP is explicitly unavailable rather than painting mouse averages onto nuclei.

The matched-partition analysis includes six WT and only three cKO mice; genotype partly confounds technical partition. These observational analyses do not establish direct binding, physical communication, anatomical layer location or a causal SATB2 mechanism. Full pySCENIC motif inference remains unavailable because both official hosts refused connections and no cache exists; no human, mm9 or mismatched database was substituted.

## Reproducibility and preserved analyses

[Executed source and configurations](../regulatory_20261026/) · [Independent hdWGCNA audit](20261028_HD_final_audit.md) · [Independent TF audit](20261028_TF_final_audit.md) · [Native hdWGCNA completion review](20261028_hdwgcna_completed.md) · [Native TF completion review](20261028_collectri_completed.md)

The original export failure is retained. Its recovery verified and copied 250 prepared/network files (812,883,283 bytes) without refitting; all six recovered enrichment TSVs exactly match the earlier computation. Original memberships, counts and gene-DE fits pass before/after hash checks. The release archive includes selected raw matrices, recoverable RDS states, native tables/figures, source/configuration and receipts; the original full H5AD and installed package libraries remain external hash-bound dependencies. The existing controller/base R software closure is not fully lockfile-backed, which limits exact cross-system reproducibility despite recorded versions and source hashes.
