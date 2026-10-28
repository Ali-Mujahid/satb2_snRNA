Takeaway: Completed: Four preserved hdWGCNA networks pass audit; CA1 associations are stronger than inhibitory leads.

This read-only audit examines the completed recovery of coexpression analysis in five frozen Seurat populations from six WT and six Satb2 cKO mice. Both pooled networks and both conditional consensus networks have complete native checkpoints, finite independent-mouse scores, and a successful final package. The recovery copied the finished scientific fits without refitting them, and its six enrichment tables exactly reproduce the tables computed before the earlier SVG export failure. The most consistent observational genotype associations are increased pooled CA1 black, decreased CA1 purple, and increased superficial-program CA1 brown scores. Inhibitory results remain exploratory because no same-partition inhibitory association passes FDR 0.05, the discovery network is dominated by C1, and C8 retains mixed molecular identity. Coexpression, enrichment, and overlap with earlier DEGs organize hypotheses; they do not establish direct SATB2 regulation, physical communication, or causal pathway activation.

## 1. Scientific scope and independent units

The exact selected populations are CA1 deep-associated expression program, 8,750 nuclei; CA1 superficial-associated expression program, 11,706; inhibitory C1, 1,718; C6, 501; and C8, 468. Thus CA1 contains 20,456 nuclei and the three inhibitory populations contain 2,687, totaling 23,143. The CA1 labels describe molecular programs rather than measured anatomical positions. These findings apply to these five frozen endpoints and should not be generalized to every CA1 or inhibitory nucleus in the broader atlas.

The fitted networks use the eligible measured RNA background, rather than only DEGs or a DEG-selected feature universe. The separate ligand–receptor display retains its strict endpoint-specific DEG-only rule. Earlier original DESeq2 results were reused for module overlap and gene direction; no gene-level differential-expression model was refitted. Module eigengene association is a continuous-score mouse-level linear-model analysis, not DESeq2 and not a new DEG list. Nuclei and overlapping metacells are not independent genotype replicates.

## 2. Lifecycle, recovery, and provenance

Completed recovery: `/Users/mujahid357/Documents/codex_new/SATB2_REGULATORY_20261026/runs/20261026_hdwgcna_v1_0_2_export_recovery`.

Preserved source fits: `/Users/mujahid357/Documents/codex_new/SATB2_REGULATORY_20261026/runs/20261026_hdwgcna_v1_0_1`.

`results/FINAL_VALIDATION.json` has schema `satb2-regulatory-final-validation/v1` and status `COMPLETE`. It states that all four native branches are `NETWORK_COMPLETE`, original inputs/analyses are unchanged, preparation SHA-256 and native checkpoint MD5 seals are verified, no scientific fits occur in packaging, and completed fits were reused without refitting. Its TF-activity and motif-regulon flags are false because this is the hdWGCNA component, not the separate TF workflow.

I independently read all 250 records in `REUSED_NATIVE_RECEIPT.json` and recomputed SHA-256 and size for both the source file and copied recovery file. All 500 comparisons matched: 812,883,283 bytes of prepared/network outputs were copied exactly. I also compared every one of the six old and recovered enrichment TSVs byte-for-byte. All matched; only downstream graphics export and packaging were recovered. The recovery configuration SHA-256 is `7542a4ed9a0fe016232902685bb79a761f6613d8c60fdcdd95c00408bf168c61`; the original native fit configuration SHA-256 remains `ad59dd5cc71c7321ae49c52769ba7de80daaf2d005721b9dc8acc0ddc1f789bb`.

The four mouse-score tables contain 240, 180, 192, and 144 rows for CA1 pooled, IN pooled, CA1 consensus, and IN consensus respectively: 756 distinct endpoint–module–mouse keys in total. Every value is finite; every endpoint–module combination has all twelve independent mice. The recovered GO figure was visually inspected through its PNG and has readable axes, labels, legends, and declared enrichment background/FDR scope. Its native PDF parses as one page (24,065 bytes). The completed component contains 39 PDF files.

## 3. Networks, genotype associations, and robustness

| Network | Genes tested | Assigned genes | Non-grey modules | Grey genes | Metacells | Independent-mouse score rows |
|---|---:|---:|---:|---:|---:|---:|
| CA1 pooled | 10,618 | 4,898 | 10 | 5,720 (53.87%) | 1,200 | 240 |
| IN pooled | 11,065 | 2,304 | 5 | 8,761 (79.18%) | 359 | 180 |
| CA1 consensus | 10,618 | 4,176 | 8 | 6,442 (60.67%) | 1,200 | 192 |
| IN consensus | 11,065 | 2,893 | 4 | 8,172 (73.85%) | 359 | 144 |

The 1,200 CA1 metacells are balanced at 600 per exact program and 600 per genotype. Inhibitory discovery contains 259 C1, 50 C6, and 50 C8 metacells; C1 therefore contributes 72.14% of inhibitory metacells. The C8/cKO_933 stratum has thirteen nuclei, below the discovery minimum of twenty; it contributes no metacells but all thirteen nuclei remain scored in the twelve-mouse results. Inhibitory metacells comprise 174 WT and 185 cKO. Metacell overlap means their number is descriptive coverage, not a biological replicate count.

Primary signed-network soft powers passed the prespecified topology gate at CA1 power 7 and IN power 6. Consensus powers are CA1 WT 7/cKO 6 and IN WT 6/cKO 6. Grey means unassigned under these settings, not biological absence. The high grey fractions are especially material for inhibitory interpretation.

| Native global BH family | Tests | FDR <0.05 | CA1 / inhibitory hits | Independent-mouse cohort |
|---|---:|---:|---:|---|
| Pooled module × endpoint | 35 | 12 | 8 / 4 | 6 WT + 6 cKO |
| Pooled same-partition | 35 | 5 | 5 / 0 | 6 WT + 3 cKO |
| Consensus module × endpoint | 28 | 5 | 2 / 3 | 6 WT + 6 cKO |
| Consensus same-partition | 28 | 0 | 0 / 0 | 6 WT + 3 cKO |

The strongest pooled associations that also pass the prespecified same-partition sensitivity are:

| Exact endpoint / pooled module | cKO − WT eigengene difference | Primary FDR | Same-partition difference | Same-partition FDR |
|---|---:|---:|---:|---:|
| CA1 DL / black | +2.80572 | 3.699e-7 | +2.80923 | 2.913e-5 |
| CA1 SL / black | +4.03048 | 8.290e-11 | +4.00869 | 9.178e-8 |
| CA1 DL / purple | −2.07299 | 7.102e-6 | −2.03506 | 1.158e-4 |
| CA1 SL / purple | −2.35614 | 3.544e-8 | −2.31976 | 2.860e-5 |
| CA1 SL / brown | +1.52925 | 0.001100 | +1.25337 | 0.04223 |

Each of these five associations retains its sign in all twelve leave-one-mouse-out fits, and its 200-iteration mouse-bootstrap interval excludes zero. This supports prioritization within the existing data, while the module definitions and genotype tests remain data-derived and observational. Raw PC1 differences have module-specific units and cannot be ranked biologically by magnitude across unrelated modules.

Four inhibitory pooled primary associations pass global FDR: blue in C1 (+0.73001, FDR 0.01236), blue in C6 (+1.20295, 0.01214), green in C1 (+0.47201, 0.02896), and green in C6 (+0.75390, 0.01095). All retain positive sign across twelve leave-one-mouse-out fits and positive bootstrap intervals. However their same-partition FDR values are 0.1655, 0.1589, 0.2160, and 0.1589 respectively. Retain these as candidate synaptic/excitability programs; they are not validated inhibitory genotype mechanisms. Losing significance in a cohort with only three cKO mice does not itself disprove the primary estimate.

Consensus primary FDR hits are CA1 SL blue (+1.25939, 0.02206), CA1 SL pink (+0.76001, 0.02621), IN C6 brown (+1.13395, 0.02206), and IN turquoise in C1 (+0.88802) and C6 (+1.25648), both FDR 0.04619. None passes its own consensus same-partition global family. Consensus IN C6 brown is closest (same-partition FDR 0.06457). C8 has no FDR-positive primary module association. Its near-zero pooled turquoise estimate changes sign in half the leave-one-mouse-out fits (6/12 retained); green retains sign in 8/12, further limiting C8-specific claims.

Fixed-module connectivity stability is high (minimum leave-one-mouse-out kME correlations: CA1 pooled 0.95091, IN pooled 0.98276, CA1 consensus 0.95614, IN consensus 0.98801; minimum top-ten hub retention 0.8 in every network). These are fixed-module, fixed-projection diagnostics, not rediscovery after omitting a mouse and not formal WGCNA module preservation.

## 4. Completed GO and original-DEG enrichment

The six recovered TSVs are unchanged from the earlier computation. The GO family contains 72,805 tests across both lineages and all fifteen pooled non-grey modules; 993 pass global FDR <0.05 (511 CA1, 482 IN). Fourteen modules have at least one GO hit; CA1 pooled blue has none at this threshold. GO background is the fitted lineage network intersect mouse BP-annotated genes, with term sizes 10–2,000; counts include zero-overlap tests. GO annotations are GO.db 3.23.1/org.Mm.eg.db 3.23.0. These GO results describe pooled definitions only and must not be transferred to consensus modules by color name.

| Pooled module | Illustrative native enriched BP term | Overlap / BP-annotated module genes | Global FDR |
|---|---|---:|---:|
| CA1 black | synaptic signaling | 24 / 130 | 0.001179 |
| CA1 purple | cell junction assembly | 13 / 96 | 0.01574 |
| CA1 brown | tissue development | 107 / 517 | 5.646e-9 |
| CA1 yellow | chemical synaptic transmission | 67 / 457 | 7.006e-8 |
| IN blue | trans-synaptic signaling | 109 / 766 | 2.925e-12 |
| IN green | synapse organization | 16 / 78 | 0.002341 |
| IN turquoise | monoatomic ion transport | 114 / 877 | 1.830e-8 |
| IN yellow | proton transmembrane transport | 15 / 113 | 1.737e-10 |

The original-DEG overlap family contains 105 tests and 40 FDR-positive rows: 16 CA1 DL, 17 CA1 SL, seven IN C1, none C6 or C8. Tests include endpoint × pooled module × all/up/down direction, globally corrected together. Their universe is fitted network genes with finite original endpoint-specific adjusted P and unambiguous original symbol/stable-ID mapping. Eligible universes differ: DL 10,600; SL 10,605; C1 9,539; C6 2,190; C8 9,201. Thus overlap counts are restricted subsets of the original full DEG lists, rather than replacement DEG totals. The ambiguous-ID and duplicate-symbol exclusions are recorded separately and are not additive attrition totals.

The clearest directional convergence is:

- CA1 black: 78/152 module genes are original DL-up DEGs (330 up DEGs in eligible background; odds ratio 42.52; FDR 6.058e-77), and 119/152 are SL-up DEGs (778 background up DEGs; OR 53.52; FDR 4.730e-105).
- CA1 purple: 58/109 are DL-down DEGs (333 background down; OR 42.13; FDR 1.123e-57); 72/109 are SL-down DEGs (704 background down; OR 30.34; FDR 1.716e-57).
- CA1 brown: 99/575 are SL-up DEGs (OR 2.864; FDR 2.839e-15), consistent in direction with its mouse-level eigengene association.
- IN C1 blue: six original up DEGs among 678 module genes eligible for this endpoint (eleven up DEGs in background; OR 15.80; FDR 0.0001548). The total overlap is eight of 27 original C1 DEGs, FDR 0.001311.
- IN C1 yellow overlaps three up DEGs (FDR 0.0009766), but this module has mitochondrial-heavy hubs and should not be promoted as a specific regulatory mechanism.

IN C6 blue contains seven of its 22 eligible DEGs but its all-direction overlap FDR is 0.1436. C8's sole original DEG, Calb1, is in inhibitory pooled blue; the all-direction overlap FDR is 0.1873. Small DEG lists and endpoint-specific eligibility limit these tests; nonsignificant overlap does not establish absence of a biological program.

As arithmetic QA only, recomputing BH from exported P values reproduces native GO and DEG-overlap adjusted P values to maximum absolute differences 4.11e-15 and 1.17e-15 respectively. No model, enrichment, or P value was refitted by this audit. Enrichment and eigengene differences arise from the same data; agreement is contextual consistency rather than independent replication. Nested GO terms also share many genes and are not 993 independent mechanisms.

## 5. Signed ligand–receptor and SATB2 context

`STRICT_LR_AND_SATB2_MODULE_LINKS.tsv` has 26 endpoint–gene rows: 21 strict-LR rows spanning seventeen unique gene symbols, all assigned to pooled non-grey modules, plus five Satb2 rows. This is a descriptive cross-reference to the existing strict-DEG LR analysis. It does not fit new interactions, relax eligibility, or convert module edges into ligand–receptor edges.

Useful signed examples, using original DESeq2 log2FC cKO/WT:

| Exact endpoint / component | Pooled module | Original log2FC | Original FDR | Interpretation |
|---|---|---:|---:|---|
| IN C6 Cntn3 | IN blue | −0.53835 | 0.003701 | Existing inhibitory→CA1 ligand candidate decreases |
| CA1 DL Ptprg | CA1 turquoise | +0.39666 | 0.01649 | Candidate receiver RNA increases |
| CA1 SL Ptprg | CA1 turquoise | +0.49044 | 2.927e-8 | Candidate receiver RNA increases |
| CA1 DL Slit1 | CA1 blue | −0.30242 | 0.01918 | Existing CA1→inhibitory ligand candidate decreases |
| IN C6 Robo1 | IN turquoise | −0.45823 | 3.384e-5 | Candidate receiver RNA decreases |
| CA1 DL Ptn | CA1 purple | −0.64643 | 5.282e-5 | Component sits in a down-associated CA1 module |
| CA1 SL Ptn | CA1 purple | −1.01569 | 4.673e-6 | Component sits in a down-associated CA1 module |
| IN C1 Ptprz1 | IN blue | +0.59000 | 0.003824 | Receiver RNA direction opposes Ptn direction |

For example, Cntn3 decreases within C6 even though the blue module's overall C6 score increases. Module direction must never be assigned to every individual member. A ligand decrease plus receptor increase is a mixed-direction transcript observation, not proof of compensatory signaling or increased/decreased net communication. Native LIANA `lr_logfc` is its original group-versus-rest score and is not this DESeq2 genotype log2FC. Coexpression edge weight and module color are also different quantities.

Satb2 maps to CA1 pooled black and original gene-level RNA is **higher** in cKO: DL log2FC +1.34993, FDR 2.083e-48; SL +1.09451, FDR 8.636e-19. This unchanged original result requires separate exon/targeting/protein validation before interpreting the experiment as complete loss of SATB2 function. Gene-level RNA alone cannot identify residual transcripts, affected exons, mosaicism, targeting behavior, or protein depletion. Satb2 is not network-eligible in the selected inhibitory lineage; C1's original adjusted P is unavailable and C6/C8 have no corresponding fit row in this cross-reference. These unavailable entries remain missing, not zero activity.

## 6. Limits that must remain visible

Genotype is partially confounded with sequencing partition. All six WT and cKO_1139/1160/1165 share the lane6–7 partition; cKO_1138/991 are in lane8 and cKO_933 is unverified. The prespecified sensitivity keeps six WT and three cKO, so its reduced power and remaining biological heterogeneity matter. Passing it supports an observational lead; it does not prove technical independence or causal genotype effects.

Pooled CA1 green contains thirteen mitochondrial genes among 333 genes, and nine of its top-ten hubs are mitochondrial; IN yellow contains twelve among 122, likewise nine of its top-ten hubs. CA1 green's strongest enrichment is cytoplasmic translation (71/307, FDR 2.787e-50), alongside ribosomal/synaptic-translation annotations; IN yellow is enriched in proton/energy-related processes. Consider these QC/metabolic/possible cytoplasmic-carryover programs until independently examined. Their low GO FDR is not evidence of a SATB2-specific pathway. Do not retrospectively remove or relabel them to improve the story.

Sparse nuclear kME uses signed bicor with the configured individual Pearson fallback. The prior read-only detection audit showed that at least 67.99% of assigned CA1 pooled genes and 61.81% of assigned IN pooled genes have >50% zero values and therefore zero gene-vector MAD; corresponding consensus proportions are 68.37% and 66.64%. This guarantees material use of the configured fallback for nuclear hub connectivity, while it does not quantify fallback in aggregated metacell discovery. The original fit log records “50 or more warnings” without retaining their individual list. Their exact types/count cannot be reconstructed after session exit. Do not call every warning a zero-MAD warning or all reported correlations purely robust bicor. Detailed reasoning remains in `HD_NATIVE_WARNING_AND_CONSENSUS_AUDIT_20261026.md`.

Consensus uses the same mice and source metacells and is not independent replication or formal preservation. Its module colors denote different gene sets: for example IN consensus turquoise overlaps pooled blue in 735 genes, while IN consensus blue overlaps pooled turquoise in 697. CA1 consensus black predominantly overlaps pooled yellow (69/79), so it is not the primary CA1 black/Satb2 module. Module correspondence requires the exported overlap table rather than shared names. Native displayed hub networks use top-twenty-hub metacell coexpression edges, not a saved full TOM, not directed regulation, and not direct binding.

C8 remains a mixed molecular population, has only one original DEG, has no FDR-positive primary module association, and contributes one small thirteen-nucleus mouse stratum absent from discovery. It should remain visible with a review flag rather than be forced into a confident cell-type/mechanistic label. Activity UMAP is not part of this component; no per-mouse mean was painted onto nuclei as though it were a nucleus measurement.

## 7. Actionable interpretation and publication recommendation

Prioritize the pooled CA1 black-up and purple-down programs, with superficial brown-up as a secondary CA1 lead. Inspect their native mouse distributions, directional original DEGs, and exact hub members together. Black hubs include Fgf12, Calb1, Ncam2, Ntm, and Dcc; purple hubs include Astn2, Cntnap5c, Gabrg3, Kcnip4, and Ptn; brown hubs include Pde10a, Ndst4, Ntrk2, and Sorcs3. These are coherent gene-program hypotheses, not proof of physical CA1 layers or direct SATB2 targets. The unexpected original Satb2 RNA increase deserves independent experimental/transcript-structure clarification before a mechanistic conclusion.

Retain inhibitory blue/green differences as candidates with an explicit same-partition limitation; prioritize exact C1/C6 members rather than using a whole-lineage score to infer every inhibitory subtype. Use Cntn3–Ptprg and Slit1–Robo1 as previously defined signed LR leads, with module context as supporting organization only. Keep mitochondrial-heavy, grey/unassigned, and mixed C8 outcomes explicit. Publication is supported for this completed component when the explorer preserves native scores, feature definitions, endpoint coverage, missing values, separate FDR families, source receipts, and the limitations above. No extra hdWGCNA fitting is required to publish the completed outputs.

Primary evidence paths, relative to the completed recovery `results/` folder:

- `FINAL_VALIDATION.json`, `REUSED_NATIVE_RECEIPT.json`.
- `hdwgcna/{pooled,consensus}_module_genotype_{associations,same_partition}.tsv`.
- `hdwgcna/{CA1,IN}/{pooled,consensus}/module_membership.tsv`, `module_eigengenes_mouse.tsv`, `module_genotype_leave_one_mouse_out.tsv`, `module_genotype_mouse_bootstrap.tsv`, `fixed_module_connectivity_stability.tsv`, `completed.json`.
- `enrichment/GO_BP_ALL_TESTS.tsv`, `GO_BP_FDR05.tsv`, `GO_ANNOTATION_COVERAGE.tsv`, `ORIGINAL_DEG_MODULE_ENRICHMENT.tsv`, `ORIGINAL_DEG_ELIGIBILITY_COVERAGE.tsv`, `STRICT_LR_AND_SATB2_MODULE_LINKS.tsv`, and `GO_BP_dotplot.pdf`.

The audit read native outputs and performed checksums, counts, finite/unique-key checks, tabulation, and BH arithmetic QA only. It did not edit native data or source, load/modify Seurat RDS objects, refit networks/models/enrichment, change labels, or overwrite prior analyses.

### Byte-identical enrichment receipts

| Recovered TSV | SHA-256 (equal to earlier failed-export run) |
|---|---|
| GO_ANNOTATION_COVERAGE.tsv | `ebd189d9b9ac713ac0cb2ee2b01f9f50d2e8fab6c9a826bad44e863d942c4f06` |
| GO_BP_ALL_TESTS.tsv | `9b59e3cd13e4fc32d398ca5ed1d3029966b0159932c2de256c0831ecff1f18c9` |
| GO_BP_FDR05.tsv | `487b6070d8e1555f69f6788b78474a2303a8ccabbbfb7428e25a65cc87d4aa62` |
| ORIGINAL_DEG_ELIGIBILITY_COVERAGE.tsv | `a7eaa5a568e502e9ab9007c17c3538d315dad0a58ed5d94bb295758e1d74574c` |
| ORIGINAL_DEG_MODULE_ENRICHMENT.tsv | `89763bac27b5a20e5f8e76ef40f321ff6dd521249473bfdd5c4c25bd3f50556c` |
| STRICT_LR_AND_SATB2_MODULE_LINKS.tsv | `a904c24a5027f6f6a0373b1e0566725082b495844fdfc5b4dc040a9b19b9892b` |
