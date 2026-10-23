# Corrected DEG-only CA1–inhibitory communication projection

The previous default used a ligand-DEG OR receptor-DEG rule, which retained non-DEG partners and other inhibitory populations. This corrected dataset requires **every ligand AND receptor subunit to be a significant, unambiguous DEG at its own exact endpoint**. Only CA1 Deep, CA1 Superficial, inhibitory C1, C6 and C8 are eligible. Other inhibitory clusters and inhibitory total are excluded.

The corrected atlas uses only the saved DESeq2 `padj < 0.05` genes from the five requested endpoints. It projects completed native LIANA results; it does not refit expression, DESeq2 or LIANA, change expression normalization, or recalculate ranks using a restricted background. Native LIANA scores and ranks retain their original full-gene/background interpretation.

## Exact endpoint DEG counts

| Endpoint | Original DEGs | Up in cKO | Down in cKO | Unambiguous resource DEG genes |
|---|---:|---:|---:|---:|
| CA1_Deep | 742 | 360 | 382 | 106 |
| CA1_Superficial | 1825 | 928 | 897 | 199 |
| IN_C1 | 27 | 11 | 16 | 7 |
| IN_C6 | 22 | 4 | 18 | 5 |
| IN_C8 | 1 | 1 | 0 | 0 |

CA1 Deep uses `P0244_seurat_CA1_ProS_LAYER_PROGRAM`, label `Deep-associated expression program enriched`; CA1 Superficial uses `P0247_seurat_CA1_ProS_LAYER_PROGRAM`, label `Superficial-associated expression program enriched`. IN C1/C6/C8 use `P0140`, `P0145`, `P0147` Seurat inhibitory cluster endpoints, respectively. Every copied DEG file is byte-identical to the original. Superficial significant symbols Pakap and Gm35438 are ambiguous across stable IDs and excluded from primary interaction support; neither occurs in the LR resource.

## Strict population-specific interaction counts

| Native scoring condition | CA1→IN, all saved ranks | IN→CA1, all saved ranks | Combined | Combined at magnitude rank≤0.05 |
|---|---:|---:|---:|---:|
| Pooled | 6 | 10 | 16 | 8 |
| WT | 5 | 10 | 15 | 9 |
| cKO | 6 | 10 | 16 | 7 |

These are 47 saved records across the three conditions and 16 distinct population × LR keys. WT/cKO counts reflect native expression qualification and ranking; their difference is descriptive and is not evidence of differential interaction activity. At magnitude rank≤0.05, the directional counts are pooled 3/5, WT 3/6 and cKO 3/4 (CA1→IN / IN→CA1). The stricter replicated/context-specific criteria retain the two C6 Cntn3→CA1 Ptprg candidates in each condition.

### Why C8 has zero strict interactions

C8's sole DEG is **Calb1**, up in cKO (log2FC +1.93206760368337; padj 0.0154178430986269). Calb1 is absent from the mouseconsensus ligand/receptor resource. C8 remains an explicit eligible population with tested directions and zero retained edges. This is not an omission and does not establish an absence of biological communication. Its original mixed-identity review flag remains attached.

## Gene up/down distinction

Positive original gene DESeq2 log2(cKO/WT) means **up in cKO**; negative means **down in cKO**. Both endpoint effects and adjusted P values accompany each interaction. `down→up`, for example, means the sender's ligand DEG is down in cKO and the receiver's receptor DEG is up. It does not assign an up/down activity change to the entire interaction. For a multi-subunit complex with opposing component signs, report mixed subunit directions rather than inventing a single complex effect. Every retained complex here is a single-subunit symbol.

The strongest strict candidates are IN C6 **Cntn3→Ptprg**, directed to both CA1 programs. Cntn3 is down in C6 (log2FC −0.5383457365; padj 0.0037005639), whereas Ptprg is up in Deep (+0.3966628983; padj 0.0164890846) and Superficial (+0.4904369824; padj 2.9270547465×10⁻8). Both endpoints have significant gene effects and all six WT and six cKO mice support the expression gate. The opposite signs still do not establish net receptor activity.

All pooled strict candidates are listed below. Duplicate reverse Nrxn/Nlgn resource orientations do not represent independent physiological pathways.

| Source | Target | Ligand | Receptor | Gene effect pair | Ligand log2(cKO/WT) | Receptor log2(cKO/WT) | Native magnitude rank |
|---|---|---|---|---|---:|---:|---:|
| IN_C1 | CA1_Superficial | Nlgn1 | Nrxn1 | up→down | +0.41412 | -0.21320 | 1.04686e-05 |
| CA1_Superficial | IN_C1 | Nrxn1 | Nlgn1 | down→up | -0.21320 | +0.41412 | 1.04686e-05 |
| IN_C6 | CA1_Superficial | Cntn3 | Ptprg | down→up | -0.53835 | +0.49044 | 8.29423e-05 |
| IN_C6 | CA1_Deep | Cntn3 | Ptprg | down→up | -0.53835 | +0.39666 | 8.86983e-05 |
| CA1_Deep | IN_C6 | Slit1 | Robo1 | down→down | -0.30242 | -0.45823 | 0.0074407 |
| CA1_Superficial | IN_C1 | App | Grm7 | up→up | +0.10733 | +0.35512 | 0.00748914 |
| IN_C6 | CA1_Deep | Slit2 | Dcc | down→up | -1.49306 | +0.28712 | 0.0321868 |
| IN_C6 | CA1_Superficial | Slit2 | App | down→up | -1.49306 | +0.10733 | 0.0496269 |
| IN_C6 | CA1_Superficial | Slit2 | Robo1 | down→down | -1.49306 | -0.35684 | 0.0536709 |
| IN_C1 | CA1_Superficial | Fstl5 | Chl1 | up→up | +0.64933 | +0.15579 | 0.0979569 |
| CA1_Deep | IN_C1 | Ptn | Ptprz1 | down→up | -0.64643 | +0.59000 | 0.202279 |
| CA1_Superficial | IN_C1 | Ptn | Ptprz1 | down→up | -1.01569 | +0.59000 | 0.382691 |
| CA1_Superficial | IN_C1 | Penk | Grm7 | up→up | +0.80508 | +0.35512 | 0.636421 |
| IN_C1 | CA1_Superficial | Fstl5 | Cntn5 | up→down | +0.64933 | -0.94853 | 1 |
| IN_C6 | CA1_Deep | Slit2 | Plxna1 | down→down | -1.49306 | -0.55460 | 1 |
| IN_C6 | CA1_Superficial | Slit2 | Plxna1 | down→down | -1.49306 | -0.71521 | 1 |

## Statistical meanings and limits

- Gene DESeq2 `log2FoldChange` is the signed cKO/WT component effect from six independent mice per genotype.
- Native LIANA `lr_logfc` is a source/target group-versus-rest specificity score, **not** a genotype fold change.
- Native `magnitude_rank` and `specificity_rank` are context priority ranks, not genotype P values or FDRs. Smaller ranks are higher priority.
- The LIANA CellPhoneDB scorer estimates cell-label specificity inside a scoring condition; it is not a standalone CellPhoneDB differential-genotype test.
- No genotype interaction-activity P value was fitted. Transcript expression and curated priors do not establish protein localization, contacts, activation/inhibition, functional synapses or anatomical layers.
- Every original biological source file, the original v3.5 payload, and all prior plots remain preserved. Version 3.6 is an explicitly labeled corrected selection.

## Downloads and reproducibility

The condition/direction `rank_aggregate.tsv` files keep exactly the requested seven native score columns: `source, target, ligand_complex, receptor_complex, lr_logfc, lr_gmeans, magnitude_rank`. Separate signed extended tables add ligand and receptor gene log2FC, adjusted P, up/down and ordered effect pairs. `DESeq2/{population}/DEG_FDR05.txt` keeps the original seven gene columns and bytes; `UP_in_cKO.txt` and `DOWN_in_cKO.txt` are endpoint-matched subsets. Significant stable-ID audit rows are also supplied.

`DEG_ONLY_ENDPOINT_AUDIT.json` contains the independent source-hash and endpoint reconstruction checks. `MANIFEST.json` gives compact transport hashes, source hashes, sign counts, rules and download inventory. The read-only script `scripts/audit_deg_only_communication.py` reproduces this projection using completed native outputs; it performs no scientific fits.
