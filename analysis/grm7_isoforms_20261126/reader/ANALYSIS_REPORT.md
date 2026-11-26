Takeaway: Grm7 RNA is detected in all 12 mice, but full-length isoform abundance is not reliably resolved.

This analysis asked whether Grm7 expression in P60 Satb2 WT and cKO mice could be resolved into particular RNA isoforms in superficial and deep CA1 and inhibitory neurons. The corrected analysis completed for all twelve mice, using the frozen Seurat nucleus assignments and the mm10-2020-A annotation named in the original Cell Ranger configurations. Five transcripts are annotated, including two complete A-like and B-like receptor models, two incomplete coding models and one nonsense-mediated-decay model. Abundant gene-level RNA and transcript-associated exon segments were detected, but uniquely assignable splice-supported counts were very small. A B-type terminal-exon entry event was observed in ten superficial, eight deep and five inhibitory mouse samples, whereas the A-specific terminal junction had only two, zero and two supporting UMIs across those respective populations. These results support Grm7 transcription and some alternative-splicing features, but do not establish accurate full-length isoform abundances, isoform switching after SATB2 loss, or changes in presynaptic protein function.

## Scientific context and question
The study contains six WT and six Satb2 cKO mice at P60. Original Seurat superficial CA1 membership is C2+C3 (10,867 nuclei); deep CA1 is C1+C5 (8,750 nuclei). C8 is not included in superficial CA1. The complete CA1/ProS set has 27,703 nuclei; the original twenty inhibitory clusters have 8,911 nuclei. These sets were kept unchanged. Inhibitory C1 and parent C6 are separately viewable, without introducing new subtype assignments. Mouse, rather than nucleus or read, is the biological replicate.

## Question, lifecycle, and endpoint
The corrected workflow completed with exit code zero. Registry run: bed9be77-5caf-40be-8633-50eaf1162396. Native run root: /Users/mujahid357/Documents/codex_new/SATB2_GRM7_ISOFORMS_20261126/run_v2. The first run remains preserved and is scientifically superseded because its transcript-compatibility check mishandled coordinate types after JSON transfer; it must not be used to infer absent isoforms. The corrected version passed coordinate-boundary, serialization-roundtrip and complete synthetic SAM/UMI counter tests. Independent review confirmed count conservation, normalization arithmetic and unchanged gene/exon/junction support across both runs.

## Key findings

### Annotation in the actual mapping reference

| Transcript | Stable transcript ID | GTF classification | Receptor interpretation |
|---|---|---|---|
| Grm7-201 | ENSMUST00000071076 | Protein coding, 11 exons | Complete 922-aa B-like C-terminal model; GTF TSL5/not_organism_supported qualification retained |
| Grm7-202 | ENSMUST00000172951 | Protein coding, 10 exons | Complete 915-aa A-like C-terminal model |
| Grm7-203 | ENSMUST00000173001 | Protein coding, 2 exons | Incomplete coding start; not a complete receptor model |
| Grm7-204 | ENSMUST00000174018 | Nonsense-mediated decay, 10 exons | Predicted NMD RNA, not a normal full-length receptor |
| Grm7-205 | ENSMUST00000174310 | Protein coding, 3 exons | Incomplete coding start; contains a B-like tail but not a third complete receptor |

The A-like tail ends in SYNNLVI; the B-like tail ends in TIPPTV. Assignments were derived from the reference CDS sequences with the GTF reading phases respected. Newer databases may annotate additional transcripts; their counts cannot be substituted into this historical mapping reference without remapping and compatibility checks.

### Gene-level RNA support

| Population | Nuclei across 12 mice | BAM-derived Grm7 UMIs | Mice with gene support |
|---|---:|---:|---:|
| Superficial CA1, C2+C3 | 10,867 | 273,028 | 12/12 |
| Deep CA1, C1+C5 | 8,750 | 192,940 | 12/12 |
| Inhibitory total | 8,911 | 79,697 | 12/12 |
| Inhibitory C1 | 1,718 | 5,506 | 12/12 |
| Inhibitory C6, parent SST | 501 | 7,284 | 12/12 |

The BAM-derived gene counts exactly matched the original source matrix in 324 of 411 sampled mouse/population combinations. The maximum absolute difference was ten UMIs under the stated alignment/tag filters. Raw, per-1,000-nuclei, per-million-library-UMI and per-1,000-Grm7-UMI values are supplied. These descriptive normalizations are not new differential-expression models. Aggregate and constituent-cluster rows overlap and must not be pooled as independent observations.

### Review correction to library normalization

A subsequent comparison with Grm1 identified that the native roster used rendered metadata `nCount_RNA`, which in this saved object is an HVG-only subtotal. It was incorrectly labeled as the full library denominator in the initial per-million export. Reviewed deliverables now use `original_nCount_RNA`, which agrees exactly with `representation_full_library_UMIs`; 72 mouse/population library totals were independently verified against sums of the original all-gene pseudobulk matrices. The correction changes only `total_gene_UMIs` and `per_million_library_UMIs`. All raw counts, per-nucleus values, gene-relative values, feature assignments and junction conclusions are unchanged. The immutable native run is preserved; use the corrected deliverables rather than its original per-million fields. See `NORMALIZATION_CORRECTION.json`, `library_denominator_audit.tsv` and `VALIDATION.json`.

### How much isoform-specific support is present?

The following are conservative uniquely compatible, diagnostic, splice-supported CB/UB counts summed across twelve mice; they are evidence lower bounds, not full isoform abundance estimates.

| Transcript | Superficial CA1 | Deep CA1 | Inhibitory total |
|---|---:|---:|---:|
| Grm7-201 | 1 | 1 | 0 |
| Grm7-202 | 2 | 0 | 2 |
| Grm7-203 | 0 | 0 | 0 |
| Grm7-204 | 6 | 3 | 3 |
| Grm7-205 | 1 | 0 | 1 |

Thus four of the five models have at least some uniquely assignable splice evidence in the requested populations, but the counts are too sparse for reliable abundance or genotype comparisons. None meets the descriptive replicated-support flag of at least three conservative UMIs in at least two mice within a population. This flag does not establish absence when unmet. All per-animal values, rather than only thresholded detections, are retained.

Diagnostic exon/junction-feature support for Grm7-202-associated regions was much greater: 1,302 UMIs in superficial CA1, 884 in deep CA1 and 578 in inhibitory neurons; it occurred in all twelve mice in each population. Exon-only signal can also originate from unprocessed RNA or incompletely annotated transcripts, so these numbers must not be reported as complete mGlu7A transcript counts. Likewise, feature support for the NMD model does not establish increased NMD activity or NMD-transcript abundance.

### Terminal-splicing events

| Event | Superficial CA1 UMIs (mice) | Deep CA1 UMIs (mice) | Inhibitory UMIs (mice) |
|---|---:|---:|---:|
| A-type terminal junction, chr6:111495900–111566017 | 2 (1/12) | 0 (0/12) | 2 (2/12) |
| B-type 92-nt exon-entry junction, chr6:111495900–111501509 | 18 (10/12) | 12 (8/12) | 6 (5/12) |
| Long B-model terminal exit, chr6:111501602–111566017 | 1 (1/12) | 2 (2/12) | 0 (0/12) |
| Partial B-tail terminal exit, chr6:111501602–111502146 | 1 (1/12) | 0 (0/12) | 1 (1/12) |

Coordinates describe the 1-based inclusive CIGAR intron gap. B-exon entry is shared by Grm7-201 and Grm7-205. It demonstrates support for that local splicing event, not unique identification of the complete B receptor transcript. Counts from different junctions can share molecules and are not additive transcript counts. The data are insufficient for a convincing WT–cKO A/B switching claim.

## Interpretation
The existing gene-level Grm7 signal is corroborated by independent BAM feature counting. It cannot be equated with a specific receptor isoform or increased presynaptic receptor abundance. Short-read nuclear libraries contain intronic, partially processed and transcript-shared RNA; the recorded Cell Ranger configurations include introns. The most defensible result is strong Grm7 gene expression with limited terminal-splicing evidence and substantial unresolved full-transcript identity. The A/B event plots are therefore more informative than assigning all gene UMIs to five named transcript bars. No read-based result here establishes an APP–GRM7 or PENK–GRM7 direct-binding mechanism.

## Artifacts
Results: isoform_evidence_per_mouse.tsv; terminal_splice_event_counts_per_mouse.tsv; gene, exon, junction, CDS and UTR count TSVs; transcripts.tsv and features.tsv. All twenty original inhibitory clusters and the CA1 clusters remain selectable or tabulated. The complete coverage grid includes 420 mouse/population combinations: 411 sampled and nine NOT_SAMPLED, whose missing counts are not converted to biological zeros.

QC and visuals: VALIDATION.json; 21 reviewed static SVG/PNG/PDF figures; interactive index.html with population, evidence, normalization, heatmap/dotplot, color, font, dot-size and column-width controls, plus plot/data downloads. Reviewed protein cartoons show only the two complete receptor models; the partial B-tail transcript is explicitly qualified. PNGs were rendered directly from SVG to avoid a high-DPI browser screenshot artifact.

The offline explorer passed Chrome/Playwright checks for all nine evidence views, all 35 selectable populations, missing-population display, appearance controls and SVG/PNG/PDF/TSV/JSON downloads with no JavaScript errors. Downloaded tables and plot-state data were checked against the displayed selection. Both an exported interactive PDF and a static schematic PDF were independently rendered and visually inspected for legibility and complete mouse labels. The browser checklist is saved in EXPLORER_QA.json. The deliverables package contains the reviewed figures, tables, report and corrected source; BAMs and the nucleus-barcode roster are not copied into it.

Provenance: COMPLETE.json; native frozen workflow/configuration; twelve sample JSONs with BAM header/index and called-barcode checksums; matching reference GTF; reference_run_binding.json. The GTF SHA256 is 3a98532da5df438d646e80cf303cbdb01a04bc2edf71e91b9aaa46a9de8952d6. All original run configurations name mm10-2020-A and include introns. The cKO_933 configuration points to its historical VSC reference path; the original reference files at that historical path were not independently rehashed. Its compatible mm10 BAM header, recovered BAM/index, called barcodes and original source matrix were checked, and the Cell Ranger 10.0 versus 10.1 distinction remains documented.

## Limitations and blockers
No access blocker remains and all twelve mice are included. Sparse diagnostic junction support prevents accurate full-length isoform quantification and inferential differential-isoform testing. Zero support does not establish absence. A/B-associated feature counts cannot be used as protein abundance, trafficking or inhibitory-current measurements. The controller reused an existing environment without a complete locked dependency closure; executable versions and source hashes are recorded. Unannotated junctions are provided as descriptive candidates, not reconstructed new transcripts.

## Recommended next step
Validate the shared upstream exon to A-specific terminal junction and the B-type 92-nt exon-entry junction using targeted RT-PCR or junction-specific RT-qPCR on independent CA1-layer/interneuron material, with amplicon sequencing. If full transcript identities and relative isoform abundances are the objective, use targeted long-read cDNA sequencing; only then connect isoform-specific protein localization and presynaptic release measurements to the electrophysiology study.

## Literature context
Corti et al. 1998, alternative mGlu7/8 RNA forms: https://pubmed.ncbi.nlm.nih.gov/9875342/
Shigemoto et al. 1997, hippocampal presynaptic localization: https://pmc.ncbi.nlm.nih.gov/articles/PMC6573434/
