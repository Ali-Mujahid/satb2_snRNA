Takeaway: Grm1 alpha/beta RNA evidence is recovered in all twelve mice, but terminal-junction support is sparse.

This analysis asked which Grm1 alpha/beta-associated RNA features are supported in P60 Satb2 WT and cKO CA1 pyramidal and inhibitory populations. The corrected workflow completed for all six WT and six cKO mice using the same frozen Seurat memberships, BAMs, alignment/tag filters and matching mm10-2020-A reference as the Grm7 analysis. Sequence verification maps Grm1-201 to the 1199-amino-acid alpha protein and Grm1-202/-203 to identical 906-amino-acid beta proteins, while retaining three competing noncoding models. Gene counts closely reproduce the original Seurat matrices and Grm1 gene RNA is enriched in the original SST/O-LM-like C2 population relative to either pyramidal sublayer. Conservative diagnostic-plus-spliced alpha RNA assignments recur in C2, but beta assignments and direct alpha/beta terminal-junction counts are sparse, particularly in C6. These findings support defined RNA features and motivate targeted validation; they do not measure receptor protein abundance or establish a WT/cKO isoform switch.

## Scope, execution and provenance
Primary groups: superficial CA1 C2+C3 (10,867 nuclei; C8 excluded), deep CA1 C1+C5 (8,750), all inhibitory neurons (8,911), inhibitory C2 SST/O-LM-like (1,307) and original parent C6 SST (501). The other original clusters remain intact and available. Mouse is the biological replicate. There are 35 populations and 420 mouse/population combinations, of which 411 are sampled and nine are explicitly NOT_SAMPLED rather than biological zero. No clustering, annotation, mouse inclusion or DESeq2 fit was changed.

Corrected native run: bcb2f2ec-930a-4240-9baf-6a1221c98f58, completed with return code zero. Source version: version-d5b7aca4dc9d211494f089d247747ef1; source SHA256 5861593eb4cbd19135181b409c8c3f0e430f0a4ea1f311a63f7751849c3358af. Plan checksum: 67bea7386337f138d3a3233af2f3468a9f02613340c75d9be89debfd18ba6142. Results are under /Users/mujahid357/Documents/codex_new/SATB2_GRM1_ISOFORMS_20261126/run_v2/results. The initial run stopped at an overly restrictive three-model assertion, before BAM counting; it produced no biological results and is preserved as a failed preflight, not evidence of absent isoforms.

The existing office-Mac Snakemake9.27 controller coordinated bounded read-only LiSC locus processing through Python/SAMtools. No remote analysis files, raw BAM transfers, installations or previous-result changes were required. The reused environment has executable-version/source records but no complete locked dependency closure.

## How transcript IDs are mapped to alpha and beta
The actual GTF identifies Grm1 as ENSMUSG00000019828 on the negative strand of mm10 chromosome10:10686059–11082356. Coding intervals were concatenated in genomic order, reverse-complemented and translated with the initial GTF phase selected in transcript order. Alpha has 1199 amino acids and a tail ending DYKQSSSTL; beta has 906 amino acids with the distinctive KKRQPEFSPSSQCPSAHVQL terminal sequence. The two beta models have identical translated-sequence hashes. These assignments are established from sequence and annotation rather than guessed from transcript suffixes or antibody intensity.

| RNA model | Stable transcript ID | GTF biotype | Interpretation |
|---|---|---|---|
| Grm1-201 | ENSMUST00000044306 | protein_coding | alpha, 1199 aa |
| Grm1-202 | ENSMUST00000105560 | protein_coding | beta, 906 aa |
| Grm1-203 | ENSMUST00000105561 | protein_coding | same beta protein, different UTR model |
| Grm1-204 | ENSMUST00000135120 | lncRNA | no annotated receptor CDS |
| Grm1-205 | ENSMUST00000155772 | lncRNA | no annotated receptor CDS |
| Grm1-206 | ENSMUST00000156826 | lncRNA | no annotated receptor CDS |

Thus six RNA models do not mean six receptor proteins. Grm1 alpha/beta are unrelated to the Grm7-201/-202 transcript labels or mGlu7A/B protein names.

## What Figure 4B in the supplied paper means
Ferraguti et al. (1998), J Comp Neurol400:391–407, used isoform-specific antibodies and transfected/null controls. Figure4B is the hippocampal synaptosome Western-blot panel: approximately150-kDa alpha and97-kDa beta immunoreactive products were detected, and the pan-mGluR1 antibody detected both. “4B” is a panel label, not a receptor subtype. Hippocampal synaptosomes pool subregions and cell sources; this panel does not demonstrate beta expression specifically in CA1 pyramidal cells. The paper's localization experiments placed alpha staining prominently in non-principal neurons and beta staining in CA3/dentate principal cells, without detectable CA1 beta staining under their conditions. RNA assignments in our study cannot directly confirm or refute those protein-localization measurements.

## Main observed results
Conservative assignments below are corrected CB/UB evidence counts summed across the six mice in each genotype. They require diagnostic RNA features, compatibility exclusively within the alpha or beta family, and compatible anchored splicing, with incompatible splice evidence excluded. They do NOT necessarily contain the isoform-discriminating terminal junction; for example, a model-specific UTR segment and a compatible constitutive splice can jointly support an assignment. They are lower-bound RNA evidence, not complete mature-transcript or protein abundance.

| Population | Alpha WT | Alpha cKO | Beta WT | Beta cKO |
|---|---:|---:|---:|---:|
| Superficial CA1 C2+C3 | 27 | 36 | 1 | 1 |
| Deep CA1 C1+C5 | 8 | 14 | 1 | 2 |
| All inhibitory neurons | 14 | 14 | 0 | 0 |
| Inhibitory C2 | 9 | 6 | 0 | 0 |
| Inhibitory C6 | 0 | 1 | 0 | 0 |

C2 alpha assignments occurred in5/6WT and4/6cKO mice; C6 had one assignment in one cKO mouse. Beta diagnostic exon-feature signal is not sufficient to establish processed beta transcripts, and zero conservative assignments do not establish beta absence. These descriptive totals must not be interpreted as genotype fold changes: sampled nuclei and library sizes differ. Raw and three normalized values per mouse are supplied.

Mean full-library gene CPM (WT/cKO): superficial433.8/370.8; deep140.0/133.3; inhibitory total459.2/457.5; C21364.6/1395.1; C6912.6/1000.4. C2 gene-level enrichment relative to CA1 is qualitatively consistent with SST-associated mGluR1 literature, but does not identify an alpha protein level or prove O-LM morphology.

### Direct terminal-splicing evidence
Coordinates below are 1-based inclusive CIGAR intron gaps. Grm1 is negative-strand, so the transcript direction runs from higher to lower genomic coordinates.

| Junction | Model/family | Superficial UMIs | Deep UMIs | All inhibitory UMIs | IN C2 UMIs | IN C6 UMIs |
|---|---|---:|---:|---:|---:|---:|
| chr10:10689903–10719222 | alpha-specific, skips beta exon | 0 | 1 | 2 | 2 | 0 |
| chr10:10698133–10719222 | beta 85-nt exon entry | 2 | 2 | 0 | 0 | 0 |
| chr10:10689903–10698047 | beta exon exit | 1 | 1 | 0 | 0 | 0 |

These rare junction counts establish local RNA splice-event support, not a whole beta-to-alpha ratio. Counts from adjacent junctions can share a molecule and cannot be summed into full transcript counts. The observed sparse beta signal in CA1 is not proof of detectable beta protein or evidence that the older immunohistochemistry is wrong. No differential isoform-switching hypothesis test is supported by these sparse terminal events.

## Methods and independent checks
Primary MAPQ255 alignments, exclusive GX=Grm1, corrected CB/UB, mm!=1 and xf bits2/4/32 excluded. Each feature is deduplicated by mouse/population/(CB,UB); annotated junctions require at least8 contiguous aligned bases on each side. A diagnostic bin requires at least20 aligned bases or a family/model-exclusive anchored junction. Shared compatibility upper bounds overlap and are not additive. All six RNA models compete for shared fragments. Unannotated junctions are retained as unvalidated candidates rather than reconstructed transcripts.

Normalization uses the full original source-library original_nCount_RNA, agreeing with representation_full_library_UMIs, and never the HVG-only rendered nCount_RNA. Source gene comparison:390/411 combinations match exactly; maximum absolute BAM/source difference is2UMIs. Synthetic tests cover shared beta model ambiguity, corrected UMI deduplication, short anchors, incompatible splices, tag exclusion, JSON coordinate transfer and negative-strand sequence orientation. Independent validation verifies gene/family UMI conservation, cluster aggregation, all twelve sample identities, six models, identical beta protein hashes and normalization arithmetic. All372 original cluster/mouse full-library totals also match independently summed all-gene matrices. The GTF SHA256 is3a98532da5df438d646e80cf303cbdb01a04bc2edf71e91b9aaa46a9de8952d6.

cKO_933 retains the recovered original Cell Ranger10.0 BAM versus10.1 for the other eleven mice. All prior run configurations named mm10-2020-A; the historical VSC reference-path files for933 were not independently rehashed, as documented in the Grm7 audit. Compatible BAM headers, reference chromosome length, original sample IDs, called-barcode membership and source counts were checked.

## Artifacts and review
Tables include alpha_beta_evidence_per_mouse.tsv, alpha_beta_specific_junction_counts_per_mouse.tsv, gene_UMI_counts_per_mouse.tsv, isoform_evidence_per_mouse.tsv, exon/junction/CDS/UTR feature counts, population_genotype_master_summary.tsv and exact feature/transcript definitions. The atlas offers all35 populations, twelve evidence views, per-mouse raw/full-library/nucleus/gene normalization, adjustable colors/fonts/column widths/dots, and SVG/PNG/PDF/TSV/JSON exports. There are28 static figures in SVG/PNG/PDF. Browser/export checks and independent count validation are saved separately. Plot styling does not refit an analysis or change saved values.

## Interpretation and next experiments
The biological conclusion is strong Grm1 gene RNA in SST-associated populations with limited isoform-resolving nuclear read evidence. C2 provides more repeated alpha-associated RNA support than C6, while beta-specific processed RNA support in the requested groups is very sparse. Neither the relative number of diagnostic UMIs nor their absence establishes receptor membrane abundance, trafficking or a change in inhibition. Prioritize alpha/beta-junction RT-PCR or RT-qPCR with amplicon sequencing on independent layer/interneuron material, then isoform-specific protein localization with appropriate specificity controls. Targeted long-read cDNA is the next method if complete transcript identity and alpha/beta proportions are essential. Electrophysiology can subsequently test receptor-dependent excitability/release, but its prediction cannot be assigned from total Grm1 counts alone.

## Primary literature
- Ferraguti et al.1998: https://pubmed.ncbi.nlm.nih.gov/9779943/
- Techlovska et al.2014: https://pmc.ncbi.nlm.nih.gov/articles/PMC4188797/
- Mouse alpha/beta protein-sequence context: https://pmc.ncbi.nlm.nih.gov/articles/PMC5129514/
