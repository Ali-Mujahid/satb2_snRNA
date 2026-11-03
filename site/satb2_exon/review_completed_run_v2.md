Takeaway: CA1 Satb2 exon depletion and skipping distinguish all six cKO mice from six WT mice.

This analysis compared Satb2 candidate-exon and splice-junction evidence in CA1/ProS excitatory nuclei from six wild-type and six conditional-knockout mice. All twelve sample extractions and the descriptive outputs completed, and an independent review reconciled the sample counts with the mouse-level table. Candidate-exon support was low in every cKO mouse and higher in every WT mouse, while the exon-skipping junction was detected in every cKO mouse and not observed in WT under the stated filters. This is strong genotype-associated evidence in the sequenced CA1/ProS nuclei that is consistent with the intended exon-targeting strategy. It does not prove the genomic recombination breakpoint, complete deletion in every targeted cell, the targeted construct's exact coordinate identity, or loss of functional SATB2 protein. Confirm the floxed construct coordinates and use an appropriate protein assay when making a complete-knockout claim.

## Scientific context and question

Twelve independent P60 mice, six WT and six cKO, were analyzed using their frozen Seurat original-full CA1/ProS candidate, excitatory-reference nucleus assignments. Whole-CA1/ProS membership includes the previously identity-conflicted C11 population; its exclusion is a separately preserved sensitivity. The assay is 10x single-nucleus RNA-seq, and corrected sample-specific `(CB, UB)` support is descriptive RNA evidence rather than the number of recombined cells.

## Question, lifecycle, and endpoint

Workbench registry `6978262e-6ef6-461f-8d9b-aaf43fed8806`, target `lisc-satb2-exon-local-20261103` and validated target hash `sha256:ae22509e459dcb5cda2abc99f509bf7e246d4d0851d81a0bd2ee2572d2dc1f93`, run root `/lisc/data/scratch/molecular_neurosciences/i-med/new_snRNA_seq/SATB2_CA1_EXON_RUN_20261103_v2`, approved plan checksum `sha256:ba828c46ca942cbe761062513da44380673320ee4a22f93cb1d4b194b5d808fc`. The run finished with return code 0; native Snakemake reports 18/18 steps finished and all twelve sample JSONs, 24 summary rows and requested plot files are present. Workbench's task-attempt projection labels 18 attempts aborted even though overall lifecycle is completed and the native log reports 18 finished jobs; the cause is unconfirmed and the discrepancy must not be ignored or used by itself to rerun the analysis.

## Key findings

In the all-CA1/ProS cohort the candidate exon had 0–3 distinct barcode/UMI pairs in each cKO mouse and 21–29 in each WT mouse. The Satb2-202 E3→E5 skipping junction had 21–46 distinct pairs in each cKO mouse and zero observed pairs in each WT mouse under the stated filters. Sample-specific values (candidate exon / skipping junction) were: cKO_1138 2/26, cKO_1139 3/42, cKO_1160 1/46, cKO_1165 0/35, cKO_933 2/21, cKO_991 0/27, wt_1122 29/0, wt_1123 25/0, wt_1178 24/0, wt_1220 21/0, wt_1226 29/0 and wt_1228 28/0. All mice had detectable Satb2 locus support and all 27,703 frozen nuclei matched BAM-associated called barcodes. The C11-excluded set removes 99 identity-conflicted nuclei and retains the same qualitative group pattern. Excluding cKO_933 leaves six WT and five cKO mice and the same qualitative pattern; it is an influence sensitivity, not an isolated software-version effect.

## Interpretation

The opposing per-mouse pattern supports genotype-associated exon skipping in the existing CA1/ProS cohort. It is consistent with a targeted conditional Satb2 allele and argues against a wholesale WT/cKO swap. The exon is defined by actual mm10 coordinates chr1:56898991–56899117 (1-based inclusive), not solely by an exon number; its identity relative to the exact targeting construct remains to be confirmed. No randomization or molecule-level inference was substituted for mouse-level replication, and no formal genotype P value or recombination efficiency was estimated.

## Artifacts

- Results: `<run_root>/results/mouse_evidence.tsv` and twelve `<run_root>/results/samples/<sample>.json` files.
- QC and visual reports: four PDF/PNG pairs for all-CA1/ProS and C11-excluded sets, with and without cKO_933. Initial visual inspection found overlapping text labels for near-zero mice; data remain available in the table and an improved additive visualization is needed.
- Provenance: source, approved plan and controller log remain under the run root. Local `independent_cohort_audit_20261103.md` and the independent output audit document input identity, barcode matching, counts and denominators. The source outputs honestly preserve `requires review` fields; the independent audit supports source matching without silently rewriting those native fields.
- Model synthesis: this local review.

## Limitations and blockers

The genomic construct map, complete full-BAM integrity, individual nucleus recombination, transcript reading frame and protein abundance were not measured. Junction counts are sparse UMI-supporting observations, zero means not observed under these filters, and CB/UB features are not necessarily distinct full-length molecules. cKO_933 uses an original Cell Ranger 10.0.0 run, whereas eleven samples use 10.1.0 reruns; excluding that animal cannot distinguish software version from mouse variation. CA1/ProS and layer labels do not establish anatomical position. The managed controller lacks an observed lockfile and host platform verification, though LiSC's existing environment and samtools versions were independently observed. Only targeted indexed locus queries and quickcheck, not complete BAM integrity, were performed. A user-specified mutant transcript may still produce residual/truncated protein; no protein or functional assay is present here.

## Recommended next step

Publish the reviewed mouse-level table and clear additive visuals with stated source and limitations, then confirm the exact floxed coordinates and compare WT versus cKO SATB2 protein in CA1 pyramidal neurons. Preserve original analyses, native figures and this separate run.
