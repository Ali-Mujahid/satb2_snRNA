# SATB2 CA1–inhibitory ligand–receptor results

Seurat-defined CA1–inhibitory analysis prioritizes 324 replicated, context-specific LR candidates.

Using six WT and six Satb2 cKO mice, this analysis asked which ligand–receptor pairs link the Seurat-defined CA1 programs with inhibitory neurons. Native LIANA consensus scoring with its CellPhoneDB method completed for pooled and genotype-specific data, and the final figure package passed validation. A stricter exploratory tier retains 324 population-specific directed candidates, 322 of which have CA1-only gene-level DESeq2 support. The two candidates with both-side support link the Sst-like inhibitory population's Cntn3 to Ptprg in each CA1 program, with oppositely signed component RNA changes. These results prioritize candidates for follow-up rather than establishing changes in physiological communication. The useful next step is to inspect the affected component genes and validate protein availability and colocalization.

## Scientific context and question

The existing Seurat original_full definitions were used without re-clustering or re-annotation:

| Population | Nuclei | WT nuclei | cKO nuclei |
|---|---:|---:|---:|
| CA1 deep-associated expression program | 8,750 | 4,485 | 4,265 |
| CA1 superficial-associated expression program | 11,706 | 6,278 | 5,428 |
| Inhibitory, all 20 Seurat clusters | 8,911 | 4,464 | 4,447 |
| Combined, distinct nuclei | 29,367 | 15,227 | 14,140 |

The main analysis uses inhibitory C1–C20 separately; a complementary resolution uses inhibitory total. Total and constituent clusters are never placed together in the same fit. CA1 mixed/weak and conflicting layer programs are outside the requested scope. Deep and superficial are molecular expression programs; physical layer anatomy was not measured.

The three requested views are CA1→inhibitory, inhibitory→CA1, and the union of both directed sets. The combined view does not establish reciprocal physiological signaling and is not independent replication.

## Completed analysis and methods

All six native consensus fits and 24 per-mouse geometric-mean support calculations completed. The final registered run is `d397e9d9-03bc-43b0-9bd8-901409248fd6`, status **completed**. This atlas packages the verified saved results; it performs **zero model refits**. Color, geometry, label and size controls affect visualization only. Native scores and previous DESeq2 exports retain their contents.

Method: fixed official LIANA 2.0.0, immutable commit `c59472ccc9de8360dbbf5016db75f8abde08dd3e`, native CellPhoneDB/Connectome/logFC/NATMI/SingleCellSignalR consensus, and frozen mouseconsensus. This is LIANA's CellPhoneDB scorer, rather than an independent standalone CellPhoneDB v5 run. The resource has 3,989 ordered LR pairs; 3,946 have all subunits represented among input genes. The native method also excludes unexpressed features and applies the detection gate.

Full native integer RNA counts were joined by nucleus ID to frozen Seurat labels. Raw duplicate-symbol features were summed before log1p normalization to 10,000 UMIs, using the full RNA library denominator. All 32,245 collapsed features remain available for native background-dependent scores. Input preparation conserved all 322,680,709 selected UMIs. There are 38 ambiguous symbols involving 78 source feature rows; ambiguous gene-level DE evidence does not supply primary DEG support.

Parameters: expression fraction ≥10% for every complex subunit, min_cells=5 within a scored group, 1,000 native cell-label permutations, seed=1337. All nuclei remain in pooled/genotype fits. Mouse support records distinguish eligible mice from mice passing the expression gates; missing or ineligible mice are not treated as zero expression.

## Key findings

There are **64,936 scored records** across the six fits and **226,930 per-mouse support records**. The same molecular LR pair in different source/target populations is counted as a separate record.

| Resolution | Pooled expression-qualified records | WT | cKO |
|---|---:|---:|---:|
| Inhibitory clusters | 20,471 | 19,467 | 21,680 |
| Inhibitory total | 1,108 | 1,040 | 1,170 |

These counts are expression-compatible predictions, not counts of significantly altered interactions. They should not be compared as evidence of increased/decreased communication between genotypes.

The initial pooled cluster view has 1,589 records with magnitude_rank≤0.05 and at least one matching endpoint's DESeq2 DEG: 691 CA1→IN and 898 IN→CA1. Adding ≥3 supporting mice per genotype and a reference-supported inhibitory endpoint retains 1,126. Adding native CellPhoneDB P≤0.05 retains 734; additionally requiring aggregate specificity_rank≤0.05 yields the stricter **324-record tier** displayed by the report's “Replicated + context-specific DEG candidates” option:

| Direction | Deep program | Superficial program | Total |
|---|---:|---:|---:|
| CA1→inhibitory | 41 | 71 | 112 |
| Inhibitory→CA1 | 83 | 129 | 212 |
| Both directed sets | 124 | 200 | 324 |

Of those 324, **322 have CA1-only DEG evidence**, and two have gene-level DEG evidence at both ends. No inhibitory-only DEG-supported record passes all of these stricter criteria. This tier is an exploratory prioritization, not a new FDR-controlled genotype interaction test.

The four total-inhibitory DEGs have **no overlap with this LR resource**. Consequently, inhibitory-total DEG-supported interactions receive their DEG evidence from the CA1 endpoint. Within inhibitory clusters, only C1 and C6 contain resource-overlapping DEGs. C18's original DESeq2 fit was not tested and remains missing, rather than zero DEGs.

The clearest both-end gene-evidence candidate is **C6 Cntn3→CA1 Ptprg**:

| Component and matching population | DESeq2 log2(cKO/WT) | Adjusted P |
|---|---:|---:|
| Cntn3, inhibitory C6 (Sst-reference-supported) | −0.53835 | 0.0037006 |
| Ptprg, CA1 deep program | +0.39666 | 0.016489 |
| Ptprg, CA1 superficial program | +0.49044 | 2.9271×10⁻⁸ |

Both population-specific interaction records have expression support in all six WT and all six cKO mice. Opposite component RNA changes do not identify a single net increase or decrease in interaction activity.

Other useful examples in the broader, explicitly curated 16-record table include CA1 Ntng1→C1 Lrrc4c and C2 Efna5→CA1 Epha6. Their DEG evidence is mainly CA1-side. Deep Slit1→C6 Robo1 has changes in both components in the broader CellPhoneDB-specific shortlist, but does not pass every stricter aggregate-specificity criterion. Selection policies accompany both tables so those tiers cannot be confused.

## Interpretation

The strongest supported conclusion is that CA1-program gene changes intersect a broad set of transcript-compatible curated interaction priors, while inhibitory evidence is concentrated in a small number of subclusters. C6 Cntn3→Ptprg deserves focused review because both matching component endpoints change and all mice support expression. C6 should be described as Sst-reference-supported; its finer backprojection-like name is a molecular resemblance, not measured morphology or connectivity.

High magnitude does not guarantee cell-group specificity. Some Nrxn/Nlgn endpoints pass specificity gates and others are ubiquitous. Curated reversed orientations can have identical scores and must not be treated as independent physiological directions. Native zero permutation P values mean below the 1,000-permutation empirical resolution, not an exact zero probability.

Resource inclusion also needs biological review: for example, retain Rims1→Cacna1c in complete native exports without presenting Rims1 as an established secreted extracellular ligand. All figures remain prior- and transcript-supported hypotheses. The protein LR resource incompletely represents neurotransmitter and metabolite communication.

## Interactive views, evidence and downloads

The [ligand–receptor atlas](index.html) displays all six saved fits. Choose inhibitory clusters or total, pooled/WT/cKO data, CA1→inhibitory, inhibitory→CA1 or both directed sets, and an evidence tier. Heatmap, dotplot, circos and circle views use the selected saved records. Colors, color limits, plot sizes, arrows, arcs and labels are display settings. They do not alter expression, scores, DESeq2 fits, annotations or genotype inference.

The default DEG-supported tier requires magnitude_rank≤0.05 and a matching unambiguous endpoint gene DESeq2 adjusted P<0.05. The stricter replicated/context-specific tier additionally requires specificity_rank≤0.05, native CellPhoneDB P≤0.05, at least three supporting WT and three supporting cKO mice, and reference-supported endpoint identity. All thresholds and selection rules should accompany exported figures. Changing a threshold selects saved candidates; it does not calculate a new statistical test.

Heatmap and network counts represent population-specific directed LR candidates after the current filters. A tested, eligible source/target pair can have zero retained candidates. Untested directions are missing/NA, rather than zero. A combined view is the union of directed sets, not an independent reciprocal-signaling test. The separate total-inhibitory fit cannot be combined with its constituent clusters.

Dot area can represent the minimum detection fraction across **every** complex subunit (`lr_prop`). This can differ from the native `ligand_props` or `receptor_props`, which refer to the native mean-expression-limiting subunit. Fractions slightly above one from sparse floating-point summation are clamped only for display; the original saved numbers remain intact. Mouse support is shown with eligible denominators, so unavailable mice are not treated as zero.

The [native seven-column interaction tables](data/full_result_tables/index.html) are exactly `source, target, ligand_complex, receptor_complex, lr_logfc, lr_gmeans, magnitude_rank`. Full tables and separately labeled DESeq2-supported tables are available for each resolution, condition and direction. The [original gene DESeq2 exports](data/existing_DESeq2_exports/index.html) keep their standard seven-column format and gene symbols unchanged. Extended native evidence retains component stable IDs, signed cKO/WT effects, adjusted P values, identity flags and mouse denominators.

Native download links point to the immutable atlas-v3.5.0 repository tag. They require internet. The complete portable results ZIP attached to the GitHub release contains all native tables and the original 258 primary plus 18 supplemental precomputed PDF/PNG/SVG panels for offline use. Interactive selected tables and modified chart exports are generated directly by the atlas.

## Score meanings

| Field | Interpretation |
|---|---|
| lr_logfc | Native cell-group-versus-rest specificity; **not** cKO/WT fold change |
| lr_gmeans | Geometric mean of endpoint mean log-normalized expression with native complex bottleneck rules |
| magnitude_rank | Lower means higher consensus magnitude priority; **not** genotype FDR |
| specificity_rank | Native aggregate cell-context specificity rank; **not** genotype evidence |
| cellphone_pvals | Native LIANA CellPhoneDB cell-label specificity P; zero is below the 1,000-permutation empirical resolution |
| deg_supported | At least one unambiguous matching sender-ligand or receiver-receptor gene has saved DESeq2 adjusted P<0.05 |
| lr_prop | Minimum detection fraction across all ligand and receptor complex subunits |
| sample_support / sample_eligible | Independent mice passing expression gates / mice eligible for that population pair |

Genotype score/rank differences are descriptive only. No replicate-aware differential-interaction model was fitted, so **no interaction is called significantly altered between WT and cKO**. Existing DESeq2 support is reused gene-level evidence, not an independent validation of an interaction.

Verified sources include the [native result counts](https://raw.githubusercontent.com/Ali-Mujahid/satb2_snRNA/atlas-v3.5.0/communication_data/full_result_tables/RESULT_COUNTS.tsv), [extended interaction evidence](https://raw.githubusercontent.com/Ali-Mujahid/satb2_snRNA/atlas-v3.5.0/communication_data/full_result_tables/interaction_evidence.tsv.gz), [scientific review](data/reviews/SCIENTIFIC_RESULT_REVIEW.md), [strong 324-candidate table](data/reviews/strong_consensus_specificity_candidates.tsv), [curated 16-candidate selection policy](data/reviews/curated_16_selection_policy.json), [native final validation](data/reviews/NATIVE_FINAL_VALIDATION.json), [dataset manifest](data/manifest.json) and [transport checksums](data/TRANSPORT_RECEIPT.json). Independent checks covered 154 input checks, 43 native scientific checks, visual review of all chart types and existence of all 276 panel asset sets. Interactive transport independently reconstructs all 64,936 native records and exact component evidence.

## Scientific limitations

Nuclear RNA does not measure secreted protein, receptor surface localization, cell contact, receptor activation/inhibition or functional synaptic communication. No replicate-aware differential-interaction model was fitted, so **no interaction is called significantly altered between WT and cKO**. Existing DESeq2 support is reused gene-level evidence, not an independent validation. Mixed/uncertain groups remain explicit, and C20 has an immune-associated identity-review flag.

## Recommended next step

Review the both-end Cntn3/Ptprg candidate and the annotated source/target gene tables first. A subsequent differential-interaction model or protein/spatial validation would answer a different question and is not silently added to this run.

Method sources: [official LIANA 2.0.0 release](https://github.com/scverse/liana/releases/tag/V2.0.0), [LIANA documentation](https://liana.readthedocs.io/en/stable/), and [Rims1 reference annotation](https://www.ncbi.nlm.nih.gov/datasets/gene/116837/). Project-specific findings come from the native saved outputs and matching original DESeq2 endpoints, not the example screenshots.
