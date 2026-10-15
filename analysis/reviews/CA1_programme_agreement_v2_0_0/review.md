Takeaway: CA1 programme names differ across graph partitions and operational thresholds; anatomy remains undetermined.

This review uses completed native original-full WT programme scores and final nucleus labels. It performs only metadata arithmetic and visual sensitivity; primary programme labels, counts and DE memberships are unchanged. Same-name overlaps are separate from diagnostic Hungarian label renaming and cannot be substituted by differently named programme pairs.

```text
       branch  leading_WT_score_threshold         programme  Scanpy_nuclei  Seurat_nuclei  same_named_overlap_nuclei  Scanpy_overlap_fraction  Seurat_overlap_fraction  same_named_Jaccard  reciprocal_overlap_ge0p5_and_n_ge20  physical_layer_assignment_validated
original_full                         0.3              Deep           6297           8750                       6017                 0.955534                 0.687657            0.666334                                 True                                False
original_full                         0.3       Superficial           5862          11706                       5258                 0.896963                 0.449171            0.427132                                False                                False
original_full                         0.3        Mixed/weak          15451           7148                       6509                 0.421267                 0.910604            0.404537                                False                                False
original_full                         0.3 Identity conflict             93             99                         93                 1.000000                 0.939394            0.939394                                 True                                False
```

Fine-reference evidence in the native cluster review refers to the shared CA1-ProS atlas subclass, not validation of separate fine graph states, atlas supertypes or physical layers. Supplied VEH programs are associated expression evidence, with Satb2 excluded. The leading WT-score rule remains operational: borderline .200176/.308745 superficial scores make Scanpy2/Seurat2 assignments threshold-sensitive.

Any native cross-tool effect row pairing different programme names is population-definition sensitivity only. A claim about deep-versus-deep or superficial-versus-superficial must use this same-name table, reciprocal membership support, actual genotype model availability and count/gene-universe bindings. Fine programme models lack core-total partition/LOO sensitivities; neither matching effect signs nor atlas confidence establishes anatomical layer loss.

The PDF presents exact original-cohort semantic confusion matrices at thresholds .2/.3/.4 and programme size sensitivity. Fixed guards retain a score gap of .5, at least three program genes detectable in 10% of WT nuclei, at least20 WT nuclei/three represented WT mice, and explicit conflicting identity groups. The primary .3 reconstruction matched every actual final nucleus programme label. TSVs also retain paired raw/corrected common cohorts, same-name overlap fractions, ARI/NMI and diagnostic-only Hungarian renaming fractions. All physical-layer validation flags remain false.
