# Independent Seurat runtime/result review

Completed native run: `runs/20261013_v1_0_5`, run ID
`9702cbef-00e6-4f64-aa37-e9b62357bd71`.

Review checkpoint: **13 October 2026, 14:26 UTC / 16:26 Vienna**. This is a
read-only review of actually published artifacts. Review supplements use
existing coordinates, metadata and fitted result tables. They do not refit
models or modify frozen workflow, source inputs or native run outputs.

## Current scope and conclusions

All **13 Seurat representation datasets** are complete and reviewed: global,
EX, IN and CA1/ProS-candidate log representations in each of the three count
branches, plus original-full global SCTransform. Each includes unintegrated
and Harmony PCA/graph/cluster/UMAP results. The original-full roster has
108,062 nuclei; paired-original and paired-CellBender each have the same
104,298 usable nuclei. Exact IDs, sample/condition identities, selected-source
full-library totals, stored HVG counts, feature panels and finite coordinate
axes pass the independent integrity checks recorded below. Thirteen explicit
QC metadata sidecars distinguish full-gene measurements from Seurat's retained
HVG-assay `nCount_RNA`/`nFeature_RNA`.

All three **marker stages are complete** (original-full at 14:13:23 UTC),
and the 12-row sampling/source-label manifest plus actual marker output review
are published and verified. Fresh paired representation clusters describe graph sensitivity.
Paired marker and genotype-DE models instead use frozen original-full cluster
memberships restricted to the paired IDs. Reused MapMyCells predictions are
source-verified original annotations, not freshly remapped labels. SCTransform
is a representation/model/feature-panel sensitivity with 2,731 shared modeled
features; this run has no separate SCT corrected-count marker or genotype-DE
analysis.

All three **pseudobulk producers and the native workflow are complete**:
4,616 endpoints comprise 2,410 tested and 2,206 explicitly not tested, with
zero error or captured warning rows. Exact mouse identities, endpoint-specific
inclusion, replication/nucleus bounds, full metadata design ranks/residual df
and gene/P/FDR denominator guards pass for all tested models. The source-bound
navigation preserves every original population, including zero or ineligible
paired definitions. It provides 3,426 PRIMARY rows with actual audit/result
paths and explicit NA status.

The final native comparator reports 25 datasets, 222 partition comparisons,
no review flags and exact independent R/Python aggregation for all 77 covered
populations per branch. That count scope covers the four Atlas totals and
Scanpy Harmony families. R-only Seurat/unintegrated families have native
conservation checks in addition to original-membership binding; they are not
claimed to have independent Python count validation.

The four main Atlas totals show closely matched primary DESeq2 poscounts/Wald
and edgeR TMM/robust quasi-likelihood effects (Spearman 0.9979–0.9995), with
different BH lists. DESeq2 partition/leave-one-mouse-out sensitivities preserve
estimable strong-hit directions while weak/all-gene rankings change more.
The final same-roster CellBender sensitivity has effect Spearman 0.8895–0.9601;
raw BH-hit directions persist but magnitudes, uncertainty and FDR availability
change. Layer-associated genes show heterogeneous expression responses rather
than a uniform sublayer reassignment. The 571 overlapping definitions per branch
remain exploratory outside the main totals and do not represent independent
biological hypotheses.

Interpretation must retain three material limits. First, Harmony strongly
changes the fine EX/CA1-ProS condition-associated partitions; coarse global
atlas purity cannot establish preservation of genotype effects. The global
original log, paired-original log and SCT Harmony fits reached the ten-iteration
cap without reported convergence; finite completed coordinates do not remove
that limitation. Second, primary known/main partition contrasts are supported
mainly by the shared main sequencing partition (six WT and three cKO), so
near-identical adjusted/main effects are not independent confirmations.
Third, positive gene-body Satb2 RNA effects do not test exon-4 retention,
protein or function, and do not reopen the study's settled genotype validation.

The detailed representation and completed core pseudobulk evidence follows. Runtime
repair history is retained as background in the appendix; v1.0.5 completed
actual-data outputs supersede the earlier failures.

## Verified prerequisites

- The active v1.0.5 frozen Seurat script was independently hash-checked against
  the reviewed current source; both SHA256 values are
  `0ee1920c6df4ca406e4482269dbd671d1e166a443d964041fbea9f1ca56ba1a4`.
  The frozen pseudobulk module also matches the current reviewed source,
  SHA256 `319f799d69f3f84d837c2261bd7cd62934e235d7eb8a95e94da5cff69e1ad126`.
  The preceding importer-only repair retained source SHA256
  `6411768a18e5e35a7a1121c5b47b464266bc388662293fed45c1af2e1244ea00`.
  Earlier v1.0.1/v1.0.2 source receipts retain the pre-repair script SHA256
  `4e76b4260c9f3457d48376d56c072b466569d5c4803ccb5647f289c6368abb4d`.
- The executed preparation receipts report 12 source mice, 32,285 stable genes,
  108,062 original-full nuclei and 104,298 usable paired nuclei.
- All twelve recorded original MapMyCells query checks report exact original
  integer count equality, matching genes/nucleus IDs, and zero unequal sparse
  entries. Annotation reuse is explicitly recorded; no new mapping is claimed.
- Actual common-HVG TSV SHA256 values were independently recomputed and agree
  between original-full and paired-CellBender for all four populations. Each
  recorded panel contains 3,000 genes.
- The frozen R config uses the intended full-library normalization, original
  library metadata, explicit Harmony 2 settings and original-full label/cluster
  sources for subsequent paired-count inference. Genotype is not a Harmony
  nuisance variable.

| Frozen population | Original-full nuclei | Paired nuclei |
| --- | ---: | ---: |
| All nuclei | 108,062 | 104,298 |
| High-support EX | 60,985 | 60,849 |
| High-support IN | 8,911 | 8,903 |
| High-support CA1/ProS candidates | 27,703 | 27,672 |

These are declared source-cohort counts, not fresh cluster numbers or anatomical
population validation. Most paired attrition lies outside the stated neuronal
subsets; its interpretation requires the per-cell QC/annotation review.

The four common Atlas totals are the main population summaries in the
pseudobulk review. Wider overlapping Seurat/Scanpy cluster families at the
configured primary resolution remain exploratory, with endpoint-specific
BH rather than adjustment across all 571 population hypotheses. Resolution
0.2/0.6 representations are additional clustering sensitivities, not added
primary model endpoints by default.

SCT is a representation/model/feature-panel sensitivity only in this run.
There is no separate SCT corrected-count marker or genotype-DE analysis;
`do.correct.umi=FALSE`, and the durable object does not retain a corrected
SCT UMI assay. Native marker defaults use the Harmony resolution-0.4 global
and neuronal-subset R labels with full-source library log normalization.
Paired Seurat marker and pseudobulk models use original-full source labels
restricted to the usable paired nuclei. Fresh paired clustering outputs are
used for representation comparisons and are a different artifact scope.

## Completed paired-CellBender global representation

Both 2D UMAP exports retain the exact 104,298-nucleus metadata order and contain
only finite values. Original-library, corrected-library and HVG-library fields
are coherent; reused external broad/subclass categories and EX/IN/CA1-ProS flags
match the frozen source labels exactly. Annotation provenance explicitly says
the source-verified MapMyCells labels are reused.

| Representation at resolution 0.4 | Clusters | SNN undirected edges | Broad-label purity among high-support nuclei | Clusters meeting preliminary 20-nuclei/3-mice-per-genotype screen |
| --- | ---: | ---: | ---: | ---: |
| Unintegrated | 35 | 3,843,404 | 0.99514 | 27 / 35 |
| Harmony | 33 | 3,856,657 | 0.99447 | 26 / 33 |

The between-representation cluster ARI is 0.89145. Purity uses 98,773 externally
high-support broad predictions and is a coarse label agreement diagnostic; it
does not validate subtype, anatomy or preservation of genotype effects. The
replication screen is not a replacement for the eventual model eligibility
receipts, gene/library universe checks and influence review. Small/mixed clusters
remain visible and several lack adequate independent-mouse support.

The successful SNN and igraph objects occupy approximately 0.105 and 0.095 GiB
per branch. This confirms that the sparse conversion bypass repaired the actual
failing allocation without reducing the cohort or changing clustering settings.
UMAP calibration logs record 45 unintegrated and 10 Harmony smooth-kNN failures;
all resulting coordinates are finite. The known Harmony BLAS warning means
single-threaded operation and should remain in runtime provenance.

Evidence: `reviews/seurat_v1_0_5_paired_cellbender_umap_check.tsv`,
`seurat_v1_0_5_paired_cellbender_cluster_summary.tsv`,
`seurat_v1_0_5_paired_cellbender_cluster_support.tsv`, COMPLETE receipt and
native graph-size audits. Six original PDFs were rendered and inspected: axes,
legends and all twelve facet labels are positioned correctly, but the native
small raster points are very faint at normal page zoom. A separately authorized
three-page supplemental PDF uses the same complete exported coordinates with
larger points, WT blue/cKO orange and seeded shuffled drawing order:
`reviews/figures/seurat_v1_0_5_paired_cellbender_visible_coordinates.pdf`.
All three pages passed visual checks. The adjacent plotting script and provenance
JSON bind every input, drawing-order seed and output hash. Native artifacts and
the frozen scientific source were not changed or refitted.

### Donor-enriched provisional DG-like cluster

In this count-sensitivity stage, resolution-0.4 cluster 16 contains 2,860 nuclei
in both representations. Harmony assigns 2,847 (99.55%) the provisional
`037 DG Glut` prediction; 2,826 of these are subclass high-support. None is
`REF_CA1_PROS_HS`, while 2,829 are `REF_EX_HS`. Harmony membership is 326 WT and
2,534 cKO nuclei, including 2,020 from cKO_933 (70.63% of this cluster).
The unintegrated version has 2,850 DG predictions, 323 WT/2,537 cKO and 2,023
cKO_933 nuclei, also with zero high-support CA1/ProS flags.

This is an actual sample-enriched molecular/atlas assignment that persists after
integration. cKO_933 sequencing partition provenance remains unverified. The
labels alone do not establish physical contamination, anatomical origin or
lineage; genotype-driven reference mismatch and sampled-region composition remain
possible explanations. Original-count representations, marker support and
independent-mouse/model influence must be reviewed before a genotype claim.

All twelve mouse counts, including wt_1123's zero membership, are retained in
`reviews/seurat_v1_0_5_paired_cellbender_cluster16_mouse_counts_all12.tsv`;
the adjacent summary and reference-subclass cross-tabs preserve full predictions
and support flags. Cluster numbers are local to this stage and representation,
and must not be assumed to identify the same population in another branch/tool.

## Completed paired-CellBender neuronal subsets and integration dependence

All three source-defined subsets exactly retain their frozen external high-support
memberships. Each keeps all 3,000 requested subset features; PCA/Harmony have
30 finite dimensions and both UMAPs have two finite dimensions. Per-capture
corrected-library totals, local log aliases and merged unique-prefix cluster
metadata agree exactly. Nonmember rows retain missing subset clusters rather
than fabricated labels. The completed cohort remains fixed across representations.

| Subset | Nuclei | Unintegrated / Harmony clusters | Between-representation ARI | Preliminary replicate-screen clusters, unintegrated / Harmony |
| --- | ---: | ---: | ---: | ---: |
| EX | 60,849 | 25 / 23 | 0.63113 | 20 / 21 |
| IN | 8,903 | 17 / 19 | 0.87827 | 14 / 15 |
| CA1/ProS candidates | 27,672 | 12 / 11 | 0.38573 | 7 / 9 |

The CA1/ProS partition changes substantially with sample correction. Its
cluster-condition normalized mutual information is 0.26938 before Harmony and
0.00727 after Harmony; mutual information is 0.38453 versus 0.01015 nats.
Weighted dominant-genotype cluster fraction changes from 0.8414 to 0.5594,
and the number of nuclei in clusters containing at least 95% one genotype drops
from 13,288 to zero. EX condition NMI similarly changes from 0.0914 to 0.0133;
IN stays near 0.0006. These are descriptive association diagnostics, not
hypothesis tests or proof that biological effects rather than technical effects
were erased. They establish that the integration sensitivity matters for this
dataset and should remain visible in the final interpretation.

Evidence: `reviews/seurat_v1_0_5_paired_cellbender_subsets_check.tsv` and
`seurat_v1_0_5_paired_cellbender_subsets_condition_association.tsv`.
Original-count counterparts, independent-mouse pseudobulk on original counts,
marker support and conditioning on fixed external identities remain necessary
before interpreting condition programs. Coarse global atlas agreement alone
does not establish preservation of these finer condition-associated partitions.

The earlier Scanpy log contains tie-corrected Wilcoxon divide warnings during descriptive
marker discovery. Constant/zero-expression genes can have undefined rank-test
statistics; marker NaNs must remain visible or be explicitly audited. These
nucleus-level marker tests are not the genotype inference endpoint.

## Completed primary original-full global representation

The original-full COMPLETE receipt records `2026-10-13T10:15:49Z`.
Independent checks retain exactly 108,062 unique frozen nuclei, all twelve mouse
identities/genotypes, all 3,000 fixed features, finite 30D PCA/Harmony and finite
2D UMAPs in the exact metadata order. Each full original library agrees with
the accepted per-nucleus count metadata; the selected original total is
849,122,137 UMIs. External categories and high-support flags match their frozen
reused annotation source exactly.

At resolution 0.4 both representations have 35 clusters, with between-branch
ARI 0.87213. Weighted broad-label purity among 99,144 high-support nuclei is
0.99432 before Harmony and 0.99480 afterward. Preliminary replication screens
support 27 and 26 clusters, respectively. These remain agreement/eligibility
diagnostics rather than independent validation of genotype programs.

The provisional DG-like donor-enriched population persists on primary original
counts: unintegrated cluster 16 has 2,861 nuclei (2,848 `037 DG Glut` predictions),
and Harmony cluster 16 has 2,858 nuclei (2,842 DG predictions). Both have zero
high-support CA1/ProS flags and approximately 70.6% cKO_933 membership. This
pattern is therefore not specific to corrected CellBender counts, while the
anatomical/reference-identity caveats above still apply.

Evidence: `reviews/seurat_v1_0_5_original_full_global_integrity_check.tsv`,
`seurat_v1_0_5_original_full_global_cluster_summary.tsv` and
`seurat_v1_0_5_original_full_global_cluster_support.tsv` plus the native receipt,
coordinate, feature and count-conservation sidecars. The completed original-full
subsets, SCT sensitivity and completed mouse-level count inference are
reviewed in the following sections.

The primary global Harmony log reaches the prespecified ten-iteration maximum
without a convergence message. Installed `harmony::harmonize` checks convergence
after every iteration, including the tenth, and prints a message if satisfied;
the ten-iteration cap therefore precedes its reported convergence criterion in
this run. Paired-original global reached the same cap. These are finite
iteration-capped representations, not runtime errors or confirmed converged
solutions. The fixed matching parameter policy remains recorded, and integration
interpretation needs this limitation. Original-count pseudobulk on fixed external
identities is not computed from these corrected coordinates.

## Completed paired-original global and matched count sensitivity

The paired-original global representation retains the same exact 104,298 nuclei
as corrected paired analysis. All 3,000 requested genes, original full-library
totals and finite PCA/Harmony/UMAP dimensions pass the independent checks.
Its selected original total is 846,367,233 UMIs, versus 793,300,571 corrected
UMIs on these identical nuclei: a 6.27% count reduction. This is a count-source
difference, not an independent measurement of an ambient contamination fraction.

Both paired-original resolution-0.4 representations have 36 clusters; their
between-representation ARI is 0.86214. Raw-versus-corrected ARI on the fixed
paired cohort is 0.89013 for unintegrated clusters and 0.92710 for Harmony
clusters. Broad high-support atlas purity remains approximately 0.9942.
Preliminary replication screens support 27/36 and 26/36 clusters. Harmony
reaches the configured ten-iteration cap without convergence. The single native
singleton grouping event remains recorded and used the preserved native
Seurat postprocessing on the sparse SNN.

Evidence: `reviews/seurat_v1_0_5_paired_original_global_check.tsv`. Supplemental
coordinate-only three-page PDFs are separately labeled and visually checked:
`reviews/figures/seurat_v1_0_5_original_full_visible_coordinates.pdf` and
`reviews/figures/seurat_v1_0_5_paired_original_visible_coordinates.pdf`.
Each contains paired unintegrated/Harmony genotype plots plus all twelve mouse
facets. The new versioned plotting script preserves the earlier plot script;
adjacent JSON receipts bind exact source, drawing order and output hashes.
No embedding, cluster, source workflow or native run was modified by plotting.

## Completed primary original-full neuronal subsets

The completed primary subsets pass 219 read-only cohort/metadata/count/feature/
coordinate checks and six separately loaded sparse SNN graph checks. Every
graph has the exact nucleus axes, finite nonnegative symmetric weights and
agreement with its resource audit. The declared 3,000-gene panels, all twelve
mouse full/HVG library audits and atlas prediction/confidence fields agree
with their frozen sources. PCA/Harmony have 30 finite dimensions and both
UMAPs have two. Harmony converged after six, five and five iterations for
EX, IN and CA1/ProS respectively.

| Population | Nuclei | Unintegrated / Harmony clusters | Between-representation ARI | Preliminary replication screens |
| --- | ---: | ---: | ---: | ---: |
| EX | 60,985 | 26 / 24 | 0.61917 | 19/26 and 21/24 |
| IN | 8,911 | 18 / 20 | 0.87443 | 14/18 and 15/20 |
| CA1/ProS candidates | 27,703 | 13 / 11 | 0.44507 | 8/13 and 9/11 |

Primary CA1/ProS condition association changes from NMI 0.23176 unintegrated
to 0.00599 with Harmony; weighted majority-condition fraction changes from
0.81266 to 0.55734. Unintegrated clusters containing at least 95% one
condition hold 12,120 nuclei, versus zero after Harmony. EX NMI changes from
0.08796 to 0.01295; IN remains approximately 0.0009–0.0010. These are
descriptive sensitivities without inferential p-values. They do not determine
biological/technical cause or loss/preservation of differential expression.
The accepted atlas confidence fields remain original-count predictions,
and CA1/ProS remains a candidate molecular population rather than anatomical
CA1 validation.

Evidence: `reviews/seurat_v1_0_5_original_full_subsets_check.tsv`,
`reviews/seurat_v1_0_5_original_full_subsets_graph_check.tsv`,
`reviews/seurat_v1_0_5_original_full_subsets_summary.tsv` and
`reviews/seurat_v1_0_5_original_full_subsets_condition_association.tsv`.
The summary also compares primary clusters restricted to paired IDs against
CellBender clusters; that comparison includes representation fitting on
different source cohorts and must not be interpreted as denoising alone.
The paired-original versus paired-CellBender comparison below holds the
fitting cohort fixed.

The readable three-page coordinate-only primary CA1/ProS PDF is visually
checked, with exact 27,703-nucleus captions and all twelve mouse facets:
`reviews/figures/seurat_v1_0_5_original_full_CA1_PROS_visible_coordinates.pdf`.
Its adjacent receipt binds source and plotting hashes, reproducibly shuffled
order and PDF SHA256
`64a9bb8b03a5258bc832991cd496dc7710abea87334efcbbdc4bc21bc482a039`.

## Completed paired-original neuronal subsets

The paired-original subset stage completed at 11:06:13 UTC. Its exported
representations pass 219 read-only technical checks: exact frozen EX/IN/CA1-ProS
memberships, matched 3,000-gene feature panels, finite and correctly ordered
PCA/Harmony/UMAP coordinates, all twelve mouse library conservation audits,
matching cluster aliases/merged metadata, and missing subset labels for
nonmembers. Harmony converged after five, five and four iterations for
EX, IN and CA1/ProS respectively.

| Population | Nuclei | Unintegrated / Harmony clusters | Between-representation ARI | Raw–CellBender ARI, unintegrated / Harmony |
| --- | ---: | ---: | ---: | ---: |
| EX | 60,849 | 25 / 22 | 0.61842 | 0.75748 / 0.70470 |
| IN | 8,903 | 18 / 18 | 0.91822 | 0.93786 / 0.89348 |
| CA1/ProS candidates | 27,672 | 12 / 11 | 0.37409 | 0.63666 / 0.63780 |

The raw–CellBender comparisons use identical paired nuclei, feature panels and
representation settings. Their fine neuronal partitions are more sensitive
than the global partition. A high global ARI cannot establish robustness of
the within-CA1/ProS or excitatory subdivisions. Preliminary replication screens
pass 20/25 and 21/22 EX clusters, 14/18 and 15/18 IN clusters, and 6/12 and
9/11 CA1/ProS clusters; final gene/model eligibility remains separate.

Condition association in paired-original CA1/ProS falls from NMI 0.26298 to
0.01124 after Harmony; weighted majority-condition fraction falls from
0.83865 to 0.57835. There are 13,699 nuclei in unintegrated clusters with at
least 95% one condition, versus zero in Harmony clusters. EX NMI falls from
0.09268 to 0.01423, while IN remains approximately 0.0006–0.0007. These are
descriptive representation associations with no inferential p-values. They do
not identify a biological or technical cause, and cannot establish loss or
preservation of genotype differential expression. Unintegrated analysis and
mouse-level count inference remain necessary sensitivities.

Evidence: `reviews/seurat_v1_0_5_paired_original_subsets_check.tsv`,
`reviews/seurat_v1_0_5_paired_original_subsets_summary.tsv` and
`reviews/seurat_v1_0_5_paired_original_subsets_condition_association.tsv`.

## Completed original-full SCTransform sensitivity

The native SCT stage completed at 10:55:22 UTC. Read-only inspection of its
durable RDS passes 40 technical checks: valid Seurat/Assay5 schema, the exact
108,062 original-full nuclei, sparse original 3,000-gene RNA counts, preserved
849,122,137 full-source UMIs, and matching full/HVG count sums for every mouse.
All four stored reductions have finite values and agree with their exported
TSVs: 30 PCA dimensions, 30 Harmony dimensions and two dimensions per UMAP.
All four nearest-neighbor/SNN graph objects retain the full nucleus axes and
sparse storage. All six fresh SCT cluster columns agree with exported metadata.

The actual residual/PCA panel is the exact intersection of modeled features
across all twelve captures: 2,731 of the requested 3,000 genes. The other 269
genes lack a fitted requested residual in one or more captures; capture-level
modeled panels contain 2,834–2,943 requested genes. The coverage table agrees
with every capture panel, and no zero residuals were fabricated. SCT residuals
were used for PCA and deliberately deleted before graph construction. The
durable RNA `data` layer contains original-library log1p values, while the
stored PCA/loadings and downstream graphs/embeddings come from SCT residuals.
This object does not contain a durable SCT assay/model or corrected UMI assay.
Pseudobulk still reopens the original integer counts.

Seurat's standard `nCount_RNA` and `nFeature_RNA` in this representation object
describe its retained 3,000-HVG assay. They must not be reused as full-gene QC
measurements. Full original-count library totals are explicitly retained in
`original_nCount_RNA`; `representation_full_library_UMIs` contains full totals
for the chosen count source (original for this SCT stage). Full-gene detection
and accepted QC measurements remain in the authoritative shared input metadata.
This assay scope also applies to the matched log representation objects.

Thirteen completed-dataset QC sidecars now expose clearly named original
full-gene, selected-source full-gene and thin-assay measurements by exact
nucleus ID, sample, condition and population. The source/roster/hash-bound
index is `reviews/metadata/20261013_v1_0_5__QC_SIDECAR_INDEX.tsv`; the ten-field
semantics table is `reviews/metadata/qc_field_semantics.tsv`. All thirteen
compressed TSVs passed a second read/hash/cohort check and total about
14.4 MB. Publication used existing metadata only, without loading/reaggregating
counts or editing native objects, outputs or source.

The supplied call uses native SCT v2, minimum detection in five nuclei,
up to 5,000 training nuclei per capture, and no additional covariate regression
or corrected-UMI output. The installed Seurat 5.5.1 method's unspecified
residual clipping default is ±sqrt(capture nuclei/30), giving bounds about
±14.90 to ±18.65 on these captures. These are source/API-derived settings;
`reviews/seurat_v1_0_5_sct_default_clip_ranges.tsv` records them per mouse.

At resolution 0.4, SCT has 35 unintegrated and 31 Harmony clusters, with
between-representation ARI 0.87265. Same-representation log-versus-SCT ARI is
0.76999 unintegrated and 0.83955 Harmony. These comparisons combine
normalization and the modeled-feature intersection, so they do not isolate
a pure normalization effect. Broad high-support atlas purity remains about
0.9944 and 0.9943; preliminary replication screens pass 29/35 and 28/31
clusters. This is a descriptive representation comparison, not independent
validation of fine identities or a mouse-level genotype test. SCT Harmony
reaches the configured ten-iteration cap without convergence; finite completed
coordinates do not remove that limitation.

Evidence: `reviews/seurat_v1_0_5_sct_rds_check.tsv`,
`reviews/seurat_v1_0_5_sct_metadata_check.tsv`,
`reviews/seurat_v1_0_5_sct_cluster_summary.tsv` and
`reviews/seurat_v1_0_5_sct_cluster_support.tsv`. Independent comparator output
in `reviews/sct_normalization_feature_comparison_v1_0_5` agrees with the two
resolution-0.4 log-versus-SCT ARIs.

The separately produced three-page coordinate-only SCT PDF is visually checked
and includes exact cohort, method and 2,731/3,000-feature captions plus every
mouse facet: `reviews/figures/seurat_v1_0_5_original_full_SCT_visible_coordinates.pdf`.
Its adjacent provenance receipt binds inputs, reproducibly shuffled drawing
order, the versioned plotting script and PDF SHA256
`e242ffdf6f2bf91774c5793bb42e0cb92faf328b3ad0629e82acb7e9561d22d6`.

On the identical 108,062-nucleus primary roster (54,931 WT and 53,131 cKO),
the following resolution-0.4 associations describe condition/cluster mixing:

| Representation | Condition–cluster NMI | Weighted majority-condition fraction |
| --- | ---: | ---: |
| Log, unintegrated | 0.013385 | 0.576604 |
| Log, Harmony | 0.007698 | 0.550360 |
| SCT, unintegrated | 0.011357 | 0.561076 |
| SCT, Harmony | 0.006953 | 0.546918 |

The cohort majority-condition baseline is 0.508329. NMI uses arithmetic
normalization; the weighted fraction is the sum of each cluster's larger
condition count divided by all nuclei. Both measures are descriptive and
contain no inferential p-values. The reduced global condition association
after Harmony does not determine whether any change is biological or
technical, and does not establish loss or preservation of genotype
differential expression. Global mixing also does not establish stability
within neuronal subsets. The SCT comparison includes its different 2,731-gene
panel. The full TSV and adjacent source/hash receipt are
`reviews/seurat_v1_0_5_original_global_condition_association.tsv`.

## Completed sampled marker exploration: sources and provisional programmes

All three marker producers are complete. The 12-row source/sampling manifest
verifies exact available IDs and original labels, sample/condition identity,
actual sample-balanced selection counts, input paths and hashes. Paired-original
and paired-CellBender have identical selected IDs and frozen original-full
Harmony resolution-0.4 labels; these are distinct from fresh paired graph labels.

| Actual marker field | Original available / selected nuclei | Paired available / selected nuclei | Original source clusters |
| --- | ---: | ---: | ---: |
| log_harmony_r0p4 | 108,062 / 20,862 | 104,298 / 20,661 | 35 |
| EX_harmony_r0p4 | 60,985 / 21,543 | 60,849 / 21,478 | 24 |
| IN_harmony_r0p4 | 8,911 / 8,215 | 8,903 / 8,207 | 20 |
| CA1_PROS_CANDIDATE_harmony_r0p4 | 27,703 / 9,794 | 27,672 / 9,766 | 11 |

The original marker CSVs contain 85,801 global, 48,369 EX, 40,751 IN and
19,748 CA1/ProS rows. All twelve branch/field tables have unique cluster ×
stable-ID pairs, supported gene IDs, finite numerical fields, coherent detection
fractions, positive effects and exact Seurat full-32,285-gene Bonferroni adjusted
P values. This marker correction is separate from pseudobulk endpoint BH.
Markers use positive sampled-nucleus Wilcoxon exploration, minimum detection
0.1, log2FC threshold 0.25 and no adjusted-P gate on the complete export.
They do not provide independent-mouse WT/cKO inference. Top-seven displays rank
by detection-fraction difference, then effect size. Curated dotplot colours are
gene-wise scaled expression; dot sizes are detection fractions among selected
nuclei. The numerical sidecars retain the actual pooled selected-nucleus values.

Representative programme support is visible: IN cluster 3 has Vip/Prox1/Adarb2,
IN cluster 2 has Sst/Trpc6, and CA1/ProS cluster 1 has Col11a1/Stac/Epha10.
These are provisional molecular descriptions, rather than new definitive
cell-type or physical sublayer assignments. Native inhibitory and CA1 deep-
associated programme PDFs were rendered and visually inspected; labels, dot
legends and numerical sidecars are readable and coherent.

Not every fine candidate has equivalent identity support. IN cluster 20 has
31 nuclei and CA1/ProS cluster 11 has 99; neither has any mouse with 20 nuclei,
so neither can support a genotype model. Their marker profiles contain both
neuronal and myeloid-associated signals: IN20 Gad1/2 plus P2ry12/Inpp5d/Dock8,
and CA1-ProS11 Slc17a7/Syt1 plus P2ry12/Mertk. On the identical selected paired
IDs, P2ry12 detection remains about 77% and 63% after correction, while some
C1qa/Slc1a3 signals decrease. These observations do not determine physical
contamination, doublets or lineage. IN16 (139 nuclei) also requires provisional
interpretation: lower Gad1/2 detection and mixed reference-subclass predictions
include `045 OB-STR-CTX Inh IMN`; no anatomical origin is inferred. CA1/ProS
cluster 7 has 970 nuclei but its largest positive marker detection-fraction
difference is only 0.018. Numerous nucleus-level rank-test rows alone do not
establish a distinct layer identity or cause for that weak specificity.

Evidence: `reviews/seurat_20261013_v1_0_5_MARKER_SCOPE_AND_SOURCE_LABEL_MANIFEST.tsv`
and adjacent provenance, plus `reviews/seurat_completed_marker_output_review_v1`
(summary, original top-three rows, candidate mouse support, curated expression
and complete source hashes). Native marker tables, selected-nucleus manifests
and curated dotplot numerical values remain unchanged.

## Completed three-branch pseudobulk: status and original-full core effects

The original-full producer completed at 12:55:21 UTC, followed by paired
producers at 13:29:28 and 14:05:37 UTC. The native workflow is complete. All 571 source populations are retained, including
explicitly ineligible definitions. Final status/error/warning evidence is:

| Engine | Endpoint | Tested | Not tested |
| --- | --- | ---: | ---: |
| DESeq2 | PRIMARY | 301 | 270 |
| edgeR | PRIMARY | 301 | 270 |
| DESeq2 | KNOWN_PARTITION_ADJUSTED | 292 | 279 |
| DESeq2 | MAIN_PARTITION_ONLY | 264 | 307 |
| DESeq2 | LOO, four core totals × twelve mice | 48 | 0 |

All 2,332 published endpoints have zero error and captured warning counts.
Read-only inspection of every tested endpoint's actual MICE_ALL12 audit
confirms exact declared inclusion rules, all twelve distinct source mice,
≥20 nuclei and positive UMIs in included mice, ≥3 mice per genotype,
full-rank metadata designs, matching residual df ≥2, and valid gene/P-value/
adjusted-P/BH-hit denominator bounds. The completed native evidence is
`results/seurat/original_full/pseudobulk/PSEUDOBULK_EXECUTION_AUDIT.tsv` and
its per-endpoint audit tree. These are technical guard/status checks, not
scientific acceptance of every overlapping population.

Paired-original and paired-CellBender each complete 1,142 PRIMARY-only
endpoints: 301 tested/270 not tested per engine. The three-branch total is
2,410 tested and 2,206 not tested, with zero errors/captured warnings. Actual
all-twelve-mouse identities, inclusion masks, positive libraries, ≥20 nuclei,
replication bounds, full design ranks and exact residual df, and denominator
bounds pass for every tested endpoint. The final source-bound navigation is published after native exact R/Python
aggregation acceptance: `reviews/pseudobulk_endpoint_navigation_20261013_v1_0_5_v1`.
It preserves all 3,426 PRIMARY rows and 4,616 total endpoint statuses, with
zero failed technical guards or review flags. Untested P/FDR counts remain
NA; source indices, actual audit/result paths and all mouse eligibility
evidence are retained.

The following effect interpretation concentrates on original primary models
and DESeq2 partition/LOO sensitivities for the four main Atlas totals. Completed
paired PRIMARY count sensitivity is reviewed in the following section, with
source hashes and exact fixed-membership/gene-universe checks. The first
four R count matrices passed the parent's independent R/Python count comparison.
Twenty-eight additional read-only checks verify the full gene axes, exact
integer count columns and audited libraries for all twelve mice, positive
fitted size factors, and agreement between produced normalized tables and
the existing count/size-factor arithmetic.

All four primary populations use six WT and six cKO mice with residual
degrees of freedom ten. DESeq2 poscounts normalization with signed Wald
tests and edgeR TMM normalization with robust quasi-likelihood F tests
show high effect agreement, while their BH gene lists differ:

| Population | Fitted genes | DESeq2 non-NA adjusted P values | DESeq2 / edgeR BH hits | Shared hits | BH-list Jaccard | Effect Spearman |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| ALL | 20,049 | 19,247 | 642 / 308 | 300 | 0.46154 | 0.99941 |
| EX | 18,603 | 18,193 | 894 / 505 | 497 | 0.55100 | 0.99947 |
| IN | 14,149 | 14,143 | 6 / 2 | 2 | 0.33333 | 0.99787 |
| CA1/ProS candidates | 16,255 | 16,255 | 2,025 / 1,829 | 1,793 | 0.86997 | 0.99864 |

The corresponding DESeq2 non-NA P-value denominators are 20,025, 18,554,
14,143 and 16,255. edgeR has non-NA P/FDR values for all fitted genes in these
four populations. NAs remain unmodified; the review does not classify them
as nonsignificant tests or attribute each to a specific cause without a
stored diagnostic. Median absolute log2FC differences between primary
methods are only about 0.0021–0.0035. Every shared BH hit has the same effect
direction. These support stable estimated effects under the two statistical
methods, rather than identical significance lists or independent biological
validation.

The known/main partition and twelve leave-one-mouse-out sensitivity models
per population use DESeq2. edgeR was fitted for PRIMARY only. Across these
DESeq2 sensitivities, all estimable primary DESeq2 BH-hit effects preserve
direction. Their
median absolute magnitude change is at most 0.0594 log2 units for ALL,
0.0704 for EX, 0.0660 for IN and 0.0472 for CA1/ProS; corresponding maximum
90th-percentile changes are about 0.1684, 0.2041, 0.1399 and 0.1455.
This stability applies to the strong primary hits. All-gene effect rankings
change more: known/main partition Spearman values are approximately 0.7749
ALL, 0.7768 EX, 0.7366 IN and 0.8877 CA1/ProS. Weak or near-zero effects
must not be called uniformly stable from the strong-hit summary. LOO checks
use effects, not comparisons of BH hit counts, because LOO deliberately
uses a different independent-filtering policy and fewer mice.

Produced partition tables and numerical design-rank checks agree in all
twelve common-population/endpoint combinations. Primary uses twelve mice,
rank two and residual df ten. Known-partition adjustment excludes cKO_933,
uses eleven mice, rank three and residual df eight; lane 8 contains two cKO
and zero WT. Main-only uses six WT and three cKO, rank two and residual df
seven. The adjusted condition contrast is supported mainly by the shared
main partition. Near-identical adjusted/main effects are expected from this
support and do not constitute independent confirmations or a genotype
contrast within lane 8.

Satb2 stable ID `ENSMUSG00000038331` is verified in the produced count and
normalized tables, with 48 existing per-mouse rows across the four totals.
Gene-body RNA effects are positive in ALL/EX/CA1-ProS (primary DESeq2 log2FC
0.9592, 0.9457 and 1.1658); IN is near zero and not significant. CA1/ProS
Satb2 effects stay positive across partition models and every omission,
with LOO log2FC 1.0938–1.2316. As the prior study plan explicitly states,
genotype validation is settled. These gene-body UMIs do not assess exon-4
retention, protein or function, and do not reopen the genotype validation.
No transcript/allele mechanism is inferred from these count matrices.

The stable-ID marker table covers 49 distinct stable IDs for the requested
CA1/layer-program genes and Satb2 target, with overlaps counted once. Strong examples in
CA1/ProS are Fstl4 −0.9809, Kirrel3 +0.9780 and Ndst4 +0.7877 log2FC;
both methods support the same directions. Other curated genes are weak or
not BH significant. Heterogeneous associated-marker responses do not prove
an anatomical superficial/deep shift, cell loss or a coherent layer identity
change. The full marker rows retain NAs and the model-universe status.

Evidence: `reviews/pseudobulk_first4_atlas_model_comparison_v1.tsv`,
`reviews/pseudobulk_first4_atlas_sensitivity_v1.tsv`,
`reviews/pseudobulk_first4_atlas_marker_rows_v1.tsv`,
`reviews/pseudobulk_first4_atlas_technical_check_v1.tsv`,
`reviews/pseudobulk_first4_atlas_partition_design_check_v1.tsv` and
`reviews/pseudobulk_Satb2_existing_per_mouse_counts_v1.tsv`, with source/hash
receipts. The prior study-context statement is in
`/Users/mujahid357/Documents/Codex/2026-09-03/plugin-creator-users-mujahid357-codex-skills-3/outputs/SATB2_twelve_mouse_validation_20261010/ANALYSIS_PLAN.md:13`.

The publication-style four-panel effect PDF is rendered and visually checked:
`reviews/figures/pseudobulk_first4_atlas_DESeq2_edgeR_effect_comparison_v1.pdf`.
It plots every finite matched gene effect, marks shared BH hits, includes
gene denominators/correlations and the original-full primary contrast, and
binds exact source tables in its adjacent provenance receipt. PDF SHA256 is
`6e8d88037fcbdfa0a7af9d8e56e8453053d3cab49278b9457a499f22f76fd45c`.

The first ALL primary endpoint took 6.56 seconds from its model-settings
publication to the final endpoint audit, including fitting and exports.
Its normalized all-gene/all-twelve table is 6.13 MB logical and about
4.01 MB allocated at review after content-preserving APFS compression.
The completed source labels predict 571 populations per branch, matching
the actual native population index exactly. The three branches declare
4,616 endpoint publications including 48 LOO checks; replication/design
eligibility gives an upper bound of 1,507 successful DESeq2 normalized-table
exports before gene-universe or fit failures. Physical volume can change
under the parent's transparent compression; this single dense endpoint
does not calibrate every model's runtime or allocation.

## Completed paired count sensitivity and layer-associated expression

The final source-bound review separates original-full to paired-original
cohort restriction from paired-original to paired-CellBender correction.
All four Atlas populations retain exact original label memberships restricted
to each branch roster; paired raw/corrected memberships are identical. Frozen
original gene universes and existing full-gene/all-twelve-mouse integer tables
pass their guards. Native independent R/Python acceptance is complete for these
four totals. No count aggregation, new fit or native output edit was performed
by this review.

Cohort restriction alone barely changes primary estimates: Spearman values
are 0.99865–0.99997, and CA1/ProS DESeq2 hits change from 2,025 to 2,024.
On the identical paired roster, count correction changes ranking and
statistical support more substantially:

| Population | Raw→corrected DESeq2 / edgeR BH hits | Effect Spearman, DESeq2 / edgeR | Shared BH hits, DESeq2 / edgeR |
| --- | ---: | ---: | ---: |
| ALL | 644→652 / 308→290 | 0.95767 / 0.96010 | 598 / 278 |
| EX | 894→856 / 503→446 | 0.95690 / 0.95881 | 795 / 429 |
| IN | 6→4 / 2→1 | 0.88950 / 0.89489 | 4 / 1 |
| CA1/ProS candidates | 2,024→1,947 / 1,829→1,639 | 0.94752 / 0.94836 | 1,806 / 1,565 |

Every estimable raw paired BH-hit effect preserves its direction in the
corrected model, but this does not imply identical magnitude, uncertainty,
rank or significance. Median absolute changes among those hits are about
0.0072–0.0177 for ALL/EX/CA1, and 0.0642 DESeq2/0.2075 edgeR for the small
IN hit set. Some weaker effects and tails change more. IN Arhgap31 stays
negative and grows from −0.816 to −0.948 in edgeR while FDR changes from
0.0063 to 0.1099; Apool remains significant (−1.323 to −1.606). Four IN genes
are explicitly excluded after correction as `ALL_ZERO_AFTER_ENDPOINT_SELECTION`
within the frozen original universe; corrected IN has 14,145 modeled genes.
This does not establish that removed RNA was ambient or that a biological
population disappeared.

Adjusted-P availability must remain separate from gained/lost significance.
CA1/ProS corrected DESeq2 has 15,309 non-NA FDRs among 16,255 modeled genes;
47 raw BH hits have a corrected FDR of NA. EX has three such raw hits and
two corrected hits with raw FDR NA; ALL has eight corrected hits with raw
FDR NA. They are not tested nonsignificant counterparts. CA1 full-list
BH Jaccard is 0.83418 for DESeq2, versus 0.85269 when restricted to genes
with both FDRs available; edgeR is 0.82239. The supplement retains both
scope definitions and all finite matched effects.

Layer-associated gene responses remain heterogeneous and persist across
counts. CA1/ProS paired-original→corrected DESeq2 estimates are Ndst4
+0.7878→+0.7830, Col11a1 +0.5046→+0.4963, Asic2 −0.3633→−0.3675,
Pcdh7 −0.3261→−0.3445, Ak5 +0.3343→+0.3247, Fstl4 −0.9807→−1.0199
and Kirrel3 +0.9782→+0.9807; both primary engines support their directions
and BH significance. Thsd7b remains near zero and nonsignificant. The
original-full DESeq2 partition/LOO review supports the strong-hit directions,
with more variable weak/all-gene rankings as documented above. This is
expression evidence within a provisional molecular candidate population,
not a uniform superficial/deep programme, physical sublayer reassignment,
neuron loss or a mechanism separated from composition.

Satb2 gene-body expression stays positive (+1.1657→+1.1573 DESeq2 in paired
CA1/ProS), preserving the earlier signal. Genotype validation is settled;
these counts do not assess exon-4 retention, protein or function.

Evidence: `reviews/pseudobulk_first4_atlas_count_sensitivity_v1` includes
24 effect comparisons, all per-mouse library/cohort audits, stable-ID marker
rows, actual engine/contrast/realized-gene checks, explicit FDR-NA distinctions,
all-zero exclusions and input/output hashes. BH is within each native
population/method/endpoint; there is no correction across all 571 overlapping
population hypotheses or the sensitivity runs. Partition and LOO fits are
DESeq2 only; edgeR is PRIMARY only.

## Final interpretation limits and completed review scope

All thirteen Seurat representation datasets, all three marker producers,
three pseudobulk producers, native count/partition comparison, explicit QC
sidecars and the review supplements have completed their declared technical
checks. The native final receipt has no runtime failure or comparator review
flags; that status does not establish every candidate's biological identity.
The four Atlas totals are the main reported population summaries. Wider
cluster/subclass endpoints remain exploratory, and explicit zero/ineligible
populations are retained rather than assigned invented values.

Harmony reduces fine EX/CA1-ProS condition association substantially. Its
sample variable is nested within genotype, so mixing alone does not determine
whether biology or technical variation changed. The three iteration-capped
global representations retain their stated convergence limit; raw-count
mouse-level inference is not calculated from integrated expression.
SCT remains a model/normalization/feature-panel representation sensitivity.
Reused external annotations, mixed/weak marker candidates, nonuniform
layer-associated genes and all-cell composition prevent automatic physical
cell-type, layer or lineage conclusions. Positive gene-body Satb2 RNA does
not reopen the prior genotype validation.

## Appendix: runtime repair history and technical background

### Earlier rendering recovery and preserved attempts

The v1.0.1 attempt stopped in Scanpy original-full all-cell figure rendering
at 10:58:51 Vienna time on 13 October, before a Seurat rule ran. The traceback
identified pandas categorical `fillna('NA')` handling in `scatter_panel`, after
upstream graph/marker work. A separately registered rendering recovery reused
accepted shared-input caches; failed runtime files were preserved. Shared
preparation receipts remain bound to accepted cache/recovery provenance.

Preserved earlier native runs:

- v1.0.1 `b954c44b-60fe-4596-af9e-553bb1609f89`;
- v1.0.2 `cf1ce285-5c89-432c-8ff4-e55956505ac6`;
- v1.0.3 `c7695452-71cf-4f4b-8823-3010fe1fca7a`.

The following repairs used source/API inspection, import-only I/O validation
or synthetic technical fixtures. Scientific representations were fitted only
through the registered native recovery, and successful actual v1.0.5 results
are reviewed in the main text.

### CellBender import failure and validated repair

At 11:13:37 Vienna time on 13 October, v1.0.2
`seurat/paired_cellbender/log_global` failed on the first capture (`wt_1122`):
`An object with name droplet_latents/data does not exist in this group`.
The log/FAILED receipt show an H5 importer failure before normalization or any
Seurat fitting. Seurat's generic reader treated CellBender root latent groups
as matrices. There is no basis to interpret this as a biological failure.

The executed runtime records R 4.6.1, Seurat 5.5.1, SeuratObject 5.4.0,
Matrix 1.7.6, Harmony 2.0.5, sctransform 0.4.3, glmGamPoi 1.24.0,
uwot 0.2.5, igraph 2.3.3, DESeq2 1.52.0 and edgeR 4.10.3.
The stage-start available-memory check passed at 27.765 GiB against the 24 GiB
guard; this failure was unrelated to the memory guard.

The parent authorized an import-only source repair, preserving the frozen failed
run. The new explicit reader reconstructs only `/matrix` CSC arrays and checks
shape, unique gene/barcode axes, GEX feature types, pointers, indices and integer
nonnegative counts. It preserves source hashes and exact selected rosters.
Corrected source libraries must match frozen `cellbender_nCount_RNA`; original
libraries must match original `nCount_RNA`. hdf5r 1.3.12 is recorded as a dependency.

Actual I/O-only validation on the same `wt_1122` source files passed:

| Count source | Genes | Frozen nuclei | Imported UMIs | Exact per-nucleus library matches |
| --- | ---: | ---: | ---: | ---: |
| Original | 32,285 | 8,546 | 52,372,773 | 8,546 / 8,546 |
| Corrected paired | 32,285 | 8,184 | 48,291,948 | 8,184 / 8,184 |

Evidence: `reviews/original_import_io_check_wt_1122.tsv` and
`reviews/cellbender_import_io_check_wt_1122.tsv`. No normalization, model fit or
label remapping was performed during these tests. New source SHA256:
`6411768a18e5e35a7a1121c5b47b464266bc388662293fed45c1af2e1244ea00`.
This is import validation, not acceptance of a completed Seurat representation.

### Graph-conversion memory failure and validated sparse repair

The v1.0.3 FAILED receipt records `2026-10-13T09:27:12Z` (11:27:12 Vienna).
All twelve capture imports completed, followed by scaling/PCA and Harmony
convergence after eight iterations. The next log lines are nearest-neighbor/SNN
construction and an error while evaluating `[` argument `i`:
`vector memory limit of 64.0 Gb reached, see mem.maxVSize()`.
Harmony separately warned that BLAS core setting failed and it ran single-threaded;
this affects performance, not the selected nuisance variables or count roster.

Installed function-body inspection identifies the allocation precisely:
Seurat 5.5.1 `RunLeiden` calls igraph's weighted adjacency reader with
`mode='lower'`; igraph 2.3.3 `ensure_no_na` indexes `x[lower.tri(x)]` before
dispatching to its sparse adjacency implementation. A 104,298-square logical
index alone occupies approximately 40.52 GiB, and integer row/column/index
intermediates add further allocations. The 13.56 GiB sparse representation
planning estimate did not account for this API's dense index. It is not evidence
of a dense count matrix or a biological QC failure. The [current official igraph source](https://github.com/igraph/rigraph/blob/main/R/utils-assert-args.R)
contains a sparse-specific triangle check.

The parent authorized a source-only repair, preserving frozen v1.0.3 evidence.
The adapter extracts the lower sparse triangle with `Matrix::tril`, excludes
self edges, preserves edge weights and retains every nucleus including isolates.
It passes that igraph object to native Seurat `RunLeiden`, retaining modularity,
the seed, ten iterations and native cluster-size ordering, then calls native
`GroupSingletons` on the original sparse SNN. No package update, namespace
override, cohort/feature/PC reduction or raised vector limit is used.

The synthetic 100-vertex weighted fixture, including isolates, has an exactly
matching edge table and exactly matching native FindClusters memberships/factor
levels at resolutions 0.2, 0.4 and 0.6. A separate 104,298-vertex two-edge
synthetic conversion retains all vertices without an all-pairs dense index;
the resulting graph occupies 0.00855 GiB. These are technical synthetic checks,
not real-data clustering or evidence of biological stability.
Evidence: `reviews/sparse_leiden_adapter_synthetic_check.tsv`.
Syntax passed. Repaired source SHA256:
`ef8dfe20d015df31355a219813e8d34efa1685281b0168a35ffd22530d3b5230`.
The source also drops the now-unused dense scaled matrix before graphs, publishes
lean PCA/Harmony/count checkpoints before graphs and retains each sparse SNN plus
its edge/size audit before clustering. At this synthetic-check checkpoint, actual-data execution was still required.
The completed v1.0.5 representations reviewed above now demonstrate successful
actual sparse graph conversion and clustering with the full declared cohorts.

An additional synthetic pre-execution API check found that an empty-matrix
`SetAssayData` call does not delete an Assay5 layer: it fails the feature-overlap
guard. That cleanup call existed at the final checkpoint before this repair and
was moved earlier for the graph memory improvement. The parent authorized
replacing both cleanup sites with the official guarded NULL layer setter.
Under the exact v5 assay option, a 6-by-6 synthetic object passed deletion twice,
unchanged counts/data/nucleus IDs, validObject and named cluster/Idents metadata.
Evidence: `reviews/assay5_layer_cleanup_synthetic_check.tsv`.
The combined current source SHA256 is
`f7e807929862c108de163febf07e30fe76d7978989c8fb1e92945dc65cbfed65`.
The earlier sparse adapter source and frozen attempted runs remain preserved.
This additional check also used no real capture or scientific model fit.

### Actual v1.0.5 repaired pregraph checkpoint

The completed pregraph TSVs independently retain exactly 104,298 unique frozen
nucleus IDs, in the same order for both 30-dimensional reductions, with all
values finite. PCA retains every one of the fixed 3,000 requested genes. Every
one of the twelve count-conservation rows agrees exactly with its frozen
corrected-library sum and nucleus count. The selected full-gene corrected total
is 793,300,571 UMIs; the representation's HVG count total is 206,552,541 UMIs.
The successful lean RDS checkpoint confirms the guarded Assay5 cleanup passed
on the actual representation. Harmony converged after eight iterations.

Evidence: `reviews/seurat_v1_0_5_paired_cellbender_pregraph_check.tsv`, the active
stage's PCA/Harmony TSVs, `PCA_ACTUAL_FEATURES.tsv`,
`CAPTURE_COUNT_CONSERVATION.tsv` and `PRE_GRAPH_HVG_ONLY.rds`.
This establishes roster/count/feature integrity and finite fitted coordinates.

