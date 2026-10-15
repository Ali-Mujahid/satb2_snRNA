# Scanpy result review — observed native-run outputs

Current native run: `runs/20261013_v1_0_5/results`, registry
`9702cbef-00e6-4f64-aa37-e9b62357bd71`. **All twelve Scanpy datasets are complete and independently reviewed.**
All four populations completed in each of the original-full, paired-original
and paired-CellBender count branches.

For every completed dataset, independent review verified exact shared
nucleus/mouse/condition metadata, saved file and common-HVG SHA256 values,
saved integer count-layer preservation, finite PCA/Harmony/UMAP coordinates,
and symmetric, finite-positive connected clustering/UMAP graphs without
isolated nuclei. Balanced marker selection retains the declared per-mouse
quota policy. Across 24 exported marker tables and 39,800 rows, scores,
p-values, adjusted p-values and log-fold changes are finite. All nine actual
PDF pages per dataset were visually inspected, totaling 108 reviewed pages.
These checks establish artifact integrity. **Python Harmony actual iterations:
`NOT_RECORDED`; convergence outcome/history: `NOT_RECORDED`.** The
in-memory result history was not persisted, and verbose iteration logging
was disabled. Finite coordinates do not establish convergence or cap exhaustion.

Cluster counts below use the declared main resolution of 0.4.
`scanpy_stage_summary_v1_0_5.tsv` preserves counts, UMI totals, file/count/HVG
hashes, cached-fit status and the explicit missing Harmony history fields.

| Count branch | Population | Nuclei | Clusters, unintegrated/Harmony | Independent review |
|---|---|---:|---:|---|
|original_full|allcells|108,062|31/28|Passed; valid cached fit|
|original_full|Excitatory|60,985|20/19|Passed|
|original_full|Inhibitory|8,911|16/17|Passed|
|original_full|CA1_ProS|27,703|12/9|Passed|
|paired_original|allcells|104,298|28/30|Passed|
|paired_original|Excitatory|60,849|21/19|Passed|
|paired_original|Inhibitory|8,903|16/17|Passed|
|paired_original|CA1_ProS|27,672|12/9|Passed|
|paired_cellbender|allcells|104,298|28/27|Passed|
|paired_cellbender|Excitatory|60,849|20/19|Passed|
|paired_cellbender|Inhibitory|8,903|15/16|Passed|
|paired_cellbender|CA1_ProS|27,672|11/9|Passed; valid cached fit|

**All three independent Python count exports also completed and passed
review.** Every branch contains 77 populations with frozen original-full
labels and all twelve mouse columns/coverage rows. All 231 count-file hashes
match their receipts; coverage sums match population nucleus totals, and
actual ALL gene-by-mouse tables preserve all 32,285 stable gene IDs as
nonnegative integers. Each mouse's gene-sum exactly equals its source-count
conservation total. Branch UMI totals are 849,122,137 (original-full),
846,367,233 (paired-original) and 793,300,571 (paired-CellBender).
`python_pseudobulk_export_review_v1_0_5.json` records these checks.
R-versus-Python gene-by-mouse equality, adjusted model eligibility and final
genotype/cross-tool conclusions remain downstream checks; no model was fit
by this independent review.

Existing MapMyCells labels are common original-query annotation evidence.
The preserved preparation receipts verified exact original integer query
counts for all twelve captures; correction-specific remapping was not
performed. The broad CA1-ProS cohort includes ProS and does not establish
physical CA1 layers. Existing workbook programs and rare mixed-signature
QC/eligibility evidence are described below without hard anatomical labels
or new exclusions. The fixed primary CA1-ProS cluster-8 population retains
93 nuclei in all three count branches, but no mouse meets the twenty-nucleus
minimum, making genotype inference ineligible for that rare group.

No scientific fit or prior-result modification was executed for this review.
The following recovery/provenance sections are background. They do not change
the current completion and validation status above.

## Historical v1_0_1 preparation and rendering failure

Run: `b954c44b-60fe-4596-af9e-553bb1609f89`.
Local results: `runs/20261013_v1_0_1/results`.
Both shared preparation stages completed, but the first allcells representation
failed during PDF rendering and did not produce a completed checkpoint.

The completed original-count preparation contains 108,062 nuclei and 32,285 stable
gene IDs from all 12 captures. Its source receipts cover 849,122,137 UMIs. Each of
the 12 saved MapMyCells-query verification records reports
`VERIFIED_EXACT_ORIGINAL_INTEGER_COUNTS`, raw normalization, exact gene/nucleus
axes, and zero unequal sparse-count entries. Their nucleus totals sum to 108,062.
This supports reuse of the existing original-query annotation as common evidence;
no new MapMyCells mapping was executed.

Independent inspection of the saved HDF5 metadata found CSR encoding, shape
`108062 × 32285`, int32 count data and 275,131,404 nonzero entries. The saved
count file is 661,570,895 bytes. These are original count inputs; this stage does
not provide genotype differential-expression results.

| Fixed reference population | Nuclei | Mice | Common HVGs |
|---|---:|---:|---:|
| All cells |108,062|12|3,000|
| Excitatory, high-support reference flag |60,985|12|3,000|
| Inhibitory, high-support reference flag |8,911|12|3,000|
| CA1-ProS, high-support reference flag |27,703|12|3,000|

For every panel, the reviewer counted exactly 3,000 unique gene IDs and verified
the TSV's SHA256 against the completed-stage receipt. Each membership TSV has
the stated number of rows and nuclei from all 12 mice. CA1-ProS is a broad
reference population including ProS; 27,703 nuclei must not be described as the
previous marker-gated anatomical CA1 cohort or as anatomically validated CA1.

The paired CellBender preparation also completed: 104,298 nuclei, 32,285 genes,
and 793,300,571 corrected UMIs. Its fixed populations contain 60,849 excitatory,
8,903 inhibitory and 27,672 CA1-ProS nuclei. Independent inspection confirmed all
12 mice in each population and all four HVG TSVs byte-identical to their
original-full source panels. A matched paired-original UMI total is not yet
available here; these full-versus-paired totals alone do not isolate correction
from cohort restriction.

No API failure or traceback was present in the completed preparation logs.
The allcells Scanpy stage produced partial marker and selected-roster TSVs before
failing in `scatter_panel`: pandas rejected `fillna('NA')` on categorical labels.
Neither `scanpy.h5ad` nor `COMPLETE.json` was written. The partial tables are not
a completed analysis. Graphs, embeddings and full cluster metadata were not
available as durable artifacts for this reviewer to validate.

On the parent's instruction, the source renderer was fixed to convert labels to
object dtype before missing-value replacement. A tiny synthetic PDF test passed
for categorical labels with and without missing values. This tested rendering,
not scientific fitting. Source-level cached-preparation reuse also passed a
synthetic fixture check; the parent controls native execution/recovery. Existing
run outputs and prior source datasets were not changed.

Cross-tool consistency and biological WT/cKO conclusions remain unsupported by
these preparation outputs and the failed representation attempt. A later native
recovery run must supply completed representation artifacts before graph,
clustering and PDF quality can be reviewed.

Observed sources: `shared/original_full/COMPLETE.json`, `branch_counts.h5ad`,
the four `common_hvg_*.tsv` and `membership_*.tsv` files, and
`logs/prepare_original_full.log`.

## Native recovery v1_0_2: completed shared-cache reuse

Registry run `cf1ce285-5c89-432c-8ff4-e55956505ac6` writes to
`runs/20261013_v1_0_2/results`. Both original-full and paired CellBender shared
stages completed with new config-specific receipts. Each records the v1_0_1
source directory, `scientific_preparation_refit_performed=false` and
`raw_count_sources_reloaded=false`. The cohort sizes and all four panel hashes
remain the values reported above.

The reviewer independently recomputed each copied count-file SHA256 and current
config-file SHA256; they match the new receipts. Count-file SHA256 values are
`e1d33c0f5143789de53de56f4c2faab4e659663418a91c2f92077ec7204a7e3f`
(original full) and
`2730ad583eccfb37cc27f660d9ff24d22204eca8a0ecfdf7be1f6071be23ec22`
(paired CellBender). These new outputs establish successful validated copying,
not a repeated scientific fit. The first new allcells representation was running
when this recovery-stage entry was written; no completion claim is made here.

## Completed v1_0_2 representations reviewed

`original_full/allcells` completed in 327 seconds with 108,062 nuclei,
31 unintegrated clusters (minimum size 43) and 28 Harmony clusters (minimum size
44). The saved H5AD is 1,534,500,632 bytes. Independent read-back checks passed
its whole-file SHA256, exact authoritative nucleus roster and sample/condition
identities, finite 30-PC arrays and both finite 2D UMAP arrays. All four saved
connectivity graphs are symmetric, finite, nonnegative and connected, with zero
isolated nuclei. Their nonzero counts are 3,123,326 and 3,149,800 for clustering
(unintegrated and Harmony), and 4,957,942 and 4,983,702 for the corresponding
UMAP graphs. Streaming the saved integer count layer reproduced its receipt's
digest exactly.

The allcells descriptive-marker rosters contain 21,660 unintegrated and 23,058
Harmony nuclei. Each is a unique subset of the complete cohort, with equal
contribution from every represented mouse within its cluster and at most
100 nuclei per mouse/cluster. This verifies the bounded descriptive selection;
these cell-level marker tests do not provide WT/cKO inference. All nine PDF pages
were rendered and visually inspected. Titles, legends, mouse labels and axes are
readable, with no clipping. The unintegrated cKO_933 composition differs visibly;
Harmony mixing alone cannot establish the origin of that difference.

`paired_cellbender/CA1_ProS` completed with 27,672 nuclei, 11 unintegrated clusters
(minimum size 183) and nine Harmony clusters (minimum size 145). Independent
read-back checks also passed file hash, exact fixed paired membership and mouse/
genotype metadata, count-layer digest, finite 30-PC and 2D embeddings and four
symmetric connected graphs without isolated nuclei. All saved broad MapMyCells
labels are Excitatory. Marker rosters contain 3,833 and 7,272 nuclei respectively
and satisfy the same within-cluster mouse-balancing/cap rules. All nine PDF pages
were rendered and inspected, with readable layout.

The CA1-ProS unintegrated plots show substantial WT/cKO cluster-composition
separation, while sample-level Harmony largely smooths that separation. Because
mouse identity is nested within genotype, this visual change does not prove that
the removed structure is technical. Keep the unintegrated control and fixed
atlas-population mouse-level pseudobulk central to biological interpretation.
No anatomical CA1 sublayer assignment or cross-tool consistency conclusion
follows from these completed representations alone.

The parent reported the v1_0_2 workflow subsequently stopped in R's corrected-H5
import. The two completed Scanpy stages above remain completed artifacts. On the
parent's instruction a source-only completed-population cache was added so they
can be copied into the next native run without refitting. A tiny synthetic fixture
passed exact copy/new receipt/source-preservation checks and rejected changed
UMAP parameters and changed genotype metadata; actual reuse awaits native
execution. Existing run outputs were not edited.

## Native recovery v1_0_3: completed representation reuse

Registry run `c7695452-71cf-4f4b-8823-3010fe1fca7a` writes to
`runs/20261013_v1_0_3/results`. Both shared caches and the two completed Scanpy
population caches succeeded. The new `original_full/allcells` and
`paired_cellbender/CA1_ProS` receipts explicitly identify their v1_0_2 source
stages and record `representation_refit_performed=false`.

The reviewer independently recomputed all ten copied artifact hashes per
population and verified them against the fresh receipts, along with the current
config hashes and original completed-receipt hashes. These are byte-identical
representations; the graph/count/metadata/PDF checks above apply unchanged.
The caches did not repeat scientific fitting. Other population and cross-tool
stages remain pending when this entry was written.

The later v1_0_3 native log reports `vector memory limit of 64.0 Gb reached`
after nearest-neighbor/SNN work. The sibling methods reviewer traced the installed
Seurat/igraph conversion path to a dense lower-triangle mask for the large
sparse graph. This is separate from the earlier corrected-H5 import and Scanpy
rendering failures. The parent controls any source repair and native recovery;
the completed Scanpy artifacts above were not affected.

## Native v1_0_4 cancellation and v1_0_5 continuation

The parent intentionally canceled v1_0_4 through the native tool before the R
stage reached a known invalid Assay5 empty-matrix cleanup call. Snakemake records
SIGTERM. This was an intentional stop, not a new unexplained memory/API failure.
The parent replaced that cleanup with a tested official NULL-layer deletion and
launched v1_0_5 with both R repairs frozen.

The two Scanpy population caches completed in v1_0_5. Their full copied-artifact
hash dictionaries agree exactly with the independently verified v1_0_3 caches,
and each new config-file hash matches its fresh receipt. Both explicitly record
`representation_refit_performed=false`. The earlier fit/graph/PDF observations
therefore apply to the same artifacts. New scientific output stages were still
pending when this entry was written.

## Fresh original-full CA1-ProS representation in v1_0_5

The newly completed original-count CA1-ProS fit contains exactly 27,703 frozen
nuclei from all 12 mice, 12 unintegrated clusters (smallest 75 nuclei) and nine
Harmony clusters (smallest 93). Independent read-back passed exact membership/
mouse/genotype metadata, whole-file SHA256 and the streamed saved integer-count
layer digest. Both 30-PC arrays and both 2D UMAP arrays are finite. All four
connectivity graphs are symmetric and connected with no isolated nuclei; stored
edge counts are 781,826 and 800,636 for clustering and 1,184,268 and 1,222,740
for UMAP (unintegrated and Harmony respectively).

Balanced descriptive-marker rosters contain 4,868 and 7,176 nuclei and satisfy
the same uniqueness, subset, within-cluster mouse-balance and 100-cell cap checks.
The marker-panel expression table has 1,305 rows, finite mean expression and
detection fractions within [0,1] to 1e-12 floating-point tolerance. All nine PDF
pages were rendered and visually inspected; their titles, legends, mouse labels
and axes are readable without clipping.

The nuclei-weighted dominant-genotype fraction across original-count CA1-ProS
clusters is 83.5% unintegrated and 54.7% after sample Harmony. This descriptive
change mirrors the paired corrected-count observation. It does not establish a
biological effect or its technical origin: mouse identity is nested within
genotype, and sample integration can remove genotype-associated structure.

## Fresh paired CellBender excitatory representation in v1_0_5

The corrected-count excitatory fit completed in 201 seconds with exactly 60,849
frozen nuclei from all 12 mice, 20 unintegrated clusters (smallest 212) and
19 Harmony clusters (smallest 191). Independent read-back checks passed exact
membership/mouse/genotype metadata, whole-file hash and streamed integer-count
layer digest, finite 30-PC and 2D embedding arrays, and four symmetric connected
graphs without isolated nuclei. Clustering graphs store 1,698,090 and 1,733,356
values; UMAP graphs store 2,608,026 and 2,653,828 (unintegrated and Harmony).
All frozen broad MapMyCells labels are Excitatory.

Balanced descriptive-marker rosters contain 11,983 and 15,275 nuclei and satisfy
the established uniqueness, cohort subset, per-cluster mouse-balance and 100-cell
cap checks. The 2,755-row marker-panel expression table is finite with detection
fractions bounded to floating-point tolerance. All nine PDF pages were rendered
and visually inspected with readable layout. A distinct cKO_933-enriched
excitatory group persists in both views, including after Harmony. This observed
sample-composition difference calls for capture/marker review before genotype
interpretation; it does not indicate an artifact-integrity failure.

Read-only metadata localization of that corrected-count excitatory observation:
Harmony cluster 6 contains 2,781 nuclei, including 1,959 cKO_933 nuclei (70.4%).
Its frozen MapMyCells subclass labels are 2,780 `037 DG Glut` and one
`038 DG-PIR Ex IMN`. All 2,781 are outside `REF_CA1_PROS_HS`. This identifies the
reference-labeled broader DG/excitatory compartment affected by the conspicuous
sample-composition pattern; it does not establish independent anatomical
validation of those subclass labels.

## Fresh original-full inhibitory representation in v1_0_5

The original-count inhibitory fit completed in 64 seconds with exactly 8,911
frozen nuclei from all 12 mice, 16 unintegrated and 17 Harmony clusters; the
smallest cluster contains 30 nuclei in each view. Independent read-back checks
passed whole-file hash, authoritative membership/mouse/genotype metadata and
streamed saved integer-count-layer digest. The 30-PC and 2D embedding arrays are
finite. All four connectivity graphs are symmetric and connected, without
isolated nuclei; clustering graphs store 231,868 and 232,340 values and UMAP
graphs 361,544 and 362,290 (unintegrated and Harmony).

Both descriptive-marker rosters contain 5,321 unique cohort nuclei and satisfy
mouse balancing within each represented cluster, with maximum quotas 89 and
87 respectively. The marker-panel expression table has 2,465 finite rows with
detection fractions bounded to floating-point tolerance. All nine PDF pages
were rendered and visually inspected with readable layout. Frozen broad labels
are all Inhibitory; major reference subclasses include Sst (2,585 nuclei),
Pvalb (2,439), Vip (1,124), RHP-COA Ndnf (678), Lamp5 Lhx6 (606) and Sncg (500).
These are existing reference annotations rather than independently established
new anatomical identities.

The shared marker audit maps 145 of 147 eligible panel rows to the source gene
IDs. The only unique unmapped symbol is `AC129186.1`, present in two CA1-related
worksheet contexts. Its absence is explicit; no substitute symbol was invented.

## Qualitative identity support for original-full inhibitory clusters

The companion `original_full_IN_identity_support.tsv` joins current Harmony
cluster membership to frozen atlas subclass assignments, existing balanced
cluster-marker rankings and the existing workbook-panel expression table. It
does not fit a new model or assign hard anatomical identities. Major broad
inhibitory groups have substantial Gad1/Gad2 detection. Sst-aligned clusters 3,
8 and 9 have Sst among their ranked markers and corresponding panel expression;
Vip cluster 10 has Vip ranked first, 90.2% panel detection and 99.5% atlas
subclass agreement. Cluster 11 is also atlas Vip-aligned (99.2%) but has lower
Vip panel detection (22.6%), so the label's evidence differs between these groups.

Cluster 0 has 97.5% atlas Lamp5 agreement and Ndnf/Lamp5/Reln ranked markers.
Cluster 6 has 100% Lamp5-Lhx6 atlas agreement and Lamp5/Gad2 marker support.
Cluster 2 has 99.0% Pvalb atlas agreement; Pvalb is not among its top 100 ranked
markers, which does not demonstrate its absence. Cluster 5 has 83.9% Pvalb
atlas agreement and Pvalb in the ranked markers. Cluster 12 is an atlas
Pvalb-chandelier-associated candidate (99.2%) with Unc5b/Vipr2 among its leading
markers; this does not establish morphology or anatomical position.

Some groups need cautious mixed/uncertain interpretation. Cluster 1 contains
52.0% Sncg and 45.5% RHP-COA Ndnf atlas assignments, with Cnr1 among its leading
markers. Cck appears as a positive descriptive marker in this cluster in both
existing outputs: unintegrated rank 89, Harmony rank 99. This supports a
Cck/Cnr1-associated inhibitory group, not a hard anatomical subtype or genotype
effect. Cluster 4 is almost evenly split between Pvalb (51.6%) and Sst (48.4%)
atlas calls; graph-level mixture does not prove single-nucleus coexpression.

Rare cluster 7 has only 65 nuclei and leading Gm/uncharacterized-gene markers;
cluster 14 has 126 nuclei from six mice with mixed LSX-related atlas labels;
cluster 15 has 30 nuclei with diffuse atlas subclass assignments. Clusters 14
and 15 have no eligible balanced marker output because their selected group
sizes do not meet the declared marker threshold. Cluster 16 has 82 nuclei,
Top2a/Melk/Dtl-associated leading markers and predominantly atlas IMN labels;
its broad Gad1/Gad2 detection is weaker than in the major groups. These should
remain rare/mixed/state-associated candidates, without new exclusions or hard
CA1 interneuron labels.

The current eligible workbook panel contains Gad1, Gad2, Slc32a1, Sst and Vip,
but does not contain Pvalb, Reln, Ndnf, Lamp5, Sncg or Cck entries. Evidence for
those latter genes here comes from the already-produced ranked marker tables
and frozen atlas, not a panel expression measurement. The companion table
records dominant/secondary atlas fractions and source marker strings per cluster.

## Harmony convergence evidence and disk planning

All twelve completed Scanpy representations use the configured Harmony
maximum of ten iterations. Their actual iteration counts and
convergence status are **not recorded**: the source calls harmonypy with
`verbose=false`, logs have no iteration/convergence lines, and H5AD `uns` has
no saved objective or iteration-history fields. The in-memory harmonypy result
exposes `objective_harmony` and `kmeans_rounds`, but these were not persisted
before the result was deleted. Finite corrected coordinates verify numerical
validity, not convergence. Do not infer either convergence or cap exhaustion for
Scanpy from these artifacts or from the separately observed R convergence logs.
Settings, active runs and scientific source were left unchanged for this review.

The historical Scanpy disk estimate, made when five stages were complete and
seven remained, gave 5.44 GiB of additional H5ADs at observed analog
sizes, or 7.08 GiB with a 30% margin. Remaining metadata/marker TSVs and PDFs
add approximately 0.5 GiB. A conservative uncompressed upper allowance for
Python pseudobulk count/membership tables is 1.5 GiB across three branches;
gzip output is expected to be smaller. Plan roughly 9-10 GiB for remaining
Scanpy outputs, excluding R/SCT outputs, swap and any optional extra recovery
copy. All three shared count files are already present. The task-created
`reviews/tmp/pdfs` render scratch occupies approximately 11 MiB and is the only
small safe scratch candidate identified here; nothing was deleted. Prior
analyses, failed-run provenance and the installed runtime remain preserved.

## Descriptive ranked-marker statistic audit

`scanpy_marker_statistic_audit.tsv` covers all twenty-four ranked-marker tables
from the twelve completed Scanpy representations. Across 39,800 exported rows, there are zero
nonfinite scores, p-values, adjusted p-values or log-fold changes. No exported
row has zero detected expression in both the tested group and its reference.
These are descriptive cluster-marker statistics, not mouse-level genotype tests.

The paired-original CA1-ProS marker-selection pools contain 4,448 unintegrated
and 8,280 Harmony nuclei. Direct streaming of the saved integer counts found
8,251 and 7,460 entirely zero genes respectively among the preserved 32,285
gene columns. Installed Scanpy's tie-corrected Wilcoxon implementation divides
the centered rank sum by its corrected standard deviation; constant ranks can
have zero standard deviation and produce the logged invalid-divide warning.
Immediately afterward the library replaces NaN scores with zero and computes
a two-sided p-value of one for that score. This mechanism is consistent with
the observed all-zero genes; the warning does not identify particular genes.
Such uninformative genes are not present in the actual first-100 marker exports
examined here.

The source writer does not apply a blanket drop of rows with NA statistics.
Genuine NaN TSV fields would be written as `NA`; none were found in these exported
statistic columns. Groups with fewer than twenty selected balanced nuclei are
excluded from marker testing by the declared group threshold. All their nuclei
remain in embeddings, clustering and applicable count aggregation. Marker
warning handling is separate from saved count and embedding integrity.

## Completed paired-original CA1-ProS integrity

The paired-original CA1-ProS stage completed with 27,672 nuclei, 32,285 genes,
12 unintegrated and nine Harmony clusters. Independent checks found exact
membership, sample/condition metadata, receipt file SHA and saved integer count
layer digest, along with finite 30-PC and 2-D coordinates. Both clustering
connectivity graphs and both UMAP connectivity graphs are symmetric, finite,
positive, connected and have no isolated nuclei. Their nonzero counts are
780,840/799,874 and 1,182,010/1,221,298 respectively. The balanced descriptive
marker rosters contain 4,448/8,280 nuclei and obey the per-sample cluster cap.
The curated expression output has 1,305 rows with finite summary values. All
nine actual PDF pages were visually inspected and have readable plots/legends
without clipping or blank rendering. This establishes artifact integrity; it
does not establish agreement with Seurat or genotype effects.

## CA1-ProS identity and workbook expression programs

`CA1_ProS_program_identity_support.tsv` joins the already-produced balanced
positive ranked markers, frozen atlas supertype assignments and workbook-panel
summaries for all nine Harmony clusters in each of the three completed
CA1-ProS count-source branches. It is descriptive evidence from existing
outputs, with no new scoring, fitting, exclusion or anatomical relabeling.
The retained population is the broad atlas CA1-ProS glutamatergic population,
which includes ProS; it is not a physical CA1 layer assignment.

The original-full clusters separate coherent workbook programs. Harmony
cluster 1 has Ndst4/Col11a1/Epha10 among its leading markers and thirteen
positive deep-VEH panel genes among its first 100 markers. Cluster 3 has
Thsd7b/Pcdh7/Lingo2/Ak5 and nine superficial-VEH panel genes, while cluster 7
has C1qtnf7/Pparg-associated leading markers and eight superficial-VEH panel
genes. Cluster 2 has strong CA1 candidate marker support (including
4921539H07Rik/Chrm3/Epha3-related markers) and coherent atlas support, but few
of the layer-associated panel markers; its layer program remains unspecified.
Some other groups combine superficial/deep panel genes, state-associated
Pcsk1/Bdnf/Nnat markers or more mixed reference assignments. These are
expression programs associated with the workbook headings, not established
physical superficial/deep positions. Satb2 alone was not used for identity.

The broad Ndst4/Col11a1/Epha10 deep-associated and
Thsd7b/Pcdh7/Lingo2 superficial-associated programs recur in both paired
branches. The C1qtnf7/Pparg-associated group retains eight superficial-VEH
panel hits in all three. Numeric cluster IDs and memberships change between
independent representations, so equal numbers are not a matched-population
claim and rank changes alone do not isolate a CellBender effect.

Cck is a positive descriptive marker in the Tenm3/Ntm/Atp2b1-associated
excitatory group in each existing CA1-ProS output: original-full cluster 6
(rank 54, log-fold change 0.893), paired-original cluster 2 (rank 51, 0.902)
and paired-CellBender cluster 7 (rank 30, 1.135). This supports excitatory Cck
expression and cautions against defining inhibitory identity with Cck alone.
The quantitative matched fixed-membership Cck count analysis is separate.
The eligible workbook-panel output does not include Cck, Satb2, Calb1,
Calb2, Pcp4 or Wfs1 entries, so this review does not claim their panel
measurements. Slc17a7 has substantial detection in most of these clusters.

All three completed CA1-ProS representations also contain a rare immune
marker mixture (their Harmony cluster 8) with leading genes such as
Fli1/Ikzf1/Runx1/Inpp5d/Siglech/Dock8/Itgam rather than the workbook CA1 or
layer program. It remains a provisional mixed-signature group requiring
review. The graph-level marker mixture does not establish single-nucleus
coexpression, a doublet call or an immune anatomical identity.

## Primary CA1-ProS mixed-signature cluster 8: mouse and frozen QC evidence

The primary original-full group has 93 nuclei from all twelve mice: 43 WT and
50 cKO. Counts per mouse are WT 8/6/12/5/7/5 and cKO 12/3/11/6/3/15 in the
sample order in `primary_CA1_ProS_cluster8_mouse_coverage.tsv`. Every mouse
has fewer than twenty nuclei. Thus zero WT and zero cKO mice meet the
predeclared per-mouse minimum, and this rare group is ineligible for the
planned genotype model. This conclusion is independent of its tentative
biological identity.

Direct inspection of the actual primary TSV confirms carried, source-audited
scDblFinder aliases `score_new` and `call_new`, plus the neutral review flags
and strict-low sensitivity field. All 93 existing `call_new` values are
`singlet`; all pass primary QC and strict-low sensitivity and have
`primary_reason=retained`. The existing score median is 0.3768 (range
0.00070–0.84790). Two nuclei have `flag_high_genes_review=True`:
`wt_1178_CACCTATTCGATGCGT-1` (410,739 UMIs, 14,595 detected genes, score
0.84790) and `wt_1226_TCTGTCCGTTCTGTCC-1` (54,279 UMIs, 8,392 genes, score
0.72797). All ribosomal and robust-low review flags are false. Expected
doublet fractions are capture-level expectations, not nucleus calls. An
earlier message incorrectly claimed scores/calls were absent because it
searched only conventional column names; the direct alias audit corrects
that claim, and the companion detail TSV now preserves these exact fields.

All 93 retain high-support frozen CA1-ProS atlas subclass calls. Subclass
bootstrapping probability has median 1.0 and minimum 0.84; aggregate
probability has median 1.0 and minimum 0.68. Fine supertype/cluster agreement
is more variable. Existing singlet calls and strong broad atlas support do
not resolve the immune marker mixture or establish anatomical origin. No
nucleus was excluded, no existing call was changed and no new QC model was
fit for this review.


## Additional completed stages: light review during SCT preparation

Paired-CellBender/Inhibitory completed with 8,903 nuclei, 32,285 genes,
15 unintegrated and 16 Harmony clusters; minimum cluster sizes are 98/80.
Paired-original/allcells completed with 104,298 nuclei, 32,285 genes,
28 unintegrated and 30 Harmony clusters; minimum sizes are 43/43. Both
metadata exports match the exact shared population nucleus IDs, mouse and
condition values with zero mismatches; all twelve mice remain represented
and exported UMAP values are finite. All nine pages of each actual PDF were
rendered and visually inspected: plots and legends are readable without
blank or clipped content. The new marker tables add 7,800 finite rows to the
statistic audit. Receipt-level count preservation is recorded, but independent
full H5AD hash/count-digest/graph checks for these two are intentionally pending
until the concurrent SCT preparation releases memory. No full-count load was
performed for this light review.

Existing metadata also verifies that the corresponding immune-marker groups
are not identical memberships across representations: original-full versus
paired-original intersects at 79 of their 93/93 nuclei (Jaccard 0.738), and
original-full versus paired-CellBender intersects at 75 of 93/145 (0.460).
The paired-original group has no mouse with at least twenty nuclei and the
paired-CellBender group has only two such mice in total. Neither can reach
three eligible mice per genotype. All their frozen calls remain singlet.
These descriptive overlaps reinforce the need to compare biological programs
and actual membership, rather than assume matching numeric cluster labels.


## Deferred H5AD checks completed after SCT checkpoint

`scanpy_deferred_H5AD_checks_v1_0_5.json` records independent bounded-streaming
checks for paired-CellBender/Inhibitory and paired-original/allcells. Both
saved file and common-HVG SHA256 values match their completed receipts;
saved H5AD nucleus order matches the exported metadata. PCA/Harmony and
unintegrated/Harmony UMAP arrays are finite with the expected 30/2 dimensions.
Saved integer count-layer digests match the stage receipts, and streamed UMI
totals exactly match the branch-specific metadata totals: 79,420,310 for
paired-CellBender/Inhibitory and 846,367,233 for paired-original/allcells.

For paired-CellBender/Inhibitory, clustering graphs have 232,394/232,206
nonzero entries and UMAP graphs 359,112/359,220. For paired-original/allcells,
the corresponding values are 3,004,440/3,030,448 and 4,759,270/4,785,086.
All eight graphs are symmetric, finite-positive and connected, with no
isolated nuclei. These checks used bounded HDF5 chunks and one sparse graph
at a time, without loading a full count matrix or fitting any scientific
model. The earlier temporary deferral is resolved. All eight completed
Scanpy representation outputs currently pass independent artifact and
visual review; actual Harmony iteration/convergence history remains unrecorded.


## Original-full Excitatory completed and reviewed

The remaining original-full representation completed with 60,985 nuclei,
32,285 genes, 20 unintegrated and nineteen Harmony clusters. Saved file and
common-HVG hashes, the integer count digest, exact shared population/mouse/
condition values and saved H5AD nucleus order all pass independent checks.
Streamed count total is 674,026,171 UMIs, exactly matching branch metadata.
All PCA/Harmony and UMAP arrays are finite. Clustering graphs have
1,701,910/1,734,324 nonzero entries and UMAP graphs 2,646,344/2,687,390;
all are symmetric, finite-positive and connected without isolated nuclei.

Balanced marker rosters contain 12,351/16,115 nuclei with equal quotas for
represented mice within each cluster, capped at 100. All 3,600 exported
marker rows have finite statistics and valid p-value/detection-fraction
ranges. The curated expression panel has 2,755 finite summary rows. All nine
actual PDF pages were visually inspected and are readable without blank or
clipped content. This brings completed and independently reviewed Scanpy
representations to nine of twelve, including all original-full populations.

A visibly distinct capture-composition feature persists in these primary
counts: Harmony cluster 10 contains 2,793 nuclei, including 1,971 from
cKO_933 (70.6%). Frozen atlas subclasses are 2,784 DG Glut and nine DG-PIR
Ex IMN; all lack the CA1-ProS high-support flag. This supports a broader
excitatory capture-composition observation rather than a CA1 identity or a
replicated genotype effect. Numeric IDs differ from the corresponding
paired-CellBender representation and are not assumed to be matched populations.
The read-only reviewer helper in `reviews/audit_completed_scanpy_h5ad.py`
uses bounded HDF5 chunks and one sparse graph at a time; it performs no
scientific fit and is outside the frozen scientific workflow source.


## Paired-CellBender allcells completed and reviewed

The corrected matched allcells representation completed with 104,298 nuclei,
32,285 genes, 28 unintegrated and 27 Harmony clusters. Exact shared population/
mouse/condition metadata, saved file and common-HVG SHA256, saved H5AD nucleus
order and integer count-layer digest pass independent review. The preserved
corrected count total is 793,300,571 UMIs, exactly matching branch metadata.
PCA/Harmony and UMAP arrays are finite with the expected 30/2 dimensions.
Clustering graphs contain 2,998,756/3,024,208 nonzero entries and UMAP graphs
4,675,676/4,692,734; all are symmetric, finite-positive and connected without
isolated nuclei.

Balanced marker rosters contain 22,241/22,163 nuclei with equal quotas within
clusters for represented mice and a cap of 100. All 4,900 exported marker rows
have finite statistics and valid p-value/detection fractions; all 3,915 curated
panel rows have finite summaries. All nine actual PDF pages were inspected and
render correctly without blank or clipped content. This brings completed and
independently reviewed Scanpy representations to ten of twelve. Paired-original
Excitatory and Inhibitory are the remaining Scanpy stages. This corrected matrix
uses frozen original-query annotation evidence; no fresh correction-specific
MapMyCells remapping was performed or claimed.


## Paired-original Excitatory completed and reviewed

The matched original-count excitatory representation completed with 60,849
nuclei, 32,285 genes, 21 unintegrated and nineteen Harmony clusters. Exact
shared nucleus/mouse/condition metadata, saved file/common-HVG SHA256,
saved H5AD nucleus order and integer count digest all pass independent review.
The preserved total is 673,879,971 UMIs, matching branch metadata. All 30-PC
and 2-D coordinate arrays are finite. Clustering graphs have
1,697,862/1,730,414 nonzero entries and UMAP graphs 2,637,990/2,679,842;
all are symmetric, finite-positive and connected, with no isolated nuclei.

Balanced marker rosters contain 13,301/14,603 nuclei with equal quotas among
represented mice within clusters and a cap of 100. All 3,800 ranked-marker
rows have finite statistics and valid p-value/detection fraction ranges;
all 2,755 panel summary rows are finite. All nine actual PDF pages were
visually inspected and render correctly. Eleven of twelve Scanpy
representations are now complete and reviewed; paired-original/Inhibitory
is the remaining stage. No scientific fit was executed for this review.

## Independent Python count exports completed and reviewed

All three native integer-aggregation exports completed with 77 populations per
branch. Each preserves the frozen original-full label source: four common atlas
totals plus 28 allcells, nineteen excitatory, seventeen inhibitory and nine
CA1-ProS Harmony clusters. Independent review checked all 231 population count
file hashes against producer receipts, all-twelve coverage for every population,
and exact agreement between coverage nucleus sums and population-index totals.
The actual ALL tables contain all 32,285 unique stable gene IDs as nonnegative
integer counts with the specified twelve mouse columns. Each mouse's gene-sum
exactly matches its recorded source conservation total. Branch totals are
849,122,137, 846,367,233 and 793,300,571 UMIs for original-full, paired-original
and paired-CellBender respectively.

Evidence is saved in `python_pseudobulk_export_review_v1_0_5.json`. These are
count-aggregation/export checks, not fitted genotype models or proof of R/Python
gene-by-mouse equality. The export eligibility ledger explicitly screens an
unadjusted condition design; final adjusted model eligibility and inference
remain the downstream R statistical module's responsibility. No cell-level
genotype test or new scientific fit was performed for this review.

For the fixed original-full CA1-ProS Harmony cluster-8 count population, the
actual export ledgers preserve the same 93 nuclei in all three branches. Each
has zero eligible WT and zero eligible cKO mice and explicitly records
`FEWER_THAN_THREE_ELIGIBLE_MICE_PER_GENOTYPE`. This frozen population is distinct
from the independently reclustered paired representations whose cluster-8
memberships and sizes differ. No genotype model is supported for this frozen
rare group under the declared thresholds.


## Final paired-original Inhibitory completion and full Scanpy review

The last representation completed in 67.53 seconds with 8,903 nuclei from
all twelve mice, 32,285 genes, sixteen unintegrated and seventeen Harmony
clusters at resolution 0.4. Independent bounded-memory read-back verified
the exact shared nucleus/mouse/condition identities, saved H5AD nucleus
order, whole-file SHA256, common-HVG SHA256 and saved integer count-layer
digest. Its preserved 83,164,878 UMIs exactly match the branch metadata.
Both 30-PC arrays and both 2-D UMAP arrays are finite.

Clustering graphs have 231,644/231,804 nonzero entries and UMAP graphs
361,122/361,162 (unintegrated/Harmony). All four are symmetric,
finite-positive and connected without isolated nuclei. The balanced marker
rosters contain 5,309/5,225 nuclei, with equal quotas among represented mice
within each cluster and maximum quotas 87/64. Its 2,900 exported ranked-marker
rows have finite scores, p-values, adjusted p-values and log-fold changes;
no row has zero detection in both the group and reference. The 2,465-row
curated expression panel has finite, bounded summaries. All nine actual
PDF pages were visually inspected and have readable plots, legends and
mouse labels without blank or clipped content. No new unexpected artifact
integrity or rendering concern was found.

This completes independent artifact and visual review of all twelve
Scanpy representations and 108 PDF pages. The full marker audit covers
24 tables and 39,800 rows without nonfinite exported statistics. All three
Python aggregation exports also passed the count-conservation/receipt
checks described above. The separate final cross-tool comparator and
mouse-level statistical review establish any tool agreement and genotype
inference; artifact validity alone does not establish biological agreement.
Python Harmony actual iterations and convergence history remain
`NOT_RECORDED` for every representation. No scientific fit or source change
was performed for this final review.
