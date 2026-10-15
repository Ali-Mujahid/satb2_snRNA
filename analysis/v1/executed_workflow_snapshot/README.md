# SATB2 all-twelve Scanpy / Seurat v5 comparison 1.0.0

Native Snakemake implementation of the user-selected matrix reanalysis. Inputs
are the previously audited original integer counts/QC roster, the reviewed
12-October CellBender archive, and same-query MapMyCells annotations. The source
tree contains code only; inputs, environments, credentials, configurations and
scientific results are outside it. Original data and previous analyses are read
only. The workflow uses a new execution directory and refuses stage overwrites.

The primary cohort is all 108,062 original-QC nuclei across six WT and six cKO
mice. Two matched sensitivities compare original and CellBender counts on
104,298 supported nuclei. No zero filling, mixed count sources, extra CellBender
barcodes, sample exclusion or label forcing is permitted. Correction diagnostics
and Cck removal flags are preserved in the separate preparation audit.

Both tools use common population-specific 3,000 features, full-library log1p
normalization to 10,000, 30 PCs, graph k=20 and resolutions 0.2/0.4/0.6 (main0.4).
Harmony2 uses capture sample, theta2/sigma0.1/lambda1/100 clusters/10 iterations;
unintegrated controls are retained. Separate UMAP graphs use k30 and cosine.
Seurat's SNN/Leiden and Scanpy's graph/Leiden implementations can differ even
under these settings; software agreement is measured, not forced. Seurat SCT v2
is a separately labeled normalization sensitivity with a documented modelable
feature intersection.

Fresh global, EX, IN and provisional CA1-ProS clustering is followed by
sample-balanced exploratory marker tests. Marker subsampling never changes the
full clustering or mouse-level pseudobulk cohorts. The supplied workbook panels
are curated expression checks. Molecular layer programs do not prove anatomy.

WT/cKO expression inference uses animal-level integer pseudobulks with DESeq2
and edgeR; primary ~condition, WT baseline. Sequencing-partition sensitivities
retain the unresolved cKO_933 assignment explicitly. Shared atlas totals and
tool-specific fresh clusters are separate endpoints. Paired count comparisons
hold the original membership and gene universe fixed. Python aggregation is
checked against R on identical memberships independently of the DE model.

All scientific tasks use four threads under a four-core controller, making the
memory-intensive stages sequential. No environment installation, network call,
cloud upload, deletion of prior data, or unrelated job management occurs in the
workflow. Python dependency setup is completed separately in a new isolated venv.
Actual peak memory and model/figure acceptance require review of execution.
