# Task 4: scVelo input assessment, 28 October 2026

Status: **INPUT_BLOCKED / NOT EXECUTED**. This is a read-only inventory, not a velocity result. No biological fit, package installation, transfer, cluster submission or original-source modification was performed. The accompanying `inventory_evidence.json` records inspected headers and distinguishes current observations from historical records.

## Objective and biological design

Assess whether genuine spliced/unspliced RNA abundances can support deterministic, stochastic and dynamical scVelo analyses and conditional latent-time inference. The existing study consists of six WT and six Satb2 cKO male P60 mice, one separate 10x GEM-X 3-prime v4 capture per mouse, sampled once. These design statements come from accepted prior metadata; they were not independently re-established from raw reads here. Mice remain the biological replicates. CA1 deep/superficial labels denote molecular expression programs, not confirmed anatomical position or successive developmental stages.

## Current observed material

All 12 original filtered gene matrices exist under `/Users/mujahid357/Desktop/data/cellranger.out`. They use the ordinary 10x gene-matrix schema, with genome field `mm10`, 32,285 features and no spliced/unspliced count separation. The four raw matrices also lack that separation. Intron-inclusive gene counts cannot be renamed or decomposed into spliced/unspliced layers without read-level information.

The current original-full atlas H5AD is:

`/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/runs/20261013_v1_0_5/results/scanpy/original_full/allcells/scanpy.h5ad`

Its header reports 108,062 nuclei × 32,285 genes; metadata includes exact sample, barcode, nucleus_id, condition and sequencing partition. Original PCA/UMAP and Harmony representations are present. Its only layer is `counts`. Each of the twelve native Scanpy H5ADs (four populations × three count branches) has only that layer. The three shared branch objects and LIANA input have no layers. Twelve historical Allen-mapping objects and three SCENE objects also lack genuine velocity layers. The SCENE embedding is not a scVelo result.

No BAM, loom, GTF or FASTQ was found in bounded recursive filename searches of Desktop/data, current codex_new projects excluding package libraries, or the supplied network folder `/Volumes/Neurowissenschaften/-Ali_on_server/snr_p60`. This does not establish absence elsewhere on the network or on LiSC.

Two older Desktop `full_adata_processed.h5ad` copies are truncated (actual file lengths approximately 479 MB and 229 MB; declared HDF5 lengths approximately 6.50 GB). Their open errors are recorded in JSON. Current native atlas objects opened successfully. These old damaged files should not be substituted for the native atlas.

The existing `.venv_scanpy` provides Scanpy 1.12.4, AnnData 0.13.4, NumPy 2.5.4 and h5py 3.16.0; scVelo and velocyto were not installed there. This is an environment observation, not a recommendation to upgrade that preserved environment.

## Historical remote evidence, not a current remote inventory

The 8 September 2026 review and table at

`/Users/mujahid357/Desktop/data/snkr_cowork/seurat_out/BAM_validation_11_samples_20260908/`

describe 11 readable, coordinate-sorted Cell Ranger 10.1.0 BAMs and indexes at:

`/lisc/data/scratch/molecular_neurosciences/i-med/scdata/cellranger_out.new/<sample>/outs/possorted_genome_bam.bam`

with `.bam.bai` indexes. BAMs totaled 285,393,910,743 bytes before indexes. Historical reference information was mm10-2020-A / GRCm38, filtered GENCODE vM23, include_introns=true, with corrected CB/UB tags observed during targeted inspection. The exact current GTF path and digest are not established by this inventory. Rerun source identity must be checked against the original atlas rather than assumed from matching sample names or summary metrics.

| Sample | Original filtered matrix today | Historical BAM, 8 September | Genuine velocity layers inspected today |
|---|---|---|---|
| cKO_1138 | Present | 24.96 GB, quickcheck PASS | Absent |
| cKO_1139 | Present | 30.29 GB, PASS | Absent |
| cKO_1160 | Present | 32.16 GB, PASS | Absent |
| cKO_1165 | Present | 27.77 GB, PASS | Absent |
| cKO_933 | Present | Missing; replacement pending then | Absent |
| cKO_991 | Present | 24.72 GB, PASS | Absent |
| wt_1122 | Present | 23.11 GB, PASS | Absent |
| wt_1123 | Present | 24.35 GB, PASS | Absent |
| wt_1178 | Present | 24.99 GB, PASS | Absent |
| wt_1220 | Present | 24.18 GB, PASS | Absent |
| wt_1226 | Present | 25.65 GB, PASS | Absent |
| wt_1228 | Present | 23.22 GB, PASS | Absent |

The registered target `lisc-satb2-office-20261013` points to login02.lisc.univie.ac.at and workspace `/lisc/data/scratch/molecular_neurosciences/i-med/new_snRNA_seq`, with Slurm partition xeon_0384. Purpose-built `inspect_compute_target` returned **SSH access to the configured host was denied**. No terminal-SSH workaround or change to access controls was attempted. Therefore all current remote existence, quota, runtime and cKO_933 availability remain unknown.

## Resource observations and supportable execution conditions

The Mac filesystem had about 6.1 GiB free; the mounted network share had about 3.6 TiB free. These observations do not measure LiSC scratch quota. Local BAM copies/sorts are unsuitable under the observed free space. Velocyto may create a barcode-sorted BAM beside its input, so a future additive workflow must explicitly stage/manage these intermediates in new scratch space without writing into original output folders. Estimate requested resources from actual verified inputs and a controlled pilot; the documentation's example runtimes are not a project ETA.

For scale only, a dense float32 108,062 × 2,000 matrix is about 0.81 GiB, and the two moment matrices, velocities, likelihood fit state, neighbors and copies add to that. This arithmetic is not a measured memory requirement. Keep counts sparse, use independent lineage runs and controlled threading, and preserve native failures rather than expanding settings until an attractive map appears.

## Conditional handoff and evidence contract

1. Restore the existing LiSC connection, freshly verify all 12 BAM/indexes (especially cKO_933), sample-plus-barcode provenance, actual matching GTF/reference and tag compatibility. Do not silently replace original atlas counts with a newer Cell Ranger rerun.
2. Use a separately frozen environment and native Workbench plan for additive spliced/unspliced/ambiguous counting. Preserve raw layers and exact gene and sample/barcode axes. Join to the frozen QC/annotation roster with explicit matched, unmatched and excluded counts per mouse. No invented layers, DEG-only RNA background or reuse of normalized/CellBender values as splicing classes.
3. Examine nuclear splicing fractions, detection/phase coverage, batch and genotype structure before fitting. Analyze CA1/ProS and inhibitory identities separately; do not infer a CA1-to-interneuron lineage or force distinct mature neuronal identities into a single path.
4. Preserve distinct deterministic, stochastic and dynamical results, model-specific neighbors/graphs and diagnostics. Dynamical fitting requires recover_dynamics. Latent time is conditional on supported kinetics and defensible root/end-state evidence. A weak or unidentifiable result is a reportable outcome.
5. Required outputs: per-mouse layer QC and axis receipts; source/config/reference/environment hashes; sparse layered H5ADs; phase portraits; fitted-parameter/likelihood/coverage tables; confidence and high-dimensional directional agreement; per-mouse/partition support; original-coordinate stream/quiver maps; explicit failed/unavailable fields; documented root and root-sensitivity analysis if applicable. Existing atlas cosmetics must not alter scientific coordinates.

## Claim limits

Inference from this design: a single P60 mature-neuron snapshot does not by itself identify developmental time or transitions between deep and superficial anatomical layers. An arrow field or latent-time gradient is not proof of differentiation, elapsed hours, or cKO-to-WT progression. Genotype must not define the earlier/root state. Changes in nuclear processing or genotype-specific kinetic regimes require inspection rather than assuming a common kinetic model. Treat any supported directed dynamics as exploratory and retain independent-mouse evidence for genotype comparisons.

One immediate user-owned prerequisite is re-establishing the existing LiSC connection. The known BAM root can then be checked before asking for a new data upload or a different compute setup.

## Primary method links, accessed 28 October 2026

- [scVelo input requirements](https://scvelo.readthedocs.io/en/stable/getting_started.html)
- [scVelo velocity modes](https://scvelo.readthedocs.io/en/stable/scvelo.tl.velocity.html)
- [scVelo dynamical modeling](https://scvelo.readthedocs.io/en/stable/DynamicalModeling.html)
- [scVelo latent-time API and root arguments](https://scvelo.readthedocs.io/en/stable/scvelo.tl.latent_time.html)
- [scVelo phase portraits and preprocessing guidance](https://scvelo.readthedocs.io/en/stable/VelocityBasics.html)
- [velocyto native BAM, barcode, GTF and sorting requirements](https://velocyto.org/velocyto.py/tutorial/cli.html)
