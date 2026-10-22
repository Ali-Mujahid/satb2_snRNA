# SATB2 Seurat-defined CA1–inhibitory LIANA 1.0.0 wrapper

This new local Snakemake workflow implements the user-selected LIANA analysis
with its native CellPhoneDB scorer. It reads the completed Seurat labels and
existing DESeq2 evidence without changing or rerunning the original analyses.
The exact scientific design is in the sibling task `design/ANALYSIS_PLAN.md`.

Native entrypoint: `workflow/Snakefile`. Required absolute JSON configuration is
described in `config/README.md`. Preparation streams sparse original counts,
scores native LIANA 2.0.0 and mouseconsensus, exports complete and DEG-supported
tables, then renders all three directed views and validates the result package.
No installs, raw-read processing, new QC/clustering/DE, external messages,
network resource retrieval or previous-output overwrites occur during the run.

The user selected reuse of the existing Snakemake 9.27.0 project controller.
Task Python is an independently prepared isolated arm64 Python3.12 environment;
the full dependency artifact hashes and LIANA source checksum are retained in
the task `environment/` directory. No Conda/Docker deployment is required.
The native runtime contract is in `runtime_contract.json`.

All heavy rules request four cores so they execute sequentially. Matrices stay
sparse and the full collapsed feature universe is retained for native background scores. A new
run_root is mandatory. Native readiness/plan approval governs scientific
execution. The output is exploratory transcript-compatible communication;
WT/cKO interaction-score comparisons are descriptive, not genotype tests.
