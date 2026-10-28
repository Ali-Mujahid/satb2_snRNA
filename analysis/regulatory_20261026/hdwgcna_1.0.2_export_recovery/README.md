# SATB2 hdWGCNA export recovery 1.0.2

The execution entrypoint for this version is `workflow/RecoverSnakefile`. It
copies the completed1.0.1 preparation and four network checkpoints into a fresh
root under strict SHA256/MD5 verification, runs the unchanged enrichment
specification with explicit SVG/Cairo PDF devices, and packages the results.
It invokes no network fitting rule and preserves both native and recovery
configuration identities. See RECOVERY_PLAN.md. The original full-workflow source
and the byte-identical native fitting script remain included for reproducibility.

## Preserved original scientific workflow

Native entrypoint: `workflow/Snakefile`. Implements the user-selected hdWGCNA analysis using package 0.4.12, commit `e3344d1f7bbac4264adf94f7aa31e0802fa8282d`, and WGCNA 1.74. It does not rerun the original nuclei clustering, annotation, DESeq2 fits or ligand–receptor scoring. The five frozen Seurat populations and raw integer counts are content-hash bound and checked against all-gene mouse pseudobulks.

Stages: sparse exact input extraction → separate pooled CA1/IN signed bicor hdWGCNA networks → conditional WT/cKO consensus sensitivity → positive-RNA-oriented PC1 module eigengenes → independent-mouse summaries/continuous-score associations and leave-one-mouse-out stability → offline GO and matching original DEG enrichment → source-preserving export/validation.

The reviewed scientific plan is `ANALYSIS_PLAN.md`. Read `config/README.md`, `docs/PREPARE_INPUTS_CONTRACT.md` and `docs/HDWGCNA_CONTRACT.md` before binding inputs. Packages are preinstalled in a separate project library; this workflow contains no installation hooks or network downloads. It uses the existing project Snakemake 9.27.0 controller explicitly selected by the user. Its Host PATH provenance is a reproducibility limitation; package runtime and downloaded artifacts are separately recorded. R is x86_64 under Rosetta on an arm64 Mac.

No scientific execution is implied by this source. Run only through the NGS Analysis Workbench immutable plan and native approval. Use a fresh absolute run root. No source analysis may be overwritten. A poor topology fit produces an honest unsupported branch status rather than an automatic power fallback. Conditional consensus failure does not invalidate an otherwise verified primary network, but must remain visible.

Published module weights are undirected co-expression correlations, not direct TF targets or physical ligand–receptor evidence. Continuous eigengene associations are not gene DEGs. Curated TF activity and motif-supported regulon inference are separate methods and are not claimed by this workflow.

External packages retain their licenses: hdWGCNA is GPL-3; WGCNA is GPL-2-or-later. Original project wrappers and source provenance are supplied for reproducibility; no third-party license is replaced. No externally downloaded tutorial data or figures are included.

Revision1.0.1 declares the preparation completion receipt only, avoiding engine-created nested directories before the strict empty-output guard. Matrices remain required, content-hash-verified artifacts in that receipt. No scientific method change.
