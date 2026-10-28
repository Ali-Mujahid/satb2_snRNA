# SATB2 CollecTRI TF activity 1.0.0

Run native decoupler 2.3.0 ULM on the five exact frozen Seurat endpoints, followed by independent-mouse continuous-score associations. This is a separate workflow from hdWGCNA and pySCENIC. No installation or network calls occur during execution. Required inputs, eligibility and preserved native units are documented in scripts/run_tf_activity.py and the scientific plan.

Entrypoint: workflow/Snakefile. Configuration: absolute JSON with `schema=satb2-tf-activity/v1`, a fresh run_root, preinstalled python_executable and Rscript_executable, and hash-bound source files. Snakemake writes only new run outputs. The final validation receipt is the completion target; no failed or unvalidated outputs qualify for publication.

Outputs: results/tf_activity native scores, gene-set FDR, target coverage, curated prior edges and original Wald enrichment; results/tf_mouse continuous-score effects, global BH, sensitivity and stability; results/atlas offline transport and results/FINAL_VALIDATION.json. Gene-set FDR and genotype association FDR are distinct. Existing gene DESeq2 fits are reused without refits. The full eligible RNA background, not just DEGs, is required for ULM.

Runtime: isolated Python3.12/decoupler2.3.0 and project R4.6.1 library. User-selected existing Snakemake9.27 controller. Exact package/artifact and input hashes are bound in the configuration; the existing controller and base R closure are not fully lockfile-backed. Native tasks are sequential within this workflow, at most two requested threads and 4GB nominal task memory. A separate authorized hdWGCNA workflow may be running concurrently on this 64GB host.

pySCENIC is explicitly blocked on unavailable official mouse ranking/motif references. No curated prior network or expression association is called a motif-supported inferred regulon, direct binding or causal regulation.
