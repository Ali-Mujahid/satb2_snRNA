Takeaway: Failed: SATB2 hdWGCNA preparation rejected engine-created output folders before any network analysis.

The follow-up asks whether CA1 layer-program and inhibitory neurons contain reproducible co-expression modules related to Satb2 genotype and the existing DEG-only interaction candidates. The native workflow controller started, but its first preparation step stopped at the output-directory guard. Snakemake had created nested parent folders for the declared matrix outputs, and the guard interpreted those folders as existing analysis content. No counts were extracted, metacells constructed, network fitted, or enrichment and figure outputs generated. This operational failure provides no biological findings about co-expression or genotype. The smallest justified recovery is to declare the preparation completion receipt as the stage output, retain the strict no-overwrite guard, and submit a new immutable plan in a fresh run directory.

## Scientific context and question

Mouse hippocampal snRNA-seq; six independent WT and six Satb2 cKO mice. Frozen Seurat populations are CA1 deep/superficial-associated expression programs and inhibitory C1/C6/C8. Full eligible RNA is required for co-expression discovery; the existing LR endpoint selection remains strict DEG-only. No new gene DESeq2 model is requested.

## Question, lifecycle, and endpoint

Registry run b721fb33-f2f5-4829-b24f-8b0a72d0877e; status failed; controller returncode1. The only observed scientific process attempt was prepare, which failed before loading study inputs. No subsequent processes completed. Target local, run root /Users/mujahid357/Documents/codex_new/SATB2_REGULATORY_20261026/runs/20261026_hdwgcna_v1. Approved plan ngs-plan-6c1e92f27eb0c87e; workflow source SHA256 a3eed4d29857e43a458b3fe879440ce2e7db83c6aec3675a8f158d1818b2c673.

## Key findings

The task log logs/prepare.log contains ValueError: Refusing to overwrite nonempty output directory at prepare_inputs.py line371. The preparation completion receipt and matrix outputs are absent after Snakemake removed its failed-job outputs. Controller/native task logs establish execution failure only; no scientific measurement or QC result was produced.

## Interpretation

The guard behaved conservatively. The wrapper declared every nested matrix output, causing engine parent-directory creation before a script that expected an empty output root. Declaring only the stage's content-hash completion receipt avoids that mismatch while preserving the original guard and later verification of every generated matrix. This changes workflow bookkeeping, not biological methods.

## Artifacts

Results: unavailable. QC/visual reports: unavailable. Provenance: copied immutable workflow, approved configuration, native controller log and task log under the registered root. Model synthesis: this failure review. No previous scientific inputs were modified by the failed task; its guard preceded extraction and any scientific output write.

## Limitations and blockers

No module, hub, eigengene, genotype, enrichment, TF-activity or regulon conclusion is supported. Earlier isolated package probes and synthetic code tests are runtime/development evidence, not study results. TF method selection remains separate and pending.

## Recommended next step

Version the native output declaration, preserve this failed run and logs, and use the Workbench readiness/immutable-plan/native-execution path for a fresh recovery run. Do not launch fitting directly from the terminal or reuse an old analysis directory.
