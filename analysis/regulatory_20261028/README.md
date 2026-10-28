# Additive native mouse pySCENIC branch

The user-selected full pySCENIC branch was launched on 28 October 2026 after the official matching mouse cisTarget resources became reachable. This source snapshot records the executed workflow, immutable configuration, plan/run receipts and verified reference metadata. A source snapshot and a running workflow are not evidence of completed regulon inference. Actual results will be added after terminal completion, native validation and independent review.

The frozen cohort contains CA1 deep-associated and superficial-associated molecular programs plus inhibitory C1/C6/C8, with 23,143 nuclei and all six WT/six cKO mice. Native GRNBoost2, positive modules, matching mm10 v10-cluster motif pruning and actual nucleus AUCell are performed separately in CA1 and inhibitory lineages. Continuous independent-mouse mean-AUC models use their own global BH correction and matched-partition sensitivity. Original DESeq2 gene results and strict DEG-only communication outputs are unchanged. Regulatory inference uses the full eligible genotype-neutral RNA background.

Native identity:

- Workflow: `satb2-pyscenic-mm10-20261028`, version `version-9463ef292fcc37085be876e0eb20068c`.
- Source SHA256: `57028aaf028ab1f0866a4115b2d2409a5f0b4ea5eed03dd3d18bd2cc72d4d52e`.
- Registry run: `5f292759-e4b6-4e2a-be49-b9c45afa9db7`.
- Plan: `ngs-plan-dd5207c61fb0690c`, checksum `dd5207c61fb0690ceb4a223e8b6f2e36f2a321ca52c86ad40d5b248390208c24`.

Read [the scientific plan](pyscenic_source/ANALYSIS_PLAN.md), [native producer contract](pyscenic_source/docs/NATIVE_PRODUCER_CONTRACT.md) and [mouse/transport contract](pyscenic_source/docs/REGULON_MOUSE_AND_TRANSPORT_CONTRACT.md). The two ranking files and large MGI annotation are reusable external references, with exact official HTTPS URLs, lengths, local SHA256 and available official SHA1 evidence in [READY_MANIFEST.json](reference_metadata/READY_MANIFEST.json). The current unchanged official mouse list includes Satb2; a reference annotation does not establish an inferred regulon. Installed runtime packages are not bundled.

The earlier [v3.7 regulatory results](https://github.com/Ali-Mujahid/satb2_snRNA/releases/tag/atlas-v3.7.0) and historical blocked-reference evidence remain preserved. New motif results require a separate additive atlas release after validation.
