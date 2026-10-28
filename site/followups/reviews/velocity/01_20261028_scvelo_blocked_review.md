Takeaway: Blocked: SATB2 RNA velocity needs verified spliced/unspliced counts and restored LiSC access.

The requested RNA-velocity analysis would test whether splicing dynamics support directed expression changes in the WT and Satb2 cKO nuclei. The inspected study objects contain ordinary RNA counts but lack the separate spliced and unspliced counts needed for scVelo. Historical records describe eleven Cell Ranger BAMs on LiSC, but the current connection was denied and the twelfth sample's read-level input remains unverified. No velocity counting, model fitting, velocity map or latent-time analysis has run. Consequently, neither temporal ordering nor genotype differences in velocity can be inferred from the available material. The next step is to restore the existing LiSC connection and verify all twelve matching BAMs and references before designing an additive counting run.

## Scientific context and question

The existing study has six independent WT and six Satb2 cKO P60 mice. The requested methods are deterministic, stochastic and dynamical scVelo, with latent time conditional on adequate kinetic evidence. Existing CA1 deep/superficial labels describe molecular programs rather than experimentally established anatomical positions or successive temporal stages. CA1 and inhibitory identities must be analyzed separately. Genotype and sequencing partition are partly confounded, so any subsequent genotype comparison needs explicit partition checks.

## Question, lifecycle, and endpoint

**Input blocked; no scientific workflow was executed or registered for this task.** No velocity environment was installed and no count layers were generated. The read-only assessment inspected original matrices, H5AD headers, bounded local/network filename searches, historical remote BAM-validation records and the registered LiSC target's returned accessibility result. It did not alter the original data or attempt an alternative SSH route.

## Key findings

- All twelve current filtered Cell Ranger matrices are ordinary gene-count matrices with 32,285 features. The four supplied raw matrices also lack splicing-class separation.
- Each of the twelve inspected native Scanpy H5ADs has a `counts` layer only. Shared branch and LIANA inputs have no genuine velocity layers. Intron-inclusive RNA counts do not provide independent spliced and unspliced abundances.
- A historical 8 September audit described eleven indexed BAMs totaling 285,393,910,743 bytes under `/lisc/data/scratch/molecular_neurosciences/i-med/scdata/cellranger_out.new/<sample>/outs/`. Historical `cKO_933` input was missing. Current BAM, index, reference and GTF identities are unknown.
- Inspection of registered target `lisc-satb2-office-20261013` returned **SSH access to the configured host was denied**. Its configured alias is `codex-satb2-lisc-office-20261013`, host `login02.lisc.univie.ac.at`.
- No BAM, loom, GTF or FASTQ was found in the bounded local and specified network searches. This is not proof that the files are absent elsewhere.

## Interpretation

The blocker is missing verified splicing inputs and remote access, not a failed biological velocity fit. scVelo explicitly requires separate spliced and unspliced count matrices ([official input documentation](https://scvelo.readthedocs.io/en/stable/getting_started.html)). Its dynamical model estimates kinetics from their joint phase relationships ([official model description](https://scvelo.readthedocs.io/en/latest/about.html)). Total RNA, expression-based pseudotime or artificially rearranged UMAP coordinates cannot substitute for that evidence. Even after counting succeeds, a single mature P60 snapshot may yield weak or unidentifiable kinetics; unavailable dynamics must remain a valid outcome.

## Artifacts

**Results:** No velocity matrices, fitted rates, velocity vectors, kinetic likelihoods or latent times were generated.

**QC and visual reports:** No velocity QC or figures exist. The read-only input assessment is at `/Users/mujahid357/Documents/codex_new/SATB2_FOLLOWUPS_20261028/velocity_readiness/STARTING_POINT_ASSESSMENT.md`.

**Provenance:** Observed HDF5 headers, original paths, historical records and accessibility evidence are recorded in `/Users/mujahid357/Documents/codex_new/SATB2_FOLLOWUPS_20261028/velocity_readiness/inventory_evidence.json`. Historical evidence is explicitly distinguished from current verification.

**Model synthesis:** This review records why the available data cannot yet support the requested biological endpoint. Because no run exists, no Workbench run-review submission applies.

## Limitations and blockers

Current LiSC files and quota cannot be inspected until access is restored. The historical BAM/reference metadata must be freshly checked against the exact source counts, samples and barcodes; matching sample names alone are insufficient. The `cKO_933` replacement location is unknown. An incomplete eleven-mouse analysis would change the requested design and has not been substituted. No conclusions about absence of dynamics or equivalence between genotypes are supported.

## Recommended next step

Restore the existing LiSC connection; then verify all twelve BAM/index pairs, corrected barcode/UMI tags, matching GTF/reference and exact sample/barcode joins. With valid inputs, prepare a separate immutable native Workbench plan for spliced/unspliced/ambiguous counting on LiSC scratch. Inspect layer coverage before approving model-specific velocity branches, and report latent time only if kinetics and root evidence justify it.
