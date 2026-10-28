# Independent IN native GRNBoost2 checkpoint audit

Takeaway: The sealed IN GRNBoost2 checkpoint passes read-only numeric, identity and eligibility checks.

The native GRN checkpoint completed at 2026-10-28T09:32:42.589565+01:00 (Europe/Vienna). This bounded audit checks completed preparation and GRN files only. It does not fit a network, prune motifs, score activity, test genotype effects, or edit native files. The results describe an unsigned expression-prediction adjacency at this stage, not motif-supported regulons or measured transcription-factor binding.

| Check | Observed result |
|---|---:|
| Exact original nuclei | 2,687 |
| Exact endpoint–mouse groups | 36 |
| Independent mice | 6 WT + 6 cKO |
| Unique, unambiguous, nonconstant eligible RNA genes | 11,051 |
| Symbols in frozen current official TF list | 1,860 |
| Eligible official input TFs | 973 |
| Active sources in native GRN | 973 |
| RNA target genes represented | 11,051 |
| Native edges streamed | 1,697,355 |
| Finite positive importance edges | 1,697,355 |
| Zero / negative importance edges | 0 / 0 |
| Exact duplicate TF–target pairs | 0 |
| Native importance range | 3.478327375e-18 to 131.9736783 |
| Native GRN selfgene edges | 0 |
| Satb2 source / target edges | 0 / 0 |
| C8/cKO933 selected nuclei retained | 13 |

Every recorded native preparation/GRN output, frozen producer source and configuration matched its size/SHA-256 receipt. Relevant original count/gene/cell, mouse and official-TF inputs were separately hash-checked. Prepared nucleus IDs and ordering, endpoint identities and mouse genotypes equal the authoritative original prepared metadata exactly. All exact endpoints retain all twelve mice and their original nucleus counts. Eligible expression gene symbols are unique and match the genotype-blind >=5% detection/valid-symbol/nonconstant audit; TF inputs equal their intersection with the exact official list. No DEG or HVG target selection was introduced.

Satb2 is in the current official list but its IN native candidate row has `eligible_RNA=False`, reason `DETECTED_IN_FEWER_THAN_5_PERCENT_IN_EVERY_ENDPOINT`, and `forced_into_TF_list=False`. There are 0 Satb2 source edges in this checkpoint. The fixed RNA gate is distinct from biological absence of SATB2 function. Inhibitory C8 remains a mixed molecular population; the thirteen C8/cKO933 nuclei are retained.

All 1,697,355 compressed adjacency rows were streamed without loading an expression matrix or a scientific object. Source/target universes, finite numeric weights, deterministic sorting, and exact pair uniqueness passed. The exact duplicate check used a 1,344,078-byte bitset rather than a multi-million-row dataframe. The native edge file was SHA-256 checked before and after streaming.

Positive feature importance is not an activation sign. These raw edges cannot yet be called motif-supported regulons, direct binding, or causal regulatory interactions. No genotype-association claim follows from this pooled GRN. Further interpretation awaits completed native motif pruning, actual-nucleus AUCell, independent-mouse models and final validation. The historical analyses and every native seal remain unchanged.

Machine-readable evidence: [20261028_IN_GRN_checkpoint_independent_QA.json](/Users/mujahid357/Documents/codex_new/SATB2_REGULATORY_20261026/reviews/20261028_IN_GRN_checkpoint_independent_QA.json). Source checkpoint: [COMPLETE.json](/Users/mujahid357/Documents/codex_new/SATB2_REGULATORY_20261026/runs/20261028_pyscenic_v1/results/pyscenic_native/IN/grn/COMPLETE.json). The top ten raw importance rows are retained in JSON for audit only, without mechanistic ranking claims.
