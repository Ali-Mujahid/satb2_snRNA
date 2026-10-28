# Independent read-only audit of the additive pySCENIC branch

The audit is outside frozen scientific source, under `project/reviews/`. It never modifies a native run, imports pySCENIC to rescore data, loads Python/R scientific objects, fits a model, or changes any P value/FDR used in results. Native per-nucleus AUC means, effects and BH arithmetic are independently reproduced only to compare exported numbers. CSV/gzip streams, eight-MiB hash blocks and endpoint–mouse running sums avoid holding expression/AUC matrices in memory.

`audit_pyscenic_readonly.py` checks completed seals only. Default mode requires final validation; `--allow-partial` returns only completed checkpoint summaries/hashes, with `PARTIAL_CHECKPOINT_AUDIT` status. It never polls or calls the Workbench controller. Run a partial check only when a useful stage completion is observed, and a final audit when root reports completion. The output path must lie outside the entire native run and must be new.

Example final invocation:

```sh
/Users/mujahid357/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3 \
 /Users/mujahid357/Documents/codex_new/SATB2_REGULATORY_20261026/reviews/audit_pyscenic_readonly.py \
 --config /Users/mujahid357/Documents/codex_new/SATB2_REGULATORY_20261026/runs/20261028_pyscenic_v1/config/pyscenic.v1.json \
 --results /Users/mujahid357/Documents/codex_new/SATB2_REGULATORY_20261026/runs/20261028_pyscenic_v1/results \
 --output /Users/mujahid357/Documents/codex_new/SATB2_REGULATORY_20261026/reviews/20261028_pySCENIC_independent_final_QA.json
```

The exact copied config basename should be confirmed from the native run directory. `--allow-partial` can be appended with a separately named output, never overwriting final QA.

## Checks required for final completion

1. **Unchanged evidence:** all completed stage seals/config/source bytes, every final native/transport receipt, and all original configured input/reference records match SHA-256/size. Native producer and downstream package final schemas/statuses match the additive pySCENIC branch. No prior DESeq2/LR/hDWGCNA/CollecTRI result is edited.
2. **Exact cells:** all23,143 selected nuclei retained, distinct across lineages; CA1 8,750 DL +11,706 SL, inhibitory C1 1,718/C6 501/C8 468. Raw native AUCell row order equals matching metadata and the full original prepared lineage order. Exact mouse/endpoint/genotype identities match; the thirteen C8/cKO933 nuclei remain scored.
3. **Actual measured AUCell:** raw matrix columns equal eligible lineage-distinct native positive regulons. AUC is finite and in[0,1]; legitimate zeros remain zeros. Every exported mean equals the streamed mean across **all exact nuclei** belonging to its endpoint/mouse, using compensated sums. This also checks supported Satb2 self-excluded raw AUC. No cluster-mean substitution or nucleus pseudo-replication.
4. **Native target eligibility:** unique native targets, positive unmodified gene2weight, explicit TF-selfgene origin/flag, full genotype-neutral background. Each eligible set has at least10 measured native targets; low-target sets remain unavailable. Candidate and final native contexts are activating/positive, retained motif NES≥3 with annotation; context is motif-resource support rather than measured binding or proven activation. Actual native ctx outputs remain available.
5. **Independent-mouse effects:** all eligible endpoint–regulon means cover six WT/six cKO; primary native effects equal those mouse mean differences. Same-partition six WT/three cKO is selected by exact partition metadata and has matching axes/effects. LOO omits each of twelve actual mice once, with correct five/six remaining genotype counts; bootstrap summaries declare200 within-genotype independent-mouse iterations. No new model is fit by QA.
6. **Multiplicity families:** recompute arithmetic BH from exported native P values, jointly over all tested primary endpoint×lineage-specific regulons; same-partition and Satb2 no-self each have distinct global families. Untested/constant cases retain missing P/FDR and are not given zeros. Native aliases padj/association_fdr agree. These are continuous-score association FDRs, not gene-DESeq2 or motif significance.
7. **SATB2 honesty:** actual current official TF-list membership and eligible RNA are distinct from discovered motif regulons. No candidate is forced; absent motif support has no activity proxy. A supported no-self sensitivity removes only the native own-gene member, with unchanged background and≥10 remaining targets; all five unavailable/supported states are explicit. Original full regulons and primary BH remain unchanged.
8. **Interpretation/visual review:** root or reviewer examines native AUC distributions, top TF/regulon effects and coverage, same-partition/sign stability, native motif evidence, at least one actual PDF and browser exports. Gene2weight/arrows are candidate orientation, not binding strength/causal regulation. Layer programs are not measured anatomical layers. C8 mixed identity, partial partition confounding, pooled-discovery reuse of the same cohort and cross-lineage AUC non-calibration remain visible.

The machine-readable JSON provides numeric errors/counts and explicit limitations, but it is not a substitute for a model-authored scientific interpretation/review. Results should only be published after the native run has terminal completion and root independently verifies native final gates.

## Synthetic verification

`test_audit_pyscenic_readonly.py` uses entirely synthetic tiny tables to test BH arithmetic, missing FDR on untested rows, actual nucleus means/zeros, axis/mean corruption detection, and honest empty-regulon axes. It performs no study analysis or scientific workflow execution. Four groups pass before the initial bounded checkpoint audit.
