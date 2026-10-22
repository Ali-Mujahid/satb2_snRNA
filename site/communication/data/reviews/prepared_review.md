# Prepared-input acceptance review — 22 October 2026

PASS. No blocking finding. This review read existing sources and native prepared
outputs only. It did not modify outputs, execute preparation, fit a model, or
change any scientific result. Validation used calculations separate from the
preparation completion receipt.

- Actual `lr_expression.h5ad`: **29,367 unique nuclei × 32,245 collapsed RNA symbols**.
- Exact frozen Seurat membership, row-to-native-ID mapping, group, mouse and genotype matched independently.
- Groups: CA1 Deep 8,750; CA1 Superficial 11,706; IN 8,911 across C1–C20. IN_Total is only a supplemental metadata row, not duplicated cells.
- Complete normalized RNA background preserved: **1,597 resource symbols + 30,648 nonresource symbols**.
- All 32,285 native stable feature IDs remain mapped. **78 source rows / 38 ambiguous symbols** are flagged; they are excluded from primary symbol-based DEG support.
- Full selected RNA library: **322,680,709 UMIs**; post-collapse all-gene UMIs equal the same total. Every per-cell denominator equals native `nCount_RNA`.
- Reviewed **96,150,440 stored expression values**. Maximum full-row inverse-log CP10k sum deviation: **0.001174601**; maximum inverse-normalization distance to an integer: **0.00076422406 UMI**. Every row reconstructs its exact integer total; recovered grand total **322,680,709**.
- One nucleus per each of 22 groups also passed a complete raw-count→symbol-collapse→log1p(CP10k) comparison, maximum per-value tolerance 1e-6.
- Resource coverage independently matches **3,946 complete interactions / 3,989 resource rows**; **43** have unavailable subunits. Original resource pair order is preserved.
- Every output hash matches its receipt. Every recorded source, reference, endpoint and configuration hash still equals its before/after digest.
- Per-mouse counts and group totals agree with H5AD. C18 remains NOT_TESTED; C20 conflicting identity remains visible.

## Interpretation constraints retained

CA1 groups are molecular expression programs, not established anatomy. All IN
fine annotations retain review flags. Native ranking and CellPhoneDB cell-label
permutation results are not replicate-aware genotype interaction tests. The
existing DESeq2 endpoints remain reused evidence; C18 missing DE is not zero.

Validated checks (including individual file hashes): 154.
