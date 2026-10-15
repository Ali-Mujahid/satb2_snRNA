The completed marker/reference evidence supports useful parent-level and “-like” display options for every uncertain original cluster. These are nonbinding editor choices: the scientific annotations, retained nuclei, frozen cluster/category memberships and statistical fits are unchanged.

Across the eight overlapping original views (Scanpy/Seurat × All/EX/IN/CA1-ProS), 70 of163 cluster records fail the existing conservative fine-label guard. This comprises61 retained records and9 excluded/audit-only records.66 have explicit mixed, low-reference or conflicting-signature status, while four have other limitations such as absent representative support or an emptied main mask. After the suggestions there are still70 scientifically uncertain records and66 explicit mixed/low/conflicting records. Every affected record now has a provisional named parent/resemblance option; zero records were promoted to a validated identity. These are overlapping cluster records, not70 independent cell types or disjoint biological populations.

| View | Uncertain before | Uncertain after provisional options |
|---|---:|---:|
| Scanpy All |16|16|
| Scanpy EX |7|7|
| Scanpy IN |6|6|
| Scanpy CA1-ProS |1|1|
| Seurat All |24|24|
| Seurat EX |8|8|
| Seurat IN |6|6|
| Seurat CA1-ProS |2|2|

Parent proposals prioritize native canonical marker support, reference composition and the original supported representative evidence. Quantitative pooled and WT-only detection/mean/specificity/mouse-support values are saved with stable gene IDs. Mere gene detection is not lineage specificity: the separate specificity-qualified marker sets require fraction detected≥0.20 and a detected-fraction contrast≥0.10 versus the existing rest population. The source marker panel has no duplicated symbols, and the builder asserts that no stable gene IDs are silently collapsed. These descriptive criteria are neither calibrated cell-type probabilities nor new inferential tests.

Useful examples include CA3/CA2-reference-like glutamatergic mixtures, Lamp5/Lhx6-reference-like GABAergic mixtures, astrocyte-like RNA with additional myeloid-associated signals, and vascular-like RNA associated with endothelial/pericyte references. Native conflict flags remain visible even when an RNA parent is recognizable. Broad immune genes can be compatible with a microglial-like parent; a conflict flag does not by itself prove a mixture. Aggregate detection cannot distinguish cells of different types, same-nucleus mixed RNA, ambient RNA or multiplets, so none of those explanations is asserted.

For IN, Sncg/Ndnf and Pvalb/Sst proposals explicitly describe reference mixtures. In Scanpy IN1, WT Sncg detection is approximately5.0% versus Ndnf29.3%; in IN4, WT Pvalb is9.8% versus Sst48.8%. A “Sncg/Pvalb-reference-like” name must not be read as a measured high-expression subtype in those records. LSX-associated is atlas provenance, not demonstrated anatomical origin; “Cholinergic neurons”, basket, bistratified and OLM worksheet headings remain expression programmes rather than forced identities. Rare IN/CA1 records remain available; low model coverage and marker identity confidence are separate questions.

CA1 identity and expression programme are kept in separate fields. All20 original CA1 clusters have programme display options;11 remain mixed/weak or identity-conflicting before and after. Nine have an operational deep/superficial-associated WT RNA programme, but every physical-layer assignment remains UNDETERMINED. Seurat CA1 cluster7 has a CA1/ProS-reference-like parent with weak local marker contrast and no supported representative candidates. Scanpy8/Seurat11 remain neuronal–myeloid/astroglial mixed-RNA records. The WT-score rules are operational and threshold-sensitive; display aliases do not change their primary programme memberships or create an anatomical layer assignment.

The JSON schema is `satb2-annotation-suggestions-2`. `suggestions[]` contains70 identity proposals; `by_key` maps `tool::population::original_cluster_id` to array index. `programme_suggestions[]` contains20 independent CA1 programme options, with `programme_by_key` lookup. `categories[]` provides grouping, and `count_summary` preserves exact before/after counts. Editor integration should offer explicit optional display-label choices, retain native labels/evidence alongside them, and carry aliases across paired branches using the ORIGINAL frozen cluster ID. Fresh paired cluster numbers cannot inherit these aliases by numeric coincidence. Excluded/empty All/EX records stay excluded and their reference-only options remain audit-only. No proposal changes a category union or licenses a different cached DE fit.

The evidence table is `cluster_label_suggestions.tsv`, editor metadata is `editor_label_suggestions.json`, and source/hash/count provenance is `REVIEW_PROVENANCE.json`. The five peer-reviewed IN/CA1 examples passed conservative wording checks. No source RNA extraction/aggregation, new scientific fit, native output modification or v1 replacement occurred.
