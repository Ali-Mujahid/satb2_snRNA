# SATB2 — separate reviewed follow-ups

These analyses preserve all previous results and original DESeq2 fits. Six WT and six cKO mice are the biological units. The common sequencing-partition analysis contains six WT and three cKO mice; it is a sensitivity comparison within the same study, not independent replication. Partial genotype/partition confounding and uncertain identities remain limitations.

## 1. Relative proportions: scCODA

The strict same-model numerical check completed all 16 models. Four pass unchanged sampling diagnostics and12 remain inconclusive (81 retained divergences). Both excitatory common-partition models pass; Scanpy inhibitory common-partition and Seurat inhibitory full-cohort also pass. No broad-population or CA1-program fit qualifies for inference. There is no further automatic tuning.

|Mixed excitatory target|Saved flagged reference|Native PIP|Native compositional log2 fold change|
|---|---|---:|---:|
|Scanpy C9|C14|0.993333|+0.646680|
|Seurat C4|C21|0.995958|+0.625121|

These same-data mixed targets largely overlap (common-partition Jaccard 0.9735), but their identities and references remain uncertain. The coefficients concern relative recovered-nucleus proportions; they do not prove a named cell type expanded or absolute tissue abundance increased. The two qualified inhibitory models select no effects, which does not prove equivalence. Native expected Bayesian FDR 0.05 is per model, distinct from gene-DE statistics. The initial family, with0 of 16 passes, remains fully available alongside the strict check.

## 2. Both-end DEG-only CellChat

Nine endpoint contexts across six molecular pairs and five pathways pass the original matching-endpoint Seurat DESeq2 DEG rule, including every required ligand and receptor subunit. Three primary triMean full-cohort associations have lower cKO scores: CA1 deep-associated Slit1 to inhibitory C6 Robo1, C6 Tenm3 to CA1 deep-associated Adgrl2, and C6 Tenm3 to CA1 superficial-associated Adgrl2. Each full-cohort BH value is 0.0064935; each corresponding primary common-partition BH value is 0.0714286, above 0.05. None of nine primary common-partition tests qualifies.

The separately specified truncated-mean sensitivity has seven full-cohort and five common-partition associations. It does not independently confirm primary results. Both methods preserve the three full-cohort directions. Every retained pathway context contains one eligible pair, so pathway and LR results duplicate evidence. There is no retained C8 edge, which does not establish absence of biological communication. Same-data DEG selection makes subsequent score associations exploratory. RNA scores do not measure protein binding, pathway activation or synaptic function. Original component RNA up/down signs remain separate from interaction-score differences.

CA1→inhibitory, inhibitory→CA1 and both-direction views each have saved dot, tile, circos and circle figures plus mouse-level plots, PDFs and native tables. Readable CellChat layouts correct margins, headers and vertical jitter only; original figures remain archived with explicit display warnings.

## 3. Monocle3 expression-state graphs

|Scope|Nuclei|Graph-associated genes|Modules|
|---|---:|---:|---:|
|CA1 deep/superficial programs|20,456|9,056|47|
|Inhibitory C1|1,718|3,669|19|
|Inhibitory C6|501|1,439|8|
|Inhibitory C8|468|2,002|11|

Each scope contains all 12 mice. CA1 has two graph components; the smaller 847-nucleus component is mostly superficial-associated and includes all 12 mice. Primary graphs are unrooted. Genotype-blind arbitrary-root distances and alternate references are provided, but changing reference materially changes rankings. A mature P60 single snapshot supplies no independent biological temporal root. These graphs are not a demonstrated WT-to-cKO or deep-to-superficial transition. Native graph_test screens cell-level autocorrelation, not independent-mouse genotype DE. Mitochondrial/astroglial-associated modules require quality/background interpretation. Full graph-gene/module tables, native CDS objects and support indexes are retained.

## 4. scVelo input blocked — no fit

Inspected objects lack genuine separate spliced/unspliced layers. Historical LiSC records describe11 BAMs; cKO_933 inputs and current access remain unverified. Intron-inclusive counts are not unspliced measurements. No sample was silently omitted and no synthetic vectors or latent times were produced. Restore the existing LiSC connection and verify all 12matching BAMs/indexes/barcodes plus GTF/reference before a separately reviewed counting/velocity plan.

## Files and use

The complete flat offline ZIP opens through index.html without Python/R when kept intact. The thin online reader loads some exact verified media from the release. Separate native archives preserve both sampler families, CellChat objects, Monocle CDSs/indexes, selected integer counts, original DESeq2 tables, source/configuration and all reviews. Native objects require compatible recorded software. The analysis and display reviewers' exact receipts remain authoritative. The independently running pySCENIC analysis is separate; none of these follow-ups establishes completed motif-supported regulons.
