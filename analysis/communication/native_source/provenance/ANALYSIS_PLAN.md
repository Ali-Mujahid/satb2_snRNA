# Seurat-defined CA1–inhibitory ligand–receptor analysis

Version 1.0.0; 22 October 2026. New analysis only; prior analyses are read-only.

## Objective and experimental unit

Predict expression-compatible protein ligand–receptor interactions between the
existing Seurat CA1 deep-associated and superficial-associated expression programs
and inhibitory populations. Deliver CA1→IN, IN→CA1, and their directed union,
with dot plots, tile/heatmaps, circos ribbons, circle networks and native LIANA
rank-aggregate tables. These three views use the same scored universe; the union
is not a third independent experiment and does not turn a directed interaction
into a reciprocal interaction.

The experiment has six independent WT and six independent Satb2 cKO mice.
Nuclei are observations within mice, not independent genotype replicates.
The requested Seurat-defined original_full memberships are frozen: 8,750 CA1
deep-program, 11,706 CA1 superficial-program and 8,911 inhibitory nuclei.
CA1 mixed/weak and conflicting programs are outside this requested scope.
All 20 inhibitory clusters remain available, including provisional/mixed labels.
IN total is a separate complementary resolution, not combined with its own
constituent clusters in the same scoring fit.

## Evidence-backed starting point

The `.audit/input_audit.md` identifies native integer full-gene RNA counts and
exact Seurat labels. The count H5AD is a shared storage container, not use of
Scanpy clustering. No HVG-only, Harmony, SCT or scaled matrix supplies expression.
Counts, annotations, sample metadata, endpoint maps and existing DESeq2 tables
are bound by SHA256 in `config/analysis.json`.

Existing PRIMARY mouse-level DESeq2 endpoints supply gene evidence: CA1 deep
P0244 (742 padj<0.05 genes), superficial P0247 (1,825), inhibitory total P0139
(4), and the matching IN cluster endpoints P0140–P0159. C18 was not tested and
must retain missing/untested status, rather than zero DEGs. No DESeq2 fits are
repeated. No edgeR, cell-level differential tests, new annotations, embeddings,
QC decisions or downsampling are included.

## Methods and endpoints

Use the fixed official LIANA 2.0.0 release at immutable source commit
`c59472ccc9de8360dbbf5016db75f8abde08dd3e`; the older PyPI 1.10.0 package has a
native duplicate-score ranking issue and is unsuitable for default aggregation.
Use LIANA's native implementations of CellPhoneDB, Connectome, logFC, NATMI and
SingleCellSignalR with the frozen mouseconsensus resource. This is LIANA's
CellPhoneDB implementation, not a claim of running standalone human CellPhoneDB
v5. Preserve all resource subunits, interaction provenance, installed versions
and downloaded source hashes. No uppercasing of mouse genes or unsupported
one-to-one human orthology assumptions.

Sum raw count feature rows sharing the same gene symbol, retain a stable-ID
mapping and ambiguity flags, normalize each nucleus to 10,000 total RNA UMIs
using its full-gene library, then log1p. Retain the full collapsed normalized feature universe for native LIANA
background-dependent scores; LIANA restricts features to the resource internally. Every complex subunit must
pass native expression/detection gates. Use expr_prop=0.10, min_cells=5,
1,000 native cell-label permutations, seed=1337, four threads. A subpopulation
with fewer than five nuclei within a mouse is unavailable for that mouse's
support estimate; those nuclei are not deleted from pooled/genotype fits.

Run the same native consensus for pooled, WT and cKO populations at IN-cluster
and IN-total resolutions, restricting eligible ordered pairs to CA1↔IN.
Append native geometric_mean on the exact unique source/target/ligand/receptor
key. Export the requested seven columns without renaming score meanings:
source, target, ligand_complex, receptor_complex, lr_logfc, lr_gmeans,
magnitude_rank. Save full native scores and extended evidence separately.

`lr_logfc` measures native source/target group-vs-rest specificity, not cKO/WT
change. `lr_gmeans` uses native mean log-normalized expression with complex
bottleneck handling. Lower magnitude_rank means stronger consensus priority;
it is not genotype FDR. Native CellPhoneDB permutation P values describe
within-dataset cell-label specificity and are not a replicate-aware genotype
test. WT/cKO differences in score or rank are descriptive.

For each mouse separately, compute native geometric_mean with the same expression
gates and pairs. Report both eligible mice and expression-supported mice by
genotype. Missing/ineligible mice are not zero expression. This guards against
one mouse driving a pooled candidate without inventing a differential
interaction P value.

Preserve full expression-qualified candidates, then produce a separately labeled
DESeq2-supported view: a ligand subunit is padj<0.05 in the sender's matching
endpoint OR a receptor subunit is padj<0.05 in the receiver's matching endpoint.
Keep component gene IDs, signed cKO/WT LFCs, padj, and tested/untested states.
Unresolved duplicate-symbol DE evidence is flagged and excluded from primary
DEG support, rather than selecting the smallest padj. Requiring both sides to
be DE is an additional descriptive flag, not the primary criterion.

## Evidence contract and validity gates

- Exact target membership and 12-mouse roster; no overlapping CA1/IN nuclei.
- Original source hashes unchanged before/after analysis; outputs only under a
  new run root created by native approved execution.
- Resource contents and installed LIANA version verified; no network during run.
- Sparse normalization; full-gene library denominators retained; no imputation.
- One-to-one native score joins; valid ranks and finite expression scores.
- Matching endpoint only; native DESeq2 seven-column exports stay unchanged.
- All group counts, mouse coverage, excluded resource components, ambiguous
  symbols, low-support and untested groups exported.
- Every plotted edge/row has an exact downloadable source-data table; readable
  top-candidate displays do not replace complete result exports.
- Independent review of methods, native outputs, figure rendering and report
  links before interpreting results.

## Claim boundary

These are predicted, transcript-compatible protein interactions and candidate
mechanisms. snRNA-seq does not directly measure secreted protein, surface
receptor availability, spatial contact, direction of physiological signalling,
activation/inhibition, synaptic connectivity or anatomical layers. Standard
protein LR resources incompletely cover neurotransmitter/metabolite biology.
Expression absence is therefore not proof of absent communication. DEG support
is reused gene-level evidence, not independent validation of an interaction.
No interaction is called significantly altered between genotypes without a
separately specified replicate-aware differential-interaction model.

## Implementation

A new local Snakemake wrapper runs preparation, six consensus fits plus
mouse-level expression support, table annotation and visualization sequentially.
The user selected reuse of the existing project controller. A clean separate
task environment and frozen reference cache are prepared outside execution.
The immutable workflow and concrete configuration are checked and planned in
NGS Analysis Workbench, then execution is requested through its native approval.
