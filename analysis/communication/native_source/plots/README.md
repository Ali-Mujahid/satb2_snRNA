# SATB2 CA1–inhibitory ligand–receptor figures

This renderer reads the completed LIANA evidence table. It does not fit LIANA,
CellPhoneDB, DESeq2, or any other model. It does not change an existing analysis.
Use a new output directory for each completed rendering; the script refuses to
overwrite a finished report.

```bash
python workflow_source/plots/plot_ccc.py \
  --results /absolute/new_run/results/tables/interaction_evidence.tsv.gz \
  --outdir /absolute/new_run/reports/communication_atlas
```

Dependencies: Python, NumPy, pandas and matplotlib from the isolated LIANA
environment. No JavaScript packages, R packages, remote assets or web server are
needed. Open `index.html` locally. Keep its `figures`, `data` and `README.md`
siblings together when copying the report to another computer.

The required `condition_pair_eligibility.tsv` is read from beside the evidence
table (or pass `--pair-eligibility /absolute/file.tsv`). It explicitly records
all evaluated pairs, including pairs returning no ligand–receptor rows. This
prevents untested pairs from being mistaken for measured zero counts.

## Populations and biological scope

- Primary: CA1 Deep and CA1 Superficial expression programs versus the Seurat
  inhibitory clusters.
- Supplementary: the same CA1 programs versus inhibitory total.
- Each scope includes CA1 → inhibitory, inhibitory → CA1, and both directions.
- WT, cKO, and pooled expression contexts are shown separately.
- Deep and Superficial denote molecular programs. Physical layers were not
  measured. Group annotations and identity-review flags remain in the exported
  data; an asterisk identifies a provisional or mixed group on the graphs.

These are transcript-supported *candidate* interactions. Nuclear RNA does not
establish the availability of secreted ligand or receptor protein, spatial
contact, functional connectivity, signaling activity, activation, or inhibition.
Sender/receiver direction is not an activating/inhibitory functional assignment.

## Evidence and ranking

The default display requires `magnitude_rank <= 0.05` and endpoint-specific
DESeq2 support. Lower LIANA magnitude rank is higher priority. Rank is **not a
genotype-comparison P value**. CellPhoneDB permutation P values, when present in
the full evidence table, are also not tests of the WT–cKO difference.

The `DEG_supported` label is defined upstream: at least one ligand subunit has
DESeq2 support in its matched source population, **or** at least one receptor
subunit has support in its matched target population. This prioritizes a
candidate; it does not constitute a statistically significant differential
interaction. Both significant and nonsignificant gene-level endpoint evidence
must remain available in the upstream audit. The renderer accepts the supplied
boolean and never substitutes a union of DEGs across unrelated populations.

The third view, `replicated_specific` (**Replicated + context-specific DEG
candidates**), is an additional exploratory priority set. It requires all of:

- endpoint-specific DESeq2 support;
- magnitude rank ≤ 0.05 (or a stricter configured display cutoff);
- LIANA specificity rank ≤ 0.05;
- native CellPhoneDB within-condition permutation P ≤ 0.05;
- at least 3 supporting independent WT mice **and** 3 supporting cKO mice;
- every inhibitory endpoint has the exact review flag `REFERENCE_SUPPORTED`.

The CA1 program endpoint is retained under its existing expression-program
definition. This view does not resolve uncertain identities, fit any additional
model, estimate a differential interaction P value, or control genotype FDR.
The specificity and CellPhoneDB criteria describe the expression context; they
are **not** WT–cKO tests. The original `DEG_supported` and `all_candidates`
views remain unchanged, retain all 20 inhibitory groups including mixed/review
groups, and remain available in the selector. The original DEG-supported view
stays the default. The stringent view can be empty for inhibitory total if its
combined identity-review field is not exactly `REFERENCE_SUPPORTED`; this is a
selection exclusion, not evidence of absent communication. Every rule is also
recorded in the manifest.

An additional context view shows every expression-supported candidate passing
the same rank criterion. Candidates missing or below this criterion in one
condition must not be called gains/losses of biological communication: ranks
depend on the expression context, cell abundance and the tested interaction set.
No differential interaction model is fitted by this renderer.

Native LIANA `lr_logfc` is a group-versus-rest expression score. It is **not**
the DESeq2 cKO/WT log2 fold change. Do not interpret its sign as a genotype effect.

## Figures

All figure types are supplied as PDF, PNG and editable-text SVG. Every plot has
an adjacent `_SourceData.tsv` with exactly its plotted values, plus a manifest
record with its selection rules. The interface links to these files directly.

### Dotplot

Rows are the top 24 ligand–receptor identities ranked by their best magnitude
rank across all three conditions within the selected scope/direction/evidence
view. This is a descriptive prioritization, not a genotype test. Ligand–receptor
order and cell-pair order are fixed across conditions. Cell pairs are split into
pages of at most 20 columns; every retained pair is available in a later page.

Color is magnitude rank, on the same `[0, rank cutoff]` scale in every condition.
Dot area encodes the **minimum detection fraction across every subunit of both
the ligand complex in its source and receptor complex in its target**, with a
common size scale. The preferred source fields are
`ligand_all_subunit_min_prop` and `receptor_all_subunit_min_prop` (and their
minimum is the true complex-level `lr_prop`). Native `ligand_props` and
`receptor_props` remain unchanged in the audit: those refer to subunits selected
by the native mean-expression bottleneck and may miss a less detected subunit.
If the full-complex fields are absent, the renderer falls back explicitly to
native endpoint-subunit fractions, changes the legend, and records this weaker
definition in the plot data and manifest. A tiny base area keeps present rows
visible; the exact formula and fractions are exported. Missing
dots are absent or below the criterion, not measured expression zeros. This
selection is limited for legibility; full qualifying and unfiltered tables
remain available.

### Source–target tileplot

Each tile counts distinct predicted ligand–receptor identities for its tested
source/target direction under the selected criteria. Counts are not nuclei,
synapses, or statistically significant genotype changes. Same-family and
unrequested directions are gray/NA. A numeric zero means no candidate passed
the criteria for that tested direction. All source/target groups remain on the
axes, including zero-edge groups. WT, cKO and pooled use a shared count scale.

### Circos / chord diagram

These are actual closed, finite-width ribbons between cell-group sectors.
Blue ribbons originate in CA1; orange ribbons originate in inhibitory neurons.
Tips mark receptor recipients. Node order is fixed and every group is shown.
Equal sector sizes improve legibility. Ribbon width is proportional to the
candidate count **within its endpoint sector**; widths across different sectors
do not share an absolute count scale. The source data contain exact counts and
arc coordinates. This graph aggregates ligand–receptor identities at the group
level; use the dotplots and tables to identify the underlying genes/complexes.

### Directed circle network

Arrows run from ligand source to receptor recipient. Width is candidate count;
node area is the sum of incoming and outgoing candidate counts. Node positions,
width scales and size scales are fixed across WT, cKO and pooled contexts.
Isolated groups are retained. Arrow colors identify source family, not
activation/inhibition. Exact per-edge counts and node positions/sizes are
exported separately.

## Tables

- `data/interaction_evidence.tsv.gz`: full completed evidence table, copied
  without altering its bytes if its input was gzip compressed.
- `data/node_annotations.tsv`: complete group names, supplied annotations and
  identity-review flags.
- `data/condition_pair_eligibility.tsv`: explicitly evaluated directions and
  source/target nucleus counts for every condition, including zero-return pairs.
- `data/plot_manifest.json`: figure selection, all condition/direction totals,
  plotting parameters, and input SHA-256 for reproducibility.
- Each figure folder has `qualifying_interactions.tsv.gz`, including the full
  supplied evidence columns, and `source_target_counts.tsv`.
- Each folder also has `liana_rank_aggregate_7columns.tsv`, with exactly the
  requested columns: `source`, `target`, `ligand_complex`, `receptor_complex`,
  `lr_logfc`, `lr_gmeans`, and `magnitude_rank`. Condition, scope, direction and
  selection are encoded by its folder; use the full table for provenance.

## Input contract

One row per scope/condition/source/target/ligand/receptor. Required columns:
`scope`, `condition`, `source`, `target`, `ligand_complex`, `receptor_complex`,
`lr_logfc`, `lr_gmeans`, `magnitude_rank`, `direction`, `deg_supported`,
`ligand_props`, `receptor_props`. The renderer preserves optional evidence
columns including sample support and eligible-mouse denominators. Fractions
must be 0–1, ranks 0–1 or NA, and DEG support explicitly boolean.

Recognized scopes are `IN_clusters` and `IN_total`; condition labels are `WT`,
`cKO` and `pooled`; directions are `CA1_to_IN` or `IN_to_CA1`. Canonical groups
are `CA1_Deep`, `CA1_Superficial`, `IN_C1` … `IN_C20`, and `IN_Total`. Optional
`source_label`, `target_label`, `source_review`, `target_review` improve annotation
reporting. `REFERENCE_SUPPORTED` does not receive a review asterisk; any other
nonempty substantive review flag is retained and marked.

The output contains static precomputed figures. HTML controls select among
saved panels; they do not recalculate ranks, expression, detection thresholds,
complex aggregation, or differential statistics.
