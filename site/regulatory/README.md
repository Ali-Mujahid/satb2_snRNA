# Regulatory / independent-mouse explorer data contract

This additive UI displays completed, validated saved results and performs no
scientific modeling. With no `data.js` it displays **Awaiting verified analysis
outputs**. Synthetic tests live outside `site/`. The existing nuclei and strict
DEG-only LR explorers retain their original scientific inputs and results.

The SATB2 activity button selects the saved curated activity scores across all
compatible endpoints and both genotypes, with native units and automatic axes.
The full-prior primary analysis and the self-target-excluded sensitivity remain
separate. Network edge colors represent saved RNA-correlation signs or curated
regulatory-prior signs; genotype effects are reported in the saved WT–cKO tables.

## Loading and identifiers

Create `data.js` only from completed, validated saved outputs. Assign either
`window.SATB2_REGULATORY_DATA = payload` or
`window.SATB2_REGULATORY_DATA_GZIP = gzip_base64_encoded_JSON`.
The latter uses the existing local `../vendor/fflate.min.js` library. No network
API or runtime analysis is called. The single-file offline publication builder must
inline this payload, CSS, app and existing vendor libraries before distributing it.

IDs are exact opaque strings. The adapter may use endpoint IDs `CA1_DL`, `CA1_SL`,
`IN_C1`, `IN_C6`, `IN_C8`; retain the full original Seurat / DESeq2 endpoint ID in
`source_endpoint_id` and its exact membership hash. Endpoint `family` is `CA1` or
`IN`. Do not silently map or combine different memberships.

```typescript
interface Payload {
  schema: 'satb2-regulatory-atlas/v1';
  analysis_id: string;
  source_fingerprint: string;
  status: 'complete' | 'partial' | 'awaiting_analysis';
  created_at?: string;
  methods: Method[];
  endpoints: Endpoint[];
  samples: Sample[];
  features: Feature[];
  values: MouseValue[];
  effects?: SavedEffect[];
  hubs?: Hub[];
  nodes?: NetworkNode[];
  edges?: NetworkEdge[];
  embeddings?: Embedding[];
  coverage?: Coverage[];
  provenance?: object;
}
interface Method {
  id: string;
  label: string;
  family: 'hdWGCNA' | 'curated_tf_activity' | 'motif_regulon' | 'mouse_validation';
  status: 'complete' | 'partial' | 'awaiting_analysis' | 'failed' | 'not_available';
  score_label: string;     // e.g. the actual native module eigengene score name
  units: string;
  activity_unit: 'mouse' | 'nucleus' | 'mouse_and_nucleus';
  endpoint_ids?: string[]; // REQUIRED for shared, pooled-lineage network records
  network_edge_basis?: string; // e.g. weighted coexpression / curated TF-target
  interpretation?: string;
  eligibility_policy?: string;
  resources?: object;
  source_files?: object[]; // paths + SHA256, code, version, parameters, run IDs
}
interface Endpoint {
  id: string;
  label: string;
  family: 'CA1' | 'IN';
  source_endpoint_id: string;
  membership_sha256: string;
  identity_warning?: string;
}
interface Sample { id: string; genotype: 'WT' | 'cKO'; }
interface Feature {
  id: string; label: string; method_id: string;
  kind: 'module' | 'tf' | 'regulon' | 'gene' | 'lr_score';
  endpoint_id?: string; // include if separately fitted modules are endpoint-specific
  color?: string; units?: string; source_id?: string;
}
interface MouseValue {
  method_id: string; endpoint_id: string; feature_id: string; sample_id: string;
  value: number | null;
  eligible: boolean;
  missing_reason?: string;
  n_nuclei?: number | null;
  coverage_n?: number | null; coverage_total?: number | null;
  source_fit_id?: string; normalization?: string;
  // Additional original source fields are preserved in CSV / JSON.
}
interface SavedEffect {
  method_id: string; endpoint_id: string; feature_id: string;
  contrast: string; effect_name: string; effect: number | null;
  statistic?: number | null; pvalue?: number | null; padj?: number | null;
  tested: boolean; n_WT?: number | null; n_cKO?: number | null;
  reason?: string; source_fit_id?: string;
}
interface Hub {
  method_id: string; endpoint_id: string | null; feature_id: string; // feature = module
  gene_id: string; gene_symbol: string;
  kME: number | null; hub_rank: number | null;
  // Preserve native module membership, connectivity, q-values and source fields.
}
interface NetworkNode {
  id: string; // globally unique, including method + endpoint + module as needed
  method_id: string; endpoint_id: string | null; feature_id: string;
  label: string; gene_id?: string;
  kind: 'module_gene' | 'tf' | 'target';
  color?: string; kME?: number | null; hub_rank?: number | null;
}
interface NetworkEdge {
  method_id: string; endpoint_id: string | null; feature_id: string;
  source: string; target: string; // exact globally unique node IDs
  weight: number; sign?: 'positive' | 'negative' | 'unknown';
  directed: boolean; // FALSE for hdWGCNA coexpression; no arrows invented
  edge_basis: string; // weighted coexpression / curated prior / motif support
  motif_evidence?: string; source_id?: string;
}
interface Embedding {
  id: string; label: string; method_id: string; endpoint_id: string;
  basis: string;
  coordinates_kind: 'original_source';
  coordinates_sha256: string;
  nuclei: { id: string; sample_id: string; x: number; y: number; }[];
  activity: Record<string, (number | null)[]>; // feature_id -> exact nucleus axis
  source_score_description: string;
}
interface Coverage {
  method_id: string; endpoint_id: string; sample_id?: string; feature_id?: string;
  eligible: boolean; reason?: string;
  n_nuclei?: number | null; n_metacells?: number | null;
  measured_targets?: number | null; required_targets?: number | null;
  // Additional supplied coverage / warnings / model metadata stay exportable.
}
```

## Minimum hdWGCNA adapter

Use one method ID for each separately fitted lineage/network when appropriate.
Each module is a `Feature(kind='module')`. Supply one `MouseValue` per independent
mouse × exact endpoint × module, including missing / ineligible mice as explicit
rows with `value=null` and reasons. Do not package each metacell or nucleus as an
independent biological replicate. Separately fitted modules require distinct IDs;
similar colors or names do not imply equivalence.

For a network fitted once across a lineage, use `endpoint_id=null` on its hubs,
nodes and edges, and declare all compatible endpoint IDs in `Method.endpoint_ids`.
Shared module features omit `endpoint_id`. These network records are retained
once when any compatible endpoint is selected, never duplicated for each program
or cluster. The UI labels them **Pooled lineage network; endpoint selection does
not refit edges**. Mouse values and saved effects still require exact endpoint IDs.

Module membership / kME / hub rank tables map to `hubs`. Export the plotted hub
network's actual weighted edges as `edges` with `directed=false` and matching
nodes; a hub list alone is insufficient to invent connections. Native module
eigengene sign orientation and adjustment / harmonization choices belong in
`Method.interpretation` / `resources`. kME is module membership, not a TF activity
score, and neither is gene RNA expression.

Mouse-level model results map to `effects` without recomputing tests in JavaScript.
Use the actual contrast, effect units, cohort, FDR correction and model status.
If inference is unavailable, leave `effects=[]`; descriptive mouse plots remain
available. The genotype display filter never changes these original model rows.

## Display and export rules

- Unique `values` key: method + endpoint + feature + independent mouse. Unknown
  references, duplicate keys, non-finite values, unaligned vectors and artificial
  UMAP coordinates are rejected. Numeric missing values use JSON `null`.
- Native values are preserved. For mouse plots, display z = (value − mean) /
  sample SD, with reference fixed to **all eligible finite mice for that exact
  method + endpoint + feature**, including both genotypes. Fewer than two values
  or zero SD gives unavailable z, not zero. The transform is descriptive only.
- Optional nucleus UMAP display z uses all finite scores on that exact saved
  source axis, including both genotypes, and is labeled separately. No mouse
  inference is obtained from nucleus scores. Never substitute TF RNA expression
  for an unavailable activity vector.
- Heatmaps / dotplots retain missing cells. Boxes use eligible mice; finite
  ineligible points can optionally be shown open and never enter box summaries.
- CSV contains every filtered source row plus explicitly named display values;
  plot limits only reduce drawing. Settings and plot JSON include source identity.
  PDF / PNG reproduce the current appearance; SVG retains vector geometry.
- Curated activity and motif-supported regulons are separate method families.
  A missing TF method stays unavailable. The strict DEG-only communication atlas
  remains an independent, unchanged source with its existing three directions.

## Developer checks

Synthetic fixtures are under `tests/fixtures/regulatory_synthetic.cjs` and must
never be copied to `site/` or published as biological results. Run
`node tests/regulatory.cjs` for schema, eligibility, display, source-preservation
and plot/export-state checks. Browser acceptance must use the authorized live
HTTPS site after publication, not a workaround for blocked local-file browsing.
