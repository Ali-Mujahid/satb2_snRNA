# Ligand–receptor atlas data contract

The communication explorer packages completed LIANA 2.0.0 results from registered run `d397e9d9-03bc-43b0-9bd8-901409248fd6`. Packaging reads saved outputs. It performs no expression transformation, annotation, DESeq2 fit or ligand–receptor model fit.

## Offline interactive transport

`site/communication/data/data.js` assigns `window.CCC_DATA_GZIP`. Decode base64, gunzip with the already bundled `fflate`, decode UTF-8, then parse JSON. There is no fetch requirement for this dataset, so saved Pages files can also work from an offline folder.

The decoded schema is `satb2-liana-atlas/v1`:

- `columns`: ordered names for the 48 scalar values in each row array; original native scalar columns are retained.
- `datasets['IN_clusters|pooled']` (and WT/cKO or IN_total): `rows`, `eligible_pairs`, `scope`, `condition`.
- `nodes`: exact population ID, display label, native annotation and identity-review string, DESeq2 endpoint/test status and saved DEG count. IN C18 remains untested/missing.
- `components[population_id][complex_name]`: exact native array of per-subunit matching gene DESeq2 evidence, including stable IDs, signed log2(cKO/WT), adjusted P, ambiguity and test status. Use source + ligand_complex for ligand evidence and target + receptor_complex for receptor evidence. Arrays are verified invariant across the six fitted contexts.
- `provenance`: immutable native method/source identity, normalization, cohort, parameters and score meanings.

Scope and condition are dataset metadata, rather than repeated inside each row. The original two JSON component columns are replaced by the exact dictionary lookup. This removes repetition without discarding endpoint evidence. Numeric values retain the original floating-point precision; missing values become null. Original native tables preserve their exact scientific bytes and textual precision.

## Dataset counts and validation

| Resolution | pooled | WT | cKO |
|---|---:|---:|---:|
| IN_clusters | 20,471 | 19,467 | 21,680 |
| IN_total | 1,108 | 1,040 | 1,170 |

There are 64,936 exact unique `(scope, condition, source, target, ligand_complex, receptor_complex)` keys. Each cluster-context dataset has 80 eligible directed pairs and each total-context dataset has four.

At magnitude_rank≤0.05, the pooled cluster DEG-supported subset has 1,589 records. Its stricter reference-supported, replicated and context-specific subset has 324 records: 112 CA1→IN and 212 IN→CA1, with two records having gene-level DEG evidence at both ends. Packaging asserts these counts against the unchanged native evidence.

The strict tier additionally requires specificity_rank≤0.05, cellphone_pvals≤0.05, sample_support_WT≥3, sample_support_cKO≥3 and both endpoint review strings `REFERENCE_SUPPORTED`. This saved policy gives an empty strict IN_total view because heterogeneous IN_Total is not an atomic reference-supported cluster. Do not silently override that result.

`lr_prop` is the minimum detection proportion across every subunit of both complexes. The transport preserves tiny sparse-summation overshoots (maximum 1.0000000000000093). Clamp to [0,1] only when calculating plot areas. Native `ligand_props` and `receptor_props` use mean-expression bottleneck subunits and can differ from the all-subunit minimum.

`eligible_pairs` is essential: a tested eligible pair with no retained records can be shown as zero; an untested direction is NA. IN→IN and CA1→CA1 were not fitted here. IN_Total and its constituent IN clusters belong to different models and must never be combined as one fit.

## Meaning of statistical fields

`lr_logfc` measures native cell-group-versus-rest specificity. It is not cKO/WT log2 fold change. `magnitude_rank` and `specificity_rank` are native consensus priorities, not genotype P values or FDR. `cellphone_pvals` are within-context cell-label specificity probabilities from LIANA's native CellPhoneDB scorer, not a standalone CellPhoneDB v5 run and not genotype tests. A native zero permutation P means below empirical 1/1000 resolution.

Only matching existing mouse-level gene DESeq2 effects use signed log2(cKO/WT). No replicate-aware differential-interaction model was fitted. WT/cKO differences in counts, ranks or scores are descriptive. Reusing gene DESeq2 evidence does not independently validate an interaction. Nuclear RNA and curated priors do not measure contact, surface protein, physiological direction or activation/inhibition. Molecular deep/superficial programs are not measured anatomy.

## Native downloads and source provenance

Native `communication_data/full_result_tables` and `communication_data/existing_DESeq2_exports` live outside the Pages tree because the existing atlas is close to its Pages size ceiling. Pages download indexes link to the immutable `atlas-v3.5.0` repository tag. Those links require internet; the full portable release ZIP contains all native tables and 276 original PDF/PNG/SVG panel sets for offline use.

Every requested native interaction table has exactly `source, target, ligand_complex, receptor_complex, lr_logfc, lr_gmeans, magnitude_rank`. Matching gene DESeq2 exports retain their original standard seven columns and symbols. Extended evidence, all eligible/supporting mouse records and component stable-ID audits are also preserved.

The native full-gene component JSON is transported as `.json.gz`; gunzip restores its original 76,645,727 bytes exactly. `site/communication/data/TRANSPORT_RECEIPT.json` records source and transport byte lengths, SHA256 hashes and the only adapted native download-index link. All scientific native files are unchanged. A compact `manifest.json` records the interactive JSON/gzip/script hashes, schema and validation counts.

Executed scientific source/configuration, original model-source provenance, runtime lock/audit files and the frozen mouse resource are copied under `analysis/communication`. Input hashes identify original expression and Seurat membership inputs; the large expression matrices are not duplicated into this publication.

## Regenerate the transport

Use Python 3.9 or later; packaging uses only its standard library:

```bash
python3 scripts/package_communication_atlas.py \
  --source /absolute/path/to/SATB2_communication_report \
  --repo /absolute/path/to/satb2_snRNA_github
```

This command rewrites only communication publication transport, copied native downloads, provenance and interpretation artifacts. It never invokes scientific tools. The local project used its existing Python launcher for this packaging command. The default immutable download tag is specified by `IMMUTABLE_DATA_BASE` in the packaging source.
