# Grm7 RNA and isoform-feature evidence: twelve P60 mice

Open **index.html** in Chrome or Safari after extracting the complete folder. It works offline. Do not open it inside the ZIP archive.

The analysis includes six WT and six Satb2 cKO mice. Superficial CA1 is C2+C3; deep CA1 is C1+C5. C8 is excluded from the superficial group. The frozen Seurat assignments, original matrices, clusters and DESeq2 analyses were preserved.

## Start here

- `ANALYSIS_REPORT.md`: findings, limitations, methods and next experiments.
- `index.html`: per-mouse heatmaps/dotplots; select population, evidence type and normalization. Appearance controls are under “Appearance and feature selection.” SVG, PNG, PDF, TSV and plot-state JSON downloads are available.
- `figures/`: 21 figures, each supplied as SVG, PNG and PDF. `01` shows transcript structures; `02` shows gene counts; `03` conservative transcript-specific splicing; `04` diagnostic features; `05` exons; `06` junctions; `07` ambiguity; `08` presynaptic receptor schematic; `09` A/B terminal splicing.

## Read the counts correctly

Five transcripts are annotated in the actual mm10-2020-A mapping reference. Only Grm7-201 and Grm7-202 are complete receptor models. Grm7-203 and Grm7-205 have incomplete coding starts; Grm7-204 is annotated as nonsense-mediated decay. Transcript-specific splice support is sparse: these data do not reliably quantify complete isoform abundances or demonstrate a WT/cKO isoform switch.

`gene_UMI_counts_per_mouse.tsv` contains BAM-derived gene support and separately defined evidence categories. `isoform_evidence_per_mouse.tsv` separates diagnostic features, conservative unique splice support and overlapping compatible upper bounds. An exon-only diagnostic signal does not establish a full-length mature transcript. `terminal_splice_event_counts_per_mouse.tsv` counts four defined A/B terminal events directly. The B-type exon-entry event is shared by Grm7-201 and Grm7-205.

Exon, junction, CDS, UTR and disjoint exonic-bin files provide feature-level counts. Features can overlap within a molecule and must not be summed to reconstruct gene totals. Counts use corrected barcode/UMI pairs. Raw counts and three descriptive normalizations are supplied; these are not new DESeq2 fits.

The coverage grid has 420 mouse/population combinations. Nine combinations contain no sampled nuclei and are marked `NOT_SAMPLED`; they are not biological zeros. Zero splice support in a sampled population also does not prove absence of a transcript.

## Provenance

`VALIDATION.json` records the independent count checks. `EXPLORER_QA.json` records browser/download checks. `COMPLETE.json` records the corrected native run and source hashes. `reference_run_binding.json` records the mapping-reference checks and historical cKO_933 reference-path limitation. `Grm7.reference.gtf`, `annotation.json`, `transcripts.tsv` and `features.tsv` define the models and features used.

The reviewed per-million values use the full original library (`original_nCount_RNA`), verified against all-gene matrix sums. The native run initially used the saved object's HVG-only `nCount_RNA` subtotal. `NORMALIZATION_CORRECTION.json` and `library_denominator_audit.tsv` document this export correction; raw and per-nucleus counts are unchanged. Prefer these reviewed exports over the native run's original per-million field.

`source/` contains the frozen corrected counting workflow, configuration, tests and reviewed export code. Input BAMs and the full nucleus-barcode roster are not included. These remain at their original locations; scripts do not embed credentials. The first, superseded transcript-compatibility run is intentionally excluded from this package.

No differential-isoform hypothesis test, full-length transcript reconstruction, protein measurement or direct ligand–receptor binding assay was performed.
