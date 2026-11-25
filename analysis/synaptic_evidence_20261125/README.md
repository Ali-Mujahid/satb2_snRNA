# SATB2 CA1 ligand–receptor candidates and synaptic references

The Pages entry is `site/synaptic_evidence/`. Machine-readable results are in `analysis/synaptic_evidence_20261125/results/`. The standalone HTML includes the data, plot code and PDF renderer. The atlas loads a checksum-verified copy from its tagged GitHub release revision to respect the existing Pages size budget. The DEG inputs are the original Seurat mouse pseudobulk results for deep C1+C5 and inhibitory C1–C20, plus the superficial C2+C3 refit after removal of CA1 C8. C18 has no eligible original model. The prior pooled LIANA file for these exact populations supplies the candidates. Each molecule passes adjusted P < 0.05 at its own endpoint. WT and cKO mouse-support fractions require ≥5 nuclei of each endpoint and ≥5% detection of each component in its own population, within mouse. Colors in the Connectome-seq and candidate RNA plots use the local genotype contrast.

## Inputs and execution

Inputs named below are stored under the workstation analysis directory `SATB2_CA1_SYNAPTIC_PROTEOMICS_20261125/sources/`. Obtain the source files from the authors. The output manifest records checksums for each source and each fitted model table used.

* [SynDive Supplementary Data 2](https://media.springernature.com/full/springer-static/esm/art%3A10.1038%2Fs41467-025-63119-5/MediaObjects/41467_2025_63119_MOESM4_ESM.xlsx). Both strata protein sheets are used. Split grouped symbols on semicolons. Plot the median of finite reported values from single-gene protein groups with unambiguous symbol-to-ID mapping. Author flags and individual groups are retained. FASS compares each stratum against the other strata. Tissue compares it against mean log₂ abundance across all four. These measurements are from the reference animals.
* [Synaptome.db 0.99.17](https://bioconductor.posit.co/packages/3.22/data/annotation/html/synaptome.db.html) and [database archive](https://datashare.ed.ac.uk/download/DS_10283_4757.zip). Extract the official R sources into `sources/synaptome.db` and the SQLite file (`AH107282`) into `sources`. The R driver uses the supplied `getGeneInfoByName` and `getPPIbyIDs` against the extracted database. The Python script validates the matched study and PPI keys against those native results and queries the full DEG set. Proteomic entries use the source’s mouse-symbol mappings and retain experimental species and region. Pair assay records use the database PPI type, method and PMID. PPI anatomical scope is unrestricted.
* Save the G2C [Bayés mouse–human PSD](https://www.genes2cognition.org/publications/human-mouse-psp/) and [Fernández PSD95](https://www.genes2cognition.org/publications/tap-psd-95/) pages as `human-mouse-psp.html` and `tap-psd-95.html`. The identity mapping uses stable mouse Ensembl IDs or reported mouse symbols and indexed MGI identifiers. Table S1 species calls and S2d/S1 complex membership flags preserve listed, ambiguous and unmarked results. Coverage is bounded by these publications. The human column follows the authors’ cross-species assignments.
* [Connectome-seq Extended Data Fig. 9b](https://www.nature.com/articles/s41592-026-03026-9/figures/15) provides the 30 genes, transcribed in displayed order from the [source image](https://media.springernature.com/full/springer-static/esm/art%3A10.1038%2Fs41592-026-03026-9/MediaObjects/41592_2026_3026_Fig15_ESM.jpg) saved as `connect_fig15.jpg`. The published comparison is between Purkinje clusters with different fractions of pons-connected cells. Local gene-set tests use the uniquely named genes with finite per-gene adjusted P in each population. The universe, locally testable panel, DEG counts, local direction, enrichment P and Benjamini–Hochberg adjusted P across all fitted populations are in the test table. Shared genes can occur in multiple endpoints. Protein preparations can also recur between reference databases.
* Local model paths and expected hashes are read from `SATB2_CA1_IN_EXPANDED_20261102/workflow/config/analysis.json`. The exact C8-excluded superficial model and strict RNA candidates come from `SATB2_CA1_EPHYS_PLANNING_20261125`.

From the workspace root, using Python with openpyxl and lxml, and R with DBI, RSQLite, dplyr and dbplyr installed, run

```bash
project_dir="$PWD/SATB2_CA1_SYNAPTIC_PROTEOMICS_20261125"
analysis_dir="$PWD/satb2_snRNA_github/analysis/synaptic_evidence_20261125"
(cd "$project_dir" && Rscript "$analysis_dir/synaptome.R")
python3 "$analysis_dir/analyze.py" --workspace "$PWD" --sources "$project_dir/sources" --output "$analysis_dir/results"
python3 "$analysis_dir/build_reader.py"
```

A missing reference estimate or database entry leaves detection unresolved. Co-membership indicates a shared preparation. RNA changes, protein presence and a molecular interaction assay each address different biological questions. These nuclear RNA data have no paired synaptic barcodes.
