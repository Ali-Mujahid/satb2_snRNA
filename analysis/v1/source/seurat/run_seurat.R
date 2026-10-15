#!/usr/bin/env Rscript
# All-twelve original-count comparison. Authoring/parse inspection alone does not run it.
# Usage: Rscript --vanilla run_seurat.R CONFIG.json STAGE
# STAGE: log_global | sct_global | annotate_subsets | markers | pseudobulk
# All outputs go to a fresh stage directory; prior analyses are read-only inputs.

fail <- function(...) stop(paste0(..., collapse = ""), call. = FALSE)
`%or%` <- function(x, y) if (is.null(x)) y else x
read_tsv <- function(p) utils::read.delim(p, check.names = FALSE, stringsAsFactors = FALSE)
truth <- function(x) !is.na(x) & tolower(as.character(x)) %in% c("true", "1", "yes")
utc <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
safe <- function(x) gsub("[^A-Za-z0-9_-]", "_", x)
write_tsv <- function(x, p) {
  if (file.exists(p)) fail("Refusing to overwrite ", p)
  utils::write.table(x, p, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}
save_new <- function(x, p) {
  if (file.exists(p)) fail("Refusing to overwrite ", p)
  saveRDS(x, p, compress = FALSE)
}
new_dir <- function(p) {
  if (dir.exists(p)) {
    # Snakemake can create the parent of COMPLETE.json before invoking a rule.
    # Only that genuinely empty directory is reusable; hidden files, failed
    # receipts and partial outputs all count as content and remain protected.
    if (length(list.files(p, all.files = TRUE, no.. = TRUE))) fail("Output directory contains prior/partial evidence; use a new run root: ", p)
  } else if (!dir.create(p, recursive = TRUE)) fail("Cannot create ", p)
  p
}
check_unique <- function(x, what) {
  if (anyNA(x) || any(!nzchar(x)) || anyDuplicated(x)) fail("Invalid unique IDs: ", what)
}
roster <- c(wt_1122 = "WT", wt_1123 = "WT", wt_1178 = "WT", wt_1220 = "WT",
  wt_1226 = "WT", wt_1228 = "WT", cKO_1138 = "cKO", cKO_1139 = "cKO",
  cKO_1160 = "cKO", cKO_1165 = "cKO", cKO_933 = "cKO", cKO_991 = "cKO")

snn_to_igraph_sparse <- function(snn) {
  if (!inherits(snn, "sparseMatrix") || nrow(snn) != ncol(snn)) fail("SNN must be square and sparse")
  if (!identical(rownames(snn), colnames(snn))) fail("SNN nucleus axes differ")
  check_unique(colnames(snn), "SNN nuclei")
  if (anyNA(snn@x) || any(!is.finite(snn@x)) || any(snn@x < 0)) fail("Invalid SNN edge weights")
  # Seurat 5.5.1 RunLeiden's igraph path requests mode='lower'. In igraph
  # 2.3.3, ensure_no_na first indexes x[lower.tri(x)], allocating a dense
  # N-by-N logical mask despite a sparse SNN. Extract the identical weighted
  # lower triangle through Matrix instead; remove self edges and retain every
  # vertex, including isolates. Memory is proportional to stored SNN edges.
  edges <- Matrix::summary(Matrix::tril(snn, k = -1L))
  edges <- edges[edges$x != 0, , drop = FALSE]
  graph <- igraph::make_empty_graph(n = ncol(snn), directed = FALSE)
  if (nrow(edges)) graph <- igraph::add_edges(graph,
    edges = as.vector(t(as.matrix(edges[, c("i", "j"), drop = FALSE]))),
    attr = list(weight = as.numeric(edges$x)))
  igraph::set_vertex_attr(graph, "name", value = colnames(snn))
}

sparse_leiden_membership <- function(graph, snn, resolution, random_seed, n_iter = 10L) {
  if (igraph::vcount(graph) != ncol(snn) ||
      !identical(igraph::vertex_attr(graph, "name"), colnames(snn))) fail("Leiden vertex roster changed")
  # Passing an igraph object bypasses only the unsafe adjacency conversion.
  # Native Seurat RunLeiden still performs seeded modularity Leiden and its
  # cluster size ordering. Native GroupSingletons preserves FindClusters'
  # default postprocessing on the original sparse SNN.
  ids <- Seurat::RunLeiden(object = graph, leiden_method = "igraph",
    partition.type = "RBConfigurationVertexPartition", leiden_objective_function = "modularity",
    resolution.parameter = resolution, random.seed = random_seed, n.iter = n_iter)
  names(ids) <- colnames(snn)
  ids <- get("GroupSingletons", envir = asNamespace("Seurat"))(
    ids = ids, SNN = snn, group.singletons = TRUE, verbose = TRUE)
  levels <- unique(as.character(ids))
  is_int <- grepl("^[0-9]+$", levels)
  levels <- c(levels[is_int][order(as.integer(levels[is_int]))], sort(levels[!is_int]))
  factor(ids, levels = levels)
}

drop_scaled_layer <- function(obj) {
  if (!inherits(obj[["RNA"]], "Assay5")) fail("Expected Seurat v5 RNA assay for layer cleanup")
  # Assay5's official NULL setter deletes the layer and its feature/cell maps.
  # An empty matrix is not a deletion and fails the feature-overlap guard.
  if ("scale.data" %in% SeuratObject::Layers(obj[["RNA"]])) obj[["RNA"]]$scale.data <- NULL
  if ("scale.data" %in% SeuratObject::Layers(obj[["RNA"]])) fail("Scaled layer was not removed")
  methods::validObject(obj)
  obj
}

load_contract <- function(cfg) {
  required <- c("sample_manifest_tsv", "roster_tsv", "gene_manifest_tsv",
    "fixed_features_tsv", "output_root", "resource_budget_gib")
  if (!all(required %in% names(cfg))) fail("Incomplete configuration")
  sm <- read_tsv(cfg$sample_manifest_tsv)
  if (!all(c("sample", "condition", "sequencing_partition") %in% names(sm))) fail("Sample manifest schema")
  check_unique(sm$sample, "samples")
  if (!setequal(sm$sample, names(roster))) fail("All twelve confirmed mice required")
  sm <- sm[match(names(roster), sm$sample), , drop = FALSE]
  if (any(sm$condition != unname(roster)) || anyNA(sm$sequencing_partition) ||
      any(!nzchar(sm$sequencing_partition))) fail("Sample genotype/partition contract")
  md <- read_tsv(cfg$roster_tsv)
  needed <- c("nucleus_id", "barcode", "sample", "condition", "nCount_RNA", "passes_primary_QC", "call_new")
  if (!all(needed %in% names(md))) fail("Frozen roster schema")
  check_unique(md$nucleus_id, "nuclei")
  if (!setequal(unique(md$sample), names(roster)) || anyNA(md$condition) ||
      any(md$condition != unname(roster[md$sample]))) fail("Frozen animal roster mismatch")
  if (any(!truth(md$passes_primary_QC)) || any(md$call_new != "singlet")) fail("Roster includes nuclei outside agreed primary QC/singlets")
  if (nrow(md) != as.integer(cfg$expected_nuclei %or% 108062L)) fail("Frozen nucleus count changed")
  if (any(md$nucleus_id != paste0(md$sample, "_", md$barcode))) fail("sample_barcode identity contract")
  md$sequencing_partition <- sm$sequencing_partition[match(md$sample, sm$sample)]
  old <- grepl("^(historical_|legacy_|old_)", names(md)) |
    names(md) %in% c("cluster", "seurat_clusters", "EX_cluster", "IN_cluster", "CA1_cluster", "layer", "layer_label")
  md <- md[, !old, drop = FALSE]
  rownames(md) <- md$nucleus_id
  genes <- read_tsv(cfg$gene_manifest_tsv)
  if (!all(c("gene_id", "gene_symbol") %in% names(genes))) fail("Gene manifest schema")
  check_unique(genes$gene_id, "stable genes")
  if (any(grepl("[_|]", genes$gene_id))) fail("Stable gene IDs would be renamed by Seurat")
  genes$symbol_missing <- is.na(genes$gene_symbol) | !nzchar(genes$gene_symbol)
  genes$display_symbol <- ifelse(genes$symbol_missing, genes$gene_id, genes$gene_symbol)
  hv <- read_tsv(cfg$fixed_features_tsv)$gene_id
  check_unique(hv, "common HVGs")
  if (length(hv) != as.integer(cfg$nfeatures %or% 3000L) || any(!hv %in% genes$gene_id)) fail("Common HVG contract")
  count_source <- cfg$count_source %or% "original_full"
  if (!count_source %in% c("original_full", "paired_original", "paired_cellbender")) fail("Unrecognized count_source")
  list(samples = sm, metadata = md, genes = genes, features = hv, count_source = count_source)
}

# Read the explicit modern 10x /matrix group. CellBender adds root-level latent
# groups; they are not expression matrices and must never be enumerated as such.
read_10x_matrix_csc <- function(path) {
  if (!requireNamespace("hdf5r", quietly = TRUE)) fail("hdf5r is required for explicit CSC import")
  if (!requireNamespace("Matrix", quietly = TRUE)) fail("Matrix is required for explicit CSC import")
  f <- hdf5r::H5File$new(path, mode = "r")
  on.exit(f$close_all(), add = TRUE)
  required <- c("matrix/data", "matrix/indices", "matrix/indptr", "matrix/shape",
    "matrix/barcodes", "matrix/features/id", "matrix/features/feature_type")
  if (!all(vapply(required, f$exists, logical(1)))) fail("Source lacks the explicit 10x /matrix CSC schema: ", path)
  shape <- as.numeric(f[["matrix/shape"]]$read())
  if (length(shape) != 2L || any(!is.finite(shape)) || any(shape <= 0) ||
      any(shape != floor(shape)) || any(shape > .Machine$integer.max)) fail("Invalid CSC dimensions: ", path)
  ids <- as.character(f[["matrix/features/id"]]$read())
  barcodes <- as.character(f[["matrix/barcodes"]]$read())
  kinds <- as.character(f[["matrix/features/feature_type"]]$read())
  check_unique(ids, "H5 stable feature IDs"); check_unique(barcodes, "H5 barcodes")
  if (length(ids) != shape[1L] || length(barcodes) != shape[2L] ||
      length(kinds) != shape[1L] || anyNA(kinds) || any(kinds != "Gene Expression")) fail("CSC feature/barcode/GEX axis mismatch: ", path)
  values <- as.numeric(f[["matrix/data"]]$read())
  indices <- as.numeric(f[["matrix/indices"]]$read())
  indptr <- as.numeric(f[["matrix/indptr"]]$read())
  if (length(indices) != length(values) || length(indptr) != shape[2L] + 1L ||
      anyNA(values) || any(!is.finite(values)) || any(values < 0) || any(values != floor(values)) ||
      anyNA(indices) || any(indices < 0) || any(indices >= shape[1L]) || any(indices != floor(indices)) ||
      anyNA(indptr) || indptr[1L] != 0 || tail(indptr, 1L) != length(values) ||
      any(diff(indptr) < 0) || any(indptr != floor(indptr)) || any(indptr > .Machine$integer.max)) fail("Invalid nonnegative integer CSC source arrays: ", path)
  x <- methods::new("dgCMatrix", i = as.integer(indices), p = as.integer(indptr),
    x = values, Dim = as.integer(shape), Dimnames = list(ids, barcodes))
  if (!isTRUE(methods::validObject(x))) fail("Invalid imported CSC matrix: ", path)
  x
}

# Stream one complete declared capture at a time, then select exact frozen nuclei.
read_capture <- function(sample, contract, ids = NULL) {
  row <- contract$samples[match(sample, contract$samples$sample), , drop = FALSE]
  md <- contract$metadata[contract$metadata$sample == sample, , drop = FALSE]
  if (!is.null(ids)) md <- md[md$nucleus_id %in% ids, , drop = FALSE]
  if (!nrow(md)) return(NULL)
  if ("filtered_h5" %in% names(row) && !is.na(row$filtered_h5) && nzchar(row$filtered_h5)) {
    sha_field <- intersect(c("source_h5_sha256", "filtered_h5_sha256"), names(row))
    if (!length(sha_field) || is.na(row[[sha_field[1L]]]) || !grepl("^[0-9a-fA-F]{64}$", row[[sha_field[1L]]])) fail("Frozen source H5 SHA256 required: ", sample)
    source_path <- row$filtered_h5
    actual_sha256 <- digest::digest(file = source_path, algo = "sha256")
    if (!identical(tolower(row[[sha_field[1L]]]), actual_sha256)) fail("Source H5 changed since freeze: ", sample)
    x <- read_10x_matrix_csc(row$filtered_h5)
    j <- match(md$barcode, colnames(x))
  } else if ("counts_rds" %in% names(row) && !is.na(row$counts_rds) && nzchar(row$counts_rds)) {
    if (!"counts_rds_sha256" %in% names(row) || is.na(row$counts_rds_sha256)) fail("Frozen count RDS SHA256 required: ", sample)
    source_path <- row$counts_rds
    actual_sha256 <- digest::digest(file = source_path, algo = "sha256")
    if (!identical(tolower(row$counts_rds_sha256), actual_sha256)) fail("Count RDS changed since freeze: ", sample)
    x <- readRDS(row$counts_rds)
    j <- match(md$nucleus_id, colnames(x))
  } else fail("No original H5 or accepted count RDS for ", sample)
  if (!inherits(x, "sparseMatrix")) fail("Original counts must be sparse: ", sample)
  if (anyNA(j) || !setequal(rownames(x), contract$genes$gene_id)) fail("Feature/barcode source mismatch: ", sample)
  x <- methods::as(x[contract$genes$gene_id, j, drop = FALSE], "dgCMatrix")
  colnames(x) <- md$nucleus_id
  if (anyNA(x@x) || any(!is.finite(x@x)) || any(x@x < 0) || any(x@x != floor(x@x))) fail("Noninteger count source; no rounding authorized: ", sample)
  lib <- as.numeric(Matrix::colSums(x))
  if (any(lib <= 0)) fail("Zero full count library in audited usable roster: ", sample)
  if (contract$count_source != "paired_cellbender" && !identical(lib, as.numeric(md$nCount_RNA))) fail("Full original library totals mismatch: ", sample)
  if (contract$count_source == "paired_cellbender") {
    if (!"cellbender_nCount_RNA" %in% names(md) || anyNA(md$cellbender_nCount_RNA) ||
        !identical(lib, as.numeric(md$cellbender_nCount_RNA))) fail("Corrected library totals disagree with frozen CellBender audit: ", sample)
  }
  md$representation_full_library_UMIs <- lib
  if (contract$count_source == "paired_cellbender") md$corrected_nCount_RNA <- lib
  list(counts = x, metadata = md, library = lib, source_path = source_path, source_sha256 = actual_sha256)
}

log_full_library <- function(x, libraries, target) {
  y <- x %*% Matrix::Diagonal(x = target / libraries)
  y <- methods::as(y, "dgCMatrix"); dimnames(y) <- dimnames(x)
  y@x <- log1p(y@x); y
}

read_labels <- function(cfg, contract) {
  if (is.null(cfg$labels_tsv)) fail("External MapMyCells label TSV required for annotation")
  z <- read_tsv(cfg$labels_tsv)
  if (!"nucleus_id" %in% names(z)) fail("Label schema lacks nucleus_id")
  check_unique(z$nucleus_id, "external label nuclei")
  if (contract$count_source == "original_full") {
    if (!setequal(z$nucleus_id, contract$metadata$nucleus_id)) fail("External labels must cover EXACT frozen primary roster")
  } else if (any(!contract$metadata$nucleus_id %in% z$nucleus_id)) fail("Fixed original external labels do not cover paired roster")
  z <- z[match(contract$metadata$nucleus_id, z$nucleus_id), , drop = FALSE]
  allowed <- grepl("^(reference_|REF_|mapmycells_)", names(z)) |
    names(z) %in% c("analysis_label", "CA1_supported", "annotation_status")
  z <- z[, c("nucleus_id", names(z)[allowed]), drop = FALSE]
  rownames(z) <- z$nucleus_id
  z
}

add_labels <- function(obj, labels) {
  fields <- setdiff(names(labels), "nucleus_id")
  obj <- SeuratObject::AddMetaData(obj, labels[colnames(obj), fields, drop = FALSE])
  obj$annotation_provenance <- "source_verified_existing_MapMyCells_reused"
  obj
}

runtime_record <- function(cfg, out, stage) {
  pkgs <- c("Seurat", "SeuratObject", "Matrix", "hdf5r", "harmony", "sctransform", "glmGamPoi", "uwot", "igraph", "DESeq2", "edgeR")
  write_tsv(data.frame(package = c("R", pkgs), version = c(as.character(getRversion()),
    vapply(pkgs, function(p) if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p)) else "ABSENT", character(1)))),
    file.path(out, "PACKAGE_VERSIONS.tsv"))
  writeLines(capture.output(utils::sessionInfo()), file.path(out, "SESSION_INFO.txt"))
  jsonlite::write_json(list(stage = stage, started_utc = utc(), config = cfg,
    status = "EXECUTION_IN_PROGRESS", scientific_acceptance = "NOT_ASSESSED"), file.path(out, "START.json"), pretty = TRUE, auto_unbox = TRUE)
}

memory_guard <- function(cfg) {
  if (is.null(cfg$memory_probe_python) || !nzchar(cfg$memory_probe_python)) fail("memory_probe_python is required for available-RAM guard")
  command <- "import psutil; print(psutil.virtual_memory().available/(1024**3))"
  answer <- system2(cfg$memory_probe_python, args = c("-c", shQuote(command)), stdout = TRUE, stderr = TRUE)
  available <- suppressWarnings(as.numeric(answer))
  if (length(available) != 1L || !is.finite(available)) fail("Available-memory probe failed; no scientific fitting started")
  minimum <- as.numeric(cfg$minimum_available_gib %or% 24)
  if (available < minimum) fail("Available RAM ", round(available, 2), " GiB below stage guard ", minimum, " GiB")
  data.frame(probed_utc = utc(), available_gib = available, minimum_available_gib = minimum,
    resource_budget_gib = cfg$resource_budget_gib, probe = "psutil.virtual_memory.available; stage-start guard")
}

plot_representation <- function(obj, out, prefix) {
  groups <- intersect(c("condition", "sample", "sequencing_partition", "reference_broad_prediction", "reference_subclass_name"), names(obj[[]]))
  for (reduction in intersect(c("umap_unintegrated", "umap_harmony"), SeuratObject::Reductions(obj))) {
    fields <- c(groups[!groups %in% "reference_subclass_name"], paste0(prefix, "_", sub("umap_", "", reduction), "_r", gsub("\\.", "p", as.character(CFG$primary_resolution %or% .4))))
    fields <- fields[fields %in% names(obj[[]])]
    for (field in fields) {
      p <- Seurat::DimPlot(obj, reduction = reduction, group.by = field, raster = TRUE, raster.dpi = c(1600, 1600), pt.size = .2)
      ggplot2::ggsave(file.path(out, paste0(reduction, "_", safe(field), ".pdf")), p, width = 10, height = 8, limitsize = FALSE)
    }
    p <- Seurat::DimPlot(obj, reduction = reduction, group.by = "condition", split.by = "sample", ncol = 4,
      raster = TRUE, raster.dpi = c(800, 800), pt.size = .2)
    ggplot2::ggsave(file.path(out, paste0(reduction, "_all12_mouse_facets.pdf")), p, width = 16, height = 11, limitsize = FALSE)
  }
}

representation <- function(contract, cfg, out, prefix = "log", ids = NULL, sct = FALSE) {
  if (sct && contract$count_source != "original_full") fail("SCT sensitivity is original_full only")
  features <- contract$features
  if (!is.null(ids)) {
    if (length(ids) < 100L) fail("Insufficient nuclei for subset representation")
    if (any(!ids %in% contract$metadata$nucleus_id)) fail("Unknown subset nuclei")
  }
  n_expected <- if (is.null(ids)) nrow(contract$metadata) else length(ids)
  one_dense_gib <- as.double(length(features)) * n_expected * 8 / 1024^3
  preliminary_peak <- 6 * one_dense_gib + 6
  write_tsv(data.frame(expected_nuclei = n_expected, nfeatures = length(features),
    one_dense_scale_gib = one_dense_gib, preliminary_peak_gib = preliminary_peak,
    resource_budget_gib = cfg$resource_budget_gib,
    caveat = "pre-fit conservative estimate; actual sparse count sizes refine it after loading"), file.path(out, "PRE_FIT_RESOURCE_PLAN.tsv"))
  if (preliminary_peak > cfg$resource_budget_gib) fail("Pre-fit memory plan exceeds target budget")
  parts <- list(); data_parts <- list(); residual_parts <- list(); modeled_features <- list(); meta <- list(); source_audit <- list()
  for (sample in contract$samples$sample) {
    message("Reading ", contract$count_source, " capture ", sample)
    cap <- read_capture(sample, contract, ids)
    if (is.null(cap)) next
    h <- cap$counts[features, , drop = FALSE]
    d <- log_full_library(h, cap$library, as.numeric(cfg$normalization_target %or% 10000))
    if (sct) {
      # Fit on the full capture; request only features with enough observed nuclei
      # for the declared SCT min_cells=5 model policy. The shared intersection is
      # determined after all captures, never filled with invented zero residuals.
      detection <- Matrix::rowSums(cap$counts[features, , drop = FALSE] > 0)
      requested <- features[detection >= 5L]
      write_tsv(data.frame(gene_id = features, detected_nuclei = as.integer(detection),
        requested_for_SCT_residual = features %in% requested,
        reason = ifelse(features %in% requested, "AT_LEAST_FIVE_DETECTED_NUCLEI", "FEWER_THAN_FIVE_DETECTED_NUCLEI")),
        file.path(out, paste0("SCT_REQUESTED_FEATURES_", sample, ".tsv")))
      temp <- Seurat::CreateSeuratObject(cap$counts, min.cells = 0, min.features = 0)
      temp <- Seurat::SCTransform(temp, vst.flavor = "v2", residual.features = requested, min_cells = 5L,
        variable.features.n = length(features), ncells = min(5000L, ncol(temp)),
        vars.to.regress = NULL, do.correct.umi = FALSE, conserve.memory = TRUE,
        return.only.var.genes = FALSE, seed.use = cfg$seed, verbose = TRUE)
      resid <- SeuratObject::GetAssayData(temp, assay = "SCT", layer = "scale.data")
      available <- intersect(features, rownames(resid))
      modeled_features[[sample]] <- available
      residual_parts[[sample]] <- resid[available, colnames(h), drop = FALSE]
      rm(temp, resid); gc()
    }
    parts[[sample]] <- h; data_parts[[sample]] <- d; meta[[sample]] <- cap$metadata
    source_audit[[sample]] <- data.frame(sample = sample, nuclei = ncol(h), full_genes = nrow(cap$counts),
      count_source = contract$count_source, source_full_UMIs = sum(cap$library), HVG_UMIs = sum(h), full_capture_sparse_gib = as.numeric(object.size(cap$counts)) / 1024^3)
    source_audit[[sample]]$source_path <- cap$source_path
    source_audit[[sample]]$source_sha256 <- cap$source_sha256
    rm(cap, h, d); gc()
  }
  x <- do.call(cbind, parts); normalized <- do.call(cbind, data_parts)
  md <- do.call(rbind, meta); rownames(md) <- md$nucleus_id
  md <- md[match(colnames(x), md$nucleus_id), , drop = FALSE]
  rm(parts, data_parts, meta); gc()
  if (is.null(ids) && !setequal(colnames(x), contract$metadata$nucleus_id)) fail("Incomplete global roster")
  dense_gib <- as.double(nrow(x)) * ncol(x) * 8 / 1024^3
  estimate <- 5 * as.numeric(object.size(x)) / 1024^3 + 4 * dense_gib + 2
  write_tsv(data.frame(n_nuclei = ncol(x), n_features = nrow(x), sparse_HVG_counts_gib = as.numeric(object.size(x)) / 1024^3,
    one_dense_scale_matrix_gib = dense_gib, planning_peak_gib = estimate, declared_budget_gib = cfg$resource_budget_gib,
    caveat = "conservative planning estimate, not measured RSS; capture loading/SCT peaks logged separately"), file.path(out, "RESOURCE_ESTIMATE.tsv"))
  if (estimate > cfg$resource_budget_gib) fail("Planned representation exceeds memory budget; choose a larger target, preserve all nuclei")
  obj <- Seurat::CreateSeuratObject(counts = x, meta.data = md, min.cells = 0, min.features = 0)
  obj$original_nCount_RNA <- md$nCount_RNA
  obj$HVG_nCount_RNA <- as.numeric(Matrix::colSums(x))
  obj <- SeuratObject::SetAssayData(obj, assay = "RNA", layer = "data", new.data = normalized)
  SeuratObject::VariableFeatures(obj) <- features
  rm(x, normalized); gc()
  if (sct) {
    common_modeled <- Reduce(intersect, modeled_features)
    coverage <- data.frame(gene_id = features, n_captures_with_modeled_residual = vapply(features, function(id) sum(vapply(modeled_features, function(z) id %in% z, logical(1))), integer(1)))
    coverage$included_in_SCT_representation <- coverage$gene_id %in% common_modeled
    coverage$reason <- ifelse(coverage$included_in_SCT_representation, "MODELED_IN_ALL12", "NO_FITTED_RESIDUAL_IN_ONE_OR_MORE_CAPTURES; low_detection_or_model_eligibility")
    write_tsv(coverage, file.path(out, "SCT_MODELED_FEATURE_COVERAGE.tsv"))
    for (sample in names(modeled_features)) write_tsv(data.frame(gene_id = modeled_features[[sample]]), file.path(out, paste0("SCT_MODELED_FEATURES_", sample, ".tsv")))
    features <- features[features %in% common_modeled]
    if (length(features) <= as.integer(cfg$npcs %or% 30L)) fail("Insufficient shared modeled SCT features; no zero residuals fabricated")
    residuals <- do.call(cbind, lapply(residual_parts, function(x) x[features, , drop = FALSE])); residuals <- residuals[, colnames(obj), drop = FALSE]
    obj <- SeuratObject::SetAssayData(obj, assay = "RNA", layer = "scale.data", new.data = residuals)
    SeuratObject::VariableFeatures(obj) <- features
    rm(residual_parts, residuals); gc()
  } else obj <- Seurat::ScaleData(obj, features = features, vars.to.regress = NULL,
    do.center = TRUE, do.scale = TRUE, scale.max = as.numeric(cfg$scale_max %or% 10), verbose = TRUE)
  npcs <- as.integer(cfg$npcs %or% 30L)
  obj <- Seurat::RunPCA(obj, features = features, npcs = npcs, seed.use = cfg$seed, weight.by.var = TRUE, verbose = TRUE)
  if (ncol(SeuratObject::Embeddings(obj, "pca")) != npcs) fail("Requested PC count not supported")
  write_tsv(data.frame(gene_id = rownames(SeuratObject::Loadings(obj[["pca"]]))),
    file.path(out, "PCA_ACTUAL_FEATURES.tsv"))
  set.seed(cfg$seed)
  obj <- harmony::RunHarmony(obj, group.by.vars = "sample", reduction.use = "pca", dims.use = seq_len(npcs),
    reduction.save = "harmony", project.dim = FALSE, theta = as.numeric(cfg$harmony_theta %or% 2),
    sigma = as.numeric(cfg$harmony_sigma %or% .1), lambda = as.numeric(cfg$harmony_lambda %or% 1),
    nclust = as.integer(cfg$harmony_nclust %or% 100L), max_iter = as.integer(cfg$harmony_max_iter %or% 10L),
    ncores = as.integer(cfg$harmony_ncores %or% 4L), verbose = TRUE)
  # PCA/Harmony are sufficient for every graph/UMAP step. Drop the now-unused
  # dense scaled/residual matrix before graph allocation and publish a lean
  # checkpoint so a later graph/API failure leaves inspectable fitted evidence.
  obj <- drop_scaled_layer(obj)
  for (red in c("pca", "harmony")) write_tsv(data.frame(nucleus_id = colnames(obj),
    SeuratObject::Embeddings(obj, red), check.names = FALSE), file.path(out, paste0(red, ".tsv")))
  write_tsv(do.call(rbind, source_audit), file.path(out, "CAPTURE_COUNT_CONSERVATION.tsv"))
  save_new(obj, file.path(out, "PRE_GRAPH_HVG_ONLY.rds"))
  write_tsv(data.frame(implementation = "sparse_lower_edges_to_igraph_then_native_Seurat_RunLeiden_GroupSingletons",
    objective = "modularity", seed = cfg$seed, n_iter = 10L, self_edges = FALSE,
    singleton_grouping = TRUE, semantic_change = "none; bypass dense lower.tri adjacency check"),
    file.path(out, "CLUSTERING_IMPLEMENTATION.tsv"))
  resolutions <- as.numeric(cfg$resolutions %or% c(.2, .4, .6))
  for (branch in c("unintegrated", "harmony")) {
    red <- if (branch == "unintegrated") "pca" else "harmony"
    message("Building sparse ", branch, " kNN/SNN graph")
    obj <- Seurat::FindNeighbors(obj, reduction = red, dims = seq_len(npcs), k.param = as.integer(cfg$k_neighbors %or% 20L),
      nn.method = "rann", nn.eps = 0, l2.norm = FALSE, prune.SNN = 1/15,
      graph.name = c(paste0("nn_", branch), paste0("snn_", branch)), verbose = TRUE)
    snn <- obj[[paste0("snn_", branch)]]
    graph <- snn_to_igraph_sparse(snn)
    write_tsv(data.frame(branch = branch, n_vertices = igraph::vcount(graph),
      n_undirected_edges = igraph::ecount(graph), snn_stored_values = length(snn@x),
      snn_gib = as.numeric(object.size(snn))/1024^3, graph_gib = as.numeric(object.size(graph))/1024^3),
      file.path(out, paste0("GRAPH_RESOURCE_", branch, ".tsv")))
    save_new(snn, file.path(out, paste0("SNN_", branch, ".rds")))
    for (r in resolutions) {
      key <- paste0(prefix, "_", branch, "_r", gsub("\\.", "p", as.character(r)))
      message("Native Seurat modularity Leiden on sparse ", branch, " graph, resolution ", r)
      ids <- sparse_leiden_membership(graph, snn, resolution = r, random_seed = cfg$seed, n_iter = 10L)
      obj[[key]] <- ids
      SeuratObject::Idents(obj) <- ids
      obj[["seurat_clusters"]] <- SeuratObject::Idents(obj)
    }
    rm(graph, snn); gc()
    obj <- Seurat::RunUMAP(obj, reduction = red, dims = seq_len(npcs), reduction.name = paste0("umap_", branch),
      reduction.key = if (branch == "harmony") "UMAPH_" else "UMAPU_", umap.method = "uwot",
      n.neighbors = as.integer(cfg$umap_neighbors %or% 30L), min.dist = as.numeric(cfg$umap_min_dist %or% .3),
      metric = "cosine", seed.use = cfg$seed, verbose = TRUE)
  }
  if (!is.null(cfg$labels_tsv)) obj <- add_labels(obj, read_labels(cfg, contract))
  obj$representation_method <- if (sct) "per_capture_SCT_v2_shared_modeled_feature_intersection_sensitivity" else "full_library_log1p_fixed_common_features"
  obj$count_source <- contract$count_source
  obj$representation_count_scope <- "HVG_ONLY; inference reopens all source genes"
  if (!prefix %in% c("log", "sct")) {
    for (field in grep(paste0("^", prefix, "_(unintegrated|harmony)_r"), names(obj[[]]), value = TRUE)) {
      alias <- sub(paste0("^", prefix, "_"), "log_", field)
      obj[[alias]] <- obj[[field, drop = TRUE]]
    }
  }
  for (red in setdiff(SeuratObject::Reductions(obj), c("pca", "harmony"))) write_tsv(data.frame(nucleus_id = colnames(obj), SeuratObject::Embeddings(obj, red), check.names = FALSE), file.path(out, paste0(red, ".tsv")))
  write_tsv(data.frame(nucleus_id = colnames(obj), obj[[]], check.names = FALSE)[, !duplicated(c("nucleus_id", names(obj[[]]))), drop = FALSE], file.path(out, "nucleus_metadata.tsv"))
  write_tsv(data.frame(gene_id = features), file.path(out, "fixed_features.tsv"))
  plot_representation(obj, out, prefix)
  # PCA already exists; remove dense residual/scaled matrix from durable checkpoint.
  obj <- drop_scaled_layer(obj)
  save_new(obj, file.path(out, "representation_HVG_only.rds"))
  invisible(obj)
}

annotate_subsets <- function(contract, cfg, out) {
  labels <- read_labels(cfg, contract)
  subsets <- list(EX = if ("REF_EX_HS" %in% names(labels)) truth(labels$REF_EX_HS) else rep(FALSE, nrow(labels)),
    IN = if ("REF_IN_HS" %in% names(labels)) truth(labels$REF_IN_HS) else rep(FALSE, nrow(labels)))
  if ("CA1_supported" %in% names(labels)) subsets$CA1_SUPPORTED <- truth(labels$CA1_supported)
  else if ("REF_CA1_PROS_HS" %in% names(labels)) subsets$CA1_PROS_CANDIDATE <- truth(labels$REF_CA1_PROS_HS)
  if (any(subsets$EX & subsets$IN)) fail("External EX and IN high-support memberships overlap")
  if (length(subsets) < 3L) fail("No supported CA1 or explicitly candidate CA1/ProS membership supplied")
  metadata <- contract$metadata
  for (key in names(subsets)) {
    ids <- labels$nucleus_id[subsets[[key]]]
    folder <- new_dir(file.path(out, key))
    write_tsv(data.frame(nucleus_id = ids), file.path(folder, "subset_membership.tsv"))
    coverage <- table(factor(metadata$sample[match(ids, metadata$nucleus_id)], levels = names(roster)))
    write_tsv(data.frame(sample = names(roster), condition = unname(roster), n_nuclei = as.integer(coverage)), file.path(folder, "ALL12_COVERAGE.tsv"))
    if (length(ids) < 100L) {
      writeLines("NOT_RUN: fewer than 100 supported nuclei; no label or membership forced.", file.path(folder, "NOT_RUN.txt")); next
    }
    sub_contract <- contract
    if (!is.null(cfg$subset_features_tsv[[key]])) {
      features <- read_tsv(cfg$subset_features_tsv[[key]])$gene_id
      check_unique(features, paste0(key, " features"))
      if (any(!features %in% contract$genes$gene_id)) fail("Unknown subset feature IDs")
      sub_contract$features <- features
    }
    obj <- representation(sub_contract, cfg, folder, prefix = key, ids = ids)
    fields <- grep(paste0("^", key, "_"), names(obj[[]]), value = TRUE)
    for (field in fields) {
      metadata[[field]] <- NA_character_
      metadata[[field]][match(colnames(obj), metadata$nucleus_id)] <- as.character(obj[[field, drop = TRUE]])
    }
    rm(obj); gc()
  }
  write_tsv(metadata, file.path(out, "subset_cluster_metadata.tsv"))
}

# Marker exploration uses a recorded balanced subsample. It never changes the
# all-nucleus clustering roster, labels or animal-level inferential count source.
marker_sampling <- function(md, field, cfg) {
  if (!field %in% names(md)) fail("Marker grouping column unavailable: ", field)
  keep <- !is.na(md[[field]]) & nzchar(as.character(md[[field]]))
  z <- md[keep, , drop = FALSE]
  groups <- sort(unique(as.character(z[[field]])))
  cap <- min(as.integer(cfg$marker_max_cells_per_cluster %or% 600L),
    max(1L, floor(as.integer(cfg$marker_max_total_cells %or% 30000L) / length(groups))))
  set.seed(cfg$seed)
  chosen <- unlist(lapply(groups, function(group) {
    pool <- z[z[[field]] == group, , drop = FALSE]
    by_sample <- split(pool$nucleus_id, pool$sample)
    per_sample <- max(1L, floor(cap / length(by_sample)))
    unlist(lapply(by_sample, function(ids) {
      ids <- sort(ids)
      if (length(ids) > per_sample) ids[sample.int(length(ids), per_sample)] else ids
    }), use.names = FALSE)
  }), use.names = FALSE)
  audit <- z[, c("nucleus_id", "sample", "condition", field), drop = FALSE]
  audit$selected_for_marker_exploration <- audit$nucleus_id %in% chosen
  audit$sampling_policy <- paste0("seeded sample-balanced <=", cap, " nuclei/cluster; marker exploration only")
  list(ids = chosen, audit = audit)
}

full_gene_marker_object <- function(ids, md, contract, cfg) {
  counts <- list(); normalized <- list()
  for (sample in contract$samples$sample) {
    cap <- read_capture(sample, contract, ids)
    if (is.null(cap)) next
    counts[[sample]] <- cap$counts
    normalized[[sample]] <- log_full_library(cap$counts, cap$library, as.numeric(cfg$normalization_target %or% 10000))
    rm(cap); gc()
  }
  x <- do.call(cbind, counts); d <- do.call(cbind, normalized)
  rm(counts, normalized); gc()
  estimate <- 5 * as.numeric(object.size(x)) / 1024^3 + 2
  if (estimate > cfg$resource_budget_gib) fail("Marker sparse-object planning estimate exceeds declared budget")
  obj <- Seurat::CreateSeuratObject(x, min.cells = 0, min.features = 0, meta.data = md[colnames(x), , drop = FALSE])
  obj <- SeuratObject::SetAssayData(obj, assay = "RNA", layer = "data", new.data = d)
  rm(x, d); gc(); obj
}

combine_cluster_metadata <- function(contract, cfg) {
  md <- contract$metadata
  # Registries come from full-original source metadata, before selecting paired
  # nuclei. A real original cluster can therefore remain a declared endpoint
  # even when its paired membership is empty; no label or nucleus is invented.
  level_registry <- list(); registry_sources <- list()
  register <- function(field, values, path) {
    text <- as.character(values)
    text <- text[!is.na(text) & nzchar(text)]
    level_registry[[field]] <<- sort(unique(text))
    registry_sources[[field]] <<- path
  }
  source_root <- cfg$cluster_metadata_root %or% cfg$output_root
  global_path <- file.path(source_root, "log_global", "nucleus_metadata.tsv")
  global <- read_tsv(global_path)
  fields <- grep("^log_.*_r", names(global), value = TRUE)
  if (contract$count_source == "original_full" && !setequal(global$nucleus_id, md$nucleus_id)) fail("Global cluster metadata roster changed")
  if (any(!md$nucleus_id %in% global$nucleus_id)) fail("Fixed original global clusters do not cover paired roster")
  for (field in fields) {
    register(field, global[[field]], global_path)
    md[[field]] <- global[[field]][match(md$nucleus_id, global$nucleus_id)]
  }
  # Raw external label levels are preserved as source-defined names. The usual
  # high-support membership filter still applies, including to zero endpoints.
  for (field in grep("^reference_.*(prediction|name|label)$", names(global), value = TRUE)) register(field, global[[field]], global_path)
  sub_path <- file.path(source_root, "annotate_subsets", "subset_cluster_metadata.tsv")
  if (file.exists(sub_path)) {
    sub <- read_tsv(sub_path)
    if (contract$count_source == "original_full" && !setequal(sub$nucleus_id, md$nucleus_id)) fail("Subset cluster metadata roster changed")
    if (any(!md$nucleus_id %in% sub$nucleus_id)) fail("Fixed original subset clusters do not cover paired roster")
    fields <- grep("^(EX|IN|CA1_SUPPORTED|CA1_PROS_CANDIDATE)_.*_r", names(sub), value = TRUE)
    for (field in fields) {
      register(field, sub[[field]], sub_path)
      md[[field]] <- sub[[field]][match(md$nucleus_id, sub$nucleus_id)]
    }
  }
  labels <- read_labels(cfg, contract)
  for (field in setdiff(names(labels), "nucleus_id")) md[[field]] <- labels[[field]][match(md$nucleus_id, labels$nucleus_id)]
  if (!is.null(cfg$scanpy_metadata_paths)) {
    for (population in names(cfg$scanpy_metadata_paths)) {
      tab <- read_tsv(cfg$scanpy_metadata_paths[[population]])
      if (!all(c("nucleus_id", "cluster_harmony", "cluster_unintegrated") %in% names(tab))) fail("Scanpy metadata cluster schema: ", population)
      check_unique(tab$nucleus_id, paste0("Scanpy ", population))
      if (population == "allcells" && any(!md$nucleus_id %in% tab$nucleus_id)) fail("Scanpy allcells roster incomplete")
      for (branch in c("harmony", "unintegrated")) {
        field <- paste0("SCANPY_", toupper(safe(population)), "_", branch)
        register(field, tab[[paste0("cluster_", branch)]], cfg$scanpy_metadata_paths[[population]])
        md[[field]] <- as.character(tab[[paste0("cluster_", branch)]][match(md$nucleus_id, tab$nucleus_id)])
      }
    }
  }
  rownames(md) <- md$nucleus_id
  attr(md, "population_level_registry") <- level_registry
  attr(md, "population_level_registry_sources") <- registry_sources
  md
}

marker_tables_and_pdfs <- function(contract, cfg, out) {
  md <- combine_cluster_metadata(contract, cfg)
  resolution <- gsub("\\.", "p", as.character(cfg$primary_resolution %or% .4))
  fields <- cfg$marker_cluster_fields %or% grep(paste0("_harmony_r", resolution, "$"), names(md), value = TRUE)
  for (field in fields) {
    folder <- new_dir(file.path(out, safe(field)))
    sampled <- marker_sampling(md, field, cfg)
    write_tsv(sampled$audit, file.path(folder, "MARKER_SUBSAMPLING_MANIFEST.tsv"))
    if (length(unique(md[sampled$ids, field])) < 2L) {
      writeLines("NOT_RUN: fewer than two represented clusters.", file.path(folder, "NOT_RUN.txt")); next
    }
    obj <- full_gene_marker_object(sampled$ids, md, contract, cfg)
    SeuratObject::Idents(obj) <- as.character(obj[[field, drop = TRUE]])
    markers <- Seurat::FindAllMarkers(obj, assay = "RNA", only.pos = TRUE, test.use = "wilcox", slot = "data",
      min.pct = as.numeric(cfg$marker_min_pct %or% .1), logfc.threshold = as.numeric(cfg$marker_logfc_threshold %or% .25),
      return.thresh = 1, densify = FALSE, random.seed = cfg$seed, verbose = TRUE)
    if (!nrow(markers)) { writeLines("No exploratory marker rows passed declared gates.", file.path(folder, "EMPTY.txt")); rm(obj); gc(); next }
    markers$gene_id <- as.character(markers$gene)
    markers$gene_symbol <- contract$genes$display_symbol[match(markers$gene_id, contract$genes$gene_id)]
    markers$marker_inference_scope <- "sampled_nucleus_rank_exploration; not WT/cKO animal inference"
    markers$specificity_pct_difference <- markers$pct.1 - markers$pct.2
    markers <- markers[order(markers$cluster, -markers$specificity_pct_difference, -markers$avg_log2FC), , drop = FALSE]
    write_tsv(markers, file.path(folder, "ALL_EXPLORATORY_MARKERS.tsv"))
    utils::write.csv(markers, file.path(folder, "ALL_EXPLORATORY_MARKERS.csv"), row.names = FALSE)
    top <- do.call(rbind, lapply(split(markers, markers$cluster), function(z) head(z, 7L)))
    write_tsv(top, file.path(folder, "TOP7_SPECIFIC_MARKERS.tsv"))
    for (group in unique(top$cluster)) {
      ids <- unique(top$gene_id[top$cluster == group])
      p <- Seurat::DotPlot(obj, assay = "RNA", features = ids, group.by = field, scale = TRUE) +
        ggplot2::ggtitle(paste(field, "cluster", group, "marker expression; exploratory subsample")) +
        ggplot2::scale_x_discrete(labels = setNames(contract$genes$display_symbol[match(ids, contract$genes$gene_id)], ids)) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
      write_tsv(p$data, file.path(folder, paste0("cluster_", safe(group), "_DOTPLOT_NUMERICAL_VALUES.tsv")))
      ggplot2::ggsave(file.path(folder, paste0("cluster_", safe(group), "_TOP7_DOTPLOT.pdf")), p, width = 10, height = max(5, min(16, length(unique(obj[[field, drop = TRUE]])) * .25 + 3)))
    }
    if (!is.null(cfg$marker_tsv) && file.exists(cfg$marker_tsv)) {
      panel <- read_tsv(cfg$marker_tsv)
      if (!all(c("set", "gene_symbol") %in% names(panel))) fail("marker_tsv requires set and gene_symbol")
      for (set in unique(panel$set)) {
        ids <- contract$genes$gene_id[contract$genes$gene_symbol %in% panel$gene_symbol[panel$set == set]]
        if (!length(ids)) next
        p <- Seurat::DotPlot(obj, features = ids, assay = "RNA", group.by = field) + ggplot2::ggtitle(paste(set, "curated program; not a label gate")) +
          ggplot2::scale_x_discrete(labels = setNames(contract$genes$display_symbol[match(ids, contract$genes$gene_id)], ids)) +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1))
        write_tsv(p$data, file.path(folder, paste0("curated_", safe(set), "_NUMERICAL_VALUES.tsv")))
        ggplot2::ggsave(file.path(folder, paste0("curated_", safe(set), ".pdf")), p, width = max(10, min(24, length(ids) * .25)), height = 10)
      }
    }
    rm(obj); gc()
  }
}

population_manifest <- function(md, cfg) {
  membership <- list(); audit <- list(); index <- 0L
  add <- function(family, label, ids) {
    index <<- index + 1L
    key <- paste0(sprintf("P%04d_", index), safe(family))
    membership[[key]] <<- ids
    audit[[key]] <<- data.frame(population_key = key, family = family, label = label,
      n_nuclei = length(ids), source_label_branch = if (startsWith(family, "ATLAS_") || !is.null(cfg$cluster_metadata_root)) "original_full" else (cfg$count_source %or% "original_full"),
      label_interpretation = if (grepl("CA1_PROS", family)) "CA1/ProS_candidate_not_anatomical_CA1" else "declared_analysis_identity")
  }
  add("ATLAS_TOTAL_ALL", "ALL", md$nucleus_id)
  flag_map <- c(EX = "REF_EX_HS", IN = "REF_IN_HS", CA1_PROS = "REF_CA1_PROS_HS", CA1_SUPPORTED = "CA1_supported")
  for (unit in names(flag_map)) {
    flag <- flag_map[[unit]]
    if (flag %in% names(md)) add(paste0("ATLAS_TOTAL_", unit), unit, md$nucleus_id[truth(md[[flag]])])
  }
  resolution <- gsub("\\.", "p", as.character(cfg$primary_resolution %or% .4))
  fields <- cfg$pseudobulk_label_fields %or% c(intersect(c("reference_broad_prediction", "reference_subclass_name"), names(md)),
    grep(paste0("_(unintegrated|harmony)_r", resolution, "$"), names(md), value = TRUE), grep("^SCANPY_", names(md), value = TRUE))
  for (field in fields) {
    if (!field %in% names(md)) fail("Pseudobulk grouping unavailable: ", field)
    supported <- rep(TRUE, nrow(md))
    if (field == "reference_subclass_name" && "reference_subclass_high_support" %in% names(md)) supported <- truth(md$reference_subclass_high_support)
    if (field == "reference_broad_prediction" && "reference_broad_high_support" %in% names(md)) supported <- truth(md$reference_broad_high_support)
    valid <- supported & !is.na(md[[field]]) & nzchar(as.character(md[[field]]))
    family <- if (startsWith(field, "SCANPY_")) field else if (startsWith(field, "reference_")) paste0("ATLAS_", field) else paste0("SEURAT_", field)
    registry <- attr(md, "population_level_registry")[[field]]
    levels <- if (!is.null(registry)) registry else sort(unique(as.character(md[[field]][valid])))
    if (any(!as.character(md[[field]][valid]) %in% levels)) fail("Selected membership has a label outside its original source registry: ", field)
    for (level in levels) add(family, level, md$nucleus_id[valid & as.character(md[[field]]) == level])
  }
  if (!length(membership)) fail("No declared pseudobulk populations")
  list(memberships = membership, audit = do.call(rbind, audit))
}

# Sparse source counts are streamed once per capture. Only gene-by-animal
# aggregates are dense; no all-gene all-nucleus count object is assembled.
streamed_pseudobulk <- function(contract, cfg, out, script_dir) {
  source(file.path(script_dir, "pseudobulk_core.R"), local = TRUE)
  md <- combine_cluster_metadata(contract, cfg)
  registry <- attr(md, "population_level_registry")
  registry_sources <- attr(md, "population_level_registry_sources")
  registry_table <- do.call(rbind, lapply(names(registry), function(field) data.frame(
    source_column = rep(field, length(registry[[field]])),
    original_source_level = registry[[field]],
    source_metadata_path = rep(registry_sources[[field]], length(registry[[field]])))))
  if (!is.null(registry_table)) write_tsv(registry_table, file.path(out, "SOURCE_POPULATION_LEVEL_REGISTRY.tsv"))
  populations <- population_manifest(md, cfg)
  write_tsv(populations$audit, file.path(out, "POPULATION_INDEX.tsv"))
  members <- populations$memberships
  memberships <- do.call(rbind, lapply(names(members), function(key) data.frame(population_key = rep(key, length(members[[key]])), nucleus_id = members[[key]])))
  write_tsv(memberships, file.path(out, "POPULATION_MEMBERSHIPS.tsv"))
  genes <- contract$genes
  original_id_order <- genes$gene_id
  mouse <- contract$samples[, c("sample", "condition", "sequencing_partition"), drop = FALSE]
  rownames(mouse) <- mouse$sample
  expected <- satb2_pb_expected_partitions()[mouse$sample]
  if (any(mouse$sequencing_partition != unname(expected))) fail("Partition provenance changed; review the copied statistical-core partition contract before fitting")
  aggregates <- lapply(members, function(x) matrix(0, nrow(genes), nrow(mouse), dimnames = list(original_id_order, mouse$sample)))
  n_nuclei <- matrix(0L, length(members), nrow(mouse), dimnames = list(names(members), mouse$sample))
  conservation <- list()
  for (sample in mouse$sample) {
    cap <- read_capture(sample, contract)
    cap_ids <- colnames(cap$counts)
    ii <- list(); jj <- list()
    for (p in seq_along(members)) {
      ids <- intersect(members[[p]], cap_ids)
      n_nuclei[p, sample] <- length(ids)
      ii[[p]] <- match(ids, cap_ids); jj[[p]] <- rep.int(p, length(ids))
    }
    assignment <- Matrix::sparseMatrix(i = unlist(ii), j = unlist(jj), x = 1,
      dims = c(length(cap_ids), length(members)))
    y <- as.matrix(cap$counts %*% assignment)
    if (any(!is.finite(y)) || any(y < 0) || any(y != floor(y)) || any(y > .Machine$integer.max)) fail("Invalid integer pseudobulk aggregate: ", sample)
    for (p in seq_along(members)) {
      aggregates[[p]][, sample] <- y[, p]
      expected_total <- sum(cap$library[cap_ids %in% members[[p]]])
      if (sum(y[, p]) != expected_total) fail("Population count conservation failed: ", sample, "/", names(members)[p])
    }
    conservation[[sample]] <- data.frame(sample = sample, count_source = contract$count_source,
      source_nuclei = length(cap_ids), source_total_UMIs = sum(cap$library), aggregate_memberships_nuclei = sum(n_nuclei[, sample]),
      aggregate_memberships_UMIs = sum(y), population_sums_conserved = TRUE,
      note = "memberships can overlap; summed population totals are not disjoint source totals")
    rm(cap, y, assignment); gc()
  }
  write_tsv(do.call(rbind, conservation), file.path(out, "STREAMED_COUNT_CONSERVATION.tsv"))
  if (identical(contract$count_source, "paired_cellbender")) {
    if (is.null(cfg$prefilter_universe_root)) fail("CellBender count sensitivity requires frozen original-count gene universes")
    count_kind <- "CellBender_corrected_integer_count_sensitivity_fixed_original_labels"
  } else count_kind <- paste0("original_CellRanger_integer_UMI_", contract$count_source)
  statuses <- list()
  for (p in seq_along(members)) {
    key <- names(members)[p]
    folder <- new_dir(file.path(out, key)); auditdir <- new_dir(file.path(folder, "AUDIT"))
    y <- aggregates[[p]]; storage.mode(y) <- "integer"
    coverage <- mouse; coverage$n_nuclei <- n_nuclei[p, mouse$sample]; coverage$total_UMIs <- colSums(y)
    rownames(coverage) <- coverage$sample
    write_tsv(coverage, file.path(folder, "ALL12_POPULATION_COVERAGE.tsv"))
    write_tsv(data.frame(gene_id = genes$gene_id, gene_symbol = genes$display_symbol, y, check.names = FALSE), file.path(folder, "SOURCE_INTEGER_COUNTS_ALL_GENES_ALL12.tsv"))
    selection <- satb2_pb_select(coverage, "PRIMARY")
    prefilter <- satb2_pb_prefilter(y, selection)
    if (!is.null(cfg$prefilter_universe_root)) {
      source_index <- read_tsv(file.path(cfg$prefilter_universe_root, "POPULATION_INDEX.tsv"))
      here <- populations$audit[populations$audit$population_key == key, , drop = FALSE]
      same <- which(source_index$family == here$family & source_index$label == here$label)
      if (length(same) != 1L) fail("Fixed-universe population does not uniquely match original: ", key)
      old <- read_tsv(file.path(cfg$prefilter_universe_root, source_index$population_key[same], "AUDIT", "FROZEN_PRIMARY_GENE_PREFILTER_AUDIT.tsv"))
      universe <- old$gene_id[truth(old$passes_primary_prefilter)]
      prefilter$passes_primary_prefilter <- prefilter$gene_id %in% universe
      prefilter$universe_policy <- "frozen_from_original_count_branch"
    } else {
      universe <- prefilter$gene_id[prefilter$passes_primary_prefilter]
      prefilter$universe_policy <- "declared_primary_original_count_prefilter"
    }
    write_tsv(prefilter, file.path(auditdir, "FROZEN_PRIMARY_GENE_PREFILTER_AUDIT.tsv"))
    fits <- list(list(method = "DESeq2", endpoint = "PRIMARY"), list(method = "edgeR", endpoint = "PRIMARY"))
    if (truth(cfg$run_partition_sensitivities %or% (contract$count_source == "original_full"))) fits <- c(fits,
      list(list(method = "DESeq2", endpoint = "KNOWN_PARTITION_ADJUSTED"), list(method = "DESeq2", endpoint = "MAIN_PARTITION_ONLY")))
    loo_families <- cfg$loo_families %or% c("ATLAS_TOTAL_ALL", "ATLAS_TOTAL_EX", "ATLAS_TOTAL_IN", "ATLAS_TOTAL_CA1_PROS")
    aliases <- c(TOTAL_ALL = "ATLAS_TOTAL_ALL", TOTAL_REF_EX_HS = "ATLAS_TOTAL_EX", TOTAL_REF_IN_HS = "ATLAS_TOTAL_IN", TOTAL_REF_CA1_PROS_HS = "ATLAS_TOTAL_CA1_PROS")
    remap <- loo_families %in% names(aliases); loo_families[remap] <- unname(aliases[loo_families[remap]])
    family <- populations$audit$family[match(key, populations$audit$population_key)]
    do_loo <- truth(cfg$run_loo %or% (contract$count_source == "original_full")) && family %in% loo_families
    if (do_loo) fits <- c(fits, lapply(mouse$sample, function(sample) list(method = "DESeq2", endpoint = "LOO", omit = sample)))
    primary_results <- NULL; loo_results <- list()
    for (fit in fits) {
      label <- if (fit$endpoint == "LOO") paste0("LOO_omit_", fit$omit) else fit$endpoint
      ans <- satb2_pb_run_endpoint(y, coverage, genes, universe,
        outdir = file.path(folder, fit$method, label), method = fit$method, endpoint = fit$endpoint,
        omit = fit$omit, independent_filtering = fit$endpoint != "LOO", count_kind = count_kind,
        auditdir = file.path(auditdir, fit$method, label))
      ans$audit$population_key <- key; ans$audit$count_source <- contract$count_source
      statuses[[paste(key, fit$method, label, sep = "__")]] <- ans$audit
      if (fit$method == "DESeq2" && fit$endpoint == "PRIMARY") primary_results <- ans$results
      if (fit$endpoint == "LOO" && !is.null(ans$results)) loo_results[[fit$omit]] <- ans$results
    }
    if (!is.null(primary_results)) {
      z <- primary_results
      pdf(file.path(folder, "DESeq2_PRIMARY_MA.pdf"), width = 8, height = 6)
      ok <- is.finite(z$baseMean) & z$baseMean > 0 & is.finite(z$log2FoldChange)
      significant <- !is.na(z$padj) & z$padj < .05
      if (any(ok)) {
        plot(log10(z$baseMean[ok]), z$log2FoldChange[ok], pch = 16, cex = .25,
          col = ifelse(significant[ok], "#be3c39", "#8d939b"),
          xlab = "log10 DESeq2 baseMean", ylab = "cKO / WT log2 fold change",
          main = paste(key, contract$count_source))
        abline(h = 0, col = "#444444", lty = 2)
      } else { plot.new(); title("No finite fitted effects; complete table retained") }
      dev.off()
      pdf(file.path(folder, "DESeq2_PRIMARY_VOLCANO.pdf"), width = 8, height = 6)
      ok <- is.finite(z$log2FoldChange) & is.finite(z$padj) & z$padj >= 0
      yy <- -log10(pmax(z$padj, .Machine$double.xmin))
      if (any(ok)) {
        plot(z$log2FoldChange[ok], yy[ok], pch = 16, cex = .25,
          col = ifelse(significant[ok], "#be3c39", "#8d939b"),
          xlab = "cKO / WT log2 fold change", ylab = "-log10 BH adjusted P value",
          main = paste(key, contract$count_source))
        abline(h = -log10(.05), col = "#444444", lty = 2)
        label_rows <- head(which(ok)[order(z$padj[ok], -abs(z$log2FoldChange[ok]))], 10L)
        text(z$log2FoldChange[label_rows], yy[label_rows],
          labels = genes$display_symbol[match(z$gene_id[label_rows], genes$gene_id)], pos = 4, cex = .5)
      } else { plot.new(); title("No finite BH adjusted P values; complete table retained") }
      dev.off()
    }
    if (!is.null(primary_results) && length(loo_results)) {
      z <- primary_results
      lfc <- vapply(loo_results, function(x) x$log2FoldChange[match(z$gene_id, x$gene_id)], numeric(nrow(z)))
      if (is.null(dim(lfc))) lfc <- matrix(lfc, nrow = nrow(z))
      finite <- is.finite(lfc)
      influence <- data.frame(gene_id = z$gene_id, primary_log2FC = z$log2FoldChange, primary_padj = z$padj,
        n_LOO_estimable = rowSums(finite), n_LOO_same_direction = rowSums(finite & sign(lfc) == sign(z$log2FoldChange), na.rm = TRUE),
        min_LOO_log2FC = apply(lfc, 1L, function(x) if (any(is.finite(x))) min(x[is.finite(x)]) else NA_real_),
        max_LOO_log2FC = apply(lfc, 1L, function(x) if (any(is.finite(x))) max(x[is.finite(x)]) else NA_real_))
      write_tsv(influence, file.path(auditdir, "LOO_EFFECT_INFLUENCE_ONLY.tsv"))
    }
    rm(y); gc()
  }
  audit <- do.call(rbind, statuses)
  write_tsv(audit, file.path(out, "PSEUDOBULK_EXECUTION_AUDIT.tsv"))
  if (any(audit$status == "ERROR")) fail("Pseudobulk endpoint errors recorded; results require review")
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 2L) fail("Usage: Rscript --vanilla run_seurat.R CONFIG.json STAGE")
  CFG <<- jsonlite::fromJSON(args[1L], simplifyVector = TRUE)
  CFG$seed <- as.integer(CFG$seed %or% 20261013L)
  stage <- args[2L]
  if (!stage %in% c("log_global", "sct_global", "annotate_subsets", "markers", "pseudobulk")) fail("Unknown stage")
  needed <- c("Seurat", "SeuratObject", "Matrix", "hdf5r", "harmony", "igraph", "jsonlite", "ggplot2", "scattermore", "future", "digest")
  if (stage == "sct_global") needed <- c(needed, "sctransform", "glmGamPoi")
  if (stage == "pseudobulk") needed <- c(needed, "DESeq2", "edgeR", "statmod")
  missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) fail("Missing runtime packages; no installation performed: ", paste(missing, collapse = ", "))
  if (utils::packageVersion("Seurat") < "5.0.0") fail("Seurat v5 required")
  options(Seurat.object.assay.version = "v5", future.globals.maxSize = as.numeric(CFG$resource_budget_gib) * 1024^3)
  future::plan("sequential")
  set.seed(CFG$seed)
  memory <- memory_guard(CFG)
  contract <- load_contract(CFG)
  out <- new_dir(file.path(CFG$output_root, stage))
  file_arg <- grep("^--file=", commandArgs(), value = TRUE)
  script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[1L]), mustWork = TRUE))
  runtime_record(CFG, out, stage)
  write_tsv(memory, file.path(out, "STAGE_START_AVAILABLE_MEMORY.tsv"))
  write_tsv(data.frame(gene_id = contract$genes$gene_id, gene_symbol = contract$genes$gene_symbol), file.path(out, "GENE_ID_SYMBOL_MANIFEST.tsv"))
  input_files <- unique(c(CFG$sample_manifest_tsv, CFG$roster_tsv, CFG$gene_manifest_tsv, CFG$fixed_features_tsv, CFG$labels_tsv))
  input_files <- input_files[!is.na(input_files) & nzchar(input_files)]
  write_tsv(data.frame(path = input_files, sha256 = vapply(input_files, function(p) digest::digest(file = p, algo = "sha256"), character(1))), file.path(out, "CONFIG_INPUT_SHA256.tsv"))
  tryCatch({
    switch(stage,
      log_global = representation(contract, CFG, out),
      sct_global = representation(contract, CFG, out, prefix = "sct", sct = TRUE),
      annotate_subsets = annotate_subsets(contract, CFG, out),
      markers = marker_tables_and_pdfs(contract, CFG, out),
      pseudobulk = streamed_pseudobulk(contract, CFG, out, script_dir))
    jsonlite::write_json(list(status = "EXECUTED_REVIEW_REQUIRED", stage = stage, count_source = contract$count_source,
      completed_utc = utc(), source_nuclei = nrow(contract$metadata), scientific_acceptance = "NOT_ASSESSED"),
      file.path(out, "COMPLETE.json"), pretty = TRUE, auto_unbox = TRUE)
  }, error = function(e) {
    jsonlite::write_json(list(status = "FAILED_REVIEW_REQUIRED", stage = stage, message = conditionMessage(e), failed_utc = utc()),
      file.path(out, "FAILED.json"), pretty = TRUE, auto_unbox = TRUE)
    fail(conditionMessage(e))
  })
}

if (sys.nframe() == 0L) main()
