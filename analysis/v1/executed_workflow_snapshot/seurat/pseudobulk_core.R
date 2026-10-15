# Twelve independent mice: reusable original-UMI pseudobulk module.
# Source defines functions only. Call run_mouse_pseudobulk(...) explicitly.
# No installer, reference label forcing, marker gate, or automatic scientific acceptance.
# Dependencies: Matrix, DESeq2, edgeR, statmod, S4Vectors, SummarizedExperiment.
# Readiness tested: R 4.6.1; Matrix 1.7-6; DESeq2 1.52.0; edgeR 4.10.3; statmod 1.5.2.
# Unit tests use fabricated data; runtime versions are recorded on actual execution.
# Caller must provide ORIGINAL gene_id-by-nucleus integer UMIs. Integer checks cannot
# distinguish original counts from rounded corrected counts; provenance is the caller's duty.
# DESeq2 exports ALL_RESULTS.txt / DEG_FDR05.txt with seven standard fields.
# edgeR exports native fields, not a DESeq2-shaped analogy. Internal harmonized
# fields are comparison-only. Stable-ID/status evidence lives in the AUDIT tree.
# No directional-split or NOT_TESTED marker files are created.
# Caller must freeze reviewed count/metadata/population source hashes upstream.

satb2_pb_roster <- function() {
  c(wt_1122 = "WT", wt_1123 = "WT", wt_1178 = "WT", wt_1220 = "WT",
    wt_1226 = "WT", wt_1228 = "WT", cKO_1138 = "cKO", cKO_1139 = "cKO",
    cKO_1160 = "cKO", cKO_1165 = "cKO", cKO_933 = "cKO", cKO_991 = "cKO")
}

satb2_pb_verified_partitions <- function() c("23HJT3LT4_lanes6_7", "23HJTCLT4_lane8")
satb2_pb_main_partition <- function() "23HJT3LT4_lanes6_7"
satb2_pb_expected_partitions <- function() {
  p <- setNames(rep(satb2_pb_main_partition(), length(satb2_pb_roster())), names(satb2_pb_roster()))
  p[c("cKO_1138", "cKO_991")] <- "23HJTCLT4_lane8"
  p["cKO_933"] <- "unverified"
  p
}
satb2_pb_contrast <- function() c("condition", "cKO", "WT")
satb2_pb_result_columns <- function() c("gene_symbol", "baseMean", "log2FoldChange",
                                      "lfcSE", "stat", "pvalue", "padj")
satb2_pb_fail <- function(...) stop(paste0(..., collapse = ""), call. = FALSE)
satb2_pb_utc <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")

satb2_pb_write_table <- function(x, path) {
  if (file.exists(path)) satb2_pb_fail("Refusing to overwrite: ", path)
  tmp <- paste0(path, ".partial")
  if (file.exists(tmp)) satb2_pb_fail("Existing partial output: ", tmp)
  utils::write.table(x, tmp, sep = "\t", quote = FALSE, row.names = FALSE,
                     col.names = TRUE, na = "NA", fileEncoding = "UTF-8")
  if (!file.rename(tmp, path)) satb2_pb_fail("Cannot publish: ", path)
}

satb2_pb_write_text <- function(x, path) {
  if (file.exists(path)) satb2_pb_fail("Refusing to overwrite: ", path)
  tmp <- paste0(path, ".partial")
  if (file.exists(tmp)) satb2_pb_fail("Existing partial output: ", tmp)
  writeLines(x, tmp, useBytes = TRUE)
  if (!file.rename(tmp, path)) satb2_pb_fail("Cannot publish: ", path)
}

satb2_pb_validate_input <- function(counts, metadata, genes, populations) {
  if (!requireNamespace("Matrix", quietly = TRUE)) satb2_pb_fail("Matrix is required.")
  if (!inherits(counts, "sparseMatrix")) satb2_pb_fail("counts must be a sparse Matrix of original gene UMIs.")
  if (!isTRUE(methods::validObject(counts))) satb2_pb_fail("Invalid sparse count object.")
  counts <- methods::as(counts, "dgCMatrix")
  if (!nrow(counts) || !ncol(counts) || anyNA(counts@x) ||
      any(!is.finite(counts@x)) || any(counts@x < 0) || any(counts@x != floor(counts@x))) {
    satb2_pb_fail("Counts must be finite nonnegative whole-number original UMIs.")
  }
  for (z in list(rownames(counts), colnames(counts))) {
    if (is.null(z) || anyNA(z) || any(!nzchar(z)) || anyDuplicated(z) || any(grepl("[\t\r\n]", z)))
      satb2_pb_fail("Count axes must have unique, nonmissing IDs without tab/newline characters.")
  }
  if (any(Matrix::colSums(counts) <= 0)) satb2_pb_fail("Input QC nuclei must have positive original UMI libraries.")
  if (!is.data.frame(metadata) || !all(c("sample", "condition", "sequencing_partition") %in% names(metadata)))
    satb2_pb_fail("metadata requires nucleus ID rownames, sample, condition and explicit sequencing_partition.")
  metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)
  if (!"nucleus_id" %in% names(metadata)) metadata$nucleus_id <- rownames(metadata)
  if (anyNA(metadata$nucleus_id) || anyDuplicated(metadata$nucleus_id) ||
      !identical(rownames(metadata), as.character(metadata$nucleus_id)) ||
      !setequal(colnames(counts), as.character(metadata$nucleus_id)))
    satb2_pb_fail("Unique metadata rownames must equal nucleus_id and exactly match count columns.")
  metadata <- metadata[match(colnames(counts), metadata$nucleus_id), , drop = FALSE]
  metadata$sample <- as.character(metadata$sample)
  metadata$condition <- as.character(metadata$condition)
  metadata$sequencing_partition <- as.character(metadata$sequencing_partition)
  roster <- satb2_pb_roster(); samples <- names(roster)
  if (anyNA(metadata$sample) || !setequal(unique(metadata$sample), samples))
    satb2_pb_fail("The source must contain all twelve confirmed mice; no global sample exclusions.")
  if (anyNA(metadata$condition) || any(metadata$condition != unname(roster[metadata$sample])))
    satb2_pb_fail("Genotype disagrees with the confirmed six-WT/six-cKO roster.")
  if (anyNA(metadata$sequencing_partition) || any(!nzchar(trimws(metadata$sequencing_partition))))
    satb2_pb_fail("Unknown partitions must be explicit (e.g. unverified), never blank/NA.")
  metadata$sequencing_partition <- trimws(metadata$sequencing_partition)
  mouse <- data.frame(sample = samples, condition = unname(roster), stringsAsFactors = FALSE)
  mouse$sequencing_partition <- vapply(samples, function(s) {
    p <- unique(metadata$sequencing_partition[metadata$sample == s])
    if (length(p) != 1L) satb2_pb_fail("Inconsistent partition within mouse ", s)
    p
  }, character(1))
  if (any(mouse$sequencing_partition != unname(satb2_pb_expected_partitions()[samples])))
    satb2_pb_fail("Sequencing partitions disagree with the reviewed twelve-mouse assignment; cKO_933 must remain unverified.")
  if (!is.data.frame(genes) || !all(c("gene_id", "gene_symbol") %in% names(genes)))
    satb2_pb_fail("genes requires gene_id and gene_symbol.")
  genes <- as.data.frame(genes[, c("gene_id", "gene_symbol"), drop = FALSE], stringsAsFactors = FALSE)
  genes$gene_id <- as.character(genes$gene_id); genes$gene_symbol <- as.character(genes$gene_symbol)
  if (anyNA(genes$gene_id) || any(!nzchar(genes$gene_id)) || anyDuplicated(genes$gene_id) ||
      !setequal(genes$gene_id, rownames(counts))) satb2_pb_fail("Stable gene IDs must be unique and exactly match counts.")
  genes <- genes[match(rownames(counts), genes$gene_id), , drop = FALSE]; rownames(genes) <- NULL
  if (any(grepl("[\t\r\n]", genes$gene_symbol), na.rm = TRUE)) satb2_pb_fail("Gene symbols contain tab/newline characters.")
  genes$symbol_missing <- is.na(genes$gene_symbol) | !nzchar(trimws(genes$gene_symbol))
  genes$display_symbol <- ifelse(genes$symbol_missing, genes$gene_id, genes$gene_symbol)
  if (!is.list(populations) || !length(populations) || is.null(names(populations)) ||
      anyNA(names(populations)) || anyDuplicated(names(populations)) ||
      any(!grepl("^[A-Za-z0-9][A-Za-z0-9_-]*$", names(populations))))
    satb2_pb_fail("populations needs unique safe directory keys and lists of nucleus IDs.")
  for (name in names(populations)) {
    ids <- populations[[name]]
    if (!is.character(ids) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids) ||
        any(!ids %in% metadata$nucleus_id)) satb2_pb_fail("Invalid or duplicate population membership: ", name)
  }
  list(counts = counts, metadata = metadata, genes = genes, populations = populations, mouse = mouse)
}

satb2_pb_aggregate <- function(counts, metadata, ids, mouse) {
  j <- match(ids, metadata$nucleus_id); samples <- mouse$sample
  if (length(j)) {
    a <- Matrix::sparseMatrix(i = seq_along(j), j = match(metadata$sample[j], samples), x = 1,
                             dims = c(length(j), length(samples)))
    y <- as.matrix(counts[, j, drop = FALSE] %*% a)
    expected_gene <- as.numeric(Matrix::rowSums(counts[, j, drop = FALSE]))
  } else {
    y <- matrix(0, nrow(counts), length(samples)); expected_gene <- rep(0, nrow(counts))
  }
  dimnames(y) <- list(rownames(counts), samples)
  if (anyNA(y) || any(!is.finite(y)) || any(y < 0) || any(y != floor(y)) ||
      any(y > .Machine$integer.max) || !identical(as.numeric(rowSums(y)), expected_gene))
    satb2_pb_fail("Pseudobulk integer/count-conservation validation failed.")
  storage.mode(y) <- "integer"
  md <- mouse
  md$n_nuclei <- as.integer(table(factor(metadata$sample[j], levels = samples)))
  md$total_UMIs <- as.numeric(colSums(y))
  if (any((md$n_nuclei == 0) != (md$total_UMIs == 0))) satb2_pb_fail("Absent population and zero library do not agree.")
  rownames(md) <- samples
  list(counts = y, metadata = md, conserved_total_UMIs = sum(y))
}

satb2_pb_select <- function(md, endpoint, omit = NULL) {
  audit <- md
  eligible <- audit$n_nuclei >= 20L & audit$total_UMIs > 0
  reason <- ifelse(audit$n_nuclei == 0L, "ABSENT_POPULATION",
                   ifelse(audit$n_nuclei < 20L, "FEWER_THAN_20_NUCLEI", "INCLUDED"))
  verified <- audit$sequencing_partition %in% satb2_pb_verified_partitions()
  audit$partition_provenance_verified <- verified
  if (endpoint == "KNOWN_PARTITION_ADJUSTED") {
    reason[!verified] <- paste0(reason[!verified], ";UNVERIFIED_PARTITION")
    eligible <- eligible & verified
  } else if (endpoint == "MAIN_PARTITION_ONLY") {
    same <- audit$sequencing_partition == satb2_pb_main_partition()
    reason[!same] <- paste0(reason[!same], ";OUTSIDE_PRESPECIFIED_MAIN_PARTITION")
    eligible <- eligible & same
  } else if (!endpoint %in% c("PRIMARY", "LOO")) satb2_pb_fail("Unknown endpoint: ", endpoint)
  redundant <- FALSE
  if (endpoint == "LOO") {
    if (length(omit) != 1L || !omit %in% audit$sample) satb2_pb_fail("LOO omission must name one roster mouse.")
    redundant <- !eligible[match(omit, audit$sample)]
    eligible[audit$sample == omit] <- FALSE
    reason[audit$sample == omit] <- "PRESPECIFIED_LEAVE_ONE_MOUSE_OUT_OMISSION"
  }
  audit$included_in_model <- eligible
  audit$inclusion_or_exclusion_reason <- reason
  selected <- audit[eligible, , drop = FALSE]
  selected$condition <- factor(selected$condition, levels = c("WT", "cKO"))
  selected$sequencing_partition <- droplevels(factor(selected$sequencing_partition))
  nwt <- sum(selected$condition == "WT"); nko <- sum(selected$condition == "cKO")
  formula <- if (endpoint == "KNOWN_PARTITION_ADJUSTED") ~ sequencing_partition + condition else ~ condition
  skip <- ""; mm <- NULL; residual_df <- NA_integer_
  if (redundant) skip <- "OMITTED_MOUSE_ALREADY_INELIGIBLE_NO_MODEL_CHANGE"
  else if (min(nwt, nko) < 3L) skip <- "FEWER_THAN_THREE_ELIGIBLE_MICE_PER_GENOTYPE"
  else if (endpoint == "KNOWN_PARTITION_ADJUSTED" && nlevels(selected$sequencing_partition) < 2L)
    skip <- "FEWER_THAN_TWO_VERIFIED_PARTITIONS"
  else {
    mm <- stats::model.matrix(formula, selected)
    residual_df <- nrow(mm) - qr(mm)$rank
    if (qr(mm)$rank != ncol(mm)) skip <- "DESIGN_NOT_FULL_RANK"
    else if (residual_df < 2L) skip <- "FEWER_THAN_TWO_RESIDUAL_DEGREES_OF_FREEDOM"
  }
  list(audit = audit, metadata = selected, included = eligible, formula = formula,
       design = paste(deparse(formula), collapse = ""), model_matrix = mm,
       n_WT = nwt, n_cKO = nko, residual_df = residual_df, skip = skip,
       omitted_mouse = if (is.null(omit)) NA_character_ else omit)
}

satb2_pb_prefilter <- function(y, selection) {
  xx <- y[, selection$included, drop = FALSE]
  threshold <- min(selection$n_WT, selection$n_cKO)
  n <- rowSums(xx >= 10L)
  data.frame(gene_id = rownames(y), total_UMIs_all12 = rowSums(y),
             n_included_mice_with_ge10_UMIs = n, required_mice = threshold,
             passes_primary_prefilter = threshold > 0L & n >= threshold, stringsAsFactors = FALSE)
}

satb2_pb_prepare_endpoint <- function(y, selection, universe) {
  audit <- data.frame(gene_id = rownames(y), in_frozen_primary_universe = rownames(y) %in% universe,
                      total_UMIs_endpoint = rowSums(y[, selection$included, drop = FALSE]), stringsAsFactors = FALSE)
  audit$included_gene_in_endpoint <- audit$in_frozen_primary_universe & audit$total_UMIs_endpoint > 0
  audit$reason <- ifelse(!audit$in_frozen_primary_universe, "OUTSIDE_FROZEN_PRIMARY_UNIVERSE",
                         ifelse(audit$total_UMIs_endpoint == 0, "ALL_ZERO_AFTER_ENDPOINT_SELECTION", "INCLUDED"))
  x <- y[audit$included_gene_in_endpoint, selection$included, drop = FALSE]
  skip <- selection$skip
  if (!nzchar(skip) && !nrow(x)) skip <- "NO_NONZERO_GENES_IN_FROZEN_PRIMARY_UNIVERSE"
  if (!nzchar(skip) && any(colSums(x) <= 0)) skip <- "ZERO_INCLUDED_MOUSE_LIBRARY_AFTER_FIXED_GENE_SELECTION"
  list(counts = x, audit = audit, skip = skip)
}

satb2_pb_export_results <- function(stable, genes, outdir, auditdir = file.path(outdir, "AUDIT")) {
  if (!dir.exists(auditdir) && !dir.create(auditdir, recursive = TRUE)) satb2_pb_fail("Cannot create audit directory.")
  wanted <- c("gene_id", satb2_pb_result_columns()[-1L])
  if (!is.data.frame(stable) || !identical(names(stable), wanted) || anyDuplicated(stable$gene_id) ||
      any(!stable$gene_id %in% genes$gene_id)) satb2_pb_fail("Invalid stable-ID result schema.")
  gg <- genes[match(stable$gene_id, genes$gene_id), , drop = FALSE]
  readable <- data.frame(gene_symbol = gg$display_symbol, stable[, wanted[-1L], drop = FALSE],
                         stringsAsFactors = FALSE, check.names = FALSE)
  if (!identical(names(readable), satb2_pb_result_columns())) satb2_pb_fail("Seven-column export schema changed.")
  sig <- !is.na(stable$padj) & stable$padj < .05
  masks <- list(ALL_RESULTS = rep(TRUE, nrow(stable)), DEG_FDR05 = sig)
  for (name in names(masks)) {
    ii <- which(masks[[name]])
    satb2_pb_write_table(readable[ii, , drop = FALSE], file.path(outdir, paste0(name, ".txt")))
    satb2_pb_write_table(data.frame(export_row = seq_along(ii), stable[ii, , drop = FALSE],
                                   source_gene_symbol = gg$gene_symbol[ii],
                                   display_symbol = gg$display_symbol[ii],
                                   gene_symbol_missing_used_gene_id = gg$symbol_missing[ii], check.names = FALSE),
                         file.path(auditdir, paste0(name, "_STABLE_ID_AUDIT.tsv")))
  }
  invisible(readable)
}

satb2_pb_export_edger <- function(native, genes, outdir, auditdir = file.path(outdir, "AUDIT")) {
  wanted <- c("gene_id", "logFC", "logCPM", "F", "PValue", "FDR")
  if (!is.data.frame(native) || !identical(names(native), wanted) || anyDuplicated(native$gene_id) ||
      any(!native$gene_id %in% genes$gene_id)) satb2_pb_fail("Invalid native edgeR QL result schema.")
  if (!dir.exists(auditdir) && !dir.create(auditdir, recursive = TRUE)) satb2_pb_fail("Cannot create audit directory.")
  gg <- genes[match(native$gene_id, genes$gene_id), , drop = FALSE]
  readable <- data.frame(gene_symbol = gg$display_symbol, native[, wanted[-1L], drop = FALSE],
                         stringsAsFactors = FALSE, check.names = FALSE)
  masks <- list(ALL_RESULTS = rep(TRUE, nrow(native)), DEG_FDR05 = !is.na(native$FDR) & native$FDR < .05)
  for (name in names(masks)) {
    ii <- which(masks[[name]])
    satb2_pb_write_table(readable[ii, , drop = FALSE], file.path(outdir, paste0(name, ".txt")))
    satb2_pb_write_table(data.frame(export_row = seq_along(ii), native[ii, , drop = FALSE],
                                   source_gene_symbol = gg$gene_symbol[ii], display_symbol = gg$display_symbol[ii],
                                   gene_symbol_missing_used_gene_id = gg$symbol_missing[ii], check.names = FALSE),
                         file.path(auditdir, paste0(name, "_STABLE_ID_AUDIT.tsv")))
  }
  invisible(readable)
}

satb2_pb_fit_deseq2 <- function(x, selection, independent_filtering) {
  warnings <- character()
  collect <- function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") }
  dds <- withCallingHandlers({
    ds <- DESeq2::DESeqDataSetFromMatrix(x, selection$metadata, selection$formula)
    DESeq2::DESeq(ds, test = "Wald", sfType = "poscounts", fitType = "parametric",
                  betaPrior = FALSE, minReplicatesForReplace = Inf, parallel = FALSE, quiet = TRUE)
  }, warning = collect)
  rr <- withCallingHandlers(DESeq2::results(dds, contrast = satb2_pb_contrast(), alpha = .05,
                                           independentFiltering = independent_filtering,
                                           lfcThreshold = 0, pAdjustMethod = "BH", parallel = FALSE), warning = collect)
  stable <- data.frame(gene_id = rownames(rr), as.data.frame(rr)[, satb2_pb_result_columns()[-1L], drop = FALSE],
                       stringsAsFactors = FALSE, check.names = FALSE)
  sf <- DESeq2::sizeFactors(dds)
  if (anyNA(sf) || any(!is.finite(sf)) || any(sf <= 0)) satb2_pb_fail("Invalid DESeq2 size factors.")
  list(results = stable, size_factors = setNames(sf, colnames(dds)), warnings = unique(warnings),
       fit_type_actual = attr(DESeq2::dispersionFunction(dds), "fitType"))
}

satb2_pb_fit_edger <- function(x, selection) {
  warnings <- character()
  collect <- function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") }
  native <- withCallingHandlers({
    y <- edgeR::DGEList(counts = x)
    y <- edgeR::normLibSizes(y, method = "TMM")
    y <- edgeR::estimateDisp(y, selection$model_matrix, robust = TRUE)
    fit <- edgeR::glmQLFit(y, selection$model_matrix, robust = TRUE)
    k <- match("conditioncKO", colnames(selection$model_matrix))
    if (is.na(k)) satb2_pb_fail("Missing cKO-versus-WT model coefficient.")
    tt <- edgeR::topTags(edgeR::glmQLFTest(fit, coef = k), n = Inf,
                        adjust.method = "BH", sort.by = "none")$table
    data.frame(gene_id = rownames(tt), tt, row.names = NULL, check.names = FALSE)
  }, warning = collect)
  stable <- data.frame(gene_id = native$gene_id, baseMean = NA_real_, log2FoldChange = native$logFC,
                       lfcSE = NA_real_, stat = native$F, pvalue = native$PValue, padj = native$FDR,
                       stringsAsFactors = FALSE, check.names = FALSE)
  list(results = stable, native = native, warnings = unique(warnings))
}

satb2_pb_run_endpoint <- function(y, md, genes, universe, outdir, method, endpoint,
                                  omit = NULL, independent_filtering = TRUE,
                                  auditdir = file.path(outdir, "AUDIT"), count_kind = "original_CellRanger_integer_UMI") {
  if (!dir.create(outdir, recursive = TRUE, showWarnings = FALSE)) satb2_pb_fail("Endpoint output already exists/unwritable: ", outdir)
  if (!dir.create(auditdir, recursive = TRUE, showWarnings = FALSE)) satb2_pb_fail("Endpoint audit already exists/unwritable: ", auditdir)
  selection <- satb2_pb_select(md, endpoint, omit)
  prepared <- satb2_pb_prepare_endpoint(y, selection, universe)
  satb2_pb_write_table(selection$audit, file.path(auditdir, "MICE_ALL12_AUDIT.tsv"))
  partition_replication <- as.data.frame(table(selection$metadata$sequencing_partition,
                                              selection$metadata$condition), stringsAsFactors = FALSE)
  names(partition_replication) <- c("sequencing_partition", "condition", "n_included_mice")
  satb2_pb_write_table(partition_replication, file.path(auditdir, "PARTITION_BY_GENOTYPE_REPLICATION.tsv"))
  satb2_pb_write_table(prepared$audit, file.path(auditdir, "GENE_UNIVERSE_AUDIT.tsv"))
  settings <- c(method = method, endpoint = endpoint, design = selection$design, numerator = "cKO", denominator = "WT",
                biological_unit = "independent_mouse", counts = count_kind,
                min_nuclei_per_mouse = "20", min_mice_per_genotype = "3", min_residual_df = "2",
                gene_universe = "fixed_primary_prefilter_genes; endpoint_all_zero_genes_recorded_and_removed",
                p_adjustment = "BH within this endpoint only", alpha = "0.05", seed_policy = "caller_recorded",
                omitted_mouse = if (is.null(omit)) "none" else omit,
                independentFiltering = if (method == "DESeq2") as.character(independent_filtering) else "not_applicable",
                cooksCutoff = if (method == "DESeq2") "DESeq2_default; NA_retained" else "not_applicable",
                automatic_count_replacement = "disabled",
                normalization = if (method == "DESeq2") "poscounts" else "TMM",
                statistic = if (method == "DESeq2") "signed_Wald" else "unsigned_quasi_likelihood_F",
                export_fields = if (method == "DESeq2") "gene_symbol_baseMean_log2FoldChange_lfcSE_stat_pvalue_padj" else "gene_symbol_logFC_logCPM_F_PValue_FDR",
                robust_QL = if (method == "edgeR") "TRUE" else "not_applicable")
  satb2_pb_write_table(data.frame(setting = names(settings), value = unname(settings)), file.path(auditdir, "MODEL_SETTINGS.tsv"))
  row <- data.frame(method = method, endpoint = endpoint, omitted_mouse = selection$omitted_mouse,
                    status = "NOT_TESTED", reason = prepared$skip, n_WT = selection$n_WT, n_cKO = selection$n_cKO,
                    n_nuclei = sum(selection$metadata$n_nuclei), residual_df = selection$residual_df,
                    n_frozen_primary_genes = length(universe), n_endpoint_genes = nrow(prepared$counts),
                    n_pvalue = NA_integer_, n_padj = NA_integer_, DEG_FDR05 = NA_integer_,
                    warning_count = 0L, design = selection$design, output_directory = outdir, stringsAsFactors = FALSE)
  fit <- NULL
  if (!nzchar(prepared$skip)) {
    fit <- tryCatch(if (method == "DESeq2") satb2_pb_fit_deseq2(prepared$counts, selection, independent_filtering)
                    else if (method == "edgeR") satb2_pb_fit_edger(prepared$counts, selection)
                    else satb2_pb_fail("Unknown method."), error = function(e) e)
    if (inherits(fit, "error")) {
      row$status <- "ERROR"; row$reason <- conditionMessage(fit)
      satb2_pb_write_text(row$reason, file.path(auditdir, "ERROR.txt")); fit <- NULL
    } else {
      if (method == "DESeq2") satb2_pb_export_results(fit$results, genes, outdir, auditdir)
      else satb2_pb_export_edger(fit$native, genes, outdir, auditdir)
      row$status <- "TESTED_REVIEW_REQUIRED"; row$reason <- ""
      row$n_pvalue <- sum(!is.na(fit$results$pvalue)); row$n_padj <- sum(!is.na(fit$results$padj))
      row$DEG_FDR05 <- sum(!is.na(fit$results$padj) & fit$results$padj < .05)
      row$warning_count <- length(fit$warnings)
      satb2_pb_write_text(if (length(fit$warnings)) fit$warnings else "No captured fitting warnings.", file.path(auditdir, "MODEL_WARNINGS.txt"))
      if (method == "DESeq2") {
        normalized <- matrix(NA_real_, nrow(y), ncol(y), dimnames = dimnames(y))
        sf <- fit$size_factors
        normalized[, names(sf)] <- sweep(y[, names(sf), drop = FALSE], 2L, sf, "/")
        satb2_pb_write_table(data.frame(gene_id = genes$gene_id, gene_symbol = genes$display_symbol,
                                       normalized, check.names = FALSE), file.path(outdir, "NORMALIZED_UMI_ALL_GENES_ALL12.tsv"))
        satb2_pb_write_table(data.frame(sample = names(sf), fitted_size_factor = unname(sf)), file.path(auditdir, "FITTED_SIZE_FACTORS.tsv"))
        satb2_pb_write_text(paste0("Requested parametric; actual fit type: ", paste(fit$fit_type_actual, collapse = ",")),
                           file.path(auditdir, "DISPERSION_FIT_TYPE.txt"))
      }
    }
  }
  satb2_pb_write_table(row, file.path(auditdir, "EXECUTION_AUDIT.tsv"))
  list(audit = row, results = if (is.null(fit)) NULL else fit$results)
}

run_mouse_pseudobulk <- function(counts, metadata, genes, populations, outdir, seed = 20261010) {
  input <- satb2_pb_validate_input(counts, metadata, genes, populations)
  needed <- c("Matrix", "DESeq2", "edgeR", "statmod", "S4Vectors", "SummarizedExperiment")
  absent <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(absent)) satb2_pb_fail("Required packages unavailable; no automatic installation: ", paste(absent, collapse = ", "))
  if (length(seed) != 1L || !is.numeric(seed) || is.na(seed) || !is.finite(seed) || seed != floor(seed) || seed < 0 || seed > .Machine$integer.max)
    satb2_pb_fail("seed must be one finite nonnegative integer.")
  if (length(outdir) != 1L || is.na(outdir) || !nzchar(outdir)) satb2_pb_fail("outdir must be a fresh dedicated directory.")
  if (dir.exists(outdir) && length(list.files(outdir, all.files = TRUE, no.. = TRUE)))
    satb2_pb_fail("Refusing to mix with or overwrite existing outputs: ", outdir)
  if (!dir.exists(outdir) && !dir.create(outdir, recursive = TRUE)) satb2_pb_fail("Cannot create output directory.")
  outdir <- normalizePath(outdir, mustWork = TRUE)
  audit_root <- file.path(outdir, "AUDIT")
  if (!dir.create(audit_root)) satb2_pb_fail("Cannot create central audit directory.")
  lock <- file.path(outdir, ".execution.lock")
  if (!dir.create(lock, showWarnings = FALSE)) satb2_pb_fail("Execution lock exists; inspect ownership before restart.")
  token <- paste(Sys.getpid(), satb2_pb_utc(), sep = "_")
  writeLines(token, file.path(lock, "owner.txt"))
  on.exit({
    owned <- tryCatch(identical(readLines(file.path(lock, "owner.txt"), warn = FALSE), token), error = function(e) FALSE)
    if (owned) unlink(lock, recursive = TRUE)
  }, add = TRUE)
  satb2_pb_write_table(data.frame(package = c("R", needed),
                                 version = c(as.character(getRversion()), vapply(needed, function(p) as.character(utils::packageVersion(p)), character(1)))),
                       file.path(audit_root, "R_PACKAGE_VERSIONS.tsv"))
  satb2_pb_write_text(capture.output(utils::sessionInfo()), file.path(audit_root, "SESSION_INFO.txt"))
  satb2_pb_write_table(input$mouse, file.path(audit_root, "CONFIRMED_ALL12_MOUSE_ROSTER.tsv"))
  satb2_pb_write_table(input$genes, file.path(audit_root, "STABLE_GENE_ID_SYMBOL_AUDIT.tsv"))
  duplicate <- !input$genes$symbol_missing & (duplicated(input$genes$gene_symbol) | duplicated(input$genes$gene_symbol, fromLast = TRUE))
  satb2_pb_write_table(input$genes[duplicate, , drop = FALSE], file.path(audit_root, "DUPLICATE_SYMBOL_STABLE_IDS.tsv"))
  satb2_pb_write_text(c("Execution evidence requires scientific review; this module does not declare project completion.",
                        "Each nucleus belongs to its confirmed independent mouse; nuclei are never biological replicates.",
                        "Population lists may overlap/nest. No anatomical layer labels or marker gates are imposed.",
                        "Original all-12 raw zero columns identify absent populations; normalization/model ineligibility is exported as NA/audit status.",
                        "Primary DESeq2 and edgeR share a frozen input gene universe, but their realized tested/BH families may differ.",
                        "Each endpoint has separate BH adjustment. There is no pooled-method or across-population significance claim.",
                        "LOO is influence evidence and does not authorize data-dependent mouse exclusion.",
                        "Readable DESeq2 exports contain its standard seven columns; edgeR exports its native fields.",
                        "No directional-split or NOT_TESTED marker files are produced; status evidence remains in AUDIT.",
                        "Caller must freeze reviewed original-count, metadata and population source identities before execution.",
                        paste0("Seed recorded: ", seed, "; these count-model fits do not require random sampling.")), file.path(audit_root, "METHODS_AND_LIMITS.txt"))
  audits <- list(); coverage <- list(); memberships <- list(); conservation <- list()
  for (pop in names(input$populations)) {
    dir <- file.path(outdir, pop)
    if (!dir.create(dir)) satb2_pb_fail("Cannot create population directory: ", pop)
    ids <- input$populations[[pop]]
    pb <- satb2_pb_aggregate(input$counts, input$metadata, ids, input$mouse)
    y <- pb$counts; md <- pb$metadata
    coverage[[pop]] <- data.frame(population = pop, md, row.names = NULL)
    mm <- input$metadata[match(ids, input$metadata$nucleus_id), c("nucleus_id", "sample", "condition"), drop = FALSE]
    memberships[[pop]] <- data.frame(population = rep(pop, nrow(mm)), mm, row.names = NULL)
    conservation[[pop]] <- data.frame(population = pop, n_nuclei = length(ids), total_original_UMIs = pb$conserved_total_UMIs,
                                      per_gene_mouse_sum_conserved = TRUE)
    satb2_pb_write_table(data.frame(gene_id = input$genes$gene_id, gene_symbol = input$genes$display_symbol,
                                   y, check.names = FALSE), file.path(dir, "ORIGINAL_RAW_UMI_ALL_GENES_ALL12.tsv"))
    satb2_pb_write_table(md, file.path(dir, "ALL12_POPULATION_COVERAGE.tsv"))
    primary_selection <- satb2_pb_select(md, "PRIMARY")
    pf <- satb2_pb_prefilter(y, primary_selection)
    pop_audit <- file.path(audit_root, pop)
    if (!dir.create(pop_audit)) satb2_pb_fail("Cannot create population audit directory.")
    satb2_pb_write_table(pf, file.path(pop_audit, "FROZEN_PRIMARY_GENE_PREFILTER_AUDIT.tsv"))
    universe <- pf$gene_id[pf$passes_primary_prefilter]
    run <- function(method, endpoint, omit = NULL, independent_filtering = TRUE) {
      label <- if (endpoint == "LOO") paste0("LOO_omit_", omit) else endpoint
      ans <- satb2_pb_run_endpoint(y, md, input$genes, universe, file.path(dir, method, label),
                                   method, endpoint, omit, independent_filtering,
                                   auditdir = file.path(pop_audit, method, label))
      ans$audit$population <- pop
      audits[[paste(pop, method, label, sep = "__")]] <<- ans$audit
      ans
    }
    primary <- run("DESeq2", "PRIMARY")
    run("edgeR", "PRIMARY")
    run("DESeq2", "KNOWN_PARTITION_ADJUSTED")
    run("DESeq2", "MAIN_PARTITION_ONLY")
    loo <- list()
    for (omit in md$sample) {
      ans <- run("DESeq2", "LOO", omit, independent_filtering = FALSE)
      if (!is.null(ans$results)) {
        z <- ans$results
        z$omitted_mouse <- omit
        loo[[omit]] <- z
      }
    }
    if (length(loo) && !is.null(primary$results)) {
      base <- primary$results; lfc <- vapply(loo, function(z) z$log2FoldChange[match(base$gene_id, z$gene_id)], numeric(nrow(base)))
      if (is.null(dim(lfc))) lfc <- matrix(lfc, nrow = nrow(base))
      finite <- is.finite(lfc)
      st <- data.frame(gene_id = base$gene_id, gene_symbol = input$genes$display_symbol[match(base$gene_id, input$genes$gene_id)],
                        primary_log2FC = base$log2FoldChange, primary_padj = base$padj,
                        n_LOO_estimable = rowSums(finite),
                        n_LOO_same_direction = rowSums(finite & sign(lfc) == sign(base$log2FoldChange), na.rm = TRUE),
                        min_LOO_log2FC = apply(lfc, 1L, function(z) if (any(is.finite(z))) min(z[is.finite(z)]) else NA_real_),
                        max_LOO_log2FC = apply(lfc, 1L, function(z) if (any(is.finite(z))) max(z[is.finite(z)]) else NA_real_),
                        stringsAsFactors = FALSE)
      satb2_pb_write_table(st, file.path(pop_audit, "LOO_EFFECT_INFLUENCE_ONLY.tsv"))
    }
  }
  audit <- do.call(rbind, audits); rownames(audit) <- NULL
  satb2_pb_write_table(audit, file.path(audit_root, "EXECUTION_AUDIT_REPORT.tsv"))
  satb2_pb_write_table(do.call(rbind, coverage), file.path(outdir, "ALL_POPULATION_ALL12_COVERAGE.tsv"))
  satb2_pb_write_table(do.call(rbind, memberships), file.path(audit_root, "POPULATION_MEMBERSHIP_ID_AUDIT.tsv"))
  satb2_pb_write_table(do.call(rbind, conservation), file.path(audit_root, "ORIGINAL_COUNT_CONSERVATION_AUDIT.tsv"))
  report <- list(status = if (any(audit$status == "ERROR")) "EXECUTION_FINISHED_WITH_ERRORS_REVIEW_REQUIRED" else "EXECUTED_REVIEW_REQUIRED",
                  endpoint_errors = sum(audit$status == "ERROR"),
                  endpoints_not_tested = sum(audit$status == "NOT_TESTED"), audit = audit,
                  source_n_nuclei = ncol(input$counts), source_n_genes = nrow(input$counts),
                  n_populations = length(input$populations), seed = seed, completed_execution_utc = satb2_pb_utc(),
                  outdir = outdir, scientific_acceptance = "NOT_ASSESSED")
  saveRDS(report, file.path(audit_root, "EXECUTION_REPORT.rds"), compress = TRUE)
  report
}
