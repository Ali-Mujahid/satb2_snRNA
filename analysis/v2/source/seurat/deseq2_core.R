# DESeq2-only 14 October revision. Source defines functions; never fits on source().
# Mouse pseudobulk; WT baseline and cKO/WT signed Wald effects.
# Stable IDs define modeled rows in AUDIT; readable exports use gene symbols.
# Missing symbols are NA. Duplicate symbols remain distinct modeled rows.
# IN: every nonempty mouse aggregate, >=5 counts in >=3 eligible mice. Others: >=20 nuclei/mouse,
# >=10 counts in min(eligible WT, eligible cKO) mice. All models require >=3/group,
# full-rank design and >=2 residual df. BH is within each endpoint.

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
satb2_pb_result_columns <- function() c("gene", "baseMean", "log2FoldChange",
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
  genes$display_symbol <- ifelse(genes$symbol_missing, NA_character_, genes$gene_symbol)
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

satb2_pb_select <- function(md, endpoint, omit = NULL, min_nuclei = 20L) {
  audit <- md
  if (length(min_nuclei) != 1L || !min_nuclei %in% c(1L, 20L)) satb2_pb_fail("Unknown declared nuclei policy")
  eligible <- audit$n_nuclei >= min_nuclei & audit$total_UMIs > 0
  reason <- ifelse(audit$n_nuclei == 0L, "ABSENT_POPULATION",
                   ifelse(audit$n_nuclei < min_nuclei, paste0("FEWER_THAN_", min_nuclei, "_NUCLEI"), "INCLUDED"))
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

satb2_pb_prefilter <- function(y, selection, inhibitory = FALSE) {
  xx <- y[, selection$included, drop = FALSE]
  count_cutoff <- if (inhibitory) 5L else 10L
  required <- if (inhibitory) 3L else min(selection$n_WT, selection$n_cKO)
  n <- rowSums(xx >= count_cutoff)
  data.frame(gene_id = rownames(y), total_UMIs_all12 = rowSums(y),
    count_cutoff_per_mouse = count_cutoff, n_included_mice_above_count_cutoff = n,
    required_mice = required, passes_primary_prefilter = required > 0L & n >= required,
    stringsAsFactors = FALSE)
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
  readable <- data.frame(gene = gg$display_symbol, stable[, wanted[-1L], drop = FALSE],
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
                                   gene_symbol_missing_exported_NA = gg$symbol_missing[ii], check.names = FALSE),
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

# Cache reuse is mathematical reuse of an unchanged endpoint, not a fresh fit.
# Every selected count, gene, mouse, formula, contrast, package version and
# DESeq2 result policy must agree. Candidate paths never authorize reuse alone.
satb2_pb_cache <- function(candidate, y, selection, prepared, universe, independent_filtering, membership_sha256 = NULL) {
  no <- function(reason) list(accepted = FALSE, reason = reason)
  if (is.null(candidate)) return(no("NO_CANDIDATE"))
  need <- c("counts", "gene_prefilter", "gene_universe", "mice", "settings", "results", "execution", "versions", "session_info", "statistical_core", "size_factors", "fit_type", "warnings", "membership_metadata", "branch_roster")
  if (!is.list(candidate) || !all(need %in% names(candidate))) return(no("INCOMPLETE_CANDIDATE_CONTRACT"))
  if (any(!vapply(candidate[need], file.exists, logical(1)))) return(no("CANDIDATE_FILE_MISSING"))
  if (is.null(membership_sha256) || is.null(candidate$membership_sha256) ||
      !identical(as.character(candidate$membership_sha256), as.character(membership_sha256))) return(no("SOURCE_POPULATION_MEMBERSHIP_DIFFERENT"))
  read <- function(p) utils::read.delim(p, check.names = FALSE, stringsAsFactors = FALSE)
  truth <- function(x) !is.na(x) & tolower(as.character(x)) %in% c("true", "1", "yes")
  tryCatch({
    old <- read(candidate$counts)
    if (!all(c("gene_id", colnames(y)) %in% names(old)) || anyDuplicated(old$gene_id) ||
        !identical(as.character(old$gene_id), rownames(y))) return(no("COUNT_GENE_AXIS_DIFFERENT"))
    old_y <- as.matrix(old[, colnames(y), drop = FALSE]); storage.mode(old_y) <- "double"
    if (!identical(as.numeric(old_y), as.numeric(y))) return(no("ALL12_COUNTS_DIFFERENT"))
    pf <- read(candidate$gene_prefilter)
    if (!all(c("gene_id", "passes_primary_prefilter") %in% names(pf)) ||
        !identical(sort(as.character(pf$gene_id[truth(pf$passes_primary_prefilter)])), sort(universe))) return(no("FROZEN_UNIVERSE_DIFFERENT"))
    gu <- read(candidate$gene_universe)
    if (!all(c("gene_id", "included_gene_in_endpoint") %in% names(gu)) ||
        !identical(as.character(gu$gene_id[truth(gu$included_gene_in_endpoint)]), rownames(prepared$counts))) return(no("REALIZED_GENE_UNIVERSE_DIFFERENT"))
    mice <- read(candidate$mice)
    if (!identical(as.character(mice$sample), as.character(selection$audit$sample)) ||
        !identical(truth(mice$included_in_model), selection$included) ||
        !identical(as.character(mice$condition), as.character(selection$audit$condition)) ||
        !identical(as.character(mice$sequencing_partition), as.character(selection$audit$sequencing_partition))) return(no("MOUSE_SELECTION_OR_METADATA_DIFFERENT"))
    ss <- read(candidate$settings); settings <- setNames(as.character(ss$value), ss$setting)
    expect <- c(method = "DESeq2", design = selection$design, numerator = "cKO", denominator = "WT",
      biological_unit = "independent_mouse", independentFiltering = as.character(independent_filtering),
      cooksCutoff = "DESeq2_default; NA_retained", automatic_count_replacement = "disabled",
      normalization = "poscounts", statistic = "signed_Wald", alpha = "0.05")
    if (anyNA(settings[names(expect)]) || any(settings[names(expect)] != expect)) return(no("MODEL_SETTINGS_DIFFERENT"))
    vv <- read(candidate$versions)
    if (!all(c("package", "version") %in% names(vv))) return(no("VERSION_SCHEMA_DIFFERENT"))
    pkgs <- c("R", "DESeq2", "S4Vectors", "SummarizedExperiment", "Matrix")
    actual <- c(as.character(getRversion()), vapply(pkgs[-1L], function(p) as.character(utils::packageVersion(p)), character(1)))
    prior <- as.character(vv$version[match(pkgs, vv$package)])
    session <- paste(readLines(candidate$session_info, warn = FALSE), collapse = " ")
    for (i in which(is.na(prior))) {
      token <- regmatches(session, regexpr(paste0(pkgs[i], "_[0-9][0-9.-]*"), session))
      if (length(token) && nzchar(token)) prior[i] <- sub(paste0("^", pkgs[i], "_"), "", token)
    }
    if (anyNA(prior) || any(prior != actual)) return(no("MODEL_RUNTIME_VERSION_DIFFERENT"))
    # This specific audited core binds legacy settings omitted from its compact
    # settings table: parametric fit request, Wald, betaPrior=FALSE, replacement
    # disabled, explicit contrast, and parallel=FALSE. Unknown source versions
    # are recomputed through the native revision instead of assumed equivalent.
    approved <- c("319f799d69f3f84d837c2261bd7cd62934e235d7eb8a95e94da5cff69e1ad126",
      getOption("satb2.revision.statistical_core_sha256", character()))
    if (!digest::digest(file = candidate$statistical_core, algo = "sha256") %in% approved)
      return(no("UNRECOGNIZED_AUDITED_LEGACY_STATISTICAL_CORE"))
    audit <- read(candidate$execution)
    if (nrow(audit) != 1L || audit$status != "TESTED_REVIEW_REQUIRED" || audit$warning_count != 0L ||
        audit$n_WT != selection$n_WT || audit$n_cKO != selection$n_cKO ||
        audit$residual_df != selection$residual_df || audit$design != selection$design) return(no("PRIOR_EXECUTION_GUARD_DIFFERENT"))
    z <- read(candidate$results)
    wanted <- c("gene_id", satb2_pb_result_columns()[-1L])
    if (!all(wanted %in% names(z)) || anyDuplicated(z$gene_id) ||
        !identical(as.character(z$gene_id), rownames(prepared$counts))) return(no("RESULT_GENE_AXIS_DIFFERENT"))
    z <- z[, wanted, drop = FALSE]
    if (sum(!is.na(z$pvalue)) != audit$n_pvalue || sum(!is.na(z$padj)) != audit$n_padj ||
        sum(!is.na(z$padj) & z$padj < .05) != audit$DEG_FDR05) return(no("PRIOR_RESULT_DENOMINATOR_DIFFERENT"))
    sf <- read(candidate$size_factors)
    if (!identical(as.character(sf$sample), colnames(prepared$counts)) ||
        any(!is.finite(sf$fitted_size_factor)) || any(sf$fitted_size_factor <= 0)) return(no("SIZE_FACTOR_GUARD_DIFFERENT"))
    list(accepted = TRUE, reason = "EXACT_UNCHANGED_COUNTS_GENE_UNIVERSE_DESIGN_RUNTIME_AND_RESULT_POLICY",
      fit = list(results = z, size_factors = setNames(sf$fitted_size_factor, sf$sample), warnings = character(),
        fit_type_actual = paste(readLines(candidate$fit_type, warn = FALSE), collapse = " ")),
      receipt = data.frame(role = need, path = unlist(candidate[need], use.names = FALSE),
        sha256 = vapply(candidate[need], function(p) digest::digest(file = p, algo = "sha256"), character(1)),
        stringsAsFactors = FALSE))
  }, error = function(e) no(paste0("CACHE_GUARD_ERROR:", conditionMessage(e))))
}

satb2_pb_run_endpoint <- function(y, md, genes, universe, outdir, endpoint = "PRIMARY",
    omit = NULL, independent_filtering = TRUE, count_kind, min_nuclei = 20L,
    gene_count_cutoff = 10L, gene_required_mice = NA_integer_, cache = NULL, membership_sha256 = NULL) {
  if (dir.exists(outdir) && length(list.files(outdir, all.files = TRUE, no.. = TRUE)))
    satb2_pb_fail("Refusing existing endpoint content: ", outdir)
  if (!dir.exists(outdir) && !dir.create(outdir, recursive = TRUE)) satb2_pb_fail("Cannot create endpoint: ", outdir)
  auditdir <- file.path(outdir, "AUDIT")
  if (!dir.create(auditdir)) satb2_pb_fail("Cannot create endpoint audit: ", auditdir)
  selection <- satb2_pb_select(md, endpoint, omit, min_nuclei)
  prepared <- satb2_pb_prepare_endpoint(y, selection, universe)
  satb2_pb_write_table(selection$audit, file.path(auditdir, "MICE_ALL12_AUDIT.tsv"))
  replication <- as.data.frame(table(selection$metadata$sequencing_partition, selection$metadata$condition), stringsAsFactors = FALSE)
  names(replication) <- c("sequencing_partition", "condition", "n_included_mice")
  satb2_pb_write_table(replication, file.path(auditdir, "PARTITION_BY_GENOTYPE_REPLICATION.tsv"))
  satb2_pb_write_table(prepared$audit, file.path(auditdir, "GENE_UNIVERSE_AUDIT.tsv"))
  settings <- c(method = "DESeq2", endpoint = endpoint, design = selection$design,
    numerator = "cKO", denominator = "WT", biological_unit = "independent_mouse", counts = count_kind,
    min_nuclei_per_mouse = as.character(min_nuclei), min_mice_per_genotype = "3", min_residual_df = "2",
    primary_gene_count_cutoff = as.character(gene_count_cutoff),
    primary_gene_required_mice = as.character(gene_required_mice),
    gene_universe = "fixed_original_primary_prefilter; endpoint_all_zero_genes_audited_and_removed",
    p_adjustment = "BH within this endpoint only", alpha = "0.05",
    omitted_mouse = if (is.null(omit)) "none" else omit,
    independentFiltering = as.character(independent_filtering), cooksCutoff = "DESeq2_default; NA_retained",
    automatic_count_replacement = "disabled", normalization = "poscounts", statistic = "signed_Wald",
    test = "Wald", fitType_requested = "parametric", betaPrior = "FALSE", parallel = "FALSE",
    export_fields = paste(satb2_pb_result_columns(), collapse = ","),
    missing_symbol_policy = "NA; stable IDs only in AUDIT", duplicate_symbol_policy = "separate modeled stable-ID rows retained")
  satb2_pb_write_table(data.frame(setting = names(settings), value = unname(settings)), file.path(auditdir, "MODEL_SETTINGS.tsv"))
  row <- data.frame(method = "DESeq2", endpoint = endpoint, omitted_mouse = selection$omitted_mouse,
    status = "NOT_TESTED", reason = prepared$skip, n_WT = selection$n_WT, n_cKO = selection$n_cKO,
    n_nuclei = sum(selection$metadata$n_nuclei), residual_df = selection$residual_df,
    n_frozen_primary_genes = length(universe), n_endpoint_genes = nrow(prepared$counts),
    n_pvalue = NA_integer_, n_padj = NA_integer_, DEG_FDR05 = NA_integer_, warning_count = 0L,
    design = selection$design, execution_origin = "NOT_FITTED", cache_reason = "NOT_EVALUATED",
    output_directory = outdir, stringsAsFactors = FALSE)
  row$min_included_mouse_nuclei <- if (nrow(selection$metadata)) min(selection$metadata$n_nuclei) else NA_integer_
  row$max_included_mouse_nuclei <- if (nrow(selection$metadata)) max(selection$metadata$n_nuclei) else NA_integer_
  row$coverage_flag <- if (min_nuclei == 1L && nrow(selection$metadata) && any(selection$metadata$n_nuclei < 5L))
    "LOW_NUCLEUS_COVERAGE_EXPLORATORY" else "STANDARD_DECLARED_COVERAGE"
  fit <- NULL
  if (!nzchar(prepared$skip)) {
    reused <- satb2_pb_cache(cache, y, selection, prepared, universe, independent_filtering, membership_sha256)
    row$cache_reason <- reused$reason
    if (isTRUE(reused$accepted)) {
      fit <- reused$fit; row$execution_origin <- "REUSED_EXACT_UNCHANGED_DESEQ2_ENDPOINT"
      satb2_pb_write_table(reused$receipt, file.path(auditdir, "REUSED_ENDPOINT_SOURCE_HASHES.tsv"))
    } else {
      row$execution_origin <- "FRESH_NATIVE_DESEQ2_FIT"
      fit <- tryCatch(satb2_pb_fit_deseq2(prepared$counts, selection, independent_filtering), error = function(e) e)
    }
    if (inherits(fit, "error")) {
      row$status <- "ERROR"; row$reason <- conditionMessage(fit)
      satb2_pb_write_text(row$reason, file.path(auditdir, "ERROR.txt")); fit <- NULL
    } else {
      satb2_pb_export_results(fit$results, genes, outdir, auditdir)
      row$status <- "TESTED_REVIEW_REQUIRED"; row$reason <- ""
      row$n_pvalue <- sum(!is.na(fit$results$pvalue)); row$n_padj <- sum(!is.na(fit$results$padj))
      row$DEG_FDR05 <- sum(!is.na(fit$results$padj) & fit$results$padj < .05)
      row$warning_count <- length(fit$warnings)
      satb2_pb_write_text(if (length(fit$warnings)) fit$warnings else "No captured fitting warnings.", file.path(auditdir, "MODEL_WARNINGS.txt"))
      satb2_pb_write_table(data.frame(sample = names(fit$size_factors), fitted_size_factor = unname(fit$size_factors)), file.path(auditdir, "FITTED_SIZE_FACTORS.tsv"))
      satb2_pb_write_text(paste(fit$fit_type_actual, collapse = ","), file.path(auditdir, "DISPERSION_FIT_TYPE.txt"))
    }
  }
  # Publication barrier: no consumer treats an endpoint as complete before this.
  satb2_pb_write_table(row, file.path(auditdir, "EXECUTION_AUDIT.tsv"))
  list(audit = row, results = if (is.null(fit)) NULL else fit$results)
}

