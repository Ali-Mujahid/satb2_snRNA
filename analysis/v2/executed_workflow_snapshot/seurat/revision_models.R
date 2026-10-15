# DESeq2-only endpoint orchestration; source() performs no fitting.
legacy_metadata_cache <- new.env(parent = emptyenv())
legacy_narrow_metadata <- function(path, fields) {
  key <- paste(path, paste(fields, collapse = "|"), sep = "\r")
  if (exists(key, envir = legacy_metadata_cache, inherits = FALSE)) return(get(key, envir = legacy_metadata_cache))
  header <- read.delim(path, nrows = 0L, check.names = FALSE)
  required_columns(header, fields, "legacy membership evidence")
  classes <- rep("NULL", ncol(header)); classes[match(fields, names(header))] <- "character"
  tab <- read.delim(path, colClasses = classes, check.names = FALSE, stringsAsFactors = FALSE)
  assign(key, tab, envir = legacy_metadata_cache); tab
}
legacy_candidate <- function(cfg, branch, row, endpoint_label) {
  root <- cfg$deseq2$cache_pseudobulk_root_by_branch[[branch]]
  if (is.null(root)) root <- file.path(cfg$previous_run_root, "results", "seurat", branch, "pseudobulk")
  if (!file.exists(file.path(root, "POPULATION_INDEX.tsv"))) return(NULL)
  index <- read_tsv(file.path(root, "POPULATION_INDEX.tsv"))
  pop_alias <- c(ALL = "allcells", ALLCELLS = "allcells", EX = "Excitatory", EXCITATORY = "Excitatory", IN = "Inhibitory", INHIBITORY = "Inhibitory", CAI = "CA1_ProS", CA1_PROS = "CA1_ProS")
  pop_name <- unname(pop_alias[toupper(row$population)])
  if (is.na(pop_name) || !row$analysis_kind %in% c("total", "cluster")) return(NULL)
  tool <- tolower(row$tool)
  if (row$analysis_kind == "total") {
    family <- paste0("ATLAS_TOTAL_", c(allcells = "ALL", Excitatory = "EX", Inhibitory = "IN", CA1_ProS = "CA1_PROS")[[pop_name]])
    source <- file.path(cfg$previous_run_root, "results", "shared", branch, "nucleus_metadata.tsv")
    flag <- c(allcells = "nucleus_id", Excitatory = "REF_EX_HS", Inhibitory = "REF_IN_HS", CA1_ProS = "REF_CA1_PROS_HS")[[pop_name]]
    metadata <- legacy_narrow_metadata(source, unique(c("nucleus_id", flag)))
    ids <- if (pop_name == "allcells") metadata$nucleus_id else metadata$nucleus_id[truth(metadata[[flag]])]
    hit <- which(index$family == family)
    source_roster <- source
  } else {
    field <- if (tool == "scanpy") "cluster_harmony" else
      c(allcells = "log_harmony_r0p4", Excitatory = "EX_harmony_r0p4", Inhibitory = "IN_harmony_r0p4", CA1_ProS = "CA1_PROS_CANDIDATE_harmony_r0p4")[[pop_name]]
    family <- if (tool == "scanpy") paste0("SCANPY_", toupper(pop_name), "_harmony") else paste0("SEURAT_", field)
    source <- if (tool == "scanpy") file.path(cfg$previous_run_root, "results", "scanpy", "original_full", pop_name, "nucleus_metadata.tsv") else
      if (pop_name == "allcells") file.path(cfg$previous_run_root, "results", "seurat", "original_full", "log_global", "nucleus_metadata.tsv") else
        file.path(cfg$previous_run_root, "results", "seurat", "original_full", "annotate_subsets", c(Excitatory = "EX", Inhibitory = "IN", CA1_ProS = "CA1_PROS_CANDIDATE")[[pop_name]], "nucleus_metadata.tsv")
    if (!file.exists(source)) return(NULL)
    metadata <- legacy_narrow_metadata(source, c("nucleus_id", field))
    source_roster <- file.path(cfg$previous_run_root, "results", "shared", branch, "nucleus_metadata.tsv")
    rr <- legacy_narrow_metadata(source_roster, "nucleus_id")
    ids <- metadata$nucleus_id[metadata[[field]] == as.character(row$analysis_label) & metadata$nucleus_id %in% rr$nucleus_id]
    hit <- which(index$family == family & as.character(index$label) == as.character(row$analysis_label))
  }
  if (length(hit) != 1L) return(NULL)
  pop_key <- index$population_key[hit]
  pop <- file.path(root, pop_key); ad <- file.path(pop, "AUDIT", "DESeq2", endpoint_label)
  list(counts = file.path(pop, "SOURCE_INTEGER_COUNTS_ALL_GENES_ALL12.tsv"),
    gene_prefilter = file.path(pop, "AUDIT", "FROZEN_PRIMARY_GENE_PREFILTER_AUDIT.tsv"),
    gene_universe = file.path(ad, "GENE_UNIVERSE_AUDIT.tsv"), mice = file.path(ad, "MICE_ALL12_AUDIT.tsv"),
    settings = file.path(ad, "MODEL_SETTINGS.tsv"), results = file.path(ad, "ALL_RESULTS_STABLE_ID_AUDIT.tsv"),
    execution = file.path(ad, "EXECUTION_AUDIT.tsv"), versions = file.path(root, "PACKAGE_VERSIONS.tsv"),
    session_info = file.path(root, "SESSION_INFO.txt"),
    statistical_core = cfg$deseq2$legacy_statistical_core %or% file.path(cfg$study_root, "source", "seurat", "pseudobulk_core.R"),
    size_factors = file.path(ad, "FITTED_SIZE_FACTORS.tsv"), fit_type = file.path(ad, "DISPERSION_FIT_TYPE.txt"),
    warnings = file.path(ad, "MODEL_WARNINGS.txt"), membership_metadata = source, branch_roster = source_roster,
    membership_sha256 = member_sha(as.character(ids)))
}
revision_candidate <- function(pop_counts, pop_audit, endpoint, stage_root, source_core, membership_path, index_path, membership_sha256) {
  ad <- file.path(endpoint, "AUDIT")
  list(counts = pop_counts, gene_prefilter = file.path(pop_audit, "FROZEN_PRIMARY_GENE_PREFILTER_AUDIT.tsv"),
    gene_universe = file.path(ad, "GENE_UNIVERSE_AUDIT.tsv"), mice = file.path(ad, "MICE_ALL12_AUDIT.tsv"),
    settings = file.path(ad, "MODEL_SETTINGS.tsv"), results = file.path(ad, "ALL_RESULTS_STABLE_ID_AUDIT.tsv"),
    execution = file.path(ad, "EXECUTION_AUDIT.tsv"), versions = file.path(stage_root, "PACKAGE_VERSIONS.tsv"),
    session_info = file.path(stage_root, "SESSION_INFO.txt"), statistical_core = source_core,
    size_factors = file.path(ad, "FITTED_SIZE_FACTORS.tsv"), fit_type = file.path(ad, "DISPERSION_FIT_TYPE.txt"),
    warnings = file.path(ad, "MODEL_WARNINGS.txt"), membership_metadata = membership_path,
    branch_roster = index_path, membership_sha256 = membership_sha256)
}
population_identity <- function(row) paste(row$tool, row$population, row$analysis_kind, row$category_family, row$analysis_label, sep = "\r")

verify_count_acceptance <- function(sc, branch, finalized) {
  acceptance <- sc$count_acceptance_json_by_branch[[branch]]
  if (is.null(acceptance) || !file.exists(acceptance)) fail("Independent count acceptance receipt required before model fitting")
  ac <- jsonlite::fromJSON(acceptance, simplifyVector = TRUE)
  key <- sc$count_acceptance_bool_key %or% "exact_integer_counts_PASS"
  if (!isTRUE(ac[[key]])) fail("Independent count acceptance did not pass")
  expected_config <- getOption("satb2.revision.master_config_sha256", NULL)
  if (is.null(expected_config) || is.null(ac$config_sha256) || !identical(as.character(ac$config_sha256), as.character(expected_config)))
    fail("Independent count acceptance is not bound to the current immutable master")
  bindings <- ac$validated_sources
  if (is.null(bindings)) fail("Independent count acceptance lacks validated source hashes")
  for (path in c(file.path(finalized, "POPULATION_INDEX.tsv"), file.path(finalized, "POPULATION_ATOMS.tsv"))) {
    recorded <- bindings[[normalizePath(path, mustWork = TRUE)]]
    if (is.null(recorded)) fail("Independent count acceptance lacks the consumed branch population source: ", path)
    hash_guard(path, as.character(recorded), "accepted finalized branch definition")
  }
  invisible(TRUE)
}

run_revision_models <- function(cfg, branch, out, source_core) {
  sc <- cfg$deseq2; finalized <- sc$finalized_root_by_branch[[branch]]
  if (!file.exists(file.path(finalized, "COMPLETE.json"))) fail("Final population producer incomplete")
  verify_count_acceptance(sc, branch, finalized)
  index <- read_tsv(file.path(finalized, "POPULATION_INDEX.tsv")); genes <- read_tsv(file.path(finalized, "GENE_ID_SYMBOL_MANIFEST.tsv"))
  genes$gene_id <- as.character(genes$gene_id); genes$gene_symbol <- as.character(genes$gene_symbol)
  genes$symbol_missing <- is.na(genes$gene_symbol) | !nzchar(trimws(genes$gene_symbol))
  genes$display_symbol <- ifelse(genes$symbol_missing, NA_character_, genes$gene_symbol)
  required_columns(index, c("population_key", "tool", "population", "analysis_kind", "analysis_label", "category_family", "model_enabled", "core_total", "min_nuclei_per_mouse", "counts_file", "counts_sha256"), "model populations")
  samples <- names(satb2_pb_roster()); statuses <- list(); same_run <- list(); universes <- list()
  for (i in seq_len(nrow(index))) {
    row <- index[i, , drop = FALSE]; pop <- row$population_key
    folder <- new_dir(file.path(out, row$tool, pop)); pop_audit <- new_dir(file.path(folder, "AUDIT"))
    hash_guard(row$counts_file, row$counts_sha256, paste0("finalized count array/", pop))
    y <- read_integer_counts(row$counts_file, genes, samples)
    md <- read_tsv(file.path(dirname(row$counts_file), "ALL12_POPULATION_COVERAGE.tsv"))
    required_columns(md, c("sample", "condition", "sequencing_partition", "n_nuclei", "total_UMIs"), "model coverage")
    if (!identical(as.character(md$sample), samples) || any(md$condition != unname(satb2_pb_roster())) ||
        any(md$sequencing_partition != unname(satb2_pb_expected_partitions())) ||
        !identical(as.numeric(colSums(y)), as.numeric(md$total_UMIs))) fail("All12 model coverage/count roster mismatch")
    rownames(md) <- samples
    inhibitory <- toupper(row$population) %in% c("IN", "INHIBITORY")
    min_nuclei <- if (inhibitory) 1L else 20L
    if (row$min_nuclei_per_mouse != min_nuclei) fail("Model nuclei policy changed")
    selection <- satb2_pb_select(md, "PRIMARY", min_nuclei = min_nuclei)
    pf <- satb2_pb_prefilter(y, selection, inhibitory)
    universe <- pf$gene_id[pf$passes_primary_prefilter]; policy <- "declared_primary_original_count_prefilter"
    universe_root <- sc$prefilter_universe_root_by_branch[[branch]]
    if (branch == "paired_cellbender" && is.null(universe_root)) fail("Corrected sensitivity requires paired-original frozen universe")
    if (!is.null(universe_root)) {
      old_index <- read_tsv(file.path(universe_root, "POPULATION_INDEX.tsv"))
      hit <- which(vapply(seq_len(nrow(old_index)), function(j) population_identity(old_index[j, , drop = FALSE]), character(1)) == population_identity(row))
      if (length(hit) != 1L) fail("Frozen universe has no unique matching population identity")
      old_pop <- old_index[hit, , drop = FALSE]
      old_pf <- read_tsv(file.path(universe_root, old_pop$tool, old_pop$population_key, "AUDIT", "FROZEN_PRIMARY_GENE_PREFILTER_AUDIT.tsv"))
      required_columns(old_pf, c("gene_id", "passes_primary_prefilter"), "source frozen universe")
      universe <- as.character(old_pf$gene_id[truth(old_pf$passes_primary_prefilter)])
      if (any(!universe %in% genes$gene_id)) fail("Frozen universe includes unsupported genes")
      pf$passes_primary_prefilter <- pf$gene_id %in% universe; policy <- "frozen_from_declared_original_count_branch"
      write_tsv(data.frame(path = file.path(universe_root, old_pop$tool, old_pop$population_key, "AUDIT", "FROZEN_PRIMARY_GENE_PREFILTER_AUDIT.tsv"),
        sha256 = sha_file(file.path(universe_root, old_pop$tool, old_pop$population_key, "AUDIT", "FROZEN_PRIMARY_GENE_PREFILTER_AUDIT.tsv"))), file.path(pop_audit, "FROZEN_UNIVERSE_SOURCE_HASH.tsv"))
    }
    pf$universe_policy <- policy
    write_tsv(pf, file.path(pop_audit, "FROZEN_PRIMARY_GENE_PREFILTER_AUDIT.tsv"))
    if (!truth(row$model_enabled)) {
      status <- data.frame(method = "DESeq2", endpoint = "PRIMARY", omitted_mouse = NA_character_, status = "NOT_TESTED",
        reason = if (row$n_nuclei == 0L) "EMPTY_ORIGINAL_SOURCE_CLUSTER_AFTER_MASK_OR_PAIRING_RETAINED_IN_AUDIT" else "EXCLUDED_ORIGINAL_ALL_EX_CLUSTER_BELOW100_AFTER_SHARED_INTERSECTION", n_WT = selection$n_WT, n_cKO = selection$n_cKO,
        n_nuclei = sum(md$n_nuclei), residual_df = selection$residual_df, n_frozen_primary_genes = length(universe),
        n_endpoint_genes = NA_integer_, n_pvalue = NA_integer_, n_padj = NA_integer_, DEG_FDR05 = NA_integer_, warning_count = 0L,
        design = selection$design, execution_origin = "NOT_FITTED", cache_reason = "MODEL_DISABLED_BY_PRESPECIFIED_MASK", output_directory = folder)
      status$min_included_mouse_nuclei <- if (nrow(selection$metadata)) min(selection$metadata$n_nuclei) else NA_integer_
      status$max_included_mouse_nuclei <- if (nrow(selection$metadata)) max(selection$metadata$n_nuclei) else NA_integer_
      status$coverage_flag <- "EXCLUDED_BY_PRESPECIFIED_ATOMIC_MASK"
      ad <- new_dir(file.path(folder, "PRIMARY", "AUDIT")); write_tsv(status, file.path(ad, "EXECUTION_AUDIT.tsv")); write_tsv(selection$audit, file.path(ad, "MICE_ALL12_AUDIT.tsv"))
      status$population_key <- pop; status$tool <- row$tool; status$population <- row$population; status$analysis_kind <- row$analysis_kind; status$analysis_label <- row$analysis_label; status$branch <- branch
      statuses[[paste(pop, "PRIMARY")]] <- status
      next
    }
    endpoints <- list(list(endpoint = "PRIMARY"))
    if (branch == "original_full" && truth(row$core_total) && truth(sc$run_original_core_sensitivities %or% TRUE)) {
      endpoints <- c(endpoints, list(list(endpoint = "KNOWN_PARTITION_ADJUSTED"), list(endpoint = "MAIN_PARTITION_ONLY")),
        lapply(samples, function(sample) list(endpoint = "LOO", omit = sample)))
    }
    primary <- NULL; loo <- list()
    for (spec in endpoints) {
      label <- if (spec$endpoint == "LOO") paste0("LOO_omit_", spec$omit) else spec$endpoint
      selected <- satb2_pb_select(md, spec$endpoint, spec$omit, min_nuclei)
      pp <- satb2_pb_prepare_endpoint(y, selected, universe)
      fingerprint <- digest::digest(list(counts = pp$counts, genes = universe, membership_sha256 = row$membership_sha256,
        mice = selected$metadata[, c("sample", "condition", "sequencing_partition")],
        design = selected$design, independent_filtering = spec$endpoint != "LOO", version = as.character(utils::packageVersion("DESeq2"))), algo = "sha256")
      candidate <- if (nzchar(pp$skip)) NULL else same_run[[fingerprint]] %or% legacy_candidate(cfg, branch, row, label)
      ep <- file.path(folder, label)
      ans <- satb2_pb_run_endpoint(y, md, genes, universe, ep, endpoint = spec$endpoint,
        omit = spec$omit, independent_filtering = spec$endpoint != "LOO", min_nuclei = min_nuclei,
        gene_count_cutoff = if (inhibitory) 5L else 10L,
        gene_required_mice = if (inhibitory) 3L else min(selection$n_WT, selection$n_cKO),
        count_kind = if (branch == "paired_cellbender") "CellBender_integer_count_sensitivity_fixed_original_cluster_union" else paste0("original_CellRanger_integer_UMI_", branch),
        cache = candidate, membership_sha256 = row$membership_sha256)
      ans$audit$population_key <- pop; ans$audit$tool <- row$tool; ans$audit$population <- row$population; ans$audit$analysis_kind <- row$analysis_kind; ans$audit$analysis_label <- row$analysis_label; ans$audit$branch <- branch
      statuses[[paste(pop, label)]] <- ans$audit
      if (!is.null(ans$results)) same_run[[fingerprint]] <- revision_candidate(row$counts_file, pop_audit, ep, out, source_core,
        file.path(finalized, "POPULATION_MEMBERSHIPS.tsv.gz"), file.path(finalized, "POPULATION_INDEX.tsv"), row$membership_sha256)
      if (spec$endpoint == "PRIMARY") primary <- ans$results
      if (spec$endpoint == "LOO" && !is.null(ans$results)) loo[[spec$omit]] <- ans$results
    }
    if (!is.null(primary) && length(loo)) {
      lfc <- vapply(loo, function(z) z$log2FoldChange[match(primary$gene_id, z$gene_id)], numeric(nrow(primary)))
      if (is.null(dim(lfc))) lfc <- matrix(lfc, nrow = nrow(primary))
      finite <- is.finite(lfc)
      write_tsv(data.frame(gene_id = primary$gene_id, primary_log2FoldChange = primary$log2FoldChange, primary_padj = primary$padj,
        n_LOO_estimable = rowSums(finite), n_LOO_same_direction = rowSums(finite & sign(lfc) == sign(primary$log2FoldChange), na.rm = TRUE),
        min_LOO_log2FoldChange = apply(lfc, 1L, function(x) if (any(is.finite(x))) min(x[is.finite(x)]) else NA_real_),
        max_LOO_log2FoldChange = apply(lfc, 1L, function(x) if (any(is.finite(x))) max(x[is.finite(x)]) else NA_real_)), file.path(pop_audit, "LOO_EFFECT_INFLUENCE_ONLY.tsv"))
    }
    rm(y); gc()
  }
  audit <- do.call(rbind, statuses); rownames(audit) <- NULL
  write_tsv(audit, file.path(out, "PSEUDOBULK_EXECUTION_AUDIT.tsv")); write_tsv(index, file.path(out, "POPULATION_INDEX.tsv"))
  write_tsv(genes[, c("gene_id", "gene_symbol")], file.path(out, "AUDIT_GENE_ID_SYMBOL_MANIFEST.tsv"))
  if (any(audit$status == "ERROR")) fail("Native DESeq2 endpoint errors recorded; inspect retained audit")
  invisible(list(n_endpoints = nrow(audit), n_tested = sum(audit$status == "TESTED_REVIEW_REQUIRED"),
    n_not_tested = sum(audit$status == "NOT_TESTED"), n_reused = sum(audit$execution_origin == "REUSED_EXACT_UNCHANGED_DESEQ2_ENDPOINT"),
    n_fresh = sum(audit$execution_origin == "FRESH_NATIVE_DESEQ2_FIT")))
}
