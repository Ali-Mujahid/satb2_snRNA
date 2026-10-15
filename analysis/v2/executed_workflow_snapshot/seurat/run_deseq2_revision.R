#!/usr/bin/env Rscript
# Native-workflow entry point. Authoring/parsing/source() is not execution.
# Rscript --vanilla run_deseq2_revision.R --config MASTER.json --stage
#   prepare|finalize_populations|deseq2 --branch BRANCH --outdir NEW_STAGE_DIR
# Existing outputs, failed evidence and previous analyses are never overwritten.

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 8L || !identical(args[c(1L, 3L, 5L, 7L)], c("--config", "--stage", "--branch", "--outdir")))
    stop("Usage: Rscript --vanilla run_deseq2_revision.R --config MASTER.json --stage STAGE --branch BRANCH --outdir NEW_STAGE_DIR", call. = FALSE)
  stage <- args[4L]; branch <- args[6L]
  if (!stage %in% c("prepare", "finalize_populations", "deseq2") || !branch %in% c("original_full", "paired_original", "paired_cellbender"))
    stop("Unknown declared stage/branch", call. = FALSE)
  needed <- c("Matrix", "hdf5r", "jsonlite", "digest")
  if (stage == "deseq2") needed <- c(needed, "DESeq2", "S4Vectors", "SummarizedExperiment")
  absent <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(absent)) stop(paste("Required installed packages absent; no installation:", paste(absent, collapse = ",")), call. = FALSE)
  file_arg <- grep("^--file=", commandArgs(), value = TRUE)
  script <- normalizePath(sub("^--file=", "", file_arg[1L]), mustWork = TRUE)
  dir <- dirname(script); core <- file.path(dir, "deseq2_core.R")
  source(core, local = TRUE); source(file.path(dir, "revision_io.R"), local = TRUE)
  source(file.path(dir, "revision_models.R"), local = TRUE)
  options(satb2.revision.statistical_core_sha256 = sha_file(core))
  cfg <- jsonlite::fromJSON(args[2L], simplifyVector = FALSE)
  options(satb2.revision.master_config_sha256 = sha_file(args[2L]))
  if (is.null(cfg$deseq2)) fail("Master must contain deseq2 contract")
  probe <- cfg$deseq2$memory_probe_python %or% cfg$memory_probe_python
  if (is.null(probe) || !file.exists(probe)) fail("Configured memory-probe Python is required")
  answer <- system2(probe, args = c("-c", shQuote("import psutil; print(psutil.virtual_memory().available/(1024**3))")), stdout = TRUE, stderr = TRUE)
  available <- suppressWarnings(as.numeric(answer)); minimum <- as.numeric(cfg$deseq2$minimum_available_gib %or% 8)
  if (length(available) != 1L || !is.finite(available) || available < minimum) fail("Available-memory guard failed before source reads/models")
  out <- new_dir(args[8L]); out <- normalizePath(out, mustWork = TRUE)
  write_tsv(data.frame(package = c("R", needed), version = c(as.character(getRversion()), vapply(needed, function(p) as.character(utils::packageVersion(p)), character(1)))), file.path(out, "PACKAGE_VERSIONS.tsv"))
  satb2_pb_write_text(capture.output(utils::sessionInfo()), file.path(out, "SESSION_INFO.txt"))
  write_tsv(data.frame(probed_utc = satb2_pb_utc(), available_gib = available, minimum_available_gib = minimum), file.path(out, "STAGE_START_AVAILABLE_MEMORY.tsv"))
  sources <- c(script, core, file.path(dir, "revision_io.R"), file.path(dir, "revision_models.R"), normalizePath(args[2L], mustWork = TRUE))
  write_tsv(data.frame(path = sources, sha256 = vapply(sources, sha_file, character(1))), file.path(out, "SOURCE_AND_MASTER_HASHES.tsv"))
  # Frozen source metadata hashes can be independently supplied by the master.
  if (!is.null(cfg$deseq2$input_hashes_tsv)) {
    h <- read_tsv(cfg$deseq2$input_hashes_tsv); required_columns(h, c("path", "sha256"), "frozen inputs")
    for (i in seq_len(nrow(h))) hash_guard(h$path[i], h$sha256[i], "frozen master input")
    write_tsv(h, file.path(out, "FROZEN_INPUT_HASH_GUARDS.tsv"))
  }
  jsonlite::write_json(list(status = "EXECUTION_IN_PROGRESS", stage = stage, branch = branch,
    started_utc = satb2_pb_utc(), config_sha256 = sha_file(args[2L]), method = "DESeq2_only", scientific_acceptance = "NOT_ASSESSED"), file.path(out, "START.json"), pretty = TRUE, auto_unbox = TRUE)
  tryCatch({
    result <- switch(stage,
      prepare = prepare_revision(cfg, branch, out),
      finalize_populations = finalize_populations(cfg, branch, out),
      deseq2 = run_revision_models(cfg, branch, out, core))
    published <- list.files(out, recursive = TRUE, all.files = TRUE, no.. = TRUE, full.names = FALSE)
    published <- published[!grepl("\\.partial$", published)]
    output_hashes <- setNames(as.list(vapply(file.path(out, published), sha_file, character(1))), published)
    jsonlite::write_json(c(list(status = "EXECUTED_REVIEW_REQUIRED", stage = stage, branch = branch,
      completed_utc = satb2_pb_utc(), method = "DESeq2_only", scientific_acceptance = "NOT_ASSESSED",
      scientific_fitting_performed = stage == "deseq2" && (result$n_fresh %or% 0L) > 0L,
      raw_nucleus_count_aggregation_performed = stage == "prepare",
      output_sha256 = output_hashes), result),
      file.path(out, "COMPLETE.json"), pretty = TRUE, auto_unbox = TRUE)
  }, error = function(e) {
    jsonlite::write_json(list(status = "FAILED_REVIEW_REQUIRED", stage = stage, branch = branch,
      error = conditionMessage(e), failed_utc = satb2_pb_utc()), file.path(out, "FAILED.json"), pretty = TRUE, auto_unbox = TRUE)
    fail(conditionMessage(e))
  })
}
if (sys.nframe() == 0L) main()
