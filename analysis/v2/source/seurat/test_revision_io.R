# Fabricated IO/aggregation fixtures only. Does not fit a model or read study H5s.
# Rscript --vanilla test_revision_io.R SOURCE_DIRECTORY
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Provide source directory")
source(file.path(args[1L], "deseq2_core.R")); source(file.path(args[1L], "revision_io.R"))
source(file.path(args[1L], "revision_models.R"))
fixture <- tempfile("satb2_revision_SYNTHETIC_IO_"); dir.create(fixture)
samples <- names(satb2_pb_roster()); genes <- data.frame(gene_id = paste0("FAKE", 1:4), gene_symbol = paste0("FakeSymbol", 1:4))
x <- Matrix::Matrix(matrix(c(10,2,0,5,4,0,3,3), 4, 2), sparse = TRUE)
x <- methods::as(x, "dgCMatrix"); bars <- c("FAKEBAR1-1", "FAKEBAR2-1")
source_rows <- list(); md <- list()
for (sample in samples) {
  p <- file.path(fixture, paste0(sample, ".h5")); f <- hdf5r::H5File$new(p, "w")
  m <- f$create_group("matrix"); fe <- m$create_group("features")
  m[["data"]] <- x@x; m[["indices"]] <- x@i; m[["indptr"]] <- x@p
  m[["shape"]] <- x@Dim; m[["barcodes"]] <- bars; fe[["id"]] <- genes$gene_id
  fe[["feature_type"]] <- rep("Gene Expression", 4L); f$create_group("droplet_latents")
  f$close_all()
  source_rows[[sample]] <- data.frame(sample = sample, condition = unname(satb2_pb_roster()[sample]),
    sequencing_partition = unname(satb2_pb_expected_partitions()[sample]), filtered_h5 = p, source_h5_sha256 = sha_file(p))
  md[[sample]] <- data.frame(nucleus_id = paste0(sample, "_", bars), barcode = bars, sample = sample,
    condition = unname(satb2_pb_roster()[sample]), nCount_RNA = as.numeric(Matrix::colSums(x)))
}
sm <- do.call(rbind, source_rows); metadata <- do.call(rbind, md); rownames(metadata) <- NULL
write_tsv(sm, file.path(fixture, "samples.tsv")); write_tsv(metadata, file.path(fixture, "roster.tsv")); write_tsv(genes, file.path(fixture, "genes.tsv"))
write_tsv(genes[1:2, ], file.path(fixture, "panel.tsv"))
ii <- list(); membership <- list()
for (tool in c("Scanpy", "Seurat")) for (cluster in c("0", "1")) {
  k <- paste(tool, "IN", cluster, sep = "_"); ids <- metadata$nucleus_id[metadata$barcode == bars[as.integer(cluster) + 1L]]
  ii[[k]] <- data.frame(population_key = k, tool = tool, population = "IN", cluster_id = cluster, n_nuclei = length(ids),
    membership_sha256 = member_sha(ids), model_enabled = TRUE)
  membership[[k]] <- data.frame(population_key = k, nucleus_id = ids)
}
write_tsv(do.call(rbind, ii), file.path(fixture, "atoms.tsv")); write_tsv(do.call(rbind, membership), file.path(fixture, "members.tsv.gz"))
cfg <- list(plot_gene_panel_tsv = file.path(fixture, "panel.tsv"), deseq2 = list(
  sample_manifest_tsv_by_branch = list(original_full = file.path(fixture, "samples.tsv")),
  roster_tsv_by_branch = list(original_full = file.path(fixture, "roster.tsv")), gene_manifest_tsv = file.path(fixture, "genes.tsv"),
  population_index_tsv_by_branch = list(original_full = file.path(fixture, "atoms.tsv")),
  population_memberships_tsv_by_branch = list(original_full = file.path(fixture, "members.tsv.gz"))))
prepared <- new_dir(file.path(fixture, "prepared")); prepare_revision(cfg, "original_full", prepared)
jsonlite::write_json(list(status = "SYNTHETIC_FIXTURE_COMPLETE"), file.path(prepared, "COMPLETE.json"))
cfg$deseq2$prepared_root_by_branch <- list(original_full = prepared)
categories <- data.frame(tool = c("Scanpy", "Seurat"), population = "IN", category_family = "FAKE_CATEGORY", label = "FakeProgram",
  member_cluster_ids = '["0","1"]', cluster_scope = "original_frozen", source_branch = "original_full")
write_tsv(categories, file.path(fixture, "categories.tsv")); cfg$deseq2$category_definition_tsv <- file.path(fixture, "categories.tsv")
final <- new_dir(file.path(fixture, "finalized")); finalize_populations(cfg, "original_full", final)
stats <- read_tsv(file.path(prepared, "gene_cluster_mouse_aggregates.tsv.gz"))
stopifnot(nrow(stats) == 4L * 12L * 2L, all(stats$n_nuclei == 1L), all(stats$cluster_scope == "original_frozen"))
a <- stats[stats$tool == "Scanpy" & stats$cluster_id == "0" & stats$gene_id == "FAKE1", ]
stopifnot(all(a$sum_counts == 10), all(a$n_detected == 1), max(abs(a$sum_log1p_normalized - log1p(10000*10/17))) < 1e-12)
fi <- read_tsv(file.path(final, "POPULATION_INDEX.tsv")); stopifnot(nrow(fi) == 8L, all(fi$model_enabled), all(fi$min_nuclei_per_mouse == 1L))
total <- fi[fi$analysis_kind == "total", ]
for (p in total$counts_file) {y <- read_integer_counts(p, genes, samples); stopifnot(all(y[, 1] == as.numeric(Matrix::rowSums(x))), all(colSums(y) == 27L))}
bad <- cfg; cc <- categories; cc$member_cluster_ids[1] <- '["0","0"]'; write_tsv(cc, file.path(fixture, "bad_categories.tsv")); bad$deseq2$category_definition_tsv <- file.path(fixture, "bad_categories.tsv")
rejection <- tryCatch({finalize_populations(bad, "original_full", new_dir(file.path(fixture, "bad_final"))); FALSE}, error = function(e) grepl("Repeated annotation cluster", conditionMessage(e)))
stopifnot(rejection)
cat("SYNTHETIC explicit CSC/extra-root import, all12 aggregation, exact library conservation, full-library panel normalization, IN1 coverage, disjoint union, duplicate rejection PASS; NO fits.\n")
cat("Fixture:", fixture, "\n")
