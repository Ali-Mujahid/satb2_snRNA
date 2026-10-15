# Bounded count/panel preparation and audited atomic-population sums.
# Functions only. Scientific execution belongs to the registered native workflow.
`%or%` <- function(x, y) if (is.null(x)) y else x
fail <- satb2_pb_fail
truth <- function(x) !is.na(x) & tolower(as.character(x)) %in% c("true", "1", "yes")
read_tsv <- function(p) utils::read.delim(p, check.names = FALSE, stringsAsFactors = FALSE)
check_unique <- function(x, what) {
  if (anyNA(x) || any(!nzchar(x)) || anyDuplicated(x) || any(grepl("[\t\r\n]", x))) fail("Invalid IDs: ", what)
}
new_dir <- function(p) {
  if (dir.exists(p) && length(list.files(p, all.files = TRUE, no.. = TRUE))) fail("Existing output content: ", p)
  if (!dir.exists(p) && !dir.create(p, recursive = TRUE)) fail("Cannot create: ", p)
  p
}
write_tsv <- function(x, p) {
  if (file.exists(p) || file.exists(paste0(p, ".partial"))) fail("Existing output/partial: ", p)
  tmp <- paste0(p, ".partial")
  if (grepl("\\.gz$", p)) {
    con <- gzfile(tmp, "wt", compression = 6L)
    tryCatch(utils::write.table(x, con, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"), finally = close(con))
  } else utils::write.table(x, tmp, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
  if (!file.rename(tmp, p)) fail("Cannot publish: ", p)
}
sha_file <- function(p) digest::digest(file = p, algo = "sha256")
member_sha <- function(ids) digest::digest(if (length(ids)) paste0(paste(sort(ids, method = "radix"), collapse = "\n"), "\n") else "", algo = "sha256", serialize = FALSE)
required_columns <- function(tab, wanted, name) if (!all(wanted %in% names(tab))) fail("Missing ", name, " fields: ", paste(setdiff(wanted, names(tab)), collapse = ","))
hash_guard <- function(path, expected, role) {
  if (length(path) != 1L || is.na(path) || !file.exists(path) || length(expected) != 1L || is.na(expected) ||
      !grepl("^[0-9a-fA-F]{64}$", expected) || sha_file(path) != tolower(expected)) fail("SHA256/input guard failed: ", role)
}

load_count_contract <- function(cfg, branch) {
  sc <- cfg$deseq2
  sm <- read_tsv(sc$sample_manifest_tsv_by_branch[[branch]])
  md <- read_tsv(sc$roster_tsv_by_branch[[branch]])
  genes <- read_tsv(sc$gene_manifest_tsv)
  required_columns(sm, c("sample", "condition", "sequencing_partition", "filtered_h5", "source_h5_sha256"), "capture manifest")
  required_columns(md, c("nucleus_id", "barcode", "sample", "condition", "nCount_RNA"), "frozen count roster")
  required_columns(genes, c("gene_id", "gene_symbol"), "genes")
  check_unique(sm$sample, "captures"); check_unique(md$nucleus_id, "nuclei"); check_unique(genes$gene_id, "genes")
  roster <- satb2_pb_roster(); samples <- names(roster)
  if (!setequal(sm$sample, samples) || !setequal(unique(md$sample), samples) ||
      anyNA(md$condition) || any(md$condition != unname(roster[md$sample]))) fail("All12 independent-mouse roster mismatch")
  sm <- sm[match(samples, sm$sample), , drop = FALSE]
  if (any(sm$condition != unname(roster)) || any(sm$sequencing_partition != unname(satb2_pb_expected_partitions()))) fail("Reviewed sample/partition assignment changed")
  if (any(md$nucleus_id != paste0(md$sample, "_", md$barcode))) fail("sample_barcode identity mismatch")
  md$sequencing_partition <- sm$sequencing_partition[match(md$sample, sm$sample)]
  if (branch == "paired_cellbender") required_columns(md, "cellbender_nCount_RNA", "corrected libraries")
  genes$gene_id <- as.character(genes$gene_id); genes$gene_symbol <- as.character(genes$gene_symbol)
  genes$symbol_missing <- is.na(genes$gene_symbol) | !nzchar(trimws(genes$gene_symbol))
  genes$display_symbol <- ifelse(genes$symbol_missing, NA_character_, genes$gene_symbol)
  list(samples = sm, metadata = md, genes = genes, count_source = branch)
}

# Explicit /matrix CSC importer avoids CellBender's nonmatrix latent root groups.
read_10x_matrix_csc <- function(path) {
  f <- hdf5r::H5File$new(path, mode = "r"); on.exit(f$close_all(), add = TRUE)
  required <- c("matrix/data", "matrix/indices", "matrix/indptr", "matrix/shape", "matrix/barcodes", "matrix/features/id", "matrix/features/feature_type")
  if (!all(vapply(required, f$exists, logical(1)))) fail("Explicit10x CSC schema absent: ", path)
  shape <- as.numeric(f[["matrix/shape"]]$read()); ids <- as.character(f[["matrix/features/id"]]$read())
  bars <- as.character(f[["matrix/barcodes"]]$read()); kinds <- as.character(f[["matrix/features/feature_type"]]$read())
  check_unique(ids, "H5 genes"); check_unique(bars, "H5 barcodes")
  if (length(shape) != 2L || anyNA(shape) || any(!is.finite(shape)) || any(shape <= 0) || any(shape != floor(shape)) ||
      any(shape > .Machine$integer.max) || length(ids) != shape[1L] || length(bars) != shape[2L] ||
      length(kinds) != shape[1L] || anyNA(kinds) || any(kinds != "Gene Expression")) fail("H5 axis/GEX mismatch")
  xx <- as.numeric(f[["matrix/data"]]$read()); ii <- as.numeric(f[["matrix/indices"]]$read()); pp <- as.numeric(f[["matrix/indptr"]]$read())
  if (length(ii) != length(xx) || length(pp) != shape[2L] + 1L || anyNA(xx) || any(!is.finite(xx)) || any(xx < 0) || any(xx != floor(xx)) ||
      anyNA(ii) || any(ii < 0 | ii >= shape[1L]) || any(ii != floor(ii)) || anyNA(pp) || pp[1L] != 0 ||
      tail(pp, 1L) != length(xx) || any(diff(pp) < 0) || any(pp != floor(pp)) || any(pp > .Machine$integer.max)) fail("Invalid integer CSC arrays")
  x <- methods::new("dgCMatrix", i = as.integer(ii), p = as.integer(pp), x = xx,
    Dim = as.integer(shape), Dimnames = list(ids, bars))
  methods::validObject(x); x
}
read_capture <- function(sample, contract) {
  row <- contract$samples[match(sample, contract$samples$sample), , drop = FALSE]
  md <- contract$metadata[contract$metadata$sample == sample, , drop = FALSE]
  hash_guard(row$filtered_h5, row$source_h5_sha256, paste0("capture/", sample))
  x <- read_10x_matrix_csc(row$filtered_h5); j <- match(md$barcode, colnames(x))
  if (anyNA(j) || !setequal(rownames(x), contract$genes$gene_id)) fail("Count gene/barcode identity mismatch: ", sample)
  x <- methods::as(x[contract$genes$gene_id, j, drop = FALSE], "dgCMatrix"); colnames(x) <- md$nucleus_id
  lib <- as.numeric(Matrix::colSums(x))
  expected <- if (contract$count_source == "paired_cellbender") md$cellbender_nCount_RNA else md$nCount_RNA
  if (any(lib <= 0) || anyNA(expected) || !identical(lib, as.numeric(expected))) fail("Chosen-source full-library audit mismatch: ", sample)
  list(counts = x, metadata = md, library = lib)
}
read_integer_counts <- function(path, genes, samples) {
  z <- read_tsv(path); required_columns(z, c("gene_id", samples), "gene×12 counts")
  check_unique(as.character(z$gene_id), "count stable genes")
  if (!identical(as.character(z$gene_id), genes$gene_id)) fail("Stable gene order mismatch: ", path)
  y <- as.matrix(z[, samples, drop = FALSE]); storage.mode(y) <- "double"
  if (anyNA(y) || any(!is.finite(y)) || any(y < 0) || any(y != floor(y)) || any(y > .Machine$integer.max)) fail("Noninteger/overflow aggregate: ", path)
  rownames(y) <- genes$gene_id; storage.mode(y) <- "integer"; y
}
load_atomic <- function(cfg, branch, contract) {
  sc <- cfg$deseq2
  index <- read_tsv(sc$population_index_tsv_by_branch[[branch]])
  mm <- read_tsv(sc$population_memberships_tsv_by_branch[[branch]])
  if ("branch" %in% names(index)) index <- index[index$branch == branch, , drop = FALSE]
  if ("branch" %in% names(mm)) mm <- mm[mm$branch == branch, , drop = FALSE]
  required_columns(index, c("population_key", "tool", "population", "cluster_id", "n_nuclei", "membership_sha256"), "atomic registry")
  required_columns(mm, c("population_key", "nucleus_id"), "atomic memberships")
  check_unique(index$population_key, "atomic keys")
  if (any(!grepl("^[A-Za-z0-9][A-Za-z0-9_-]*$", index$population_key))) fail("Unsafe atomic directory key")
  if (any(!mm$population_key %in% index$population_key) || any(!mm$nucleus_id %in% contract$metadata$nucleus_id) ||
      anyDuplicated(paste(mm$population_key, mm$nucleus_id, sep = "\r"))) fail("Atomic memberships changed/duplicated")
  members <- setNames(lapply(index$population_key, function(k) sort(as.character(mm$nucleus_id[mm$population_key == k]), method = "radix")), index$population_key)
  if (any(vapply(members, length, integer(1)) != as.integer(index$n_nuclei)) ||
      any(vapply(members, member_sha, character(1)) != tolower(index$membership_sha256))) fail("Atomic size/membership SHA mismatch")
  # Atomic clusters within one tool/population must form a disjoint partition.
  for (f in unique(paste(index$tool, index$population, sep = "\r"))) {
    ids <- unlist(members[paste(index$tool, index$population, sep = "\r") == f], use.names = FALSE)
    if (anyDuplicated(ids)) fail("Overlapping atomic clusters inside tool/population")
  }
  list(index = index, members = members, memberships = mm)
}

prepare_revision <- function(cfg, branch, out) {
  contract <- load_count_contract(cfg, branch); aa <- load_atomic(cfg, branch, contract)
  index <- aa$index; members <- aa$members; genes <- contract$genes; mouse <- contract$samples[, c("sample", "condition", "sequencing_partition")]
  samples <- mouse$sample; panel <- read_tsv(cfg$plot_gene_panel_tsv %or% cfg$deseq2$plot_gene_panel_tsv)
  required_columns(panel, c("gene_id", "gene_symbol"), "plot panel"); check_unique(panel$gene_id, "panel genes")
  if (any(!panel$gene_id %in% genes$gene_id)) fail("Plot panel includes unsupported stable genes")
  panel <- genes[match(panel$gene_id, genes$gene_id), c("gene_id", "gene_symbol"), drop = FALSE]
  n <- length(members); aggregates <- vector("list", n); names(aggregates) <- names(members); cached <- rep(FALSE, n)
  for (p in seq_len(n)) {
    source <- if ("source_counts_file" %in% names(index)) index$source_counts_file[p] else NA_character_
    if (!is.na(source) && nzchar(source)) {
      if (!"source_counts_sha256" %in% names(index)) fail("Cached count SHA required")
      exact <- all(c("source_membership_sha256", "source_count_branch") %in% names(index)) &&
        !is.na(index$source_membership_sha256[p]) && !is.na(index$source_count_branch[p]) &&
        index$source_membership_sha256[p] == index$membership_sha256[p] && index$source_count_branch[p] == branch
      if (exact) {
        hash_guard(source, index$source_counts_sha256[p], paste0("cached atomic counts/", index$population_key[p]))
        aggregates[[p]] <- read_integer_counts(source, genes, samples); cached[p] <- TRUE
      } else aggregates[[p]] <- matrix(0L, nrow(genes), length(samples), dimnames = list(genes$gene_id, samples))
    } else aggregates[[p]] <- matrix(0L, nrow(genes), length(samples), dimnames = list(genes$gene_id, samples))
  }
  coverage <- list(); stats <- list(); conservation <- list()
  for (sample in samples) {
    cap <- read_capture(sample, contract); ids <- colnames(cap$counts)
    ii <- lapply(members, function(x) match(intersect(x, ids), ids)); nn <- vapply(ii, length, integer(1))
    assignment <- Matrix::sparseMatrix(i = unlist(ii, use.names = FALSE), j = rep.int(seq_len(n), nn), x = 1,
      dims = c(length(ids), n))
    expected <- as.numeric(crossprod(assignment, cap$library))
    fresh <- which(!cached)
    if (length(fresh)) {
      y <- as.matrix(cap$counts %*% assignment[, fresh, drop = FALSE])
      for (j in seq_along(fresh)) aggregates[[fresh[j]]][, sample] <- y[, j]
    }
    actual <- vapply(aggregates, function(y) sum(y[, sample]), numeric(1))
    if (!identical(as.numeric(actual), expected)) fail("Atomic full-library conservation mismatch: ", sample)
    px <- Matrix::drop0(cap$counts[panel$gene_id, , drop = FALSE])
    sums <- as.matrix(px %*% assignment); detected <- px; detected@x[] <- 1
    det <- as.matrix(detected %*% assignment)
    norm <- px %*% Matrix::Diagonal(x = 10000 / cap$library); norm <- methods::as(norm, "dgCMatrix"); norm@x <- log1p(norm@x)
    logs <- as.matrix(norm %*% assignment)
    for (p in seq_len(n)) {
      keys <- data.frame(tool = index$tool[p], branch = branch, population = index$population[p],
        cluster_scope = "original_frozen", cluster_id = as.character(index$cluster_id[p]), sample = sample,
        condition = unname(satb2_pb_roster()[sample]), n_nuclei = nn[p], stringsAsFactors = FALSE)
      if ("model_enabled" %in% names(index)) keys$model_enabled <- truth(index$model_enabled[p])
      coverage[[paste(sample, p)]] <- keys
      z <- keys[rep(1L, nrow(panel)), , drop = FALSE]
      z$gene_id <- panel$gene_id; z$gene_symbol <- panel$gene_symbol
      z$sum_counts <- sums[, p]; z$n_detected <- det[, p]; z$sum_log1p_normalized <- logs[, p]
      stats[[paste(sample, p)]] <- z
    }
    conservation[[sample]] <- data.frame(sample = sample, branch = branch, source_nuclei = length(ids), source_UMIs = sum(cap$library),
      atomic_membership_UMIs = sum(actual), exact_population_library_conservation = TRUE,
      note = "tool/population partitions may overlap across tools or nested biological populations")
    rm(cap, px, norm, detected, assignment, sums, det, logs); gc()
  }
  for (p in seq_len(n)) {
    folder <- new_dir(file.path(out, index$population_key[p])); audit <- new_dir(file.path(folder, "AUDIT"))
    y <- aggregates[[p]]
    if (anyNA(y) || any(!is.finite(y)) || any(y < 0) || any(y != floor(y)) || any(y > .Machine$integer.max)) fail("Invalid final atomic counts")
    storage.mode(y) <- "integer"
    c <- mouse; c$n_nuclei <- vapply(samples, function(s) sum(contract$metadata$sample[match(members[[p]], contract$metadata$nucleus_id)] == s), integer(1))
    c$total_UMIs <- as.numeric(colSums(y))
    if (any((c$n_nuclei == 0L) != (c$total_UMIs == 0))) fail("Absent/zero atomic population discrepancy")
    counts_file <- file.path(audit, "SOURCE_INTEGER_COUNTS_ALL_GENES_ALL12.tsv.gz")
    write_tsv(data.frame(gene_id = genes$gene_id, gene_symbol = genes$gene_symbol, y, check.names = FALSE), counts_file)
    write_tsv(c, file.path(audit, "ALL12_POPULATION_COVERAGE.tsv"))
    index$counts_file[p] <- counts_file; index$counts_sha256[p] <- sha_file(counts_file)
    index$count_execution_origin[p] <- if (cached[p]) "REUSED_EXACT_BOUND_ATOMIC_COUNTS" else "FRESH_NATIVE_ATOMIC_AGGREGATE"
  }
  write_tsv(index, file.path(out, "POPULATION_INDEX.tsv")); write_tsv(aa$memberships, file.path(out, "POPULATION_MEMBERSHIPS.tsv.gz"))
  write_tsv(genes[, c("gene_id", "gene_symbol")], file.path(out, "GENE_ID_SYMBOL_MANIFEST.tsv"))
  write_tsv(do.call(rbind, coverage), file.path(out, "cluster_mouse_counts.tsv"))
  write_tsv(do.call(rbind, stats), file.path(out, "gene_cluster_mouse_aggregates.tsv.gz"))
  write_tsv(do.call(rbind, conservation), file.path(out, "STREAMED_COUNT_CONSERVATION.tsv"))
  invisible(list(n_atomic_populations = n, n_nuclei = nrow(contract$metadata), n_panel_genes = nrow(panel)))
}

finalize_populations <- function(cfg, branch, out) {
  sc <- cfg$deseq2; prepared <- sc$prepared_root_by_branch[[branch]]
  if (!file.exists(file.path(prepared, "COMPLETE.json"))) fail("Atomic producer incomplete")
  index <- read_tsv(file.path(prepared, "POPULATION_INDEX.tsv"))
  mm <- read_tsv(file.path(prepared, "POPULATION_MEMBERSHIPS.tsv.gz"))
  genes <- read_tsv(file.path(prepared, "GENE_ID_SYMBOL_MANIFEST.tsv")); genes$gene_id <- as.character(genes$gene_id)
  required_columns(index, c("population_key", "tool", "population", "cluster_id", "n_nuclei", "membership_sha256", "counts_file", "counts_sha256"), "prepared atoms")
  check_unique(index$population_key, "prepared atoms")
  samples <- names(satb2_pb_roster()); arrays <- setNames(vector("list", nrow(index)), index$population_key)
  atomic_coverage <- setNames(vector("list", nrow(index)), index$population_key)
  for (p in seq_len(nrow(index))) {
    hash_guard(index$counts_file[p], index$counts_sha256[p], "prepared atomic array")
    arrays[[p]] <- read_integer_counts(index$counts_file[p], genes, samples)
    atomic_coverage[[p]] <- read_tsv(file.path(dirname(index$counts_file[p]), "ALL12_POPULATION_COVERAGE.tsv"))
    ids <- as.character(mm$nucleus_id[mm$population_key == index$population_key[p]])
    if (length(ids) != index$n_nuclei[p] || member_sha(ids) != index$membership_sha256[p]) fail("Prepared atom roster changed")
  }
  definitions <- list(); atomic_sets <- list(); registry <- list()
  add <- function(tool, pop, kind, label, family, atoms, enabled = TRUE, cache_key = NA_character_) {
    ordinal <- length(definitions) + 1L
    key <- sprintf("P%04d_%s_%s_%s", ordinal, gsub("[^A-Za-z0-9_-]", "_", tool), gsub("[^A-Za-z0-9_-]", "_", pop), toupper(kind))
    if (anyDuplicated(atoms) || any(!atoms %in% index$population_key)) fail("Population union contains duplicate/unknown atoms")
    if (length(atoms) && any(index$tool[match(atoms, index$population_key)] != tool | index$population[match(atoms, index$population_key)] != pop))
      fail("Population union crossed tools/biological source populations")
    ids <- as.character(mm$nucleus_id[mm$population_key %in% atoms])
    if (anyDuplicated(ids)) fail("Population atomic union is not disjoint")
    definitions[[key]] <<- data.frame(population_key = key, tool = tool, population = pop,
      analysis_kind = kind, analysis_label = as.character(label), category_family = family,
      model_enabled = enabled, core_total = kind == "total", n_nuclei = length(ids),
      membership_sha256 = member_sha(ids), source_label_branch = "original_full", branch = branch,
      min_nuclei_per_mouse = if (toupper(pop) %in% c("IN", "INHIBITORY")) 1L else 20L,
      cache_population_key = cache_key, stringsAsFactors = FALSE)
    atomic_sets[[key]] <<- atoms
    registry[[key]] <<- data.frame(population_key = rep(key, length(ids)), nucleus_id = sort(ids, method = "radix"), stringsAsFactors = FALSE)
    key
  }
  # Masks are cohort-specific: original sizes for original_full, common paired
  # sizes for both paired branches. Small CAI/IN atoms are retained for modeling
  # subject to independent-mouse eligibility; no overall100 cut applies there.
  groups <- unique(index[, c("tool", "population"), drop = FALSE])
  keep <- if ("model_enabled" %in% names(index)) truth(index$model_enabled) else
    !(toupper(index$population) %in% c("ALL", "ALLCELLS", "EX", "EXCITATORY") & index$n_nuclei < 100L)
  expected <- index$n_nuclei > 0L & !(toupper(index$population) %in% c("ALL", "ALLCELLS", "EX", "EXCITATORY") & index$n_nuclei < 100L)
  if (!identical(keep, expected)) fail("Atomic mask violates declared All/EX100 and CAI/IN no100 policy")
  # No iterative filtering is redone here: initial annotation owns the reviewed
  # shared intersection. Verify the two tools' retained total rosters agree.
  for (pop in unique(index$population)) {
    tools <- unique(index$tool[index$population == pop]); rosters <- lapply(tools, function(tool)
      sort(as.character(mm$nucleus_id[mm$population_key %in% index$population_key[index$tool == tool & index$population == pop & keep]]), method = "radix"))
    if (length(rosters) > 1L && any(!vapply(rosters[-1L], identical, logical(1), rosters[[1L]]))) fail("Tool total rosters differ after shared mask: ", pop)
  }
  for (g in seq_len(nrow(groups))) {
    j <- which(index$tool == groups$tool[g] & index$population == groups$population[g])
    retained <- index$population_key[j[keep[j]]]
    add(groups$tool[g], groups$population[g], "total", "TOTAL_FILTERED", "TOTAL", retained)
    for (p in j) {
      cache_key <- if ("cache_population_key" %in% names(index)) index$cache_population_key[p] else NA_character_
      add(index$tool[p], index$population[p], "cluster", index$cluster_id[p], "CLUSTER", index$population_key[p], keep[p], cache_key)
    }
  }
  definition_path <- sc$category_definition_tsv %or% cfg$category_definition_tsv
  if (!is.null(definition_path)) {
    cats <- read_tsv(definition_path)
    required_columns(cats, c("tool", "population", "category_family", "label", "member_cluster_ids", "cluster_scope", "source_branch"), "annotation categories")
    if ("branch" %in% names(cats)) cats <- cats[cats$branch == branch, , drop = FALSE]
    if (any(cats$cluster_scope != "original_frozen") || any(cats$source_branch != "original_full")) fail("Only original frozen cluster unions enter inference")
    seen <- list()
    for (i in seq_len(nrow(cats))) {
      cl <- as.character(jsonlite::fromJSON(cats$member_cluster_ids[i]))
      if (anyDuplicated(cl)) fail("Repeated annotation cluster")
      eligible <- which(index$tool == cats$tool[i] & index$population == cats$population[i])
      match_j <- eligible[match(cl, as.character(index$cluster_id[eligible]))]
      if (anyNA(match_j)) fail("Annotation names cluster outside its tool/population registry")
      family <- paste(cats$tool[i], cats$population[i], cats$category_family[i], sep = "\r")
      if (any(index$population_key[match_j] %in% seen[[family]])) fail("Category family assigns one atom to multiple labels")
      seen[[family]] <- c(seen[[family]], index$population_key[match_j])
      atoms <- index$population_key[match_j[keep[match_j]]]
      kind <- if ("analysis_kind" %in% names(cats)) as.character(cats$analysis_kind[i]) else
        if (cats$category_family[i] == "CA1_LAYER_EXPRESSION_PROGRAM") "layer_program" else
          if (cats$category_family[i] == "IN_EXPRESSION_PROGRAM") "expression_program" else "category"
      add(cats$tool[i], cats$population[i], kind, cats$label[i], cats$category_family[i], atoms)
    }
    membership_path <- sc$category_membership_tsv %or% cfg$category_membership_tsv
    if (!is.null(membership_path)) {
      cm <- read_tsv(membership_path)
      required_columns(cm, c("tool", "branch", "population", "category_family", "label", "nucleus_id"), "category memberships")
      cm <- cm[cm$branch == branch, , drop = FALSE]
      for (key in names(definitions)) {
        d <- definitions[[key]]
        if (d$analysis_kind %in% c("cluster", "total")) next
        given <- as.character(cm$nucleus_id[cm$tool == d$tool & cm$population == d$population &
          cm$category_family == d$category_family & cm$label == d$analysis_label])
        if (anyDuplicated(given) || !identical(sort(given, method = "radix"), registry[[key]]$nucleus_id))
          fail("Reviewed category membership does not equal exact disjoint atomic union")
      }
      write_tsv(data.frame(path = membership_path, sha256 = sha_file(membership_path)), file.path(out, "CATEGORY_MEMBERSHIP_SOURCE_HASH.tsv"))
    }
    write_tsv(cats, file.path(out, "SOURCE_CATEGORY_DEFINITION.tsv"))
  }
  pi <- do.call(rbind, definitions); rownames(pi) <- NULL
  manifest <- list(); counts_hash <- list()
  for (key in names(definitions)) {
    atoms <- atomic_sets[[key]]; y <- matrix(0, nrow(genes), length(samples), dimnames = list(genes$gene_id, samples))
    c <- data.frame(sample = samples, condition = unname(satb2_pb_roster()), sequencing_partition = unname(satb2_pb_expected_partitions()), n_nuclei = 0L, total_UMIs = 0)
    for (atom in atoms) {
      y <- y + arrays[[atom]]
      cv <- atomic_coverage[[atom]]; cv <- cv[match(samples, cv$sample), , drop = FALSE]
      if (anyNA(cv$n_nuclei)) fail("Atomic all12 coverage changed")
      c$n_nuclei <- c$n_nuclei + as.integer(cv$n_nuclei); c$total_UMIs <- c$total_UMIs + as.numeric(cv$total_UMIs)
    }
    if (any(y > .Machine$integer.max) || !identical(as.numeric(colSums(y)), as.numeric(c$total_UMIs)) ||
        sum(c$n_nuclei) != pi$n_nuclei[match(key, pi$population_key)]) fail("Atomic union count/nucleus conservation failed")
    storage.mode(y) <- "integer"
    folder <- new_dir(file.path(out, key)); audit <- new_dir(file.path(folder, "AUDIT"))
    p <- file.path(audit, "SOURCE_INTEGER_COUNTS_ALL_GENES_ALL12.tsv.gz")
    write_tsv(data.frame(gene_id = genes$gene_id, gene_symbol = genes$gene_symbol, y, check.names = FALSE), p)
    write_tsv(c, file.path(audit, "ALL12_POPULATION_COVERAGE.tsv"))
    pi$counts_file[match(key, pi$population_key)] <- p; pi$counts_sha256[match(key, pi$population_key)] <- sha_file(p)
    manifest[[key]] <- data.frame(population_key = rep(key, length(atoms)), atomic_population_key = atoms)
  }
  write_tsv(pi, file.path(out, "POPULATION_INDEX.tsv")); write_tsv(do.call(rbind, registry), file.path(out, "POPULATION_MEMBERSHIPS.tsv.gz"))
  write_tsv(do.call(rbind, manifest), file.path(out, "POPULATION_ATOMS.tsv")); write_tsv(genes, file.path(out, "GENE_ID_SYMBOL_MANIFEST.tsv"))
  write_tsv(index[, intersect(c("population_key", "tool", "population", "cluster_id", "n_nuclei", "membership_sha256", "model_enabled"), names(index)), drop = FALSE], file.path(out, "ATOMIC_MASK_REGISTRY.tsv"))
  invisible(list(n_populations = nrow(pi), n_enabled = sum(pi$model_enabled), n_atomic = nrow(index)))
}
