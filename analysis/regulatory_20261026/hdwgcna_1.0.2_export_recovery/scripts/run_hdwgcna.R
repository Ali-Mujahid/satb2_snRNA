#!/usr/bin/env Rscript
# Drafted against hdWGCNA 0.4.12, e3344d1f7bbac4264adf94f7aa31e0802fa8282d.
# Scientific execution belongs to the approval-bound Snakemake workflow.
# This file never installs packages, changes original analyses, or runs DESeq2.

`%||%` <- function(x, y) if (is.null(x)) y else x

arguments <- function(x) {
  ans <- list(); i <- 1L
  while (i <= length(x)) {
    if (!grepl('^--', x[i]) || i == length(x)) stop('Expected --key value arguments')
    ans[[sub('^--', '', x[i])]] <- x[i + 1L]; i <- i + 2L
  }
  for (key in c('config', 'input-dir', 'outdir')) if (is.null(ans[[key]])) stop('Missing --', key)
  ans
}

defaults <- function() list(
  execute = FALSE, seed = 2601026L, threads = 4L,
  expression_fraction = 0.05, powers = c(1:10, seq(12, 30, 2)),
  topology_fit_min = 0.8, mean_connectivity_min = 5,
  min_metacells = 30L, min_module_size = 50L, merge_cut_height = 0.2,
  deep_split = 4L, max_block_size = 3000L, max_p_outliers = 0.05,
  neighbor_hvgs = 2000L, neighbor_pcs = 30L, max_iter = 5000L,
  n_hubs = 20L, n_hub_plot_modules = 12L, bootstrap_iterations = 200L,
  consensus = TRUE, consensus_min_metacells_per_condition = 30L,
  consensus_min_mice_per_condition = 4L, consensus_min_metacells_per_mouse = 2L,
  same_partition_values = c('23HJT3LT4_lanes6_7'), same_lane_values = c('6', '7', 'L006', 'L007'),
  exploratory_power = NULL,
  CA1 = list(k = 25L, max_shared = 10L, min_cells = 50L, target_metacells = 50L),
  IN = list(k = 10L, max_shared = 3L, min_cells = 20L, target_metacells = 25L)
)

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, '.writing'); con <- if (grepl('[.]gz$', path)) gzfile(tmp, 'wt') else file(tmp, 'wt')
  tryCatch(write.table(x, con, sep = '\t', quote = FALSE, row.names = FALSE, na = 'NA'), finally = close(con))
  if (!file.rename(tmp, path)) stop('Cannot atomically write ', path)
}
write_json <- function(x, path) {
  tmp <- paste0(path, '.writing'); jsonlite::write_json(x, tmp, auto_unbox = TRUE, pretty = TRUE, null = 'null', na = 'null')
  if (!file.rename(tmp, path)) stop('Cannot write ', path)
}
save_rds <- function(x, path) {
  tmp <- paste0(path, '.writing'); saveRDS(x, tmp, compress = 'gzip')
  if (!file.rename(tmp, path)) stop('Cannot write ', path)
}
md5 <- function(paths) setNames(as.character(tools::md5sum(paths)), basename(paths))
fingerprint <- function(x) digest::digest(x, algo = 'sha256', serialize = TRUE)

ensure_identity <- function(outdir, expected) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(outdir, 'run_identity.json')
  if (file.exists(path)) {
    old <- jsonlite::read_json(path, simplifyVector = TRUE)
    if (!identical(old$fingerprint, expected)) stop('Output directory belongs to different inputs/settings/code; choose a NEW run: ', outdir)
  } else {
    if (length(list.files(outdir, all.files = TRUE, no.. = TRUE))) stop('Nonempty output lacks a verifiable run identity; choose a NEW run: ', outdir)
    write_json(list(fingerprint = expected, created_utc = format(Sys.time(), tz = 'UTC', usetz = TRUE)), path)
  }
}

checkpoint_valid <- function(outdir, expected) {
  p <- file.path(outdir, 'completed.json'); if (!file.exists(p)) return(FALSE)
  ck <- jsonlite::read_json(p, simplifyVector = TRUE)
  if (!identical(ck$fingerprint, expected)) stop('Completed output has different inputs/settings; choose a NEW output directory: ', outdir)
  f <- file.path(outdir, names(ck$files))
  if (any(!file.exists(f))) stop('Completed checkpoint has missing files; preserve it and choose a NEW run: ', outdir)
  if (!identical(unname(as.character(tools::md5sum(f))), unname(as.character(ck$files)))) stop('Completed checkpoint checksum failed; preserve it and choose a NEW run: ', outdir)
  TRUE
}
complete <- function(outdir, expected, details) {
  f <- list.files(outdir, recursive = TRUE, full.names = TRUE)
  f <- f[!grepl('completed[.]json$|[.]writing$|[.]log$', f)]
  hashes <- as.character(tools::md5sum(f)); names(hashes) <- substring(f, nchar(outdir) + 2L)
  write_json(c(list(fingerprint = expected, completed_utc = format(Sys.time(), tz = 'UTC', usetz = TRUE), files = as.list(hashes)), details), file.path(outdir, 'completed.json'))
}

save_plot <- function(p, stem, width = 10, height = 7) {
  dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(paste0(stem, '.png'), p, width = width, height = height, dpi = 180, bg = 'white', limitsize = FALSE)
  ggplot2::ggsave(paste0(stem, '.pdf'), p, width = width, height = height, device = grDevices::cairo_pdf, bg = 'white', limitsize = FALSE)
  ggplot2::ggsave(paste0(stem, '.svg'), p, width = width, height = height, device = grDevices::svg, bg = 'white', limitsize = FALSE)
}

read_input <- function(input_dir, lineage) {
  d <- file.path(input_dir, lineage)
  files <- file.path(d, c('counts.mtx.gz', 'genes.tsv', 'cells.tsv'))
  if (any(!file.exists(files))) stop('Missing lineage inputs: ', paste(files[!file.exists(files)], collapse = ', '))
  genes <- read.delim(files[2], check.names = FALSE, stringsAsFactors = FALSE)
  cells <- read.delim(files[3], check.names = FALSE, stringsAsFactors = FALSE)
  required <- c('barcode', 'nucleus_id', 'sample_id', 'condition', 'endpoint', 'library_umi_full_genes')
  if (!all(required %in% names(cells))) stop('Required cells.tsv columns absent: ', paste(setdiff(required, names(cells)), collapse = ','))
  if (!all(c('feature_id', 'gene_symbol') %in% names(genes))) stop('Required gene identifiers absent')
  con <- gzfile(files[1], 'rt'); counts <- tryCatch(as(Matrix::readMM(con), 'CsparseMatrix'), finally = close(con))
  if (!identical(dim(counts), c(nrow(genes), nrow(cells)))) stop('MTX dimensions disagree with row-ordered metadata')
  if (anyDuplicated(genes$feature_id) || anyDuplicated(cells$nucleus_id)) stop('Nonunique feature/nucleus key')
  if (any(counts@x < 0) || any(abs(counts@x - round(counts@x)) > 1e-8)) stop('Expected nonnegative integer raw UMIs')
  expected <- if (lineage == 'CA1') c('CA1_Deep', 'CA1_Superficial') else c('IN_C1', 'IN_C6', 'IN_C8')
  if (!setequal(unique(cells$endpoint), expected)) stop('Unexpected exact endpoint memberships')
  cells$condition <- ifelse(tolower(cells$condition) == 'wt', 'WT', ifelse(tolower(cells$condition) %in% c('cko', 'ko'), 'cKO', NA_character_))
  if (anyNA(cells$condition)) stop('Unknown condition encoding')
  if (any(!is.finite(cells$library_umi_full_genes) | cells$library_umi_full_genes <= 0)) stop('Invalid original full-library denominator')
  if (any(Matrix::colSums(counts) > cells$library_umi_full_genes + 1e-6)) stop('Collapsed counts exceed original full library')
  # Seurat disallows underscores/pipes. Keep an explicit reversible feature-key map;
  # output symbols are the original supplied symbols, never silently renamed.
  genes$analysis_key <- make.unique(gsub('[_|]', '-', genes$feature_id), sep = '.duplicate')
  rownames(counts) <- genes$analysis_key; colnames(counts) <- cells$nucleus_id
  rownames(cells) <- cells$nucleus_id
  for (s in unique(cells$sample_id)) if (length(unique(cells$condition[cells$sample_id == s])) != 1L) stop('Mouse appears in multiple genotypes')
  list(counts = counts, genes = genes, cells = cells, files = files, endpoints = expected)
}

log_normalize <- function(counts, original_library) {
  x <- counts %*% Matrix::Diagonal(x = 10000 / original_library)
  dimnames(x) <- dimnames(counts)
  x@x <- log1p(x@x); x
}
make_seurat <- function(inp) {
  obj <- Seurat::CreateSeuratObject(inp$counts, meta.data = inp$cells, min.cells = 0, min.features = 0, assay = 'RNA')
  if (!identical(rownames(obj), rownames(inp$counts))) stop('Unexpected Seurat feature renaming')
  obj <- SeuratObject::SetAssayData(obj, assay = 'RNA', layer = 'data', new.data = log_normalize(inp$counts, inp$cells$library_umi_full_genes))
  Seurat::Idents(obj) <- 'endpoint'; obj
}

select_genes <- function(inp, fraction) {
  fractions <- sapply(inp$endpoints, function(ep) Matrix::rowMeans(inp$counts[, inp$cells$endpoint == ep, drop = FALSE] > 0))
  missing_symbol <- grepl('^__gene_id__', inp$genes$feature_id) | is.na(inp$genes$gene_symbol) | !nzchar(inp$genes$gene_symbol)
  if ('symbol_missing' %in% names(inp$genes)) missing_symbol <- missing_symbol | tolower(as.character(inp$genes$symbol_missing)) %in% c('true', '1')
  detected <- apply(fractions, 1L, max) >= fraction
  eligible <- detected & !missing_symbol
  list(keys = rownames(inp$counts)[eligible], table = cbind(inp$genes, as.data.frame(fractions),
    passes_detection_fraction = detected, excluded_missing_symbol = missing_symbol, eligible = eligible))
}

metacells <- function(obj, inp, par, cfg, outdir, name) {
  # HVGs are used ONLY for the neighborhood PCA, not for network-gene selection.
  obj <- Seurat::FindVariableFeatures(obj, nfeatures = cfg$neighbor_hvgs, verbose = FALSE)
  obj <- Seurat::ScaleData(obj, features = Seurat::VariableFeatures(obj), verbose = FALSE)
  npcs <- min(cfg$neighbor_pcs, length(Seurat::VariableFeatures(obj)) - 1L, ncol(obj) - 1L)
  obj <- Seurat::RunPCA(obj, features = Seurat::VariableFeatures(obj), npcs = npcs, verbose = FALSE, seed.use = cfg$seed)
  grouping <- split(seq_len(nrow(inp$cells)), paste(inp$cells$sample_id, inp$cells$endpoint, sep = '::'))
  pieces <- list(); coverage <- list()
  for (group in names(grouping)) {
    idx <- grouping[[group]]; meta <- inp$cells[idx[1L], , drop = FALSE]
    reason <- if (length(idx) < par$min_cells) 'below_min_cells' else NA_character_
    piece <- NULL
    if (is.na(reason)) {
      # Separate official calls isolate failed small groups: v0.4.12's global
      # wrapper fails to remove exactly one NULL result from its list.
      attempt <- tryCatch({
        cur <- hdWGCNA::MetacellsByGroups(obj, group.by = c('sample_id', 'endpoint'), ident.group = 'endpoint',
          cells.use = inp$cells$nucleus_id[idx], k = par$k, max_shared = par$max_shared,
          min_cells = par$min_cells, target_metacells = par$target_metacells,
          max_iter = cfg$max_iter, reduction = 'pca', dims = seq_len(npcs), assay = 'RNA', layer = 'counts', mode = 'sum', wgcna_name = name)
        hdWGCNA::GetMetacellObject(cur, wgcna_name = name)
      }, error = function(e) { reason <<- paste0('native_metacell_failure: ', conditionMessage(e)); NULL })
      if (inherits(attempt, 'Seurat') && ncol(attempt) >= 2L) piece <- attempt else if (is.na(reason)) reason <- 'fewer_than_two_metacells'
    }
    coverage[[group]] <- data.frame(sample_id = meta$sample_id, condition = meta$condition, endpoint = meta$endpoint,
      n_nuclei = length(idx), n_metacells = if (is.null(piece)) 0L else ncol(piece), discovery_status = if (is.null(piece)) reason else 'included',
      all_nuclei_scored = TRUE)
    if (!is.null(piece)) pieces[[group]] <- piece
  }
  coverage <- do.call(rbind, coverage); rownames(coverage) <- NULL
  write_tsv(coverage, file.path(outdir, 'metacell_coverage.tsv'))
  if (!length(pieces)) stop('No valid metacell groups')
  mc <- if (length(pieces) == 1L) pieces[[1L]] else merge(pieces[[1L]], y = pieces[-1L])
  mc <- SeuratObject::JoinLayers(mc, assay = 'RNA')
  mc$condition <- rep(NA_character_, ncol(mc))
  members <- strsplit(as.character(mc$cells_merged), ',', fixed = TRUE); names(members) <- colnames(mc)
  records <- vector('list', length(members)); den <- numeric(length(members))
  for (j in seq_along(members)) {
    ix <- match(members[[j]], inp$cells$nucleus_id)
    if (anyNA(ix) || length(ix) != par$k || anyDuplicated(ix)) stop('Invalid native metacell membership')
    meta <- inp$cells[ix, , drop = FALSE]
    if (length(unique(meta$sample_id)) != 1L || length(unique(meta$endpoint)) != 1L || length(unique(meta$condition)) != 1L) stop('Mixed mouse/endpoint/genotype metacell')
    mc$condition[j] <- meta$condition[1L]; mc$sample_id[j] <- meta$sample_id[1L]; mc$endpoint[j] <- meta$endpoint[1L]
    den[j] <- sum(meta$library_umi_full_genes)
    records[[j]] <- data.frame(metacell_id = colnames(mc)[j], nucleus_id = meta$nucleus_id, sample_id = meta$sample_id,
      endpoint = meta$endpoint, condition = meta$condition)
  }
  mc$library_umi_full_genes <- den
  mc <- SeuratObject::SetAssayData(mc, assay = 'RNA', layer = 'data',
    new.data = log_normalize(SeuratObject::LayerData(mc, assay = 'RNA', layer = 'counts'), den))
  write_tsv(do.call(rbind, records), file.path(outdir, 'metacell_membership.tsv.gz'))
  write_tsv(cbind(metacell_id = colnames(mc), mc[[]]), file.path(outdir, 'metacell_metadata.tsv'))
  save_plot(ggplot2::ggplot(coverage, ggplot2::aes(sample_id, endpoint, fill = n_metacells)) + ggplot2::geom_tile(color = 'white') +
    ggplot2::geom_text(ggplot2::aes(label = n_metacells), size = 3) + ggplot2::scale_fill_gradient(low = '#edf4f8', high = '#285d86') +
    ggplot2::labs(title = paste(name, 'metacell discovery coverage'), x = 'Independent mouse', y = NULL, fill = 'Metacells') +
    ggplot2::theme_minimal(base_size = 11) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1)), file.path(outdir, 'figures/metacell_coverage'), height = 4)
  obj <- hdWGCNA::SetMetacellObject(obj, mc, wgcna_name = name)
  list(obj = obj, mc = mc, members = members, coverage = coverage)
}

power_choice <- function(tab, cfg) {
  good <- is.finite(tab$SFT.R.sq) & tab$SFT.R.sq >= cfg$topology_fit_min & is.finite(tab$slope) & tab$slope < 0 &
    is.finite(tab$mean.k.) & tab$mean.k. >= cfg$mean_connectivity_min
  if (any(good)) return(list(power = min(tab$Power[good]), status = 'PASSED_PREDEFINED_TOPOLOGY_GATE'))
  if (!is.null(cfg$exploratory_power)) {
    p <- as.numeric(cfg$exploratory_power)
    if (length(p) != 1L || !p %in% tab$Power) stop('Exploratory power must be one explicitly tested power')
    return(list(power = p, status = 'EXPLORATORY_EXPLICIT_POWER_POOR_TOPOLOGY'))
  }
  list(power = NULL, status = 'STOPPED_POOR_TOPOLOGY_NO_AUTOMATIC_FALLBACK')
}
power_outputs <- function(tab, cfg, outdir) {
  write_tsv(tab, file.path(outdir, 'soft_power_diagnostics.tsv'))
  if (!'group' %in% names(tab)) tab$group <- 'pooled'
  long <- do.call(rbind, lapply(c('SFT.R.sq', 'slope', 'mean.k.'), function(key) data.frame(Power = tab$Power, group = tab$group, metric = key, value = tab[[key]])))
  p <- ggplot2::ggplot(long, ggplot2::aes(Power, value, color = group)) + ggplot2::geom_line() + ggplot2::geom_point() +
    ggplot2::facet_wrap(~metric, scales = 'free_y', nrow = 1) + ggplot2::labs(title = 'Signed bicor soft-power diagnostics',
      subtitle = paste0('Gate: fit >= ', cfg$topology_fit_min, '; slope < 0; mean connectivity >= ', cfg$mean_connectivity_min), y = NULL) + ggplot2::theme_minimal(base_size = 11)
  save_plot(p, file.path(outdir, 'figures/soft_power_diagnostics'), width = 13, height = 4)
}

fit_mouse <- function(x, value = 'eigengene_mean') {
  x$condition <- factor(x$condition, levels = c('WT', 'cKO'))
  nwt <- sum(x$condition == 'WT'); nko <- sum(x$condition == 'cKO')
  empty <- data.frame(n_WT = nwt, n_cKO = nko, effect_cKO_minus_WT = NA_real_, SE = NA_real_, statistic = NA_real_, pvalue = NA_real_, CI_low = NA_real_, CI_high = NA_real_, model_status = 'insufficient_or_constant')
  if (nwt < 3L || nko < 3L || !all(is.finite(x[[value]])) || stats::sd(x[[value]]) <= 0) return(empty)
  fit <- stats::lm(stats::reformulate('condition', response = value), data = x)
  sm <- summary(fit)$coefficients; if (!'conditioncKO' %in% rownames(sm)) return(empty)
  ci <- suppressMessages(stats::confint(fit, 'conditioncKO'))
  data.frame(n_WT = nwt, n_cKO = nko, effect_cKO_minus_WT = sm['conditioncKO', 1L], SE = sm['conditioncKO', 2L],
    statistic = sm['conditioncKO', 3L], pvalue = sm['conditioncKO', 4L], CI_low = ci[1L], CI_high = ci[2L], model_status = 'mouse_level_lm')
}
mouse_effects <- function(means, cfg, lineage, network) {
  groups <- split(means, interaction(means$endpoint, means$module, drop = TRUE))
  primary <- list(); sensitivity <- list(); loo <- list(); boots <- list()
  for (key in names(groups)) {
    x <- groups[[key]]; info <- data.frame(lineage = lineage, network = network, endpoint = x$endpoint[1], module = x$module[1])
    primary[[key]] <- cbind(info, fit_mouse(x))
    use <- rep(FALSE, nrow(x))
    if ('sequencing_partition' %in% names(x)) use <- use | x$sequencing_partition %in% cfg$same_partition_values
    if ('lane' %in% names(x)) use <- use | as.character(x$lane) %in% as.character(cfg$same_lane_values)
    sensitivity[[key]] <- cbind(info, partition = 'prespecified_lane6_7', fit_mouse(x[use, , drop = FALSE]))
    loo[[key]] <- do.call(rbind, lapply(unique(x$sample_id), function(s) cbind(info, omitted_mouse = s, fit_mouse(x[x$sample_id != s, , drop = FALSE]))))
    wt <- x$eigengene_mean[x$condition == 'WT']; ko <- x$eigengene_mean[x$condition == 'cKO']
    if (length(wt) >= 3L && length(ko) >= 3L && cfg$bootstrap_iterations > 0) {
      draws <- replicate(cfg$bootstrap_iterations, mean(sample(ko, length(ko), TRUE)) - mean(sample(wt, length(wt), TRUE)))
      qs <- stats::quantile(draws, c(.025, .975), names = FALSE)
      boots[[key]] <- cbind(info, bootstrap_iterations = cfg$bootstrap_iterations, CI_low = qs[1], CI_high = qs[2], fraction_positive = mean(draws > 0), resampling_unit = 'independent_mouse_within_genotype')
    }
  }
  list(primary = do.call(rbind, primary), sensitivity = do.call(rbind, sensitivity), loo = do.call(rbind, loo), bootstrap = if (length(boots)) do.call(rbind, boots) else data.frame())
}

fixed_module_stability <- function(expr, me, mc_meta, modules, cfg) {
  out <- list(); j <- 0L
  for (mod in setdiff(unique(as.character(modules$module)), 'grey')) {
    genes <- intersect(modules$gene_name[modules$module == mod], colnames(expr))
    if (length(genes) < 3L || !mod %in% colnames(me)) next
    full <- as.numeric(WGCNA::bicor(as.matrix(expr[, genes, drop = FALSE]), me[, mod], maxPOutliers = cfg$max_p_outliers, use = 'p', pearsonFallback = 'individual'))
    for (s in unique(mc_meta$sample_id)) {
      keep <- mc_meta$sample_id != s; if (sum(keep) < 10L) next
      reduced <- as.numeric(WGCNA::bicor(as.matrix(expr[keep, genes, drop = FALSE]), me[keep, mod], maxPOutliers = cfg$max_p_outliers, use = 'p', pearsonFallback = 'individual'))
      nh <- min(10L, length(genes)); a <- order(full, decreasing = TRUE)[seq_len(nh)]; b <- order(reduced, decreasing = TRUE)[seq_len(nh)]
      j <- j + 1L; out[[j]] <- data.frame(module = mod, omitted_mouse = s, n_metacells_remaining = sum(keep), n_genes = length(genes),
        fixed_kME_correlation = suppressWarnings(stats::cor(full, reduced, use = 'pairwise.complete.obs')), top10_hub_fraction_retained = length(intersect(a, b)) / nh,
        method = 'fixed_module_fixed_nuclear_PC1_projection_connectivity_LOMO', formal_module_preservation = FALSE)
    }
  }
  if (length(out)) do.call(rbind, out) else data.frame()
}

network_outputs <- function(obj, inp, cfg, outdir, name, lineage, network, mc, members, expr) {
  all_modules <- hdWGCNA::GetModules(obj, wgcna_name = name)
  assigned <- droplevels(all_modules[as.character(all_modules$module) != 'grey', , drop = FALSE])
  if (!nrow(assigned)) stop('No non-grey modules identified; no module association can be calculated')
  # Native code orients raw PC1 to correlate positively with mean RNA. No Harmony
  # on sample IDs or genotype is applied, and no genotype covariate is regressed.
  obj <- hdWGCNA::ModuleEigengenes(obj, modules = assigned, group.by.vars = NULL, vars.to.regress = NULL,
    pc_dim = 1, assay = 'RNA', npcs = min(10L, min(table(assigned$module)) - 1L), seed.use = cfg$seed, wgcna_name = name, verbose = FALSE)
  obj <- hdWGCNA::ModuleConnectivity(obj, group.by = 'endpoint', group_name = inp$endpoints,
    corFnc = 'bicor', corOptions = paste0("use='p', maxPOutliers=", cfg$max_p_outliers, ", pearsonFallback='individual'"),
    harmonized = FALSE, assay = 'RNA', layer = 'data', sparse = FALSE, reassign_modules = FALSE, wgcna_name = name)
  km <- hdWGCNA::GetModules(obj, wgcna_name = name)
  for (col in setdiff(names(km), names(all_modules))) all_modules[[col]] <- km[[col]][match(all_modules$gene_name, km$gene_name)]
  all_modules$kME <- vapply(seq_len(nrow(all_modules)), function(i) {
    key <- paste0('kME_', as.character(all_modules$module[i])); if (key %in% names(all_modules)) as.numeric(all_modules[[key]][i]) else NA_real_
  }, numeric(1))
  gi <- match(all_modules$gene_name, inp$genes$analysis_key)
  all_modules$gene_symbol <- inp$genes$gene_symbol[gi]; all_modules$gene_ids <- inp$genes$gene_ids[gi]
  all_modules$color_hex <- vapply(as.character(all_modules$color), function(col) {v<-grDevices::col2rgb(col)[,1];grDevices::rgb(v[1],v[2],v[3],maxColorValue=255)}, character(1))
  all_modules$assigned <- as.character(all_modules$module) != 'grey'
  write_tsv(all_modules, file.path(outdir, 'modules.tsv'))
  write_tsv(all_modules, file.path(outdir, 'module_membership.tsv'))
  hubs <- all_modules[all_modules$assigned & is.finite(all_modules$kME), , drop = FALSE]
  hubs <- hubs[order(as.character(hubs$module), -hubs$kME), , drop = FALSE]
  hubs$hub_rank <- ave(hubs$kME, as.character(hubs$module), FUN = seq_along)
  write_tsv(hubs, file.path(outdir, 'all_gene_hub_ranks.tsv'))
  top <- hubs[hubs$hub_rank <= cfg$n_hubs, , drop = FALSE]
  write_tsv(top, file.path(outdir, 'top_hubs.tsv'))
  write_tsv(top, file.path(outdir, 'hub_genes.tsv'))
  me <- as.matrix(hdWGCNA::GetMEs(obj, harmonized = FALSE, wgcna_name = name))
  me <- me[inp$cells$nucleus_id, , drop = FALSE]
  rna <- SeuratObject::LayerData(obj, assay = 'RNA', layer = 'data')
  module_rna <- do.call(cbind, lapply(colnames(me), function(mod) Matrix::colMeans(rna[all_modules$gene_name[as.character(all_modules$module) == mod], , drop = FALSE])))
  dimnames(module_rna) <- list(inp$cells$nucleus_id, colnames(me))
  orientation <- vapply(colnames(me), function(mod) stats::cor(me[, mod], module_rna[, mod]), numeric(1))
  if (any(orientation < -1e-8, na.rm = TRUE)) stop('Native eigengene orientation disagrees with positive mean RNA')
  write_tsv(data.frame(module = names(orientation), correlation_with_mean_RNA = orientation, harmonized = FALSE, sign = 'positive_mean_RNA'), file.path(outdir, 'eigengene_orientation.tsv'))
  write_tsv(cbind(inp$cells, as.data.frame(me)), file.path(outdir, 'nucleus_eigengenes.tsv.gz'))
  group <- split(seq_len(nrow(inp$cells)), paste(inp$cells$sample_id, inp$cells$endpoint, sep = '::'))
  means <- do.call(rbind, lapply(group, function(ix) {
    meta <- inp$cells[ix[1], intersect(c('sample_id', 'condition', 'endpoint', 'lane', 'sequencing_partition'), names(inp$cells)), drop = FALSE]
    cbind(meta[rep(1L, ncol(me)), , drop = FALSE], module = colnames(me), n_nuclei = length(ix),
      eigengene_mean = as.numeric(colMeans(me[ix, , drop = FALSE])), module_RNA_mean = as.numeric(colMeans(module_rna[ix, , drop = FALSE])))
  })); rownames(means) <- NULL
  means$sample <- means$sample_id; means$value <- means$eigengene_mean
  write_tsv(means, file.path(outdir, 'mouse_endpoint_eigengenes.tsv'))
  write_tsv(means, file.path(outdir, 'module_eigengenes_mouse.tsv'))
  effects <- mouse_effects(means, cfg, lineage, network)
  write_tsv(effects$primary, file.path(outdir, 'module_genotype_associations_unadjusted.tsv'))
  write_tsv(effects$primary, file.path(outdir, 'effects.tsv'))
  write_tsv(effects$sensitivity, file.path(outdir, 'module_genotype_same_partition_unadjusted.tsv'))
  write_tsv(effects$loo, file.path(outdir, 'module_genotype_leave_one_mouse_out.tsv'))
  write_tsv(effects$bootstrap, file.path(outdir, 'module_genotype_mouse_bootstrap.tsv'))
  mc_me <- do.call(rbind, lapply(members, function(ids) colMeans(me[ids, , drop = FALSE])))
  dimnames(mc_me) <- list(names(members), colnames(me))
  ex <- expr[rownames(mc_me), , drop = FALSE]; mc_meta <- mc[[]][rownames(mc_me), , drop = FALSE]
  stability <- fixed_module_stability(ex, mc_me, mc_meta, all_modules, cfg)
  write_tsv(stability, file.path(outdir, 'fixed_module_connectivity_stability.tsv'))
  # Bounded top-hub co-expression edge list: all eligible genes remain in discovery.
  edges <- list()
  for (mod in unique(as.character(top$module))) {
    keys <- intersect(top$gene_name[top$module == mod], colnames(ex)); if (length(keys) < 2L) next
    cor <- WGCNA::bicor(as.matrix(ex[, keys, drop = FALSE]), maxPOutliers = cfg$max_p_outliers, use = 'p', pearsonFallback = 'individual')
    ix <- which(upper.tri(cor), arr.ind = TRUE)
    edges[[mod]] <- data.frame(module = mod, gene1 = inp$genes$gene_symbol[match(keys[ix[, 1]], inp$genes$analysis_key)],
      gene2 = inp$genes$gene_symbol[match(keys[ix[, 2]], inp$genes$analysis_key)], bicor = cor[ix], edge_type = 'top_hub_RNA_coexpression_not_direct_binding')
  }
  write_tsv(if (length(edges)) do.call(rbind, edges) else data.frame(), file.path(outdir, 'top_hub_coexpression_edges.tsv'))
  write_tsv(if (length(edges)) do.call(rbind, edges) else data.frame(), file.path(outdir, 'network_edges.tsv'))
  label <- paste(lineage, network)
  heat <- means; heat$mouse_endpoint <- paste(heat$sample_id, heat$endpoint, sep = ' / ')
  heat$display_z <- ave(heat$eigengene_mean, heat$module, FUN = function(x) if (stats::sd(x) > 0) as.numeric(scale(x)) else rep(0, length(x)))
  save_plot(ggplot2::ggplot(heat, ggplot2::aes(mouse_endpoint, module, fill = display_z)) + ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(low = '#397eb5', mid = 'white', high = '#c63d4f') + ggplot2::labs(title = paste(label, 'mouse-level eigengenes'), subtitle = 'Display z-score only; raw positive-RNA-oriented PC1 retained in tables', x = NULL, y = NULL, fill = 'ME z') +
    ggplot2::theme_minimal(base_size = 10) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)), file.path(outdir, 'figures/mouse_eigengene_heatmap'), width = 14, height = max(5, min(20, ncol(me) * .23)))
  mr <- aggregate(module_RNA_mean ~ endpoint + condition + module, means, mean)
  save_plot(ggplot2::ggplot(mr, ggplot2::aes(paste(endpoint, condition), module, color = module_RNA_mean, size = module_RNA_mean)) + ggplot2::geom_point() + ggplot2::scale_color_gradient(low = '#d8e6f1', high = '#c63d4f') +
    ggplot2::labs(title = paste(label, 'module RNA'), subtitle = 'Equal-mouse averages of mean log1p(CP10k) across module genes', x = NULL, y = NULL, color = 'Mean RNA', size = 'Mean RNA') + ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1)), file.path(outdir, 'figures/module_RNA_dotplot'), height = max(5, min(20, ncol(me) * .23)))
  plot_mods <- names(sort(table(hubs$module), decreasing = TRUE))[seq_len(min(cfg$n_hub_plot_modules, length(unique(hubs$module))))]
  hp <- top[as.character(top$module) %in% plot_mods & top$hub_rank <= 10L, , drop = FALSE]
  hp$gene_panel <- paste(hp$module, hp$gene_symbol, sep = ': ')
  save_plot(ggplot2::ggplot(hp, ggplot2::aes(kME, reorder(gene_panel, kME), fill = color)) + ggplot2::geom_col() + ggplot2::scale_fill_identity() + ggplot2::facet_wrap(~module, scales = 'free_y', ncol = 3) +
    ggplot2::labs(title = paste(label, 'top RNA co-expression hubs'), subtitle = 'Bicor kME; this is not TF binding evidence', x = 'Module membership correlation (kME)', y = NULL) + ggplot2::theme_minimal(base_size = 9), file.path(outdir, 'figures/top_hubs'), height = max(6, ceiling(length(plot_mods) / 3) * 3))
  net <- hdWGCNA::GetNetworkData(obj, wgcna_name = name)
  for (b in seq_along(net$dendrograms)) {
    draw <- function() {
      ids <- if (!is.null(net$blockGenes)) net$blockGenes[[b]] else seq_len(nrow(all_modules))
      colors <- setNames(as.character(all_modules$color), all_modules$gene_name)[colnames(expr)[ids]]
      WGCNA::plotDendroAndColors(net$dendrograms[[b]], colors, groupLabels = 'Modules', dendroLabels = FALSE, hang = .03, main = paste(label, 'block', b))
    }
    stem <- file.path(outdir, 'figures', paste0('gene_dendrogram_block_', b))
    grDevices::cairo_pdf(paste0(stem, '.pdf'), width = 12, height = 5); tryCatch(draw(), finally = grDevices::dev.off())
    grDevices::svg(paste0(stem, '.svg'), width = 12, height = 5); tryCatch(draw(), finally = grDevices::dev.off())
    grDevices::png(paste0(stem, '.png'), width = 2160, height = 900, res = 180); tryCatch(draw(), finally = grDevices::dev.off())
  }
  # Keep recoverable scientific state without duplicating the full nucleus counts.
  state <- list(modules = all_modules, MEs = me, mouse_means = means, effects = effects, metacell_object = mc,
    metacell_members = members, datExpr = expr, fixed_metacell_MEs = mc_me, network = net,
    eligible_genes = hdWGCNA::GetWGCNAGenes(obj, wgcna_name = name), settings = cfg,
    interpretation = 'Exploratory RNA coexpression; no causal binding or formal preservation inference; no genotype DEGs generated')
  save_rds(state, file.path(outdir, 'analysis_state.rds'))
  write_json(list(status = 'NETWORK_COMPLETE', lineage = lineage, network = network, n_nuclei_scored = nrow(me), n_metacells = ncol(mc),
    n_tested_genes = ncol(expr), n_modules = ncol(me), n_grey_genes = sum(!all_modules$assigned),
    TOM_saved = FALSE, block_size = cfg$max_block_size, harmonized_eigengenes = FALSE,
    confounding = 'Genotype partially confounded with technical partition; mouse-level association is exploratory',
    preservation = 'Fixed-module leave-one-mouse-out connectivity only; formal WGCNA preservation not computed'), file.path(outdir, 'status.json'))
  invisible(state)
}

construct <- function(obj, cfg, outdir, name, power, consensus = FALSE) {
  previous <- getwd(); setwd(outdir); on.exit(setwd(previous), add = TRUE)
  obj <- hdWGCNA::ConstructNetwork(obj, soft_power = power, consensus = consensus, tom_outdir = 'TOM', tom_name = name,
    overwrite_tom = FALSE, store_tom_in_seurat = FALSE, wgcna_name = name,
    maxBlockSize = cfg$max_block_size, randomSeed = cfg$seed, corType = 'bicor', networkType = 'signed', TOMType = 'signed',
    maxPOutliers = cfg$max_p_outliers, pearsonFallback = 'individual', minModuleSize = cfg$min_module_size,
    mergeCutHeight = cfg$merge_cut_height, deepSplit = cfg$deep_split, pamStage = FALSE,
    saveConsensusTOMs = FALSE, saveIndividualTOMs = FALSE, useDiskCache = FALSE, verbose = 3)
  # v0.4.12 records a nonexistent TOM path when saving is disabled. Clear that
  # misleading pointer so later connectivity code cannot mistake it for a TOM.
  net <- hdWGCNA::GetNetworkData(obj, wgcna_name = name); net$TOMFiles <- character()
  hdWGCNA::SetNetworkData(obj, net, wgcna_name = name)
}

run_primary <- function(lineage, input_dir, rootout, cfg, script_hash) {
  outdir <- file.path(rootout, lineage, 'pooled'); dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  inp <- read_input(input_dir, lineage); fp <- fingerprint(list(inputs = md5(inp$files), config = cfg, script = script_hash, lineage = lineage, network = 'pooled'))
  ensure_identity(outdir, fp)
  if (checkpoint_valid(outdir, fp)) { message('Verified completed pooled checkpoint: ', lineage); return(invisible(TRUE)) }
  name <- paste0(lineage, '_pooled'); par <- cfg[[lineage]]
  write_json(list(input_md5 = as.list(md5(inp$files)), config = cfg, script_md5 = script_hash, upstream_sha = 'e3344d1f7bbac4264adf94f7aa31e0802fa8282d'), file.path(outdir, 'provenance.json'))
  selected <- select_genes(inp, cfg$expression_fraction); write_tsv(selected$table, file.path(outdir, 'network_gene_selection.tsv'))
  if (length(selected$keys) < 500L) stop('Fewer than 500 eligible full-RNA genes')
  obj <- make_seurat(inp); obj <- hdWGCNA::SetupForWGCNA(obj, features = selected$keys, wgcna_name = name)
  mc_result <- metacells(obj, inp, par, cfg, outdir, name); obj <- mc_result$obj; mc_result$obj <- NULL
  if (ncol(mc_result$mc) < cfg$min_metacells) stop('Too few metacells for discovery')
  obj <- hdWGCNA::SetDatExpr(obj, group_name = inp$endpoints, group.by = 'endpoint', use_metacells = TRUE, assay = 'RNA', layer = 'data', wgcna_name = name)
  expr <- hdWGCNA::GetDatExpr(obj, wgcna_name = name)
  selected$table$retained_after_goodGenes <- selected$table$analysis_key %in% colnames(expr)
  write_tsv(selected$table, file.path(outdir, 'network_genes.tsv'))
  write_tsv(data.frame(analysis_key = colnames(expr), retained_after_goodGenes = TRUE), file.path(outdir, 'network_genes_tested.tsv'))
  obj <- hdWGCNA::TestSoftPowers(obj, powers = cfg$powers, networkType = 'signed', corFnc = 'bicor',
    corOptions = list(use = 'p', maxPOutliers = cfg$max_p_outliers, pearsonFallback = 'individual'), blockSize = cfg$max_block_size, wgcna_name = name)
  pt <- hdWGCNA::GetPowerTable(obj, wgcna_name = name); power_outputs(pt, cfg, outdir); power <- power_choice(pt, cfg)
  write_json(power, file.path(outdir, 'selected_power.json'))
  if (is.null(power$power)) { write_json(list(status = power$status, lineage = lineage, modules_fitted = FALSE), file.path(outdir, 'status.json')); return(invisible(FALSE)) }
  obj <- construct(obj, cfg, outdir, name, power$power)
  network_outputs(obj, inp, cfg, outdir, name, lineage, 'pooled', mc_result$mc, mc_result$members, expr)
  complete(outdir, fp, list(lineage = lineage, network = 'pooled', topology_status = power$status))
  invisible(TRUE)
}

run_consensus <- function(lineage, input_dir, rootout, cfg, script_hash) {
  pooled <- file.path(rootout, lineage, 'pooled'); outdir <- file.path(rootout, lineage, 'consensus')
  if (!file.exists(file.path(pooled, 'completed.json'))) {
    write_json(list(status = 'NOT_RUN_PRIMARY_NETWORK_INCOMPLETE'), file.path(rootout, paste0(lineage, '_consensus_not_run.json'))); return(invisible(FALSE))
  }
  state_path <- file.path(pooled, 'analysis_state.rds')
  input_files <- file.path(input_dir, lineage, c('counts.mtx.gz', 'genes.tsv', 'cells.tsv'))
  fp <- fingerprint(list(inputs = md5(input_files), pooled_state = md5(state_path), config = cfg, script = script_hash, network = 'consensus'))
  ensure_identity(outdir, fp)
  if (checkpoint_valid(outdir, fp)) { message('Verified completed consensus checkpoint: ', lineage); return(invisible(TRUE)) }
  state <- readRDS(state_path); mc <- state$metacell_object
  support <- do.call(rbind, lapply(c('WT', 'cKO'), function(g) {
    ct <- table(mc$sample_id[mc$condition == g]); data.frame(condition = g, n_metacells = sum(ct),
      n_eligible_mice = sum(ct >= cfg$consensus_min_metacells_per_mouse))
  }))
  write_tsv(support, file.path(outdir, 'consensus_coverage.tsv'))
  if (any(support$n_metacells < cfg$consensus_min_metacells_per_condition | support$n_eligible_mice < cfg$consensus_min_mice_per_condition)) {
    write_json(list(status = 'NOT_RUN_INSUFFICIENT_CONSENSUS_COVERAGE', minimum_mc = cfg$consensus_min_metacells_per_condition,
      minimum_mice = cfg$consensus_min_mice_per_condition), file.path(outdir, 'status.json')); return(invisible(FALSE))
  }
  inp <- read_input(input_dir, lineage)
  name <- paste0(lineage, '_consensus'); obj <- make_seurat(inp)
  obj <- hdWGCNA::SetupForWGCNA(obj, features = state$eligible_genes, wgcna_name = name)
  obj <- hdWGCNA::SetMetacellObject(obj, mc, wgcna_name = name)
  obj <- hdWGCNA::SetMultiExpr(obj, group_name = inp$endpoints, group.by = 'endpoint', multi.group.by = 'condition', multi_groups = c('WT', 'cKO'),
    use_metacells = TRUE, assay = 'RNA', layer = 'data', wgcna_name = name)
  obj <- hdWGCNA::TestSoftPowersConsensus(obj, powers = cfg$powers, networkType = 'signed', corFnc = 'bicor', setDatExpr = FALSE,
    corOptions = list(use = 'p', maxPOutliers = cfg$max_p_outliers, pearsonFallback = 'individual'), blockSize = cfg$max_block_size, wgcna_name = name)
  pt <- hdWGCNA::GetPowerTable(obj, wgcna_name = name); power_outputs(pt, cfg, outdir)
  multi <- hdWGCNA::GetMultiExpr(obj, wgcna_name = name)
  genes <- read.delim(file.path(pooled, 'network_genes.tsv'), stringsAsFactors = FALSE, check.names = FALSE)
  genes$retained_after_goodGenes <- genes$analysis_key %in% colnames(multi[[1L]]$data)
  write_tsv(genes, file.path(outdir, 'network_genes.tsv'))
  choices <- lapply(names(multi), function(g) power_choice(pt[pt$group == g, , drop = FALSE], cfg)); names(choices) <- names(multi)
  write_json(choices, file.path(outdir, 'selected_power.json'))
  if (any(vapply(choices, function(x) is.null(x$power), logical(1)))) {
    write_json(list(status = 'NOT_RUN_POOR_CONSENSUS_TOPOLOGY', modules_fitted = FALSE), file.path(outdir, 'status.json')); return(invisible(FALSE))
  }
  obj <- construct(obj, cfg, outdir, name, vapply(choices, function(x) x$power, numeric(1)), consensus = TRUE)
  expr <- do.call(rbind, lapply(multi, function(x) as.data.frame(x$data)))
  # rbind prefixes set names; preserve original metacell IDs by matching data rows.
  rownames(expr) <- unlist(lapply(multi, function(x) rownames(x$data)), use.names = FALSE)
  expr <- expr[colnames(mc), , drop = FALSE]
  network_outputs(obj, inp, cfg, outdir, name, lineage, 'consensus', mc, state$metacell_members, expr)
  pmod <- state$modules; cmod <- read.delim(file.path(outdir, 'modules.tsv'), stringsAsFactors = FALSE)
  overlap <- merge(pmod[, c('gene_name', 'module')], cmod[, c('gene_name', 'module')], by = 'gene_name', suffixes = c('_pooled', '_consensus'))
  write_tsv(as.data.frame(table(overlap$module_pooled, overlap$module_consensus)), file.path(outdir, 'pooled_consensus_module_overlap.tsv'))
  complete(outdir, fp, list(lineage = lineage, network = 'consensus', primary_analysis = FALSE))
  invisible(TRUE)
}

final_associations <- function(outdir, cfg) {
  for (network in c('pooled', 'consensus')) for (kind in c('module_genotype_associations', 'module_genotype_same_partition')) {
    paths <- file.path(outdir, c('CA1', 'IN'), network, paste0(kind, '_unadjusted.tsv'))
    paths <- paths[file.exists(paths)]
    if (!length(paths)) next
    combined <- do.call(rbind, lapply(paths, read.delim, check.names = FALSE, stringsAsFactors = FALSE))
    combined$padj <- stats::p.adjust(combined$pvalue, method = 'BH')
    combined$fdr_family <- paste(network, kind, 'all tested CA1+IN module x endpoint associations')
    write_tsv(combined, file.path(outdir, paste0(network, '_', kind, '.tsv')))
    good <- combined[is.finite(combined$effect_cKO_minus_WT), , drop = FALSE]
    if (!nrow(good)) next
    good$label <- paste(good$lineage, good$endpoint, good$module, sep = ' / ')
    good$significant <- !is.na(good$padj) & good$padj < .05
    p <- ggplot2::ggplot(good, ggplot2::aes(effect_cKO_minus_WT, reorder(label, effect_cKO_minus_WT), color = significant)) +
      ggplot2::geom_vline(xintercept = 0, color = '#bfc7cc') + ggplot2::geom_segment(ggplot2::aes(x = CI_low, xend = CI_high, yend = label), linewidth = .5) + ggplot2::geom_point(size = 1.7) +
      ggplot2::scale_color_manual(values = c('FALSE' = '#8295a0', 'TRUE' = '#c63d4f')) + ggplot2::labs(title = paste(network, 'mouse-level module associations'),
        subtitle = 'Exploratory continuous PC1 contrasts; technical partition partly confounds genotype; these are not DEGs', x = 'cKO - WT mean module eigengene (95% CI)', y = NULL) + ggplot2::theme_minimal(base_size = 9)
    save_plot(p, file.path(outdir, 'figures', paste0(network, '_', kind)), width = 12, height = max(6, min(30, nrow(good) * .17)))
  }
}

main <- function() {
  args <- arguments(commandArgs(trailingOnly = TRUE)); raw <- jsonlite::read_json(args$config, simplifyVector = TRUE)
  provided <- raw$hdwgcna %||% list()
  unknown <- setdiff(names(provided), c(names(defaults()), 'enabled', 'version', 'commit', 'network_type', 'cor_type'))
  if (length(unknown)) stop('Unrecognized hdwgcna settings (no silent aliases): ', paste(unknown, collapse = ', '))
  if (!identical(provided$network_type %||% 'signed', 'signed') || !identical(provided$cor_type %||% 'bicor', 'bicor')) stop('Approved network must be signed bicor')
  if (!identical(provided$version %||% '0.4.12', '0.4.12') || !identical(provided$commit %||% 'e3344d1f7bbac4264adf94f7aa31e0802fa8282d', 'e3344d1f7bbac4264adf94f7aa31e0802fa8282d')) stop('Config upstream version/commit differs from pinned source')
  cfg <- utils::modifyList(defaults(), provided)
  if (!isTRUE(cfg$execute)) stop('hdwgcna.execute must be true in the immutable approval-bound workflow configuration')
  if (as.integer(cfg$max_block_size) > 4000L) stop('max_block_size must remain <=4000 for this approved resource plan')
  required <- c('Seurat', 'SeuratObject', 'Matrix', 'hdWGCNA', 'WGCNA', 'ggplot2', 'jsonlite', 'digest')
  for (p in required) if (!requireNamespace(p, quietly = TRUE)) stop('Missing package: ', p, '; install through the separate managed environment setup')
  if (as.character(utils::packageVersion('hdWGCNA')) != '0.4.12') stop('Pinned hdWGCNA 0.4.12 required')
  desc <- utils::packageDescription('hdWGCNA'); sha <- desc$RemoteSha %||% desc$GithubSHA1
  if (!is.null(sha) && !startsWith('e3344d1f7bbac4264adf94f7aa31e0802fa8282d', sha)) stop('hdWGCNA source SHA differs from approved pin')
  cfg$runtime_versions <- setNames(lapply(required, function(p) as.character(utils::packageVersion(p))), required)
  outdir <- normalizePath(args$outdir, mustWork = FALSE)
  set.seed(cfg$seed); options(future.globals.maxSize = 8 * 1024^3)
  WGCNA::enableWGCNAThreads(nThreads = max(2L, as.integer(cfg$threads)))
  ggplot2::theme_set(ggplot2::theme_minimal(base_size = 11))
  script_arg <- grep('^--file=', commandArgs(), value = TRUE)
  script_hash <- if (length(script_arg)) unname(md5(sub('^--file=', '', script_arg[1L]))) else 'unavailable'
  inputs <- unlist(lapply(c('CA1', 'IN'), function(lin) file.path(args[['input-dir']], lin, c('counts.mtx.gz', 'genes.tsv', 'cells.tsv'))))
  if (any(!file.exists(inputs))) stop('Missing prepared lineage inputs')
  ensure_identity(outdir, fingerprint(list(input_hashes = as.character(tools::md5sum(inputs)), config = cfg, script = script_hash)))
  capture.output(sessionInfo(), file = file.path(outdir, 'sessionInfo.txt'))
  write_json(list(config = cfg, upstream_sha = 'e3344d1f7bbac4264adf94f7aa31e0802fa8282d', installed_remote_sha = sha,
    script_md5 = script_hash, all_gene_DESeq2_results_reused_without_refit = TRUE, new_gene_DEGs_generated = FALSE), file.path(outdir, 'run_parameters.json'))
  errors <- list(); optional_errors <- list()
  # Finish both primary networks before attempting either conditional consensus.
  for (lineage in c('CA1', 'IN')) {
    tryCatch(run_primary(lineage, args[['input-dir']], outdir, cfg, script_hash), error = function(e) {
      errors[[paste(lineage, 'pooled')]] <<- conditionMessage(e)
      d <- file.path(outdir, lineage, 'pooled'); dir.create(d, recursive = TRUE, showWarnings = FALSE)
      write_json(list(status = 'FAILED', message = conditionMessage(e)), file.path(d, 'failure.json'))
    }); gc()
  }
  if (isTRUE(cfg$consensus)) for (lineage in c('CA1', 'IN')) {
    tryCatch(run_consensus(lineage, args[['input-dir']], outdir, cfg, script_hash), error = function(e) {
      optional_errors[[paste(lineage, 'consensus')]] <<- conditionMessage(e)
      d <- file.path(outdir, lineage, 'consensus'); dir.create(d, recursive = TRUE, showWarnings = FALSE)
      write_json(list(status = 'FAILED_OPTIONAL_CONSENSUS', message = conditionMessage(e)), file.path(d, 'failure.json'))
    }); gc()
  }
  final_associations(outdir, cfg)
  final_status <- if (length(errors)) 'FINISHED_WITH_PRIMARY_FAILURES' else if (length(optional_errors)) 'FINISHED_WITH_OPTIONAL_FAILURES' else 'FINISHED_CHECK_INDIVIDUAL_TOPOLOGY_STATUSES'
  write_json(list(status = final_status, primary_failures = errors, optional_consensus_failures = optional_errors,
    note = 'No formal module preservation; no new DESeq2 fits; mouse-level module associations are exploratory.'), file.path(outdir, 'workflow_status.json'))
  if (length(errors)) stop('One or more network branches failed; see workflow_status.json')
  if (length(optional_errors)) warning('Optional consensus failed; verified primary results are retained. See workflow_status.json')
}

if (sys.nframe() == 0L) main()
