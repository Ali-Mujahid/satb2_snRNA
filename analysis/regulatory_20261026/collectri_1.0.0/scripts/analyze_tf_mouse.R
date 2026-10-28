#!/usr/bin/env Rscript
# Continuous mouse-level CollecTRI ULM associations; never a gene-DE analysis.
# Execute only through the reviewed native workflow. Import/source for synthetic
# tests does not execute main(). No package installation or online requests.

`%||%` <- function(x, y) if (is.null(x)) y else x
ENDPOINTS <- c('CA1_Deep', 'CA1_Superficial', 'IN_C1', 'IN_C6', 'IN_C8')
PRIMARY_MODE <- 'mouse_log1p_original_sizefactor_expression'
NOSELF_MODE <- 'mouse_log1p_original_sizefactor_expression_no_Satb2_self_target'

policy_defaults <- function() list(
  seed = 1337L, bootstrap_iterations = 200L,
  same_partition = '23HJT3LT4_lanes6_7',
  expected_primary_WT = 6L, expected_primary_cKO = 6L,
  expected_partition_WT = 6L, expected_partition_cKO = 3L,
  variation_relative_tolerance = 1e-12, residual_relative_tolerance = 1e-12,
  heatmap_top_tfs = 30L, boxplot_top_tfs = 5L, volcano_labels_per_endpoint = 3L,
  satb2_tf = 'Satb2', require_satb2_no_self = TRUE, plots = TRUE
)
need <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
sha256 <- function(path) digest::digest(file = path, algo = 'sha256', serialize = FALSE)
receipt <- function(path) {
  need(file.exists(path) && !dir.exists(path) && !nzchar(Sys.readlink(path)), paste('Invalid regular input file:', path))
  list(path = normalizePath(path, mustWork = TRUE), bytes = unname(file.info(path)$size), sha256 = sha256(path))
}
check_record <- function(record, label) {
  need(is.list(record) && all(c('path', 'bytes', 'sha256') %in% names(record)), paste('Missing bound input record:', label))
  got <- receipt(record$path)
  need(identical(got$sha256, record$sha256) && isTRUE(got$bytes == record$bytes), paste('Changed input:', label))
  got
}
json_write <- function(x, path) {
  need(!file.exists(path), paste('Refusing to overwrite:', path))
  jsonlite::write_json(x, path, pretty = TRUE, auto_unbox = TRUE, null = 'null', na = 'null', digits = NA)
}
write_tsv <- function(x, path) {
  need(!file.exists(path), paste('Refusing to overwrite:', path))
  write.table(x, path, sep = '\t', quote = FALSE, row.names = FALSE, na = 'NA')
}
read_tsv <- function(path, columns) {
  x <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, na.strings = c('NA', 'NaN', ''))
  need(!anyDuplicated(names(x)) && all(columns %in% names(x)), paste('Missing/duplicate columns in', path))
  x
}
as_flag <- function(x) {
  y <- tolower(as.character(x)); need(all(y %in% c('true', 'false', '1', '0')), 'Invalid Boolean coverage/eligibility value')
  y %in% c('true', '1')
}
args_parse <- function(x) {
  need(length(x) >= 6L && length(x) %% 2L == 0L, 'Expected --key value arguments')
  keys <- sub('^--', '', x[seq.int(1L, length(x), 2L)])
  out <- as.list(x[seq.int(2L, length(x), 2L)]); names(out) <- keys
  need(!anyDuplicated(keys) && setequal(keys, c('config', 'input-dir', 'outdir')), 'Use --config ABS --input-dir ABS --outdir ABS')
  need(all(vapply(out, function(p) grepl('^/', p), logical(1))), 'Paths must be absolute')
  out
}

sample_partition <- function(cfg, roster) {
  candidate <- list(samples = check_record(cfg$samples, 'samples'))
  for (name in names(cfg$additional_hash_inputs %||% list())) if (grepl('mouse.*coverage|sample.*coverage', name, ignore.case = TRUE)) {
    candidate[[name]] <- check_record(cfg$additional_hash_inputs[[name]], name)
  }
  for (name in names(candidate)) {
    tab <- read_tsv(candidate[[name]]$path, 'sample')
    partition_col <- intersect(c('sequencing_partition', 'partition'), names(tab))
    if (!length(partition_col)) next
    tab <- tab[tab$sample %in% roster$sample, , drop = FALSE]
    values <- lapply(roster$sample, function(s) unique(as.character(tab[[partition_col[1L]]][tab$sample == s])))
    if (any(lengths(values) != 1L) || anyNA(unlist(values)) || any(!nzchar(unlist(values)))) next
    return(list(partition = unlist(values, use.names = FALSE), source = candidate[[name]], column = partition_col[1L]))
  }
  stop('No unambiguous sample-level partition in bound samples/coverage metadata', call. = FALSE)
}

model_score <- function(x, cfg, expected_WT = NULL, expected_cKO = NULL) {
  nwt <- sum(x$condition == 'WT'); nko <- sum(x$condition == 'cKO')
  wt <- x$escore[x$condition == 'WT']; ko <- x$escore[x$condition == 'cKO']
  base <- data.frame(n_WT = nwt, n_cKO = nko,
    mean_WT = if (length(wt)) mean(wt) else NA_real_, mean_cKO = if (length(ko)) mean(ko) else NA_real_,
    effect_cKO_minus_WT = if (length(wt) && length(ko)) mean(ko) - mean(wt) else NA_real_,
    SE = NA_real_, statistic = NA_real_, pvalue = NA_real_, CI_low = NA_real_, CI_high = NA_real_,
    residual_df = NA_integer_, tested = FALSE, model_status = 'NOT_TESTED')
  reject <- function(reason) {base$model_status <- reason; base}
  if (!is.null(expected_WT) && (nwt != expected_WT || nko != expected_cKO)) return(reject('UNEXPECTED_MOUSE_COVERAGE'))
  if (nwt < 3L || nko < 3L) return(reject('FEWER_THAN_3_INDEPENDENT_MICE_IN_ONE_GENOTYPE'))
  if (anyDuplicated(x$sample)) return(reject('DUPLICATE_MOUSE'))
  if (any(!is.finite(x$escore))) return(reject('NONFINITE_SCORE'))
  scale <- max(1, abs(x$escore))
  if (diff(range(x$escore)) <= cfg$variation_relative_tolerance * scale) return(reject('CONSTANT_OR_NUMERICALLY_CONSTANT_SCORE'))
  rss <- sum((wt - mean(wt))^2) + sum((ko - mean(ko))^2)
  if (rss <= nrow(x) * (cfg$residual_relative_tolerance * scale)^2) return(reject('ZERO_OR_NUMERICALLY_ZERO_WITHIN_GROUP_VARIANCE'))
  x$condition <- factor(x$condition, levels = c('WT', 'cKO'))
  fit <- stats::lm(escore ~ condition, data = x)
  cf <- summary(fit)$coefficients['conditioncKO', ]
  if (fit$rank != 2L || fit$df.residual != nrow(x) - 2L || any(!is.finite(cf))) return(reject('NONESTIMABLE_LINEAR_MODEL'))
  delta <- stats::qt(.975, fit$df.residual) * cf[2L]
  base$effect_cKO_minus_WT <- unname(cf[1L]); base$SE <- unname(cf[2L]); base$statistic <- unname(cf[3L]); base$pvalue <- unname(cf[4L])
  base$CI_low <- unname(cf[1L] - delta); base$CI_high <- unname(cf[1L] + delta); base$residual_df <- fit$df.residual
  base$tested <- TRUE; base$model_status <- 'MOUSE_LEVEL_CONTINUOUS_LM'; base
}

adjust_family <- function(x, family) {
  x$association_fdr <- NA_real_
  ix <- which(x$tested & is.finite(x$pvalue))
  if (length(ix)) x$association_fdr[ix] <- stats::p.adjust(x$pvalue[ix], method = 'BH')
  x$padj <- x$association_fdr; x$fdr_family <- family
  x$gene_set_fdr_used_for_filtering <- FALSE
  x
}

validate_native_scores <- function(scores, coverage, roster, mode) {
  need(all(scores$mode == mode), 'Supplementary and primary score modes must not be mixed')
  need(all(scores$endpoint %in% ENDPOINTS) && !anyDuplicated(scores[c('endpoint', 'TF', 'sample')]), 'Duplicate/unknown endpoint-TF-mouse score')
  need(all(scores$sample %in% roster$sample), 'Unknown mouse in score table')
  need(all(scores$condition == roster$condition[match(scores$sample, roster$sample)]), 'Mouse genotype differs from frozen roster')
  need(all(is.finite(scores$escore)) && all(is.finite(scores$gene_set_fdr)) && all(scores$gene_set_fdr >= 0 & scores$gene_set_fdr <= 1), 'Nonfinite native score/gene-set FDR')
  need(all(scores$score_unit == 'native_ULM_t_value'), 'Unexpected native score unit')
  need(all(as_flag(scores$eligible)) && all(scores$n_used_targets >= 5L) && all(scores$n_background_genes >= 6L), 'Invalid native TF score eligibility')
  need(!anyDuplicated(coverage[c('endpoint', 'TF')]) && setequal(unique(coverage$endpoint), ENDPOINTS), 'Coverage must identify all five exact endpoints without duplicate TF rows')
  for (i in seq_len(nrow(coverage))) {
    c <- coverage[i, ]; x <- scores[scores$endpoint == c$endpoint & scores$TF == c$TF, , drop = FALSE]
    eligible <- as_flag(c$eligible_tmin5) && c$n_used_targets >= 5L
    if (eligible) {
      need(nrow(x) == 12L && setequal(x$sample, roster$sample), paste('Eligible TF lacks all12 mice:', c$endpoint, c$TF))
      need(length(unique(x$n_used_targets)) == 1L && x$n_used_targets[1L] == c$n_used_targets &&
        length(unique(x$n_background_genes)) == 1L && x$n_background_genes[1L] == c$n_background_genes,
        paste('Score/coverage disagree:', c$endpoint, c$TF))
    } else need(nrow(x) == 0L, paste('Ineligible prior TF unexpectedly scored:', c$endpoint, c$TF))
  }
  k <- paste(coverage$endpoint, coverage$TF, sep = '::')
  need(all(paste(scores$endpoint, scores$TF, sep = '::') %in% k), 'Native scored TF missing from all-TF coverage')
  invisible(TRUE)
}

analyze_variant <- function(scores, coverage, roster, cfg, variant) {
  primary <- list(); sensitive <- list(); loo <- list(); boots <- list(); stability <- list()
  for (i in seq_len(nrow(coverage))) {
    cov <- coverage[i, , drop = FALSE]
    info <- data.frame(endpoint = cov$endpoint, endpoint_id = cov$endpoint_id, source_endpoint_id = cov$source_endpoint_id,
      TF = cov$TF, network_variant = variant, n_used_targets = cov$n_used_targets, n_background_genes = cov$n_background_genes,
      score_unit = 'native_ULM_t_value')
    x <- scores[scores$endpoint == cov$endpoint & scores$TF == cov$TF, , drop = FALSE]
    info$n_scored_mice <- nrow(x); info$n_expected_mice <- 12L
    if (!nrow(x)) {
      missing <- model_score(data.frame(sample = character(), condition = character(), escore = numeric()), cfg)
      missing$model_status <- paste0('NOT_MODELED_', cov$reason)
      primary[[i]] <- cbind(info, missing)
      sensitive[[i]] <- cbind(info, partition = cfg$same_partition, missing)
      next
    }
    fit <- model_score(x, cfg, cfg$expected_primary_WT, cfg$expected_primary_cKO)
    primary[[i]] <- cbind(info, fit)
    same <- x[x$sequencing_partition == cfg$same_partition, , drop = FALSE]
    sensitive[[i]] <- cbind(info, partition = cfg$same_partition,
      model_score(same, cfg, cfg$expected_partition_WT, cfg$expected_partition_cKO))
    if (!fit$tested) next
    omitted <- lapply(roster$sample, function(s) cbind(info, omitted_mouse = s, model_score(x[x$sample != s, , drop = FALSE], cfg)))
    cur_loo <- do.call(rbind, omitted); loo[[length(loo) + 1L]] <- cur_loo
    valid <- cur_loo$tested & is.finite(cur_loo$effect_cKO_minus_WT)
    base_sign <- sign(fit$effect_cKO_minus_WT)
    stability[[length(stability) + 1L]] <- cbind(info,
      primary_effect = fit$effect_cKO_minus_WT, n_loo_models = nrow(cur_loo), n_loo_tested = sum(valid),
      n_loo_same_sign = sum(sign(cur_loo$effect_cKO_minus_WT[valid]) == base_sign),
      loo_same_sign_fraction = if (any(valid)) mean(sign(cur_loo$effect_cKO_minus_WT[valid]) == base_sign) else NA_real_,
      loo_effect_min = if (any(valid)) min(cur_loo$effect_cKO_minus_WT[valid]) else NA_real_,
      loo_effect_max = if (any(valid)) max(cur_loo$effect_cKO_minus_WT[valid]) else NA_real_)
    wt <- x$escore[x$condition == 'WT']; ko <- x$escore[x$condition == 'cKO']
    draws <- replicate(cfg$bootstrap_iterations, mean(sample(ko, length(ko), replace = TRUE)) - mean(sample(wt, length(wt), replace = TRUE)))
    ci <- stats::quantile(draws, c(.025, .975), names = FALSE, type = 7)
    boots[[length(boots) + 1L]] <- cbind(info, n_WT = length(wt), n_cKO = length(ko), bootstrap_iterations = cfg$bootstrap_iterations,
      primary_effect = fit$effect_cKO_minus_WT, bootstrap_CI_low = ci[1L], bootstrap_CI_high = ci[2L],
      fraction_positive = mean(draws > 0), fraction_negative = mean(draws < 0),
      resampling_unit = 'independent_mouse_with_replacement_within_each_genotype',
      interval_type = 'exploratory_percentile_95_percent_not_multiplicity_adjusted')
  }
  family <- if (variant == 'full_CollecTRI_prior') 'BH_all_tested_TFs_x_five_endpoints_full_prior' else 'BH_tested_Satb2_no_self_target_x_five_endpoints'
  combine <- function(x, empty) if (length(x)) do.call(rbind, x) else empty
  model_empty <- data.frame(endpoint = character(), TF = character(), omitted_mouse = character(), tested = logical(), effect_cKO_minus_WT = numeric(), pvalue = numeric())
  list(primary = adjust_family(do.call(rbind, primary), family),
    sensitivity = adjust_family(do.call(rbind, sensitive), paste0(family, '_same_partition_sensitivity')),
    loo = combine(loo, model_empty),
    stability = combine(stability, data.frame(endpoint = character(), TF = character(), n_loo_tested = integer(), loo_same_sign_fraction = numeric())),
    bootstrap = combine(boots, data.frame(endpoint = character(), TF = character(), bootstrap_iterations = integer(), bootstrap_CI_low = numeric(), bootstrap_CI_high = numeric())))
}

save_plot <- function(p, stem, width = 11, height = 7) {
  dir.create(dirname(stem), showWarnings = FALSE, recursive = TRUE)
  for (ext in c('png', 'pdf', 'svg')) {
    path <- paste0(stem, '.', ext); need(!file.exists(path), paste('Existing plot:', path))
    device <- switch(ext, png = 'png', pdf = grDevices::cairo_pdf, svg = grDevices::svg)
    ggplot2::ggsave(path, p, width = width, height = height, dpi = 180, device = device, bg = 'white', limitsize = FALSE)
  }
}
short_endpoint <- function(x) c(CA1_Deep = 'CA1 DL', CA1_Superficial = 'CA1 SL', IN_C1 = 'IN C1', IN_C6 = 'IN C6', IN_C8 = 'IN C8')[x]
rank_effects <- function(x) x[order(is.na(x$association_fdr), x$association_fdr, is.na(x$pvalue), x$pvalue, -abs(x$effect_cKO_minus_WT), x$TF), , drop = FALSE]

volcano <- function(effects, title, cfg) {
  x <- effects[effects$tested, , drop = FALSE]
  if (!nrow(x)) return(ggplot2::ggplot() + ggplot2::annotate('text', x = 0, y = 0, label = 'No eligible continuous-score tests') + ggplot2::labs(title = title) + ggplot2::theme_void())
  x$endpoint_label <- factor(short_endpoint(x$endpoint), levels = short_endpoint(ENDPOINTS))
  x$minus_log10_association_fdr <- -log10(pmax(x$association_fdr, .Machine$double.xmin))
  x$direction <- ifelse(x$association_fdr < .05, ifelse(x$effect_cKO_minus_WT > 0, 'Higher in cKO', 'Lower in cKO'), 'Below association threshold')
  lab <- do.call(rbind, lapply(split(x, x$endpoint), function(z) head(rank_effects(z), cfg$volcano_labels_per_endpoint)))
  p <- ggplot2::ggplot(x, ggplot2::aes(effect_cKO_minus_WT, minus_log10_association_fdr, color = direction)) +
    ggplot2::geom_point(size = 1.3, alpha = .85) + ggplot2::geom_hline(yintercept = -log10(.05), linetype = 2, color = '#aab5bb') +
    ggplot2::geom_vline(xintercept = 0, color = '#dce1e4') + ggplot2::facet_wrap(~endpoint_label, scales = 'free_x', ncol = 3, drop = FALSE) +
    ggplot2::scale_color_manual(values = c('Higher in cKO' = '#c63d4f', 'Lower in cKO' = '#397eb5', 'Below association threshold' = '#c3cbd0')) +
    ggplot2::labs(title = title, subtitle = 'Independent-mouse continuous ULM score association; these are not gene DEGs',
      x = 'Mean ULM score difference (cKO - WT)', y = '-log10 genotype-association BH FDR', color = NULL) + ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = 'bottom')
  p + ggrepel::geom_text_repel(data = lab, ggplot2::aes(label = TF), size = 2.6, max.overlaps = Inf, seed = cfg$seed, min.segment.length = Inf, show.legend = FALSE)
}

make_plots <- function(result, scores, cfg, outdir, no_self = NULL) {
  effects <- result$primary; tested <- effects[effects$tested, , drop = FALSE]
  score_selection <- list()
  if (nrow(tested)) {
    order <- rank_effects(tested); tfs <- head(unique(order$TF), cfg$heatmap_top_tfs)
    heat <- effects[effects$TF %in% tfs, , drop = FALSE]
    heat$TF <- factor(heat$TF, levels = rev(tfs)); heat$endpoint_label <- factor(short_endpoint(heat$endpoint), levels = short_endpoint(ENDPOINTS))
    heat$display_effect <- ifelse(heat$tested, heat$effect_cKO_minus_WT, NA_real_)
    limit <- max(abs(heat$display_effect), na.rm = TRUE); if (!is.finite(limit) || limit <= 0) limit <- 1
    heat$star <- ifelse(!is.na(heat$association_fdr) & heat$association_fdr < .05, '*', '')
    p <- ggplot2::ggplot(heat, ggplot2::aes(endpoint_label, TF, fill = display_effect)) + ggplot2::geom_tile(color = 'white') +
      ggplot2::geom_text(ggplot2::aes(label = star), size = 4) + ggplot2::scale_fill_gradient2(low = '#397eb5', mid = 'white', high = '#c63d4f', limits = c(-limit, limit), na.value = '#e4e8eb') +
      ggplot2::labs(title = 'CollecTRI TF score differences across the five populations',
        subtitle = 'Top TFs by global association FDR; * association FDR < 0.05; gray = not tested', x = NULL, y = NULL, fill = 'cKO - WT\nULM score') +
      ggplot2::theme_minimal(base_size = 11)
    save_plot(p, file.path(outdir, 'figures/TF_effect_heatmap'), width = 9, height = max(5, length(tfs) * .25))
    score_selection[[1L]] <- data.frame(plot = 'TF_effect_heatmap', endpoint = 'all_five', TF = tfs, selection = 'minimum_global_association_FDR_then_P_then_absolute_effect')
  } else save_plot(ggplot2::ggplot() + ggplot2::annotate('text', 0, 0, label = 'No eligible continuous-score tests') + ggplot2::theme_void(), file.path(outdir, 'figures/TF_effect_heatmap'))
  save_plot(volcano(effects, 'TF activity: full-prior mouse associations', cfg), file.path(outdir, 'figures/TF_association_volcano'), width = 12, height = 8)
  save_plot(volcano(result$sensitivity, 'TF activity: same-partition sensitivity', cfg), file.path(outdir, 'figures/TF_same_partition_volcano'), width = 12, height = 8)
  satb2_status <- list()
  for (ep in ENDPOINTS) {
    ep_effect <- effects[effects$endpoint == ep, , drop = FALSE]
    top <- head(rank_effects(ep_effect[ep_effect$tested & ep_effect$TF != cfg$satb2_tf, , drop = FALSE])$TF, cfg$boxplot_top_tfs)
    selection <- unique(c(cfg$satb2_tf, top)); x <- scores[scores$endpoint == ep & scores$TF %in% selection, , drop = FALSE]
    x$TF <- factor(x$TF, levels = selection); x$condition <- factor(x$condition, levels = c('WT', 'cKO'))
    missing <- selection[!selection %in% as.character(x$TF)]
    labels <- data.frame(TF = factor(selection, levels = selection), label = vapply(selection, function(tf) {
      e <- ep_effect[ep_effect$TF == tf, , drop = FALSE]
      if (!nrow(e)) return('Unavailable in curated prior')
      if (!e$tested) return(e$model_status)
      paste0('Difference=', signif(e$effect_cKO_minus_WT, 3), '; association q=', format(e$association_fdr, digits = 2, scientific = TRUE))
    }, character(1)))
    p <- ggplot2::ggplot(x, ggplot2::aes(condition, escore, fill = condition)) +
      ggplot2::geom_boxplot(width = .55, outlier.shape = NA, alpha = .25) +
      ggplot2::geom_point(ggplot2::aes(color = condition), position = ggplot2::position_jitter(width = .10, height = 0, seed = cfg$seed), size = 2) +
      ggplot2::facet_wrap(~TF, ncol = 3, scales = 'free_y', drop = FALSE) +
      ggplot2::scale_fill_manual(values = c(WT = '#397eb5', cKO = '#cf8334'), drop = FALSE) +
      ggplot2::scale_color_manual(values = c(WT = '#397eb5', cKO = '#cf8334'), drop = FALSE) +
      ggplot2::geom_text(data = labels, ggplot2::aes(x = 1.5, y = Inf, label = label), inherit.aes = FALSE, vjust = 1.3, size = 2.45) +
      ggplot2::labs(title = paste(short_endpoint(ep), 'individual mouse TF scores'),
        subtitle = 'Satb2 plus five selected TFs; six WT and six cKO mice; post hoc top TF selection is exploratory',
        x = NULL, y = 'Native CollecTRI ULM t-value') + ggplot2::theme_minimal(base_size = 11) + ggplot2::theme(legend.position = 'none')
    if (length(missing)) p <- p + ggplot2::geom_text(data = data.frame(TF = factor(missing, levels = selection), message = 'No modeled activity\nNo RNA proxy substituted'),
      ggplot2::aes(x = 1.5, y = 0, label = message), inherit.aes = FALSE, size = 3)
    save_plot(p, file.path(outdir, 'figures', paste0(ep, '_Satb2_top5_mouse_scores')), width = 12, height = 7)
    score_selection[[length(score_selection) + 1L]] <- data.frame(plot = paste0(ep, '_Satb2_top5_mouse_scores'), endpoint = ep, TF = selection,
      selection = ifelse(selection == cfg$satb2_tf, 'prespecified_Satb2', 'post_hoc_global_association_FDR_then_P_then_absolute_effect'))
    row <- ep_effect[ep_effect$TF == cfg$satb2_tf, , drop = FALSE]
    satb2_status[[ep]] <- data.frame(endpoint = ep, TF = cfg$satb2_tf, present_in_prior = nrow(row) > 0,
      activity_modeled = any(scores$endpoint == ep & scores$TF == cfg$satb2_tf),
      association_tested = nrow(row) > 0 && row$tested, status = if (nrow(row)) row$model_status else 'TF_ABSENT_FROM_COLLECTRI_PRIOR')
  }
  if (!is.null(no_self)) {
    joined <- rbind(transform(effects[effects$TF == cfg$satb2_tf, , drop = FALSE], panel = 'Full curated prior'),
      transform(no_self$primary, panel = 'Satb2 self-target excluded'))
    p <- ggplot2::ggplot(joined[joined$tested, , drop = FALSE], ggplot2::aes(effect_cKO_minus_WT, short_endpoint(endpoint), color = panel)) +
      ggplot2::geom_vline(xintercept = 0, color = '#c3cbd0') + ggplot2::geom_errorbar(ggplot2::aes(xmin = CI_low, xmax = CI_high), orientation = 'y', position = ggplot2::position_dodge(width = .4), width = .18) +
      ggplot2::geom_point(position = ggplot2::position_dodge(width = .4), size = 2.3) + ggplot2::labs(title = 'Satb2 score association: self-target sensitivity',
        subtitle = 'Same eligible RNA background; curated target-set activity is not independent proof of SATB2 binding',
        x = 'cKO - WT ULM score difference (mouse-level 95% CI)', y = NULL, color = NULL) + ggplot2::theme_minimal(base_size = 11) + ggplot2::theme(legend.position = 'bottom')
    save_plot(p, file.path(outdir, 'figures/Satb2_self_target_sensitivity'), width = 11, height = 5)
  }
  write_tsv(do.call(rbind, score_selection), file.path(outdir, 'TF_PLOT_SELECTION.tsv'))
  write_tsv(do.call(rbind, satb2_status), file.path(outdir, 'SATB2_STATUS.tsv'))
}

write_variant <- function(result, prefix, outdir) {
  write_tsv(result$primary, file.path(outdir, paste0(prefix, 'ASSOCIATIONS.tsv')))
  write_tsv(result$sensitivity, file.path(outdir, paste0(prefix, 'SAME_PARTITION_ASSOCIATIONS.tsv')))
  write_tsv(result$loo, file.path(outdir, paste0(prefix, 'LOO_MODELS.tsv')))
  write_tsv(result$stability, file.path(outdir, paste0(prefix, 'LOO_STABILITY.tsv')))
  write_tsv(result$bootstrap, file.path(outdir, paste0(prefix, 'MOUSE_BOOTSTRAP.tsv')))
  write_tsv(result$primary[, intersect(c('endpoint', 'endpoint_id', 'source_endpoint_id', 'TF', 'network_variant', 'n_used_targets', 'n_background_genes', 'n_scored_mice', 'n_expected_mice', 'n_WT', 'n_cKO', 'tested', 'model_status'), names(result$primary)), drop = FALSE],
    file.path(outdir, paste0(prefix, 'MODEL_ELIGIBILITY.tsv')))
}

main <- function() {
  options(digits = 17)
  args <- args_parse(commandArgs(trailingOnly = TRUE))
  for (pkg in c('jsonlite', 'digest', 'ggplot2', 'ggrepel')) need(requireNamespace(pkg, quietly = TRUE), paste('Missing managed R dependency:', pkg))
  config_receipt <- receipt(args$config); cfg <- jsonlite::read_json(args$config, simplifyVector = FALSE)
  need(identical(cfg$schema, 'satb2-tf-activity/v1'), 'Unexpected TF config schema')
  base_policy <- policy_defaults(); base_policy$seed <- cfg$random_seed %||% cfg$seed %||% base_policy$seed
  provided_policy <- cfg$mouse_associations %||% list()
  need(all(names(provided_policy) %in% names(base_policy)), paste('Unknown mouse-association policy:', paste(setdiff(names(provided_policy), names(base_policy)), collapse = ', ')))
  policy <- utils::modifyList(base_policy, provided_policy)
  need(policy$bootstrap_iterations == 200L && policy$expected_primary_WT == 6L && policy$expected_primary_cKO == 6L,
    'Protocol fixes 200 independent-mouse bootstraps and all12 primary mice')
  need(policy$expected_partition_WT == 6L && policy$expected_partition_cKO == 3L && identical(policy$same_partition, '23HJT3LT4_lanes6_7'), 'Protocol fixes the known same-partition 6WT/3cKO sensitivity')
  need(policy$heatmap_top_tfs >= 1L && policy$heatmap_top_tfs <= 60L && policy$boxplot_top_tfs == 5L, 'Invalid bounded plot policy')
  output <- args$outdir
  need(!file.exists(output) || (dir.exists(output) && !nzchar(Sys.readlink(output)) && !length(list.files(output, all.files = TRUE, no.. = TRUE))),
    'Output must be nonexistent or strictly empty; never reuse a nonempty directory')
  native <- normalizePath(args[['input-dir']], mustWork = TRUE)
  completion_path <- file.path(native, 'COMPLETE.json'); native_complete <- jsonlite::read_json(completion_path, simplifyVector = FALSE)
  need(identical(native_complete$schema, 'satb2-tf-activity-complete/v1') && identical(native_complete$status, 'COMPLETE') &&
    identical(native_complete$config$sha256, config_receipt$sha256) && identical(native_complete$analysis_id, cfg$analysis_id), 'Native activity step incomplete or config/analysis differs')
  required_native <- c('TF_ACTIVITY_PER_MOUSE.tsv', 'TF_TARGET_COVERAGE.tsv', 'PROVENANCE.json')
  if (isTRUE(policy$require_satb2_no_self)) required_native <- c(required_native, 'SATB2_NO_SELF_TARGET_ACTIVITY_PER_MOUSE.tsv', 'SATB2_NO_SELF_TARGET_COVERAGE.tsv')
  need(all(required_native %in% names(native_complete$output_files)), 'Required native tables absent from completed producer receipt')
  inputs <- list(native_complete = receipt(completion_path), native_provenance = receipt(file.path(native, 'PROVENANCE.json')), samples = check_record(cfg$samples, 'samples'))
  # Verify the producer's complete native evidence; never accept partial score files.
  for (rel in names(native_complete$output_files)) {
    record <- native_complete$output_files[[rel]]
    need(identical(normalizePath(record$path, mustWork = TRUE), normalizePath(file.path(native, rel), mustWork = TRUE)), 'Native output receipt points outside its declared table')
    inputs[[paste0('native/', rel)]] <- check_record(record, paste('native', rel))
  }
  for (name in names(cfg$additional_hash_inputs %||% list())) inputs[[paste0('additional/', name)]] <- check_record(cfg$additional_hash_inputs[[name]], name)
  source_path <- sub('^--file=', '', grep('^--file=', commandArgs(), value = TRUE)[1L]); source_receipt <- receipt(source_path)
  roster <- read_tsv(inputs$samples$path, c('sample', 'condition'))
  need(nrow(roster) == 12L && !anyDuplicated(roster$sample) && sum(roster$condition == 'WT') == 6L && sum(roster$condition == 'cKO') == 6L, 'Expected exactly12 independent mice, six of each genotype')
  partition <- sample_partition(cfg, roster); roster$sequencing_partition <- partition$partition
  inputs$partition_metadata <- partition$source
  required <- c('endpoint', 'endpoint_id', 'source_endpoint_id', 'mode', 'TF', 'sample', 'escore', 'gene_set_fdr', 'gene_set_fdr_family',
    'n_used_targets', 'n_background_genes', 'score_unit', 'condition', 'n_nuclei', 'complete_gene_UMIs', 'eligible')
  scores <- read_tsv(file.path(native, 'TF_ACTIVITY_PER_MOUSE.tsv'), required)
  coverage <- read_tsv(file.path(native, 'TF_TARGET_COVERAGE.tsv'), c('endpoint', 'endpoint_id', 'source_endpoint_id', 'mode', 'TF', 'eligible_tmin5', 'n_used_targets', 'n_background_genes', 'reason'))
  coverage <- coverage[coverage$mode == PRIMARY_MODE, , drop = FALSE]
  validate_native_scores(scores, coverage, roster, PRIMARY_MODE)
  scores$sequencing_partition <- roster$sequencing_partition[match(scores$sample, roster$sample)]
  no_scores <- no_cov <- NULL
  if (isTRUE(policy$require_satb2_no_self) || file.exists(file.path(native, 'SATB2_NO_SELF_TARGET_ACTIVITY_PER_MOUSE.tsv'))) {
    no_scores <- read_tsv(file.path(native, 'SATB2_NO_SELF_TARGET_ACTIVITY_PER_MOUSE.tsv'), c(required, 'network_variant'))
    no_cov <- read_tsv(file.path(native, 'SATB2_NO_SELF_TARGET_COVERAGE.tsv'), c('endpoint', 'endpoint_id', 'source_endpoint_id', 'mode', 'TF', 'eligible_tmin5', 'n_used_targets', 'n_background_genes', 'reason'))
    no_cov <- no_cov[no_cov$mode == NOSELF_MODE, , drop = FALSE]
    need(all(no_scores$TF == policy$satb2_tf) && all(no_cov$TF == policy$satb2_tf), 'No-self sensitivity must contain only the prespecified Satb2 TF')
    validate_native_scores(no_scores, no_cov, roster, NOSELF_MODE)
    no_scores$sequencing_partition <- roster$sequencing_partition[match(no_scores$sample, roster$sample)]
  }
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  tryCatch({
    set.seed(policy$seed)
    result <- analyze_variant(scores, coverage, roster, policy, 'full_CollecTRI_prior')
    no_self <- if (!is.null(no_scores)) analyze_variant(no_scores, no_cov, roster, policy, 'SATB2_only_self_target_excluded') else NULL
    write_variant(result, 'TF_', output)
    if (!is.null(no_self)) write_variant(no_self, 'SATB2_NO_SELF_TARGET_', output)
    write_tsv(scores, file.path(output, 'TF_MOUSE_SCORES.tsv'))
    if (!is.null(no_scores)) write_tsv(no_scores, file.path(output, 'SATB2_NO_SELF_TARGET_MOUSE_SCORES.tsv'))
    write_tsv(roster[, c('sample', 'condition', 'sequencing_partition')], file.path(output, 'MOUSE_PARTITION_AUDIT.tsv'))
    if (isTRUE(policy$plots)) make_plots(result, scores, policy, output, no_self)
    capture.output(sessionInfo(), file = file.path(output, 'sessionInfo.txt'))
    for (name in names(inputs)) need(identical(receipt(inputs[[name]]$path), inputs[[name]]), paste('Input changed during association analysis:', name))
    need(identical(receipt(args$config), config_receipt) && identical(receipt(source_path), source_receipt), 'Config or source changed during analysis')
    model_counts <- data.frame(endpoint = ENDPOINTS,
      n_tested = vapply(ENDPOINTS, function(ep) sum(result$primary$endpoint == ep & result$primary$tested), integer(1)),
      n_association_fdr_lt05 = vapply(ENDPOINTS, function(ep) sum(result$primary$endpoint == ep & !is.na(result$primary$association_fdr) & result$primary$association_fdr < .05), integer(1)))
    write_tsv(model_counts, file.path(output, 'MODEL_SUMMARY.tsv'))
    json_write(list(schema = 'satb2-tf-mouse-provenance/v1', analysis_id = cfg$analysis_id, config = config_receipt, source = source_receipt, inputs = inputs, policy = policy,
      model = 'native_ULM_t_score ~ condition; unweighted continuous mouse-level linear model',
      primary_FDR = 'BH jointly across every tested TF x five endpoint primary score association',
      sensitivity_FDR = 'separate BH for same-partition sensitivity; separate families for Satb2 no-self-target',
      native_gene_set_fdr = 'Preserved within-mouse prior-target enrichment FDR; never used to filter or label genotype associations',
      constant_outcomes = 'No tests for globally constant or numerically zero within-genotype residual variance; effect means may remain descriptive',
      resampling_unit = 'independent_mouse_within_genotype; 200 percentile bootstraps; leave-one-mouse-out summaries are sensitivity checks',
      confounding = 'Genotype partly confounds technical partition. All associations and selected plots are exploratory, not causally separated technical effects.',
      interpretation = 'Signed curated-prior target-set activity is not direct binding or independent proof of SATB2 activity; self-target exclusion is a sensitivity, not validation.',
      original_gene_DESeq2_fits_changed = FALSE, new_gene_DE_analysis = FALSE, motif_regulon_inference = FALSE,
      software = lapply(c('ggplot2', 'ggrepel', 'jsonlite', 'digest'), function(p) list(package = p, version = as.character(utils::packageVersion(p))))), file.path(output, 'PROVENANCE.json'))
    paths <- list.files(output, recursive = TRUE, full.names = TRUE); paths <- paths[!dir.exists(paths)]
    receipts <- lapply(paths, function(p) {r <- receipt(p); r$path <- NULL; r}); names(receipts) <- substring(paths, nchar(output) + 2L)
    json_write(list(schema = 'satb2-tf-mouse-complete/v1', status = 'COMPLETE', analysis_id = cfg$analysis_id,
      finished_UTC = format(Sys.time(), tz = 'UTC', usetz = TRUE), config = config_receipt, source = source_receipt, inputs = inputs, outputs = receipts,
      checks = list(all12_mouse_roster = TRUE, original_native_scores_preserved = TRUE, global_primary_BH = TRUE,
        separate_sensitivity_BH = TRUE, native_gene_set_FDR_not_association_FDR = TRUE, no_nucleus_pseudoreplication = TRUE,
        no_new_gene_DESeq2 = TRUE, inputs_config_source_SHA256_unchanged = TRUE)), file.path(output, 'COMPLETE.json'))
  }, error = function(e) {
    if (!file.exists(file.path(output, 'FAILED.json'))) json_write(list(status = 'FAILED', message = conditionMessage(e)), file.path(output, 'FAILED.json'))
    stop(e)
  })
}

if (sys.nframe() == 0L) main()
