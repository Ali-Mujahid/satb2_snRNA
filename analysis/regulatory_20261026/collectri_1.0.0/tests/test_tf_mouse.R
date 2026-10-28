# Synthetic fixtures only. This file never reads study scores or writes plots.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
env <- new.env(); sys.source(args[1L], envir = env)
cfg <- env$policy_defaults()
roster <- data.frame(sample = paste0('fixture_', seq_len(12L)), condition = rep(c('WT', 'cKO'), each = 6L),
  sequencing_partition = c(rep(cfg$same_partition, 9L), rep('other', 3L)))
x <- data.frame(sample = roster$sample, condition = roster$condition, escore = c(1:6, 3:8))
fit <- env$model_score(x, cfg, 6L, 6L)
stopifnot(fit$tested, fit$model_status == 'MOUSE_LEVEL_CONTINUOUS_LM', abs(fit$effect_cKO_minus_WT - 2) < 1e-12, fit$residual_df == 10L)
x$escore <- rep(1, 12L)
stopifnot(!env$model_score(x, cfg)$tested, env$model_score(x, cfg)$model_status == 'CONSTANT_OR_NUMERICALLY_CONSTANT_SCORE')
x$escore <- rep(c(0, 1), each = 6L)
stopifnot(!env$model_score(x, cfg)$tested, env$model_score(x, cfg)$model_status == 'ZERO_OR_NUMERICALLY_ZERO_WITHIN_GROUP_VARIANCE')
rows <- lapply(seq_along(env$ENDPOINTS), function(i) data.frame(endpoint = env$ENDPOINTS[i], endpoint_id = paste0('fixture_', i), source_endpoint_id = paste0('fixture_source_', i),
  mode = env$PRIMARY_MODE, TF = 'Satb2', sample = roster$sample, escore = c(1:6, 1:6 + (i - 1) * .4), gene_set_fdr = .99,
  gene_set_fdr_family = 'native_fixture_only', n_used_targets = 5L, n_background_genes = 100L,
  score_unit = 'native_ULM_t_value', condition = roster$condition, n_nuclei = 10L, complete_gene_UMIs = 100L,
  eligible = TRUE, sequencing_partition = roster$sequencing_partition))
scores <- do.call(rbind, rows)
cov <- data.frame(endpoint = env$ENDPOINTS, endpoint_id = paste0('fixture_', 1:5), source_endpoint_id = paste0('fixture_source_', 1:5),
  mode = env$PRIMARY_MODE, TF = 'Satb2', eligible_tmin5 = TRUE, n_used_targets = 5L, n_background_genes = 100L, reason = 'MODELED')
env$validate_native_scores(scores, cov, roster, env$PRIMARY_MODE)
set.seed(cfg$seed)
result <- env$analyze_variant(scores, cov, roster, cfg, 'full_CollecTRI_prior')
stopifnot(nrow(result$primary) == 5L, all(result$primary$tested), all(result$primary$n_WT == 6L), all(result$primary$n_cKO == 6L),
  all(result$primary$n_scored_mice == 12L), all(result$sensitivity$n_WT == 6L), all(result$sensitivity$n_cKO == 3L),
  nrow(result$loo) == 60L, all(result$stability$n_loo_models == 12L), all(result$bootstrap$bootstrap_iterations == 200L),
  identical(result$primary$association_fdr, stats::p.adjust(result$primary$pvalue, 'BH')),
  identical(result$primary$association_fdr, result$primary$padj), all(!result$primary$gene_set_fdr_used_for_filtering))
# Native target-set FDR does not gate or alter the continuous mouse contrast.
y <- scores[scores$endpoint == env$ENDPOINTS[2], ]; a <- env$model_score(y, cfg)
y$gene_set_fdr <- .00001; b <- env$model_score(y, cfg)
stopifnot(identical(a, b))
bad <- scores[-1L, ]
stopifnot(inherits(tryCatch(env$validate_native_scores(bad, cov, roster, env$PRIMARY_MODE), error = identity), 'error'))
# A no-self endpoint with insufficient targets is explicitly unavailable, not tested.
no_cov <- cov; no_cov$mode <- env$NOSELF_MODE; no_cov$eligible_tmin5[5] <- FALSE; no_cov$n_used_targets[5] <- 0L; no_cov$reason[5] <- 'FEWER_THAN_5_ELIGIBLE_TARGETS'
no_scores <- scores[scores$endpoint != env$ENDPOINTS[5], ]; no_scores$mode <- env$NOSELF_MODE
env$validate_native_scores(no_scores, no_cov, roster, env$NOSELF_MODE)
set.seed(cfg$seed)
no_result <- env$analyze_variant(no_scores, no_cov, roster, cfg, 'SATB2_only_self_target_excluded')
stopifnot(sum(no_result$primary$tested) == 4L, !no_result$primary$tested[5L], no_result$primary$n_scored_mice[5L] == 0L,
  is.na(no_result$primary$association_fdr[5L]), all(grepl('Satb2_no_self_target', no_result$primary$fdr_family)))
cat('Synthetic continuous-score, coverage, BH-family, LOO, bootstrap and unavailable-target tests passed. No study data read.\n')
