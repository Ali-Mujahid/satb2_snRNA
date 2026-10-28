#!/usr/bin/env Rscript
# Positive motif-supported regulon AUCell: independent-mouse continuous associations.
# Scientific execution only via the immutable native workflow. Sourcing is inert.
# No count-model/DESeq2/edgeR, motif inference, nucleus pseudo-replication, online requests.

`%||%` <- function(x, y) if (is.null(x)) y else x
ENDPOINTS <- c('CA1_Deep', 'CA1_Superficial', 'IN_C1', 'IN_C6', 'IN_C8')

policy_defaults <- function() list(
  seed = 1337L, bootstrap_iterations = 200L,
  same_partition = '23HJT3LT4_lanes6_7',
  expected_primary_WT = 6L, expected_primary_cKO = 6L,
  expected_partition_WT = 6L, expected_partition_cKO = 3L,
  variation_relative_tolerance = 1e-12, residual_relative_tolerance = 1e-12,
  heatmap_top_regulons = 30L, boxplot_top_regulons = 6L, volcano_labels_per_endpoint = 3L, plots = TRUE
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
  x$association_fdr <- rep(NA_real_,nrow(x))
  ix <- which(x$tested & is.finite(x$pvalue))
  if (length(ix)) x$association_fdr[ix] <- stats::p.adjust(x$pvalue[ix], method = 'BH')
  x$padj <- x$association_fdr; x$fdr_family <- rep(family,nrow(x))
  x$activity_threshold_used_for_filtering <- rep(FALSE,nrow(x))
  x
}


lineage_for <- function(endpoint) ifelse(startsWith(as.character(endpoint), 'CA1_'), 'CA1', 'IN')
endpoint_id <- function(endpoint) unname(c(CA1_Deep='CA1_DL', CA1_Superficial='CA1_SL', IN_C1='IN_C1', IN_C6='IN_C6', IN_C8='IN_C8')[endpoint])
source_endpoint_id <- function(endpoint) unname(c(CA1_Deep='P0244_seurat_CA1_ProS_LAYER_PROGRAM', CA1_Superficial='P0247_seurat_CA1_ProS_LAYER_PROGRAM',
  IN_C1='P0140_seurat_Inhibitory_CLUSTER', IN_C6='P0145_seurat_Inhibitory_CLUSTER', IN_C8='P0147_seurat_Inhibitory_CLUSTER')[endpoint])
short_endpoint <- function(x) unname(c(CA1_Deep='CA1 DL', CA1_Superficial='CA1 SL', IN_C1='IN C1', IN_C6='IN C6', IN_C8='IN C8')[x])

validate_native_scores <- function(scores, coverage, roster, lineage_status, supplement=FALSE) {
  need(setequal(lineage_status$lineage, c('CA1','IN')) && nrow(lineage_status)==2L && !anyDuplicated(lineage_status$lineage), 'Need two unambiguous lineage statuses')
  need(all(lineage_status$status %in% c('COMPLETE_WITH_REGULONS','EMPTY_NO_SUPPORTED_REGULONS')), 'Native lineage is not complete')
  need(all(scores$endpoint %in% ENDPOINTS) && all(coverage$endpoint %in% ENDPOINTS), 'Unknown endpoint')
  need(!anyDuplicated(scores[c('endpoint','regulon_id','sample')]) && !anyDuplicated(coverage[c('endpoint','regulon_id')]), 'Duplicate regulon keys')
  need(all(scores$sample %in% roster$sample), 'Unknown mouse in native score table')
  need(all(scores$condition == roster$condition[match(scores$sample,roster$sample)]), 'Native mouse genotype changed')
  for (x in list(scores,coverage)) {
    need(all(x$lineage == lineage_for(x$endpoint)) && all(startsWith(as.character(x$regulon_id),paste0(x$lineage,'::'))), 'Regulon lineage/ID mismatch')
    need(all(x$endpoint_id == endpoint_id(x$endpoint)) && all(x$source_endpoint_id == source_endpoint_id(x$endpoint)), 'Native endpoint identity mismatch')
    need(all(nzchar(x$TF)) && all(nzchar(x$native_name)) && all(x$n_targets>=0) && all(x$n_background_genes>=x$n_targets), 'Missing native regulon identity/coverage')
    need(all(grepl('\\(\\+\\)$',x$native_name)), 'Only native positive motif-supported regulons are permitted')
  }
  need(all(is.finite(scores$aucell_mean)) && all(scores$aucell_mean>=0 & scores$aucell_mean<=1), 'Invalid native AUCell AUC')
  need(all(scores$score_unit=='native_AUCell_AUC') && all(as_flag(scores$eligible)), 'Invalid native score eligibility/unit')
  need(all(scores$n_nuclei>0) && all(scores$n_nuclei==as.integer(scores$n_nuclei)), 'Invalid nucleus coverage')
  if(!supplement)need(all(coverage$n_targets>0),'Native discovered regulons must retain positive target count')
  if(!supplement)for (lin in c('CA1','IN')) {
    status <- lineage_status[lineage_status$lineage==lin,,drop=FALSE]
    c <- coverage[coverage$lineage==lin,,drop=FALSE]; s <- scores[scores$lineage==lin,,drop=FALSE]
    regs <- unique(c$regulon_id); eps <- ENDPOINTS[lineage_for(ENDPOINTS)==lin]
    need(length(regs)==status$n_regulons, 'Regulon count differs from completed native status')
    if (status$status=='EMPTY_NO_SUPPORTED_REGULONS') {need(!length(regs) && !nrow(s), 'Empty native lineage has manufactured activities'); next}
    need(length(regs)>0 && setequal(c$endpoint,eps) && nrow(c)==length(regs)*length(eps), 'Each discovered lineage regulon needs each matching endpoint')
  }
  if(supplement)need(nrow(coverage)==5L && setequal(coverage$endpoint,ENDPOINTS) && all(coverage$TF=='Satb2') && all(coverage$regulon_id==paste0(coverage$lineage,'::Satb2(+)_NO_SELF_TARGET')),'No-self coverage must explicitly identify all five exact Satb2 endpoint states')
  for (i in seq_len(nrow(coverage))) {
    c <- coverage[i,,drop=FALSE]; x <- scores[scores$endpoint==c$endpoint & scores$regulon_id==c$regulon_id,,drop=FALSE]
    if (as_flag(c$eligible)) {
      need(c$n_targets>=10L,'Eligible regulon has fewer than10 measured targets')
      need(nrow(x)==12L && setequal(x$sample,roster$sample),'Eligible regulon lacks all twelve independent mice')
      for (k in c('lineage','native_name','TF','n_targets','n_background_genes')) need(length(unique(x[[k]]))==1L && x[[k]][1L]==c[[k]],paste('Native score/coverage mismatch:',k))
    } else need(nrow(x)==0L,'Unavailable regulon has manufactured activity')
  }
  need(all(paste(scores$endpoint,scores$regulon_id) %in% paste(coverage$endpoint,coverage$regulon_id)), 'Scored regulon absent from coverage')
  invisible(TRUE)
}

analyze_regulons <- function(scores,coverage,roster,cfg,family='BH_all_tested_positive_motif_regulons_x_matching_five_endpoints_both_lineages') {
  primary <- sensitive <- loo <- boots <- stability <- list()
  for (i in seq_len(nrow(coverage))) {
    cov <- coverage[i,,drop=FALSE]
    info <- cov[,c('lineage','endpoint','endpoint_id','source_endpoint_id','regulon_id','native_name','TF','n_targets','n_background_genes'),drop=FALSE]
    info$score_unit <- 'native_AUCell_AUC'; info$n_expected_mice <- 12L
    x <- scores[scores$endpoint==cov$endpoint & scores$regulon_id==cov$regulon_id,,drop=FALSE]
    x$escore <- x$aucell_mean; info$n_scored_mice <- nrow(x)
    if (!nrow(x)) {
      missing <- model_score(data.frame(sample=character(),condition=character(),escore=numeric()),cfg)
      missing$model_status <- paste0('NOT_MODELED_',cov$reason)
      primary[[i]] <- cbind(info,missing); sensitive[[i]] <- cbind(info,partition=cfg$same_partition,missing); next
    }
    fit <- model_score(x,cfg,cfg$expected_primary_WT,cfg$expected_primary_cKO)
    primary[[i]] <- cbind(info,fit)
    same <- x[x$sequencing_partition==cfg$same_partition,,drop=FALSE]
    sensitive[[i]] <- cbind(info,partition=cfg$same_partition,model_score(same,cfg,cfg$expected_partition_WT,cfg$expected_partition_cKO))
    if (!fit$tested) next
    cur <- do.call(rbind,lapply(roster$sample,function(s) cbind(info,omitted_mouse=s,model_score(x[x$sample!=s,,drop=FALSE],cfg))))
    loo[[length(loo)+1L]] <- cur; valid <- cur$tested & is.finite(cur$effect_cKO_minus_WT)
    stability[[length(stability)+1L]] <- cbind(info,primary_effect=fit$effect_cKO_minus_WT,n_loo_models=nrow(cur),n_loo_tested=sum(valid),
      n_loo_same_sign=sum(sign(cur$effect_cKO_minus_WT[valid])==sign(fit$effect_cKO_minus_WT)),
      loo_same_sign_fraction=if(any(valid))mean(sign(cur$effect_cKO_minus_WT[valid])==sign(fit$effect_cKO_minus_WT)) else NA_real_,
      loo_effect_min=if(any(valid))min(cur$effect_cKO_minus_WT[valid]) else NA_real_,loo_effect_max=if(any(valid))max(cur$effect_cKO_minus_WT[valid]) else NA_real_)
    wt<-x$escore[x$condition=='WT'];ko<-x$escore[x$condition=='cKO']
    draws<-replicate(cfg$bootstrap_iterations,mean(sample(ko,length(ko),replace=TRUE))-mean(sample(wt,length(wt),replace=TRUE)))
    ci<-stats::quantile(draws,c(.025,.975),names=FALSE,type=7)
    boots[[length(boots)+1L]]<-cbind(info,n_WT=length(wt),n_cKO=length(ko),bootstrap_iterations=cfg$bootstrap_iterations,
      primary_effect=fit$effect_cKO_minus_WT,bootstrap_CI_low=ci[1L],bootstrap_CI_high=ci[2L],fraction_positive=mean(draws>0),fraction_negative=mean(draws<0),
      resampling_unit='independent_mouse_with_replacement_within_each_genotype',interval_type='exploratory_percentile_95_percent_not_multiplicity_adjusted')
  }
  model_empty<-cbind(coverage[FALSE,c('lineage','endpoint','endpoint_id','source_endpoint_id','regulon_id','native_name','TF','n_targets','n_background_genes'),drop=FALSE],
    score_unit=character(),n_expected_mice=integer(),n_scored_mice=integer(),model_score(data.frame(sample=character(),condition=character(),escore=numeric()),cfg)[FALSE,,drop=FALSE])
  combine<-function(x,empty) if(length(x))do.call(rbind,x) else empty
  list(primary=adjust_family(combine(primary,model_empty),family),sensitivity=adjust_family(combine(sensitive,cbind(model_empty,partition=character())),paste0(family,'_same_partition_sensitivity')),
    loo=combine(loo,cbind(model_empty,omitted_mouse=character())),
    stability=combine(stability,data.frame(lineage=character(),endpoint=character(),regulon_id=character(),n_loo_tested=integer(),loo_same_sign_fraction=numeric())),
    bootstrap=combine(boots,data.frame(lineage=character(),endpoint=character(),regulon_id=character(),bootstrap_iterations=integer(),bootstrap_CI_low=numeric(),bootstrap_CI_high=numeric())))
}

save_plot <- function(p,stem,width=11,height=7) {
  dir.create(dirname(stem),showWarnings=FALSE,recursive=TRUE)
  for(ext in c('pdf','svg','png')) {
    path<-paste0(stem,'.',ext);need(!file.exists(path),paste('Refusing plot overwrite',path))
    ggplot2::ggsave(path,p,width=width,height=height,dpi=180,bg='white',limitsize=FALSE,
      device=switch(ext,pdf=grDevices::cairo_pdf,svg=grDevices::svg,png='png'))
  }
}
rank_effects <- function(x)x[order(is.na(x$association_fdr),x$association_fdr,is.na(x$pvalue),x$pvalue,-abs(x$effect_cKO_minus_WT),x$regulon_id),,drop=FALSE]
empty_plot <- function(label)ggplot2::ggplot()+ggplot2::annotate('text',0,0,label=label)+ggplot2::theme_void()
make_plots <- function(result,scores,cfg,outdir,no_self=NULL) {
  effects<-result$primary;tested<-effects[effects$tested,,drop=FALSE]; selections<-list()
  order<-rank_effects(tested); regs<-head(unique(order$regulon_id),cfg$heatmap_top_regulons)
  if(length(regs)) {
    heat<-effects[effects$regulon_id %in% regs,,drop=FALSE];heat$regulon_id<-factor(heat$regulon_id,levels=rev(regs));heat$endpoint_label<-factor(short_endpoint(heat$endpoint),levels=short_endpoint(ENDPOINTS))
    heat$display_effect<-ifelse(heat$tested,heat$effect_cKO_minus_WT,NA_real_);lim<-max(abs(heat$display_effect),na.rm=TRUE);if(!is.finite(lim)||lim<=0)lim<-1
    heat$star<-ifelse(!is.na(heat$association_fdr)&heat$association_fdr<.05,'*','')
    plot<-ggplot2::ggplot(heat,ggplot2::aes(endpoint_label,regulon_id,fill=display_effect))+ggplot2::geom_tile(color='white')+ggplot2::geom_text(ggplot2::aes(label=star),size=4)+
      ggplot2::scale_fill_gradient2(low='#397eb5',mid='white',high='#c63d4f',limits=c(-lim,lim),na.value='#e4e8eb')+
      ggplot2::labs(title='Motif-supported regulon AUCell differences',subtitle='Post hoc top regulons by global association FDR; * q<0.05; lineage-specific definitions',x=NULL,y=NULL,fill='cKO - WT\nmean AUC')+ggplot2::theme_minimal(base_size=10)
    selections[[1L]]<-data.frame(plot='REGULON_effect_heatmap',endpoint='matching_lineage',regulon_id=regs,selection='post_hoc_global_association_FDR_then_P_then_absolute_effect')
  } else plot<-empty_plot('No eligible nonconstant regulon association tests')
  save_plot(plot,file.path(outdir,'figures/REGULON_effect_heatmap'),9,max(5,length(regs)*.25))
  for(variant in c('primary','sensitivity')) {
    x<-result[[variant]];x<-x[x$tested,,drop=FALSE]
    if(nrow(x)) {
      x$endpoint_label<-factor(short_endpoint(x$endpoint),levels=short_endpoint(ENDPOINTS));x$y<--log10(pmax(x$association_fdr,.Machine$double.xmin));x$direction<-ifelse(x$association_fdr<.05,ifelse(x$effect_cKO_minus_WT>0,'Higher in cKO','Lower in cKO'),'Below threshold')
      lab<-do.call(rbind,lapply(split(x,x$endpoint),function(z)head(rank_effects(z),cfg$volcano_labels_per_endpoint)))
      plot<-ggplot2::ggplot(x,ggplot2::aes(effect_cKO_minus_WT,y,color=direction))+ggplot2::geom_point(size=1.4,alpha=.8)+ggplot2::geom_hline(yintercept=-log10(.05),color='#b8c2c9',linetype=2)+ggplot2::geom_vline(xintercept=0,color='#dce1e4')+
        ggplot2::facet_wrap(~endpoint_label,ncol=3,scales='free_x')+ggplot2::scale_color_manual(values=c('Higher in cKO'='#c63d4f','Lower in cKO'='#397eb5','Below threshold'='#c3cbd0'))+
        ggrepel::geom_text_repel(data=lab,ggplot2::aes(label=native_name),size=2.6,max.overlaps=Inf,min.segment.length=Inf,seed=cfg$seed,show.legend=FALSE)+
        ggplot2::labs(title=paste('Regulon independent-mouse associations:',variant),subtitle='Continuous mean AUCell AUC; these are not gene DEGs',x='Mean AUCell difference (cKO - WT)',y='-log10 genotype-association BH FDR',color=NULL)+ggplot2::theme_minimal(base_size=11)+ggplot2::theme(legend.position='bottom')
    } else plot<-empty_plot('No eligible nonconstant regulon association tests')
    save_plot(plot,file.path(outdir,paste0('figures/REGULON_',variant,'_association_volcano')),12,8)
  }
  for(ep in ENDPOINTS) {
    e<-rank_effects(effects[effects$endpoint==ep & effects$tested,,drop=FALSE]);regs<-head(e$regulon_id,cfg$boxplot_top_regulons)
    x<-scores[scores$endpoint==ep & scores$regulon_id %in% regs,,drop=FALSE]
    if(nrow(x)) {
      x$regulon_id<-factor(x$regulon_id,levels=regs);x$condition<-factor(x$condition,levels=c('WT','cKO'))
      labs<-e[e$regulon_id %in% regs,,drop=FALSE];labs$regulon_id<-factor(labs$regulon_id,levels=regs);labs$label<-paste0('Difference=',signif(labs$effect_cKO_minus_WT,3),'; q=',format(labs$association_fdr,digits=2,scientific=TRUE))
      plot<-ggplot2::ggplot(x,ggplot2::aes(condition,aucell_mean,fill=condition))+ggplot2::geom_boxplot(width=.55,alpha=.25,outlier.shape=NA)+ggplot2::geom_point(ggplot2::aes(color=condition),position=ggplot2::position_jitter(width=.1,height=0,seed=cfg$seed),size=2)+
        ggplot2::facet_wrap(~regulon_id,ncol=3,scales='free_y')+ggplot2::scale_fill_manual(values=c(WT='#397eb5',cKO='#cf8334'))+ggplot2::scale_color_manual(values=c(WT='#397eb5',cKO='#cf8334'))+
        ggplot2::geom_text(data=labs,ggplot2::aes(x=1.5,y=Inf,label=label),inherit.aes=FALSE,vjust=1.3,size=2.4)+ggplot2::labs(title=paste(short_endpoint(ep),'mouse regulon scores'),subtitle='One point per independent mouse; top six selected post hoc; exploratory',x=NULL,y='Native mean AUCell AUC')+ggplot2::theme_minimal(base_size=10)+ggplot2::theme(legend.position='none')
      selections[[length(selections)+1L]]<-data.frame(plot=paste0(ep,'_top6_regulons'),endpoint=ep,regulon_id=regs,selection='post_hoc_global_association_FDR_then_P_then_absolute_effect')
    } else plot<-empty_plot(paste(short_endpoint(ep),'No eligible nonconstant regulon associations'))
    save_plot(plot,file.path(outdir,paste0('figures/',ep,'_top6_regulons')),12,7)
  }
  if(!is.null(no_self)) {
    sat<-result$primary[result$primary$TF=='Satb2',,drop=FALSE];sat$panel<-rep('Native full regulon',nrow(sat));ns<-no_self$primary;ns$panel<-rep('Self target excluded',nrow(ns));joined<-rbind(sat,ns);jj<-joined[joined$tested,,drop=FALSE]
    if(nrow(jj))plot<-ggplot2::ggplot(jj,ggplot2::aes(effect_cKO_minus_WT,short_endpoint(endpoint),color=panel))+ggplot2::geom_vline(xintercept=0,color='#c3cbd0')+ggplot2::geom_errorbar(ggplot2::aes(xmin=CI_low,xmax=CI_high),orientation='y',position=ggplot2::position_dodge(width=.4),width=.18)+ggplot2::geom_point(position=ggplot2::position_dodge(width=.4),size=2.3)+ggplot2::labs(title='Satb2 motif regulon: self-target sensitivity',subtitle='Native gene2weight may include TF own gene; this sensitivity is not independent binding evidence',x='cKO - WT mean AUCell AUC (mouse-level 95% CI)',y=NULL,color=NULL)+ggplot2::theme_minimal(base_size=11)+ggplot2::theme(legend.position='bottom') else plot<-empty_plot('No modelable motif-supported Satb2/self-target sensitivity; no proxy substituted')
    save_plot(plot,file.path(outdir,'figures/Satb2_regulon_self_target_sensitivity'),11,5)
  }
  empty<-data.frame(plot=character(),endpoint=character(),regulon_id=character(),selection=character())
  write_tsv(if(length(selections))do.call(rbind,selections) else empty,file.path(outdir,'REGULON_PLOT_SELECTION.tsv'))
}

main <- function() {
  options(digits=17);args<-args_parse(commandArgs(trailingOnly=TRUE))
  for(pkg in c('jsonlite','digest','ggplot2','ggrepel'))need(requireNamespace(pkg,quietly=TRUE),paste('Missing managed dependency:',pkg))
  config_receipt<-receipt(args$config);cfg<-jsonlite::read_json(args$config,simplifyVector=FALSE);need(identical(cfg$schema,'satb2-pyscenic/v1'),'Unexpected pySCENIC config schema');need(cfg$min_regulon_targets==10L,'Protocol fixes ten measured regulon targets')
  policy<-policy_defaults();policy$seed<-cfg$random_seed %||% cfg$seed %||% policy$seed;provided<-cfg$mouse_associations %||% list();need(all(names(provided)%in%names(policy)),'Unknown mouse-association policy');policy<-utils::modifyList(policy,provided)
  need(policy$bootstrap_iterations==200L && policy$expected_primary_WT==6L && policy$expected_primary_cKO==6L,'Protocol fixes 200 mouse bootstraps and 6WT/6cKO')
  need(policy$expected_partition_WT==6L && policy$expected_partition_cKO==3L && identical(policy$same_partition,'23HJT3LT4_lanes6_7'),'Protocol fixes known same-partition sensitivity')
  need(policy$heatmap_top_regulons>=1L && policy$heatmap_top_regulons<=60L && policy$boxplot_top_regulons==6L,'Invalid bounded plot policy')
  out<-args$outdir;need(!file.exists(out)||(dir.exists(out)&&!nzchar(Sys.readlink(out))&&!length(list.files(out,all.files=TRUE,no..=TRUE))),'Refusing nonempty or nonregular output directory')
  native<-normalizePath(args[['input-dir']],mustWork=TRUE);seal<-file.path(native,'COMPLETE.json');completed<-jsonlite::read_json(seal,simplifyVector=FALSE)
  need(identical(completed$schema,'satb2-pyscenic-native-complete/v1')&&identical(completed$status,'COMPLETE')&&identical(completed$config$sha256,config_receipt$sha256)&&identical(completed$analysis_id,cfg$analysis_id),'Incomplete native regulon producer or changed config/analysis')
  required<-c('REGULON_ACTIVITY_PER_MOUSE.tsv','REGULON_COVERAGE.tsv','REGULON_TARGETS.tsv','LINEAGE_STATUS.tsv','PROVENANCE.json','SATB2_NO_SELF_TARGET_ACTIVITY_PER_MOUSE.tsv','SATB2_NO_SELF_TARGET_COVERAGE.tsv')
  need(all(required%in%names(completed$output_files)),'Required native regulon inputs are not sealed')
  inputs<-list(native_complete=receipt(seal),samples=check_record(cfg$samples,'samples'))
  for(rel in names(completed$output_files)) {
    rec<-completed$output_files[[rel]];expected<-normalizePath(file.path(native,rel),mustWork=TRUE)
    need(startsWith(expected,paste0(native,'/')) && identical(normalizePath(rec$path,mustWork=TRUE),expected),'Native output receipt escapes declared root')
    inputs[[paste0('native/',rel)]]<-check_record(rec,paste('native',rel))
  }
  found<-list.files(native,recursive=TRUE,full.names=FALSE,all.files=TRUE,no..=TRUE);found<-found[!dir.exists(file.path(native,found))]
  need(setequal(found,c(names(completed$output_files),'COMPLETE.json')),'Unsealed native producer files')
  for(name in names(cfg$additional_hash_inputs %||% list()))inputs[[paste0('additional/',name)]]<-check_record(cfg$additional_hash_inputs[[name]],name)
  src<-sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1L]);src_receipt<-receipt(src)
  roster<-read_tsv(inputs$samples$path,c('sample','condition'));need(nrow(roster)==12L&&!anyDuplicated(roster$sample)&&sum(roster$condition=='WT')==6L&&sum(roster$condition=='cKO')==6L,'Expected exactly12 independent mice, six WT/six cKO')
  partition<-sample_partition(cfg,roster);roster$sequencing_partition<-partition$partition;inputs$partition_metadata<-partition$source
  columns<-c('lineage','endpoint','endpoint_id','source_endpoint_id','regulon_id','native_name','TF','n_targets','n_background_genes','eligible')
  scores<-read_tsv(file.path(native,'REGULON_ACTIVITY_PER_MOUSE.tsv'),c(columns,'sample','condition','n_nuclei','aucell_mean','score_unit'))
  coverage<-read_tsv(file.path(native,'REGULON_COVERAGE.tsv'),c(columns,'reason'));statuses<-read_tsv(file.path(native,'LINEAGE_STATUS.tsv'),c('lineage','status','n_regulons','n_nuclei','n_genes'))
  no_scores<-read_tsv(file.path(native,'SATB2_NO_SELF_TARGET_ACTIVITY_PER_MOUSE.tsv'),c(columns,'sample','condition','n_nuclei','aucell_mean','score_unit','network_variant'));no_cov<-read_tsv(file.path(native,'SATB2_NO_SELF_TARGET_COVERAGE.tsv'),c(columns,'reason','network_variant'))
  validate_native_scores(no_scores,no_cov,roster,statuses,supplement=TRUE);no_scores$sequencing_partition<-roster$sequencing_partition[match(no_scores$sample,roster$sample)]
  validate_native_scores(scores,coverage,roster,statuses);scores$sequencing_partition<-roster$sequencing_partition[match(scores$sample,roster$sample)]
  dir.create(out,showWarnings=FALSE,recursive=TRUE)
  tryCatch({
    set.seed(policy$seed);result<-analyze_regulons(scores,coverage,roster,policy);no_self<-analyze_regulons(no_scores,no_cov,roster,policy,'BH_available_Satb2_motif_no_self_target_x_five_endpoints')
    for(pair in list(c('primary','ASSOCIATIONS'),c('sensitivity','SAME_PARTITION_ASSOCIATIONS'),c('loo','LOO_MODELS'),c('stability','LOO_STABILITY'),c('bootstrap','MOUSE_BOOTSTRAP')))write_tsv(result[[pair[1L]]],file.path(out,paste0('REGULON_',pair[2L],'.tsv')))
    for(pair in list(c('primary','ASSOCIATIONS'),c('sensitivity','SAME_PARTITION_ASSOCIATIONS'),c('loo','LOO_MODELS'),c('stability','LOO_STABILITY'),c('bootstrap','MOUSE_BOOTSTRAP')))write_tsv(no_self[[pair[1L]]],file.path(out,paste0('SATB2_NO_SELF_TARGET_',pair[2L],'.tsv')))
    write_tsv(no_scores,file.path(out,'SATB2_NO_SELF_TARGET_MOUSE_SCORES.tsv'))
    eligible_cols<-c('lineage','endpoint','endpoint_id','source_endpoint_id','regulon_id','native_name','TF','n_targets','n_background_genes','n_scored_mice','n_expected_mice','n_WT','n_cKO','tested','model_status')
    write_tsv(result$primary[,eligible_cols,drop=FALSE],file.path(out,'REGULON_MODEL_ELIGIBILITY.tsv'));write_tsv(scores,file.path(out,'REGULON_MOUSE_SCORES.tsv'));write_tsv(roster[,c('sample','condition','sequencing_partition')],file.path(out,'MOUSE_PARTITION_AUDIT.tsv'))
    if(isTRUE(policy$plots))make_plots(result,scores,policy,out,no_self)
    summary<-data.frame(endpoint=ENDPOINTS,n_tested=vapply(ENDPOINTS,function(ep)sum(result$primary$endpoint==ep & result$primary$tested),integer(1)),n_association_fdr_lt05=vapply(ENDPOINTS,function(ep)sum(result$primary$endpoint==ep & !is.na(result$primary$association_fdr) & result$primary$association_fdr<.05),integer(1)))
    sat_status<-do.call(rbind,lapply(ENDPOINTS,function(ep){x<-result$primary[result$primary$endpoint==ep & result$primary$TF=='Satb2',,drop=FALSE];n<-no_self$primary[no_self$primary$endpoint==ep,,drop=FALSE];data.frame(endpoint=ep,lineage=lineage_for(ep),TF='Satb2',native_positive_regulon_present=nrow(x)>0L,native_regulon_ids=paste(x$regulon_id,collapse=';'),n_targets=if(nrow(x))x$n_targets[1L] else NA_integer_,native_AUCell_available=any(scores$endpoint==ep & scores$TF=='Satb2'),primary_association_tested=any(x$tested),primary_status=if(nrow(x))paste(x$model_status,collapse=';') else 'NO_MOTIF_SUPPORTED_NATIVE_SATB2_REGULON',self_excluded_AUCell_available=any(no_scores$endpoint==ep),self_excluded_association_tested=any(n$tested),self_excluded_status=paste(n$model_status,collapse=';'))}))
    write_tsv(sat_status,file.path(out,'SATB2_STATUS.tsv'))
    write_tsv(summary,file.path(out,'MODEL_SUMMARY.tsv'));capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
    for(name in names(inputs))need(identical(receipt(inputs[[name]]$path),inputs[[name]]),paste('Input changed:',name))
    need(identical(receipt(args$config),config_receipt)&&identical(receipt(src),src_receipt),'Source/config changed during analysis')
    json_write(list(schema='satb2-regulon-mouse-provenance/v1',analysis_id=cfg$analysis_id,config=config_receipt,source=src_receipt,inputs=inputs,policy=policy,
      model='mean_native_nucleus_AUCell_AUC ~ condition; unweighted independent-mouse continuous linear model',
      primary_FDR='BH jointly across every tested lineage-specific positive motif-supported regulon x its matching endpoints, both lineages',
      sensitivity_FDR='Separate global BH across same-partition 6WT/3cKO associations; separate Satb2 no-self family across available five endpoints',native_activity='Unthresholded nucleus-level native AUCell scores averaged per exact endpoint and mouse by native producer',
      interpretation='Motif-supported expression-derived candidate regulons; enrichment is not direct TF binding or causal regulation; AUC is not TF RNA expression',
      confounding='Genotype partly confounds technical partition; analyses exploratory; C8 mixed identity remains under review',
      original_DESeq2_and_LR_unchanged=TRUE,no_nucleus_pseudoreplication=TRUE,no_activity_threshold_filter=TRUE,positive_motif_regulons_only=TRUE,min_measured_regulon_targets=10L,Satb2_self_target_sensitivity='Separate native AUCell with identical expression ranks/background and only selfgene removed when supported; original full regulons unchanged',
      software=lapply(c('ggplot2','ggrepel','jsonlite','digest'),function(p)list(package=p,version=as.character(utils::packageVersion(p))))),file.path(out,'PROVENANCE.json'))
    paths<-list.files(out,recursive=TRUE,full.names=TRUE);paths<-paths[!dir.exists(paths)];receipts<-lapply(paths,function(p){x<-receipt(p);x$path<-NULL;x});names(receipts)<-substring(paths,nchar(out)+2L)
    json_write(list(schema='satb2-regulon-mouse-complete/v1',status='COMPLETE',analysis_id=cfg$analysis_id,config=config_receipt,source=src_receipt,inputs=inputs,outputs=receipts,
      finished_UTC=format(Sys.time(),tz='UTC',usetz=TRUE),checks=list(all12_mouse_roster=TRUE,native_AUCell_means_preserved=TRUE,global_primary_BH_both_lineages=TRUE,separate_sensitivity_BH=TRUE,
      no_activity_threshold_filter=TRUE,no_nucleus_pseudoreplication=TRUE,no_new_gene_DESeq2=TRUE,positive_motif_regulons_only=TRUE,inputs_config_source_SHA256_unchanged=TRUE)),file.path(out,'COMPLETE.json'))
  },error=function(e){if(!file.exists(file.path(out,'FAILED.json')))json_write(list(status='FAILED',message=conditionMessage(e)),file.path(out,'FAILED.json'));stop(e)})
}
if(sys.nframe()==0L)main()
