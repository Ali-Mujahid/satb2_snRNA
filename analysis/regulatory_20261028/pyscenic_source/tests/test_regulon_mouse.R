#!/usr/bin/env Rscript
# Synthetic scores only. Sourcing the analysis script is inert.
args<-commandArgs(trailingOnly=TRUE);stopifnot(length(args)==1L)
source(args[1L])
expect_error<-function(expr){x<-try(force(expr),silent=TRUE);stopifnot(inherits(x,'try-error'))}
roster<-data.frame(sample=c(paste0('wt_fake',1:6),paste0('ko_fake',1:6)),condition=rep(c('WT','cKO'),each=6),sequencing_partition=c(rep('23HJT3LT4_lanes6_7',6),rep('23HJT3LT4_lanes6_7',3),rep('other',3)))
cover<-do.call(rbind,lapply(ENDPOINTS,function(ep){lin<-lineage_for(ep);data.frame(lineage=lin,endpoint=ep,endpoint_id=endpoint_id(ep),source_endpoint_id=source_endpoint_id(ep),regulon_id=paste0(lin,'::FakeTF(+)'),native_name='FakeTF(+)',TF='FakeTF',n_targets=10L,n_background_genes=100L,eligible=TRUE,reason='ELIGIBLE')}))
scores<-do.call(rbind,lapply(seq_len(nrow(cover)),function(i){c<-cover[i,,drop=FALSE];x<-c[rep(1L,12),];x$sample<-roster$sample;x$condition<-roster$condition;x$n_nuclei<-20L;x$aucell_mean<-.2+rep(c(.001,.003,.006,.009,.014,.021),2)+ifelse(x$condition=='cKO',.05,0);x$score_unit<-'native_AUCell_AUC';x$sequencing_partition<-roster$sequencing_partition;x}))
status<-data.frame(lineage=c('CA1','IN'),status='COMPLETE_WITH_REGULONS',n_regulons=1L,n_nuclei=c(20456L,2687L),n_genes=100L)
validate_native_scores(scores,cover,roster,status)
p<-policy_defaults();p$plots<-FALSE;set.seed(1337);a<-analyze_regulons(scores,cover,roster,p)
stopifnot(nrow(a$primary)==5L,all(a$primary$tested),all(abs(a$primary$effect_cKO_minus_WT-.05)<1e-12),all(a$primary$n_WT==6L),all(a$primary$n_cKO==6L),all(a$sensitivity$n_WT==6L),all(a$sensitivity$n_cKO==3L))
stopifnot(max(abs(a$primary$association_fdr-p.adjust(a$primary$pvalue,'BH')))<1e-14,nrow(a$loo)==60L,nrow(a$bootstrap)==5L,all(a$bootstrap$bootstrap_iterations==200L),all(a$stability$n_loo_same_sign==12L))
x<-data.frame(sample=roster$sample,condition=roster$condition,escore=rep(.2,12));stopifnot(!model_score(x,p,6,6)$tested)
x$escore<-ifelse(x$condition=='WT',.2,.3);stopifnot(model_score(x,p,6,6)$model_status=='ZERO_OR_NUMERICALLY_ZERO_WITHIN_GROUP_VARIANCE')
x$escore<-scores$aucell_mean[1:12];stopifnot(!model_score(x[-1,],p,6,6)$tested)
bad<-scores;bad$regulon_id[1]<-'IN::FakeTF(+)';expect_error(validate_native_scores(bad,cover,roster,status))
bad<-scores[-1,];expect_error(validate_native_scores(bad,cover,roster,status))
bad<-scores;bad$score_unit[1]<-'RNA';expect_error(validate_native_scores(bad,cover,roster,status))
bad<-scores;bad$aucell_mean[1]<-NA;expect_error(validate_native_scores(bad,cover,roster,status))
bad<-cover;bad$n_targets[1]<-9L;expect_error(validate_native_scores(scores,bad,roster,status))
# Zero-regulon lineages remain explicit; all-zero/empty activity is never imputed.
e<-cover[FALSE,];ss<-scores[FALSE,];zero<-status;zero$status<-'EMPTY_NO_SUPPORTED_REGULONS';zero$n_regulons<-0L
validate_native_scores(ss,e,roster,zero);empty<-analyze_regulons(ss,e,roster,p);stopifnot(nrow(empty$primary)==0L,nrow(empty$sensitivity)==0L)
# Supplement unsupported in all five populations is a legitimate missing result.
ncov<-cover;ncov$regulon_id<-paste0(ncov$lineage,'::Satb2(+)_NO_SELF_TARGET');ncov$native_name<-'Satb2(+)';ncov$TF<-'Satb2';ncov$n_targets<-0L;ncov$eligible<-FALSE;ncov$reason<-'NO_SUPPORTED_NATIVE_SATB2_REGULON';ncov$network_variant<-'SATB2_native_regulon_self_target_excluded'
validate_native_scores(ss,ncov,roster,zero,supplement=TRUE);n<-analyze_regulons(ss,ncov,roster,p,'BH_Satb2_supplement');stopifnot(nrow(n$primary)==5L,!any(n$primary$tested),all(is.na(n$primary$association_fdr)),all(n$primary$fdr_family=='BH_Satb2_supplement'))
# Family separation: adding a supplementary P value cannot enter primary BH.
stopifnot(!identical(unique(a$primary$fdr_family),unique(n$primary$fdr_family)))
cat('PASS: native identities, all12 mice, min10 targets, mean-difference units, global BH, separate sensitivity, LOO/200bootstrap, constants, no-self and empty lineages\n')
