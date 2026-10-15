library(hdf5r)
library(jsonlite)
archive <- "/Users/mujahid357/Desktop/SATB2_LiSC_archive_transfer_20261012.BmzvXa"
out <- "/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/audit_cellbender"
dir.create(out, recursive=TRUE, showWarnings=FALSE)
qcbase <- file.path(archive,"office_analysis/runs/00_source_qc_audit_v1/results/source_qc")
samples <- sort(list.dirs(file.path(archive,"cellbender_out"),recursive=FALSE,full.names=FALSE))
marker_symbols <- c("Gad1","Gad2","Slc17a7","Slc6a1","Pvalb","Sst","Vip","Reln","Cck","Lamp5","Ndnf","Satb2","Wfs1","Ndst4","Col11a1","Asic2","Calb1","Sox5","Hs6st3","Iqgap2","Ryr3","Epha6","Aqp4","Gfap","Plp1","Mbp","Malat1","Rbfox3","Tbr1")
all_summary <- list(); all_markers <- list(); all_curve <- list(); all_ambient <- list()
scan_matrix <- function(f, ids, marker_rows, paired_barcodes) {
  bc <- f[["matrix/barcodes"]][]
  p <- f[["matrix/indptr"]][]
  pair <- bc %in% paired_barcodes
  n <- length(bc); nnz <- tail(p,1)
  totals <- numeric(n)
  msum <- numeric(length(marker_rows)); mdetection <- integer(length(marker_rows))
  nonnegative <- TRUE; integer_counts <- TRUE
  row_to_marker <- integer(length(ids)); row_to_marker[marker_rows] <- seq_along(marker_rows)
  if(nnz>0) for(start in seq(1,nnz,by=1000000)) {
    end <- min(start+999999,nnz)
    d <- f[["matrix/data"]][start:end]
    ix <- f[["matrix/indices"]][start:end]+1L
    col <- findInterval(seq.int(start-1,end-1),p[-length(p)])
    nonnegative <- nonnegative && all(is.finite(d) & d>=0)
    integer_counts <- integer_counts && all(d==floor(d))
    cs <- rowsum(as.matrix(as.numeric(d)),col,reorder=FALSE)
    ci <- as.integer(rownames(cs)); totals[ci] <- totals[ci]+cs[,1]
    m <- row_to_marker[ix]
    keep <- m>0L & pair[col] & d>0
    if(any(keep)) {
      ms <- rowsum(as.matrix(as.numeric(d[keep])),m[keep],reorder=FALSE)
      mi <- as.integer(rownames(ms)); msum[mi] <- msum[mi]+ms[,1]
      mdetection <- mdetection+tabulate(m[keep],nbins=length(marker_rows))
    }
  }
  list(barcodes=bc,totals=totals,detected_genes=diff(p),marker_sum=msum,marker_detection=mdetection,nonnegative=nonnegative,integer_counts=integer_counts)
}
for(sample in samples) {
  cat("Diagnostic audit",sample,"\n")
  cbpath <- file.path(archive,"cellbender_out",sample,paste0(sample,"_filtered.h5"))
  fullpath <- file.path(archive,"cellbender_out",sample,paste0(sample,".h5"))
  originalpath <- file.path(archive,"cellranger_out",sample,"outs/filtered_feature_bc_matrix.h5")
  q <- read.delim(file.path(qcbase,paste0(sample,"_nucleus_QC.tsv")),check.names=FALSE)
  q <- q[q$passes_primary_QC,]
  full <- H5File$new(fullpath,mode="r")
  full_bc <- full[["matrix/barcodes"]][]
  full_p <- full[["matrix/indptr"]][]
  qi <- match(q$barcode,full_bc)
  stopifnot(!anyNA(qi))
  full_nnz_qc <- full_p[qi+1]-full_p[qi]
  analyzed <- full[["metadata/barcodes_analyzed"]][]
  q$cellbender_analyzed <- q$barcode %in% analyzed
  q$cellbender_nonzero_full <- full_nnz_qc>0
  full$close_all(); rm(full,full_bc,full_p,qi,full_nnz_qc,analyzed); gc(verbose=FALSE)
  cb <- H5File$new(cbpath,mode="r")
  ids <- cb[["matrix/features/id"]][]
  symbols <- cb[["matrix/features/name"]][]
  marker_rows <- match(marker_symbols,symbols)
  stopifnot(!anyNA(marker_rows))
  paired <- q$barcode[q$cellbender_analyzed & q$cellbender_nonzero_full]
  cbs <- scan_matrix(cb,ids,marker_rows,paired)
  original <- H5File$new(originalpath,mode="r")
  stopifnot(identical(ids,original[["matrix/features/id"]][]))
  ori <- scan_matrix(original,ids,marker_rows,paired)
  original$close_all(); rm(original)
  qm <- match(q$barcode,cbs$barcodes)
  q$cellbender_filtered_present <- !is.na(qm)
  q$cellbender_nCount_RNA <- ifelse(is.na(qm),0,cbs$totals[qm])
  q$cellbender_nFeature_RNA <- ifelse(is.na(qm),0,cbs$detected_genes[qm])
  q$paired_sensitivity_usable <- q$cellbender_analyzed & q$cellbender_nonzero_full & q$cellbender_nCount_RNA>0
  q$umi_fraction_removed <- 1-q$cellbender_nCount_RNA/q$nCount_RNA
  stopifnot(all(q$cellbender_nonzero_full == (q$cellbender_nCount_RNA>0)))
  write.table(q,file.path(out,paste0(sample,"_original_cellbender_primary_paired_roster.tsv")),sep="\t",quote=FALSE,row.names=FALSE)
  record <- list(sample=sample,full_corrected_path=fullpath,filtered_corrected_path=cbpath,original_filtered_path=originalpath,primary_QC_nuclei=nrow(q),paired_supported=sum(q$paired_sensitivity_usable),zero_corrected=sum(!q$cellbender_nonzero_full),not_analyzed=sum(!q$cellbender_analyzed),corrected_filtered_barcodes=length(cbs$barcodes),count_nonnegative=cbs$nonnegative,count_integer=cbs$integer_counts,total_original_UMIs_paired=sum(q$nCount_RNA[q$paired_sensitivity_usable]),total_corrected_UMIs_paired=sum(q$cellbender_nCount_RNA[q$paired_sensitivity_usable]),median_fraction_removed_paired=median(q$umi_fraction_removed[q$paired_sensitivity_usable]),q95_fraction_removed_paired=as.numeric(quantile(q$umi_fraction_removed[q$paired_sensitivity_usable],.95)),max_fraction_removed_paired=max(q$umi_fraction_removed[q$paired_sensitivity_usable]),zeroed_original_median_UMI=median(q$nCount_RNA[!q$cellbender_nonzero_full]),zeroed_original_median_genes=median(q$nFeature_RNA[!q$cellbender_nonzero_full]))
  elbo_fields <- c("train","test")
  for(kind in elbo_fields) {
    ep <- cb[[paste0("metadata/learning_curve_",kind,"_epoch")]][]
    el <- cb[[paste0("metadata/learning_curve_",kind,"_elbo")]][]
    all_curve[[length(all_curve)+1]] <- data.frame(sample=sample,split=kind,epoch=ep,ELBO=el)
    record[[paste0(kind,"_finite")]] <- all(is.finite(el))
    record[[paste0(kind,"_last_epoch")]] <- tail(ep,1)
    record[[paste0(kind,"_last_ELBO")]] <- tail(el,1)
    record[[paste0(kind,"_last20_slope_per_epoch")]] <- unname(coef(lm(tail(el,20)~tail(ep,20)))[2])
    record[[paste0(kind,"_last20_sd")]] <- sd(tail(el,20))
  }
  record$estimator <- cb[["metadata/estimator"]][]
  record$target_fpr <- cb[["metadata/target_false_positive_rate"]][]
  ambient <- cb[["global_latents/ambient_expression"]][]
  ao <- head(order(ambient,decreasing=TRUE),20)
  all_ambient[[length(all_ambient)+1]] <- data.frame(sample=sample,rank=seq_along(ao),gene_id=ids[ao],gene_symbol=symbols[ao],ambient_probability=ambient[ao])
  all_markers[[length(all_markers)+1]] <- data.frame(sample=sample,gene_id=ids[marker_rows],gene_symbol=marker_symbols,paired_nuclei=length(paired),original_UMIs=ori$marker_sum,corrected_UMIs=cbs$marker_sum,fraction_removed=ifelse(ori$marker_sum>0,1-cbs$marker_sum/ori$marker_sum,NA),original_detected_nuclei=ori$marker_detection,corrected_detected_nuclei=cbs$marker_detection)
  all_summary[[length(all_summary)+1]] <- record
  cb$close_all()
  rm(cb,cbs,ori,q,ids,symbols,marker_rows,paired,qm,record,ambient,ao)
  gc(verbose=FALSE)
}
write_json(all_summary,file.path(out,"local_cellbender_diagnostics.json"),pretty=TRUE,auto_unbox=TRUE)
write.table(do.call(rbind,lapply(all_summary,as.data.frame)),file.path(out,"local_cellbender_diagnostics.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
write.table(do.call(rbind,all_markers),file.path(out,"marker_removal_paired.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
write.table(do.call(rbind,all_curve),file.path(out,"learning_curves.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
write.table(do.call(rbind,all_ambient),file.path(out,"top_ambient_genes.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
cat("DONE paired:",sum(vapply(all_summary,function(x)x$paired_supported,numeric(1))),"\n")
