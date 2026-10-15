library(hdf5r)
archive <- "/Users/mujahid357/Desktop/SATB2_LiSC_archive_transfer_20261012.BmzvXa"
out <- "/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/audit_cellbender"
membership <- read.delim(file.path(out,"identity_marker_diagnostic_membership.tsv"),check.names=FALSE)
markers <- c("Gad1","Gad2","Slc17a7","Pvalb","Sst","Vip","Reln","Cck","Ndnf","Satb2","Wfs1","Ndst4","Col11a1","Asic2","Calb1","Sox5","Hs6st3","Iqgap2","Ryr3","Epha6")
records <- list()
for(sample in sort(unique(membership$sample))) {
  cat("Within reference marker removal",sample,"\n")
  mm <- membership[membership$sample==sample & membership$paired_sensitivity_usable,]
  IN <- mm$barcode[mm$ReferenceHighSupport_IN]
  CA1 <- mm$barcode[mm$ReferenceHighSupport_CA1_ProS]
  for(kind in c("original","corrected")) {
    path <- if(kind=="original") file.path(archive,"cellranger_out",sample,"outs/filtered_feature_bc_matrix.h5") else file.path(archive,"cellbender_out",sample,paste0(sample,"_filtered.h5"))
    f <- H5File$new(path,mode="r")
    names <- f[["matrix/features/name"]][]
    ids <- f[["matrix/features/id"]][]
    mr <- match(markers,names);stopifnot(!anyNA(mr))
    lookup <- integer(length(names));lookup[mr] <- seq_along(markers)
    bc <- f[["matrix/barcodes"]][];p <- f[["matrix/indptr"]][]
    inc <- bc %in% IN;ca1c <- bc %in% CA1
    sums <- matrix(0,nrow=length(markers),ncol=2)
    detects <- matrix(0L,nrow=length(markers),ncol=2)
    for(start in seq(1,tail(p,1),by=1000000)) {
      end <- min(start+999999,tail(p,1))
      ix <- f[["matrix/indices"]][start:end]+1L
      m <- lookup[ix];pos <- which(m>0)
      if(!length(pos)) next
      d <- f[["matrix/data"]][start:end][pos]
      col <- findInterval(start+pos-2,p[-length(p)])
      m <- m[pos]
      for(g in 1:2) {
        keep <- (if(g==1) inc[col] else ca1c[col]) & d>0
        if(any(keep)) {
          sm <- rowsum(as.matrix(as.numeric(d[keep])),m[keep],reorder=FALSE)
          si <- as.integer(rownames(sm));sums[si,g] <- sums[si,g]+sm[,1]
          detects[,g] <- detects[,g]+tabulate(m[keep],nbins=length(markers))
        }
      }
    }
    for(g in 1:2) records[[length(records)+1]] <- data.frame(sample=sample,frozen_reference_population=if(g==1)"ReferenceHighSupport_IN" else "ReferenceHighSupport_CA1_ProS",count_kind=kind,gene_id=ids[mr],gene_symbol=markers,paired_nuclei=if(g==1)length(IN)else length(CA1),UMIs=sums[,g],detected_nuclei=detects[,g])
    f$close_all();rm(f,names,ids,mr,lookup,bc,p,inc,ca1c,sums,detects);gc(verbose=FALSE)
  }
}
write.table(do.call(rbind,records),file.path(out,"marker_removal_within_frozen_reference.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
