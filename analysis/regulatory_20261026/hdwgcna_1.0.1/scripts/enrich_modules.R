#!/usr/bin/env Rscript
# Offline enrichment of completed hdWGCNA modules; called only by the native workflow.
suppressPackageStartupMessages({library(jsonlite);library(AnnotationDbi);library(org.Mm.eg.db);library(GO.db);library(ggplot2)})
args <- commandArgs(trailingOnly=TRUE)
arg <- function(key) {i<-match(key,args);if(is.na(i)||i==length(args))stop("Missing ",key);args[[i+1]]}
cfg<-fromJSON(arg("--config"),simplifyVector=FALSE)
res<-normalizePath(arg("--results"),mustWork=TRUE)
out<-arg("--outdir")
if(dir.exists(out)&&length(list.files(out,all.files=TRUE,no..=TRUE)))stop("Refusing nonempty enrichment output")
dir.create(out,recursive=TRUE)
read_tsv<-function(p)read.delim(p,check.names=FALSE,stringsAsFactors=FALSE)
write_tsv<-function(x,p)write.table(x,p,sep="\t",quote=FALSE,row.names=FALSE,na="NA")
all_results<-list(); coverage<-list(); overlaps<-list(); lr_links<-list(); de_coverage<-list(); skipped<-character()
endpoint_lineage<-c(CA1_Deep="CA1",CA1_Superficial="CA1",IN_C1="IN",IN_C6="IN",IN_C8="IN")
endpoints<-cfg$existing_deseq2
lr<-read_tsv(cfg$strict_lr_pooled$path)
for(lineage in c("CA1","IN")){
  folder<-file.path(res,lineage,"pooled")
  if(!file.exists(file.path(folder,"module_membership.tsv"))){skipped<-c(skipped,lineage);next}
  mods<-read_tsv(file.path(folder,"module_membership.tsv"))
  gene_col<-intersect(c("gene_symbol","gene_name","gene"),names(mods))[[1]]
  stopifnot("module"%in%names(mods),!anyDuplicated(mods[[gene_col]]))
  names(mods)[names(mods)==gene_col]<-"gene_symbol"
  bg<-unique(mods$gene_symbol)
  modules<-setdiff(unique(mods$module),c("grey","gray","unassigned","0"))
  valid_keys<-intersect(bg,keys(org.Mm.eg.db,keytype="SYMBOL"))
  annotations<-if(length(valid_keys))suppressMessages(AnnotationDbi::select(org.Mm.eg.db,keys=valid_keys,keytype="SYMBOL",columns=c("GOALL","ONTOLOGYALL")))else data.frame(SYMBOL=character(),GOALL=character(),ONTOLOGYALL=character())
  annotations<-unique(annotations[!is.na(annotations$GOALL)&annotations$ONTOLOGYALL=="BP",c("SYMBOL","GOALL")])
  annotated_bg<-intersect(bg,annotations$SYMBOL)
  termgenes<-split(annotations$SYMBOL,annotations$GOALL)
  termgenes<-lapply(termgenes,unique)
  termgenes<-termgenes[vapply(termgenes,length,integer(1))>=10&vapply(termgenes,length,integer(1))<=2000]
  for(module in modules){
    g<-intersect(mods$gene_symbol[mods$module==module],annotated_bg)
    coverage[[length(coverage)+1]]<-data.frame(lineage,module,n_module_genes=sum(mods$module==module),n_module_BP_annotated=length(g),n_network_genes=length(bg),n_network_BP_annotated=length(annotated_bg),n_terms_tested=length(termgenes))
    if(!length(g))next
    for(term in names(termgenes)){
      tg<-intersect(termgenes[[term]],annotated_bg);hit<-intersect(g,tg)
      a<-length(hit);b<-length(g)-a;c<-length(tg)-a;d<-length(annotated_bg)-a-b-c
      ft<-fisher.test(matrix(c(a,b,c,d),nrow=2,byrow=TRUE),alternative="greater")
      all_results[[length(all_results)+1]]<-data.frame(lineage,module,GO_ID=term,overlap=a,module_annotated_n=length(g),term_background_n=length(tg),background_n=length(annotated_bg),odds_ratio=unname(ft$estimate),pvalue=ft$p.value,overlap_symbols=paste(sort(hit),collapse=";"))
    }
  }
  for(ep in names(endpoint_lineage)[endpoint_lineage==lineage]){
    de<-read_tsv(endpoints[[ep]]$path)
    gene<-intersect(c("gene","gene_symbol","symbol"),names(de))[[1]]
    stopifnot(all(c("log2FoldChange","padj")%in%names(de)))
    duplicated_symbol<-duplicated(de[[gene]])|duplicated(de[[gene]],fromLast=TRUE)
    multi_stable_id<-if("gene_ids"%in%names(mods))grepl(";",mods$gene_ids,fixed=TRUE)else rep(TRUE,nrow(mods))
    unambiguous_network<-mods$gene_symbol[!multi_stable_id]
    tested<-is.finite(de$padj)&is.finite(de$log2FoldChange)&!duplicated_symbol&de[[gene]]%in%unambiguous_network
    endpoint_bg<-intersect(bg,de[[gene]][tested])
    pass<-tested&de$padj<0.05
    de_coverage[[length(de_coverage)+1]]<-data.frame(lineage,endpoint=ep,n_network_genes=length(bg),n_finite_unambiguous_DESeq2_background=length(endpoint_bg),excluded_multi_stable_ID_features=sum(multi_stable_id),excluded_duplicate_fit_symbols=sum(duplicated_symbol),background_policy="Network genes intersect endpoint-matched finite adjusted-P eligible genes with unique original symbol and one stable ID")
    selected<-list(all=unique(de[[gene]][pass]),up=unique(de[[gene]][pass&de$log2FoldChange>0]),down=unique(de[[gene]][pass&de$log2FoldChange<0]))
    for(direction in names(selected))for(module in modules){
      mg<-intersect(mods$gene_symbol[mods$module==module],endpoint_bg);dg<-intersect(selected[[direction]],endpoint_bg);hit<-intersect(mg,dg)
      if(!length(mg)||!length(endpoint_bg))next
      a<-length(hit);b<-length(mg)-a;c<-length(dg)-a;d<-length(endpoint_bg)-a-b-c
      ft<-fisher.test(matrix(c(a,b,c,d),nrow=2,byrow=TRUE),alternative="greater")
      overlaps[[length(overlaps)+1]]<-data.frame(lineage,endpoint=ep,module,direction,overlap=a,module_tested_n=length(mg),DEG_background_n=length(dg),endpoint_eligible_network_background_n=length(endpoint_bg),odds_ratio=unname(ft$estimate),pvalue=ft$p.value,overlap_symbols=paste(sort(hit),collapse=";"),source_endpoint=endpoints[[ep]]$source_endpoint_id)
    }
    ep_lr<-unique(c(lr$ligand_complex[lr$source==ep],lr$receptor_complex[lr$target==ep]))
    # Strict existing source uses only complete single-symbol pairs at present;
    # split complexes for a transparent per-component module cross-reference.
    components<-unique(unlist(strsplit(ep_lr,"_",fixed=TRUE)))
    components<-unique(c("Satb2",components))
    for(g in components){
      mm<-mods[mods$gene_symbol==g,,drop=FALSE]
      dd<-de[de[[gene]]==g,,drop=FALSE]
      lr_links[[length(lr_links)+1]]<-data.frame(lineage,endpoint=ep,gene_symbol=g,module=if(nrow(mm))as.character(mm$module[[1]])else NA_character_,network_eligible=nrow(mm)>0,DESeq2_log2FC=if(nrow(dd)==1)dd$log2FoldChange[[1]]else NA_real_,DESeq2_padj=if(nrow(dd)==1)dd$padj[[1]]else NA_real_,deseq_rows=nrow(dd),strict_LR_component=g%in%components[components!="Satb2"])
    }
  }
}
go<-if(length(all_results))do.call(rbind,all_results)else data.frame(lineage=character(),module=character(),GO_ID=character(),overlap=integer(),module_annotated_n=integer(),term_background_n=integer(),background_n=integer(),odds_ratio=numeric(),pvalue=numeric(),overlap_symbols=character(),padj=numeric(),term=character())
if(nrow(go)){
  go$padj<-p.adjust(go$pvalue,"BH")
  descriptions<-suppressMessages(AnnotationDbi::select(GO.db,keys=unique(go$GO_ID),keytype="GOID",columns="TERM"))
  go$term<-descriptions$TERM[match(go$GO_ID,descriptions$GOID)]
  go<-go[order(go$padj,go$pvalue),]
}
ov<-if(length(overlaps))do.call(rbind,overlaps)else data.frame(lineage=character(),endpoint=character(),module=character(),direction=character(),overlap=integer(),module_tested_n=integer(),DEG_background_n=integer(),endpoint_eligible_network_background_n=integer(),odds_ratio=numeric(),pvalue=numeric(),overlap_symbols=character(),source_endpoint=character(),padj=numeric())
if(nrow(ov)){ov$padj<-p.adjust(ov$pvalue,"BH");ov<-ov[order(ov$padj,ov$pvalue),]}
write_tsv(go,file.path(out,"GO_BP_ALL_TESTS.tsv"));write_tsv(if(nrow(go))go[go$padj<.05,,drop=FALSE]else go,file.path(out,"GO_BP_FDR05.tsv"))
write_tsv(if(length(coverage))do.call(rbind,coverage)else data.frame(),file.path(out,"GO_ANNOTATION_COVERAGE.tsv"))
write_tsv(ov,file.path(out,"ORIGINAL_DEG_MODULE_ENRICHMENT.tsv"));write_tsv(if(length(lr_links))do.call(rbind,lr_links)else data.frame(),file.path(out,"STRICT_LR_AND_SATB2_MODULE_LINKS.tsv"))
write_tsv(if(length(de_coverage))do.call(rbind,de_coverage)else data.frame(),file.path(out,"ORIGINAL_DEG_ELIGIBILITY_COVERAGE.tsv"))
if(nrow(go)){
  significant<-go[go$padj<.05&go$overlap>0,,drop=FALSE]
  if(nrow(significant)){
    top<-do.call(rbind,lapply(split(significant,paste(significant$lineage,significant$module)),head,3))
    top<-head(top[order(top$padj),],60)
    top$label<-paste(top$lineage,top$module,sep=" · ")
    p<-ggplot(top,aes(x=label,y=reorder(term,-log10(padj)),size=overlap,color=-log10(pmax(padj,1e-300))))+geom_point()+scale_color_viridis_c()+labs(x="Fitted coexpression module",y=NULL,color="−log10 FDR",size="Overlap genes",title="GO biological-process enrichment",subtitle="Measured network genes with BP annotation; BH across both lineage/module/term test families")+theme_bw(base_size=10)+theme(axis.text.x=element_text(angle=60,hjust=1))
    ggsave(file.path(out,"GO_BP_dotplot.pdf"),p,width=13,height=max(7,min(18,nrow(top)*.24)))
    ggsave(file.path(out,"GO_BP_dotplot.svg"),p,width=13,height=max(7,min(18,nrow(top)*.24)))
    ggsave(file.path(out,"GO_BP_dotplot.png"),p,width=13,height=max(7,min(18,nrow(top)*.24)),dpi=180)
  }
}
write_json(list(status=if(length(skipped))"PARTIAL_NETWORK_UNAVAILABLE"else"COMPLETE",skipped_lineages=skipped,GO_db=as.character(packageVersion("GO.db")),org_Mm_eg_db=as.character(packageVersion("org.Mm.eg.db")),GO_background="Fitted lineage network gene symbols with mouse GO BP annotation",GO_term_sizes=c(10,2000),GO_FDR_family="All tested lineage × non-grey module × eligible BP term",DEG_background="Fitted network genes intersect matching endpoint finite adjusted-P gene universe, excluding multi-ID collapsed features and duplicate fit symbols",DEG_FDR_family="All tested lineage × matching endpoint × nonempty eligible module × direction",caveat="Over-representation is exploratory gene-set association, not direct binding or causal pathway activity"),file.path(out,"COMPLETE.json"),auto_unbox=TRUE,pretty=TRUE)
