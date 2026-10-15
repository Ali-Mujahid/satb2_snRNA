args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==4L, args[1]=='--config',args[3]=='--outdir')
cfg <- jsonlite::fromJSON(args[2],simplifyVector=TRUE)
out <- args[4]
if(dir.exists(out)&&length(list.files(out))) stop('No overwrite')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
for(branch in cfg$branches) for(repr in c('unintegrated','harmony')) {
 p<-file.path(cfg$previous_run_root,'results/seurat',branch,'annotate_subsets/CA1_PROS_CANDIDATE',paste0('SNN_',repr,'.rds'))
 g<-readRDS(p);Matrix::writeMM(g,file.path(out,paste0(branch,'_',repr,'.mtx')))
 writeLines(rownames(g),file.path(out,paste0(branch,'_',repr,'_ids.txt')))
 edges<-Matrix::summary(Matrix::tril(g,k=-1L));edges<-edges[edges$x!=0,,drop=FALSE]
 graph<-igraph::make_empty_graph(n=ncol(g),directed=FALSE)
 if(nrow(edges)) graph<-igraph::add_edges(graph,as.vector(t(as.matrix(edges[,c('i','j')]))),attr=list(weight=as.numeric(edges$x)))
 graph<-igraph::set_vertex_attr(graph,'name',value=colnames(g))
 for(iter in c(2L,10L,30L)) {
  ids<-Seurat::RunLeiden(object=graph,leiden_method='igraph',partition.type='RBConfigurationVertexPartition',leiden_objective_function='modularity',resolution.parameter=.4,random.seed=20261013,n.iter=iter)
  names(ids)<-colnames(g)
  ids<-get('GroupSingletons',envir=asNamespace('Seurat'))(ids=ids,SNN=g,group.singletons=TRUE,verbose=FALSE)
  write.table(data.frame(nucleus_id=colnames(g),cluster=as.character(ids)),file.path(out,paste0(branch,'_',repr,'_leiden',iter,'.tsv')),sep='\t',row.names=FALSE,quote=FALSE)
 }
}
jsonlite::write_json(list(status='COMPLETE',stage='fixed_CA1_graph_native_Seurat_iteration_control',n_iter=c(2L,10L,30L),seed=20261013,primary_labels_changed=FALSE),file.path(out,'COMPLETE.json'),auto_unbox=TRUE,pretty=TRUE)
