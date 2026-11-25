# Run from the local analysis folder containing sources/.
library(DBI)
library(RSQLite)
library(dplyr)
source("sources/synaptome.db/R/geneInfo.R")
source("sources/synaptome.db/R/ppiTables.R")
db <- dbConnect(SQLite(),Sys.glob("sources/*.db.sqlite"),flags=SQLITE_RO)
get_dbconn <- function() db
pairs <- read.delim("../SATB2_CA1_EPHYS_PLANNING_20261125/strict_DEG_ligand_receptor_pairs.tsv")
genes <- sort(unique(c(pairs$ligand,pairs$receptor)))
hits <- getGeneInfoByName(genes)
hits <- hits[hits$MouseName %in% genes,]
write.table(hits,"native_synaptome_R_lookup.tsv",quote=FALSE,sep="\t",row.names=FALSE)
edges <- getPPIbyIDs(unique(hits$GeneID),type="limited")
write.table(edges,"native_synaptome_R_ppi.tsv",quote=FALSE,sep="\t",row.names=FALSE)
dbDisconnect(db)
