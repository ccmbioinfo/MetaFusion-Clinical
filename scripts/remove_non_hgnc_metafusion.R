if (!requireNamespace("limma", quietly = TRUE))
  BiocManager::install("limma")
if (!requireNamespace("AnnotationDbi", quietly = TRUE))
  BiocManager::install("AnnotationDbi")
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
  BiocManager::install("org.Hs.eg.db")

if (!requireNamespace("dplyr", quietly = TRUE))
  BiocManager::install("dplyr")
if (!requireNamespace("stringr", quietly = TRUE))
  install.packages("stringr")

library(limma)
library(dplyr)
library(stringr)

#COMMAND LINE ARGUMENTS
args = commandArgs(trailingOnly=TRUE)
#INPUT/OUTPUT FILES
cluster <- strsplit(grep('--cluster*', args, value = TRUE), split = '=')[[1]][[2]]
cluster.non_hgnc <- strsplit(grep('--out*', args, value = TRUE), split = '=')[[1]][[2]]
print(cluster)
print(cluster.non_hgnc)

# Metafusion total result file
cluster <- read.csv(file=cluster, sep="\t", header=T)
colnames(cluster)[colnames(cluster) == "X.cluster_type"] <- "cluster_type"


#REMOVE NON-HGNC SYMBOLS FROM CLUSTER

remove_non_HGNC <- function(genes){
  # input: genes (comma-separated string of gene names)
  genes<- unlist(str_split(genes, pattern=","))
  # check whether an HGNC symbol remains
  genes.hgnc <- genes[!is.na(limma::alias2SymbolTable(genes))]
  genes.hgnc <- paste(genes.hgnc, collapse=",")
  if (nchar(genes.hgnc) ){
    return(genes.hgnc)
  } else {
    return(paste(genes, collapse=","))
  }
}
cluster$gene1 <- unlist(lapply(cluster$gene1, remove_non_HGNC))
cluster$gene2 <- unlist(lapply(cluster$gene2, remove_non_HGNC))



# write cluster to file
colnames(cluster)[colnames(cluster) == "cluster_type"] <- "#cluster_type"
write.table(cluster , file=cluster.non_hgnc, quote=FALSE, sep='\t', append=F, col.names = T, row.names=F)



