#!/usr/bin/env Rscript

#Rscript annotate_metafusion_output_with_pfam.R --cluster=/Users/mapostolides/MetaFusion/ISOFORMS/trusight_2019_2020_2020b_all_samples.star_fusion_frame.Jan-8-2021/final.isoform_merge.per_sample.cluster.CANCER_FUSIONS --out=/Users/mapostolides/Sandbox/metafusion.ann.cluster --gene_symbols=/Users/mapostolides/Drive/dowling_lab/MH/hgTables.gene_symbol.ENSG.ENST.ENSP.Nov-13-2020.tsv --pfam_gene_ann=/Users/mapostolides/Drive/fusion-pipeline/trusight/pfam_annotation/domains_per_gene.Rds --pfam_names_descriptions=/Users/mapostolides/Drive/fusion-pipeline/trusight/pfam_annotation/pdb_pfam_mapping.txt

# if (!requireNamespace("limma", quietly = TRUE))
#   BiocManager::install("limma")
# if (!requireNamespace("dplyr", quietly = TRUE))
#   install.packages("dplyr", ,repos = "http://cran.us.r-project.org")
# if (!requireNamespace("stringr", quietly = TRUE))
#   install.packages("stringr", ,repos = "http://cran.us.r-project.org")
#library(limma,quietly = TRUE)
#library(dplyr,quietly = TRUE)
#library(stringr,quietly = TRUE)
suppressMessages(suppressWarnings(library(limma,warn.conflicts = F, quietly = T)))
suppressMessages(suppressWarnings(library(dplyr,warn.conflicts = F, quietly = T)))
suppressMessages(suppressWarnings(library(stringr,warn.conflicts = F, quietly = T)))
#COMMAND LINE ARGUMENTS

args = commandArgs(trailingOnly=TRUE)
#INPUT/OUTPUT FILES
metafusion.cluster <- strsplit(grep('--cluster*', args, value = TRUE), split = '=')[[1]][[2]]
metafusion.cluster.ann <- strsplit(grep('--out*', args, value = TRUE), split = '=')[[1]][[2]]
#print(metafusion.cluster)
#print(metafusion.cluster.ann)

#ANNOTATION FILES
#--gene_symbols=/Users/mapostolides/Drive/dowling_lab/MH/hgTables.gene_symbol.ENSG.ENST.ENSP.Nov-13-2020.tsv --pfam_gene_ann=/Users/mapostolides/Drive/fusion-pipeline/trusight/pfam_annotation/domains_per_gene.Rds --pfam_names_descriptions=/Users/mapostolides/Drive/fusion-pipeline/trusight/pfam_annotation/pdb_pfam_mapping.txt
#ensg.enst.ensp.gene_symbol.file <-"/Users/mapostolides/Drive/dowling_lab/MH/hgTables.gene_symbol.ENSG.ENST.ENSP.Nov-13-2020.tsv"
ensg.enst.ensp.gene_symbol.file  <- strsplit(grep('--gene_symbols*', args, value = TRUE), split = '=')[[1]][[2]]
#domains_per_gene <- readRDS(file="/Users/mapostolides/Drive/fusion-pipeline/trusight/pfam_annotation/domains_per_gene.Rds")
domains_per_gene <-  strsplit(grep('--pfam_gene_ann*', args, value = TRUE), split = '=')[[1]][[2]]
domains_per_gene <- readRDS(domains_per_gene)

#pfam.names.descriptions <- read.csv("/Users/mapostolides/Drive/fusion-pipeline/trusight/pfam_annotation/pdb_pfam_mapping.txt", sep="\t")
pfam.names.descriptions <- strsplit(grep('--pfam_names_descriptions*', args, value = TRUE), split = '=')[[1]][[2]]
pfam.names.descriptions <- read.csv(pfam.names.descriptions, sep="\t")


#metafusion.cluster <- read.csv(file="/Users/mapostolides/MetaFusion/ISOFORMS/trusight_2019_2020_2020b_all_samples.star_fusion_frame.Jan-8-2021/final.isoform_merge.per_sample.cluster.CANCER_FUSIONS", sep="\t", header=T)
metafusion.cluster <- read.csv(file=metafusion.cluster, sep="\t", header=T)



# HELPER FUNCTIONS
#alias2SymbolTable_null <-function(name){
#        symbol<-alias2SymbolTable(name)
#        if(is.na(symbol) ){
#                symbol<- name
#        }
#        return(symbol)
#}
annotate_pfam_IDs_with_names <- function(domains) {
  # takes as input a comma-separated list of PFAM IDs (e.g "PF00664,PF00005")
  # returns annotations with names (e.g. "PF00664|ABC_membrane,PF00005|ABC_tran")
  domains_vect <- c()
  for (domain in unlist(str_split(domains, ","))){
    #print(domain)
    #print(pfam.names.descriptions[which(pfam.names.descriptions$PFAM_ACC == domain),]$PFAM_Name)
    name <-pfam.names.descriptions[which(pfam.names.descriptions$PFAM_ACC == domain),]$PFAM_Name
    
    if( length(name) > 0) {
      domains_vect <- c(domains_vect, paste0(domain, "|", name))
    } else
      domains_vect <- c(domains_vect, domain)
  }
  domains_vect <- paste0(domains_vect, collapse=",")
  #domains_vect <- paste0(sort(domains_vect), collapse=",")
  return(domains_vect)
}

# LOAD REFERENCE FILES

# LOAD ENSG TO GENE_NAME FILE
ensg.enst.ensp.gene_symbol <- read.csv(ensg.enst.ensp.gene_symbol.file, sep="\t")
#rename columns
cols <- c("gene_symbol", "gene_name", "tx_name", "protein_name")
colnames(ensg.enst.ensp.gene_symbol) <- cols

# Load file containing per-gene PFAM domains with ENSG as key

# Load file contain domain names and descriptions
# Remove trailing ".12" numbers from PFAM_ACC
pfam.names.descriptions$PFAM_ACC<-unlist(lapply(str_split(pfam.names.descriptions$PFAM_ACC, pattern="\\."), "[[", 1))
# subset df
pfam.names.descriptions <- unique(pfam.names.descriptions[c("PFAM_ACC", "PFAM_Name", "PFAM_desc")])

#LOAD CLUSTER FILE
#print(colnames(metafusion.cluster))

colnames(metafusion.cluster)[colnames(metafusion.cluster) == "X.gene1"] <- "gene1"

#convert MetaFusion breakpoitns from 0-->1 based
# metafusion.cluster$breakpoint_1 <- metafusion.cluster$breakpoint_1 + 1
# metafusion.cluster$breakpoint_2 <- metafusion.cluster$breakpoint_2 + 1

# rename genes
#metafusion.cluster$gene1 <- unlist(lapply(FUN=alias2SymbolTable_null, metafusion.cluster$gene1))
#metafusion.cluster$gene2 <- unlist(lapply(FUN=alias2SymbolTable_null, metafusion.cluster$gene2))
metafusion.cluster$gene1 <- as.character(metafusion.cluster$gene1)
metafusion.cluster$gene2 <- as.character(metafusion.cluster$gene2) 

metafusion.cluster.full <-metafusion.cluster
#metafusion.cluster<-metafusion.cluster.full[1:4,]

# ANNOTATE FUSIONS WITH PFAM_ID INFO
##NOTE:
# columns of domains_per_gene " "contig"  "start"   "end"     "pfam_id"", 
# "start" will always be smaller than "end". Strand is not considered when assigning
# start and end, always left --> right

metafusion.cluster$domains_kept_gene1 <- "NA"
metafusion.cluster$domains_removed_gene1 <- "NA"
metafusion.cluster$domains_kept_gene2 <- "NA"
metafusion.cluster$domains_removed_gene2 <- "NA"

#START FOR LOOP
for (i in 1:nrow(metafusion.cluster) ){
  fusion <- metafusion.cluster[i,]
  tryCatch(
  {
  domains_kept_gene1 <- ""
  domains_removed_gene1 <- ""
  process_gene1 <- TRUE
  # Get ENSG IDs for fusion
  gene_name1<- unique(ensg.enst.ensp.gene_symbol[ensg.enst.ensp.gene_symbol$gene_symbol == fusion$gene1,]$gene_name)[1]
  # Check whether to process genes
  if(length(gene_name1) ==0 | is.na(fusion$strand1)) { process_gene1 <- FALSE }

  # GENE1
  #check if gene1 ID exists
  if(process_gene1 ) {
    if(gene_name1 %in% names(domains_per_gene)) {

    #convert end from char-->int to integer
    domains_per_gene[[gene_name1]]$end <- as.integer(domains_per_gene[[gene_name1]]$end)

    domains_gene1.df <- domains_per_gene[[gene_name1]]
      if (fusion$strand1 == "+") {# kept are upstream (smaller loci), removed are downstream (larger loci)
        domains_kept_gene1 <- domains_gene1.df[domains_gene1.df$end < fusion$breakpoint_1  ,]
        domains_kept_gene1 <- unique(domains_kept_gene1$pfam_id)
        domains_removed_gene1 <- domains_gene1.df[domains_gene1.df$start > fusion$breakpoint_1  ,]
        domains_removed_gene1 <- unique(domains_removed_gene1$pfam_id)

      } else if (fusion$strand1 == "-"){ #kept are upstream (larger loci), removed are downstream (smaller loci)
        domains_kept_gene1 <- domains_gene1.df[domains_gene1.df$start >  fusion$breakpoint_1  ,]
        domains_kept_gene1 <- unique(domains_kept_gene1$pfam_id)
        domains_removed_gene1 <- domains_gene1.df[domains_gene1.df$end <fusion$breakpoint_1  ,]
        domains_removed_gene1 <- unique(domains_removed_gene1$pfam_id)

      } else {
        print("Ambiguous strand")
      }

    }
  }
  },
  error=function(cond) {
    #message(paste0("Annotation failed for fusion", paste0(gene_name1, "\t", gene_name2)))
    message(paste0("PFAM Annotation failed for gene1 ",fusion$gene1, ":", gene_name1,  " of fusion: ", paste0(fusion$gene1, ":", gene_name1, "\\t", fusion$gene2, ":", gene_name2)))
    #metafusion.cluster[i,] <- fusion
    return(NA)
  })
  
  # Make empty domain fields "NA"
  fusion$domains_kept_gene1 <- dplyr::na_if(paste0(domains_kept_gene1 , collapse=","), "")
  fusion$domains_removed_gene1 <- dplyr::na_if(paste0(domains_removed_gene1, collapse=","), "")
  tryCatch(
    {
  domains_kept_gene2 <- ""
  domains_removed_gene2 <- ""
  process_gene2 <- TRUE
  gene_name2<- unique(ensg.enst.ensp.gene_symbol[ensg.enst.ensp.gene_symbol$gene_symbol == fusion$gene2,]$gene_name)[1]
  if(length(gene_name2) ==0| is.na(fusion$strand2)) { process_gene2 <- FALSE }

  #GENE2
  if(process_gene2 ) {
    if(gene_name2 %in% names(domains_per_gene)) {
     #print(paste0("gene_name2: ", ))
    #convert end from char-->int to integer
    #print(as.integer(domains_per_gene[[gene_name2]]$end))
    domains_per_gene[[gene_name2]]$end <- as.integer(domains_per_gene[[gene_name2]]$end)
     #print(paste0("gene_name2: TEST", gene_name2))

    domains_gene2.df <- domains_per_gene[[gene_name2]]
      if (fusion$strand2 == "-") {# kept are downstream (smaller loci), removed are upstream (larger loci)
        domains_kept_gene2 <- domains_gene2.df[domains_gene2.df$end < fusion$breakpoint_2  ,]
        domains_kept_gene2 <- unique(domains_kept_gene2$pfam_id)
        domains_removed_gene2 <- domains_gene2.df[domains_gene2.df$start > fusion$breakpoint_2  ,]
        domains_removed_gene2 <- unique(domains_removed_gene2$pfam_id)

      }else if (fusion$strand2 == "+"){# Kept are downstream (larger loci), removed are upstream (smaller loci)
        domains_kept_gene2 <- domains_gene2.df[domains_gene2.df$start >  fusion$breakpoint_2  ,]
        domains_kept_gene2 <- unique(domains_kept_gene2$pfam_id)
        domains_removed_gene2 <- domains_gene2.df[domains_gene2.df$end <  fusion$breakpoint_2  ,]
        domains_removed_gene2 <- unique(domains_removed_gene2$pfam_id)

      } else {
        print("Ambiguous strand")
      }
    }
  }
  # Make empty domain fields "NA"
  fusion$domains_kept_gene2 <- dplyr::na_if(paste0(domains_kept_gene2 , collapse=","), "")
  fusion$domains_removed_gene2 <- dplyr::na_if(paste0(domains_removed_gene2, collapse=","), "")

  metafusion.cluster[i,] <- fusion
  
},
        error=function(cond) {
            #message(paste0("Annotation failed for fusion", paste0(gene_name1, "\t", gene_name2)))
            message(paste0("PFAM Annotation failed for gene1 ",fusion$gene1, ":", gene_name1,  " of fusion: ", paste0(fusion$gene1, ":", gene_name1, "\\t", fusion$gene2, ":", gene_name2)))
            #metafusion.cluster[i,] <- fusion
            return(NA)
        })
  
  #END FOR LOOP
}

# Add domain names to PFAM IDs
metafusion.cluster$domains_kept_gene1 <- unlist(lapply(metafusion.cluster$domains_kept_gene1, FUN=annotate_pfam_IDs_with_names))
metafusion.cluster$domains_kept_gene2 <- unlist(lapply(metafusion.cluster$domains_kept_gene2, FUN=annotate_pfam_IDs_with_names))
metafusion.cluster$domains_removed_gene1 <- unlist(lapply(metafusion.cluster$domains_removed_gene1, FUN=annotate_pfam_IDs_with_names))
metafusion.cluster$domains_removed_gene2 <- unlist(lapply(metafusion.cluster$domains_removed_gene2, FUN=annotate_pfam_IDs_with_names))

# write cluster to file
colnames(metafusion.cluster)[colnames(metafusion.cluster) == "X.cluster_type"] <- "#cluster_type"
write.table(metafusion.cluster , file=metafusion.cluster.ann, quote=FALSE, sep='\t', append=F, col.names = T, row.names=F)



