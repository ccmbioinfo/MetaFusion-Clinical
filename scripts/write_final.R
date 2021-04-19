#!/usr/bin/env Rscript
#Rscript write_final.R --cluster=/Users/mapostolides/MetaFusion.clinical/CLINICAL_RUNS/trusight_sample_subset_test.Mar-16-2021/final.n2.cluster

suppressMessages(suppressWarnings(library(writexl,warn.conflicts = F, quietly = T)))
suppressMessages(suppressWarnings(library(stringr,warn.conflicts = F, quietly = T)))
suppressMessages(suppressWarnings(library(dplyr,warn.conflicts = F, quietly = T)))
#COMMAND LINE ARGUMENTS

args = commandArgs(trailingOnly=TRUE)
#INPUT/OUTPUT FILES
cluster <- strsplit(grep('--cluster*', args, value = TRUE), split = '=')[[1]][[2]]
#cluster <- "/Users/mapostolides/MetaFusion.clinical/CLINICAL_RUNS/trusight_sample_subset_test.Mar-16-2021/final.n2.cluster"

print(cluster)

metafusion.cluster <- read.csv(cluster, sep="\t", header=T)
colnames(metafusion.cluster)[colnames(metafusion.cluster) == "X.cluster_type"] <- "cluster_type"
colnames(metafusion.cluster)

#rearr_cols <- c( "gene1", "gene2", "num_tools", "max_split_cnt", "max_span_cnt", "frame", "cancer_db_hits", "samples",  "chr1", "breakpoint_1", "chr2", "breakpoint_2",  "inferred_fusion_type", "disease", "tools", "rna_type1", "rna_type2", "strand1", "strand2", "clinical_samples", "num_clinical_samples", "prev_samples", "num_prev_samples", "junction_sequence", "domains_kept_gene1", "domains_removed_gene1", "domains_kept_gene2", "domains_removed_gene2", "gene1_on_bnd", "gene1_close_to_bnd", "gene2_on_bnd", "gene2_close_to_bnd", "exon1", "exon2", "sample_type", "fusion_IDs", "cluster_type")
rearr_cols <- c( "gene1", "gene2", "num_tools", "max_split_cnt", "max_span_cnt", "frame", "cancer_db_hits", "samples",  "chr1", "breakpoint_1", "chr2", "breakpoint_2",  "inferred_fusion_type", "disease", "tools", "rna_type1", "rna_type2", "strand1", "strand2", "clinical_samples", "num_clinical_samples", "prev_samples", "num_prev_samples", "junction_sequence", "domains_kept_gene1", "domains_removed_gene1", "domains_kept_gene2", "domains_removed_gene2", "gene1_on_bnd", "gene1_close_to_bnd", "gene2_on_bnd", "gene2_close_to_bnd", "exon1", "exon2", "sample_type", "fusion_IDs")
metafusion.cluster <- metafusion.cluster[ ,rearr_cols]

# ADD IN LAST FILTERS HERE IF NEEDED
# Remove non-clinical low-coverage fusions
cluster.filt <- metafusion.cluster %>% filter( !((max_split_cnt + max_span_cnt) < 3 & num_clinical_samples == 0) )
# Remove non-clinical truncated fusions
cluster.filt <- cluster.filt %>% filter( !(grepl("Truncated", inferred_fusion_type) & num_clinical_samples == 0) )
# Remove non-clinical fusions which are highly recurrent
cluster.filt <- cluster.filt %>% filter( !(num_clinical_samples == 0 &  num_prev_samples >= 5) )

# write cluster to file
#colnames(metafusion.cluster)[colnames(metafusion.cluster) == "cluster_type"] <- "#cluster_type"
#colnames(cluster.filt)[colnames(cluster.filt) == "cluster_type"] <- "#cluster_type"

# Write filtered output
print(paste0(cluster, ".filt"))
print(paste0(cluster, ".filt.xlsx"))
write_xlsx(cluster.filt, path = paste0(cluster, ".filt.xlsx"))
write.table(cluster.filt , file=paste0(cluster, ".filt"), quote=FALSE, sep='\t', append=F, col.names = T, row.names=F, na = "NA")

# Write regular output
print(paste0(cluster, ".xlsx"))
write_xlsx(metafusion.cluster, path = paste0(cluster, ".xlsx"))



