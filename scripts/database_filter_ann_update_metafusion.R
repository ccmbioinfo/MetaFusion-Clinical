#library(RSQLite)
#library(dplyr)

suppressMessages(suppressWarnings(library(dplyr,warn.conflicts = F, quietly = T)))
suppressMessages(suppressWarnings(library(RSQLite,warn.conflicts = F, quietly = T)))

#COMMAND LINE ARGUMENTS
args = commandArgs(trailingOnly=TRUE)
#INPUT/OUTPUT FILES
cluster <- strsplit(grep('--cluster*', args, value = TRUE), split = '=')[[1]][[2]]
cluster.filt <- strsplit(grep('--out*', args, value = TRUE), split = '=')[[1]][[2]]
update.hist <- as.integer(strsplit(grep('--update_hist*', args, value = TRUE), split = '=')[[1]][[2]])
database.path <- strsplit(grep('--database*', args, value = TRUE), split = '=')[[1]][[2]]
#print(database.path)

#print(cluster)
#print(cluster.filt)
print(update.hist)
#q()

# Initialize database:
#conn <- dbConnect(RSQLite::SQLite(), "/MetaFusion/reference_files/historical_fusions.db")
conn <- dbConnect(RSQLite::SQLite(), database.path) 

#conn <- dbConnect(RSQLite::SQLite(), "/Users/mapostolides/MetaFusion.clinical/reference_files/historical_fusions.db")
#conn <- dbConnect(RSQLite::SQLite(), "/hpf/largeprojects/ccmbio/mapostolides/MetaFusion.clinical/reference_files/historical_fusions.db")


false_positives <- dbGetQuery(conn, "SELECT * from false_positives")
historical_fusions <- dbGetQuery(conn, "SELECT * from historical_fusions")
# rename cols
colnames(false_positives) <- paste("FP_", colnames(false_positives), sep="")

# Metafusion total result file
cluster.orig <- cluster <- read.csv(file=cluster, sep="\t", header=T)
colnames(cluster)[colnames(cluster) == "X.cluster_type"] <- "cluster_type"

#Remove false positives
cluster <- anti_join(cluster, false_positives, by=c("gene1"="FP_gene1", "gene2"="FP_gene2"), keep=TRUE)
cluster <- anti_join(cluster, false_positives, by=c("gene2"="FP_gene1", "gene1"="FP_gene2"), keep=TRUE)
# clear false_positives
#dbClearResult(false_positives)

print(paste0("Cluster rows: ",nrow(cluster.orig), " filtered rows: ", nrow(cluster)))

query.clin <- 1
if (query.clin == 1){

  # GENERATE CLINICAL KEEP LIST
  print("GENERATE CLINICAL KEEP LIST")
  # Set "prev_samples" and "Num_prev_samples" fields to "NA"
  cluster$clinical_samples <- "NA"
  cluster$num_clinical_samples <- "NA"
  
  for(i in 1:nrow(cluster)) {
    fusion <- cluster[i, ]
  
    #genes query: Run only if the fusion is not an ITD (gene1 == gene2)
    if (fusion$gene1 == fusion$gene2) {
      query.genes <- paste0("SELECT * FROM clinical_fusions WHERE (gene1='",fusion$gene1,"' AND gene2='", fusion$gene2, "')")
      query.genes.rev <- paste0("SELECT * FROM clinical_fusions WHERE (gene1='",fusion$gene2,"' AND gene2='", fusion$gene1, "')")
      df.genes <- dbGetQuery(conn, query.genes )
      df.genes.rev <- dbGetQuery(conn, query.genes.rev )
      df.genes <-rbind(df.genes, df.genes.rev)
    }
  
    # breakpoints query
    slop <- 10
    query.bps <- paste0("SELECT * FROM clinical_fusions WHERE (ref_pos1 > ",fusion$breakpoint_1-slop," AND ref_pos1 < ",     fusion$breakpoint_1+slop, " AND ref_pos2 > ",fusion$breakpoint_2-slop," AND ref_pos2 < ", fusion$breakpoint_2+slop, " AND chr1='",fusion$chr1,"' AND chr2='", fusion$chr2, "'", ") ")
    query.bps.rev <- paste0("SELECT * FROM clinical_fusions WHERE (ref_pos1 > ",fusion$breakpoint_2-slop," AND ref_pos1 < ",     fusion$breakpoint_2+slop, " AND ref_pos2 > ",fusion$breakpoint_1-slop," AND ref_pos2 < ", fusion$breakpoint_1+slop, " AND chr1='",fusion$chr2,"' AND chr2='", fusion$chr1, "'", ") ")
    df.bps <- dbGetQuery(conn, query.bps )
    df.bps.rev <- dbGetQuery(conn, query.bps.rev )
    df.bps <- rbind(df.bps, df.bps.rev)
 
    # Remove duplicates. If an ITD, we don't have gene query.
    if (fusion$gene1 == fusion$gene2) {
      df.matches <- unique(rbind(df.genes, df.bps))
    } else {
      df.matches <- unique(df.bps)
    }

    # remove spaces from sample names
    df.matches$sample <- stringr::str_replace(pattern="\\s+", replacement="-", df.matches$sample)
  
    # change to vector, since randomly becoming factor
    fusion$samples <- as.vector(fusion$samples)
  
    samples <- df.matches$sample
	# remove spaces from sample name
    samples <- gsub(pattern=" +", replacement="", samples)

    # Add "clinical_samples" and "num_clinical_samples" fields
    fusion$clinical_samples <- dplyr::na_if(paste0(unique(samples )   , collapse=","), "")
    fusion$num_clinical_samples <- length(unique(samples))
  
    #Update fusion in cluster
    cluster[i, ] <- fusion
  
  }

}

query.hist <- 1
if (query.hist == 1){

  # QUERY DATABASE FOR PREVIOUSLY SEEN FUSIONS
  print("QUERY DATABASE FOR PREVIOUSLY SEEN FUSIONS")
  # Set "prev_samples" and "num_prev_samples" fields to "NA"
  cluster$prev_samples <- "NA"
  cluster$num_prev_samples <- "NA"
  
  for(i in 1:nrow(cluster)) {
    #i=43
    fusion <- cluster[i, ]
  
    # genes query
    query.genes <- paste0("SELECT * FROM historical_fusions WHERE (gene1='",fusion$gene1,"' AND gene2='", fusion$gene2, "')")
    query.genes.rev <- paste0("SELECT * FROM historical_fusions WHERE (gene1='",fusion$gene2,"' AND gene2='", fusion$gene1, "')")
    df.genes <-rbind(dbGetQuery(conn, query.genes ), dbGetQuery(conn, query.genes.rev ))
  
    # breakpoints query
    slop <- 10
    query.bps <- paste0("SELECT * FROM historical_fusions WHERE (breakpoint_1 > ",fusion$breakpoint_1-slop," AND breakpoint_1 < ",     fusion$breakpoint_1+slop, " AND breakpoint_2 > ",fusion$breakpoint_2-slop," AND breakpoint_2 < ", fusion$breakpoint_2+slop," AND chr1='",fusion$chr1,",' AND chr2='", fusion$chr2, "'",  ") ")
    query.bps.rev <- paste0("SELECT * FROM historical_fusions WHERE (breakpoint_1 > ",fusion$breakpoint_2-slop," AND breakpoint_1 < ",     fusion$breakpoint_2+slop, " AND breakpoint_2 > ",fusion$breakpoint_1-slop," AND breakpoint_2 < ", fusion$breakpoint_1+slop," AND chr1='",fusion$chr2,",' AND chr2='", fusion$chr1, "'", ") " )
    df.bps <- rbind(dbGetQuery(conn, query.bps ), dbGetQuery(conn, query.bps.rev ))
  
    #matches are hits based on either gene names, or where both breakpoints are within 10bp
    df.matches <- unique(rbind(df.genes, df.bps))
  
    # change to vector, since randomly becoming factor
    fusion$samples <- as.vector(fusion$samples)
    # remove "self" sample from previously seen
    samples <- df.matches$samples[!df.matches$samples %in% c(fusion$samples)]
	# remove random garbage after SR-1234-GBM format
	# remove spaces 
	samples <- gsub(pattern=" +", replacement="", samples)
    # Add "prev_samples" and "num_prev_samples" fields
    fusion$prev_samples <- dplyr::na_if(paste0(unique(samples )   , collapse=","), "")
    fusion$num_prev_samples <- length(unique(samples))
  
  
    cluster[i, ] <- fusion
  }
}

# UPDATE HISTORICAL FUSIONS TABLE
if (update.hist == 1){
  cluster.update <- cluster[!names(cluster) %in% c("prev_samples", "num_prev_samples", "clinical_samples", "num_clinical_samples", "junction_sequence")]
  # Create unique index (one time operation)
  #dbSendQuery(conn, "CREATE UNIQUE INDEX unique_idx ON historical_fusions(gene1,gene2,tools,samples,chr1,breakpoint_1,chr2,breakpoint_2);" )

  # Clear the stage, add current run to stage
  print("Clear the stage, add current run to stage")
  del_res <- dbSendQuery(conn,"delete from stage;")
  print(del_res)
  dbClearResult(del_res)
  print("dbAppendTable(conn, stage, cluster.update)")
  dbAppendTable(conn, "stage", cluster.update)
  #append_res <- dbAppendTable(conn, "stage", cluster.update)
  #dbClearResult(append_res)
  #stage <- RSQLite::dbGetQuery(conn, "SELECT * FROM stage;")
  # Only insert cluster entries not seen before, based on fields in unique_idx index
  #dbSendQuery(conn, "insert or ignore into historical_fusions select * from stage;")
  print("dbSendQuery(conn, insert or ignore into historical_fusions select * from stage;)")
  insert_res <- dbSendQuery(conn, "insert or ignore into historical_fusions select * from stage;")
  print(insert_res)
  dbClearResult(insert_res)
}

#q()
# WRITE TO FILE AND DISCONNECT
colnames(cluster)[colnames(cluster) == "cluster_type"] <- "#cluster_type"
write.table(cluster , file=cluster.filt, quote=FALSE, sep='\t', append=F, col.names = T, row.names=F, na = "NA")
dbDisconnect(conn)
