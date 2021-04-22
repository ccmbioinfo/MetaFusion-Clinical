#Rscript update_fusion_database_tables.R --table=clinical_fusions --operation=add --excel=/Users/mapostolides/Sandbox/test_xlsx_metafusion/final.n2.cluster.xlsx --database=/Users/mapostolides/MetaFusion.clinical/reference_files/historical_fusions.db

library(RSQLite)
library(dplyr)
library(stringr)
library(writexl)

#COMMAND LINE ARGUMENTS
args = commandArgs(trailingOnly=TRUE)
#INPUT/OUTPUT FILES
cluster.xlsx <- strsplit(grep('--excel*', args, value = TRUE), split = '=')[[1]][[2]]
database.path <- strsplit(grep('--database*', args, value = TRUE), split = '=')[[1]][[2]]
table <- strsplit(grep('--table*', args, value = TRUE), split = '=')[[1]][[2]]
operation <- strsplit(grep('--operation*', args, value = TRUE), split = '=')[[1]][[2]]

print(cluster.xlsx)
print(database.path)
print(table)
print(operation )

# Get and format curent time
current_time <- Sys.time()
current_time <- str_replace_all(current_time, c(" " = "." , "," = "" ))

# Open database connection
#conn <- dbConnect(RSQLite::SQLite(), "/Users/mapostolides/MetaFusion.clinical/reference_files/historical_fusions.db")
#conn <- dbConnect(RSQLite::SQLite(), "/hpf/largeprojects/ccmbio/mapostolides/MetaFusion.clinical/reference_files/historical_fusions.db")
conn <- dbConnect(RSQLite::SQLite(), database.path)
print("connection to database established")
clinical_fusions <- dbGetQuery(conn, "SELECT * from clinical_fusions")
false_positives <- dbGetQuery(conn, "SELECT * from false_positives")



if (operation == "remove"){
  print("may need to implement this later. Not sure. Right now 'remove' operation does not work.")
}

if(table == "clinical_fusions"){
  # Read in excel
  cluster <- readxl::read_xlsx(cluster.xlsx)
  # Update clinical_fusions table
  if(operation == "add"){
    print("updating clinical_fusions table")
    # Subset and format cluster to make compatible with clincal_fusions table
    clinical.colnames <- c("comment", "sample", "fusion", "gene1", "gene2", "chr1", "ref_pos1", 
    "strand1", "chr2", "ref_pos2", "strand2", "split_cnt", "span_cnt")
    cluster.subset <- cluster[c("tools", "samples",  "gene1",  "gene1", "gene2",  "chr1", "breakpoint_1", 
                            "strand1","chr2","breakpoint_2", "strand2",  "max_split_cnt", "max_span_cnt"  ) ]
    names(cluster.subset) <- clinical.colnames
    cluster.subset$fusion <- paste0(cluster.subset$gene1, "--", cluster.subset$gene2)
    
    # Clear the stage
    print("Clear the stage")
    del_res <- dbSendQuery(conn,"delete from clinical_stage;")
    print(del_res)
    dbClearResult(del_res)
    
    # Append to stage
    print("dbAppendTable(conn, clinical_stage, cluster.subset)")
    dbAppendTable(conn, "clinical_stage", cluster.subset)
    dbGetQuery(conn, "SELECT * from clinical_stage")
  
    print("dbSendQuery(conn, insert or ignore into clinical_fusions select * from clinical_stage;)")
    insert_res <- dbSendQuery(conn, "insert or ignore into clinical_fusions select * from clinical_stage;")
    print(insert_res)
    dbClearResult(insert_res)
    
    # assign updated
    clinical_fusions.updated <- dbGetQuery(conn, "SELECT * from clinical_fusions")
    
    # WRITE EXCEL FILES
    write_xlsx(clinical_fusions, path = paste0("clinical_fusions.", current_time, ".xlsx"))
    write_xlsx(clinical_fusions.updated, path = paste0("clinical_fusions.updated.", current_time, ".xlsx"))
    
  } else if (operation == "remove"){
    print("remove not yet implemented")
  } else if (operation == "view") {
    write_xlsx(clinical_fusions, path = paste0("clinical_fusions.view.", current_time, ".xlsx"))
  }
  
# UPDATE FALSE POSITIVE TABLE
} else if(table == "false_positives"){
  # Read in excel
  FP.cluster <- readxl::read_xlsx(cluster.xlsx)
  
  # Subset and format cluster to make compatible with false_positives table
  FP.cluster.subset <- cluster[c("gene1", "gene2") ]
  
  # Clear the stage
  print("Clear the FP_stage")
  del_res <- dbSendQuery(conn,"delete from FP_stage;")
  print(del_res)
  dbClearResult(del_res)
  
  # Append to stage
  print("dbAppendTable(conn, FP_stage, FP.cluster.subset)")
  dbAppendTable(conn, "FP_stage", FP.cluster.subset)
  FP_stage <- dbGetQuery(conn, "SELECT * from FP_stage")
  
  # Update false_positives table
  if(operation == "add"){
    print("dbSendQuery(conn, insert or ignore into false_positives select * from FP_stage;)")
    insert_res <- dbSendQuery(conn, "insert or ignore into false_positives select * from FP_stage;")
    print(insert_res)
    dbClearResult(insert_res)
    
    # assign updated
    false_positives.updated <- dbGetQuery(conn, "SELECT * from false_positives")
    
    # WRITE EXCEL FILES
    write_xlsx(false_positives, path = paste0("false_positives.", current_time, ".xlsx"))
    write_xlsx(false_positives.updated, path = paste0("false_positives.updated.", current_time, ".xlsx"))
    
  } else if (operation == "remove"){
    print("remove not yet implemented")
  } else if (operation == "view") {
    write_xlsx(false_positives, path = paste0("false_positives.view.", current_time, ".xlsx"))
  }
  
} else if(table == "historical_fusions"){
  historical_fusions <- dbGetQuery(conn, "SELECT * from historical_fusions")
  if(operation == "view"){
    print("viewing historical_fusions table")
    write_xlsx(historical_fusions, path = paste0("historical_fusions.view.", current_time, ".xlsx"))
    
  } else if(operation == "add"){
    print("add to historical_fusions table is done by running MetaFusion, and update is automatic. 
          This is not the correct way to update the historical_fusions table")
  } else if (operation == "remove"){
    print("remove not yet implemented. It also probably should not be implemented in this way")
  }
}
  
#  DISCONNECT
dbDisconnect(conn)


# write updated table to xlsx
# if(table == "clinical_fusions"){
#   write_xlsx(clinical_fusions, path = paste0("clinical_fusions.", current_time, ".xlsx"))
#   write_xlsx(clinical_fusions.updated, path = paste0("clinical_fusions.updated.", current_time, ".xlsx"))
# } else if(table == "false_positives"){
#   write_xlsx(false_positives, path = paste0("false_positives.", current_time, ".xlsx"))
#   write_xlsx(false_positives.updated, path = paste0("false_positives.updated.", current_time, ".xlsx"))
#   } #else if (table == "historical_fusions")

# Code to format FP and clinical table and stage
format <- 0 
if (format){
  #FALSE POSITIVES
  # create FP_stage
  res <- dbWriteTable(conn, "FP_stage", false_positives, overwrite=T)
  # Create unique index (one time operation)
  res<- dbSendQuery(conn, "CREATE UNIQUE INDEX unique_idx_FP ON false_positives(gene1,gene2);" )
  dbClearResult(res)
  print("created FP_stage and created unique index for table false_positives")
  #CLINICAL
  # re-write clinical_fusions table with split_span_cnt
  clinical.df<- readRDS( file="/hpf/largeprojects/ccmbio/mapostolides/Sandbox/update_db_test/test_rob_update/clinical.df.all_samples.splt_spn_cnt.Mar-19-2021.Rds")
  res <- dbWriteTable(conn, "clinical_fusions", clinical.df, overwrite=T)
  # create clinical_stage
  clinical_fusions <- dbGetQuery(conn, "SELECT * from clinical_fusions")
  res <- dbWriteTable(conn, "clinical_stage", clinical_fusions, overwrite=T)
  res<- dbSendQuery(conn, "CREATE UNIQUE INDEX unique_idx_clinical ON clinical_fusions(comment,sample,fusion,gene1,gene2,chr1,ref_pos1,strand1,chr2,ref_pos2,strand2);" )
  dbClearResult(res)
  dbDisconnect(conn)
  print("created clinical_stage and created unique index for table clinical_fusions")
  # Write clinical_fusions and false_positives tables
  current_time <- Sys.time()
  current_time <- str_replace_all(current_time, c(" " = "." , "," = "" ))
  write_xlsx(clinical_fusions, path = paste0("clinical_fusions.", current_time, ".xlsx"))
  write_xlsx(false_positives, path = paste0("false_positives.", current_time, ".xlsx"))
  q()
}



