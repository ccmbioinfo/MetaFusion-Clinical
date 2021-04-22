#!/bin/bash
# Script to update "clinical_fusions" or "false_positives" tables

# Absolute path to database 
database=/hpf/largeprojects/ccmbio/mapostolides/Sandbox/update_db_test/historical_database_from_rob.reformat.TEST_SCRAP.db

# Absolute path to excel spreadsheet with selected calls 
# Only needed if doing "add" operation. Can use placeholder path for "view"
#excel=/hpf/largeprojects/ccmbio/mapostolides/Sandbox/update_db_test/test_view/false_positives.UPDATE.xlsx
excel=/hpf/largeprojects/ccmbio/mapostolides/Sandbox/update_db_test/test_view/clinical_fusions.UPDATE.xlsx
#excel=/placeholder/file/path/nothing.xlsx

# Operation "add", "view" 
#operation=add
operation=view

# Table to view or update "clinical_fusions" or "false_positives"
#table=clinical_fusions
#table=false_positives
table=historical_fusions

# Do not modify below this line (unless you know what you're doing!)
module load Singularity
db_dir=$(dirname $database)
scripts=/hpf/largeprojects/ccmbio/mapostolides/MetaFusion.clinical/scripts
singularity exec -B $db_dir -B $PWD -B $scripts -B $excel /hpf/largeprojects/ccmbio/mapostolides/MetaFusion.clinical/MetaFusion.clinical.simg Rscript $scripts/update_fusion_database_tables.R --excel=$excel --database=$database --table=$table --operation=$operation
