#!/bin/bash
# Script to update "clinical_fusions" or "false_positives" tables

# Absolute path to database 
database=/Users/maposto/MetaFusion-Clinical/MWE/historical_database.db

# Absolute path to excel spreadsheet with selected calls 
# Only needed if doing "add" operation. Can use placeholder path for "view"
excel=/Users/maposto/MetaFusion-Clinical/MWE/run1.out/final.n2.cluster.filt.xlsx
#excel=/placeholder/file/path/nothing.xlsx

# Operation "add", "view" 
operation=add
#operation=view

# Table to view or update database tables 
table=clinical_fusions
#table=false_positives
#table=historical_fusions

# Run update script
scripts=/Users/maposto/MetaFusion-Clinical/scripts
cd  /Users/maposto/MetaFusion-Clinical/MWE

Rscript $scripts/update_fusion_database_tables.R --excel=$excel --database=$database --table=$table --operation=$operation

cd /Users/maposto/MetaFusion-Clinical/scripts
