#!/bin/bash

#DATABASE
database=/Users/maposto/MetaFusion-Clinical/MWE/historical_database.db

fusiontools=/Users/maposto/MetaFusion-Clinical/scripts
#REFERENCE FILES FILES
mkdir $runs_dir
ref_dir=/Users/maposto/MetaFusion-Clinical/reference_files
gene_bed=$ref_dir/new_bed.total.Oct-1-2020.uniq.bed
gene_info=$ref_dir/Homo_sapiens.gene_info
genome_fasta=$ref_dir/human_g1k_v37_decoy.fasta
recurrent_bedpe=$ref_dir/blocklist_breakpoints.bedpe

outdir=/Users/maposto/MetaFusion-Clinical/MWE/run2
echo generating output in $outdir
mkdir $outdir
cff=/Users/maposto/MetaFusion-Clinical/MWE/run2.cff

bash $fusiontools/MetaFusion.clinical.sh --outdir $outdir \
                 --cff $cff  \
                 --gene_bed $gene_bed \
                 --annotate_exons \
                 --fusion_annotator \
                 --genome_fasta $genome_fasta \
                 --gene_info $gene_info \
                 --num_tools=2  \
                 --per_sample \
                 --recurrent_bedpe $recurrent_bedpe \
                 --scripts $fusiontools \
                 --database $database \
                 --update_hist \
                 --ref_dir $ref_dir

