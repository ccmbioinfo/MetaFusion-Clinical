#!/bin/bash

#Change date to current date. Can also add tag to this string for multiple runs
date=Feb-17-2021
date=Feb-16-2021
date=Feb-22-2021.re-run
date=Feb-23-2021
date=Mar-16-2021
#date=Mar-29-2021
date=Mar-30-2021
date=Apr-6-2021
date=Apr-12-2021

#DATABASE
database=/MetaFusion/reference_files/historical_fusions.db
#DATASETS
brca_4=0
trusight_frame=0
sim45_sim52=0
trusight_2019_2020_2020b=0
trusight_2019_2020_2020b_2020c=0
trusight_2019_2020_2020b_2020c_dragen=0
ntrk_control=1
trusight_missed=0
trusight_sample_subset_test=0
trycatch_test=0

fusiontools=/MetaFusion/scripts
#REFERENCE FILES FILES
runs_dir=/MetaFusion/CLINICAL_RUNS
mkdir $runs_dir
gene_bed=/MetaFusion/reference_files/new_bed.total.Oct-1-2020.uniq.bed
gene_info=/MetaFusion/reference_files/Homo_sapiens.gene_info
genome_fasta=/MetaFusion/reference_files/human_g1k_v37_decoy.fasta
recurrent_bedpe=/MetaFusion/reference_files/blocklist_breakpoints.bedpe
ref_dir=/MetaFusion/reference_files


if [ $trusight_2019_2020_2020b -eq 1 ]; then
echo TRUSIGHT 2019+2020+2020b combined
#outdir=$runs_dir/trusight_2019_2020_2020b_all_samples.$date
outdir=$runs_dir/trusight_2019_2020_2020b_all_samples.star_fusion_frame.$date
echo generating output in $outdir
mkdir $outdir
#cff=/MetaFusion/test_data/trusight_cff/2019.2020.2020b.ALL_SAMPLES_ROB.merged.cff
cff=/MetaFusion/test_data/trusight_cff/2019.2020.2020b.trusight.star_fusion_frame.cff

bash MetaFusion.clinical.sh --outdir $outdir \
                 --cff $cff  \
                 --gene_bed $gene_bed \
                 --annotate_exons \
                 --fusion_annotator \
                 --genome_fasta $genome_fasta \
                 --gene_info $gene_info \
                 --num_tools=2  \
                 --per_sample \
                 --recurrent_bedpe $recurrent_bedpe \
                 --scripts $fusiontools
fi

# SIM45.SIM52.combined
if [ $sim45_sim52 -eq 1 ]; then
echo SIM45.SIM52
#outdir=$runs_dir/SIM45.SIM52.benchmark.$date
outdir=$runs_dir/SIM45.SIM52.benchmark.arriba.star_fusion.star_seqr.$date
echo generating output in $outdir
mkdir $outdir
#cff=/MetaFusion/test_data/cff/dream.sim45.sim52.cff
cff=/MetaFusion/test_data/cff_test/dream.arriba.star_fusion.star_seqr.cff
truth_fusions=/MetaFusion/test_data/truth_sets/dream.sim45.sim52.truth_set.dat

bash MetaFusion.clinical.sh --outdir $outdir \
                 --cff $cff  \
                 --gene_bed $gene_bed \
                 --annotate_exons \
                 --fusion_annotator \
                 --genome_fasta $genome_fasta \
                 --gene_info $gene_info \
                 --num_tools=2  \
                 --per_sample \
                 --recurrent_bedpe $recurrent_bedpe \
				 --truth_set $truth_fusions \
                 --scripts $fusiontools

fi


#/MetaFusion/test_data/trusight_cff/trusight.frame_info.cff
# trusight Frame Info
if [ $trusight_frame -eq 1 ]; then
echo trusight
outdir=$runs_dir/trusight.clinical.frame_info.$date
echo generating output in $outdir
mkdir $outdir
cff=/MetaFusion/test_data/trusight_cff/trusight.frame_info.cff

bash MetaFusion.clinical.sh --outdir $outdir \
                 --cff $cff  \
                 --gene_bed $gene_bed \
                 --annotate_exons \
                 --fusion_annotator \
                 --genome_fasta $genome_fasta \
                 --gene_info $gene_info \
                 --num_tools=2  \
                 --per_sample \
                 --recurrent_bedpe $recurrent_bedpe \
                 --scripts $fusiontools
fi

#BT474.KPL4.MCF7.SKBR3
if [ $brca_4 -eq 1 ]; then
echo BT474.KPL4.MCF7.SKBR3
#outdir=$runs_dir/BT474.KPL4.MCF7.SKBR3.$date
outdir=$runs_dir/BT474.KPL4.MCF7.SKBR3.7callers.chimerascan.$date
#outdir=$runs_dir/BT474.KPL4.MCF7.SKBR3.arriba.star_fusion.star_seqr.$date
#outdir=$runs_dir/BT474.KPL4.MCF7.SKBR3.arriba.star_fusion.star_seqr.Dec-10-2020
#outdir=$runs_dir/BT474.KPL4.MCF7.SKBR3.Dec-9-2020
cff=/MetaFusion/test_data/cff/BRCA.7callers.chimerascan.cff 
#cff=/MetaFusion/test_data/cff/BRCA.cff
#cff=/MetaFusion/test_data/cff_test/BRCA.arriba.star_fusion.star_seqr.cff
truth_fusions=/MetaFusion/test_data/truth_sets/BRCA.truth_set.dat

echo generating output in $outdir
bash MetaFusion.clinical.sh --outdir $outdir \
                 --cff $cff  \
                 --annotate_exons \
                 --gene_bed $gene_bed \
                 --fusion_annotator \
                 --gene_info $gene_info \
                 --genome_fasta $genome_fasta \
                 --truth_set $truth_fusions \
                 --num_tools=2 \
                 --per_sample \
                 --recurrent_bedpe $recurrent_bedpe \
                 --scripts $fusiontools
fi

# NTRK_control  
if [ $ntrk_control -eq 1 ]; then
echo NTRK_control 
outdir=$runs_dir/NTRK_control.$date
echo generating output in $outdir
truth_fusions=/MetaFusion/test_data/truth_sets/NTRK_control.truth_set.dat
cff=/MetaFusion/test_data/cff/NTRK_control.frame_info.cff

bash MetaFusion.clinical.sh --outdir $outdir \
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

fi

if [ $trusight_missed -eq 1 ]; then
echo TRUSIGHT missed 
outdir=$runs_dir/trusight_missed.$date
echo generating output in $outdir
mkdir $outdir
cff=/MetaFusion/test_data/trusight_cff/trusight.missed.cff

bash MetaFusion.clinical.sh --outdir $outdir \
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
                 --update_hist \
                 --ref_dir $ref_dir 
fi

if [ $trusight_2019_2020_2020b_2020c_dragen -eq 1 ]; then
echo TRUSIGHT 2019_2020_2020b_2020c_dragen
outdir=$runs_dir/trusight_2019_2020_2020b_2020c_dragen.$date
echo generating output in $outdir
mkdir $outdir
cff=/MetaFusion/test_data/trusight_cff/2019.2020.2020b.2020c_dragen.trusight.cff

bash MetaFusion.clinical.sh --outdir $outdir \
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
fi

if [ $trusight_2019_2020_2020b_2020c -eq 1 ]; then
echo TRUSIGHT 2019_2020_2020b_2020c 
outdir=$runs_dir/trusight_trusight_2019_2020_2020b_2020c.$date
echo generating output in $outdir
mkdir $outdir
cff=/MetaFusion/test_data/trusight_cff/2019.2020.2020b.2020c.trusight.cff

bash MetaFusion.clinical.sh --outdir $outdir \
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
                 --update_hist \
                 --ref_dir $ref_dir
fi

if [ $trusight_sample_subset_test -eq 1 ]; then
echo TRUSIGHT sample_subset_test
outdir=$runs_dir/trusight_sample_subset_test.$date
echo generating output in $outdir
mkdir $outdir
cff=/MetaFusion/test_data/trusight_cff/new_sample_subset_test.cff

bash MetaFusion.clinical.sh --outdir $outdir \
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
                 --update_hist \
                 --ref_dir $ref_dir
fi


