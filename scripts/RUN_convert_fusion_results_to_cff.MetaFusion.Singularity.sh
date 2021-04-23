#!/bin/bash

# to set
topdir=/hpf/largeprojects/ccmbio/mapostolides/PROJECTS/trusight/convert_dragen_to_cff
merged_cff_name=my_favourite_merged_cff_name.cff
sampleinfo=$topdir/sampleinfo.dragen.3LGG
outdir=$topdir/my_outdir_name
mkdir $outdir

# DON'T CHANGE BELOW THIS LINE (unless you know what you're doing :D   )
module load Singularity

metafusion_dir=/hpf/largeprojects/ccmbio/mapostolides/MetaFusion.clinical/

caller_file_dir=$topdir/caller_output_files/

dataset=trusight

tools=$(echo dragen )

for tool in ${tools[@]};do 

    raw_file_dir=$caller_file_dir/$dataset/$tool
    result_files=$(ls $raw_file_dir/*/*)
    
    echo generating cff for $tool
    for result_file in ${result_files[@]};do
    	sample=$(basename $(dirname $result_file))
    	echo $sample, $tool, $result_file, $outdir

      echo singularity exec -B $PWD -B $metafusion_dir $metafusion_dir/MetaFusion.clinical.simg python $metafusion_dir/scripts/convert_fusion_results_to_cff.py --sample $sample --sample_info_file $sampleinfo --tool $tool --fusion_result_file $result_file  --outdir $outdir
      singularity exec -B $PWD -B $metafusion_dir $metafusion_dir/MetaFusion.clinical.simg python $metafusion_dir/scripts/convert_fusion_results_to_cff.py --sample $sample --sample_info_file $sampleinfo --tool $tool --fusion_result_file $result_file  --outdir $outdir
    done

done

# Merge all .cff files into one combined file, "merged.cff"

files=$(ls $outdir/*cff | grep -v $merged_cff_name)
cat $files > $outdir/$merged_cff_name
