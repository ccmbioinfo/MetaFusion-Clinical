#!/bin/bash
#module load R/3.5.1
#module load python/2.7.11
cff=$1
outdir=$2
fusiontools_dir=$3
per_sample=$4

#Generate intersections file for breakpoints (df_bed.pair_to_pair(df_bed, slop=10, rdn=False))
if [ $per_sample -eq 1 ]; then
  fid_intersection_file=$outdir/FID.intersections.breakpoints.per_sample.tsv
  python $fusiontools_dir/intersect_breakpoints.per_sample.py $cff > $fid_intersection_file
else
  fid_intersection_file=$outdir/FID.intersections.breakpoints.tsv
  python $fusiontools_dir/intersect_breakpoints.py $cff > $fid_intersection_file 
fi

# Graph clustering
fid_clusters_file=$outdir/FID.clusters.tsv
#ls -l $fid_intersection_file
#ls -l $fusiontools_dir/cluster_intersections.local.R
Rscript $fusiontools_dir/cluster_intersections.R $fid_intersection_file $fid_clusters_file 

# generate cluster file using clustered FIDs and cff file
#echo python $fusiontools_dir/generate_cluster_file.py $cff $fid_clusters_file
python $fusiontools_dir/generate_cluster_file.py --choose_bp $cff $fid_clusters_file
