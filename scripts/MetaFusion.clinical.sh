#!/bin/bash
#STEPS
rename=1
annotate=1
annotate_exons=1
merge=1
isohunter=1
remove_non_hgnc=1
fusionannotator=1
db_filt_ann=1
RT_call_filter=1
blck_filter=1
ANC_filter=1
rank=1
add_junction_sequence=1
pfam_ann=1
final=1

# Default 0 variables
update_hist=0

# Loop through arguments and process them
while test $# -gt 0;do
    case $1 in
        -n=*|--num_tools=*)
        num_tools="${1#*=}"
        shift
        ;;
        --outdir)
        outdir="$2"
        shift 2
        ;;
        --cff)
        cff="$2"
        shift 2
        ;;
        --gene_bed)
        gene_bed="$2"
        shift 2
        ;;
        --gene_info)
        gene_info="$2"
        shift 2
        ;;
        --genome_fasta)
        genome_fasta="$2"
        shift 2
        ;;
        --ref_dir)
        reference_file_dir="$2"
        shift 2
        ;;
        --truth_set)
        truth_set="$2"
        shift 2
        ;;
        --recurrent_bedpe)
        recurrent_bedpe="$2"
        shift 2
        ;;
        --scripts)
        fusiontools="$2"
        shift 2
        ;;
        --fusion_annotator)
        FA=1
        shift
        ;;
        --annotate_exons)
        exons=1
        shift
        ;;
        --per_sample)
        per_sample=1
        shift
        ;;
        --update_hist)
        update_hist=1
        shift
        ;;
        --database)
        database="$2"
        shift 2
        ;;
        *)
        #OTHER_ARGUMENTS+=("$1")
        shift # Remove generic argument from processing
        ;;
    esac
done
#echo num_tools $num_tools 
#echo outdir $outdir 
#echo cff  $cff 
#echo gene_bed $gene_bed 
#echo gene_info $gene_info 
#echo genome_fasta $genome_fasta 
#echo truth_set $truth_set
echo UPDATE HIST: $update_hist
#METAFUSION WORKFLOW
mkdir $outdir

#Check CFF file format:
#Remove entries with nonconformming chromosome name
#Remove "." from strand field and replace with "NA"
cat $cff | awk '$1 ~ /[0-9XY]/ && $4 ~ /[0-9XY]/ ' |  awk 'BEGIN{FS=OFS="\t"} $3 !~ /^[-+]$/{$3="NA"} 1' | awk 'BEGIN{FS=OFS="\t"} $6 !~ /^[-+]$/{$6="NA"} 1'   > $outdir/$(basename $cff).reformat 
cff=$outdir/$(basename $cff).reformat
#NEED TO INSERT +/-/NA for strand, make NA if other

#Rename cff
if [ $rename -eq 1 ]; then
  echo Rename cff
  python rename_cff_file_genes.MetaFusion.py $cff $gene_info > $outdir/$(basename $cff).renamed
fi
cff=$outdir/$(basename $cff).renamed

#Annotate cff
if [ $annotate -eq 1 ]; then
  if [ $genome_fasta ]; then 
    echo Annotate cff, extract sequence surrounding breakpoint
    echo python reann_cff_fusion.py --cff $cff --gene_bed $gene_bed --ref_fa $genome_fasta \> $outdir/$(basename $cff).reann.WITH_SEQ
    python reann_cff_fusion.py --cff $cff --gene_bed $gene_bed --ref_fa $genome_fasta > $outdir/$(basename $cff).reann.WITH_SEQ
  else 
    echo Annotate cff, no extraction of sequence surrounding breakpoint
    echo python reann_cff_fusion.py --cff $cff --gene_bed $gene_bed \> $outdir/$(basename $cff).reann.NO_SEQ
    python reann_cff_fusion.py --cff $cff --gene_bed $gene_bed > $outdir/$(basename $cff).reann.NO_SEQ
  fi
fi
# Assign .cff based on SEQ or NOSEQ
if [ $genome_fasta ]; then 
  cff=$outdir/$(basename $cff).reann.WITH_SEQ
  echo cff $cff
else
  cff=$outdir/$(basename $cff).reann.NO_SEQ
  echo cff $cff
fi

if [ $annotate_exons -eq 1 ] && [ $exons -eq 1 ]; then
  echo Add adjacent exons to cff
  #echo python $fusiontools/extract_closest_exons.py $cff $gene_bed $genome_fasta  \> $outdir/$(basename $cff).exons
  #python $fusiontools/extract_closest_exons.py $cff $gene_bed $genome_fasta  > $outdir/$(basename $cff).exons
  python $fusiontools/extract_closest_exons.RNA_dev.py $cff $gene_bed $genome_fasta  > $outdir/$(basename $cff).exons
fi
# assign cff as ".exons" if --annotate_exons flag was specified
if [ $exons -eq 1 ]; then
  cff=$outdir/$(basename $cff).exons
fi

#Merge
cluster=$outdir/$(basename $cff).cluster
if [ $merge -eq 1 ]; then
  echo Merge cff by genes and breakpoints
  echo bash RUN_cluster_genes_breakpoints.sh $cff $outdir $fusiontools \> $cluster
  bash RUN_cluster_genes_breakpoints.sh $cff $outdir $fusiontools > $cluster
fi

#ISOHUNTER
if [ $isohunter -eq 1 ]; then
  echo Merge cff by breakpoints to separate isoforms 
    echo bash $fusiontools/RUN_cluster_breakpoints.sh $cff $outdir $fusiontools $per_sample \> $outdir/$(basename $cff).isoforms.per_sample.cluster 
    bash $fusiontools/RUN_cluster_breakpoints.sh $cff $outdir $fusiontools $per_sample > $outdir/$(basename $cff).isoforms.per_sample.cluster 
fi
cluster=$outdir/$(basename $cff).isoforms.per_sample.cluster

# REMOVE NON-HGNC
if [ $remove_non_hgnc -eq 1 ]; then
  echo Rscript $fusiontools/remove_non_hgnc_metafusion.R --cluster=$cluster --out=$outdir/$(basename $cluster).non_hgnc_removed 
  Rscript $fusiontools/remove_non_hgnc_metafusion.R --cluster=$cluster --out=$outdir/$(basename $cluster).non_hgnc_removed 
fi
cluster=$outdir/$(basename $cluster).non_hgnc_removed

#FUSIONANNOTATOR
reference_files=$(dirname $gene_bed)
if [ $fusionannotator -eq 1 ] && [ $FA -eq 1 ]; then
  echo Running FusionAnnotator 
  bash RUN_FusionAnnotator.sh $outdir $cluster $fusiontools $reference_files full
  echo Adding FusionAnnotator database hits to $cluster file
  FA_db_file=$outdir/cluster.preds.collected.gencode_mapped.wAnnot.CANCER_FUSIONS
  echo python add_db_hits_to_cluster.py $cluster $FA_db_file \> $outdir/$(basename $cluster).cancer_db_hits
  python add_db_hits_to_cluster.py $cluster $FA_db_file full > $outdir/$(basename $cluster).cancer_db_hits
fi
cluster=$outdir/$(basename $cluster).cancer_db_hits


# FILTER and annotate using historical database
# 0) check to see database file exists
# 1) filter out manually curated false positives
# 2) annotate with historical calls and clinical table 
# 3) Update historical_table with calls from current run
# 4) SEPARATE CLINICAL calls so they bypass filters
ls $database
if [ $? -eq 1 ]; then
    echo database file $database does not exist or was not submitted as an argument 
    exit 1
fi

if [ $db_filt_ann -eq 1 ]; then
  if [ $update_hist -eq 1 ]; then
	echo Updating historical fusion table with current run
  else
	echo NOT updating historical fusions table with current run
  fi
  echo Rscript $fusiontools/database_filter_ann_update_metafusion.R --cluster=$cluster --out=$outdir/$(basename $cluster).db_filt_ann --update_hist=$update_hist --database=$database
  Rscript $fusiontools/database_filter_ann_update_metafusion.R --cluster=$cluster --out=$outdir/$(basename $cluster).db_filt_ann --update_hist=$update_hist --database=$database
  # SEPERATE CLINICAL calls so they bypass filters
  cat $outdir/$(basename $cluster).db_filt_ann | awk '$30 > 0' > $outdir/$(basename $cluster).db_filt_ann.clinical_hits
fi
cluster_clinical=$outdir/$(basename $cluster).db_filt_ann.clinical_hits
cluster=$outdir/$(basename $cluster).db_filt_ann

#ReadThrough Callerfilter
if [ $RT_call_filter -eq 1 ]; then
  echo ReadThrough, callerfilter $num_tools
  cat $cluster | grep ReadThrough > $outdir/$(basename $cluster).ReadThrough
  cat $cluster | grep -v ReadThrough > $outdir/$(basename $cluster).RT_filter
  python callerfilter_num.py --cluster $outdir/$(basename $cluster).RT_filter --num_tools $num_tools > $outdir/$(basename $cluster).RT_filter.callerfilter.$num_tools
fi
if [ $? -eq 1 ]; then
	echo No calls pass ReadThrough and callerfilter. See the following file for unfiltered calls
	echo $cluster
	exit 1
fi
cluster_RT_call=$outdir/$(basename $cluster).RT_filter.callerfilter.$num_tools 


# Blocklist Filter
if [ $blck_filter -eq 1 ]; then
  echo blocklist filter
  #echo bash blocklist_filter_recurrent_breakpoints.sh $cff $cluster_RT_call $outdir $recurrent_bedpe
  bash blocklist_filter_recurrent_breakpoints.sh $cff $cluster_RT_call $outdir $recurrent_bedpe > $outdir/$(basename $cluster).RT_filter.callerfilter.$num_tools.blck_filter
fi
cluster=$outdir/$(basename $cluster).RT_filter.callerfilter.$num_tools.blck_filter

# Adjacent Noncoding filter 
if [ $ANC_filter -eq 1 ]; then
  echo ANC adjacent noncoding filter
  python filter_adjacent_noncoding.py $cluster > $outdir/$(basename $cluster).ANC_filter  
fi
cluster=$outdir/$(basename $cluster).ANC_filter

#Rank and add in clinical calls 
if [ $rank -eq 1 ]; then
  echo Rank  
  # add clinical calls back in which bypassed filters
  cat $cluster_clinical | head -1 > $outdir/$(basename $cluster).with_clinical
  cat $cluster_clinical $cluster | grep -v cluster_type | sort | uniq >> $outdir/$(basename $cluster).with_clinical
  # rank and remove duplicates
  python rank_cluster_file.py $outdir/$(basename $cluster).with_clinical | uniq > $outdir/$(basename $cluster).with_clinical.ranked 
fi
cluster=$outdir/$(basename $cluster).with_clinical.ranked

# ADD JUNCTION SEQUENCE TO FINAL OUTPUT FILE
if [ $add_junction_sequence -eq 1 ]; then
    echo add junction sequence to cluster file 
    echo python add_junction_seq.py --cluster $cluster --ref_fa $genome_fasta \> $outdir/$(basename $cluster).jxn_seq
    python add_junction_seq.py --cluster $cluster --ref_fa $genome_fasta > $outdir/$(basename $cluster).jxn_seq
fi
cluster=$outdir/$(basename $cluster).jxn_seq

if [ $pfam_ann -eq 1 ]; then
   gene_symbol_file=$reference_file_dir/hgTables.gene_symbol.ENSG.ENST.ENSP.Nov-13-2020.tsv
   pfam_domains_per_gene=$reference_file_dir/domains_per_gene.Rds
   pfam_names_descriptions=$reference_file_dir/pdb_pfam_mapping.txt
   echo Rscript $fusiontools/annotate_metafusion_output_with_pfam.R --cluster=$cluster --out=$outdir/$(basename $cluster).pfam_ann --gene_symbols=$gene_symbol_file --pfam_gene_ann=$pfam_domains_per_gene --pfam_names_descriptions=$pfam_names_descriptions
   Rscript $fusiontools/annotate_metafusion_output_with_pfam.R --cluster=$cluster --out=$outdir/$(basename $cluster).pfam_ann --gene_symbols=$gene_symbol_file --pfam_gene_ann=$pfam_domains_per_gene --pfam_names_descriptions=$pfam_names_descriptions
else
     echo PFAM annotation not performed
fi
cluster=$outdir/$(basename $cluster).pfam_ann

if [ $final -eq 1 ]; then
  #Generate final cluster
  final_cluster=$outdir/final.n$num_tools.cluster
  cat $cluster > $final_cluster
  # Generate final xlsx
  echo Rscript write_final.R --cluster=$final_cluster
  Rscript write_final.R --cluster=$final_cluster
fi
