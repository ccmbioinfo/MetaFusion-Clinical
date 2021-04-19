#!/usr/bin/env python
import sys
import pandas as pd
import numpy as np
import pygeneann_MetaFusion as pygeneann
import pybedtools.bedtool as bedtools
import itertools
import sequtils
import argparse

debug = 1

parser = argparse.ArgumentParser()
parser.add_argument('cff', action='store', help='CFF file')

args = parser.parse_args()
cff=args.cff
#cff="/MetaFusion/EXON_RUNS/trusight.exon_annot.Oct-20-2020/KIAA--BRAF.exons.cff"

#INTERSECT FUSIONS BY BREAKPOINTS
def intersect_fusions_by_breakpoints():
    lines=[line for line in open(cff, "r")]
    fusion=pygeneann.CffFusion(lines[0])
    header=fusion.zone1_attrs + fusion.zone2_attrs + fusion.zone3_attrs + fusion.zone4_attrs
    df_cff=pd.read_csv(cff, sep='\t', keep_default_na=False, index_col=False, names=header)

    # Combine sample name and chr to allow for only same-sample intersections
    df_cff['chr1'] = df_cff['chr1'] + "_" + df_cff['sample_name']
    df_cff['chr2'] = df_cff['chr2'] + "_" + df_cff['sample_name']
    #print(df_cff['chr1'])
    #exit(0)
    #create BedTools object with appropriate column names
    print >> sys.stderr, "create BedTools object with appropriate column names"
    df_bed=df_cff[['chr1','pos1','pos1','chr2','pos2','pos2', 'fusion_id']]
    df_bed.columns=['chr1','pos1','pos1_2','chr2','pos2','pos2_2', 'fusion_id']
    df_bed.loc[:,['pos1_2','pos2_2']] +=1
    df_bed=bedtools.BedTool.from_dataframe(df_bed)

    #Intersect fusions: NOTE: only keeps fusions that intersect
    print >> sys.stderr, "Intersect fusions: NOTE: rdn=False, keeps self-intersections"
    df_intersect=df_bed.pair_to_pair(df_bed, slop=10, rdn=False)
    df=df_intersect.to_dataframe(header=None).iloc[:,0:14]
    df.columns = ['chr1','pos1','pos1_2','chr2','pos2','pos2_2', 'fusion_id', 'chr1_1','pos1_1','pos1_2_1','chr2_1','pos2_1','pos2_2_1', 'fusion_id_lst']
    df=df[['fusion_id','fusion_id_lst']]
    #write paired F_IDs to tsv
    return df

df = intersect_fusions_by_breakpoints()
df.to_csv(sys.stdout,header=True,index=True, sep="\t")

exit(0)
#INTERSECT FUSIONS BY EXONS 
def intersect_exons_bedtools():
    lines=[line for line in open(cff, "r")]
    fusion=pygeneann.CffFusion(lines[0])
    header=fusion.zone1_attrs + fusion.zone2_attrs + fusion.zone3_attrs + fusion.zone4_attrs
    sys.stderr.write(str(header) + "\n")
    df_cff=pd.read_csv(cff, sep='\t', keep_default_na=False, index_col=False, names=header)
    
    # Split exons into 'chr1','start1','end1','chr2','start2','end2'
    # Split closest_exon1
    df_cff[['exon_chr1','exon_range1']] = df_cff['closest_exon1'].str.split(':',expand=True)
    df_cff[['exon_start1','exon_end1']] = df_cff['exon_range1'].str.split('-',expand=True)

    # Split closest_exon2
    df_cff[['exon_chr2','exon_range2']] = df_cff['closest_exon2'].str.split(':',expand=True)
    df_cff[['exon_start2','exon_end2']] = df_cff['exon_range2'].str.split('-',expand=True)

    
    #create BedTools object with appropriate column names
    print >> sys.stderr, "create BedTools object with appropriate column names"
    df_bed=df_cff[['exon_chr1','exon_start1','exon_end1','exon_chr2','exon_start2','exon_end2', 'fusion_id']]
    #df_bed.loc[:,['pos1_2','pos2_2']] +=1
    df_bed=bedtools.BedTool.from_dataframe(df_bed)
    
    #Intersect fusions: 
    print >> sys.stderr, "Intersect fusions: NOTE: rdn=False, keeps self-intersections"
    #df_intersect=df_bed.pair_to_pair(df_bed, slop=100, rdn=False)
    df_intersect=df_bed.pair_to_pair(df_bed, rdn=False)
    df=df_intersect.to_dataframe(header=None).iloc[:,0:14]
    df.columns = ['chr1','pos1','pos1_2','chr2','pos2','pos2_2', 'fusion_id', 'chr1_1','pos1_1','pos1_2_1','chr2_1','pos2_1','pos2_2_1', 'fusion_id_lst'] 
    df=df[['fusion_id','fusion_id_lst']]
    #write paired F_IDs to tsv
    return df

df = intersect_exons_bedtools()
df.to_csv(sys.stdout,header=True,index=True, sep="\t")

exit(0)
fusion_dict = intersect_fusions_by_exons(cff)
if debug: sys.stderr.write(str(fusion_dict) + "\n")
#exit(0)
#count = df.shape[0] + 1 
count = 0
#print header
print("\t" + "fusion_id" + "\t" + "fusion_id_lst")
for key in fusion_dict.keys():
    lst=fusion_dict[key]
    #print(lst)
    edges=list(itertools.permutations(lst, 2))
    #print(edges)
    for edge in edges:
        print("\t".join([str(count)] + list(edge)))
        count += 1
