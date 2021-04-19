#!/usr/bin/env python
import sys
import pandas as pd
import numpy as np
import  pygeneann_MetaFusion as pygeneann
import pybedtools.bedtool as bedtools
import itertools
import sequtils
import argparse
from collections import Counter

parser = argparse.ArgumentParser()
parser.add_argument('cff', action='store', help='CFF file')
parser.add_argument('FIDs', action='store', help='clustered FIDs file')
parser.add_argument('--choose_bp', action='store_true')

args = parser.parse_args()
cff=args.cff
FIDs=args.FIDs
choose_bp=args.choose_bp

def choose_most_frequent_breakpoint(breakpoint_list):
    """ Selects breakpoint which is most common from breakpoint_list
    returns: list of breakpoints
    """ 
    #counts = Counter(breakpoint_list) 
    most_common_breakpoint = Counter(breakpoint_list).most_common(1)[0][0]
    return [str(most_common_breakpoint)]

def choose_most_frequent_chromosome(chr_list):
    """ Selects chr which is most common from chr_list
    returns: list of chr 
    """ 
    most_common_chr = Counter(chr_list).most_common(1)[0][0]
    return [str(most_common_chr)]

def output_clustered_fusions(fusion_list, cluster_type):
        use_reference_fusion=1
        #use first fusion in fusion as reference fusion to account for flipped fusions for BP_Cluster
        if use_reference_fusion:
            try:
                # exclude defuse as ref fusion since defuse often gets the order wrong
                ref_fus = [fusion for fusion in fusion_list if fusion.tool != 'defuse'][0]
            except IndexError:
                ref_fus = fusion_list[0]
            #flip all fusions not ordered the same way as reference:
            #print("REFERENCE FUSION:")
            #print(ref_fus.tostring())
            for fusion in fusion_list:
                #print(fusion.tool,ref_fus.t_gene1,fusion.t_gene1,ref_fus.t_gene2,fusion.t_gene2,  not (ref_fus.t_gene1 == fusion.t_gene1) and ( ref_fus.t_gene1 == fusion.t_gene2 or ref_fus.t_gene2 == fusion.t_gene1 ))
                if not (ref_fus.t_gene1 == fusion.t_gene1) and ( ref_fus.t_gene1 == fusion.t_gene2 or ref_fus.t_gene2 == fusion.t_gene1 ):
                    #flip fusion:
                    fusion.chr1, fusion.chr2 = fusion.chr2, fusion.chr1
                    fusion.pos1, fusion.pos2 = fusion.pos2, fusion.pos1
                    fusion.strand1, fusion.strand2 = fusion.strand2, fusion.strand1
                    fusion.t_gene1, fusion.t_gene2 = fusion.t_gene2, fusion.t_gene1
                    fusion.t_area1, fusion.t_area2 = fusion.t_area2, fusion.t_area1
                    fusion.reann_type1,fusion.reann_type2 = fusion.reann_type2,fusion.reann_type1
                    fusion.closest_exon1,fusion.closest_exon2 = fusion.closest_exon2,fusion.closest_exon1
                    fusion.rna_type1, fusion.rna_type2 = fusion.rna_type2, fusion.rna_type1 
        #need to remove NA values from below lists, but only if there are valid gene names in them
        gene1_list_orig = [f.t_gene1 for f in fusion_list]
        gene1_list = [gene for gene in gene1_list_orig if gene !="NA"]
        if len(gene1_list) < 1:
            gene1_list = gene1_list_orig
        gene2_list_orig = [f.t_gene2 for f in fusion_list]
        gene2_list = [gene for gene in gene2_list_orig if gene !="NA"]
        if len(gene2_list) < 1:
            gene2_list = gene2_list_orig

        max_split_cnt = max([f.split_cnt for f in fusion_list])
        max_span_cnt = max([f.span_cnt for f in fusion_list])
        sample_list = [f.sample_name for f in fusion_list]
        
        fusion_IDs = [f.fusion_id for f in fusion_list]
        disease_list = [f.disease for f in fusion_list]
        tool_list = [f.tool for f in fusion_list]
        sample_type_list = [f.sample_type for f in fusion_list]
        gene1_on_bndry = "True" in [f.gene1_on_bndry for f in fusion_list]
        gene1_close_to_bndry = "True" in [f.gene1_close_to_bndry for f in fusion_list]
        gene2_on_bndry = "True" in [f.gene2_on_bndry for f in fusion_list]
        gene2_close_to_bndry = "True" in [f.gene2_close_to_bndry for f in fusion_list]
        #closest exons
        exon1 = ",".join(set([f.closest_exon1 for f in fusion_list]))
        exon2 = ",".join(set([f.closest_exon2 for f in fusion_list]))
        #categories
        category_list = [f.category for f in fusion_list]

        # PRIORITIZE CATEGORIES TO REMOVE MULTIPLE CATEGORIES, also rename categories
        if "ReadThrough" in category_list:
            category_list = ["ReadThrough"]
        elif "GeneFusion" in category_list:
            category_list = ["CodingFusion"]
        elif "TruncatedCoding" in category_list:
            category_list = ["TruncatedCoding"]
        elif "TruncatedNoncoding" in category_list:
            category_list = ["TruncatedNoncoding"]
        elif "NoDriverGene" in category_list:
            category_list = ["NoHeadGene"]

        # added 4 new fields to final .cluster file for purposes of validation 
        chr1_list = choose_most_frequent_chromosome([str(f.chr1) for f in fusion_list])
        breakpoint_1_list = [str(f.pos1) for f in fusion_list]
        if choose_bp:
            breakpoint_1_list = choose_most_frequent_breakpoint(breakpoint_1_list)
        chr2_list = choose_most_frequent_chromosome([str(f.chr2) for f in fusion_list])
        breakpoint_2_list = [str(f.pos2) for f in fusion_list]
        if choose_bp:
            breakpoint_2_list = choose_most_frequent_breakpoint(breakpoint_2_list)
        #cancer DB hits placeholder:
        cancer_db_hits = "NA"
        #RNA types
        rna_type1 = list(set([f.rna_type1 for f in fusion_list]))
        rna_type2 = list(set([f.rna_type2 for f in fusion_list]))
        # Frame info
        frame = []
        for fusion in fusion_list:
            if fusion.tool == "fusionmap": frame.append("fusionmap:" + fusion.frame)
            elif fusion.tool == "arriba": frame.append("arriba:" + fusion.frame)
            elif fusion.tool == "star_fusion": frame.append("star_fusion:" + fusion.frame)
            else: continue
        if len(frame) == 0: frame = ["NA"]
        frame.sort()
        # Strand info
        strand1 = set(fusion.strand1 for fusion in fusion_list if fusion.tool != "defuse")
        strand2 = set(fusion.strand2 for fusion in fusion_list if fusion.tool != "defuse")
        if len(strand1) == 0: strand1=set(["NA"])
        if len(strand2) == 0: strand2=set(["NA"])
        #print strand info
        #print([fusion.strand1 for fusion in fusion_list])
        #print([fusion.strand2 for fusion in fusion_list])
        
        # Main output: correct version
        print "\t".join(map(str, [cluster_type, ",".join(list(set(gene1_list))), ",".join(list(set(gene2_list))), max_split_cnt, max_span_cnt, ",".join(list(set(sample_type_list))), ",".join(list(set(disease_list))), ",".join(list(set(tool_list))), ",".join(list(set(category_list))), gene1_on_bndry, gene1_close_to_bndry, gene2_on_bndry, gene2_close_to_bndry, cancer_db_hits, ",".join(list(set(sample_list))), ",".join(list(set(chr1_list))), "|".join(list(set(breakpoint_1_list))), ",".join(list(set(chr2_list))), "|".join(list(set(breakpoint_2_list))), exon1, exon2 ,",".join(list(set(fusion_IDs))), ",".join(rna_type1), ",".join(rna_type2), len(list(set(tool_list))), ",".join(list(set(frame))), ",".join(list(strand1)), ",".join(list(strand2)) ]))
        # below removes "set" from breakpoints and callers, allowing us to see which breakpoints correspond to which callers
        #print "\t".join(map(str, [cluster_type, ",".join(list(set(gene1_list))), ",".join(list(set(gene2_list))), max_split_cnt, max_span_cnt, ",".join(list(set(sample_type_list))), ",".join(list(set(disease_list))), ",".join(list(tool_list)), ",".join(list(set(category_list))), gene1_on_bndry, gene1_close_to_bndry, gene2_on_bndry, gene2_close_to_bndry, cancer_db_hits, ",".join(list(set(sample_list))), ",".join(list(set(chr1_list))), "|".join(list(breakpoint_1_list)), ",".join(list(set(chr2_list))), "|".join(list(breakpoint_2_list)), exon1, exon2 ,",".join(list(set(fusion_IDs))), ",".join(rna_type1), ",".join(rna_type2), len(list(set(tool_list))), ",".join(list(set(frame))) ]))


# Load cff file
lines=[line for line in open(cff, "r")]
fusion=pygeneann.CffFusion(lines[0])
header=fusion.zone1_attrs + fusion.zone2_attrs + fusion.zone3_attrs + fusion.zone4_attrs
df_cff=pd.read_csv(cff, sep='\t', keep_default_na=False, index_col=False, names=header)

# load FIDs file
FID_clusters = [line for line in open(FIDs, "r")]
# output header
pygeneann.output_cluster_header()
for cluster in FID_clusters:
    if cluster.startswith('FIDs'): continue
    FID_lst = cluster.rstrip().split(",")
    df_cluster = df_cff[df_cff['fusion_id'].isin(FID_lst)].to_csv(header=None, index=False, sep="\t").rstrip()
    fusion_list = []
    for line in df_cluster.split("\n"):
        fusion=pygeneann.CffFusion(line)
        fusion_list.append(fusion)
    #output fusion list
    output_clustered_fusions(fusion_list, "Cluster")
