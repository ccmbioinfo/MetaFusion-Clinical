#!/usr/bin/env python
import sys
import  pygeneann_MetaFusion as pygeneann
import sequtils
import pysam
import argparse

# debug=1 for head, debug = 2 for tail
# set debug = 3 for STX16--RAE1
debug=0

if debug == 3:
  cff_file = "/MetaFusion/EXON_RUNS/outdir-Aug-27-2020/STX16--RAE1.ALL.cff.reann.NO_SEQ.55929088.arriba.1entry" 
  ensbed="/MetaFusion/test_data/cff_test/ens_known_genes.STX16--RAE1.renamed.bed"
  ref_fa="/MetaFusion/reference_files/human_g1k_v37_decoy.fasta"
#PARSER
else:
  parser = argparse.ArgumentParser()
  parser.add_argument('cff_file', action='store', help='CFF file, can be .cff or cff.reann')
  parser.add_argument('ensbed', action='store', help='Ensemble gene file')
  parser.add_argument('ref_fa', action='store', help='Reference genome file')
  args = parser.parse_args()
  cff_file = args.cff_file
  ensbed = args.ensbed
  ref_fa=args.ref_fa

def choose_closest_exon(exon_list, pos, strand, ori):
    #chooses the closest exon from a list of adjacent exons either upstream (head) or downstream (tail) of the fusion breakpoint
    smallest_dist=float('inf')
    slop=10 # keep "smallest_dist" the same if within slop 
    closest_exon=""
    if debug == 2 and ori == 'tail': 
        sys.stderr.write("There are " + str(len(exon_list)) + " exons in the list\n")
        sys.stderr.write(str([",".join([str(exon.chr) + ":" + str(exon.start) + "-" + str(exon.end), 
                                         exon.gene_name, exon.type]) for exon in exon_list]) + "\n")
    for exon in exon_list:
      if debug == 1 and ori=="head": sys.stderr.write(str(exon.__dict__) + "\n")
      if debug == 2 and ori=="tail": sys.stderr.write(str(exon.__dict__) + "\n")
      exon_igv=str(exon.chr) + ":" + str(exon.start) + "-" + str(exon.end)
      # get dist, exon.start will always be larger value than pos2
      if (strand=="+" and ori=="tail") or (strand=="-" and ori=="head"):
        dist=abs(exon.start - pos) 
        if debug == 1 and ori=="head": sys.stderr.write("For exon " + exon_igv + ", the value of exon.start - pos is " + str(dist) + "\n")
        if debug == 2 and ori=="tail": sys.stderr.write("For exon " + exon_igv + ", the value of exon.start - pos is " + str(dist) + "\n")
      elif (strand=="+" and ori=="head") or (strand=="-" and ori=="tail"):
        dist=abs(pos - exon.end) 
        if debug == 1 and ori=="head": sys.stderr.write("For exon " + exon_igv + ", the value of pos - exon.start is " + str(dist) + "\n")
        if debug == 2 and ori=="tail": sys.stderr.write("For exon " + exon_igv + ", the value of pos - exon.start is " + str(dist) + "\n")
      # overwrite closest_exon_lst, must be at least "slop" bp closer than "smallest_dist"
      if dist < smallest_dist - slop: 
          closest_exon_lst = [exon]
          smallest_dist = dist
      # append exon to closest_exon_lst if within "slop" bp
      elif smallest_dist - slop < dist < smallest_dist + slop:
          closest_exon_lst.append(exon)
    #filter exons by RNA type
    primary_rna_types=["lincRNA", "miRNA","mRNA", "snRNA",  "snoRNA","rRNA"]
    filtered_exon_lst=[exon for exon in closest_exon_lst if exon.rna_type in primary_rna_types]
    if len(filtered_exon_lst) > 0: closest_exon_lst=filtered_exon_lst

    # choose longest exon from closest_exon_lst
    max_exon_length = max([exon.length for exon in closest_exon_lst])
    for exon in sorted(closest_exon_lst): 
      if exon.length == max_exon_length: closest_exon = exon
    if debug and ori == 'tail': 
      sys.stderr.write("Closest exon: " + str(closest_exon.chr) + ":" + str(closest_exon.start) + "-" + str(closest_exon.end) + "\n")
      #if (strand=="+" and ori=="tail") or (strand=="-" and ori=="head"): print(closest_exon.start - pos)
      #elif (strand=="+" and ori=="head") or (strand=="-" and ori=="tail"): print(pos - exon.start)
    return closest_exon

def assign_exon_to_cff(adjacent_exons, fusion, ori):
    """
    Returns exon_str, fusion.rna_type 
    """
    if ori == "head": pos,strand = fusion.pos1, fusion.strand1
    elif ori == "tail":  pos,strand =fusion.pos2, fusion.strand2
    try:
      exon=choose_closest_exon(adjacent_exons, pos, strand, ori)
      exon_str=str(exon.chr) + ":" + str(exon.start) + "-" + str(exon.end)
    except:
      exon_str="NA"
      exon="NA"
    if exon == "NA":
        return  exon_str,"NA" 
    else:
        return  exon_str,exon.rna_type

gene_ann = pygeneann.GeneAnnotation(ensbed)
for line in open(cff_file, "r"):
    fusion = pygeneann.CffFusion(line)
    adjacent_exons1 = gene_ann.get_closest_exon_lists(fusion.chr1, fusion.pos1)
    adjacent_exons2 = gene_ann.get_closest_exon_lists(fusion.chr2, fusion.pos2)
    #exons downstream of breakpoint on + strand
    if debug == 1:
      sys.stderr.write("HEAD GENE: " + fusion.t_gene1 + "\n")
      sys.stderr.write("\t".join(["BREAKPOINT: " , str(fusion.chr1) , str(fusion.pos1) , str(fusion.strand1), "\n"])) 
      sys.stderr.write("PREV HEAD features\n")
      sys.stderr.write(str([str(exon.chr) + ":" + str(exon.start) + "-" + str(exon.end) + "|" + str(exon.rna_type) for exon in adjacent_exons1[0]]) + "\n")
      sys.stderr.write("NEXT HEAD features\n")
      sys.stderr.write(str([str(exon.chr) + ":" + str(exon.start) + "-" + str(exon.end) + "|" + str(exon.rna_type) for exon in adjacent_exons1[1]])+ "\n" )
    #HEAD LEFT OF BREAKPOINT (PREV) 
    if fusion.strand1 == '+':
      fusion.closest_exon1, fusion.rna_type1 = assign_exon_to_cff(adjacent_exons1[0], fusion, "head")
    #HEAD RIGHT OF BREAKPOINT (NEXT) 
    elif fusion.strand1 == '-':
      fusion.closest_exon1, fusion.rna_type1 = assign_exon_to_cff(adjacent_exons1[1], fusion, "head")
    else: raise Exception('Strand must be either + or -')

    if debug == 2: 
      sys.stderr.write("TAIL GENE: " + fusion.t_gene2 + "\n")
      sys.stderr.write("\t".join(["BREAKPOINT: " , str(fusion.chr2) , str(fusion.pos2) , str(fusion.strand2), "\n"])) 
      sys.stderr.write("PREV TAIL features\n")
      sys.stderr.write(str([str(exon.chr) + ":" + str(exon.start) + "-" + str(exon.end) + "|" + str(exon.rna_type) for exon in adjacent_exons2[0]]) + "\n")
      sys.stderr.write("NEXT TAIL features\n")
      sys.stderr.write(str([str(exon.chr) + ":" + str(exon.start) + "-" + str(exon.end) + "|" + str(exon.rna_type) for exon in adjacent_exons2[1]])+ "\n" )

    #TAIL LEFT OF BREAKPOINT (PREV)  
    if fusion.strand2 == '-':
      if debug == 2: sys.stderr.write("Choosing tail gene feature on '" + fusion.strand2 + "' strand\n")
      fusion.closest_exon2, fusion.rna_type2 = assign_exon_to_cff(adjacent_exons2[0], fusion, "tail")
    #TAIL RIGHT OF BREAKPOINT (NEXT)
    elif fusion.strand2 == '+':
      fusion.closest_exon2, fusion.rna_type2 = assign_exon_to_cff(adjacent_exons2[1], fusion, "tail")
    else: raise Exception('Strand must be either + or -')
    print fusion.tostring() 
