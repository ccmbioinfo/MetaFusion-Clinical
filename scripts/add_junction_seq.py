#!/usr/bin/env python
import sys
import  pygeneann_MetaFusion as pygeneann
import sequtils
import pysam
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('--cluster', action='store', help='cluster file')
parser.add_argument('--ref_fa', required=False, action='store', help='Reference genome file')

args = parser.parse_args()
cluster_file = args.cluster
# Assign reference fasta if provided by user
if args.ref_fa is not None:
  ref_fa=args.ref_fa

pygeneann.output_cluster_header()
for line in open(cluster_file, "r"):
    if line.startswith("#"):
        continue
    fusion = pygeneann.CategoryFusions(line)

    # Get fusion seq only if specified by user
    if args.ref_fa is not None: 
      #pygeneann.get_fusion_seq_cluster(fusion, ref_fa, 30)
      pygeneann.get_fusion_seq_cluster(fusion, ref_fa, 100)

    fusion.out()

