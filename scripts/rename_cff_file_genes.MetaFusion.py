#!/usr/bin/env python
import os
import sys
import pygeneann_MetaFusion as pygeneann
import pandas as pd
import sequtils
import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument('cff_file', action='store', help='CFF file before annotation. if there are multiple gene names in a field, names MUST be comma-separated lists')
parser.add_argument('gene_info_file', action='store', help='Homo_sapiens.gene_info')

args = parser.parse_args()
cff_file = args.cff_file
gene_info_file=args.gene_info_file

# tab-delimited file with a column contating all aliases for each gene
#tax_id	GeneID	Symbol	LocusTag	Synonyms	dbXrefs	chromosome	map_location	description	type_of_gene	Symbol_from_nomenclature_authority	Full_name_from_nomenclature_authority	Nomenclature_status	Other_designations	Modification_date	Feature_type
#9606	1	A1BG	-	A1B|ABG|GAB|HYST2477	MIM:138670|HGNC:HGNC:5|Ensembl:ENSG00000121410	19	19q13.43	alpha-1-B glycoprotein	protein-coding	A1BG	alpha-1-B glycoprotein	O	alpha-1B-glycoprotein|HEL-S-163pA|epididymis secretory sperm binding protein Li 163pA	20200313	-

#Open NCBI file and create df
ncbi_gene_info_file = open(gene_info_file)
df = pd.read_csv(ncbi_gene_info_file, sep='\t')

#pd.set_option('display.max_rows', None)
#pd.set_option('display.max_columns', None)

#Create one row for each alias
#https://stackoverflow.com/questions/58523316/split-rows-in-pandas-dataframe
df = pd.concat( (df.Synonyms.str.split('|', expand=True), df[['Symbol','chromosome']]), axis=1).melt(id_vars=['Symbol','chromosome'], value_name='Synomyms').dropna()

df.rename(columns = {'Synomyms': 'Alias', 'Symbol': 'HGNC_Symbol'}, inplace = True)
#Setting column values all to "str" so they can be compared with Python strings
df = df.astype(str)


def clean_weird_input(gene_name): 
    #$ cat $cluster | sed 's/(.\+)//g' | sed 's/\//,/g' 
    # remove arriba brackets, e.g. RNU6-121P(14818),AQP10(8524)
    gene_name = re.sub('\([0-9]+\)', '', gene_name) 
    # remove integrate slash-delimited paralogs
    gene_name = re.sub('/', ',', gene_name) 
    return gene_name 

def alias2hgnc(df, query, chr):
    #check to see if query is already HGNC symbol
    hgnc = df.loc[df.HGNC_Symbol == query].HGNC_Symbol.values.tolist()
    is_hgnc = True if len(hgnc) > 0 else False 
    if is_hgnc: 
        return query
    else:
        hgnc_lst = df.loc[(df['Alias'] == query) & (df['chromosome'].str.match(chr))].HGNC_Symbol.values.tolist()
        #sys.stderr.write(str(hgnc_lst) + "\n")
        return "NA" if len(hgnc_lst) == 0 else hgnc_lst[0]

def select_gene_name(gene_lst, chr):
    """takes a list of gene names and selects a single name"""
    #handle list of genes from caller
    gene_hgnc_lst = []
    for gene in gene_lst:
        hgnc_gene = alias2hgnc(df, gene, chr)
        gene_hgnc_lst.append(hgnc_gene)
    try: return [gen for gen in gene_hgnc_lst if gen != "NA"][0]
    except: return ','.join([gen for gen in gene_hgnc_lst if gen != "NA"])

for line in open(cff_file, "r"):
    fusion = pygeneann.CffFusion(line)
    #clean weird input
    fusion.t_gene1 = clean_weird_input(fusion.t_gene1) 
    fusion.t_gene2 = clean_weird_input(fusion.t_gene2) 
    # t_gene1
    t_gene1 = select_gene_name(fusion.t_gene1.split(","), fusion.chr1[3:5])# remove "chr" prefix
    if len(t_gene1) != 0: fusion.t_gene1 = t_gene1 
    # t_gene2
    t_gene2 = select_gene_name(fusion.t_gene2.split(","), fusion.chr2[3:5])
    if len(t_gene2) != 0: fusion.t_gene2 = t_gene2 

    print(fusion.tostring())    
