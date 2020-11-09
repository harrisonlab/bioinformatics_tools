#!/usr/bin/python


'''
Tool to parse gene gff files in a tab seperated format that can be
accpeted by the SMURF webserver for secondary metabolite prediction
'''

#-----------------------------------------------------
# Step 1
# Import variables & load input files
#-----------------------------------------------------

import sys,argparse

ap = argparse.ArgumentParser()
ap.add_argument('--gff',required=True,type=str,help='.gff file of gene annotations')
conf = ap.parse_args()

with open(conf.gff) as f:
    gff_lines = f.readlines()
f.close()

#-----------------------------------------------------
# Step 2
# Parse gff annotations into a tab-seperated format
# accepted by SMURF
#-----------------------------------------------------
# Col.
# 1. protein ID
# 2. chromosome/contig
# 3. 5' gene start
# 4. 3' gene stop
# 5. protein name/function/definition (if available)

for line in gff_lines:
    line = line.rstrip()
    split_line = line.split("\t")
    # print split_line
    if split_line[2] and "mRNA" in split_line[2]:
        contig = split_line[0]
        gene_start = split_line[3]
        gene_end = split_line[4]
        col_9 = split_line[8]
        gene_id = col_9.split(';')[0]
        gene_id = gene_id.replace('ID=', '')
        print "\t".join([gene_id, contig, gene_start, gene_end, col_9])
