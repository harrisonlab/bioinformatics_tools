#!/home/gomeza/miniconda3/envs/gene_pred/bin/python


'''
This tool extracts and parses regions predicted to contain secondary metabolite
synthesis genes by SMURF. These regions are output as gff features.
'''

#-----------------------------------------------------
# Step 1
# Import variables & load input files
#-----------------------------------------------------

import sys,argparse
import re
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
#from sets import Set
from collections import defaultdict
from operator import itemgetter

ap = argparse.ArgumentParser()
ap.add_argument('--smurf_clusters',required=True,type=str,help='SMURF output file containing secondary metabolite clusters')
ap.add_argument('--smurf_backbone',required=True,type=str,help='SMURF output file containing backbone genes for secondary metabolite clusters')


conf = ap.parse_args()

with open(conf.smurf_clusters) as f:
    cluster_lines = f.readlines()
    # print "monkeys"
f.close()

with open(conf.smurf_backbone) as f:
    backbone_lines = f.readlines()
    # print "monkeys"
f.close()


#-----------------------------------------------------
# Step 2
# Build a dictionary of products for backbone genes
#
#-----------------------------------------------------

product_dict = defaultdict(list)
for line in backbone_lines:
    line = line.rstrip()
    split_line = line.split()
    backbone_id = split_line[0]
    product = split_line[6]
    product_dict[backbone_id].append(product)

#-----------------------------------------------------
# Step 3
# Use the start and stop location of each feature
# labelled as a "cluster" to ubild a gff annotation
# and print to stdout
#
#-----------------------------------------------------

coordinate_set = set([])
new_cluster = False
program = 'SMURF'
contig = 'First'
i = 0
for line in cluster_lines:
    line = line.rstrip()
    if line.startswith("Cluster:"):
        if contig == 'First':
            continue
        # new_cluster = True
        # print "\n".join(cluster_lines)
        i += 1
        cluster_start = min(coordinate_set)
        cluster_end = max(coordinate_set)
        cluster_func = "".join(["ID=", "Cluster_", str(i), ";Notes=",  product, ";"])
        print ("\t".join([contig, program, product, str(cluster_start), str(cluster_end), '.', strand, '.', cluster_func]))
        cluster_lines = []
        coordinate_set = set([])
    elif line.startswith("Backbone_gene_id") or line == "":
        continue
    else:
        split_line = line.split("\t")
        backbone_id = split_line[0]
        gene_id = split_line[1]
        contig = split_line[3]
        # feature = 'SecMet'
        product = "".join(product_dict[backbone_id])
        col_6 = int(split_line[5])
        col_7 = int(split_line[6])
        if col_6 < col_7:
            gene_start = col_6
            gene_end = col_7
            strand = '+'
        elif col_6 > col_7:
            gene_end = col_6
            gene_start = col_7
            strand = '-'
        coordinate_set.add(gene_start)
        coordinate_set.add(gene_end)

# for line in smurf_lines:
#     line = line.rstrip()
