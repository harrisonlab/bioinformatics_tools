#!/home/gomeza/miniconda3/envs/gene_pred/bin/python

'''
This program extracts all the GO annotations for each gene.
'''


#-----------------------------------------------------
# Step 1
# Import variables & load input files
#-----------------------------------------------------

import sys
import argparse
import re
import numpy
#my_set = set()
#from sets import Set
from Bio import SeqIO
from collections import defaultdict
from operator import itemgetter

ap = argparse.ArgumentParser()
ap.add_argument('--interpro',required=True,type=str,help='The genome assembly')


conf = ap.parse_args()

with open(conf.interpro) as f:
    interpro_lines = f.readlines()

annot_dict = defaultdict(list)

for line in interpro_lines:
    line = line.rstrip()
    split_line = line.split("\t")
    gene_id = split_line[0]
    m = re.findall("GO:.......", line)
    if m:
        annot_dict[gene_id].extend(m)

# max_annots = 0
# for gene_id in annot_dict.keys():
#     annotation_set = set(annot_dict[gene_id])
#     #count number of terms in list
#     num_annots = len(annotation_set)
#     if num_annots >= max_annots:
#         max_annots = num_annots

# additional_columns = [""] * max_annots

for gene_id in annot_dict.keys():
    annotation_set = set(annot_dict[gene_id])
    # num_annots = len(annotation_set)
    # print gene_id + "\t" + "\t".join(annotation_set) + "\t" + "\t".join(additional_columns[num_annots:])
    print (gene_id + "\t" + ", ".join(annotation_set))

# print annot_dict