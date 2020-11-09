#!/usr/bin/python

'''
This program is used to build information on all the genes predicted in
a fusarium genome. These commands take information on location of blast
hits against reference genome & suppliment this information with information
from other files detailing functional annotations and orthology status.
'''

#-----------------------------------------------------
# Step 1
# Import variables & load input files
#-----------------------------------------------------

import sys
import argparse
import re
from sets import Set
from collections import defaultdict
from operator import itemgetter

ap = argparse.ArgumentParser()
ap.add_argument('--blast_tbl',required=True,type=str,help='The .tbl output of a swissprot blast')
ap.add_argument('--blast_db_fasta',required=True,type=str,help='The original fasta file from which the blast database was constructed')

conf = ap.parse_args()

with open(conf.blast_tbl) as f:
    blast_tbl_lines = f.readlines()

with open(conf.blast_db_fasta) as f:
    fasta_lines = f.readlines()



#-----------------------------------------------------
# Step 2
# Build a database of the original fasta files.
# This will be used to link Blast hits to function.
#-----------------------------------------------------

fasta_id_set = Set([])
fasta_dict = defaultdict(list)
for line in fasta_lines:
    line = line.rstrip()
    if line.startswith('>'):
        # print(line)
        split_line = line.split("|")
        # db_info = split_line[0]
        db_id = split_line[1]
        # print (db_id)
        functional_info = split_line[-1]
        functional_info_split = functional_info.split(" ")
        # OS_index_list = [i for i, s in enumerate(functional_info_split) if 'OS=' in s]
        # OS_index = "".join(str(x) for x in OS_index_list)
        index_list = [i for i, s in enumerate(functional_info_split) if '=' in s]
        OS_index = int(index_list[0])
        # GN_index = "".join(str(x) for x in GN_index_list[])
        GN_index = int(index_list[1])
        # if len(GN_index_list) == 0:
        #     print line
        # print (OS_index)
        gene_function = " ".join(functional_info_split[1: int(OS_index)])
        # print (gene_function)
        species_name = " ".join(functional_info_split[int(OS_index) : int(GN_index)])
        species_name = species_name.replace('OS=', '')
        # print (species_name)

        fasta_dict[db_id]=[db_id, gene_function, species_name]



#-----------------------------------------------------
# Step 2
# Build a database of the original fasta files.
# This will be used to link Blast hits to function.
#-----------------------------------------------------

gene_id_set = Set([])
# fasta_dict = defaultdict(list)
for line in blast_tbl_lines:
    line = line.rstrip()
    split_line = line.split("\t")
    gene_id = split_line[0]
    if gene_id not in gene_id_set:
        id_col = split_line[1]
        id_col_split = id_col.split("|")
        db_id = id_col_split[1]
        gene_id_set.add(gene_id)
        split_line.extend(fasta_dict[db_id])
        print ("\t".join(split_line))
