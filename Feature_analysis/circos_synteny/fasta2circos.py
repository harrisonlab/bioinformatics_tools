#!/usr/bin/python

'''
This program is used to convert fasta files into input format for circos
'''


#-----------------------------------------------------
# Step 1
# Import variables & load input files
#-----------------------------------------------------

import sys
import argparse
import re
import numpy
from sets import Set
from Bio import SeqIO
from collections import defaultdict
from operator import itemgetter

ap = argparse.ArgumentParser()
ap.add_argument('--genome',required=True,type=str,help='The genome assembly')
ap.add_argument('--contig_prefix',required=True,type=str,help='The prexix for contigs names')

conf = ap.parse_args()

#-----------------------------------------------------
# Step 2
# Identify the length and gene density of FoL chromosomes
#-----------------------------------------------------

prefix = conf.contig_prefix
i = 0

contig_length_dict = defaultdict(list)

genome_file = open(conf.genome, 'r')
for cur_record in SeqIO.parse(genome_file,"fasta"):
    seq_id = cur_record.id
    i += 1
    seq_len = len(cur_record.seq)
    outline = " ".join(["chr", "-", str(prefix) + str(seq_id), str(i), "0", str(seq_len), "chr" + str(i)])
    print(outline)
