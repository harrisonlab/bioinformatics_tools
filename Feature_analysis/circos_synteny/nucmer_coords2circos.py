#!/usr/bin/python

'''
A tool to convert NUCmer output files into genome alignment features in Circos format.
'''

import sys,argparse
import re
from collections import defaultdict

#-----------------------------------------------------
# Step 1
# Import variables & load input files
#-----------------------------------------------------


ap = argparse.ArgumentParser()
ap.add_argument('--inp_coords',required=True,type=str,help='Coordinates')
ap.add_argument('--queery_id',required=True,type=str,help='Queery ID to add to contig names in the first column')
ap.add_argument('--ref_id',required=True,type=str,help='Reference ID to add to contig names in the second column')
conf = ap.parse_args()

pos_genes = {}
tab_lines = []
hit_fasta_lines = []
miss_fasta_lines = []


with open(conf.inp_coords) as f:
    coords_lines = f.readlines()

q_id = conf.queery_id
r_id = conf.ref_id

#-----------------------------------------------------
# Step 2
# Parse input columns following the NUCmer header line
# Output column information in Circos format
#-----------------------------------------------------

start_reading = False
for line in coords_lines:
    line = line.rstrip()
    if line.startswith('[S1]'):
        start_reading = True
    elif start_reading == True:
        split_line = line.split('\t')
        r_start = split_line[0]
        r_end = split_line[1]
        q_start = split_line[2]
        q_end = split_line[3]
        r_contig = "_".join([r_id, split_line[10]])
        q_contig = "_".join([q_id, split_line[11]])
        print "\t".join([q_contig, q_start, q_end, r_contig, r_start, r_end])
