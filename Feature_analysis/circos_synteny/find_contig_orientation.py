#!/usr/bin/python

'''
This program looks at link-files used in circos plots and determines whether
contigs should be orientated in forward of reverse orientation between them
'''

#-----------------------------------------------------
# Step 1
# Import variables & load input files
#-----------------------------------------------------

import sys
import argparse
import re
from collections import defaultdict
from collections import Set

ap = argparse.ArgumentParser()
# ap.add_argument('--reference_col',required=True,type=int,choices=['1','2'],help='Column containing the columns you want other contigs to be ordered to')
ap.add_argument('--links_file',required=True,type=str,help='Input links file ordered by the reference contigs')

conf = ap.parse_args()

with open(conf.links_file) as f:
    links_lines = f.readlines()

prev_contig = ''
prev_link_start = ''
forward_list = []
reverse_list = []
seen_list = []
seen_set = set()
for line in links_lines:
    line = line.rstrip()
    split_line = line.split()
    ref_contig = split_line[3]
    this_contig = split_line[0]
    this_link_start = split_line[1]
    if this_contig in seen_set:
        continue
    if prev_contig == this_contig:
        if this_link_start >= prev_link_start:
            forward_list.append(this_contig)
        elif this_link_start <= prev_link_start:
            reverse_list.append(this_contig)
        seen_set.add(this_contig)
        seen_list.append(this_contig)
    prev_contig = this_contig
    prev_link_start = this_link_start

print "Note - This program gives indication of contig order based upon the first two links it encounters of a particular contig. Check all orientations for validity."

print "Order of all seen contigs"
print ", ".join(seen_list)
print "Same orientation as reference contigs:"
print ", ".join(forward_list)
print "Reverse orientation to reference contigs:"
print ", ".join(reverse_list)



# if conf.reference_col == 1:
# elif conf.reference_col == 2:
