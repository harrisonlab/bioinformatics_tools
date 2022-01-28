#!/usr/bin/env python

'''
Provide a measure of mean coverage in a sliding window of a specified size
'''

import sys,argparse
import numpy as np
from collections import defaultdict

#######################################
#            Import variables         #
#                                     #
#                                     #
#######################################

ap = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.ArgumentDefaultsHelpFormatter)

ap.add_argument('--cov',required=True,type=str,help='Summary of reads depth by bp as output by samtools depth and pasted into a single file using BASH paste')
# ap.add_argument('--size',required=True,type=int,help='The size of the sliding window')

conf = ap.parse_args() #sys.argv

line_list = []
total_cov = 0
prev_contig = 'first'
prev_pos = 0
i = 0
with open(conf.cov) as f:
    for line in f:
        split_line = line.split("\t")
        contig = split_line[0]
        contig_pos = split_line[1]
        cov = int(split_line[2])
        i +=1
        if contig != prev_contig:
            if prev_contig == 'first':
                total_cov += cov
                prev_contig = contig
                continue
            depth = np.divide(total_cov, i-1)
            outline = [prev_contig, prev_pos, str(depth)]
            print "\t".join(outline)
            total_cov = cov
            i = 1
        elif i < 1000:
            total_cov += cov
        elif i == 1000:
            total_cov += cov
            depth = np.divide(total_cov, i)
            outline = [contig, contig_pos, str(depth)]
            print "\t".join(outline)
            i = 0
            total_cov = 0
        prev_contig = contig
        prev_pos = contig_pos

depth = np.divide(total_cov, i)
outline = [contig, contig_pos, str(depth)]
print "\t".join(outline)
