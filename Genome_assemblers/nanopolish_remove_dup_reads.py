#!/usr/bin/python
import sys,argparse
from os import path
from collections import defaultdict



#######################################
#            Import variables         #
#                                     #
#                                     #
#######################################

ap = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
ap.add_argument('--fastq',required=True,type=str,help='fastq file of minion reads')
ap.add_argument('--out',required=True,type=str,help='output contigs file')
conf = ap.parse_args()

header_set = set()
exclude_line = 0

out_f = open(conf.out, 'w')

with open(conf.fastq) as f:
   for line in f:
       line = line.rstrip()
       if exclude_line > 0:
           exclude_line -= 1
           continue
        #    print(line + "\n")
       elif line.startswith("@"):
        #    read = line.split(" ")[0]
        #    runid = line.split(" ")[1]
        #    header = "_".join([read, runid])
            header = line
            if header in header_set:
               exclude_line = 3
               print(line + "\n")
               # next
               continue
            header_set.add(header)
       out_f.write(line + "\n")
