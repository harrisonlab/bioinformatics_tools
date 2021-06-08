#!/usr/bin/python

'''
This script will run parse the ouput file from Optimus into a csv file that can
be read into excel for analysis of potential CrisprCas sites by the end user.
'''

import os
import sys,argparse
import re
from collections import defaultdict
# import gffutils
# from collections import defaultdict
# from itertools import chain


#######################################
#           Load sys. args.           #
#                                     #
#                                     #
#######################################


ap = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
ap.add_argument('--inp',required=True,type=str,help='The output file from OPTIMus')
ap.add_argument('--out',required=True,type=str,help='The name of the tab seperated file this script will output')

conf = ap.parse_args() #sys.argv
f = open(conf.inp)
out_f = open(conf.out, 'w')

print("In file is:\t" + conf.inp)

#######################################
#           Open the input file.      #
#                                     #
#                                     #
#######################################
'''
Read lines of input until the next protospacer is found.
For each gene within that protospacer group do the following:
1) add an accession to the dictionary.
2) Collect the co-ordinates of each occurence of the protospacer sequence
3) append to the dictionary
    \tprotospacer sequence:co-ordinates:other genes in this protospacer group.
'''

global_dict = defaultdict(list)
global_genes = []

def group_by_heading( source_file ):
    buffer= []
    for line in source_file:
        line = line.rstrip()
        # print(line)
        if line.startswith( ">" ):
            # print(line)
            buffer.append( line )
        else:
            # print("Newline\n")
            if buffer: yield buffer
            buffer= [ line ]
    yield buffer

def natural_sort(l):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    return sorted(l, key = alphanum_key)


for heading_and_lines in group_by_heading(f):
    # print(heading_and_lines)
    group_dict = defaultdict(list)
    gene_list = []
    protospacer = heading_and_lines[0]
    for line in heading_and_lines[1:]:
        gene = re.sub('_coor.*$', '', line)
        coordinates = re.sub('^.*_coor-', '', line)
        group_dict[gene].append(coordinates)
        gene_list.append(gene)
    sorted_genes = natural_sort(set(gene_list))
    # print(sorted_genes)
    for gene in sorted_genes:
        filtered_genes = [ g for g in sorted_genes if not g.startswith(gene) ]
        if not filtered_genes:
            filtered_genes = ['unique']
        coor_list = group_dict[gene]
        # print(coor_list)
        coor_list.sort()
        positions = ','.join(coor_list)
        this_entry = protospacer + ":" + positions + ":" + ",".join(filtered_genes)
        global_dict[gene].append(this_entry)

    group_dict.clear()
    # print(gene + "\t" + "\t".join(global_dict[gene]) + "\n")
    # exit()



all_genes = global_dict.keys()
sorted_all_genes = natural_sort(all_genes)

print("Writing to:\t" + conf.out)
for gene in sorted_all_genes:
    out_f.write(gene + "\t" + "\t".join(global_dict[gene]) + "\n")
