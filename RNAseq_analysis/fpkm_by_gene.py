#!/usr/bin/python

'''
This tool ranks genes by FPKM as output by coverage quantification (Cufflinks -G).
The program accepts multiple cufflinks output files and determines final rank
of gene expression based upon mean expression accross these repeats.
'''

import sys,argparse
from collections import defaultdict
from sets import Set
import numpy as np

#######################################
#            Import variables         #
#                                     #
#                                     #
#######################################

ap = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.ArgumentDefaultsHelpFormatter)

ap.add_argument('--fpkm_file',required=True,type=str,nargs='+',help='One or more cufflinks gene_fpkm_tracking files')
ap.add_argument('--CAZY_headers',required=False,type=str,help='file containing headers of CAZY genes')
ap.add_argument('--effectorP_headers',required=False,type=str,help='file containing headers of effectorP genes')
ap.add_argument('--transposon_headers',required=False,type=str,help='file containing headers of identified as transposon genes')
ap.add_argument('--metabolite_headers',required=False,type=str,help='file containing headers of identified as secondary metabolite genes')
ap.add_argument('--annotation_table',required=False,type=str,help='The gene annotation file for this organism')


conf = ap.parse_args() #sys.argv
fpkm_dict = defaultdict(list)
CAZY_set = Set()
effP_set = Set()
transposon_set = Set()
metabolite_set = Set()
annotation_dict = defaultdict(list)

#######################################
#            Import variables         #
#                                     #
#                                     #
#######################################

if conf.CAZY_headers is not None:
    with open(conf.CAZY_headers) as f:
        for line in f:
            line = line.rstrip()
            if line[-2] == 'T' and line[-1].isdigit:
                CAZY_gene = line[:-2]
            else:
                CAZY_gene = line.split(".")[0]
            CAZY_set.add(CAZY_gene)

if conf.effectorP_headers is not None:
    with open(conf.effectorP_headers) as f:
        for line in f:
            line = line.rstrip()
            # print line
            if line[-2] == 'T' and line[-1].isdigit:
                effP_gene = line[:-2]
                # print effP_gene
            else:
                effP_gene = line.split(".")[0]
            effP_set.add(effP_gene)

if conf.transposon_headers is not None:
    with open(conf.transposon_headers) as f:
        for line in f:
            line = line.rstrip()
            if line[-2] == 'T' and line[-1].isdigit:
                transposon_gene = line[:-2]
            else:
                transposon_gene = line.split(".")[0]
            transposon_set.add(transposon_gene)

if conf.metabolite_headers is not None:
    with open(conf.metabolite_headers) as f:
        for line in f:
            line = line.rstrip()
            if line[-2] == 'T' and line[-1].isdigit:
                metabolite_gene = line[:-2]
            else:
                metabolite_gene = line.split(".")[0]
            metabolite_set.add(metabolite_gene)

if conf.annotation_table is not None:
    with open(conf.annotation_table) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith("query_id"):
                continue
            split_line = line.split("\t")
            transcript_id = split_line[0]
            # print transcript_id
            if transcript_id[-2] == 'T' and transcript_id[-1].isdigit:
                # print transcript_id
                gene_id = transcript_id[:-2]
                # print gene_id
            else:
                gene_id = transcript_id.split(".")[0]
            annotation_dict[gene_id].append(line)

file_list = conf.fpkm_file
for fpkm_file in file_list:
    with open(fpkm_file) as f:
        for line in f:
            if line.startswith("tracking_id"):
                continue
            line = line.rstrip()
            split_line = line.split("\t")
            gene_id = split_line[4] #4 for 5th column
            fpkm = float(split_line[9])
            fpkm_dict[gene_id].append(fpkm)
            

list_of_tupule = []
for gene in fpkm_dict.keys():
    fpkm_list = fpkm_dict[gene]
    mean_fpkm = round(np.mean(fpkm_list))
    stdev_fpkm = round(np.std(fpkm_list))
    list_of_tupule.append([gene, mean_fpkm, stdev_fpkm])

def getKey(item):
    return item[1]

i = 0
for gene, mean_fpkm, stdev_fpkm in sorted(list_of_tupule, key=getKey, reverse=True):
    i +=1

    feature_list = []
    if gene in CAZY_set:
        feature_list.append("CAZY")
    if gene in effP_set:
        feature_list.append("EffP")
    if gene in transposon_set:
        feature_list.append("Transposon")
    if gene in metabolite_set:
        feature_list.append("Antismash")
    features = "; ".join(feature_list)
    gene_annotations = "\t".join(annotation_dict[gene])
    outline = [str(i), str(gene), str(mean_fpkm), str(stdev_fpkm), features, gene_annotations]
    print "\t".join(outline)
    # print str(i) + "\t" + str(gene) + "\t" + str(mean_fpkm)
