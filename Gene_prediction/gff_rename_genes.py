#!/mnt/shared/scratch/agomez/apps/conda/envs/perly_env/bin/python

'''
Script to rename genes in gff files. Genes are renamed from g1 based on the order
of appearence.
'''

import sys,argparse
from collections import defaultdict
#from sets import Set
my_set = set()

#-----------------------------------------------------
# Step 1
# Import variables & load input files
#-----------------------------------------------------
ap = argparse.ArgumentParser()
ap.add_argument('--inp_gff',required=True,type=str,help='input gff file')
ap.add_argument('--conversion_log', required = False, type=str, default = False, help = 'File name to store converted gene names if required')
conf = ap.parse_args() #sys.argv

if conf.conversion_log:
    logfile = conf.conversion_log

with open(conf.inp_gff) as f:
    inp_lines = f.readlines()


#-----------------------------------------------------
# Step 2
# Order gff features by input contigs and gene start
#-----------------------------------------------------

contig_list = []
gene_start_dict = defaultdict(list)
features_dict = defaultdict(list)
conversion_dict = defaultdict(list)

for line in inp_lines:
    line = line.strip("\n")
    if line.startswith("#"):
        print (line)
        continue
    split_line = line.split()
    # Col9 = split_line[8]
    if split_line[2] == 'gene':
        contig = split_line[0]
        gene_start = split_line[3]
        if contig not in contig_list:
            contig_list.append(contig)
        gene_start_dict[contig].append(int(gene_start))
        key = "_".join([contig, gene_start])
    features_dict[key].append(line)

gff_lines = []
# for contig in contig_list:
for contig in sorted(contig_list, key = lambda x: (int(x.split('_')[1]))):
    gene_start_list = gene_start_dict[contig]
    for gene_start in sorted(gene_start_list):
        # print "\t".join([contig, str(gene_start)])
        key = "_".join([contig, str(gene_start)])
        gff_lines.extend(features_dict[key])
        # print("\n".join(features_dict[key]))


#-----------------------------------------------------
# Step 3
# Raname genes
#---

split_line = []
outlines = []
gene_set = set([])
gene_dict = defaultdict(list)
i = int(0)
conversion_lines = []
# for line in inp_lines:
for line in gff_lines:
    line = line.strip("\n")
    # if line.startswith("#"):
    #     print(line)
    #     continue
    split_line = line.split()
    Col9 = split_line[8]
    if split_line[2] == 'gene':
        gene_id = Col9.split("=")[1]
        gene_id = gene_id.replace(";", "")
        i += 1
        new_id = "g" + str(i)
        conversion_lines.append("\t".join([gene_id, "-", new_id]) + "\n")
    Col9 = Col9.replace(gene_id, new_id)
    split_line[8] = Col9
    print ("\t".join(split_line))

if logfile:
    f = open(logfile,"w")
    for line in conversion_lines:
        f.write(line)
    f.close()
