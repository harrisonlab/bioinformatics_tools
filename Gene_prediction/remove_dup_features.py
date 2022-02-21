#!/home/gomeza/miniconda3/envs/perly_env/bin/python

'''
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
ap.add_argument('--out_gff',required=True,type=str,help='output gff file')
conf = ap.parse_args() #sys.argv

with open(conf.inp_gff) as f:
    inp_lines = f.readlines()

out_f = conf.out_gff

# def check_in_dict(contig, cds_boudnries):
#     gene_key = "\t".join(cds_boudnries)

transcript = ""
contig = ""
cds_boudnries = []
split_line = []
outlines = []
gene_set = set()
duplicated_set = set()
gene_dict = defaultdict(list)
for line in inp_lines:
    line = line.strip("\n")
    if line.startswith("#"):
        continue
    split_line = line.split()
    if split_line[2] == 'mRNA':
        # print line
        prev_transcript = transcript
        col9 = split_line[8]
        transcript = col9.split(';')[0].replace('ID=', '')
        key = contig +  "_".join(cds_boudnries)
        if key in gene_set:
            duplicated_set.add(prev_transcript)
        gene_set.add(key)
        cds_boudnries = []
    elif split_line[2] == 'CDS':
        contig = split_line[0]
        start = split_line[3]
        end = split_line[4]
        cds_boudnries.extend([start, end])
        # gene_dict[gene_key].append(line)
        # # if gene_key in gene_set:
        # if len(gene_dict[gene_key]) > 1:
        #     print "Duplicate gene found:\t" + gene_key
        #     # print line
        #     print "\n".join(gene_dict[gene_key])
        # # else:
        #     # gene_set.add(gene_key)
        #     # gene_dict[gene_key].append(line)
prev_transcript = transcript
col9 = split_line[8]
transcript = col9.split(';')[0].replace('ID=', '')

key = contig +  "_".join(cds_boudnries)
if key in gene_set:
    duplicated_set.add(prev_transcript)
gene_set.add(key)

# print "\n".join(gene_set)
print ("Identifiied the following duplicated transcripts:")
print ("\n".join(duplicated_set))
print ("NOTE - if any of these represent the first transcript of a gene (.t1) then an entire gene may be duplicated")
print ("if so the gene feature will need to be stripped out of the gff file seperately.")


outlines = []
for line in inp_lines:
    line = line.rstrip("\n")
    if not any(gene in line for gene in duplicated_set):
        outlines.append(line)
f = open(out_f,"w")
f.write("\n".join(outlines))
