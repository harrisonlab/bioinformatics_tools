#!/usr/bin/python


'''
This tool converts gff files into a tab-based format that foramat_embly.py can
use, along with a fasta file to produce an EMBL format file of the genome and
its annotations. THis EMBL file can then be provided to antismash for prediction
of antimsash clusters. This program requires a gff file as input.
'''

#-----------------------------------------------------
# Step 1
# Import variables & load input files
#-----------------------------------------------------

import sys,argparse
import re
# from natsort import natsorted, ns
# from Bio.Seq import Seq
# from Bio.SeqFeature import SeqFeature, FeatureLocation
from sets import Set
from collections import defaultdict
from operator import itemgetter

ap = argparse.ArgumentParser()
ap.add_argument('--fasta',required=True,type=str,help='.fasta file of the assembly')
# ap.add_argument('--regions',required=True,type=str,help='.gff file of gene annotations')
ap.add_argument('--out',required=True,type=str,help='output prefix')

conf = ap.parse_args()

# outdir = conf.out


#-----------------------------------------------------
# Step X
# Define classes
#-----------------------------------------------------
def reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N':'N'}
    return ''.join([complement[base] for base in dna[::-1]])

class Fasta_obj(object):
    """A gene identified as differentially expressed in one isolate.
    Attributes:
        gene_name: A string representing the gene name.
        conditions_tested: List of conditions that the gene was tested for DE.
        conditions_positive: A binary list of integers representing whether the
            gene tested positive for each condition listed in conditions tested.
    """

    def __init__(self):
        """initiate a fasta object"""
        self.fasta_id = ''
        self.seq = ''
        self.seq_1kb_sections_F = []
        self.seq_1kb_sections_R = []
        self.telomere_hits_F = []
        self.telomere_hits_R = []


    def search_telomeres(self):
        """Reset conditions_tested to a list of conditions"""
        n = 1000
        self.seq_1kb_sections_F = [self.seq[i:i+n] for i in range(0, len(self.seq), n)]
        for chunk in self.seq_1kb_sections_F:
            # self.telomere_hits.append(len(re.findall('TTAGGG'|'AATCCC'|'CCCTAA', chunk)))
            # self.telomere_hits.append(len(re.findall('TTAGGG', chunk)))
            self.telomere_hits_F.append(len(re.findall(r'TTAGGG|CCCTAA', chunk)))
        # rev_seq = reverse_complement(self.seq)
        # self.seq_1kb_sections_R = [rev_seq[i:i+n] for i in range(0, len(rev_seq), n)]
        # for chunk in self.seq_1kb_sections_R:
        #     # self.telomere_hits.append(len(re.findall('TTAGGG'|'AATCCC'|'CCCTAA', chunk)))
        #     # self.telomere_hits.append(len(re.findall('TTAGGG', chunk)))
        #     self.telomere_hits_R.append(len(re.findall(r'TTAGGG|CCCTAA', chunk)))



#-----------------------------------------------------
# Step 1
# Prepare an individual fasta file for each assembled contig
#-----------------------------------------------------

with open(conf.fasta) as f:
    fasta_lines = f.readlines()
f.close()

seq_dict = defaultdict()

for line in fasta_lines:
    line = line.rstrip()
    if line.startswith('>'):
        ID = line.replace('>', '').split(" ")[0]
        seq_dict[ID] = Fasta_obj()
        seq_dict[ID].fasta_id = ID
    else:
        seq_dict[ID].seq += line.upper()

outlines = []
out_circos = []
for key in seq_dict.keys():
    seq_dict[key].search_telomeres()
    # hits_list = seq_dict[key].telomere_hits_F[0:10]
    hits_list = seq_dict[key].telomere_hits_F
    for i, num_hits in enumerate(hits_list):
        outlines.append("\t".join([key, 'F', str(i + 1) + '000', str(len(hits_list)) + '000', str(num_hits)]))
        if num_hits > 0:
            out_circos.append("\t".join([key, str(i) + '000', str(i + 1) + '000', str(num_hits)]))
        # print "\t".join([key, 'F', '+' + str(i + 1) + '000', str(num_hits)])
    # hits_list = seq_dict[key].telomere_hits_F[-10:]
    # for i, num_hits in enumerate(hits_list):
    #     print "\t".join([key, 'F', '-' + str(10 - i) + '000', str(num_hits)])
    # hits_list = seq_dict[key].telomere_hits_R[-5:]
    # for i, num_hits in enumerate(hits_list):
    #     print "\t".join([key, 'R', '+' + str(i + 1) + '000', str(num_hits)])

# print ID
# print line

o = open(conf.out + '.txt', 'w')
o.write("\n".join(outlines))

o = open(conf.out + '_circos.txt', 'w')
o.write("\n".join(out_circos))
