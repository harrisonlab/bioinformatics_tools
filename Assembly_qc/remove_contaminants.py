#!/home/gomeza/miniconda3/envs/Python/bin/python

'''
This tools will split or remove contigs from an assembly based upon instructions
given in a FCSreport file by NCBI.

The tab seperated co-ordinate file may look like:

SUBID           BioProject      BioSample       Organism
--------------------------------------------------------
SUB1938501      PRJNA338256     SAMN05529104    Fusarium proliferatum

***Text on errors detected***

Exclude:
Sequence name, length, apparent source

contig_705      843     mitochondrion
contig_744      722     vector/etc
contig_103      115415  Delftia acidovorans SPH-1
contig_115      92435   Delftia acidovorans SPH-1

Trim:
Sequence name, length, span(s), apparent source
contig_137      66900   1..69   adaptor:NGB00755.1
contig_447      4231    1..1283 mitochondrion
contig_116      91572   1..73850,76685..76732,77137..77170,79043..79074 Delftia acidovorans SPH-1


The important features of this file are that contig names are in the first column.
The third column contains the instruction to exclude or split the contig.
The fourth and fith columns contain the co-ordinates where a split should be
performed.
The 2nd and 6th+ columns can contain anything. This format is based upon the
format ncbi returns information on which contigs need to be excluded / split in
WGS projects they are processing.
'''

from Bio import SeqIO
import sys,argparse
import sys
from os import path
from collections import defaultdict



#######################################
#            Import variables         #
#######################################

ap = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
ap.add_argument('--inp',required=True,type=str,help='fasta file of assembled contigs')
ap.add_argument('--out',required=True,type=str,help='output contigs file')
ap.add_argument('--coord_file',required=True,type=str,help='file giving the co-ordinates for trimming / exclusion of contigs')
ap.add_argument('--keep_mitochondria',required=False,default=False,action="store_true",help='features to be retained by this script')
conf = ap.parse_args()
out_f = open(conf.out, 'w')
keep_mt = conf.keep_mitochondria

exclude_dic = {}
split_dic = {}
length_dic = defaultdict(list)
len_list = []

exclude_lines = False
trim_lines = False
comments_following = False

#######################################
#        Open co-ordinate file        #
#######################################

f1 = open(conf.coord_file)

print("Processing the instructions file:\n")
for line in f1:
    line = line.strip()
    this_rec = line.split("\t")
    rec_name = this_rec[0]
    if line.startswith("Exclude:"):
        print ("Excluding lines:")
        exclude_lines = True
    elif line.startswith("Trim:"):
        print ("Trimming lines")
        trim_lines = True
    elif line.startswith("Sequence name,"):
        comments_following = True
    elif line == "":
        exclude_lines = False
        trim_lines = False
        comments_following = False
    elif comments_following == True and exclude_lines == True:
        if this_rec[2] == "mitochondrion" and keep_mt == True:
            print ("Ignoring instructions for putative mitochondrial sequence")
        else:
            exclude_dic[rec_name] = this_rec
    elif comments_following == True and trim_lines == True:
                if this_rec[3] == "mitochondrion" and keep_mt == True:
                    print ("Ignoring instructions for putative mitochondrial sequence")
                else:
                    split_dic[rec_name] = this_rec

print("\nThe following instuctions have been logged:\n")

print("exclude dictionary contains:\t") ,
print(len(exclude_dic.keys()))
print("\n")
print("split dictionary contatins:\t") ,
print(len(split_dic.keys()))
print("\n")

#######################################
#        Identify contigs for         #
#          editing and store          #
#          for later printing         #
#######################################


with open(conf.inp) as f2:
    for rec in SeqIO.parse(f2, "fasta"):
        if rec.id in exclude_dic:
                print("Excluding record:\t" + rec.id)
                # Do nothing
        elif rec.id in split_dic:
            split_list = split_dic[rec.id][2].split(",")
            goodseq_boundary_end = []
            goodseq_boundary_start = []
            for split in split_list:
                trim_coordinates = split.split('..')
                goodseq_boundary_end.append(trim_coordinates[0])
                goodseq_boundary_start.append(trim_coordinates[1])
            goodseq_boundary_start.insert(0,int(0))
            goodseq_boundary_end.append(len(rec.seq))
            print ("Split this accession:\t") ,
            print (rec.id)
            print (split_dic[rec.id])
            for i in xrange(0,len(goodseq_boundary_start)):
                start = int(goodseq_boundary_start[i])
                end = int(goodseq_boundary_end[i]) -1 # NCBI boundaries are inclusive of the bases that they want trimmed
                seq = str(rec.seq)[start : end]
                seq = seq.strip('N')
                length = len(seq)
                length_dic[length].append(seq)
                len_list.append(length)
        else:
            unmod_seq = str(rec.seq)
            unmod_seq = unmod_seq.strip('N')
            unmod_seq_len = len(unmod_seq)
            len_list.append(len(unmod_seq))
            length_dic[unmod_seq_len].append(unmod_seq)

#######################################
#            Sort contigs             #
#######################################

sorted_len = sorted(set(len_list), reverse=True)

#######################################
#           Rename contigs            #
#######################################

i = 0
for length in sorted_len:
    if int(length) >= int(500):
        for seq in length_dic[int(length)]:
            i += 1
            out_f.write(">contig_" + str(i) + "\n" + seq + "\n")

print ("number of contigs written:\t" + str(i))
quit()
