#!/home/gomeza/miniconda3/envs/gene_pred/bin/python

'''
This tool splits a fasta file with many accessions down into a number of
smaller files, each contatining 500 accessions.
'''

import sys,argparse
from Bio import SeqIO

ap = argparse.ArgumentParser()
ap.add_argument('--inp_fasta',required=True,type=str,help='input FASTA file')
ap.add_argument('--out_dir',required=True,type=str,help='the directory for split files to be written to')
ap.add_argument('--out_base',required=True,type=str,help='the basename for output files')
conf = ap.parse_args()

f = open(conf.inp_fasta)
out_dir = conf.out_dir
n = conf.out_base

fasta_sequences = SeqIO.parse(f,'fasta')


i = 0
file_num = 500
new_file = out_dir + "/" + n + "_" + str(file_num) + ".fa"
cur_file=open(new_file, 'w')
for seq_record in fasta_sequences:
    i += 1
    if i <= 500:
        cur_file.write(seq_record.format("fasta"))
    else:
        i = 0
        file_num += 500
        cur_file.close()
        new_file = out_dir + "/" + n + "_" + str(file_num) + ".fa"
        cur_file = open(new_file, 'w')
        cur_file.write(seq_record.format("fasta"))

cur_file.close()
f.close()
