#!/home/gomeza/miniconda3/envs/perly_env/bin/python

'''
This tool uses the output of signalP4.1 to extract SigP positive accessions from
a fasta file and add the appropriate location of cleaveage sites along with
hmm scores for signal peptides.
'''

import sys,argparse
import re
#from sets import Set
my_set = set()

#-----------------------------------------------------
# Step 1
# Import variables & load input files
#-----------------------------------------------------


ap = argparse.ArgumentParser()
ap.add_argument('--inp_fasta',required=True,type=str,help='input FASTA file')
ap.add_argument('--inp_sigP',required=True,type=str,help='File containing SignalP4.1 results')
ap.add_argument('--out_tab',required=True,type=str,help='A tab seperated file of genenames and SigP results')
ap.add_argument('--out_fasta',required=True,type=str,help='Output fasta of SigP-containing proteins')
ap.add_argument('--out_neg',required=True,type=str,help='Output fasta of SigP-nedative proteins')
conf = ap.parse_args()

pos_genes = {}
tab_lines = []
hit_fasta_lines = []
miss_fasta_lines = []


with open(conf.inp_sigP) as sigP_f:
    sigP_lines = sigP_f.readlines()

with open(conf.inp_fasta) as fasta_f:
    fasta_lines = fasta_f.readlines()


#-----------------------------------------------------
# Step 2
# Identify SigP hits within sigP file
# Store key information from hits and non-hits.
#-----------------------------------------------------

for line in sigP_lines:
    if 'Name' in line:
        m = re.search(r'Name=(.*?)\sSP=\'(.*?)\'\s(.*?)D=(.*?)\sD-cutoff=(.*?)\sNetworks=(.*)', line)
        gene_name = m.group(1)
        hit = m.group(2)
        cleavage = m.group(3)
        score = m.group(4)
        cutoff = m.group(5)
        trans_mem = m.group(6)

        if hit == 'YES':
            cleave_m = re.search(r'pos\.\s(\d+)', cleavage)
            cleavage_site=cleave_m.group(1)
            headerline="\t".join([">" + gene_name, "--HMM_score", score, "--Signal_peptide_length=", cleavage_site + "\n"])
            this_tab_line="\t".join([gene_name, hit, score, cleavage_site])
            tab_lines.append(this_tab_line)
            pos_genes[gene_name] = headerline
        if hit == 'NO':
            cleavage_site = "-"
            this_tab_line="\t".join([gene_name, hit, score, cleavage_site])
            tab_lines.append(this_tab_line)


#-----------------------------------------------------
# Step 3
# For accessions in fasta file with same name as
# SigP-hits, print to hits-outfile with the accompanying
# hmm-score and cleavage site information
#-----------------------------------------------------

for line in fasta_lines:
    if '>' in line:
        printline = False
        m = re.search(r'>(.*?)\s', line)
        header = str(m.group(1))
        # print header
        if header in pos_genes:
            new_header = pos_genes[header]
            hit_fasta_lines.append(new_header)
            printline = True
        else:
            miss_fasta_lines.append(line)
    elif printline == True:
        hit_fasta_lines.append(line)
    else:
        miss_fasta_lines.append(line)

#-----------------------------------------------------
# Step 4
# Print outfiles
#-----------------------------------------------------

with open(conf.out_tab, 'w') as o:
    o.write("Query\tname\tSignal peptide probability(SignalP-HMM)\tCleavage site\n")
    for out_tab_line in tab_lines:
        o.write(">" + out_tab_line + "\n")

with open(conf.out_fasta, 'w') as o:
    for out_line in hit_fasta_lines:
        o.write(out_line)

with open(conf.out_neg, 'w') as o:
    for out_line in miss_fasta_lines:
        o.write(out_line)
