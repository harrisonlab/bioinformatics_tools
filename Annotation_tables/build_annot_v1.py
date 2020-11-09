#!/usr/bin/python

'''
This script parses information from fasta files and gff files for the location,
sequence and functional information for annotated gene models and RxLRs.
'''

# -----------------------------------------------------
# Step 1
# Import variables & load input files
# -----------------------------------------------------

import argparse
from sets import Set
from collections import defaultdict
from operator import itemgetter
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument('--gff_format', required=True, type=str,
                choices=['gff3', 'gtf'], help='Gff file format')
ap.add_argument('--gene_gff', required=True, type=str,
                help='Gff file of predicyted gene models')
ap.add_argument('--gene_fasta', required=True, type=str,
                help='amino acid sequence of predicted proteins')
#ap.add_argument('--SigP4', required=True, type=str,
                #help='fasta file of genes testing positive for signal peptide \
                #using SigP4.1')
#ap.add_argument('--trans_mem', required=True, type=str,
                #help='txt file of headers from gene testing positive for \
                #transmembrane proteins by TMHMM')
ap.add_argument('--TFs',required=True,type=str,
		        help='Tab seperated of putative transcription factors and their \
                domains as identified by interpro2TFs.py')
#ap.add_argument('--effector_total', required=True, type=str,
                #help='txt file of all transcripts considered low confidence \
                #effector')
#ap.add_argument('--CAZY_total', required=True, type=str,
                #help='txt file of all transcripts considered CAZY')
#ap.add_argument('--ortho_name', required=True, type=str,
                #help='the name used for the organism during orthology \
                #analysis')
#ap.add_argument('--ortho_file', required=True, type=str,
                #help='txt file of ortholog groups')
#ap.add_argument('--DEG_files', required=True, nargs='+', type=str,
                #help='space seperated list of files \
                #containing DEG information')
#ap.add_argument('--raw_counts', required=True, type=str,
               # help='raw count data as output from DESeq')
#ap.add_argument('--fpkm', required=True, type=str,
               # help='normalised fpkm count data as output from DESeq')
ap.add_argument('--InterPro', required=True, type=str,
                help='The Interproscan functional annotation .tsv file')
ap.add_argument('--Swissprot', required=True, type=str,
                help='A parsed table of BLAST results against the Swissprot \
                database. Note - must have been parsed with \
                swissprot_parser.py')
ap.add_argument('--Antismash', required=True, type=str,
                help='Output of Antismash parsed into a tsv file of gene names, \
                contig, secmet function and cluster ID')

conf = ap.parse_args()


with open(conf.gene_gff) as f:
    gene_lines = f.readlines()

with open(conf.gene_fasta) as f:
    prot_lines = f.readlines()

#with open(conf.SigP4) as f:
    #sigP4_lines = f.readlines()

#with open(conf.trans_mem) as f:
    #trans_mem_lines = f.readlines()

#with open(conf.effector_total) as f:
    #effector_total_lines = f.readlines()

#with open(conf.CAZY_total) as f:
    #CAZY_total_lines = f.readlines()

#with open(conf.ortho_file) as f:
    #ortho_lines = f.readlines()

#DEG_files = conf.DEG_files
#DEG_dict = defaultdict(list)
#for DEG_file in DEG_files:
    #with open(DEG_file) as f:
        #filename = DEG_file
        #DEG_lines = f.readlines()
        #for line in DEG_lines:
            #if line.startswith('baseMean'):
                #continue
            #else:
                #split_line = line.split()
                #gene_name = split_line[0]
                #base = split_line [1]
                #log_change = split_line[2]
                #P_val = split_line[6]
                #entryname = "_".join([filename, gene_name])
                #DEG_dict[entryname].extend([base, log_change, P_val])

#with open(conf.raw_counts) as f:
   # raw_count_lines = f.readlines()

#with open(conf.fpkm) as f:
   # fpkm_lines = f.readlines()

with open(conf.InterPro) as f:
    InterPro_lines = f.readlines()

with open(conf.Swissprot) as f:
    swissprot_lines = f.readlines()

with open(conf.TFs) as f:
    TF_lines = f.readlines()

with open(conf.Antismash) as f:
    antismash_lines = f.readlines()

# -----------------------------------------------------
# Load protein sequence data into a dictionary
# -----------------------------------------------------

prot_dict = defaultdict(list)
for line in prot_lines:
    line = line.rstrip()
    if line.startswith('>'):
        header = line.replace('>', '')
    else:
        prot_dict[header] += line

# -----------------------------------------------------
# Load signalP4.1 files into a set
# -----------------------------------------------------

#SigP4_set = Set()
#for line in sigP4_lines:
    #line = line.rstrip()
    #if line.startswith('>'):
        #split_line = line.split()
        #header = split_line[0].replace('>', '')
        #SigP4_set.add(header)

# -----------------------------------------------------
# Load TMHMM headers into a set
# -----------------------------------------------------

#trans_mem_set = Set()
#for line in trans_mem_lines:
    #header = line.rstrip()
    #trans_mem_set.add(header)

# -----------------------------------------------------
# Load Effector total protein into a set
# -----------------------------------------------------

#effector_total_set = Set()
#for line in effector_total_lines:
    #header = line.rstrip()
    #effector_total_set.add(header)
    # line = line.rstrip()
    # if line.startswith('>'):
    #     split_line = line.split()
    #     header = split_line[0].replace('>', '')
    #     RxLR_total_set.add(header)

# -----------------------------------------------------
# Load CAZY total proteins into a set
# -----------------------------------------------------

#CAZY_total_set = Set()
#for line in CAZY_total_lines:
    #header = line.rstrip()
    #if 'contig' in header:
        #header = header + '.t1'
    #CAZY_total_set.add(header)
     #line = line.rstrip()
     #if line.startswith('>'):
         #split_line = line.split()
         #header = split_line[0].replace('>', '')
         #CAZY_total_set.add(header)

# -----------------------------------------------------
# Store genes and their ortholog groups in a dictionary
# -----------------------------------------------------

#organism_name = conf.ortho_name
#ortho_dict = defaultdict(list)
#for line in ortho_lines:
    #line = line.rstrip()
    #split_line = line.split()
    #orthogroup = split_line[0]
    #orthogroup = orthogroup.replace(":", "")
    #genes_in_group = [x for x in split_line if organism_name in x]
    #for gene in genes_in_group:
        #gene = gene.replace(organism_name, '').replace('|', '')
        # print gene
        #ortho_dict[gene] = orthogroup

# -----------------------------------------------------
#
# Build a dictionary of raw count data
#
# -----------------------------------------------------

#raw_read_count_dict = defaultdict(list)

#line1 = raw_count_lines.pop(0)
#line1 = line1.rstrip("\n")
#count_treatment_list = line1.split("\t")
#count_treatment_list = list(filter(None, count_treatment_list))
# print count_treatment_list

#for line in raw_count_lines:
    #line = line.rstrip("\n")
    #split_line = line.split("\t")
    #transcript_id = split_line.pop(0)
    # if not len(count_treatment_list) == len(split_line):
    #     print "error"
    # print len(count_treatment_list)
    # print len(split_line)
    #for i, treatment in enumerate(count_treatment_list):
        # i = i-1
       # raw_read_count = float(split_line[i])
    # for treatment, raw_read_count in zip(count_treatment_list, split_line):
        #dict_key = "_".join([transcript_id, treatment])
        #raw_read_count_dict[dict_key].append(raw_read_count)

# -----------------------------------------------------
#
# Build a dictionary of normalised fpkm data
#
# -----------------------------------------------------

#fpkm_dict = defaultdict(list)

#line1 = fpkm_lines.pop(0)
#line1 = line1.rstrip("\n")
#fpkm_treatment_list = line1.split("\t")
#fpkm_treatment_list = list(filter(None, fpkm_treatment_list))

#for line in fpkm_lines:
    #line = line.rstrip("\n")
    #split_line = line.split("\t")
    #transcript_id = split_line.pop(0)
    #for i, treatment in enumerate(fpkm_treatment_list):
        #fpkm = float(split_line[i])
        #dict_key = "_".join([transcript_id, treatment])
        #fpkm_dict[dict_key].append(fpkm)

# -----------------------------------------------------
#
# Build a dictionary of interproscan annotations
# Annotations first need to be filtered to remove
# redundancy. This is done by first loading anntoations
# into a set.
# -----------------------------------------------------

interpro_set = set([])
interpro_dict = defaultdict(list)

for line in InterPro_lines:
    line = line.rstrip("\n")
    split_line = line.split("\t")
    interpro_columns = []
    index_list = [0, 4, 5, 11, 12]
    for x in index_list:
        if len(split_line) > x:
            interpro_columns.append(split_line[x])
    set_line = ";".join(interpro_columns)
    if set_line not in interpro_set:
        gene_id = interpro_columns[0]
        interpro_feat = ";".join(interpro_columns[1:])
        interpro_dict[gene_id].append(interpro_feat)
    interpro_set.add(set_line)


# -----------------------------------------------------
# Build a dictionary of Swissprot annotations
# -----------------------------------------------------

swissprot_dict = defaultdict(list)

for line in swissprot_lines:
    line = line.rstrip("\n")
    split_line = line.split("\t")
    gene_id = split_line[0]
    swissprot_columns = itemgetter(14, 12, 13)(split_line)

    swissprot_dict[gene_id].extend(swissprot_columns)

#-----------------------------------------------------
# Build a dictionary of TF gene homologs
#-----------------------------------------------------

TF_dict = defaultdict(list)
for line in TF_lines:
    line = line.rstrip("\n")
    split_line = line.split("\t")
    gene_id = split_line[0]
    gene_id = gene_id.replace('.p', '.t')
    TF_function = split_line[2]
    TF_dict[gene_id].append(TF_function)

#-----------------------------------------------------
# Build a dictionary of Secondary Metabolite annotations
#-----------------------------------------------------

#i = 0
antismash_dict = defaultdict(list)
for line in antismash_lines:
    #i += 1
    #cluster = "cluster_" + str(i)
    line = line.rstrip("\n")
    split_line = line.split("\t")
    gene_id = split_line[0]
    secmet_func = split_line[2]
    cluster_genes = split_line[3]
    #.split(";")
    #for gene in cluster_genes:
        #gene = re.sub("_\d$", "", gene)
    antismash_dict[gene_id].extend([secmet_func, cluster_genes])
        # print "\t".join([gene, secmet_func, cluster])

# -----------------------------------------------------
# Step 3
# Iterate through genes in file, identifying if
# they have associated information
# -----------------------------------------------------

# Print header line:
header_line = ['transcript_id']
header_line.extend(['contig', 'start', 'stop', 'strand'])
# header_line.extend(['sigP2', 'sigP4', 'phobius', 'RxLR_motif', 'RxLR_hmm',
#                     'WY_hmm', 'RxLR_total', 'CRN_LFLAK', 'CRN_DWL',
#                     'CRN_total', 'orthogroup'])
#header_line.extend(['sigP4', 'trans_mem', 'secreted', 'effector_total', 'CAZY_total','orthogroup'])
#for treatment in sorted(set(count_treatment_list)):
   # treatment = "raw_count_" + treatment
   # header_line.append(treatment)

#for treatment in sorted(set(fpkm_treatment_list)):
   # treatment = "fpkm_" + treatment
   # header_line.append(treatment)

#for DEG_file in DEG_files:
    #file_name = DEG_file.split('/')[-1]
    #header_line.append("baseMean_" + file_name)
    #header_line.append("LFC_" + file_name)
    #header_line.append("P-val_" + file_name)
header_line.append('prot_seq')
header_line.append('swissprot_org')
header_line.append('swissprot_gene')
header_line.append('swissprot_function')
header_line.append('interproscan')
header_line.append('TFs')
header_line.append('Antismash_type')
header_line.append('Antismash_cluster')
print ("\t".join(header_line))

transcript_lines = []

if conf.gff_format == 'gff3':
    for line in gene_lines:
        line = line.rstrip()
        if line.startswith('#'):
            continue
        split_line = line.split()
        if 'transcript' in split_line[2] or 'mRNA' in split_line[2]:
            transcript_lines.append("\t".join(split_line))

if conf.gff_format == 'gtf':
    prev_id = 'first'
    transcript_line = ''

    for line in gene_lines:
        line = line.rstrip()
        if line.startswith('#'):
            continue
        split_line = line.split("\t")
        if 'CDS' in split_line[2]:
            transcript_id = split_line[8]
            split_col9 = split_line[8].split(';')
            # print split_col9
            transcript_id = "".join([x for x in split_col9 if
                                    'transcript_id' in x])
            # print transcript_id
            transcript_id = transcript_id.replace(' ', '') \
                .replace('transcript_id', '').replace('"', '')
            # print transcript_id
            if transcript_id != prev_id:
                # if prev_id == 'first':
                #     continue
                transcript_lines.append("\t".join(transcript_line))
                transcript_line = split_line
                transcript_line[2] = "mRNA"
                transcript_line[8] = transcript_id
                # print split_line
                # print transcript_line
            elif split_line[6] == '+':
                transcript_line[4] = split_line[4]
            elif split_line[6] == '-':
                transcript_line[3] = split_line[3]
            # print "\t".join([prev_id, transcript_id])
            prev_id = transcript_id
            # print transcript_id
    del transcript_lines[0]

# print "\n".join(transcript_lines)

for line in transcript_lines:
    split_line = line.split("\t")
    useful_cols = [split_line[0],  split_line[3], split_line[4], split_line[6]]
    # Set defaults
    #sigP4 = ''
    #trans_mem = ''
    #effector_total = ''
    #CAZY_total = ''
    #orthogroup = ''
    prot_seq = ''
    swissprot_cols = []
    interpro_col = []
    TFs_cols= ''
    antismash_cols=[]
    # Identify gene id
    if 'ID' in split_line[8]:
        split_col9 = split_line[8].split(';')
        transcript_id = "".join([x for x in split_col9 if 'ID' in x])
        transcript_id = transcript_id.replace('ID=', '')
    else:
        transcript_id = split_line[8]
    #if transcript_id in SigP4_set:
        #sigP4 = 'Yes'
    #if transcript_id in trans_mem_set:
        #trans_mem = ''
    #else:
        #trans_mem = 'Yes'
    #if any([sigP4 == 'Yes']) and all([trans_mem == '']):
        #secreted = 'Yes'
    #else:
        #secreted = ''
    #if transcript_id in effector_total_set:
        #effector_total = 'Yes'
    #if transcript_id in CAZY_total_set:
        #CAZY_total = 'Yes'
    #if ortho_dict[transcript_id]:
        #orthogroup = ortho_dict[transcript_id]
    #DEG_out = []
    #for DEG_file in DEG_files:
        #entryname = "_".join([DEG_file, transcript_id])
        #if DEG_dict[entryname]:
            #DEG_out.append(DEG_dict[entryname][0])
            #DEG_out.append(DEG_dict[entryname][1])
            #DEG_out.append(DEG_dict[entryname][2])
        #else:
            #DEG_out.append('0')
            #DEG_out.append('0')
            #DEG_out.append('0')

    # # Add in read count data:
    #mean_count_cols = []
    #for treatment in sorted(set(count_treatment_list)):
       # dict_key = "_".join([transcript_id, treatment])
       # expression_values = raw_read_count_dict[dict_key]
        # print expression_values
       # mean_count = np.mean(expression_values)
       # mean_count = np.round_(mean_count, decimals=2)
       # mean_count_cols.append(mean_count.astype(str))
    # print mean_count_cols
   # mean_fpkm_cols = []
   # for treatment in sorted(set(fpkm_treatment_list)):
       # dict_key = "_".join([transcript_id, treatment])
        # print dict_key
       # expression_values = fpkm_dict[dict_key]
        # print expression_values
       # mean_fpkm = np.mean(expression_values)
        # print mean_fpkm
       # mean_fpkm = np.round_(mean_fpkm, decimals=0)
       # mean_fpkm_cols.append(mean_fpkm.astype(str))
        # print mean_fpkm_cols

    # # Add in Swissprot info
    if swissprot_dict[transcript_id]:
        swissprot_cols = swissprot_dict[transcript_id]
    else:
        swissprot_cols = ['.', '.', '.']
    # # Add in interproscan info
    if interpro_dict[transcript_id]:
        interpro_col = "|".join(interpro_dict[transcript_id])
    else:
        interpro_col = '.'

    # # Add TFs info
    TFs_cols = []
    if TF_dict[transcript_id]:
        TF_functions = TF_dict[transcript_id]
        TFs_cols.append(";".join(TF_functions))
    else:
        TFs_cols.append("")

    #Add Antismash info
    if antismash_dict[transcript_id]:
        antismash_cols = antismash_dict[transcript_id]
        #secmet_func = antismash_cols[0]
        #cluster = antismash_cols[1]
        #antismash_cols = [cluster, secmet_func]
    else:
        antismash_cols = ["",""]

    prot_seq = "".join(prot_dict[transcript_id])
    # outline = [transcript_id, sigP2, phobius, RxLR_motif, RxLR_hmm,
    #            WY_hmm, CRN_LFLAK, CRN_DWL, orthogroup]
    outline = [transcript_id]
    outline.extend(useful_cols)
    # outline.extend([sigP2, sigP4, phobius, RxLR_motif, RxLR_hmm,
    #                 WY_hmm, RxLR_total, CRN_LFLAK, CRN_DWL, CRN_total,
    #                 orthogroup])
    #outline.extend([sigP4])
    #outline.extend([trans_mem, secreted])
    #outline.extend([effector_total, CAZY_total])
    #outline.append(orthogroup)
    #outline.extend(mean_count_cols)
    #outline.extend(mean_fpkm_cols)
    #outline.extend(DEG_out)
    outline.append(prot_seq)
    outline.extend(swissprot_cols)
    outline.append(interpro_col)
    outline.extend(TFs_cols)
    outline.extend(antismash_cols)
    print "\t".join(outline)
    # print DEG_out
