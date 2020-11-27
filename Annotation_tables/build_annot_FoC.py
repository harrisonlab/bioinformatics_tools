#!/usr/bin/python

'''
This program is used to build information on all the genes predicted in
an annotated genome. These commands take information on location of genes
& suppliment this information with information on interproscan domains
and swissprot annotations.
'''


#-----------------------------------------------------
# Step 1
# Import variables & load input files
#-----------------------------------------------------

import sys
import argparse
import re
from sets import Set
from collections import defaultdict
from operator import itemgetter

ap = argparse.ArgumentParser()

ap.add_argument('--genome',required=True,type=str,help='A fasta file of the assembled contigs')
ap.add_argument('--genes_gff',required=True,type=str,help='A gff file of the genes')
ap.add_argument('--genes_renamed',required=True,type=str,help='A gene name conversion table')
ap.add_argument('--SigP',required=True,type=str,help='A file containing a list of signal-peptide containing genes')
ap.add_argument('--TM_list',required=True,type=str,help='txt file of headers from gene testing positive for tranmembrane proteins by TMHMM')
ap.add_argument('--GPI_list',required=True,type=str,help='txt file of headers from gene testing positive for GPI anchors as identified by GPI-SOM')
ap.add_argument('--MIMP_list',required=True,type=str,help='A file containing a list of genes within 2kb of MIMPs')
ap.add_argument('--EffP_list',required=True,type=str,help='A file containing results of effectorP')
ap.add_argument('--CAZY_list',required=True,type=str,help='A file containing results of CAZY')

ap.add_argument('--InterPro',required=True,type=str,help='The Interproscan functional annotation .tsv file')
ap.add_argument('--Swissprot',required=True,type=str,help='A parsed table of BLAST results against the Swissprot database. Note - must have been parsed with swissprot_parser.py')
ap.add_argument('--Antismash',required=True,type=str,help='Output of Antismash parsed into a tsv file of gene names, contig, secmet function and cluster ID')
# ap.add_argument('--Smurf',required=True,type=str,help='Output of Smurf parsed into a tsv file of gene names, contig, secmet function and cluster ID')
ap.add_argument('--TFs',required=True,type=str,help='Tab seperated of putative transcription factors and their domains as identified by interpro2TFs.py')
ap.add_argument('--orthogroups', required=True,type=str,help='A file containing results of orthology analysis')
ap.add_argument('--strain_id',required=True,type=str,help='The identifier of this strain as used in the orthology analysis')
ap.add_argument('--OrthoMCL_all',required=True,type=str,nargs='+',help='The identifiers of all strains used in the orthology analysis')

conf = ap.parse_args()

with open(conf.genome) as f:
    contig_lines = f.readlines()
with open(conf.genes_gff) as f:
    unordered_gff_lines = f.readlines()
with open(conf.genes_renamed) as f:
    genes_renamed_lines = f.readlines()
with open(conf.InterPro) as f:
    InterPro_lines = f.readlines()
with open(conf.Swissprot) as f:
    swissprot_lines = f.readlines()
with open(conf.Antismash) as f:
    antismash_lines = f.readlines()
# with open(conf.Smurf) as f:
#     smurf_lines = f.readlines()
with open(conf.TFs) as f:
    TF_lines = f.readlines()
with open(conf.orthogroups) as f:
    orthogroup_lines = f.readlines()

with open(conf.SigP) as f:
    sigP_lines = f.readlines()
with open(conf.TM_list) as f:
    TM_lines = f.readlines()
with open(conf.GPI_list) as f:
    gpi_lines = f.readlines()
with open(conf.MIMP_list) as f:
    mimp_lines = f.readlines()
with open(conf.EffP_list) as f:
    effectorP_lines = f.readlines()
with open(conf.CAZY_list) as f:
    cazy_lines = f.readlines()

column_list=[]

#-----------------------------------------------------
# Step 2
# Collect information on contig length from the genome
# assembly file. This can be used to determine if a
# gene has 2kb of sequence assembled up and downstream.
# This important for knock out design.
#-----------------------------------------------------

contig_len_dict = defaultdict(list)
contig_id = ""
seq_lines = ""
for line in contig_lines:
    line = line.rstrip()
    if line.startswith(">"):
        last_seq_length = len(seq_lines)
        contig_len_dict[contig_id] = len(seq_lines)
        split_line = line.split(" ")
        contig_id = split_line[0].replace(">", "")
        seq_lines = ""
    else:
        seq_lines += line

#-----------------------------------------------------
# Step 2.5
# Order gff features by input contigs and gene start
#-----------------------------------------------------

contig_list = []
gene_start_dict = defaultdict(list)
features_dict = defaultdict(list)
ordered_gff_lines = []

for line in unordered_gff_lines:
    line = line.strip("\n")
    if line.startswith("#"):
        ordered_gff_lines.append(line)
        continue
    split_line = line.split()
    if split_line[2] == 'gene':
        contig = split_line[0]
        gene_start = split_line[3]
        if contig not in contig_list:
            contig_list.append(contig)
        gene_start_dict[contig].append(int(gene_start))
        key = "_".join([contig, gene_start])
    features_dict[key].append(line)

gff_lines = []
for contig in sorted(contig_list, key = lambda x: (int(x.split('_')[1]))):
    gene_start_list = gene_start_dict[contig]
    for gene_start in sorted(gene_start_list):
        key = "_".join([contig, str(gene_start)])
        ordered_gff_lines.extend(features_dict[key])

#-----------------------------------------------------
# Step 3
# Append co-ordinates from the FoC gene gff, showing
# gene locations.
# Also identify whether there is 2kb sequence data up
# and downstream of the gene allowing design of
# knockouts
#-----------------------------------------------------

gene_id_set = Set([])
gene_id_list = []
FoC_genes_dict = defaultdict(list)
for line in ordered_gff_lines:
    if "gff-version" in line:
        continue
    if line.startswith('#'):
        continue
    line = line.rstrip()
    split_line = line.split("\t")
    if 'mRNA' in split_line[2]:
        gene_features = split_line[8].split(';')
        gene_id = gene_features[0]
        gene_id = gene_id.replace('ID=', '')
        column_list = ["", "", "", ""]
        if gene_id not in gene_id_set:
            gene_id_list.append(gene_id)
            gene_id_set.add(gene_id)
        column_list=itemgetter(0, 3, 4, 6)(split_line)
        for column in column_list:
            FoC_genes_dict[gene_id].append(column)

        contig_id = column_list[0]
        feature_start=int(column_list[1])
        feature_end=int(column_list[2])
        contig_length = contig_len_dict[contig_id]
        if (feature_start - 2000) > 0 and (feature_end +2000) < contig_length:
            FoC_genes_dict[gene_id].append("Flank")
        else:
            FoC_genes_dict[gene_id].append("")

#-----------------------------------------------------
# Step 4
# Build a dictionary of interproscan annotations
# Annotations first need to be filtered to remove
# redundancy. This is done by first loading anntoations
# into a set.
#-----------------------------------------------------

interpro_set =  Set([])
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


#-----------------------------------------------------
# Step 5
# Build a dictionary of Swissprot annotations
#-----------------------------------------------------

swissprot_dict = defaultdict(list)

for line in swissprot_lines:
    line = line.rstrip("\n")
    split_line = line.split("\t")
    gene_id = split_line[0]
    swissprot_columns = itemgetter(14, 12, 13)(split_line)

    swissprot_dict[gene_id].extend(swissprot_columns)


#-----------------------------------------------------
# Step 5
# Build a dictionary of Secondary Metabolite annotations
#-----------------------------------------------------
#
i = 0
antismash_dict = defaultdict(list)
for line in antismash_lines:
    i += 1
    cluster = "cluster_" + str(i)
    line = line.rstrip("\n")
    split_line = line.split("\t")
    secmet_func = split_line[2]
    cluster_genes = split_line[3].split(";")
    for gene in cluster_genes:
        gene = re.sub("_\d$", "", gene)
        antismash_dict[gene].extend([secmet_func, cluster])
        # print "\t".join([gene, secmet_func, cluster])
#
# smurf_dict = defaultdict(list)
# for line in smurf_lines:
#     line = line.rstrip("\n")
#     split_line = line.split("\t")
#     gene_id = split_line[0]
#     secmet_func = split_line[2]
#     cluster = "SM_" + split_line[3]
#     smurf_dict[gene_id].extend([secmet_func, cluster])


#-----------------------------------------------------
# Step 6
# Build a dictionary of vitamin TF genes
#-----------------------------------------------------

TF_dict = defaultdict(list)
for line in TF_lines:
    line = line.rstrip("\n")
    split_line = line.split("\t")
    gene_id = split_line[0]
    TF_function = split_line[2]
    TF_dict[gene_id].append(TF_function)


#-----------------------------------------------------
# Step 7
# Build a dictionary of orthogroups
#-----------------------------------------------------

strain_id = conf.strain_id + "|"
all_isolates = conf.OrthoMCL_all

orthogroup_dict = defaultdict(list)
orthogroup_content_dict = defaultdict(list)

for line in orthogroup_lines:
    line = line.rstrip("\n")
    split_line = line.split(" ")
    orthogroup_id = split_line[0].replace(":", "")
    orthogroup_contents = []
    orthogroup_content_dict.clear()
    for isolate in all_isolates:
        num_genes = line.count((isolate + "|"))
        orthogroup_contents.append(str(isolate) + "(" + str(num_genes) + ")")
        content_counts = ":".join(orthogroup_contents)
        orthogroup_content_dict[isolate] = num_genes
    for gene_id in split_line[1:]:
        content_str = ",".join(split_line[1:])
        if strain_id in gene_id:
            gene_id = gene_id.replace(strain_id, '')
            orthogroup_dict[gene_id].extend([orthogroup_id, content_counts, content_str])


#-----------------------------------------------------
# Step 7
# Build a dictionary of orthogroups
#-----------------------------------------------------

rename_dict = defaultdict(list)
for line in genes_renamed_lines:
    line = line.rstrip("\n")
    split_line = line.split()
    old_name = split_line[0].replace("BFJ63_", "")
    new_name = split_line[1].replace("BFJ63_vA", "")
    rename_dict[old_name] = new_name

#-----------------------------------------------------
# Load signalP files into a set
#-----------------------------------------------------

SigP_set = Set()
for line in sigP_lines:
    line = line.rstrip()
    if line.startswith('>'):
        split_line = line.split()
        header = split_line[0].replace('>', '')
        SigP_set.add(header)

#-----------------------------------------------------
# Load TMHMM headers into a set
#-----------------------------------------------------

TM_set = Set()
for line in TM_lines:
    line = line.rstrip()
    header = line.split()[0]
    TM_set.add(header)

#-----------------------------------------------------
# Load GPI-anchored proteins into a set
#-----------------------------------------------------

gpi_set = Set()
for line in gpi_lines:
    line = line.rstrip()
    if line.startswith('>'):
        split_line = line.split()
        header = split_line[0].replace('>', '')
        gpi_set.add(header)

#-----------------------------------------------------
# Load MIMP proteins into a set
#-----------------------------------------------------

mimp_set = Set()
for line in mimp_lines:
    header = line.rstrip()
    mimp_set.add(header)

#-----------------------------------------------------
# Load EffectorP proteins into a set
#-----------------------------------------------------

effectorP_set = Set()
for line in effectorP_lines:
    header = line.rstrip()
    effectorP_set.add(header)

#-----------------------------------------------------
# Load CAZY proteins into a set
#-----------------------------------------------------

cazy_set = Set()
for line in cazy_lines:
    header = line.rstrip()
    cazy_set.add(header)

#-----------------------------------------------------
# Step 6
# Print final table of information on query, blast
# results and genes intersecting blast results
#-----------------------------------------------------

print ("\t".join([
"ncbi_gene_id", "internal_id", "contig", "gene_start", "gene_end", "gene_strand",
"2kb_flank",
"Within_2Kb_of_MIMP",
"Signal_peptide", "Trans-membrane_domain", "GPI_anchor_site", "Secreted",
"EffectorP", "CAZY",
"Cluster_ID", "SecMet_function",
"TF",
"Swissprot_organism", "Swissprot_hit", "Swissprot_function",
"Interpro_annotations",
"orthogroup", "genes_per_isolate", "orthogroup_contents",
]))

in_cluster = False
cluster_num = 0
for gene_id in gene_id_list:
    base_id = gene_id.split('.')[0]
    new_id = gene_id.replace(base_id, "".join(rename_dict[base_id]))
    useful_columns=[new_id, gene_id]
    useful_columns.extend(FoC_genes_dict[gene_id])

    if base_id in mimp_set:
        useful_columns.append("MIMP")
    else:
        useful_columns.append("")

    if gene_id in SigP_set:
        useful_columns.append("SigP")
    else:
        useful_columns.append("")

    if gene_id in TM_set and gene_id in SigP_set:
        useful_columns.append("TM")
    else:
        useful_columns.append("")

    if gene_id in gpi_set and gene_id in SigP_set:
        useful_columns.append("GPI")
    else:
        useful_columns.append("")

    if gene_id in SigP_set and gene_id not in TM_set and gene_id not in gpi_set:
        useful_columns.append("Secreted")
    else:
        useful_columns.append("")

    if gene_id in effectorP_set:
        useful_columns.append("EffP")
    else:
        useful_columns.append("")

    if gene_id in cazy_set:
        useful_columns.append("CAZY")
    else:
        useful_columns.append("")

    if antismash_dict[base_id]:
        antismash_cols = antismash_dict[base_id]
        secmet_func = antismash_cols[0]
        cluster = antismash_cols[1]
        useful_columns.extend([cluster, secmet_func])
    else:
        useful_columns.extend(["",""])

    if TF_dict[gene_id]:
        TF_functions = TF_dict[gene_id]
        useful_columns.append(";".join(TF_functions))
    else:
        useful_columns.append("")

    if swissprot_dict[gene_id]:
        useful_columns.extend(swissprot_dict[gene_id])
    else:
        useful_columns.extend(["","",""])

    if interpro_dict[gene_id]:
        interpro_col = "|".join(interpro_dict[gene_id])
        useful_columns.append(interpro_col)
    else:
        useful_columns.append("")

    if orthogroup_dict[gene_id]:
        useful_columns.extend(orthogroup_dict[gene_id])
    else:
        useful_columns.extend(["", ""])


    print ("\t".join(useful_columns))
