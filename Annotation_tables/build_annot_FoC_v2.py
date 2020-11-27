#!/usr/bin/python

'''
This program is used to build information on all the genes predicted in
a fusarium genome. These commands take information on location of blast
hits against reference genome & suppliment this information with information
from other files detailing functional annotations and orthology status.
'''

# /home/armita/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/Fo_build_gene_annot_table.py --FoC_genes_gff gene_pred/codingquary/F.oxysporum_fsp_cepae/Fus2_edited_v2/final/final_genes_appended.gff3 --FoC_SigP gene_pred/final_genes_signalp-4.1/F.oxysporum_fsp_cepae/Fus2_edited_v2/Fus2_edited_v2_final_sp.tab --FoC_TM_list gene_pred/trans_mem/F.oxysporum_fsp_cepae/Fus2_edited_v2/Fus2_edited_v2_tmhmm_out.txt --FoC_MIMP_list analysis/mimps/F.oxysporum_fsp_cepae/Fus2_edited_v2/Fus2_edited_v2_genes_in_2kb_mimp.txt --FoC_effectorP analysis/effectorP/F.oxysporum_fsp_cepae/Fus2_edited_v2/F.oxysporum_fsp_cepae_Fus2_edited_v2_EffectorP.txt --FoC_orthogroup analysis/orthology/orthomcl/FoC_path_vs_non_path/FoC_path_vs_non_path_orthogroups.txt




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
ap.add_argument('--blast_csv',required=True,type=str,help='The blast_pipe.sh results file')
ap.add_argument('--FoL_intersected_genes',required=True,type=str,help='A bed file of FoL genes intersecting Blast hits')
ap.add_argument('--genome',required=True,type=str,help='A fasta file of the assembled contigs')
ap.add_argument('--FoC_genes_gff',required=True,type=str,help='A gff file of the genes from the FoC')
# ap.add_argument('--FoC_interescted_reblast',required=True,type=str,help='A bed file of FoC genes intersecting the location of reciprocal blast hits')
ap.add_argument('--FoC_SigP',required=True,type=str,help='A file containing a list of signal-peptide containing genes')
ap.add_argument('--FoC_TM_list',required=True,type=str,help='A file containing a list of transmembrane containing genes')
ap.add_argument('--FoC_MIMP_list',required=True,type=str,help='A file containing a list of genes within 2kb of MIMPs')
ap.add_argument('--FoC_effectorP',required=True,type=str,help='A file containing results of effectorP')
# ap.add_argument('--FoC_expression',required=True,type=str,help='A file showing cufflinks gff features intersected with Braker gene models')
ap.add_argument('--FoC_orthogroup',required=True,type=str,help='A file containing results of orthology analysis')
ap.add_argument('--OrthoMCL_id',required=True,type=str,help='The identifier of this strain as used in the orthology analysis')
ap.add_argument('--OrthoMCL_all',required=True,type=str,nargs='+',help='The identifiers of all strains used in the orthology analysis')
ap.add_argument('--OrthoMCL_path',required=True,type=str,nargs='+',help='The identifiers of pathogenic strains used in the orthology analysis')
ap.add_argument('--OrthoMCL_nonpath',required=True,type=str,nargs='+',help='The identifiers of non-pathogenic strains used in the orthology analysis')
ap.add_argument('--InterPro',required=True,type=str,help='The Interproscan functional annotation .tsv file')
ap.add_argument('--Swissprot',required=True,type=str,help='A parsed table of BLAST results against the Swissprot database. Note - must have been parsed with swissprot_parser.py')
# ap.add_argument('--DEG_Orthogroups',required=True,type=str,help='')


# ap.add_argument('--FoL_LS_genes',required=True,type=str,help='A list of FoL LS genes and the chromosomes they are present on.')


conf = ap.parse_args()


with open(conf.blast_csv) as f:
    blast_csv_lines = f.readlines()

with open(conf.FoL_intersected_genes) as f:
    FoL_intersect_lines = f.readlines()

with open(conf.genome) as f:
    contig_lines = f.readlines()

with open(conf.FoC_genes_gff) as f:
    FoC_genes_lines = f.readlines()
#
# with open(conf.FoC_interescted_reblast) as f:
#     FoC_reblast_lines = f.readlines()
#
with open(conf.FoC_SigP) as f:
    FoC_SigP_lines = f.readlines()

with open(conf.FoC_TM_list) as f:
    FoC_tmhmm_lines = f.readlines()

with open(conf.FoC_MIMP_list) as f:
    FoC_mimp_lines = f.readlines()

with open(conf.FoC_effectorP) as f:
    FoC_effectorP_lines = f.readlines()
#
# with open(conf.FoC_expression) as f:
#     FoC_expression_lines = f.readlines()

with open(conf.FoC_orthogroup) as f:
    FoC_orthogroup_lines = f.readlines()

with open(conf.InterPro) as f:
    InterPro_lines = f.readlines()

with open(conf.Swissprot) as f:
    swissprot_lines = f.readlines()

# with open(conf.DEG_Orthogroups) as f:
#     diff_ortho_lines = f.readlines()

column_list=[]


#-----------------------------------------------------
# Step 2
# Collect information on contig length from the genome
# assembly file. This can be used to determine if a
# gene has 2kb of sequence assembled up and downstream.
# This important for knock out design.
#-----------------------------------------------------

gene_id_set = Set([])
contig_len_dict = defaultdict(list)
contig_id = ""
seq_lines = ""
for line in contig_lines:
    line = line.rstrip()
    if line.startswith(">"):
        last_seq_length = len(seq_lines)
        # print contig_id
        # print last_seq_length
        contig_len_dict[contig_id] = len(seq_lines)
        split_line = line.split(" ")
        contig_id = split_line[0].replace(">", "")
        seq_lines = ""
    else:
        seq_lines += line


#-----------------------------------------------------
# Step 2
# Read Blast csv file and store in a dictionary of hits.
#-----------------------------------------------------
# The csv file may contain columns showing relationship
# between query sequences. The first line of the file
# contains column headers. The number of column headers
# containg 'Grp' can be counted and this number of
# columns skipped when reading lines.
#
# Not all columns will contain blast hits. The number
# of columns in a line must be counted to check if a
# blast hit is present.


blast_id_set = Set([])
i=0
blast_dict = defaultdict(list)
for line in blast_csv_lines:
    line = line.rstrip()
    split_line = line.split("\t")
    if "ID" in split_line[0]:
        for column in split_line:
            if "Grp" in column:
                i+=1
        hit_contig=i+4
        hit_stand=i+9
        hit_start=i+10
        hit_end=i+11
        # extract_list=itemgetter(hit_contig, hit_start, hit_end)(split_line)
        continue
    blast_id=split_line[0]
    blast_id_set.add(blast_id)
    # column_list = ["no hit", ".", ".", "."]
    column_list = ["", "", "", ""]
    if len(split_line) > hit_contig:
        column_list=itemgetter(hit_contig, hit_start, hit_end, hit_stand)(split_line)
    for column in column_list:
        blast_dict[blast_id].append(column)

#-----------------------------------------------------
# Step 3
# Build a dictionary of intersected genes
# The final column of the input file containing
# gff annotation information was split into two columns
#-----------------------------------------------------

intersect_dict = defaultdict(list)
for line in FoL_intersect_lines:
    line = line.rstrip()
    split_line = line.split("\t")
    if "gene" in split_line[11]:
        blast_id = split_line[8].strip('";')
        blast_id = blast_id.replace('ID=', '').replace('_BlastHit_1', '')
        column_list = itemgetter(12, 13, 15)(split_line)
        for column in column_list:
            intersect_dict[blast_id].append(column)

        feature_info_list = "".join(split_line[17]).split(";")
        FoL_gene_ID = ""
        FoL_feature_desc = ""
        for feature in feature_info_list:
            if "gene_id=" in feature:
                FoL_gene_ID = re.sub(r"^.*=", '', feature)
            if "description=" in feature:
                FoL_feature_desc = re.sub(r"^.*=", '', feature)
        intersect_dict[blast_id].extend([FoL_gene_ID, FoL_feature_desc])

#-----------------------------------------------------
# Step 4
# Append co-ordinates from the FoC genome, showing
# source of original Blast queries.
# Also identify whether there is 2kb sequence data up
# and downstream of the gene allowing design of
# knockouts
#-----------------------------------------------------

gene_id_set = Set([])
FoC_genes_dict = defaultdict(list)
for line in FoC_genes_lines:
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
# Step 5
# Build a dictionary of FoC genes intersecting the hit
# location of reciprocal blast queries.
#-----------------------------------------------------

# FoC_reblast_dict = defaultdict(list)
# for line in FoC_reblast_lines:
#     line = line.rstrip()
#     split_line = line.split("\t")
#     if "transcript" in split_line[11]:
#         query_id = split_line[8].strip('";')
#         query_id = query_id.replace('ID=', '').replace('_BlastHit_1', '').replace('extracted_hit_of_', '')
#         column_list = []
#         intersect_id = split_line[17]
#         if query_id in intersect_id:
#             FoC_reblast_dict[query_id]=[intersect_id, "match"]
#         elif "match" in FoC_reblast_dict[query_id]:
#             pass
#         else:
#             FoC_reblast_dict[query_id]=[intersect_id, ""]

#-----------------------------------------------------
# Steps 6 & 7
# Build a dictionary of FoC genes containing Signal
# peptides
#
# Extend this dictionary with information on FoC genes
# containing transmembrane domains.
#
# If a protein contains a signal peptide and does not
# contain a transmembrane domain, it is treated as a
# putative secreted protein.
#-----------------------------------------------------

FoC_secreted_dict = defaultdict(list)
for line in FoC_SigP_lines:
    line = line.rstrip("\n")
    split_line = line.split("\t")
    gene_id = split_line[0].replace('>', '')
    FoC_secreted_dict[gene_id] = split_line[1:]

for line in FoC_tmhmm_lines:
    line = line.rstrip("\n")
    split_line = line.split("\t")
    gene_id = split_line[0]
    num_helices = split_line[4].replace('PredHel=', '')
    if num_helices == '0':
        FoC_secreted_dict[gene_id].extend(["NO", num_helices])
        if FoC_secreted_dict[gene_id][0] == "YES":
            FoC_secreted_dict[gene_id].append("Yes")
        else:
            FoC_secreted_dict[gene_id].append("")
    else:
        FoC_secreted_dict[gene_id].extend(["YES", num_helices, ""])

#-----------------------------------------------------
# Step 8
# Build a dictionary of FoC genes within 2Kb of MIMPs
#-----------------------------------------------------

FoC_mimp_set =  Set([])

for line in FoC_mimp_lines:
    line = line.rstrip("\n")
    split_line = line.split("\t")
    gene_id = split_line[0]
    FoC_mimp_set.add(gene_id)

#-----------------------------------------------------
# Step 9
# Build a dictionary of FoC genes that have been
# predicted as effectors by effectorP
#-----------------------------------------------------

FoC_effectorP_dict = defaultdict(list)
First = True
for line in FoC_effectorP_lines:
    if First == True:
        First = False
        continue
    line = line.rstrip("\n")
    split_line = line.split("\t")
    gene_id = split_line[0]
    if gene_id in gene_id_set:
        split_line[1] = split_line[1].replace('Effector', 'Yes').replace("Non-effector", "")
        column_list = itemgetter(1, 2)(split_line)
    else:
        column_list = ["", ""]
    FoC_effectorP_dict[gene_id].extend(column_list)


#-----------------------------------------------------
# Step 10
# Build a dictionary of Fus2 genes that intersect with
# cufflinks assembled transcripts. Extract information
# on FPKM from cufflinks transcripts intersecting Braker
# gene models
#
# In some cases a gene may intersect multiple cufflinks
# Transcripts. In this case take the stats from the
# one with the longest overlap
#
#-----------------------------------------------------

# FoC_expression_dict = defaultdict(list)
# for line in FoC_expression_lines:
#     line = line.rstrip("\n")
#     split_line = line.split("\t")
#     if "transcript" in split_line[2] and "transcript" in split_line[11]:
#         gene_id = split_line[17]
#         overlap = split_line[18]
#         cufflinks_features = split_line[8].split(";")
#         cufflinks_id = str(cufflinks_features[0])
#         cufflinks_id = cufflinks_id[9:-1]
#         FPKM = str(cufflinks_features[2])
#         FPKM = FPKM[7:-1]
#         cov = str(cufflinks_features[6])
#         cov = cov[6:-1]
#         column_list = [cufflinks_id, FPKM, cov, overlap]
#         if not FoC_expression_dict[gene_id]:
#             FoC_expression_dict[gene_id].extend(column_list)
#         else:
#             prev_overlap_lgth = FoC_expression_dict[gene_id][3]
#             if overlap > prev_overlap_lgth:
#                 FoC_expression_dict[gene_id].extend(column_list)


#-----------------------------------------------------
# Step 11
# Build a dictionary of orthogroups
#-----------------------------------------------------

strain_id = conf.OrthoMCL_id
all_isolates = conf.OrthoMCL_all
path_isolates = conf.OrthoMCL_path
non_path_isolates = conf.OrthoMCL_nonpath

strain_id = strain_id + "|"

FoC_orthogroup_dict = defaultdict(list)
orthogroup_content_dict = defaultdict(list)
# path_isolates = ["Fus2", "125", "A23"]
# non_path_isolates = ["A28", "D2", "PG"]
# all_isolates = path_isolates + non_path_isolates

for line in FoC_orthogroup_lines:
    line = line.rstrip("\n")
    split_line = line.split(" ")
    orthogroup_id = split_line[0].replace(":", "")
    orthogroup_contents = []
    orthogroup_content_dict.clear()
    for isolate in all_isolates:
        num_genes = line.count((isolate + "|"))
        orthogroup_contents.append(str(isolate) + "(" + str(num_genes) + ")")
        content_str = ":".join(orthogroup_contents)
        orthogroup_content_dict[isolate] = num_genes

    path_numbers = []
    for isolate in path_isolates:
        path_numbers.append(orthogroup_content_dict[isolate])
    max_path = max(path_numbers)
    min_path = min(path_numbers)
    non_path_numbers = []
    for isolate in non_path_isolates:
        non_path_numbers.append(orthogroup_content_dict[isolate])
    max_non_path = max(non_path_numbers)
    min_non_path = min(non_path_numbers)
    if min_path > max_non_path:
        expansion_status = "pathogen_expanded"
    elif min_non_path > max_path:
        expansion_status = "non-pathogen_expanded"
    else:
        expansion_status = ""

    for gene_id in split_line[1:]:
        if strain_id in gene_id:
            gene_id = gene_id.replace(strain_id, "")
            # FoC_orthogroup_dict[gene_id] = [orthogroup_id]

            if all(x in line for x in all_isolates):
                FoC_orthogroup_dict[gene_id].extend([orthogroup_id, "all_isolates", content_str, expansion_status])
            elif all(x not in line for x in non_path_isolates) and all(x in line for x in path_isolates):
                FoC_orthogroup_dict[gene_id].extend([orthogroup_id, "path_isolates_all", content_str, expansion_status])
            elif all(x in line for x in non_path_isolates) and all(x not in line for x in path_isolates):
                FoC_orthogroup_dict[gene_id].extend([orthogroup_id, "non_path_isolates_all", content_str, expansion_status])
            elif any(x in line for x in all_isolates):
                FoC_orthogroup_dict[gene_id].extend([orthogroup_id, "some_isolates", content_str, expansion_status])

            # if all(x in line for x in non_path_isolates):

            # print gene_id



#-----------------------------------------------------
# Step 12
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
    # print split_line
    # useful_columns = (split_line[col] for col in [0, 4, 6, 12, 13] if split_line[col])
    # useful_columns = (split_line[col], col in 0, 4, 6, 12, 13 if len(split_lin) >= col)
    # useful_columns.extend() (x in line for x in 0, 4, 6, 12, 13)
    interpro_columns = []
    index_list = [0, 4, 5, 11, 12]
    for x in index_list:
        # print x
        if len(split_line) > x:
            interpro_columns.append(split_line[x])
            # print split_line[x]


    # useful_columns = itemgetter(0, 4, 6, 12, 13)(split_line)
    set_line = ";".join(interpro_columns)
    if set_line not in interpro_set:
        gene_id = interpro_columns[0]
        interpro_feat = ";".join(interpro_columns[1:])
        interpro_dict[gene_id].append(interpro_feat)
    interpro_set.add(set_line)


#-----------------------------------------------------
# Step 13
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
# Step 13
# Build a dictionary of differentially expressed
# orthogroups according in comparison of Fus2 & Fo47
# reads aligned to Fus2 genes.
#-----------------------------------------------------

# diff_ortho_dict = defaultdict(list)
#
# for line in diff_ortho_lines:
#     line = line.rstrip("\n")
#     split_line = line.split("\t")
#     if len(split_line) > 5:
#         orthogroup_id = split_line[5]
#         fus2_id = "Fus2|" + str(split_line[0])
#         base_mean = split_line[1]
#         fold_change = split_line[2]
#         significance = split_line[4]
#         expression_info = " ".join([fus2_id, base_mean,fold_change,significance])
#         diff_ortho_dict[orthogroup_id].append(expression_info)




#-----------------------------------------------------
# Step 12
# Print final table of information on query, blast
# results and genes intersecting blast results
#-----------------------------------------------------

print ("\t".join([
"query_id", "FoC_contig", "FoC_gene_start", "FoC_gene_end", "FoC_gene_strand",
"2kb_flank",
"SigP", "P-value", "cleavage_site",
"TransMem_protein", "No._helices",
"Secreted",
"MIMP_in_2Kb",
"effectorP", "P-value",
"orthogroup", "representation", "genes_per_isolate", "expansion_status",
# "DEG_in_orthogroup", "Fus2_expression_info",
"hit_FoL_contig", "hit_start", "hit_end", "hit_strand",
"FoL_gene_start", "FoL_gene_end", "FoL_strand", "FoL_gene_ID", "FoL_gene_description",
"Swissprot_organism", "Swissprot_hit", "Swissprot_function",
"Interpro_annotations"
]))
# print ("\t".join(["query_id", "FoC_contig", "FoC_gene_start", "FoC_gene_end", "FoC_gene_strand", "cufflinks_id","FPKM", "cov", "SigP", "P-value", "cleavage_site", "TransMem_protein", "No._helices", "Secreted", "MIMP_in_2Kb", "effectorP", "P-value", "hit_FoL_contig", "hit_start", "hit_end", "hit_strand", "reblast_hit", "reblast match", "FoL_gene_start", "FoL_gene_end", "FoL_strand", "FoL_gene_ID", "FoL_gene_description"]))

# for blast_id in blast_id_set:
for gene_id in gene_id_set:
    useful_columns=[gene_id]
    # useful_columns=[blast_id]
    useful_columns.extend(FoC_genes_dict[gene_id])
    # if FoC_expression_dict[gene_id]:
    #     useful_columns.extend(FoC_expression_dict[gene_id][0:3])
    # else:
    #     useful_columns.extend(["", "", ""])
    useful_columns.extend(FoC_secreted_dict[gene_id])
    mimp_col=""
    if gene_id in FoC_mimp_set:
        mimp_col="Yes"
    useful_columns.append(mimp_col)
    useful_columns.extend(FoC_effectorP_dict[gene_id])

    if FoC_orthogroup_dict[gene_id]:
        useful_columns.extend(FoC_orthogroup_dict[gene_id])
        orthogroup_id = FoC_orthogroup_dict[gene_id][0]
        # diff_expr_ortho_info = []
        # # print orthogroup_id
        # # print diff_ortho_dict[orthogroup_id]
        # diff_expr_ortho_info.append(";".join(diff_ortho_dict[orthogroup_id]))
        # if any("P<0.05" in x for x in diff_expr_ortho_info):
        #     significance = "P<0.05"
        #     # print "yes"
        # else:
        #     significance = ""
        #     # print "no"
        # useful_columns.append(significance)
        # # print diff_expr_ortho_info
        # useful_columns.extend(diff_expr_ortho_info)

    else:
        useful_columns.extend(["", "singleton", "", ""])
        # useful_columns.extend(["", ""])

    useful_columns.extend(blast_dict[gene_id])
    # if FoC_reblast_dict[blast_id]:
    #     useful_columns.extend(FoC_reblast_dict[blast_id])
    # else:
    #     useful_columns.extend(["", ""])
    if intersect_dict[gene_id]:
        intersect_cols = intersect_dict[gene_id][0:5]
        useful_columns.extend(intersect_cols)
    else:
        useful_columns.extend(["", "", "", "", ""])



    if swissprot_dict[gene_id]:
        useful_columns.extend(swissprot_dict[gene_id])
    else:
        useful_columns.extend(["","",""])

    if interpro_dict[gene_id]:
        interpro_col = "|".join(interpro_dict[gene_id])
        useful_columns.append(interpro_col)
    else:
        useful_columns.append("")


    print ("\t".join(useful_columns))
