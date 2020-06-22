#!/home/gomeza/miniconda3/envs/antismash_py27/bin/python2.7

'''
This tool extracts and parses regions predicted to contain secondary metabolite
synthesis genes by antismash VERSION 5.0. These regions are output as gff features. Also
predicted antismash ORFs within these regions are intersected with gene models
and output as additional Gff files.
'''

#-----------------------------------------------------
# Step 1
# Import variables & load input files
#-----------------------------------------------------

import sys,argparse
import re
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from sets import Set
from collections import defaultdict
from operator import itemgetter

ap = argparse.ArgumentParser()
ap.add_argument('--inp_antismash',required=True,type=str,help='.gbk file of antismash annotations')
# ap.add_argument('--gff',required=True,type=str,help='input file of gene models in gff format')
ap.add_argument('--out_prefix',required=True,type=str,help='output directory and file prefix for output gff files')


conf = ap.parse_args()
out = conf.out_prefix

#-----------------------------------------------------
# Step 2
# Use the start and stop location of each feature
# labelled as a "cluster" to ubild a gff annotation
# and print to stdout
#
#-----------------------------------------------------

cf_outlines = []
secmet_outlines = []

for gb_record in SeqIO.parse(open(conf.inp_antismash,"r"), "genbank"):
    contig_id = gb_record.id
    for feature in gb_record.features:
        if feature.type == "cand_cluster":
            feature_label = "antismash_cluster"
            feature_type = "_".join(feature.qualifiers["product"])
            feature_start = feature.location.start +1
            col_6 = "."
            feature_stop = feature.location.end +1
            col_8 = "."
            strand = str(feature.location.strand)
            strand = strand.replace("1","+").replace("0","-")
            feature_id = feature.qualifiers["candidate_cluster_number"][0]
            feature_id = feature_id.replace(' number: ', '_')
            #feature_note =  "".join(feature.qualifiers["detection_rules"][1:])
            feature_kind = "_".join(feature.qualifiers["kind"])
            col_9 = ";".join(["ID=" + contig_id + "_Cluster_" + feature_id, "Kind=" + feature_kind])#, "Notes=" + feature_note])
            outline = "\t".join([contig_id, feature_label, feature_type, str(feature_start), str(feature_stop), col_6, strand, col_8, col_9])
            if feature_type.startswith("cf_"):
                cf_outlines.append(outline)
            else:
                secmet_outlines.append(outline)

    out_secmet_gff = open(out + "_secmet_clusters.gff","w")
    out_secmet_gff.write("\n".join(str(x) for x in secmet_outlines) + "\n")
    out_secmet_gff.close()

    out_cf_gff = open(out + "_clusterfinder_clusters.gff","w")
    out_cf_gff.write("\n".join(str(x) for x in cf_outlines) + "\n")
    out_cf_gff.close()
