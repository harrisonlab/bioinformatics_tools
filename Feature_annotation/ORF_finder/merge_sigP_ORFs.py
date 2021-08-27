#!/usr/bin/python

'''
This tool will open a gffutils database and identify gene features that overlap
one another. The hmm score for each of these features will be identified. The
feature with the greatest hmm score will be retained and a new database built.
'''

import sys,argparse
import gffutils
import re
from collections import defaultdict
from itertools import chain

#######################################
#            Import variables         #
#                                     #
#                                     #
#######################################

ap = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.ArgumentDefaultsHelpFormatter)

ap.add_argument('--inp',required=True,type=str,help='databases to input')
ap.add_argument('--id',required=True,type=str,help='The ID for newly merged features will start with this string')
# ap.add_argument('--source',required=False,type=str,default='merge_db_features.py',help='A string describing the source of newly created features')
ap.add_argument('--gff',required=False,action='store_true',help='Will print the output in .gff3 format to stdout if set')
ap.add_argument('--out',required=True,type=str,help='The name of the output database.')

conf = ap.parse_args() #sys.argv

f = conf.inp
o = open(conf.out, 'w')
s = conf.gff

#######################################
#          Open the gffutils          #
#               database              #
#                                     #
#######################################

db = gffutils.FeatureDB(f)

#######################################
#        Create a function to         #
#      Serarch for overlapping        #
#             features                #
#######################################

hmm_dict = defaultdict(list)
def merge_func(func_db):
    genes = func_db.features_of_type('transcript')
    d = {}
    i = 0
    out_lines = []
    # For each gene:
    for gene in genes:
        ID = "".join(gene.attributes['ID'])
        # See if that gene has been previously merged
        if not ID in d.keys():
            # Otherwise, look for overlapping features on the same strand
            strand = gene.strand
            overlaps = func_db.region(region = gene, strand=strand, featuretype = 'transcript')
            list = []
            hmm_dict = {}
            # For each overlapping feature, add this to a dictionary of
            # all the merged features.
            for feature in overlaps:
                x = "".join(feature.attributes['ID'])
                d[x] = ""
                # Create a list of the current overlapping features
                # A list may have already been created in the notes section
                # From the first itteration of the function, if this is the
                # case, then take the list of IDs from the notes attribute
                # and add them to the new list of IDs.
                dont_append = False
                prev_notes = []
                if 'Note' in feature.attributes:
                    note_list = []
                    note_list = feature.attributes['Note']
                    for note in note_list:
                        # Search through notes and extract hmm scores.
                        # Identify the gene with the greatest hmm score.
                        if "HMM_score:" in note:
                            # print note
                            # m = re.match(r'HMM_score:(.*)\s', str(note))
                            m = re.findall(r'HMM_score:(\d*\.\d*)', str(note))
                            overlap_hmm_score = max(m)
                            # overlap_hmm_score = m.group(1)
                            # overlap_hmm_score = re.match(r'HMM_score:(.*)', str(note))
                            # print overlap_hmm_score
                            # hmm_dict[overlap_hmm_score].append(overlap_hmm_score)
                            hmm_dict[overlap_hmm_score] = x
            if dont_append == False:
                best_hmm = max(hmm_dict.keys())
                best_feature_id = hmm_dict[best_hmm]
                parent = "".join( db[best_feature_id].attributes['Parent'] )
                out_lines.append( db[parent] )
                out_lines.append( db[best_feature_id] )
    return(out_lines)


#######################################
#        Create a database of merged  #
#                features             #
#                                     #
#######################################

# Perfrom a first itteration of merging features
tmp_lines =  merge_func(db)
db2 = gffutils.create_db(
	tmp_lines,
	from_string=True,
	dbfn=':memory:',
	force=True,
	keep_order=False,
	sort_attribute_values='merge',
	merge_strategy='merge',
	id_spec=['ID']
	)

# Perform a second itteration of merging
out_lines =  merge_func(db2)
merged_db = gffutils.create_db(
	out_lines,
	from_string=True,
	dbfn=conf.out,
	force=True,
	keep_order=False,
	sort_attribute_values='merge',
	merge_strategy='merge',
	id_spec=['ID']
	)

# If the switch was set on stdin, print a gff file to stdout
if s == True:
    print("##gff-version 3")
    for line in out_lines:
        print(str(line))
