#!/usr/bin/python

'''
make a gffutils db from a gff file
'''

import sys,argparse
import gffutils
i = 0

ap = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
ap.add_argument('--inp',required=True,type=str,help='input gff file')
ap.add_argument('--db',required=True,type=str,help='output db file')
conf = ap.parse_args() #sys.argv

id_dict = {}
def transform_func(x):
   # Ensures that all features have a unique ID
    if not 'ID' in x.attributes:
        new_id = "".join(x.attributes['Parent']) + "." + str(x.featuretype)
        if id_dict.has_key (new_id):
            id_dict[new_id] += 1
        else:
            id_dict[new_id] = 1
        x.attributes['ID'] = new_id + str(id_dict[new_id])
    return x


db1 = gffutils.create_db(
    conf.inp,
    dbfn=conf.db,
    force=True,
    keep_order=True,
    merge_strategy="merge",
    transform=transform_func,
    sort_attribute_values=False,
    )
