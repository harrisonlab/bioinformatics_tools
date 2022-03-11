#!/mnt/shared/scratch/agomez/apps/conda/bin/python

'''
This tool exracts fasta accessions from a set of headers passed in a text file.
'''

import sys,argparse,re
#from sets import Set
my_set = set()

ap = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
ap.add_argument('--fasta',required=True,type=str,help='Fasta file to extract from')
ap.add_argument('--headers',required=True,type=str,help='Headers to extract, one per line.')

conf = ap.parse_args() #sys.argv
f = conf.fasta
h = conf.headers
header_set = set()

with open(conf.headers) as h:
    for line in h:
        line.rstrip()
        header_set.add(line.rstrip())

with open(conf.fasta) as f:
    for line in f:
        if '>' in line:
            print_line = False
            split_header = line.split()
            header_id = split_header[0].replace(">", "")
            if header_id in header_set:
                print( line.rstrip() )
                print_line = True
        elif print_line == True:
            print(line.rstrip())

exit()
