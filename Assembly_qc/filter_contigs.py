#!/usr/bin/python

'''
	filter_abyss.py

This script will accept a fasta file and filter accessions
if they have a sequence length longer than X, where X is
specified in the input.

Sequence data will be redirected to standard output.

As the program reads in files, it will identify those lines
beginning with '>' and store those lines as fasta headers.
It will then check the length of the subsequent line
and if >= to X, will print the header and line to screen.

Sequence data can not be wrapped over multiple lines.

'''

import sys
from os import path

# -----
# Step 1
#	Set variables
# -----

input_file = sys.argv[1]
f = open(input_file)
min_lgth = sys.argv[2]
prev_header = ""
prev_line = ""


# -----
# Step 2
#	Filter file
# -----

for line in f:
	if line.startswith('>'):
		if len(prev_line) >= int(min_lgth):
			print prev_header + "\n" + prev_line
		prev_line = ""
		prev_header = line.rstrip()
	else:
		prev_line += line.rstrip()

# -----
# Step 3
#	Tidy up
# -----

f.close()
