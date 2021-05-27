#! /usr/bin/env python

import argparse

parser=argparse.ArgumentParser(description="program that replaces fasta records")
parser.add_argument("-i", help="input fasta", type=file)
parser.add_argument("-r", help="replacement records file", type=file)
parser.add_argument("-o", help="output file")
args = parser.parse_args()
newfasta=open(args.o,'w') 

for line in args.i:
    if line.startswith('>'):
        newname=args.r.readline()
        newfasta.write(newname)
    else: 
        newfasta.write(line)