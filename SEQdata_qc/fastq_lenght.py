#!/usr/bin/env python

import sys
######################################################################################
syntax = '''
------------------------------------------------------------------------------------
Usage: python fastq_lenght.py file.fastq 
result: .txt file same name as input name plus "_lenghts.txt" 
------------------------------------------------------------------------------------
'''
######################################################################################

if len(sys.argv) != 2:
    print (syntax)
    sys.exit()

######################################################################################

fastq_file = open(sys.argv[1], 'r')
prefix = sys.argv[1].split('.')[0]
outfile = open(prefix + '_' + 'lenghts.txt', 'w')
seq = ''
for line in fastq_file:
    line = line.rstrip('\n')
    if line.startswith('@'):
        if seq:
            outfile.write(name + '\t' + str(len(seq)) + '\n')
            seq = ""
        name = line 
    else:
        seq = line 
outfile.write(name + '\t' + str(len(seq)) + '\n')
fastq_file.close()
outfile.close()

print ('\n' + '\t' + 'File: ' + prefix + '_' + 'lenghts.txt has been created...')