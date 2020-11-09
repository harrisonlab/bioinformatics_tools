#!/usr/bin/python

'''
This tool search a file output by bedtools genomecov for regions of the genome
with a low coverage in the middle of a contig and flag these as potential
missassemblies
'''

import sys,argparse
from collections import defaultdict

#######################################
#            Import variables         #
#                                     #
#                                     #
#######################################

ap = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.ArgumentDefaultsHelpFormatter)

ap.add_argument('--genomecov',required=True,type=str,help='An output file of bedtools genomecov')
ap.add_argument('--min',required=True,type=int,help='The coverage below which regions will be flagged for manual inspection')

conf = ap.parse_args() #sys.argv

min_cov = conf.min

contig_list = []
contig_len_dict = {}
regions_dict = defaultdict(list)


#######################################
#     Read the gonomecov file         #
#                                     #
#                                     #
#######################################
# Open the file & initialise variables
# For each line in file
# Identify if the current contig is different from the previous contig
# If so, store the length of the previous contig in a dictionary
# Identify if the coverage at the current position is less than X
# If so, store the position in a dictionary
#
# After the last line has been read store the length of the final contig in
# the contig length dictionary


with open(conf.genomecov) as f:
        contig = "first"
        contig_pos = 0
        for line in f:
            line = line.rstrip()
            split_line = line.split("\t")
            prev_contig = contig
            prev_pos = contig_pos
            contig = split_line[0]
            contig_pos = split_line[1]
            cov = int(split_line[2])
            if contig != prev_contig:
                contig_list.append(contig)
                contig_len_dict[prev_contig] = prev_pos
            if cov < min_cov:
                regions_dict[contig].append(contig_pos)

prev_contig = contig
prev_pos = contig_pos
contig_len_dict[prev_contig] = prev_pos

#######################################
#     Merge stored positions          #
#             into ranges             #
#                                     #
#                                     #
#######################################
# For each contig, check the dictionary for any stored positions of low coverage
# If a default number has been set (-1000) then set the minimum value of the
# range to the current number
# If a low coverage position is within 1000 bp of the previous position then
# maintain the current minimum value for the current ranges
# If a low coverage position is more than 1000 bp from the previous position then
# set the maximum value of the previous range to the previous position and
# initiate a new range.
# print all outranges
# If a contig does not have an entry in the dictionary then print the contig name
# and length

for contig in contig_list:
    if regions_dict[contig]:
        numberlist = regions_dict[contig]
        numberset = set(numberlist)

        # print numberlist
        # print numberset


        outrange = []
        outindels = []
        min_num = int(-1000)
        prev_num = int(0)
        for num in numberlist:
            down_4bp = int(num) -4
            up_4bp = int(num) +4
            # print down_2bp
            # print num
            # print up_2bp
            # quit()
            if int(num) > 4 and str(down_4bp) not in numberset and str(up_4bp) not in numberset:
                outindels.append(str(num))
                continue
            if min_num == -1000:
                min_num = num
            elif (int(num) - int(prev_num)) < 1000:
                pass
            else:
                outrange.append(str(min_num) + "-" + str(prev_num))
                min_num = num
            prev_num = num
        outrange.append(str(min_num) + "-" + str(prev_num))

        contig_len = contig_len_dict[contig]
        print contig + "\t- length: " + str(contig_len) + "\t" + ", ".join(outrange)
        print "\tindels\t" + ", ".join(outindels)
    else:
        contig_len = contig_len_dict[contig]
        print contig + "\t- length: " + str(contig_len)
