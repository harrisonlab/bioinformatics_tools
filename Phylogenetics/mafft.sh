#!/bin/bash
#SBATCH -J mafft
#SBATCH --partition=long
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=8

#### MAFFT alignment with the highest accuracy method L-INS-i (less than <200 sequences)
#### The scripts separately aligns all the FASTA sequences contained within each file in the directory.
#### OUTPUT: aligned FASTA file, with "aligned" suffix


for input in *.fasta
do
filename=$(basename "$input")
output="${filename%.*}_aligned.fasta"
mafft --localpair --maxiterate 1000 $input >$output
done

#Convert to single line FASTA for easy parsing
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' $output >temp && mv temp $output
