#!/usr/bin/env bash
#SBATCH -J RAxML
#SBATCH --partition=short
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=4

# Generate a maximum likihood phylogney from a given nucleotide alignment

Usage="RAxML.sh <alignment.fa> <outfile_prefix> <output_directory>"+
echo "$Usage"

# ---------------
# Step 1
# Collect inputs
# ---------------

AlignmentIn=$1
Prefix=$2
OutDir=$3

echo  "RAxML with the following inputs:"
echo "Alignment - $AlignmentIn"
echo "Prefix - $Prefix"
echo "OutDir - $OutDir"

CurPath=$PWD
WorkDir=$TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}

# ---------------
# Step 2
# Copy working Files
# ---------------

mkdir -p $WorkDir
cd $WorkDir
# cp $CurPath/$AlignmentIn $WorkDir/alignment.fa

cat $CurPath/$AlignmentIn \
 | sed "s/${Prefix}://g" \
 | sed "s/genome.ctg.fa://g" \
 | sed "s/_contigs_unmasked.fa//g" \
 | sed -E "s/:.*//g" \
 > $WorkDir/alignment.fa

# ---------------
# Step 3
# Run RaxML
# ---------------
# Full options can be found here:
# https://sco.h-its.org/exelixis/resource/download/NewManual.pdf
# -f a
# Tells to run in bootstrap mode and search for the best-scoring ML tree in
# one single program run
# this means variables need to be set for the random seed and number of
# bootstrap replicates: -N 1000 -x $RANDOM
Bootstraps="1000"
raxmlHPC-AVX2 -s alignment.fa -n $Prefix -m GTRGAMMA -f a -x $RANDOM -N $Bootstraps -p $RANDOM

rm alignment.fa
mkdir -p $CurPath/$OutDir
cp $WorkDir/* $CurPath/$OutDir
