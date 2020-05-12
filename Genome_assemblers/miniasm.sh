#!/usr/bin/env bash
#SBATCH -J miniasm
#SBATCH --partition=long
#SBATCH --mem=36G

# Assemble Long read data using SMRTdenovo

Usage="miniasm.sh <rawreads.fastq.qz> <outfile_prefix> <output_directory>"
echo "$Usage"

# ---------------
# Step 1
# Collect inputs
# ---------------

RawReads=$1
Prefix=$2
OutDir=$3
echo  "Running miniasm with the following inputs:"
echo "FastaIn - $FastaIn"
echo "Prefix - $Prefix"
echo "OutDir - $OutDir"

CurPath=$PWD
WorkDir=$TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}
mkdir -p $WorkDir
cd $WorkDir

Raw=$(basename $RawReads)
cp $CurPath/$RawReads $Raw

# ---------------
# Step 2
# Rename reads
# ---------------

rename.sh in=$Raw out=$WorkDir prefix=$Prefix

# ---------------
# Step 3
# Fast all-against-all overlap of raw reads
# ---------------

minimap2 -x ava-ont -t8 $Prefix.fasta $Prefix.fasta | gzip -1 > $Prefix.paf.gz

# ---------------
# Step 4
# Concatenate pieces of read sequences to generate the final sequences
# ---------------

miniasm -f $Prefix.fasta $Prefix.paf.gz > reads.gfa

# ---------------
# Step 5
# Convert gfa file to fasta file
# ---------------

awk '/^S/{print ">"$2"\n"$3}' reads.gfa | fold > $PWD/$OutDir/$Prefix.fa

#rm -r $WorkDir
