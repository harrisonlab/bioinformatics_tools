#!/usr/bin/env bash
#SBATCH -J miniasm
#SBATCH --partition=long
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=24

# Assemble Long read data using miniasm

Usage="miniasm.sh <read.fa> <outfile_prefix> <output_directory>"
echo "$Usage"

# ---------------
# Step 1
# Collect inputs
# ---------------

RawReads=$1
Prefix=$2
OutDir=$3
echo  "Running miniasm with the following inputs:"
echo "FastaIn - $RawReads"
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

rename.sh in=$Raw out="$Prefix"_rename.fasta prefix=$Prefix

cp "$Prefix"_rename.fasta $CurPath

# ---------------
# Step 3
# Fast all-against-all overlap of raw reads
# ---------------

minimap2 -x ava-ont -t8 "$Prefix"_rename.fasta "$Prefix"_rename.fasta | gzip -1 > $Prefix.paf.gz

cp $Prefix.paf.gz $CurPath

# ---------------
# Step 4
# Concatenate pieces of read sequences to generate the final sequences
# ---------------

miniasm -f "$Prefix"_rename.fasta $Prefix.paf.gz > reads.gfa

cp reads.gfa $CurPath

# ---------------
# Step 5
# Convert gfa file to fasta file
# ---------------

awk '/^S/{print ">"$2"\n"$3}' reads.gfa | fold > $Prefix.fa

cp -r $WorkDir/* $CurPath/$OutDir/.

rm -r $WorkDir
