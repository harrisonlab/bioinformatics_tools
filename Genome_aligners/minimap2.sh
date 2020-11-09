#!/bin/bash
#SBATCH -J minimap2
#SBATCH --partition=medium
#SBATCH --mem=5G
#SBATCH --cpus-per-task=20

# Align ONT reads to an assembly.

# ---------------
# Step 1
# Collect inputs
# ---------------

Assembly=$(basename $1)
Reads=$(basename $2)
OutDir=$3


CurDir=$PWD
echo  "Running minimap2 with the following inputs:"
echo "Assembly - $Assembly"
echo "Reads - $Reads"
echo "OutDir - $OutDir"

# ---------------
# Step 2
# Copy data
# ---------------

CurPath=$PWD
WorkDir=${TMPDIR}/${SLURM_JOB_ID}
mkdir -p $WorkDir
cd $WorkDir
cp $CurDir/$1 $Assembly
cp $CurDir/$2 $Reads


# ---------------
# Step 3
# Align seq reads
# ---------------
# Prepare the assembly for alignment
# Align reads against the assembly
# Convert the SAM file to BAM in preparation for sorting.
# Sort the BAM file, in preparation for SNP calling:
# Index the bam file


minimap2 -ax map-ont $Assembly $Reads > ${Assembly}_aligned.sam
samtools view --threads 24 -bS ${Assembly}_aligned.sam -o ${Assembly}_aligned.bam
samtools sort --threads 24 -o ${Assembly}_aligned_sorted.bam ${Assembly}_aligned.bam

rm $Assembly
rm $Reads
rm ${Assembly}_aligned.sam
mkdir -p $CurDir/$OutDir
cp -r $WorkDir/* $CurDir/$OutDir/.
