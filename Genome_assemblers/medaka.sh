#!/usr/bin/env bash
#SBATCH -J medaka
#SBATCH --partition=short
#SBATCH --mem-per-cpu=10G
#SBATCH --cpus-per-task=8

# create a consensus sequence from racon corrected assemblies 

Usage="medaka.sh <assembly.fa> <corrected_reads.fq.gz> <output_directory>"
echo "$Usage"

# ---------------
# Step 1
# Collect inputs
# ---------------

AssemblyIn=$1
ReadsIn=$2
OutDir=$3

echo "Assembly - $AssemblyIn"
echo "Fasta reads - $ReadsIn"
echo "OutDir - $OutDir"

CurDir=$PWD
WorkDir=$TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}
mkdir -p $WorkDir
cd $WorkDir

Assembly=$(basename $AssemblyIn)
Reads=$(basename $ReadsIn)
cp $CurDir/$AssemblyIn $Assembly
cp $CurDir/$ReadsIn $Reads

#Prefix=$(echo $Assembly | cut -f1 -d '.')

mkdir -p $CurDir/$OutDir

medaka_consensus -i $Reads -d $Assembly -t 8

cp -r $WorkDir/* $CurDir/$OutDir/.

rm -r $WorkDir
