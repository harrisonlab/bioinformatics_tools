#!/usr/bin/env bash
#SBATCH -J stringtie
#SBATCH --partition=himem
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=16

# Aligned transcripts assembly using a reference genome.

CurDir=$PWD
WorkDir=$PWD/${SLURM_JOB_USER}_${SLURM_JOBID}

AcceptedHits=$1
OutDir=$2

Usage="stringtie.sh <Aligned_reads.bam> <Out_directory>"

echo "$Usage"
echo "Aligned reads: $AcceptedHits"
echo "Output directory: $OutDir"

mkdir -p $WorkDir/out
cd $WorkDir
mkdir -p $CurDir/$OutDir

cp $CurDir/$AcceptedHits accepted_hits.bam

stringtie -o out.gtf -p 8 accepted_hits.bam

sed -i 1,2d out.gtf

cp -r $WorkDir/* $CurDir/$OutDir/.
rm -r $WorkDir