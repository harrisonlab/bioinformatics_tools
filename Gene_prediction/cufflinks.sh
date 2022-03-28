#!/usr/bin/env bash
#SBATCH -J cufflinks
#SBATCH --partition=long
#SBATCH --mem-per-cpu=6G
#SBATCH --cpus-per-task=40


# Transcript assembly and estimation of abundance

CurDir=$PWD
WorkDir=$TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}

AcceptedHits=$1
OutDir=$2

Usage="cufflinks.sh <Aligned_reads.bam> <Out_directory>"

echo "$Usage"
echo "Aligned reads: $AcceptedHits"
echo "Output directory: $OutDir"

mkdir -p $WorkDir/out
cd $WorkDir
mkdir -p $CurDir/$OutDir

cp $CurDir/$AcceptedHits accepted_hits.bam

/home/agomez/scratch/apps/prog/cufflinks/cufflinks-2.2.1.Linux_x86_64/cufflinks -o out -p 16 --max-intron-length 4000 accepted_hits.bam 2>&1 | tee out/cufflinks.log

rm accepted_hits.bam
cp -r $WorkDir/out/* $CurDir/$OutDir/.

rm -r $WorkDir
