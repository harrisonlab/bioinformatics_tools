#!/usr/bin/env bash
#SBATCH -J qualimap
#SBATCH --partition=long
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=8

bam=$1
OutDir=$2

CWD=$PWD

WorkDir=$TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}
mkdir -p $WorkDir

cp -r $bam $WorkDir
bm=$(basename "$bam")

cd $WorkDir

qualimap bamqc -bam $bm -outfile result.pdf -c

mkdir -p $CWD/$OutDir
cp -r * $CWD/$OutDir/.
rm -r $WorkDir
