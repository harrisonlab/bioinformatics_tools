#!/usr/bin/env bash
#SBATCH --partition=medium
#SBATCH --time=24:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=8

# Assess assembly quality using quast.

Assembly=$1
OutDir=$2

CurPath=$PWD
WorkDir=$TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}

echo "running quast "
quast.py --threads 8 --eukaryote -o $WorkDir $1
echo "quast run"

mkdir -p $CurPath/$OutDir
cp -r $WorkDir/* $CurPath/$OutDir/.
echo "files copied"

rm -r $WorkDir
