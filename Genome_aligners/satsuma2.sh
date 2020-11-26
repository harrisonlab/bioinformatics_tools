#!/usr/bin/env bash
#SBATCH -J satsuma2
#SBATCH --partition=long
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=8

Sequence=$(basename $1)
Query=$(basename $2)
OutDir=$3

CurDir=$PWD

WorkDir=$TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}
mkdir -p $WorkDir
cd $WorkDir
cp $CurDir/$1 $Sequence
cp $CurDir/$2 $Query

mkdir -p out

SatsumaSynteny2 -t $Sequence -q $Query -o out

cp -r out $CurDir/$OutDir/.
cd $CurDir
rm -r $WorkDir