#!/usr/bin/env bash
#SBATCH -J braker
#SBATCH --partition=medium
#SBATCH --mem-per-cpu=6G
#SBATCH --cpus-per-task=16


Assembly=$1
OutDir=$2
AcceptedHits=$3
GeneModelName=$4
CurDir=$PWD

echo "$Assembly"
echo "$OutDir"
echo "$AcceptedHits"
echo "$GeneModelName"
echo "$CurDir"

WorkDir=$CurDir/${SLURM_JOB_USER}_${SLURM_JOBID}
mkdir -p $WorkDir
cd $WorkDir

cp $CurDir/$Assembly assembly.fa
cp $CurDir/$AcceptedHits alignedRNA.bam



braker.pl \
  --GENEMARK_PATH=/home/agomez/scratch/apps/prog/gmes_linux_64 \
  --overwrite \
  --fungus \
  --gff3 \
  --softmasking on \
  --species=$GeneModelName \
  --genome="assembly.fa" \
  --bam=alignedRNA.bam
  
#--BAMTOOLS_PATH=/home/gomeza/miniconda3/envs/gene_pred/bin \
mkdir -p $CurDir/$OutDir
cp -r braker/* $CurDir/$OutDir/.

rm -r $WorkDir
