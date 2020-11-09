#!/usr/bin/env bash
#SBATCH -J effectorP
#SBATCH --partition=medium
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=8

CurPath=$PWD
InFile="$1"
BaseName="$2"
OutDir=$CurPath/"$3"

Organism=$(echo $InFile | rev | cut -d "/" -f3 | rev)
Strain=$(echo $InFile | rev | cut -d "/" -f2 | rev)
InName=$(echo $InFile | rev | cut -d "/" -f1 | rev)

WorkDir=$PWD/${SLURM_JOB_USER}_${SLURM_JOBID}
mkdir -p $WorkDir
cd $WorkDir

cp $CurPath/$InFile proteins.fa
echo "Running effectorP for: $Organism - $Strain"
EffectorP.py -o "$BaseName".txt -E "$BaseName".fa -i proteins.fa

rm proteins.fa
mkdir -p $OutDir
cp $WorkDir/* $OutDir/.
rm -r $WorkDir
