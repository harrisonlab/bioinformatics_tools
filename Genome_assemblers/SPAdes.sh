#!/bin/bash

#Assemble contigs using SPAdes

#!/usr/bin/env bash
#SBATCH -J SPAdes
#SBATCH --partition=himem
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=16

Usage="SPAdes.sh <F_read.fa> <R_read.fa> <output_directory> <correct/only-assembler> [<coverage_cutoff>]"
echo "$Usage"

R1=$1
R2=$2
OutDir=$3
Correction=$4
Cutoff='auto'
if [ $5 ]; then
  Cutoff=$5
fi

CurPath=$PWD
WorkDir=$TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}

mkdir -p $WorkDir

F_Read=$(basename $R1)
R_Read=$(basename $R2)

cp $CurPath/$R1 $WorkDir/$F_Read
cp $CurPath/$R2 $WorkDir/$R_Read


echo  "Running SPADES with the following in='$R1 $R2' $OutDir "
echo "You have set read correction to: $Correction"
echo "Coverage cutoff set to $Cutoff"

if [[ "$Correction" == 'correct' ]]; then
  spades.py -k 21,33,55,77,99,127 -m 94 --phred-offset 33 --careful -1 $WorkDir/$F_Read -2 $WorkDir/$R_Read -t 16  -o $WorkDir/. --cov-cutoff "$Cutoff"
elif [[ "$Correction" == 'only-assembler' ]]; then
  spades.py -k 21,33,55,77,99,127 -m 94 --phred-offset 33 --careful --only-assembler -1 $WorkDir/$F_Read -2 $WorkDir/$R_Read -t 16  -o $WorkDir/. --cov-cutoff "$Cutoff"
else
  echo "Please set fourth option - whether you require read correction [correct / only-assembler]"
  exit
fi

#echo "Filtering contigs smaller than 500bp"
#mkdir -p $WorkDir/filtered_contigs
#FilterDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/abyss
#$FilterDir/filter_abyss_contigs.py $WorkDir/scaffolds.fasta 500 > $WorkDir/filtered_contigs/contigs_min_500bp.fasta


rm $WorkDir/$F_Read
rm $WorkDir/$R_Read
mkdir -p $CurPath/$OutDir
cp -r $WorkDir/* $CurPath/$OutDir/.
echo "files copied"
