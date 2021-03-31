#!/bin/bash

# Assemble contigs using SPAdes

#!/usr/bin/env bash
#SBATCH -J hybridSPAdes
#SBATCH --partition=himem
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=16


Usage="hybridSPAdes.sh <longreads.fq>  <F1_read.fa> <R1_read.fa> <readtype> <output_directory> [<coverage_cutoff>]"
echo "$Usage"

N1=$1
F1=$2
R1=$3
type=$4
OutDir=$5
Cutoff='auto'
if [ $6 ]; then
  Cutoff=$6
fi

CurPath=$PWD
WorkDir=$TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}

mkdir -p $WorkDir

Long_Read1=$(basename $N1)
F1_Read=$(basename $F1)
R1_Read=$(basename $R1)

cp $CurPath/$N1 $WorkDir/$Long_Read1
cp $CurPath/$F1 $WorkDir/$F1_Read
cp $CurPath/$R1 $WorkDir/$R1_Read

echo  "Running SPADES with the following inputs:"
echo "Long reads = $N1"
echo "F1_Read = $F1"
echo "R1_Read = $R1"
echo "Output directory will be: $CurPath/$OutDir"
echo "Coverage cutoff set to $Cutoff"

if [ $type == "nanopore" ]; then
spades.py \
    -k 21,33,55,77,99,127 \
    -m 375 \
    --phred-offset 33 \
    --careful \
    --nanopore $WorkDir/$Long_Read1 \
    --pe1-1 $WorkDir/$F1_Read \
    --pe1-2 $WorkDir/$R1_Read \
    -t 24  \
    -o $WorkDir/. \
    --cov-cutoff "$Cutoff"
elif [ $type == "pacbio" ]; then
spades.py \
    -k 21,33,55,77,99,127 \
    -m 375 \
    --phred-offset 33 \
    --careful \
    --pacbio $WorkDir/$Long_Read1 \
    -1 $WorkDir/$F1_Read \
    -2 $WorkDir/$R1_Read \
    -t 24  \
    -o $WorkDir/. \
    --cov-cutoff "$Cutoff"
  else
    echo "data type not supported"
  fi
  

rm $WorkDir/$Long_Read1
rm $WorkDir/$F1_Read
rm $WorkDir/$R1_Read
mkdir -p $CurPath/$OutDir
cp -r $WorkDir/* $CurPath/$OutDir/.
echo "files copied"
