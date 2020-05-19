#!/usr/bin/env bash
#SBATCH -J braker
#SBATCH --partition=long
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=8

Usage="codingquarry.sh InGenome.fa StringTie/Cufflinks.gtf Output_directory"
echo $Usage

set -u
set -e

CurDir=$PWD
WorkDir=$PWD/${SLURM_JOB_USER}_${SLURM_JOBID}

mkdir -p $WorkDir
cd $WorkDir


Assembly=$CurDir/$1
GTF=$CurDir/$2
OutDir=$CurDir/$3
Threads=8
echo "Running CodingQuaryPM with the following options:"
echo "Assembly - $Assembly"
echo "StringTie/Cufflinks GTF file - $GTF"
echo "Output directory - $OutDir"
# CufflinksGFF3=transcripts.gff3

CufflinksGTF_to_CodingQuarryGFF3.py $GTF > transcripts.gff3
cp $Assembly assembly.fa


/home/gomeza/miniconda3/envs/antismash_py27/bin/run_CQ-PM_stranded.sh assembly.fa transcripts.gff3 $Threads 2>&1 | tee codingquaryPM_log.txt

mkdir -p $OutDir/out
mv codingquaryPM_log.txt $OutDir
mv out/* $OutDir/out/.
#rm -r $WorkDir
