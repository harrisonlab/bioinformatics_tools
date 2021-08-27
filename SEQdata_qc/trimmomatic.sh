#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 8
#$ -l virtual_free=1G

#!/usr/bin/env bash
#SBATCH -J trimmomatic
#SBATCH --partition=long
#SBATCH --mem-per-cpu=2G
#SBATCH --cpus-per-task=8

# Commands to perform adapter trimming of either RNA of DNA data using trimmomatic
# Visualises the quality of data before and after trimming using fastqc.

USAGE="trimmomatic.sh <F_reads.fq/fq.gz> <R_reads.fq/fq.gz> <Illumina_adapters.fa> <Output_directory>"
echo "$USAGE"
echo ""

F_IN=$1
R_IN=$2
ADAPTER_FILE=$3
OutDir=$4
CUR_PATH=$PWD
WORK_DIR=$TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}

ORGANISM=$(echo $F_IN | rev | cut -d "/" -f4 | rev)
STRAIN=$(echo $F_IN | rev | cut -d "/" -f3 | rev)

mkdir -p $WORK_DIR
cd $WORK_DIR

LOGFILE="$STRAIN"_trim_log

F_NO_ADAPT="$STRAIN"_no_adapt.fq
R_NO_ADAPT="$STRAIN"_no_adapt.fq
UNPAIRED_F_NO_ADAPT="$STRAIN"_F_no_adapt_unpaired.fq
UNPAIRED_R_NO_ADAPT="$STRAIN"_R_no_adapt_unpaired.fq

trimmomatic \
  PE -phred33 -trimlog "$LOGFILE"_1.txt \
  $CUR_PATH/$F_IN $CUR_PATH/$R_IN \
  $F_NO_ADAPT $UNPAIRED_F_NO_ADAPT \
  $R_NO_ADAPT $UNPAIRED_R_NO_ADAPT \
  ILLUMINACLIP:$ADAPTER_FILE:2:30:10 \
  LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:40
# Explanation of settings:	ILLUMINACLIP:$ADAPTER_FILE:<seed_mismatches>:<palindrome_clip_threshold>:<adapter_clip_threshold> MINLEN:<minmum_finaR_length>

OutDirF=$CUR_PATH/$OutDir/F
OutDirR=$CUR_PATH/$OutDir/R

mkdir -p $OutDirF
cp $F_NO_ADAPT $OutDirF/.
cp $UNPAIRED_F_NO_ADAPT $OutDirF/.

mkdir -p $OutDirR
cp $R_NO_ADAPT $OutDirR/.
cp $UNPAIRED_R_NO_ADAPT $OutDirR/.

rm -r $WORK_DIR
