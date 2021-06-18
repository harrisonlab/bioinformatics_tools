#!/usr/bin/env bash
#SBATCH -J kraken2
#SBATCH --partition=short
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=8

#Classification of DNA sequences using Kraken2

Usage="kraken.sh long_read.fastq DB Output_directory"
#echo $Usage

# ---------------
# Step 1
# Collect inputs
# ---------------

Read=$(basename $1)
Database=$2
OutDir=$3

mkdir -p $OutDir

echo $Read

# Set working directory
CurDir=$PWD
WorkDir=$TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}
mkdir -p $WorkDir
echo $WorkDir

# Copy over input files
cd $WorkDir

cp $CurDir/$1 $Read

# ---------------
# Step 2
# Release the kraken!!!
# ---------------

kraken2 --db /data/scratch/gomeza/prog/kraken2/Foxysporum \
--threads 4 \
--classified-out classifiedseqs.fq \
--unclassified-out unclassifiedseqs.fq \
--gzip-compressed $Read \
--report KrakenReport.txt
--output KrakenResults.txt

# ---------------
# Step 3
# Copy results
# ---------------

cp -r $WorkDir/* $CurDir/$OutDir/.

rm -r $WorkDir

echo "Done"