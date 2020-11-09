#!/usr/bin/env bash
#SBATCH -J kraken2
#SBATCH --partition=short
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=8

#Classification of DNA sequences using Kraken2

Usage="kraken.sh star_aligmentUnmapped.out.mate1 star_aligmentUnmapped.out.mate2 DB Output_directory"
#echo $Usage

# ---------------
# Step 1
# Collect inputs
# ---------------

ReadMate1=$(basename $1)
ReadMate2=$(basename $2)
Database=$3
OutDir=$4

mkdir -p $OutDir

echo $ReadMate1
echo $ReadMate2

# Set working directory
CurDir=$PWD
WorkDir=$TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}
mkdir -p $WorkDir
echo $WorkDir

# Copy over input files
cd $WorkDir

cp $CurDir/$1 $ReadMate1
cp $CurDir/$2 $ReadMate2

# ---------------
# Step 2
# Release the kraken!!!
# ---------------

kraken2 --db /data/scratch/gomeza/prog/kraken2/$DataBase \
--threads 4 \
--paired --classified-out cseqs#.fq $ReadMate1 $ReadMate2 \
--output KrakenResults.txt

# ---------------
# Step 3
# Copy results
# ---------------

cp -r $WorkDir/* $CurDir/$OutDir/.

rm -r $WorkDir

echo "Done"