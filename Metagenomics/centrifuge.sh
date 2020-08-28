#!/usr/bin/env bash
#SBATCH -J centrifuge
#SBATCH --partition=medium
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=8

#Classification of DNA sequences using centrifuge

Usage="centrifuge.sh Database star_aligmentUnmapped.out.mate1 star_aligmentUnmapped.out.mate2 Output_directory"
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
# Run centrifure
# ---------------

centrifuge -p 4 \
-x /data/scratch/gomeza/prog/centrifuge/$Database \
-t -1 $ReadMate1 -2 $ReadMate2 \
--phred33 \
--report-file centrifuge_report.tsv \
-S centrifuge_results.txt 

# --min-hitlen <integer> option can be used to set a minimum length of partial hits. 
# Default 22.

# ---------------
# Step 3
# Create a Kraken-style report
# ---------------

echo "Creating Kraken-style report for visualisation"
centrifuge-kreport \
-x  /data/scratch/gomeza/prog/centrifuge/$Database \
centrifuge_results.txt > centrifuge_krakened.txt

# --min-score <integer> option set minimum score for reads to be counted
# --min-length <integer> option set minimum alignment length to the read

# ---------------
# Step 4
# Copy results
# ---------------

cp -r $WorkDir/centrifuge* $CurDir/$OutDir/.

rm -r $WorkDir

echo "Done"