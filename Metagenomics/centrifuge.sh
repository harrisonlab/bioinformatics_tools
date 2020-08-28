#!/usr/bin/env bash
#SBATCH -J centrifuge
#SBATCH --partition=medium
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=8

#Classification of DNA sequences using centrifuge

Usage="centrifuge.sh Database star_aligmentUnmapped.out.mate1 star_aligmentUnmapped.out.mate2 Output_directory"
echo $Usage

# ---------------
# Step 1
# Collect inputs
# ---------------

Database=$1
ReadMate1=$(basename $2)
ReadMate2=$(basename $3)
OutDir=$4

echo $ReadMate1
echo $ReadMate2

# Set working directory
CurDir=$PWD
WorkDir=$TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}
mkdir -p $WorkDir
echo $WorkDir

# Copy over input files
cd $WorkDir

cp $CurDir/$2 $ReadMate1
cp $CurDir/$3 $ReadMate1

# ---------------
# Step 2
# Run centrifure
# ---------------

centrifuge -p 4 \
-x /data/scratch/gomeza/prog/centrifuge/$Database \
-t -q -1 $ReadMate1 -2 $ReadMate2 \
--phred33 \
--report-file centrifuge_report.tsv \
-S centrifuge_results.txt 

# ---------------
# Step 3
# Create a Kraken-style report
# ---------------

centrifuge-kreport \
-x  /data/scratch/gomeza/prog/centrifuge/$Database \
centrifuge_results.txt > centrifuge_krakened.txt

cp -r $WorkDir/* $CurDir/$OutDir/.

rm -r $WorkDir