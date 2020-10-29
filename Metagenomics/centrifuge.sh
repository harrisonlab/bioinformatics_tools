#!/usr/bin/env bash
#SBATCH -J centrifuge
#SBATCH --partition=himem
#SBATCH --mem-per-cpu=20G
#SBATCH --cpus-per-task=20

#Classification of DNA sequences using centrifuge

Usage="centrifuge.sh Database star_aligmentUnmapped.out.mate1 star_aligmentUnmapped.out.mate2 Output_directory"
#echo $Usage

# ---------------
# Step 1
# Collect inputs
# ---------------

Database=$1
OutDir=$2
Reads=$(basename $3)

echo $Reads

mkdir -p $OutDir


# Set working directory
CurDir=$PWD
WorkDir=$TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}
mkdir -p $WorkDir
echo $WorkDir

# Copy over input files
cd $WorkDir

cp $CurDir/$3 $Reads
if [ $4 ]; then
ReadMate=$(basename $4)
cp $CurDir/$4 $ReadMate
echo $ReadMate
fi

# ---------------
# Step 2
# Run centrifure
# ---------------

if [ $4 ]; then
centrifuge -p 4 \
-x /data/scratch/gomeza/prog/centrifuge/$Database \
-t -1 $Reads -2 $ReadMate \
--phred33 \
--report-file centrifuge_report.tsv \
-S centrifuge_results.txt 
else
centrifuge -p 4 \
-x /data/scratch/gomeza/prog/centrifuge/$Database \
-U $Reads \
--phred33 \
--report-file centrifuge_report.tsv \
-S centrifuge_results.txt 
fi

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