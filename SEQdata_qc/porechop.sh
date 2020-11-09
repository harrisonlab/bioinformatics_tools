#!/usr/bin/env bash
#SBATCH -J porechop
#SBATCH --partition=medium
#SBATCH --mem-per-cpu=2G
#SBATCH --cpus-per-task=16


# Wrapper to submit Porechop jobs to the cluster
# Porechop splits MinION reads at locations of adapter sequences

Usage="porechop.sh <raw_reads.fastq.gz> <output_directory>"
echo $Usage

#######  Step 1	 ########
# Initialise values
#########################

RawReads=$(basename $1)
OutDir=$2
CurDir=$PWD
WorkDir=$TMPDIR/porechop

Prefix=$(echo $RawReads | sed 's/.fq.gz//g' | sed 's/.fastq.gz//g')

echo "Raw reads - $RawReads"
echo "Output directory - $OutDir"
echo "Output files will carry the prefix - $Prefix"

#######  Step 2	 ########
# Copy data on worker node
#########################

mkdir -p $WorkDir
cd $WorkDir
cp $CurDir/$1 $RawReads

#######  Step 3	 ########
# Run Porechop
#########################

unset PYTHONPATH
porechop-runner.py -i $RawReads -o "$Prefix"_trim.fastq.gz --threads 16 > "$Prefix"_trim_log.txt

#######  Step 3	 ########
# Tidy up
#########################

mkdir -p $CurDir/$OutDir
cp "$Prefix"_trim.fastq.gz $CurDir/$OutDir/.
cp "$Prefix"_trim_log.txt $CurDir/$OutDir/.
