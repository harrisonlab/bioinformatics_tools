#!/usr/bin/env bash
#SBATCH -J mosaic
#SBATCH --partition=long
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=24

# Resolve repeats with mosaic

Usage="mosaic.sh <reads.fastq.gz> <flyegenome.fasta> <outfile_prefix> <output_directory> <genomesize>"
echo "$Usage"

# ---------------
# Step 1
# Collect inputs
# ---------------

RawReads=$1
DraftGenome=$2
Prefix=$3
OutDir=$4
Size=$5
echo  "Running mosaic with the following inputs:"
echo "Raw Reads In - $RawReads"
echo "Genome In - $DraftGenome"
echo "Prefix - $Prefix"
echo "OutDir to - $OutDir"
echo "Estimated genome size - $Size"

CurPath=$PWD
WorkDir=$TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}
mkdir -p $WorkDir
cd $WorkDir

Raw=$(basename $RawReads)
cp $CurPath/$RawReads $Reads
cp $CurPath/$DraftGenome $Genome

# ---------------
# Step 2
# Rename reads
# ---------------

/home/gomeza/miniconda3/envs/olc_assemblers/bin/rename.sh \
in=$Reads \
out="$Prefix"reads_rename.fasta \
prefix=$Prefix

# ---------------
# Step 3
# Run Mosaic
# ---------------

mosaic --reads "$Prefix"reads_rename.fasta -o $WorkDir --genome-size $Size --flye-dir /home/gomeza/miniconda3/envs/dbg_assemblers_py27/bin/flye --contigs $Genome

cp -r $WorkDir/* $CurPath/$OutDir/.

rm -r $WorkDir