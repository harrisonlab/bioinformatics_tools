#!/usr/bin/env bash
#SBATCH -J mosaic
#SBATCH --partition=long
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=24

#export LD_LIBRARY_PATH=$HOME/install/glibc-2.27/lib:$HOME/install/glibc-2.27/lib/ld-linux-x86-64.so.2:$LD_LIBRARY_PATH

LD_PRELOAD="/home/gomeza/install/glibcpath/glibc-2.27/lib/libc.so.6 /home/gomeza/install/glibcpath/glibc-2.27/lib/libpthread.so.0  /home/gomeza/install/glibcpath/glibc-2.27/lib/ld-linux-x86-64.so.2"
export LD_PRELOAD

# Resolve repeats with mosaic

Usage="mosaic.sh <reads.fastq.gz> <flyegenome.fasta> <outfile_prefix> <output_directory> <genomesize>"
echo "$Usage"

# ---------------
# Step 1
# Collect inputs
# ---------------

RawReads=$1
FlyeDir=$2
Prefix=$3
OutDir=$4

echo  "Running mosaic with the following inputs:"
echo "Raw Reads In - $RawReads"
echo "Genome In - $FlyeDir"
echo "Prefix - $Prefix"
echo "OutDir to - $OutDir"


CurPath=$PWD
WorkDir=$TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}
mkdir -p $WorkDir
cd $WorkDir

Reads=$(basename $RawReads)
Genome=$(basename $FlyeDir)
cp $CurPath/$RawReads $Reads
cp $CurPath/$FlyeDir $Genome

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

mosaic --reads "$Prefix"reads_rename.fasta -o $WorkDir  --contigs $Genome

cp -r $WorkDir/* $CurPath/$OutDir/.

rm -r $WorkDir
