#!/usr/bin/env bash
#SBATCH -J fye
#SBATCH --partition=long
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=24

# Assemble Long read data using flye

Usage="fly.sh <read.fa> <outfile_prefix> <output_directory>"
echo "$Usage"

# ---------------
# Step 1
# Collect inputs
# ---------------

RawReads=$1
Prefix=$2
OutDir=$3
Size=$4
type=$5

echo  "Running flye with the following inputs:"
echo "Raw Reads In - $RawReads"
echo "Prefix - $Prefix"
echo "OutDir - $OutDir"
echo "Estimated genome size - $Size"

CurPath=$PWD
WorkDir=$TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}
mkdir -p $WorkDir
cd $WorkDir

Raw=$(basename $RawReads)
cp $CurPath/$RawReads $Raw

# ---------------
# Step 2
# Rename reads
# ---------------

/home/gomeza/miniconda3/envs/olc_assemblers/bin/rename.sh \
in=$Raw \
out="$Prefix"reads_rename.fasta \
prefix=$Prefix


# Run Flye

  if [ $type == "pacbioraw" ]; then
   flye --pacbio-raw "$Prefix"reads_rename.fasta --out-dir $WorkDir --genome-size $Size
  elif [ $type == "pacbiocorrected" ]; then
    flye --pacbio-corr "$Prefix"reads_rename.fasta --out-dir $WorkDir --genome-size $Size
    elif [ $type == "pacbiohifi" ]; then
    flye --pacbio-hifi "$Prefix"reads_rename.fasta --out-dir $WorkDir --genome-size $Size
    elif [ $type == "nanoraw" ]; then
    flye --nano-raw "$Prefix"reads_rename.fasta --out-dir $WorkDir --genome-size $Size
    elif [ $type == "nanocorrected" ]; then
    flye --nano-corr "$Prefix"reads_rename.fasta --out-dir $WorkDir --genome-size $Size
    elif [ $type == "contigs" ]; then
    flye --subassembliesr "$Prefix"reads_rename.fasta --out-dir $WorkDir --genome-size $Size
  else
    echo "data type not supported"
  fi


cp -r $WorkDir/* $CurPath/$OutDir/.

rm -r $WorkDir