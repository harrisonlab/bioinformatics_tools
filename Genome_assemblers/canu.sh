#!/usr/bin/env bash
#SBATCH -J canu
#SBATCH --partition=long
#SBATCH --mem-per-cpu=16G
#SBATCH --cpus-per-task=12

# Assemble PacBio data using Canu

Usage="canu.sh <reads.fq> <Genome_size[e.g.45m]> <outfile_prefix> <output_directory> [<specification_file.txt>]"
echo "$Usage"

# ---------------
# Step 1
# Collect inputs
# ---------------

FastqIn=$1
Size=$2
Prefix=$3
type=$4
OutDir=$5
AdditionalCommands=""
if [ $6 ]; then
  SpecFile=$6
  AdditionalCommands="-s $SpecFile"
fi
echo  "Running Canu with the following inputs:"
echo "FastqIn - $FastqIn"
echo "Size - $Size"
echo "Prefix - $Prefix"
echo "OutDir - $OutDir"

CurPath=$PWD
WorkDir=$TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}

# ---------------
# Step 2
# Run Canu
# ---------------

mkdir -p $WorkDir
cd $WorkDir
Fastq=$(basename $FastqIn)
cp $CurPath/$FastqIn $WorkDir/$Fastq

if [ $type == "nanopore" ]; then
   canu \
  useGrid=false \
  $AdditionalCommands \
  -overlapper=mhap \
  -utgReAlign=true \
  -d $WorkDir/assembly \
  -p $Prefix genomeSize="$Size" \
  -nanopore $Fastq \
  2>&1 | tee canu_run_log.txt
  elif [ $type == "pacbio" ]; then
   canu \
  useGrid=false \
  $AdditionalCommands \
  -overlapper=mhap \
  -utgReAlign=true \
  -d $WorkDir/assembly \
  -p $Prefix genomeSize="$Size" \
  -pacbio $Fastq \
  2>&1 | tee canu_run_log.txt
  else
    echo "data type not supported"
  fi


mkdir -p $CurPath/$OutDir
cp canu_run_log.txt $CurPath/$OutDir/.
cp $WorkDir/assembly/* $CurPath/$OutDir/.
