#!/usr/bin/env bash
#SBATCH -J count_nucl
#SBATCH --partition=short
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=2


# ---------------
# Step 1
# Collect inputs
# ---------------

CurPath=$PWD

InFile=$CurPath/$1
Genome_size=$2
OutDir=$3

FileName=$(basename $InFile)

OutName=$(echo "$FileName" | sed "s/.fastq.*//g" | sed "s/.fq.*//g")

# ---------------
# Step 2
# Copy data
# ---------------

WorkDir=$PWD/${SLURM_JOB_USER}_${SLURM_JOBID}
mkdir -p "$WorkDir"
cd "$WorkDir"

cp $InFile seqreads.fastq.gz
cat seqreads.fastq.gz | gunzip -fc > "$OutName".fastq

# ---------------
# Step 3
# Calculate estimated genome coverage
# ---------------


/data/scratch/gomeza/prog/count_nucl.pl -i "$OutName".fastq -g $Genome_Size > "$OutName"_cov.txt


mkdir -p $CurPath/$OutDir
cp "$OutName"_cov.txt $CurPath/$OutDir/.

rm -r $WorkDir
