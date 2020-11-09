#!/usr/bin/env bash
#SBATCH -J count_nucl
#SBATCH --partition=short
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=2


# ---------------
# Step 1
# Collect inputs
# ---------------

Read_F=$(basename $1)
Read_R=$(basename $2)
Genome_size=$3
OutDir=$4

#DATA_TYPE=$(echo $Read_F | rev | cut -d "/" -f6 | rev)
#READ_TYPE=$(echo $Read_F | rev | cut -d "/" -f5 | rev)
#ORGANISM=$(echo $Read_F | rev | cut -d "/" -f4 | rev)
#STRAIN=$(echo $Read_F | rev | cut -d "/" -f3 | rev)

CurDir=$PWD

# ---------------
# Step 2
# Copy data
# ---------------

WorkDir=$PWD/${SLURM_JOB_USER}_${SLURM_JOBID}
mkdir -p $WorkDir
cd $WorkDir

cp $CurDir/$1 $WorkDir
cp $CurDir/$2 $WorkDir

gunzip $Read_F
gunzip $Read_R

# ---------------
# Step 3
# Calculate estimated genome coverage
# ---------------

Sub1=*R1*.fq
Sub2=*R2*.fq

/data/scratch/gomeza/prog/count_nucl.pl -i $Sub1 -i $Sub2 -g $3 > estimated_coverage.log

cp -r $WorkDir/estimated_coverage.log $CurDir/$OutDir/.
#cp -r $WorkDir/estimated_coverage.log $CurDir/$DATA_TYPE/$READ_TYPE/$ORGANISM/$STRAIN/.
rm -r $WorkDir

