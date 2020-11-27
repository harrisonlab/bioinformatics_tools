#!/usr/bin/env bash
#SBATCH -J salmon
#SBATCH --partition=long
#SBATCH --mem=10G
#SBATCH -o salmon_"%j".out


Usage='sbatch $ProgDir/salmon.sh <Transcriptome.fa> <Read_F.fq> <Read_R.fq> <OutputFi
lePath>'

echo "$Usage"

# ---------------
# Step 1
# Collect inputs
# ---------------
#export PATH=${PATH}:/home/gomeza/miniconda3/pkgs/salmon-1.1.0-hf69c8f4_0/bin/

Transcriptome=$(basename $1)
ReadF=$(basename $2)
ReadR=$(basename $3)
OutDir=$4

echo $Transcriptome
echo $ReadF
echo $ReadR

CurDir=$PWD
WorkDir=$TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}
echo "$WorkDir"
mkdir -p $WorkDir
cd $WorkDir

cp $CurDir/$1 $Transcriptome
cp $CurDir/$2 $ReadF
cp $CurDir/$3 $ReadR


# ---
# Step 2
# Index the transcriptome
# ---
# Next, we’re going to build an index on our transcriptome.
# The index is a structure that salmon uses to quasi-map RNA-seq reads during quantification.
# The index need only be constructed once per transcriptome,
# and it can then be reused to quantify many experiments.
# We use the index command of salmon to build our index:

salmon index -t $Transcriptome -i transcripts_index --keepDuplicates -k 27
# ---
# Step 3
# Quantifying transcripts
# ---
# -l        If set to A then Automatically detect the library type
# --dumpEq  Write a file in the auxiliary directory, called eq_classes.txt
#           that contains the equivalence classes and corresponding counts that
#           were computed during quasi-mapping.
# --seqBiasPassing  Enable it to learn and correct for sequence-specific biases
#                   in the input data.
# --gcBias  Learn and correct for fragment-level GC biases in the input data

salmon quant \
    -i transcripts_index \
    -l A \
    -1 $ReadF \
    -2 $ReadR \
    --validateMappings \
    -p 4 \
    --numBootstraps 1000 \
    --dumpEq \
    --seqBias \
    --gcBias \
    -o transcripts_quant


mkdir -p $CurDir/$OutDir
mv transcripts_quant/* $CurDir/$OutDir/.

rm -r $WorkDir


