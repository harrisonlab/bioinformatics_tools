#!/usr/bin/env bash
#SBATCH -J HaplotypeCaller
#SBATCH --partition=long 
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=16

# Call SNPs independently using HaplotyCaller with the GVCF option

Usage="HaplotypeCaller_GVCF.sh <Reference.fasta> <Prefix> <alignment.bam> <OutDir>"
echo "$Usage"

# ---------------
# Step 1
# Collect inputs
# ---------------

ReferenceIn=$1
Prefix=$2
AlignmentIn=$3
OutDir=$4

echo  "Running HaplotypeCaller with the following inputs:"
echo "Reference assembly is - $ReferenceIn"
echo "Prefix - $Prefix"
echo "Aligment file - $AlignmentIn"
echo "OutDir - $OutDir"

CurPath=$PWD
WorkDir=$TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}
mkdir -p $WorkDir
cd $WorkDir

cp $CurPath/$ReferenceIn ./Reference.fasta
cp $CurPath/$AlignmentIn ./Alignment.bam

# ---------------
# Step 2
# Index samples
# ---------------

samtools faidx Reference.fasta
samtools index Alignment.bam

picard CreateSequenceDictionary R=Reference.fasta O=Reference.dict 

# ---------------
# Step 3
# Run HaplotypeCaller
# ---------------

 /scratch/software/gatk4/gatk-4.2.0.0/gatk HaplotypeCaller -ERC GVCF \
     -ploidy 1 \
     -I Alignment.bam \
     -O "$Preix"_SNP_calls.g.vcf \
     -R Reference.fasta

cp $WorkDir/"$Prefix"_SNP_calls.g.vcf $OutDir
rm -r $WorkDir