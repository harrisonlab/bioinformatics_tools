#!/usr/bin/env bash
#SBATCH -J GenotypeGVCFs
#SBATCH --partition=himem
#SBATCH --mem-per-cpu=12G
#SBATCH --cpus-per-task=30

# Perform joint genotyping on one or more samples pre-called with HaplotypeCaller

Input_vcf=$1
Reference=$2
OutDir=$3

filename=$(basename "$Reference")
output=$input/"${filename%.*}_genotyped.vcf.gz"
#output2=$input/"${filename%.*}.vcf"

/scratch/software/gatk4/gatk-4.2.0.0/gatk GenotypeGVCFs \
	-R $Reference \
	-V $Input_vcf \
	-O $OutDir/$output	
 