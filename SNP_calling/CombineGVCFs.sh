#!/usr/bin/env bash
#SBATCH -J CombineGVCFs
#SBATCH --partition=himem
#SBATCH --mem-per-cpu=12G
#SBATCH --cpus-per-task=30

# Merges one or more HaplotypeCaller GVCF files into a single GVCF with appropriate annotations

input=analysis_VP/SNP_calling/HaplotypeCaller
reference=$input/R0905_good_contigs_unmasked.fa

filename=$(basename "$reference")
output=$input/"${filename%.*}_cohort.g.vcf.gz"
output2=$input/"${filename%.*}.vcf"

/scratch/software/gatk4/gatk-4.2.0.0/gatk CombineGVCFs \
	-R $reference \
	-V $input/118923_SNP_calls.g.vcf \
	-V $input/118924_SNP_calls.g.vcf \
	-V $input/226-31_SNP_calls.g.vcf \
	-V $input/227-31_SNP_calls.g.vcf \
	-V $input/Ag02_SNP_calls.g.vcf \
	-V $input/Ag04_SNP_calls.g.vcf \
	-V $input/Ag05_SNP_calls.g.vcf \
	-V $input/Ag06_SNP_calls.g.vcf \
	-V $input/Ag08_SNP_calls.g.vcf \
	-V $input/Ag09_A_SNP_calls.g.vcf \
	-V $input/Ag11_A_SNP_calls.g.vcf \
	-V $input/Ag11_B_SNP_calls.g.vcf \
	-V $input/Ag11_C_SNP_calls.g.vcf \
	-V $input/BGV344_SNP_calls.g.vcf \
	-V $input/Hg199_SNP_calls.g.vcf \
	-V $input/ND8_SNP_calls.g.vcf \
	-V $input/ND9_SNP_calls.g.vcf \
	-V $input/NMaj_SNP_calls.g.vcf \
	-V $input/OPC304_SNP_calls.g.vcf \
	-V $input/P112_SNP_calls.g.vcf \
	-V $input/R0905_SNP_calls.g.vcf \
	-V $input/R37-15_SNP_calls.g.vcf \
	-V $input/R39-15_SNP_calls.g.vcf \
	-V $input/R41-15_SNP_calls.g.vcf \
	-V $input/R42-15_SNP_calls.g.vcf \
	-V $input/R45-15_SNP_calls.g.vcf \
	-V $input/R6-17-2_SNP_calls.g.vcf \
	-V $input/R6-17-3_SNP_calls.g.vcf \
	-V $input/R68-17-C2_SNP_calls.g.vcf \
	-V $input/R68-17-C3_SNP_calls.g.vcf \
	-V $input/SVK1_SNP_calls.g.vcf \
	-V $input/SVK2_SNP_calls.g.vcf \
	-O $output	
 