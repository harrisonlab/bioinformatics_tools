#!/usr/bin/env bash
#SBATCH -J snp_calling
#SBATCH --partition=himem
#SBATCH --mem-per-cpu=12G
#SBATCH --cpus-per-task=30

# NOTE: this is a haploid organism. For diploid organism, change "ploidy" argument to 2.
# Changes required in the script:
# VARIABLES
# Reference - the genome reference used in read mapping.
# INSIDE THE GATK command:
# To specify which BAM mapping files (output from pre_SNP_calling.sh, filename ending with "_rg" -> that is, with
# read group added) are to be used in SNP calling, use the -I argument with full path to each file following after that.
# Each new BAM file has to be specified after a separate -I

input=Home/analysis_vAG/genome_alignment/bowtie/vs_R0905
reference=R0905_good/repeat_masked/filtered_contigs/R0905_good_contigs_unmasked.fa

filename=$(basename "$reference")
output=Home/analysis_vAG/SNPs/SNP_calling_R0905_VP/"${filename%.*}_temp.vcf"
output2=Home/analysis_vAG/SNPs/SNP_calling_R0905_VP/"${filename%.*}.vcf"

# gatk=/scratch/software/gatk4/gatk-4.1.9.0

     /scratch/software/gatk4/gatk-4.1.9.0/gatk HaplotypeCaller \
     -R $reference \
     -ploidy 1 \
     -I $input/*/*/Ag02/Ag02_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/*/Ag04/Ag04_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/*/Ag05/Ag05_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/*/Ag06/Ag06_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/*/Ag08/Ag08_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/*/Ag09_A/Ag09_A_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/*/Ag11_A/Ag11_A_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/*/Ag11_B/Ag11_B_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/*/Ag11_C/Ag11_C_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/*/BGV344/BGV344_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/*/Hg199/Hg199_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/*/ND8/ND8_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/*/ND9/ND9_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/*/OPC304/OPC304_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/*/P112/P112_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/*/R0905/R0905_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/*/R37-15/R37-15_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/*/R39-15/R39-15_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/*/R41-15/R41-15_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/*/R42-15/R42-15_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/*/R45-15/R45-15_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/*/R6-17-2/R6-17-2_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/*/R6-17-3/R6-17-3_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/*/R68-17-C2/R68-17-C2_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/*/R68-17-C3/R68-17-C3_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/*/SVK1/SVK1_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/*/SVK2/SVK2_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/*/118923/118923_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/*/118924/118924_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/*/226-31/226-31_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/*/227-31/227-31_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/*/*/NMaj/NMaj_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -O $output

#Break down complex SNPs into primitive ones with VariantsToAllelicPrimitives
#This tool will take an MNP (e.g. ACCCA -> TCCCG) and break it up into separate records for each component part (A-T and A->G).
#This tool modifies only bi-allelic variants.

# /scratch/software/gatk4/gatk-4.1.9.0/gatk VariantsToAllelicPrimitives \
#    -R $reference \
#    -V $output \
#    -o $output2 \
