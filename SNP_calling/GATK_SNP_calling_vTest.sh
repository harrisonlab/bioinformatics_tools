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

input=Nd
reference=R0905_contigs_unmasked.fa

filename=$(basename "$reference")
output=analysis/popgen/SNP_calling_R0905/"${filename%.*}_temp.vcf"
output2=analysis/popgen/SNP_calling_R0905/"${filename%.*}.vcf"

gatk=/scratch/software/GenomeAnalysisTK-3.6

java -jar $gatk/GenomeAnalysisTK.jar \
     -R $reference \
     -T HaplotypeCaller \
     -ploidy 1 \
     -nct 30 \
     --allow_potentially_misencoded_quality_scores \
     -I $input/Ag02/Ag02_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/Ag04/Ag04_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/Ag06/Ag06_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -o $output

#Break down complex SNPs into primitive ones with VariantsToAllelicPrimitives
#This tool will take an MNP (e.g. ACCCA -> TCCCG) and break it up into separate records for each component part (A-T and A->G).
#This tool modifies only bi-allelic variants.

java -jar $gatk/GenomeAnalysisTK.jar \
   -T VariantsToAllelicPrimitives \
   -R $reference \
   -V $output \
   -o $output2 \
