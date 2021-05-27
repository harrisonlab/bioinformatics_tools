#!/usr/bin/env bash

#Input: VCF file containing bi-allelic SNPs. Output: a neighbour joining tree, suffix: _nj.nwk, and PDF file with the tree, suffix: _nj.pdf
vcf_file=$1
ploidy=$2
vcf_tools=/scratch/software/vcftools-0.1.15/bin
scripts=/home/gomeza/git_repos/scripts/bioinformatics_tools/SNP_calling

gzip $vcf_file
zcat "$vcf_file".gz | $vcf_tools/vcf-to-tab > "${vcf_file%.vcf}.tab"
gzip -d "$vcf_file".gz
perl $scripts/vcf_tab_to_fasta_alignment.pl -i "${vcf_file%.vcf}.tab" >"${vcf_file%.vcf}.fasta"

if [ "$ploidy" = "1" ]
then
    Rscript --vanilla $scripts/nj.R "${vcf_file%.vcf}.fasta"
elif [ "$ploidy" = "2" ]
then
    seqtk=/scratch/software/seqtk/seqtk
    #Pick one random allele at heterozygous sites to obtain haploid sequences
    $seqtk randbase "${vcf_file%.vcf}.fasta" >"${vcf_file%.vcf}_haploidised.fasta"
    Rscript --vanilla $scripts/nj.R "${vcf_file%.vcf}_haploidised.fasta"
else
    exit
fi
