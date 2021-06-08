#!/usr/bin/env bash
#SBATCH -J post_snp_calling
#SBATCH --partition=long
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=4

echo "---------------------------------------------------------------"
echo "Filter the SNPs in diploid and haploid genomes. Removes indels. Can allow a certain number of individuals with N sites."
echo "Either input just the VCF file to be filtered or the VCF file with 6 additional filtering options (see below)"
echo "----------------------------------------------------------------"
echo "Arguments: "
echo "**1st argument - input VCF file"
echo "**2nd argument - Minimum QUAL (integer).
QUAL is the phred-scaled probability of a SNP occurring at this site.
NOT AT ALL informative in GATK! Default: 40"
echo "**3rd argument - Minimum MQ (integer).
MQ gives phred-scaled probability that the read is mapped to the correct location
(a low map score will occur if reads are not mapped uniquely at that site).
Default: 30"
echo "**4th argument - Minimum depth (integer) at the site per sample. Default: 10"
echo "**5th argument - Minimum GQ (integer).
GQ is the phred-scaled probability that the sample genotype being called is
correct, given that there is a SNP at that site. Default: 30"
echo "**6th argument - Exclude sites on the basis of the proportion of missing data
(defined to be between 0 and 1, where 0 allows sites that are completely missing
and 1 indicates no missing data allowed). Default: 0.95"
echo "**7th argument - Remove indels or not? Values accepted: Y or N. Default: Y,
as all the standard structure and population genetics analyses focus on SNPs."

vcf=$1

if [ -z "$2" ]
then
mq=40
qual=30
dp=10
gq=30
na=0.95
indel=Y
else
mq=$2
qual=$3
dp=$4
gq=$5
na=$6
indel=$7
fi

filename=$(basename "$vcf")
vcffilter -f "QUAL > $qual & MQ > $mq" $vcf | vcffilter -g "DP > $dp & GQ > $gq" >temp.vcf
echo $indel
if [ "$indel" = "Y" ]
then
    vcftools --vcf temp.vcf --max-missing $na --remove-indels --recode --out ${filename%.vcf}_filtered
elif [ "$indel" = "N" ]
then
    vcftools --vcf temp.vcf --max-missing $na --recode --out ${filename%.vcf}_filtered
else
    echo "Wrong last argument!" 1>&2
	exit 1
fi
mv ${filename%.vcf}_filtered.recode.vcf ${filename%.vcf}_filtered.vcf
#rm *temp*.vcf
