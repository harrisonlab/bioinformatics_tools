# SNP calling 

## Requirements

```bash
conda activate SNP_calling

conda install picard

conda install perl-vcftools-vcf
conda install vcftools
conda install vcflib
```

```bash
for Strain in Ag02 Ag04 Ag06 Ag08; do
for Reference in $(ls R0905_contigs_unmasked.fa); do
for Reads in $(ls -d ../../../qc_dna/paired/N.ditissima/$Strain); do
echo "$Organism - $Strain"
F_Read=$(ls $Reads/F/*trim.fq.gz)
R_Read=$(ls $Reads/R/*trim.fq.gz)
echo $F_Read
echo $R_Read
Organism=$(echo $Reads | rev | cut -f2 -d '/' | rev)
Strain=$(echo $Reads | rev | cut -f1 -d '/' | rev)
Prefix=$(basename $Reference | sed 's/.fasta//g' | sed 's/.fa//g' | sed 's/.fna//g')
OutDir=$Strain/vs_${Prefix}
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_aligners
sbatch $ProgDir/bowtie.sh $Reference $F_Read $R_Read $OutDir
done
done
done
```

## Rename input mapping files in each folder by prefixing with the strain ID.

```bash
  for filename in $(ls -d analysis/genome_alignment/bowtie/vs_*/$Organism/$Strain); do
  Organism=$(echo $filename | rev | cut -f3 -d '/' | rev)
  Strain=$(echo $filename | rev | cut -f2 -d '/' | rev)
  echo "$Organism - $Strain"
    mv "$filename/*_contigs_unmasked.fa_aligned.sam" "$filename/"$Strain"_unmasked.fa_aligned.sam"
    mv "$filename/*_contigs_unmasked.fa_aligned.bam" "$filename/"$Strain"_unmasked.fa_aligned.bam"
    mv "$filename/*_contigs_unmasked.fa_aligned_sorted.bam" "$filename/"$Strain"_unmasked.fa_aligned_sorted.bam"
    mv "$filename/*_contigs_unmasked.fa_aligned_sorted.bam.index" "$filename/"$Strain"_unmasked.fa_aligned_sorted.bam.index"
  done
```



## Remove multimapping reads, discordant reads. PCR and optical duplicates. 

```bash
for Strain in Ag02 Ag04 Ag06 Ag08; do
for input in Nd/$Strain/"$Strain"_unmasked.fa_aligned.sam; do
Jobs=$(squeue -u ${USER} --noheader --array | wc -l)
while [ $Jobs -gt 5 ]; do
sleep 5m
printf "."
Jobs=$(squeue -u ${USER} --noheader --array | wc -l)
done
printf "\n"
OutDir=Nd/$Strain
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/SNP_calling
sbatch $ProgDir/pre_SNP_calling.sh $input $Strain $OutDir
done
done
```

## Prepare genome reference indexes required by GATK

```bash
reference=R0905_contigs_unmasked.fa
input=./
filename=$(basename "$reference")
output="${filename%.*}.dict"
picard CreateSequenceDictionary R=$reference O=$input/$output
samtools faidx $reference
```

### Copy index file to same folder as BAM alignments

```bash
for Strain in Ag02 Ag04 Ag06; do
Index=R0905_contigs_unmasked.fa.fai
Directory=Nd/$Strain
cp $Index $Directory
done
```

## Start SNP calling with GATK
The submission script required need to be custom-prepared for each analysis, depending on what samples are being analysed.
See GATK_SNP_calling.sh

```bash
mkdir -p analysis/popgen/SNP_calling
cd analysis/popgen/SNP_calling
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/SNP_calling
sbatch $ProgDir/GATK_SNP_calling_vTest.sh
```

# SNP calling analysis


### Filter vcf outputs, only retain biallelic high-quality SNPS with no missing data for genetic analyses.

```bash
  for vcf in $(ls *_contigs_unmasked.vcf)
  do
      echo $vcf
      ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/SNP_calling
      sbatch $ProgDir/vcf_parser.sh $vcf
  done
```