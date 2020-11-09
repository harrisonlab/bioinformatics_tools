# SNP calling 

## Requirements

```bash
conda activate SNP_calling

conda install picard

conda install perl-vcftools-vcf
conda install vcftools
conda install vcflib
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
  for Strain in Strain1 Strain2; do # Replace with the strain name
    for input in analysis/genome_alignment/bowtie/vs_*/$Organism/$Strain/"$Strain"_unmasked.fa_aligned.sam; do
    Jobs=$(squeue -u ${USER} --noheader --array | wc -l)
        while [ $Jobs -gt 5 ]; do
            sleep 5m
            printf "."
            Jobs=$(squeue -u ${USER} --noheader --array | wc -l)
        done
    printf "\n"
    OutDir=analysis/genome_alignment/bowtie/$Organism/$Strain
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/SNP_calling
    sbatch $ProgDir/pre_SNP_calling.sh $input $Strain $OutDir # This will add read group and sample name to each mapped read. Preferably, use the shortest ID possible.
    done
  done
```

## Prepare genome reference indexes required by GATK

```bash
reference=repeat_masked/filtered_contigs/R0905_contigs_unmasked.fa
input=repeat_masked/filtered_contigs/
filename=$(basename "$reference")
output="${filename%.*}.dict"
picard CreateSequenceDictionary R=$reference O=$input/$output
samtools faidx $reference
```

### Copy index file to same folder as BAM alignments

```bash
for Strain in Ag02 Ag04 Ag05 Ag06 Ag08; do
  Index=repeat_masked/filtered_contigs/R0905_contigs_unmasked.fa.fai
  Directory=analysis/genome_alignment/bowtie/$Organism/$Strain
  cp $Index $Directory
done
```

## Start SNP calling with GATK
The submission script required need to be custom-prepared for each analysis, depending on what samples are being analysed.
See GATK_SNP_calling.sh

```bash
mkdir -p analysis/popgen/SNP_calling
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/SNP_calling
sbatch $ProgDir/GATK_SNP_calling.sh
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


Only one of these can be run at a time!

### General VCF stats

```bash
cd ../

vcfstats analysis/popgen/SNP_calling/R0905_contigs_unmasked.vcf > analysis/popgen/SNP_calling/R0905_contigs_unmasked.stat
perl /home/sobczm/bin/vcftools/bin/vcf-stats \
SNP_calling3/Hg199_contigs_unmasked.vcf > SNP_calling3/Hg199_contigs_unmasked.stat
perl /home/sobczm/bin/vcftools/bin/vcf-stats \
SNP_calling3/Hg199_contigs_unmasked_filtered.vcf > SNP_calling3/Hg199_contigs_unmasked_filtered.stat

perl /home/sobczm/bin/vcftools/bin/vcf-stats \
SNP_calling_R0905/R0905_good_contigs_unmasked.vcf > SNP_calling_R0905/R0905_contigs_unmasked.stat
perl /home/sobczm/bin/vcftools/bin/vcf-stats \
SNP_calling_R0905/R0905_good_contigs_unmasked_filtered.vcf > SNP_calling_R0905/R0905_contigs_unmasked_filtered.stat
```

#Calculate the index for percentage of shared SNP alleles between the individuals.

```bash
for vcf in $(ls SNP_calling3/*_filtered.vcf)
do
  scripts=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/Popgen_analysis/snp
  echo $vcf
  $scripts/similarity_percentage.py $vcf
done

for vcf in $(ls SNP_calling_R0905/*_filtered.vcf)
do
  scripts=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/Popgen_analysis/snp
  echo $vcf
  $scripts/similarity_percentage.py $vcf
done
```

#Remove monomorphic sites (minor allele count minimum 1). Argument --vcf is the filtered VCF file, and --out is the suffix to be used for the output file.

```bash
for vcf in $(ls SNP_calling3/*_filtered.vcf)
do
    echo $vcf
    out=$(basename $vcf .vcf)
    echo $out
    $vcftools/vcftools --vcf $vcf --mac 1 --recode --out SNP_calling3/$out
done

#After filtering, kept 593720 out of a possible 758934 Sites

for vcf in $(ls SNP_calling_R0905/*_filtered.vcf)
do
    echo $vcf
    out=$(basename $vcf .vcf)
    echo $out
    $vcftools/vcftools --vcf $vcf --mac 1 --recode --out SNP_calling_R0905/$out
done

#After filtering, kept 626425 out of a possible 800178 Sites
```

#Create custom SnpEff genome database
```bash
snpeff=/home/gomeza/prog/snpEff
nano $snpeff/snpEff.config
```

##Add the following lines to the section with databases:

#---
# EMR Databases
#----
# Hg199 genome
Hg199_Ref.genome: Hg199_Ref

#Collect input files

```bash
mkdir -p $snpeff/data/Hg199_Ref
cp repeat_masked/Ref_Genomes/N.ditissima/Hg199/filtered_contigs/Hg199_contigs_unmasked.fa $snpeff/data/Hg199_Ref
cp gene_pred/codingquary/Ref_Genomes/N.ditissima/Hg199/final/final_genes_appended_renamed.gff3 $snpeff/data/Hg199_Ref
```
```bash
mkdir -p $snpeff/data/R0905
cp R0905_good/repeat_masked/filtered_contigs/R0905_good_contigs_unmasked.fa $snpeff/data/R0905
cp gene_pred/codingquary/Ref_Genomes_v2/N.ditissima/R0905/final/final_genes_appended_renamed.gff3 $snpeff/data/R0905
```
#Rename input files
```bash
cd $snpeff/data/Hg199_Ref
mv final_genes_appended_renamed.gff3 genes.gff
mv Hg199_contigs_unmasked.fa sequences.fa
```
```bash
cd $snpeff/data/R0905
mv final_genes_appended_renamed.gff3 genes.gff
mv R0905_good_contigs_unmasked.fa sequences.fa
```
#Build database using GFF3 annotation

java -jar $snpeff/snpEff.jar build -gff3 -v Hg199_Ref

java -jar $snpeff/snpEff.jar build -gff3 -v R0905

#Annotate VCF files

```bash
input=/data/scratch/gomeza/analysis/popgen/
cd $input
for a in SNP_calling3/*recode.vcf
do
    echo $a
    filename=$(basename "$a")
    java -Xmx4g -jar $snpeff/snpEff.jar -v -ud 0 Hg199_Ref $a > ${filename%.vcf}_annotated.vcf
    mv snpEff_genes.txt SNP_calling3/snpEff_genes_${filename%.vcf}.txt
    mv snpEff_summary.html  SNP_calling3/snpEff_summary__${filename%.vcf}.html
done
```
```bash
input=/data/scratch/gomeza/analysis/popgen/
cd $input
for a in SNP_calling_R0905/*recode.vcf
do
    echo $a
    filename=$(basename "$a")
    java -Xmx4g -jar $snpeff/snpEff.jar -v -ud 0 R0905 $a > ${filename%.vcf}_annotated.vcf
    mv snpEff_genes.txt SNP_calling_R0905/snpEff_genes_${filename%.vcf}.txt
    mv snpEff_summary.html  SNP_calling_R0905/snpEff_summary__${filename%.vcf}.html
done
```


#Visualise the output as heatmap and clustering dendrogram

```bash
for log in $(ls analysis/popgen/SNP_calling3/*distance.log)
do
  scripts=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/Popgen_analysis/snp
  Rscript --vanilla $scripts/distance_matrix.R $log
done
```
```bash
for log in $(ls analysis/popgen/SNP_calling_R0905/*distance.log)
do
  scripts=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/Popgen_analysis/snp
  Rscript --vanilla $scripts/distance_matrix.R $log
done
```
#Carry out PCA and plot the results

```bash
for vcf in $(ls analysis/popgen/SNP_calling3/*filtered.vcf)
do
    echo $vcf
    scripts=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/Popgen_analysis/snp
    out=$(basename $vcf contigs_unmasked_filtered.vcf)
    echo $out
    Rscript --vanilla $scripts/pca.R $vcf $out
done
```
```bash
for vcf in $(ls analysis/popgen/SNP_calling_R0905/*filtered.vcf)
do
    echo $vcf
    scripts=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/Popgen_analysis/snp
    out=$(basename $vcf contigs_unmasked_filtered.vcf)
    echo $out
    Rscript --vanilla $scripts/pca.R $vcf $out
done
```
#Calculate an NJ tree based on all the SNPs. Outputs a basic display of the tree, plus a Newick file to be used for displaying the tree in FigTree and beautifying it.

```bash
cd analysis/popgen/SNP_calling3
for vcf in $(ls *filtered.vcf)
do
    echo $vcf
    scripts=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/Popgen_analysis/snp
    $scripts/nj_tree.sh $vcf 1
done
```
```bash
cd analysis/popgen/SNP_calling_R0905
for vcf in $(ls *filtered.vcf)
do
    echo $vcf
    scripts=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/Popgen_analysis/snp
    $scripts/nj_tree.sh $vcf 1
done
```
This script isn't happy with the low levels of variation I have. TODO: Investigate other tree methods - Michelle mentioned some possibilities


#Remove low coverage samples

```bash
#Remove outgroup and low-coverage isolates. Create a cut-down VCF and filter it
#The isolate Ag11_B has a 22X coverage, so it will be removed and new statistics calculated.
vcflib=/home/sobczm/bin/vcflib/bin
cd analysis/popgen/SNP_calling_R0905
mkdir NoAg11_B/
cp *.vcf NoAg11_B/
$vcflib/vcfremovesamples R0905_good_contigs_unmasked.vcf Ag11_B > R0905_good_contigs_unmasked_NoAg11B.vcf
$vcflib/vcfremovesamples R0905_good_contigs_unmasked_filtered.vcf Ag11_B > R0905_good_contigs_unmasked_NoAg11B_filtered.vcf
```

#General VCF stats (remember that vcftools needs to have the PERL library exported)

```bash
cd analysis/popgen/SNP_calling_R0905/NoAg11_B

perl /home/sobczm/bin/vcftools/bin/vcf-stats \
R0905_good_contigs_unmasked_NoAg11B.vcf > R0905_contigs_unmasked_NoAg11B.stat
perl /home/sobczm/bin/vcftools/bin/vcf-stats \
R0905_good_contigs_unmasked_NoAg11B_filtered.vcf > R0905_contigs_unmasked_filtered_NoAg11B.stat
```
#Calculate the index for percentage of shared SNP alleles between the individuals.

```bash
for vcf in $(ls R0905_good_contigs_unmasked_NoAg11B_filtered.vcf)
do
  scripts=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/Popgen_analysis/snp
  echo $vcf
  $scripts/similarity_percentage.py $vcf
done
```

#Remove monomorphic sites (minor allele count minimum 1). Argument --vcf is the filtered VCF file, and --out is the suffix to be used for the output file.

```bash
for vcf in $(ls R0905_good_contigs_unmasked_NoAg11B_filtered.vcf)
do
    echo $vcf
    out=$(basename $vcf .vcf)
    echo $out
    $vcftools/vcftools --vcf $vcf --mac 1 --recode --out $out
done

After filtering, kept 248598 out of a possible 800178 Sites
```

#Visualise the output as heatmap and clustering dendrogram

```bash
for log in $(ls *distance.log)
do
  scripts=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/Popgen_analysis/snp
  Rscript --vanilla $scripts/distance_matrix.R $log
done
```

#Carry out PCA and plot the results

```bash
for vcf in $(ls R0905_good_contigs_unmasked_NoAg11B_filtered.vcf)
do
    echo $vcf
    scripts=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/Popgen_analysis/snp
    out=$(basename $vcf contigs_unmasked_NoAg11B_filtered.vcf)
    echo $out
    Rscript --vanilla $scripts/pca.R $vcf $out
done
```



#Calculate an NJ tree based on all the SNPs. Outputs a basic display of the tree, plus a Newick file to be used for displaying the tree in FigTree and beautifying it.

```bash
  for vcf in $(ls R0905_good_contigs_unmasked_NoAg11B_filtered.vcf)
  do
      echo $vcf
      scripts=/home/gomeza/git_repos/emr_repos/scripts/neonectria_ditissima/Popgen_analysis/snp
      $scripts/nj_tree.sh $vcf 1
  done
```
This script isn't happy with the low levels of variation I have. TODO: Investigate other tree methods - Michelle mentioned some possibilities
