# SNP calling 

## Requirements

```bash
conda activate SNP_calling

conda install picard
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
cd analysis/popgen/SNP_calling
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/SNP_calling
sbatch $ProgDir/GATK_SNP_calling.sh
```

