# SNP calling 

## Requirements

```bash
conda activate SNP_calling

conda install picard
# version 0.1.18 needed
conda install samtools=0.1.18
```

## Rename input mapping files in each folder by prefixing with the strain ID


```bash
for filename in $(ls -d analysis/genome_alignment/bowtie/$Organism/$Strain/*); do
Organism=$(echo $filename | rev | cut -f3 -d '/' | rev)
Strain=$(echo $filename | rev | cut -f2 -d '/' | rev)
echo "$Organism - $Strain"
  mv "$filename/*_contigs_unmasked.fa_aligned.sam" "$filename/"$Strain"_unmasked.fa_aligned.sam"
  mv "$filename/*_contigs_unmasked.fa_aligned.bam" "$filename/"$Strain"_unmasked.fa_aligned.bam"
  mv "$filename/*_contigs_unmasked.fa_aligned_sorted.bam" "$filename/"$Strain"_unmasked.fa_aligned_sorted.bam"
  mv "$filename/*_contigs_unmasked.fa_aligned_sorted.bam.index" "$filename/"$Strain"_unmasked.fa_aligned_sorted.bam.index"
done
```

## Remove multimapping reads, discordant reads. PCR and optical duplicates, and add read group and sample name to each mapped read (preferably, the shortest ID possible)

```bash
  for Strain in Strain1 Strain2; do # Replace with the strain name
    for input in analysis/genome_alignment/bowtie/$Organism/$Strain/"$Strain"_unmasked.fa_aligned.sam; do
    Jobs=$(squeue -u ${USER} --noheader --array | wc -l)
        while [ $Jobs -gt 5 ]; do
            sleep 5m
            printf "."
            Jobs=$(squeue -u ${USER} --noheader --array | wc -l)
        done
    printf "\n"
    OutDir=analysis/genome_alignment/bowtie/$Organism/$Strain
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/SNP_calling
    sbatch $ProgDir/pre_SNP_calling.sh $input $Strain $OutDir
    done
  done
```