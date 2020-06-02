# Genome alignment tools

1. Star: Align sequence reads to the reference genome

## Star

Spliced Transcripts Alignment to a Reference. 

### Requirements

```bash
conda install star
```

### Typical run

```bash
  for Assembly in $(ls path/to/unmasked/assembly/*_contigs_unmasked.fa)
  do
    Strain=$(echo $Assembly | rev | cut -f6 -d '/' | rev) # Edit to set your ouput directory
    Organism=$(echo $Assembly | rev | cut -f7 -d '/' | rev) # Edit to set your ouput directory
    echo "$Organism - $Strain"
    for FileF in $(ls path/to/RNAseq/reads/timepoint/F/*_trim.fq.gz)
    do
    FileR=$(echo $FileF | sed 's&/F/&/R/&g'| sed 's/_1_trim/_2_trim/g')
    echo $FileF
    echo $FileR
    Timepoint=$(echo $FileF | rev | cut -d '/' -f3 | rev)
    echo "$Timepoint"
    Sample_Name=$(echo $FileF | rev | cut -d '/' -f1 | rev | sed 's/_1_trim.fq.gz//g')
    OutDir=alignment/star/$Organism/$Strain/$Timepoint/$Sample_Name
    ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_aligners
    sbatch $ProgDir/star.sh $Assembly $FileF $FileR $OutDir
  done
done
```

If multiple RNAseq samples are used, alignment outputs can be concatenated using samtools. 

```bash
  Strain="Strain name"
  Organism="Organism name"
  mkdir -p alignment/star/$Organism/$Strain/$Timepoint/concatenated # e.g.
  samtools merge -f alignment/star/$Organism/$Strain/$Timepoint/concatenated/concatenated.bam \
  alignment/star/$Organism/$Strain/$Timepoint/"Sample_1"/star_aligmentAligned.sorted.out.bam \
  alignment/star/$Organism/$Strain/$Timepoint/"Sample_2"/star_aligmentAligned.sorted.out.bam \
  alignment/star/$Organism/$Strain/$Timepoint/"Sample_3"/star_aligmentAligned.sorted.out.bam
```