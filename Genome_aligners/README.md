# Genome alignment tools

1. Star: Spliced transcripts alignment to a reference sequence

2. Bowtie2: fast and sensitive read alignment

3. Bwa: Burrows-Wheeler Aligner for mapping low-divergent sequences against a large reference genom

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


## Bowtie



### Requirements

```bash
conda activate general_tools
conda install bowtie2
# or
# Add path to profile
PATH=${PATH}:/scratch/software/bowtie2/bowtie2-2.3.5.1-linux-x86_64
```

### Typical run

```bash
  for Strain in Strain1 Strain2 Strain3; do
    for Reference in $(ls *contigs_unmasked.fa); do
      for Reads in $(ls -d qc_dna/paired/*/$Strain); do
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

```bash
  for Strain in Ag02 Ag04 Ag06 Ag08; do
  Fusarium=$(ls ../fusarium/assembly/external_group/F.oxysporum_fsp_lycopersici/4287_chromosomal/ensembl/Fusarium_oxysporum_chromosome_and_additional_contigs.fa)
  Stenotrophomonas=$(ls ../../../../../home/armita/prog/deconseq-standalone-0.4.3/database/stenotrophomonas/all_complete_GbStenotrophomonas.fasta)
  Bacteria=$(ls ../../../../../home/armita/prog/deconseq-standalone-0.4.3/database/bacteria/bacteria_split_prinseq.fasta)
  Bacillus=$(ls ../../../../../home/armita/prog/deconseq-standalone-0.4.3/database/bacillus/GbBacillus_split_prinseq.fasta)
  Paenibacillus=$(ls ../../../../../home/armita/prog/deconseq-standalone-0.4.3/database/paenibacillus/paenibacillus_split_prinseq.fasta)
    for Reference in $(ls $Fusarium $Stenotrophomonas $Bacteria $Bacillus $Paenibacillus); do
      for Reads in $(ls -d qc_dna/paired/*/$Strain); do
      echo "$Organism - $Strain"
      F_Read=$(ls $Reads/F/*trim.fq.gz)
      R_Read=$(ls $Reads/R/*trim.fq.gz)
      echo $F_Read
      echo $R_Read
      ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/assemblers/spades
      Organism=$(echo $Reads | rev | cut -f2 -d '/' | rev)
      Strain=$(echo $Reads | rev | cut -f1 -d '/' | rev)
      Prefix=$(basename $Reference | sed 's/.fasta//g' | sed 's/.fa//g' | sed 's/.fna//g')
      OutDir=Thesis/analysis/genome_alignment/bowtie/remove_contaminant/$Organism/$Strain/vs_${Prefix}
      ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/genome_alignment
      sbatch $ProgDir/bowtie/sub_bowtie.sh $Reference $F_Read $R_Read $OutDir
      done
    done
  done
```

## Bwa

### Requirements

```bash
conda activate general_tools
conda install bwa
# or
# Add path to profile
PATH=${PATH}:/scratch/software/bowtie2/bowtie2-2.3.5.1-linux-x86_64
```
### Typical run

```bash
for Reads in $(ls -d path/to/reads/$Organism/$Strain)
  do
    Strain=$(echo $Reads | rev | cut -f1 -d '/' | rev)
    Organism=$(echo $Reads | rev | cut -f2 -d '/' | rev)
      echo "$Organism - $Strain"
      FRead=$Reads/F/*.fq.gz
      RRead=$Reads/R/*.fq.gz
      Reference=path/to/genome/sequence/*.fasta
      OutDir=alignment/bwa/vs_Ref 
      mkdir -p $OutDir
      ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Genome_aligners
      sbatch $ProgDir/bwa.sh $Strain $Reference $FRead $RRead $OutDir
  done
done
```